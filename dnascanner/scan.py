#!/usr/bin/env python3
"""
DNAscanner — extract variants from a VCF, match against a local SNPedia
database, and store results in a SQLite database for report generation.

Incremental: caches VCF extraction (skipped if file unchanged),
re-joins against SNPedia when its data is updated.

Usage:
    uv run scan --vcf /path/to/file.vcf.gz --snpedia data/snpedia.db --out results.db
"""

import argparse
import gzip
import hashlib
import os
import re
import sqlite3
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

from tqdm import tqdm


# ---------------------------------------------------------------------------
# Results database schema
# ---------------------------------------------------------------------------

RESULTS_SCHEMA = """
CREATE TABLE IF NOT EXISTS vcf_variants (
    chromosome      TEXT NOT NULL,
    position        INTEGER NOT NULL,
    ref             TEXT NOT NULL,
    alt             TEXT NOT NULL,
    genotype        TEXT NOT NULL,
    quality         REAL,
    filter_status   TEXT,
    vcf_rsids       TEXT,
    PRIMARY KEY (chromosome, position, alt)
);

CREATE TABLE IF NOT EXISTS matches (
    chromosome      TEXT NOT NULL,
    position        INTEGER NOT NULL,
    ref             TEXT NOT NULL,
    alt             TEXT NOT NULL,
    genotype        TEXT NOT NULL,
    quality         REAL,
    filter_status   TEXT,
    rsid            TEXT NOT NULL,
    gene            TEXT,
    snp_summary     TEXT,
    geno_alleles    TEXT,
    geno_magnitude  REAL,
    geno_repute     TEXT,
    geno_summary    TEXT,
    match_method    TEXT,
    snpedia_rev     INTEGER,
    matched_at      TEXT,
    PRIMARY KEY (chromosome, position, alt)
);

CREATE TABLE IF NOT EXISTS scan_meta (
    key             TEXT PRIMARY KEY,
    value           TEXT
);

CREATE INDEX IF NOT EXISTS idx_matches_rsid ON matches(rsid);
CREATE INDEX IF NOT EXISTS idx_matches_magnitude ON matches(geno_magnitude);
CREATE INDEX IF NOT EXISTS idx_matches_gene ON matches(gene);
CREATE INDEX IF NOT EXISTS idx_variants_pos ON vcf_variants(chromosome, position);
"""


def init_results_db(db_path: str) -> sqlite3.Connection:
    conn = sqlite3.connect(db_path)
    conn.execute("PRAGMA journal_mode=WAL")
    conn.execute("PRAGMA synchronous=NORMAL")
    conn.executescript(RESULTS_SCHEMA)
    conn.commit()
    return conn


# ---------------------------------------------------------------------------
# VCF extraction
# ---------------------------------------------------------------------------

def normalize_chrom(chrom: str) -> str:
    """Strip 'chr' prefix and uppercase for matching against SNPedia."""
    c = chrom.upper()
    if c.startswith("CHR"):
        c = c[3:]
    if c == "M":
        c = "MT"
    return c


def extract_rsids_from_info(info: str) -> list[str]:
    """Extract rsIDs from VEP CSQ or snpEff ANN fields."""
    rsids = set()
    for m in re.finditer(r"\brs(\d+)\b", info, re.IGNORECASE):
        rsids.add(f"rs{m.group(1)}")
    return list(rsids)


def genotype_alleles(gt_str: str, ref: str, alts: list[str]) -> tuple[str, str] | None:
    """Convert VCF GT string to a pair of allele bases."""
    if "." in gt_str:
        return None
    sep = "/" if "/" in gt_str else "|"
    indices = gt_str.split(sep)
    if len(indices) != 2:
        return None
    all_alleles = [ref] + alts
    try:
        a1 = all_alleles[int(indices[0])]
        a2 = all_alleles[int(indices[1])]
    except (ValueError, IndexError):
        return None
    return (a1, a2)


def file_fingerprint(path: str) -> str:
    """Return a fingerprint string for change detection: size + mtime."""
    st = os.stat(path)
    return f"{st.st_size}:{st.st_mtime}"


def vcf_needs_extraction(results_conn: sqlite3.Connection, vcf_path: str) -> bool:
    """Check if VCF extraction can be skipped (file unchanged)."""
    row = results_conn.execute(
        "SELECT value FROM scan_meta WHERE key = 'vcf_fingerprint'"
    ).fetchone()
    if row is None:
        return True
    return row[0] != file_fingerprint(vcf_path)


def extract_vcf_variants(vcf_path: str, results_conn: sqlite3.Connection):
    """Stream VCF and extract all SNP variants into vcf_variants table."""
    file_size = os.path.getsize(vcf_path)
    is_gzipped = vcf_path.endswith(".gz")
    opener = gzip.open if is_gzipped else open

    # Clear previous extraction
    results_conn.execute("DELETE FROM vcf_variants")
    results_conn.commit()

    batch = []
    batch_size = 5000

    with opener(vcf_path, "rt", errors="replace") as f:
        pbar = tqdm(
            total=file_size,
            unit="B",
            unit_scale=True,
            unit_divisor=1024,
            desc="Extracting VCF",
            bar_format="{l_bar}{bar}| {n_fmt}/{total_fmt} [{elapsed}<{remaining}, {rate_fmt}] {postfix}",
            colour="#38bdf8",
        )
        last_raw_pos = 0
        extracted = 0
        total_variants = 0

        for line in f:
            # Update progress
            if is_gzipped and hasattr(f, "fileobj"):
                raw_pos = f.fileobj.tell()
                if raw_pos != last_raw_pos:
                    pbar.update(raw_pos - last_raw_pos)
                    last_raw_pos = raw_pos
            else:
                pbar.update(len(line))

            if line.startswith("#"):
                continue

            total_variants += 1

            fields = line.rstrip("\n").split("\t")
            if len(fields) < 10:
                continue

            chrom_raw = fields[0]
            try:
                pos = int(fields[1])
            except ValueError:
                continue
            vcf_id = fields[2]
            ref = fields[3]
            alt_str = fields[4]
            qual = fields[5]
            filt = fields[6]
            info = fields[7]
            fmt = fields[8]
            sample = fields[9]

            alts = alt_str.split(",")

            # Skip indels
            if len(ref) > 1 or any(len(a) > 1 for a in alts):
                continue

            # Parse genotype
            fmt_fields = fmt.split(":")
            sample_fields = sample.split(":")
            gt_idx = fmt_fields.index("GT") if "GT" in fmt_fields else -1
            if gt_idx < 0 or gt_idx >= len(sample_fields):
                continue

            gt_str = sample_fields[gt_idx]
            allele_pair = genotype_alleles(gt_str, ref, alts)
            if allele_pair is None:
                continue

            a1, a2 = allele_pair
            a1_s, a2_s = sorted([a1.upper(), a2.upper()])
            genotype = f"({a1_s};{a2_s})"

            chrom_norm = normalize_chrom(chrom_raw)

            try:
                qual_f = float(qual) if qual != "." else None
            except ValueError:
                qual_f = None

            # Extract rsIDs from ID column and INFO
            rsids = []
            if vcf_id != "." and vcf_id.lower().startswith("rs"):
                rsids.append(vcf_id.lower())
            rsids.extend(r.lower() for r in extract_rsids_from_info(info))
            rsid_str = ",".join(rsids) if rsids else None

            batch.append((chrom_norm, pos, ref, alt_str, genotype, qual_f, filt, rsid_str))
            extracted += 1

            if len(batch) >= batch_size:
                results_conn.executemany(
                    "INSERT OR REPLACE INTO vcf_variants VALUES (?, ?, ?, ?, ?, ?, ?, ?)",
                    batch,
                )
                results_conn.commit()
                batch.clear()
                pbar.set_postfix_str(f"extracted={extracted:,}", refresh=False)

        # Flush remaining
        if batch:
            results_conn.executemany(
                "INSERT OR REPLACE INTO vcf_variants VALUES (?, ?, ?, ?, ?, ?, ?, ?)",
                batch,
            )

        pbar.close()

    # Save fingerprint and stats
    now = datetime.now(timezone.utc).isoformat()
    meta = [
        ("vcf_fingerprint", file_fingerprint(vcf_path)),
        ("vcf_path", vcf_path),
        ("vcf_total_variants", str(total_variants)),
        ("vcf_extracted_snps", str(extracted)),
        ("vcf_extracted_at", now),
    ]
    results_conn.executemany(
        "INSERT OR REPLACE INTO scan_meta (key, value) VALUES (?, ?)", meta
    )
    results_conn.commit()

    tqdm.write(f"  {total_variants:,} total variants, {extracted:,} SNPs extracted")


# ---------------------------------------------------------------------------
# SNPedia matching
# ---------------------------------------------------------------------------

def snpedia_fingerprint(snpedia_conn: sqlite3.Connection) -> str:
    """Fingerprint based on row counts and max revision IDs."""
    snp_count = snpedia_conn.execute("SELECT COUNT(*) FROM snps").fetchone()[0]
    geno_count = snpedia_conn.execute("SELECT COUNT(*) FROM genotypes").fetchone()[0]
    snp_max_rev = snpedia_conn.execute("SELECT MAX(revision_id) FROM snps").fetchone()[0] or 0
    geno_max_rev = snpedia_conn.execute("SELECT MAX(revision_id) FROM genotypes").fetchone()[0] or 0
    return f"snps={snp_count}:{snp_max_rev},genos={geno_count}:{geno_max_rev}"


def matching_needs_update(results_conn: sqlite3.Connection, snpedia_conn: sqlite3.Connection) -> bool:
    """Check if matches need recomputation (SNPedia data changed)."""
    row = results_conn.execute(
        "SELECT value FROM scan_meta WHERE key = 'snpedia_fingerprint'"
    ).fetchone()
    if row is None:
        return True
    return row[0] != snpedia_fingerprint(snpedia_conn)


def match_variants(results_conn: sqlite3.Connection, snpedia_conn: sqlite3.Connection):
    """Join extracted VCF variants against SNPedia and populate matches table."""

    # Load SNPedia lookups into memory
    tqdm.write("  Loading SNPedia data...")

    # Position lookup: (chrom, pos) -> (rsid, gene, summary, revision_id)
    pos_lookup: dict[tuple[str, int], tuple] = {}
    for row in snpedia_conn.execute("""
        SELECT rsid, gene, chromosome, position, summary, revision_id
        FROM snps WHERE chromosome IS NOT NULL AND position IS NOT NULL
    """):
        rsid, gene, chrom, pos, summary, rev = row
        pos_lookup[(chrom.strip().upper(), pos)] = (rsid, gene, summary, rev)

    # rsID lookup
    rsid_lookup: dict[str, tuple] = {}
    for row in snpedia_conn.execute("""
        SELECT rsid, gene, chromosome, position, summary, revision_id FROM snps
    """):
        rsid, gene, chrom, pos, summary, rev = row
        rsid_lookup[rsid.lower()] = (rsid, gene, summary, rev)

    # Genotype lookup: (rsid_lower, alleles_upper) -> (magnitude, repute, summary)
    geno_lookup: dict[tuple[str, str], tuple] = {}
    for row in snpedia_conn.execute("""
        SELECT rsid, alleles, magnitude, repute, summary
        FROM genotypes WHERE rsid IS NOT NULL AND alleles IS NOT NULL
    """):
        rsid, alleles, mag, repute, summary = row
        key = (rsid.lower(), alleles.upper())
        geno_lookup[key] = (mag or 0.0, repute, summary)
        # Also store flipped
        m = re.match(r"\((.);(.)\)", alleles)
        if m:
            flipped = f"({m.group(2)};{m.group(1)})"
            flipped_key = (rsid.lower(), flipped.upper())
            if flipped_key not in geno_lookup:
                geno_lookup[flipped_key] = (mag or 0.0, repute, summary)

    tqdm.write(f"    {len(pos_lookup):,} chr:pos, {len(rsid_lookup):,} rsIDs, {len(geno_lookup):,} genotypes")

    # Fetch all extracted variants
    variant_count = results_conn.execute("SELECT COUNT(*) FROM vcf_variants").fetchone()[0]

    # Clear previous matches
    results_conn.execute("DELETE FROM matches")

    cursor = results_conn.execute("SELECT * FROM vcf_variants")
    now = datetime.now(timezone.utc).isoformat()
    batch = []
    stats = {"matched": 0, "by_pos": 0, "by_rsid": 0, "with_geno": 0}

    pbar = tqdm(
        total=variant_count,
        desc="Matching SNPs",
        unit=" variants",
        colour="#fbbf24",
        bar_format="{l_bar}{bar}| {n_fmt}/{total_fmt} [{elapsed}<{remaining}, {rate_fmt}] {postfix}",
    )

    for chrom, pos, ref, alt, genotype, qual, filt, vcf_rsids in cursor:
        pbar.update(1)

        # Try position match
        snp = pos_lookup.get((chrom, pos))
        match_method = "pos"

        # Try rsID match
        if not snp and vcf_rsids:
            for rid in vcf_rsids.split(","):
                rid = rid.strip().lower()
                if rid in rsid_lookup:
                    snp = rsid_lookup[rid]
                    match_method = "rsid"
                    break

        if not snp:
            continue

        rsid, gene, snp_summary, snp_rev = snp
        rsid_lower = rsid.lower()

        # Genotype lookup
        geno = geno_lookup.get((rsid_lower, genotype.upper()))
        if not geno:
            # Try unsorted order from VCF
            geno = geno_lookup.get((rsid_lower, genotype.upper()))

        geno_mag, geno_repute, geno_summary = geno if geno else (None, None, None)

        batch.append((
            chrom, pos, ref, alt, genotype, qual, filt,
            rsid, gene, snp_summary,
            genotype, geno_mag, geno_repute, geno_summary,
            match_method, snp_rev, now,
        ))

        stats["matched"] += 1
        if match_method == "pos":
            stats["by_pos"] += 1
        else:
            stats["by_rsid"] += 1
        if geno_mag is not None:
            stats["with_geno"] += 1

        if len(batch) >= 5000:
            results_conn.executemany(
                "INSERT OR REPLACE INTO matches VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
                batch,
            )
            batch.clear()
            pbar.set_postfix_str(
                f"hits={stats['matched']:,} pos={stats['by_pos']:,} rsid={stats['by_rsid']:,}",
                refresh=False,
            )

    if batch:
        results_conn.executemany(
            "INSERT OR REPLACE INTO matches VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
            batch,
        )

    pbar.close()

    # Save fingerprint
    fp = snpedia_fingerprint(snpedia_conn)
    results_conn.executemany("INSERT OR REPLACE INTO scan_meta (key, value) VALUES (?, ?)", [
        ("snpedia_fingerprint", fp),
        ("matched_at", now),
        ("match_total", str(stats["matched"])),
        ("match_by_pos", str(stats["by_pos"])),
        ("match_by_rsid", str(stats["by_rsid"])),
        ("match_with_geno", str(stats["with_geno"])),
    ])
    results_conn.commit()

    notable = results_conn.execute(
        "SELECT COUNT(*) FROM matches WHERE geno_magnitude >= 2"
    ).fetchone()[0]

    tqdm.write(f"\n  Results:")
    tqdm.write(f"    Matched:            {stats['matched']:,} SNPs")
    tqdm.write(f"      by chr:pos:       {stats['by_pos']:,}")
    tqdm.write(f"      by rsID:          {stats['by_rsid']:,}")
    tqdm.write(f"    With genotype info:  {stats['with_geno']:,}")
    tqdm.write(f"    Notable (mag >= 2):  {notable:,}")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Scan a VCF against SNPedia and store results in a SQLite database.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="Example:\n  uv run scan --vcf my_dna.vcf.gz --snpedia data/snpedia.db --out results.db",
    )
    parser.add_argument("--vcf", required=True, help="Path to VCF file (.vcf or .vcf.gz)")
    parser.add_argument("--snpedia", default="data/snpedia.db", help="SNPedia SQLite database (default: data/snpedia.db)")
    parser.add_argument("--out", default="results.db", help="Output results SQLite database (default: results.db)")
    args = parser.parse_args()

    if not Path(args.vcf).exists():
        print(f"Error: VCF file not found: {args.vcf}", file=sys.stderr)
        sys.exit(1)
    if not Path(args.snpedia).exists():
        print(f"Error: SNPedia database not found: {args.snpedia}", file=sys.stderr)
        sys.exit(1)

    vcf_size = os.path.getsize(args.vcf)
    snpedia_size = os.path.getsize(args.snpedia)

    print()
    print("  DNAscanner")
    print(f"  VCF:     {args.vcf} ({vcf_size / 1024 / 1024:.0f} MB)")
    print(f"  SNPedia: {args.snpedia} ({snpedia_size / 1024 / 1024:.0f} MB)")
    print(f"  Output:  {args.out}")
    print()

    results_conn = init_results_db(args.out)
    snpedia_conn = sqlite3.connect(args.snpedia)

    t0 = time.time()

    # Step 1: Extract VCF variants (skip if file unchanged)
    if vcf_needs_extraction(results_conn, args.vcf):
        print("  Step 1: Extracting variants from VCF...")
        extract_vcf_variants(args.vcf, results_conn)
    else:
        cached_count = results_conn.execute("SELECT COUNT(*) FROM vcf_variants").fetchone()[0]
        print(f"  Step 1: VCF unchanged, using cached extraction ({cached_count:,} SNPs)")

    # Step 2: Match against SNPedia (skip if SNPedia unchanged)
    if matching_needs_update(results_conn, snpedia_conn):
        print("\n  Step 2: Matching against SNPedia...")
        match_variants(results_conn, snpedia_conn)
    else:
        cached_matches = results_conn.execute("SELECT COUNT(*) FROM matches").fetchone()[0]
        print(f"  Step 2: SNPedia unchanged, using cached matches ({cached_matches:,} matches)")

    snpedia_conn.close()

    elapsed = time.time() - t0
    print(f"\n  Completed in {elapsed:.1f}s")

    # Print top notable findings
    notable = results_conn.execute("""
        SELECT geno_magnitude, geno_repute, rsid, gene, genotype, geno_summary, snp_summary
        FROM matches
        WHERE geno_magnitude >= 2
        ORDER BY geno_magnitude DESC
        LIMIT 20
    """).fetchall()

    if notable:
        print(f"\n  Top findings (magnitude >= 2):")
        print(f"  {'Mag':>5}  {'Repute':<6}  {'rsID':<12}  {'Gene':<8}  {'Geno':<7}  Summary")
        print(f"  {'---':>5}  {'------':<6}  {'----':<12}  {'----':<8}  {'----':<7}  -------")
        for mag, repute, rsid, gene, geno, geno_sum, snp_sum in notable:
            mag_s = f"{mag:.1f}" if mag else "-"
            repute_s = (repute or "-")[:6]
            gene_s = (gene or "-")[:8]
            summary = (geno_sum or snp_sum or "-")[:50]
            print(f"  {mag_s:>5}  {repute_s:<6}  {rsid:<12}  {gene_s:<8}  {geno:<7}  {summary}")

    results_conn.close()
    print()


if __name__ == "__main__":
    main()
