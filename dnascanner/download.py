#!/usr/bin/env python3
"""
SNPedia bulk downloader.

Downloads SNP, genotype, and genoset data from SNPedia's MediaWiki API
into a local SQLite database. Supports incremental updates — only fetches
pages that are new or have been revised since the last download.

Usage:
    uv run download [--db data/snpedia.db] [--delay 1.0] [--batch 50]
"""

import argparse
import random
import re
import sqlite3
import time
from datetime import datetime, timezone

import requests
from tqdm import tqdm

API_URL = "https://bots.snpedia.com/api.php"
USER_AGENT = "dnascanner/1.0 (local mirror; contact: github.com/viberaven/dnascanner)"
MAX_RETRIES = 6
MAX_BACKOFF = 30  # cap retry wait at 30s — 502s are transient, no point waiting longer

# Persistent session — reuses connections and cookies (Incapsula CDN likes this)
session = requests.Session()
session.headers["User-Agent"] = USER_AGENT


# ---------------------------------------------------------------------------
# Database
# ---------------------------------------------------------------------------

SCHEMA = """
CREATE TABLE IF NOT EXISTS snps (
    rsid            TEXT PRIMARY KEY,
    gene            TEXT,
    chromosome      TEXT,
    position        INTEGER,
    orientation     TEXT,
    stabilized_orientation TEXT,
    gmaf            REAL,
    summary         TEXT,
    raw_wikitext    TEXT,
    revision_id     INTEGER,
    revision_ts     TEXT,
    imported_at     TEXT
);

CREATE TABLE IF NOT EXISTS genotypes (
    page_title      TEXT PRIMARY KEY,
    rsid            TEXT,
    alleles         TEXT,
    magnitude       REAL,
    repute          TEXT,
    summary         TEXT,
    raw_wikitext    TEXT,
    revision_id     INTEGER,
    revision_ts     TEXT,
    imported_at     TEXT
);

CREATE TABLE IF NOT EXISTS genosets (
    page_title      TEXT PRIMARY KEY,
    magnitude       REAL,
    repute          TEXT,
    summary         TEXT,
    criteria        TEXT,
    raw_wikitext    TEXT,
    revision_id     INTEGER,
    revision_ts     TEXT,
    imported_at     TEXT
);

CREATE TABLE IF NOT EXISTS download_meta (
    key             TEXT PRIMARY KEY,
    value           TEXT
);

CREATE INDEX IF NOT EXISTS idx_genotypes_rsid ON genotypes(rsid);
CREATE INDEX IF NOT EXISTS idx_snps_gene ON snps(gene);
CREATE INDEX IF NOT EXISTS idx_snps_chromosome ON snps(chromosome);
"""


def init_db(db_path: str) -> sqlite3.Connection:
    conn = sqlite3.connect(db_path)
    conn.execute("PRAGMA journal_mode=WAL")
    conn.execute("PRAGMA synchronous=NORMAL")
    conn.executescript(SCHEMA)
    conn.commit()
    return conn


# ---------------------------------------------------------------------------
# MediaWiki API helpers
# ---------------------------------------------------------------------------

def api_request(params: dict, retries: int = MAX_RETRIES) -> dict:
    """Make a GET request to the SNPedia MediaWiki API with retry logic.

    Uses capped exponential backoff with jitter and a persistent session
    to stay friendly with the Incapsula CDN.
    """
    params["format"] = "json"

    for attempt in range(retries):
        try:
            resp = session.get(API_URL, params=params, timeout=30)
            resp.raise_for_status()
            return resp.json()
        except (requests.RequestException, ValueError) as e:
            # Capped backoff: 3s, 6s, 12s, 24s, 30s, 30s — plus jitter
            base_wait = min(3 * (2 ** attempt), MAX_BACKOFF)
            jitter = random.uniform(0.5, 1.5)
            wait = base_wait * jitter

            if attempt < retries - 1:
                tqdm.write(f"  [retry {attempt+1}/{retries}] {e} — waiting {wait:.0f}s")
                time.sleep(wait)
            else:
                raise RuntimeError(f"API request failed after {retries} attempts: {e}") from e


def fetch_category_members(category: str, delay: float) -> list[dict]:
    """Fetch all members of a category, handling continuation."""
    members = []
    params = {
        "action": "query",
        "list": "categorymembers",
        "cmtitle": f"Category:{category}",
        "cmlimit": "500",
        "cmtype": "page",
    }

    pbar = tqdm(desc=f"Listing {category}", unit=" pages", colour="#38bdf8",
                bar_format="{l_bar}{bar}| {n_fmt} pages [{elapsed}, {rate_fmt}]")

    while True:
        data = api_request(params)
        batch = data["query"]["categorymembers"]
        members.extend(batch)
        pbar.update(len(batch))

        if "continue" not in data:
            break
        params["cmcontinue"] = data["continue"]["cmcontinue"]
        # Category listing is heavy; use 2x delay
        time.sleep(delay * 2)

    pbar.close()
    return members


def fetch_page_revisions(titles: list[str]) -> dict:
    """Fetch latest revision ID and timestamp for a batch of pages.

    Returns {title: {"revid": int, "timestamp": str}} for each page found.
    """
    data = api_request({
        "action": "query",
        "titles": "|".join(titles),
        "prop": "revisions",
        "rvprop": "ids|timestamp",
    })
    result = {}
    for page in data["query"]["pages"].values():
        if "revisions" in page:
            rev = page["revisions"][0]
            result[page["title"]] = {
                "revid": rev["revid"],
                "timestamp": rev["timestamp"],
            }
    return result


def fetch_page_contents(titles: list[str]) -> dict:
    """Fetch full wikitext content + revision info for a batch of pages.

    Returns {title: {"revid": int, "timestamp": str, "content": str}}.
    """
    data = api_request({
        "action": "query",
        "titles": "|".join(titles),
        "prop": "revisions",
        "rvprop": "ids|timestamp|content",
        "rvslots": "main",
    })
    result = {}
    for page in data["query"]["pages"].values():
        if "revisions" in page:
            rev = page["revisions"][0]
            content = rev.get("slots", {}).get("main", {}).get("*", "")
            # Fallback for legacy format
            if not content:
                content = rev.get("*", "")
            result[page["title"]] = {
                "revid": rev["revid"],
                "timestamp": rev["timestamp"],
                "content": content,
            }
    return result


# ---------------------------------------------------------------------------
# Wiki template parser
# ---------------------------------------------------------------------------

def parse_template_params(wikitext: str, template_name: str) -> dict:
    """Extract key=value parameters from a named wiki template."""
    pattern = r"\{\{" + re.escape(template_name) + r"\s*\n(.*?)\}\}"
    match = re.search(pattern, wikitext, re.DOTALL | re.IGNORECASE)
    if not match:
        return {}

    params = {}
    for line in match.group(1).split("\n"):
        line = line.strip()
        if line.startswith("|"):
            line = line[1:]
        if "=" in line:
            key, _, value = line.partition("=")
            params[key.strip().lower()] = value.strip()
    return params


def parse_snp(title: str, wikitext: str) -> dict:
    """Parse a SNP page into structured fields."""
    p = parse_template_params(wikitext, "Rsnum")
    return {
        "rsid": title.lower(),
        "gene": p.get("gene") or p.get("gene_s") or None,
        "chromosome": p.get("chromosome"),
        "position": _int_or_none(p.get("position")),
        "orientation": p.get("orientation"),
        "stabilized_orientation": p.get("stabilizedorientation"),
        "gmaf": _float_or_none(p.get("gmaf")),
        "summary": p.get("summary"),
    }


def parse_genotype(title: str, wikitext: str) -> dict:
    """Parse a genotype page into structured fields."""
    # Genotype pages use {{Genotype}} template
    p = parse_template_params(wikitext, "Genotype")
    # Title format is like "Rs53576(A;G)"
    rsid = None
    alleles = None
    m = re.match(r"((?:Rs|rs|i)\d+)\(([^)]+)\)", title)
    if m:
        rsid = m.group(1).lower()
        alleles = f"({m.group(2)})"
    return {
        "rsid": rsid,
        "alleles": alleles,
        "magnitude": _float_or_none(p.get("magnitude")),
        "repute": p.get("repute") or None,
        "summary": p.get("summary"),
    }


def parse_genoset(title: str, wikitext: str) -> dict:
    """Parse a genoset page into structured fields."""
    p = parse_template_params(wikitext, "Genoset")
    return {
        "magnitude": _float_or_none(p.get("magnitude")),
        "repute": p.get("repute") or None,
        "summary": p.get("summary"),
        "criteria": p.get("criteria"),
    }


def _int_or_none(val):
    if val is None:
        return None
    try:
        return int(val)
    except (ValueError, TypeError):
        return None


def _float_or_none(val):
    if val is None:
        return None
    try:
        return float(val)
    except (ValueError, TypeError):
        return None


# ---------------------------------------------------------------------------
# Incremental download logic
# ---------------------------------------------------------------------------

def get_stored_revisions(conn: sqlite3.Connection, table: str, titles: list[str]) -> dict:
    """Get {title: revision_id} for pages already in the database."""
    if table == "snps":
        key_col = "rsid"
        # SNP titles are stored lowercased
        placeholders = ",".join("?" for _ in titles)
        rows = conn.execute(
            f"SELECT {key_col}, revision_id FROM {table} WHERE {key_col} IN ({placeholders})",
            [t.lower() for t in titles],
        ).fetchall()
        return {r[0]: r[1] for r in rows}
    else:
        placeholders = ",".join("?" for _ in titles)
        rows = conn.execute(
            f"SELECT page_title, revision_id FROM {table} WHERE page_title IN ({placeholders})",
            titles,
        ).fetchall()
        return {r[0]: r[1] for r in rows}


def upsert_snp(conn: sqlite3.Connection, title: str, content: str, revid: int, rev_ts: str):
    now = datetime.now(timezone.utc).isoformat()
    parsed = parse_snp(title, content)
    conn.execute("""
        INSERT INTO snps (rsid, gene, chromosome, position, orientation,
                          stabilized_orientation, gmaf, summary, raw_wikitext,
                          revision_id, revision_ts, imported_at)
        VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
        ON CONFLICT(rsid) DO UPDATE SET
            gene=excluded.gene, chromosome=excluded.chromosome,
            position=excluded.position, orientation=excluded.orientation,
            stabilized_orientation=excluded.stabilized_orientation,
            gmaf=excluded.gmaf, summary=excluded.summary,
            raw_wikitext=excluded.raw_wikitext,
            revision_id=excluded.revision_id, revision_ts=excluded.revision_ts,
            imported_at=excluded.imported_at
    """, (
        parsed["rsid"], parsed["gene"], parsed["chromosome"],
        parsed["position"], parsed["orientation"],
        parsed["stabilized_orientation"], parsed["gmaf"],
        parsed["summary"], content, revid, rev_ts, now,
    ))


def upsert_genotype(conn: sqlite3.Connection, title: str, content: str, revid: int, rev_ts: str):
    now = datetime.now(timezone.utc).isoformat()
    parsed = parse_genotype(title, content)
    conn.execute("""
        INSERT INTO genosets (page_title, magnitude, repute, summary,
                              criteria, raw_wikitext, revision_id, revision_ts, imported_at)
        VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)
        ON CONFLICT(page_title) DO UPDATE SET
            magnitude=excluded.magnitude, repute=excluded.repute,
            summary=excluded.summary, criteria=excluded.criteria,
            raw_wikitext=excluded.raw_wikitext,
            revision_id=excluded.revision_id, revision_ts=excluded.revision_ts,
            imported_at=excluded.imported_at
    """, (title, parsed["magnitude"], parsed["repute"], parsed["summary"],
          parsed.get("criteria"), content, revid, rev_ts, now))


def upsert_genotype_row(conn: sqlite3.Connection, title: str, content: str, revid: int, rev_ts: str):
    now = datetime.now(timezone.utc).isoformat()
    parsed = parse_genotype(title, content)
    conn.execute("""
        INSERT INTO genotypes (page_title, rsid, alleles, magnitude, repute,
                               summary, raw_wikitext, revision_id, revision_ts, imported_at)
        VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
        ON CONFLICT(page_title) DO UPDATE SET
            rsid=excluded.rsid, alleles=excluded.alleles,
            magnitude=excluded.magnitude, repute=excluded.repute,
            summary=excluded.summary, raw_wikitext=excluded.raw_wikitext,
            revision_id=excluded.revision_id, revision_ts=excluded.revision_ts,
            imported_at=excluded.imported_at
    """, (title, parsed["rsid"], parsed["alleles"], parsed["magnitude"],
          parsed["repute"], parsed["summary"], content, revid, rev_ts, now))


def upsert_genoset(conn: sqlite3.Connection, title: str, content: str, revid: int, rev_ts: str):
    now = datetime.now(timezone.utc).isoformat()
    parsed = parse_genoset(title, content)
    conn.execute("""
        INSERT INTO genosets (page_title, magnitude, repute, summary,
                              criteria, raw_wikitext, revision_id, revision_ts, imported_at)
        VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)
        ON CONFLICT(page_title) DO UPDATE SET
            magnitude=excluded.magnitude, repute=excluded.repute,
            summary=excluded.summary, criteria=excluded.criteria,
            raw_wikitext=excluded.raw_wikitext,
            revision_id=excluded.revision_id, revision_ts=excluded.revision_ts,
            imported_at=excluded.imported_at
    """, (title, parsed["magnitude"], parsed["repute"], parsed["summary"],
          parsed.get("criteria"), content, revid, rev_ts, now))


# ---------------------------------------------------------------------------
# Main download orchestration
# ---------------------------------------------------------------------------

CATEGORIES = [
    ("Is_a_snp", "snps"),
    ("Is_a_genotype", "genotypes"),
    ("Is_a_genoset", "genosets"),
]


def download_category(
    conn: sqlite3.Connection,
    category: str,
    table: str,
    batch_size: int,
    delay: float,
    refresh: bool = False,
):
    """Download all pages in a category, skipping unchanged ones."""
    print(f"\n{'='*60}")
    print(f"Category: {category} -> table: {table}")
    if refresh:
        print(f"  Mode: refresh (checking remote revisions)")
    print(f"{'='*60}")

    # Step 1: Get all page titles in this category
    members = fetch_category_members(category, delay)
    total = len(members)

    # Step 2: Process in batches with adaptive delay
    stats = {"new": 0, "updated": 0, "skipped": 0, "errors": 0}
    titles = [m["title"] for m in members]

    # Adaptive delay: starts at base, increases on 502s, decreases on success streaks
    current_delay = delay
    success_streak = 0

    pbar = tqdm(
        total=total,
        desc=f"Downloading {table}",
        unit=" pages",
        colour="#4ade80",
        bar_format="{l_bar}{bar}| {n_fmt}/{total_fmt} [{elapsed}<{remaining}, {rate_fmt}] {postfix}",
    )

    for i in range(0, total, batch_size):
        batch_titles = titles[i : i + batch_size]

        # Check which pages we already have
        stored = get_stored_revisions(conn, table, batch_titles)

        if not refresh:
            # Skip entire batch if all pages exist locally
            needs_download = []
            for title in batch_titles:
                lookup_key = title.lower() if table == "snps" else title
                if lookup_key in stored:
                    stats["skipped"] += 1
                else:
                    needs_download.append(title)

            if not needs_download:
                pbar.update(len(batch_titles))
                pbar.set_postfix_str(
                    f"new={stats['new']:,} upd={stats['updated']:,} "
                    f"skip={stats['skipped']:,} err={stats['errors']:,}  "
                    f"delay={current_delay:.1f}s",
                    refresh=False,
                )
                continue
        else:
            # Fetch remote revisions and compare
            time.sleep(current_delay)
            try:
                remote_revs = fetch_page_revisions(batch_titles)
                success_streak += 1
            except RuntimeError as e:
                tqdm.write(f"  ERROR fetching revisions: {e}")
                stats["errors"] += len(batch_titles)
                pbar.update(len(batch_titles))
                current_delay = min(current_delay * 2, delay * 10)
                success_streak = 0
                tqdm.write(f"  Delay increased to {current_delay:.1f}s")
                continue

            needs_download = []
            for title in batch_titles:
                if title not in remote_revs:
                    continue
                remote_rev = remote_revs[title]["revid"]
                lookup_key = title.lower() if table == "snps" else title
                stored_rev = stored.get(lookup_key)
                if stored_rev is None:
                    needs_download.append(title)
                elif stored_rev != remote_rev:
                    needs_download.append(title)
                else:
                    stats["skipped"] += 1

            if not needs_download:
                pbar.update(len(batch_titles))
                pbar.set_postfix_str(
                    f"new={stats['new']:,} upd={stats['updated']:,} "
                    f"skip={stats['skipped']:,} err={stats['errors']:,}  "
                    f"delay={current_delay:.1f}s",
                    refresh=False,
                )
                continue

        # Download full content for pages that need it
        time.sleep(current_delay)
        try:
            contents = fetch_page_contents(needs_download)
            success_streak += 1
        except RuntimeError as e:
            tqdm.write(f"  ERROR fetching content: {e}")
            stats["errors"] += len(needs_download)
            pbar.update(len(batch_titles))
            current_delay = min(current_delay * 2, delay * 10)
            success_streak = 0
            tqdm.write(f"  Delay increased to {current_delay:.1f}s")
            continue

        # Ease back toward base delay after 10 consecutive successes
        if success_streak >= 10 and current_delay > delay:
            current_delay = max(current_delay * 0.8, delay)
            success_streak = 0

        # Upsert into the database
        for title, page_data in contents.items():
            lookup_key = title.lower() if table == "snps" else title
            was_stored = lookup_key in stored

            try:
                if table == "snps":
                    upsert_snp(conn, title, page_data["content"],
                               page_data["revid"], page_data["timestamp"])
                elif table == "genotypes":
                    upsert_genotype_row(conn, title, page_data["content"],
                                        page_data["revid"], page_data["timestamp"])
                elif table == "genosets":
                    upsert_genoset(conn, title, page_data["content"],
                                   page_data["revid"], page_data["timestamp"])

                if was_stored:
                    stats["updated"] += 1
                else:
                    stats["new"] += 1
            except Exception as e:
                tqdm.write(f"  Error processing {title}: {e}")
                stats["errors"] += 1

        conn.commit()
        pbar.update(len(batch_titles))
        pbar.set_postfix_str(
            f"new={stats['new']:,} upd={stats['updated']:,} "
            f"skip={stats['skipped']:,} err={stats['errors']:,}  "
            f"delay={current_delay:.1f}s",
            refresh=False,
        )

    pbar.close()
    conn.commit()

    tqdm.write(f"  Done: {stats['new']:,} new, {stats['updated']:,} updated, "
               f"{stats['skipped']:,} skipped, {stats['errors']:,} errors")
    return stats


# ---------------------------------------------------------------------------
# Pre-built database download
# ---------------------------------------------------------------------------

PREBUILT_URL = "https://data.viberaven.com/dnascanner/snpedia.db"


def download_prebuilt(db_path: str) -> bool:
    """Try to download a pre-built SNPedia database. Returns True on success."""
    import shutil
    from pathlib import Path

    print(f"  Fetching pre-built database from {PREBUILT_URL}...")

    try:
        resp = session.get(PREBUILT_URL, stream=True, timeout=30)
        resp.raise_for_status()
    except requests.RequestException as e:
        print(f"  Pre-built download failed: {e}")
        return False

    total = int(resp.headers.get("content-length", 0))
    tmp_path = db_path + ".download"

    try:
        with open(tmp_path, "wb") as f:
            with tqdm(
                total=total or None,
                unit="B",
                unit_scale=True,
                unit_divisor=1024,
                desc="Downloading",
                colour="#38bdf8",
            ) as pbar:
                for chunk in resp.iter_content(chunk_size=256 * 1024):
                    f.write(chunk)
                    pbar.update(len(chunk))

        # Verify it's a valid SQLite database
        conn = sqlite3.connect(tmp_path)
        count = conn.execute("SELECT COUNT(*) FROM snps").fetchone()[0]
        conn.close()

        if count == 0:
            print(f"  Pre-built database is empty, ignoring")
            Path(tmp_path).unlink(missing_ok=True)
            return False

        # Replace the target file
        shutil.move(tmp_path, db_path)
        print(f"  Downloaded pre-built database ({count:,} SNPs)")
        return True

    except Exception as e:
        print(f"  Pre-built download failed: {e}")
        Path(tmp_path).unlink(missing_ok=True)
        return False


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description="Download SNPedia into a local SQLite database.")
    parser.add_argument("--db", default="data/snpedia.db", help="SQLite database path (default: data/snpedia.db)")
    parser.add_argument("--delay", type=float, default=1.0, help="Delay between API requests in seconds (default: 1.0)")
    parser.add_argument("--batch", type=int, default=50, help="Pages per API batch (default: 50, max 50)")
    parser.add_argument("--category", choices=["snps", "genotypes", "genosets", "all"], default="all",
                        help="Which category to download (default: all)")
    parser.add_argument("--refresh", action="store_true",
                        help="Check remote revisions and update existing pages if newer (slower, makes extra API calls)")
    parser.add_argument("--no-prebuilt", action="store_true",
                        help="Skip pre-built database download, mirror directly from SNPedia API")
    args = parser.parse_args()

    batch_size = min(args.batch, 50)  # MediaWiki API limit for multi-title queries

    # Ensure parent directory exists
    from pathlib import Path
    Path(args.db).parent.mkdir(parents=True, exist_ok=True)

    print(f"SNPedia Downloader")
    print(f"  Database: {args.db}")

    # Try pre-built download first
    db_exists = Path(args.db).exists()
    need_mirror = True

    if not args.no_prebuilt:
        if not db_exists:
            # No local DB — try pre-built first
            if download_prebuilt(args.db):
                need_mirror = False
        elif args.refresh:
            # Refresh mode — fetch latest pre-built, then refresh from API
            download_prebuilt(args.db)
            # Still need to mirror for refresh
            need_mirror = True

    if not need_mirror and not args.refresh:
        # Pre-built download succeeded, no refresh needed
        conn = sqlite3.connect(args.db)
        for table in ["snps", "genotypes", "genosets"]:
            count = conn.execute(f"SELECT COUNT(*) FROM {table}").fetchone()[0]
            print(f"  {table}: {count:,} rows")
        conn.close()
        return

    # Mirror from SNPedia API
    print(f"  Delay: {args.delay}s between requests")
    print(f"  Batch size: {batch_size}")
    print(f"  API: {API_URL}")

    conn = init_db(args.db)

    start_time = time.time()
    total_stats = {"new": 0, "updated": 0, "skipped": 0, "errors": 0}

    categories_to_run = CATEGORIES if args.category == "all" else [
        c for c in CATEGORIES if c[1] == args.category
    ]

    for category, table in categories_to_run:
        stats = download_category(conn, category, table, batch_size, args.delay, args.refresh)
        for k in total_stats:
            total_stats[k] += stats[k]

    # Record completion time
    now = datetime.now(timezone.utc).isoformat()
    conn.execute(
        "INSERT OR REPLACE INTO download_meta (key, value) VALUES (?, ?)",
        ("last_download", now),
    )
    conn.commit()

    elapsed = time.time() - start_time
    print(f"\n{'='*60}")
    print(f"Completed in {elapsed/60:.1f} minutes")
    print(f"  Total: {total_stats['new']:,} new, {total_stats['updated']:,} updated, "
          f"{total_stats['skipped']:,} skipped, {total_stats['errors']:,} errors")

    # Print DB stats
    for table in ["snps", "genotypes", "genosets"]:
        count = conn.execute(f"SELECT COUNT(*) FROM {table}").fetchone()[0]
        print(f"  {table}: {count:,} rows")

    conn.close()


if __name__ == "__main__":
    main()
