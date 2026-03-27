#!/usr/bin/env python3
"""
Report generator — reads scan results from SQLite and produces an
interactive HTML report and/or a Markdown report.

Usage:
    uv run report
    uv run report --format md
    uv run report --format all
"""

import argparse
import html
import os
import sqlite3
import sys
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path

from tqdm import tqdm


# ---------------------------------------------------------------------------
# Shared data loading
# ---------------------------------------------------------------------------

@dataclass
class ReportData:
    vcf_name: str
    total_variants: str
    total_matches: int
    with_genotype: int
    notable: int
    bad: int
    good: int
    rows: list[tuple]  # (chrom, pos, rsid, gene, geno, mag, repute, geno_sum, snp_sum)


def get_meta(conn: sqlite3.Connection, key: str, default: str = "") -> str:
    row = conn.execute("SELECT value FROM scan_meta WHERE key = ?", (key,)).fetchone()
    return row[0] if row else default


def load_report_data(db_path: str) -> ReportData:
    conn = sqlite3.connect(db_path)

    total_matches = conn.execute("SELECT COUNT(*) FROM matches").fetchone()[0]
    with_genotype = conn.execute("SELECT COUNT(*) FROM matches WHERE geno_magnitude IS NOT NULL").fetchone()[0]
    notable = conn.execute("SELECT COUNT(*) FROM matches WHERE geno_magnitude >= 2").fetchone()[0]
    bad = conn.execute("SELECT COUNT(*) FROM matches WHERE LOWER(geno_repute) = 'bad'").fetchone()[0]
    good = conn.execute("SELECT COUNT(*) FROM matches WHERE LOWER(geno_repute) = 'good'").fetchone()[0]

    vcf_path = get_meta(conn, "vcf_path", "unknown")
    total_variants = get_meta(conn, "vcf_total_variants", "?")

    rows = conn.execute("""
        SELECT chromosome, position, rsid, gene, genotype,
               geno_magnitude, geno_repute, geno_summary, snp_summary
        FROM matches
        ORDER BY COALESCE(geno_magnitude, -1) DESC
    """).fetchall()

    conn.close()

    return ReportData(
        vcf_name=Path(vcf_path).name,
        total_variants=f"{int(total_variants):,}" if total_variants.isdigit() else total_variants,
        total_matches=total_matches,
        with_genotype=with_genotype,
        notable=notable,
        bad=bad,
        good=good,
        rows=rows,
    )


# ---------------------------------------------------------------------------
# Markdown report
# ---------------------------------------------------------------------------

def generate_markdown(data: ReportData, out_path: str):
    lines = []
    lines.append("# DNAscanner Report")
    lines.append("")
    lines.append(f"Generated {datetime.now().strftime('%Y-%m-%d %H:%M')} | VCF: {data.vcf_name} | {data.total_variants} variants scanned")
    lines.append("")
    lines.append("## Summary")
    lines.append("")
    lines.append(f"- **SNPedia matches:** {data.total_matches:,}")
    lines.append(f"- **With genotype info:** {data.with_genotype:,}")
    lines.append(f"- **Notable (mag >= 2):** {data.notable:,}")
    lines.append(f"- **Bad repute:** {data.bad:,}")
    lines.append(f"- **Good repute:** {data.good:,}")
    lines.append("")

    notable_rows = [r for r in data.rows if (r[5] or 0) >= 2]
    mild_rows = [r for r in data.rows if 1 <= (r[5] or 0) < 2]
    other_rows = [r for r in data.rows if (r[5] or 0) < 1]

    def md_table_row(row, include_location=False):
        chrom, pos, rsid, gene, geno, mag, repute, geno_sum, snp_sum = row
        mag_str = f"{mag:.1f}" if mag is not None else "-"
        repute_str = (repute or "-").strip()
        gene_str = gene or "-"
        summary = geno_sum or snp_sum or "-"
        snpedia_url = f"https://www.snpedia.com/index.php/{rsid}"
        base = f"| {mag_str} | {repute_str} | [{rsid}]({snpedia_url}) | {gene_str} | `{geno}` |"
        if include_location:
            return f"{base} {chrom} | {pos:,} | {summary} |"
        return f"{base} {summary} |"

    if notable_rows:
        lines.append("## Notable findings (magnitude >= 2)")
        lines.append("")
        lines.append("| Mag | Repute | rsID | Gene | Genotype | Summary |")
        lines.append("|----:|--------|------|------|----------|---------|")
        for row in tqdm(notable_rows, desc="Building markdown (notable)", unit=" rows", colour="#4ade80"):
            lines.append(md_table_row(row))
        lines.append("")

    if mild_rows:
        lines.append("## Mildly interesting (magnitude 1-2)")
        lines.append("")
        lines.append("| Mag | Repute | rsID | Gene | Genotype | Summary |")
        lines.append("|----:|--------|------|------|----------|---------|")
        for row in tqdm(mild_rows, desc="Building markdown (mild)", unit=" rows", colour="#4ade80"):
            lines.append(md_table_row(row))
        lines.append("")

    if other_rows:
        lines.append("## Other matches")
        lines.append("")
        lines.append("| Mag | Repute | rsID | Gene | Genotype | Chr | Position | Summary |")
        lines.append("|----:|--------|------|------|----------|-----|----------|---------|")
        for row in tqdm(other_rows, desc="Building markdown (other)", unit=" rows", colour="#4ade80"):
            lines.append(md_table_row(row, include_location=True))

    lines.append("")
    lines.append("---")
    lines.append("")
    lines.append("**Disclaimer:** This report is for informational and educational purposes only. "
                 "It is NOT medical advice. SNPedia data is community-curated and may contain errors. "
                 "Genetic variants interact in complex ways — a single SNP rarely determines an outcome. "
                 "Always consult a qualified healthcare professional or genetic counselor for medical decisions.")
    lines.append("")

    with open(out_path, "w") as f:
        f.write("\n".join(lines))

    out_size = os.path.getsize(out_path)
    print(f"\n  Markdown: {out_path} ({out_size / 1024:.0f} KB)")


# ---------------------------------------------------------------------------
# HTML report
# ---------------------------------------------------------------------------

REPORT_TEMPLATE = """\
<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>DNAscanner Report</title>
<style>
  :root {{
    --bg: #0f172a;
    --surface: #1e293b;
    --border: #334155;
    --text: #e2e8f0;
    --text-muted: #94a3b8;
    --accent: #38bdf8;
    --good: #4ade80;
    --bad: #f87171;
    --neutral: #fbbf24;
  }}
  * {{ margin: 0; padding: 0; box-sizing: border-box; }}
  body {{
    font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', system-ui, sans-serif;
    background: var(--bg);
    color: var(--text);
    line-height: 1.6;
    padding: 2rem;
  }}
  .container {{ max-width: 1200px; margin: 0 auto; }}
  h1 {{
    font-size: 1.8rem;
    margin-bottom: 0.5rem;
    color: var(--accent);
  }}
  .meta {{
    color: var(--text-muted);
    font-size: 0.85rem;
    margin-bottom: 2rem;
    padding-bottom: 1rem;
    border-bottom: 1px solid var(--border);
  }}
  .stats {{
    display: grid;
    grid-template-columns: repeat(auto-fit, minmax(180px, 1fr));
    gap: 1rem;
    margin-bottom: 2rem;
  }}
  .stat-card {{
    background: var(--surface);
    border: 1px solid var(--border);
    border-radius: 8px;
    padding: 1rem 1.2rem;
  }}
  .stat-card .value {{
    font-size: 1.8rem;
    font-weight: 700;
    color: var(--accent);
  }}
  .stat-card .label {{
    color: var(--text-muted);
    font-size: 0.8rem;
    text-transform: uppercase;
    letter-spacing: 0.05em;
  }}
  .filter-group {{
    display: flex;
    gap: 0.5rem;
    align-items: center;
    flex-wrap: wrap;
  }}
  .filter-row {{
    display: flex;
    gap: 1.5rem;
    margin-bottom: 0.75rem;
    align-items: center;
    flex-wrap: wrap;
  }}
  .filter-label {{
    color: var(--text-muted);
    font-size: 0.75rem;
    text-transform: uppercase;
    letter-spacing: 0.05em;
    min-width: 5rem;
  }}
  .filter-group button {{
    background: var(--surface);
    color: var(--text-muted);
    border: 1px solid var(--border);
    border-radius: 6px;
    padding: 0.4rem 0.8rem;
    cursor: pointer;
    font-size: 0.85rem;
    transition: all 0.15s;
  }}
  .filter-group button:hover, .filter-group button.active {{
    background: var(--accent);
    color: var(--bg);
    border-color: var(--accent);
  }}
  input[type="search"] {{
    background: var(--surface);
    color: var(--text);
    border: 1px solid var(--border);
    border-radius: 6px;
    padding: 0.5rem 0.8rem;
    font-size: 0.9rem;
    width: 100%;
    max-width: 400px;
    margin-bottom: 1rem;
  }}
  input[type="search"]::placeholder {{ color: var(--text-muted); }}
  table {{
    width: 100%;
    border-collapse: collapse;
    font-size: 0.85rem;
  }}
  th {{
    background: var(--surface);
    color: var(--text-muted);
    text-transform: uppercase;
    font-size: 0.75rem;
    letter-spacing: 0.05em;
    padding: 0.7rem 0.8rem;
    text-align: left;
    position: sticky;
    top: 0;
    cursor: pointer;
    user-select: none;
    border-bottom: 2px solid var(--border);
  }}
  th:hover {{ color: var(--accent); }}
  td {{
    padding: 0.6rem 0.8rem;
    border-bottom: 1px solid var(--border);
    vertical-align: top;
  }}
  tr:hover {{ background: rgba(56, 189, 248, 0.05); }}
  .mag {{
    display: inline-block;
    min-width: 2rem;
    text-align: center;
    padding: 0.15rem 0.4rem;
    border-radius: 4px;
    font-weight: 600;
    font-size: 0.8rem;
  }}
  .mag-high {{ background: rgba(248, 113, 113, 0.2); color: var(--bad); }}
  .mag-med {{ background: rgba(251, 191, 36, 0.2); color: var(--neutral); }}
  .mag-low {{ background: rgba(148, 163, 184, 0.15); color: var(--text-muted); }}
  .repute-good {{ color: var(--good); }}
  .repute-bad {{ color: var(--bad); }}
  .genotype {{
    font-family: 'SF Mono', 'Fira Code', 'Consolas', monospace;
    font-weight: 600;
    color: var(--accent);
  }}
  .rsid a {{
    color: var(--accent);
    text-decoration: none;
  }}
  .rsid a:hover {{ text-decoration: underline; }}
  .gene {{ color: var(--neutral); font-weight: 500; }}
  .no-data {{ color: var(--text-muted); font-style: italic; }}
  .disclaimer {{
    margin-top: 3rem;
    padding: 1.2rem;
    background: rgba(248, 113, 113, 0.08);
    border: 1px solid rgba(248, 113, 113, 0.25);
    border-radius: 8px;
    color: var(--text-muted);
    font-size: 0.8rem;
  }}
  .disclaimer strong {{ color: var(--bad); }}
</style>
</head>
<body>
<div class="container">
  <h1>DNAscanner Report</h1>
  <div class="meta">
    Generated {generated_at} | VCF: {vcf_name} | {total_variants} variants scanned
  </div>

  <div class="stats">
    <div class="stat-card">
      <div class="value">{total_matches}</div>
      <div class="label">SNPedia matches</div>
    </div>
    <div class="stat-card">
      <div class="value">{with_genotype}</div>
      <div class="label">With genotype info</div>
    </div>
    <div class="stat-card">
      <div class="value">{notable_count}</div>
      <div class="label">Notable (mag &ge; 2)</div>
    </div>
    <div class="stat-card">
      <div class="value">{bad_count}</div>
      <div class="label">Bad repute</div>
    </div>
    <div class="stat-card">
      <div class="value">{good_count}</div>
      <div class="label">Good repute</div>
    </div>
  </div>

  <div class="filter-row">
    <span class="filter-label">Repute</span>
    <div class="filter-group" id="repute-filters">
      <button class="active" onclick="setFilter('repute','all',this)">All</button>
      <button onclick="setFilter('repute','bad',this)">Bad</button>
      <button onclick="setFilter('repute','good',this)">Good</button>
    </div>
  </div>
  <div class="filter-row">
    <span class="filter-label">Magnitude</span>
    <div class="filter-group" id="mag-filters">
      <button onclick="setFilter('mag','all',this)">All</button>
      <button class="active" onclick="setFilter('mag','interesting',this)">Notable + Mildly interesting (mag&ge;1)</button>
      <button onclick="setFilter('mag','notable',this)">Notable (mag&ge;2)</button>
    </div>
  </div>

  <input type="search" id="search" placeholder="Search by rsID, gene, or summary..." oninput="doSearch(this.value)">

  <table id="results">
    <thead>
      <tr>
        <th onclick="sortTable(0)">Mag</th>
        <th onclick="sortTable(1)">Repute</th>
        <th onclick="sortTable(2)">rsID</th>
        <th onclick="sortTable(3)">Gene</th>
        <th onclick="sortTable(4)">Your Genotype</th>
        <th onclick="sortTable(5)">Summary</th>
        <th onclick="sortTable(6)">Chr</th>
        <th onclick="sortTable(7)">Position</th>
      </tr>
    </thead>
    <tbody>
{table_rows}
    </tbody>
  </table>

  <div class="disclaimer">
    <strong>Disclaimer:</strong> This report is for informational and educational purposes only.
    It is NOT medical advice. SNPedia data is community-curated and may contain errors.
    Genetic variants interact in complex ways and a single SNP rarely determines an outcome.
    Always consult a qualified healthcare professional or genetic counselor for medical decisions.
  </div>
</div>

<script>
let activeRepute = 'all';
let activeMag = 'interesting';

function setFilter(group, value, btn) {{
  btn.closest('.filter-group').querySelectorAll('button').forEach(b => b.classList.remove('active'));
  btn.classList.add('active');
  if (group === 'repute') activeRepute = value;
  else if (group === 'mag') activeMag = value;
  applyFilters();
}}

function applyFilters() {{
  document.querySelectorAll('#results tbody tr').forEach(tr => {{
    const repute = tr.dataset.repute || '';
    const mag = parseFloat(tr.dataset.mag) || 0;
    let show = true;
    if (activeRepute === 'bad') show = repute === 'bad';
    else if (activeRepute === 'good') show = repute === 'good';
    if (show && activeMag === 'interesting') show = mag >= 1;
    else if (show && activeMag === 'notable') show = mag >= 2;
    tr.style.display = show ? '' : 'none';
  }});
}}

function doSearch(q) {{
  q = q.toLowerCase();
  document.querySelectorAll('#results tbody tr').forEach(tr => {{
    tr.style.display = tr.textContent.toLowerCase().includes(q) ? '' : 'none';
  }});
}}

let sortDir = {{}};
function sortTable(col) {{
  const tbody = document.querySelector('#results tbody');
  const rows = Array.from(tbody.querySelectorAll('tr'));
  sortDir[col] = !(sortDir[col] || false);
  const dir = sortDir[col] ? 1 : -1;
  rows.sort((a, b) => {{
    let va = a.children[col].dataset.sort || a.children[col].textContent;
    let vb = b.children[col].dataset.sort || b.children[col].textContent;
    const na = parseFloat(va), nb = parseFloat(vb);
    if (!isNaN(na) && !isNaN(nb)) return (na - nb) * dir;
    return va.localeCompare(vb) * dir;
  }});
  rows.forEach(r => tbody.appendChild(r));
}}
// Sort by magnitude descending on load, then apply default filters
sortDir[0] = true;
sortTable(0);
applyFilters();
</script>
</body>
</html>
"""


def esc(text: str | None) -> str:
    if text is None:
        return ""
    return html.escape(str(text))


def mag_class(mag: float | None) -> str:
    if mag is None:
        return "mag-low"
    if mag >= 3:
        return "mag-high"
    if mag >= 1.5:
        return "mag-med"
    return "mag-low"


def repute_class(repute: str | None) -> str:
    if not repute:
        return ""
    r = repute.strip().lower()
    if r == "bad":
        return "repute-bad"
    if r == "good":
        return "repute-good"
    return ""


def generate_html(data: ReportData, out_path: str):
    rows_html = []
    for chrom, pos, rsid, gene, geno, mag, repute, geno_sum, snp_sum in tqdm(
        data.rows, desc="Building HTML", unit=" rows", colour="#4ade80"
    ):
        mag_str = f"{mag:.1f}" if mag is not None else ""
        repute_lower = (repute or "").strip().lower()
        summary = geno_sum or snp_sum or ""
        snpedia_url = f"https://www.snpedia.com/index.php/{esc(rsid)}"

        row_html = (
            f'      <tr data-repute="{esc(repute_lower)}" data-mag="{mag or 0}">'
            f'<td class="mag {mag_class(mag)}" data-sort="{mag if mag is not None else -1}">'
            f'{esc(mag_str) or "<span class=no-data>-</span>"}</td>'
            f'<td class="{repute_class(repute)}">'
            f'{esc(repute_lower) or "<span class=no-data>-</span>"}</td>'
            f'<td class="rsid"><a href="{snpedia_url}" target="_blank">{esc(rsid)}</a></td>'
            f'<td class="gene">{esc(gene) or "<span class=no-data>-</span>"}</td>'
            f'<td class="genotype">{esc(geno)}</td>'
            f'<td>{esc(summary) or "<span class=no-data>-</span>"}</td>'
            f'<td>{esc(chrom)}</td>'
            f'<td data-sort="{pos}">{pos:,}</td>'
            f'</tr>'
        )
        rows_html.append(row_html)

    report_html = REPORT_TEMPLATE.format(
        generated_at=datetime.now().strftime("%Y-%m-%d %H:%M"),
        vcf_name=esc(data.vcf_name),
        total_variants=data.total_variants,
        total_matches=f"{data.total_matches:,}",
        with_genotype=f"{data.with_genotype:,}",
        notable_count=f"{data.notable:,}",
        bad_count=f"{data.bad:,}",
        good_count=f"{data.good:,}",
        table_rows="\n".join(rows_html),
    )

    with open(out_path, "w") as f:
        f.write(report_html)

    out_size = os.path.getsize(out_path)
    print(f"\n  HTML: {out_path} ({out_size / 1024:.0f} KB)")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Generate a report from SNP scan results.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="Examples:\n"
               "  uv run report --db results.db\n"
               "  uv run report --db results.db --format md --out report.md\n"
               "  uv run report --db results.db --format all --out report",
    )
    parser.add_argument("--db", default="results.db", help="Scan results SQLite database (default: results.db)")
    parser.add_argument("--out", default="report", help="Output path (default: report). Extension added automatically unless --format is html/md")
    parser.add_argument("--format", choices=["html", "md", "all"], default="html",
                        help="Output format: html, md, or all (default: html)")
    args = parser.parse_args()

    if not Path(args.db).exists():
        print(f"Error: Results database not found: {args.db}", file=sys.stderr)
        print("Run `uv run scan` first to generate it.", file=sys.stderr)
        sys.exit(1)

    print()
    print("  DNAscanner Report Generator")
    print(f"  Input:  {args.db}")
    print(f"  Format: {args.format}")
    print()

    data = load_report_data(args.db)

    if args.format == "html":
        out = args.out if args.out.endswith(".html") else f"{args.out}.html"
        generate_html(data, out)
    elif args.format == "md":
        out = args.out if args.out.endswith(".md") else f"{args.out}.md"
        generate_markdown(data, out)
    elif args.format == "all":
        base = args.out.removesuffix(".html").removesuffix(".md")
        generate_html(data, f"{base}.html")
        generate_markdown(data, f"{base}.md")

    print(f"\n  {data.total_matches:,} matches | {data.notable:,} notable | {data.bad:,} bad | {data.good:,} good")
    print()


if __name__ == "__main__":
    main()
