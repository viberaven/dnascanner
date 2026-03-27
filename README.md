# DNAscanner

Scan your whole-genome sequencing VCF against a local mirror of [SNPedia](https://www.snpedia.com/) and generate an interactive HTML report of notable genetic findings.

Everything runs locally — no data leaves your machine.

## How it works

The pipeline has three stages, each producing a SQLite database that feeds the next:

```
uv run download          uv run scan                  uv run report

 SNPedia API        your VCF + data/snpedia.db         results.db
     |                       |                             |
     v                       v                             v
 data/snpedia.db  ───>  results.db  ──────────────>  report.html
```

1. **`uv run download`** pulls SNPedia's wiki pages via the MediaWiki API into `data/snpedia.db`. Parses wiki templates into structured columns (gene, chromosome, position, magnitude, repute, summary) while keeping raw wikitext. Incremental — re-runs only fetch pages with new revisions.

2. **`uv run scan`** extracts SNP variants from your VCF into `results.db`, then matches them against `data/snpedia.db` by chromosome:position (primary) and rsID from VEP/CSQ annotations (fallback). Incremental — skips VCF extraction if the file hasn't changed, skips matching if SNPedia data hasn't changed.

3. **`uv run report`** reads `results.db` and generates a self-contained interactive HTML report. Can be re-run anytime to regenerate the report without re-scanning.

## Requirements

- Python 3.10+
- [uv](https://docs.astral.sh/uv/) (recommended) or pip

## Setup

```bash
git clone <repo-url> && cd dnascanner
uv sync
```

## Usage

### Step 1: Download SNPedia

```bash
# Download everything (~112K SNPs, ~105K genotypes, ~283 genosets)
# Takes 3-5 hours due to API rate limiting
uv run download

# Or download just one category
uv run download --category snps

# Adjust delay if hitting too many 502s
uv run download --delay 3

# Resume an interrupted download — skips pages already in DB
uv run download

# Check for updated pages on SNPedia and re-download them
uv run download --refresh
```

| Flag | Default | Description |
|------|---------|-------------|
| `--db` | `data/snpedia.db` | SQLite database path |
| `--delay` | `1.0` | Seconds between API requests |
| `--batch` | `50` | Pages per API batch (max 50) |
| `--category` | `all` | `snps`, `genotypes`, `genosets`, or `all` |
| `--refresh` | off | Check remote revisions and update stale pages |

### Step 2: Scan your VCF

```bash
uv run scan --vcf /path/to/your.vcf.gz
```

| Flag | Default | Description |
|------|---------|-------------|
| `--vcf` | *(required)* | Path to VCF file (`.vcf` or `.vcf.gz`) |
| `--snpedia` | `data/snpedia.db` | SNPedia SQLite database |
| `--out` | `results.db` | Output results database |

On first run, extracts all SNP variants from the VCF and matches them against SNPedia. On subsequent runs:

- **VCF unchanged** — skips extraction entirely, uses cached variants
- **SNPedia updated** (new download) — re-runs matching against updated data
- **Both unchanged** — instant, nothing to do

### Step 3: Generate report

```bash
uv run report

# Markdown output
uv run report --format md

# Both HTML and Markdown
uv run report --format all
```

| Flag | Default | Description |
|------|---------|-------------|
| `--db` | `results.db` | Scan results database |
| `--out` | `report` | Output path (extension added automatically) |
| `--format` | `html` | `html`, `md`, or `all` |

Can be re-run anytime — useful after updating SNPedia data and re-scanning.

### Typical workflow

```bash
# First time: full pipeline
uv run download
uv run scan --vcf /mnt/dna/my_genome.vcf.gz
uv run report

# Later: update SNPedia and regenerate
uv run download                                     # only fetches changes
uv run scan --vcf /mnt/dna/my_genome.vcf.gz         # only re-matches
uv run report                                        # regenerate HTML
```

## Supported VCF formats

- GRCh38 (hg38) assembly with `chr` prefix
- DeepVariant, GATK HaplotypeCaller, or similar callers
- snpEff `ANN` and/or VEP `CSQ` annotations for rsID extraction
- Gzipped (`.vcf.gz`) or plain (`.vcf`)
- Single-sample VCFs

## The report

The HTML report includes:

- Summary cards: total matches, genotype info, notable findings, good/bad repute counts
- Interactive table: sortable by any column, filterable by repute (good/bad/notable), full-text search
- Each rsID links directly to its SNPedia page
- Findings are ranked by magnitude (see below)

Self-contained file — open in any browser, no external dependencies.

### Magnitude

Magnitude is SNPedia's crowd-sourced "importance score" for a genotype, ranging from 0 to 10. It is not a scientific metric — it's an editorial judgment by SNPedia contributors about how noteworthy a finding is. The **repute** field (good/bad) indicates the direction.

| Magnitude | Meaning | Examples |
|-----------|---------|----------|
| 0 | Common, no known significance | Most of your genome |
| 1-2 | Mildly interesting | Earwax type, eye color |
| 2-3 | Notable | Caffeine metabolism, lactose tolerance |
| 3-4 | Significant | Drug response, carrier status |
| 4+ | Medically relevant | Disease risk, pharmacogenomics |
| 6+ | High impact | BRCA mutations, APOE-e4 homozygous |
| 10 | Usually data quality flags | Not actual findings |

The report sorts by magnitude descending so the most interesting findings appear first.

## Database schemas

### `data/snpedia.db` — SNPedia mirror

| Table | Description |
|-------|-------------|
| `snps` | One row per SNP: rsid, gene, chromosome, position, orientation, gmaf, summary, raw wikitext, revision tracking |
| `genotypes` | One row per SNP+allele combo: rsid, alleles, magnitude, repute, summary |
| `genosets` | Combinations of multiple SNPs: magnitude, repute, summary, criteria |

### `results.db` — Scan results

| Table | Description |
|-------|-------------|
| `vcf_variants` | Extracted SNP variants from VCF: chromosome, position, ref, alt, genotype, quality |
| `matches` | Variants matched to SNPedia: all variant fields + rsid, gene, magnitude, repute, summary |
| `scan_meta` | Fingerprints for incremental logic, scan timestamps, stats |

## Matching logic

1. Normalize VCF chromosome (`chr1` -> `1`, `chrM` -> `MT`)
2. Look up `(chromosome, position)` in the SNPedia SNP table
3. If no match, extract rsIDs from the VCF `ID` column and VEP `CSQ` / snpEff `ANN` info fields, then look up by rsID
4. Convert VCF genotype (`0/1` with REF=A, ALT=G) to SNPedia format `(A;G)`
5. Look up genotype magnitude and repute, trying both allele orderings

## Disclaimer

This tool is for informational and educational purposes only. It is not medical advice. SNPedia data is community-curated and may contain errors. Genetic variants interact in complex ways — a single SNP rarely determines an outcome. Always consult a qualified healthcare professional or genetic counselor for medical decisions.

## License

The source code of this project is licensed under the [MIT License](LICENSE).

The `data/snpedia.db` database contains content from [SNPedia](https://www.snpedia.com/), which is available under the [Creative Commons Attribution-Noncommercial-Share Alike 3.0 United States License](https://creativecommons.org/licenses/by-nc-sa/3.0/us/) (CC BY-NC-SA 3.0 US). If you redistribute `data/snpedia.db`, it must be shared under the same CC BY-NC-SA 3.0 US license with attribution to SNPedia. The non-commercial clause means you cannot sell or use the database in commercial products without separate permission from SNPedia.
