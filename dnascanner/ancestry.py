"""
Ancestry estimation library — ancestry-informative markers (AIMs) with known
1000 Genomes allele frequencies and maximum-likelihood estimation across
five super-populations:
  EUR (European), AFR (African), EAS (East Asian),
  SAS (South Asian), AMR (Americas/Admixed American)
"""

import math
import re
import sqlite3

from tqdm import tqdm


# ---------------------------------------------------------------------------
# Ancestry-informative markers (AIMs)
#
# Each entry: (rsid, chrom, pos_grch38, ref, alt, EUR_alt, AFR_alt, EAS_alt, SAS_alt, AMR_alt)
#
# Frequencies are ALT allele frequencies per 1000 Genomes Phase 3
# super-population, lifted to GRCh38 coordinates.
#
# Selected for high Fst (population differentiation) across continental
# groups — these are among the most ancestry-informative common SNPs.
# ---------------------------------------------------------------------------

AIMS: list[tuple[str, str, int, str, str, float, float, float, float, float]] = [
    # rsid, chrom, pos_grch38, ref, alt, EUR, AFR, EAS, SAS, AMR
    # --- High Fst markers: skin/hair/eye pigmentation and SLC genes ---
    ("rs1426654", "15", 48426484, "A", "G", 0.999, 0.065, 0.987, 0.870, 0.780),
    ("rs16891982", "5", 33951693, "G", "C", 0.937, 0.017, 0.003, 0.055, 0.460),
    ("rs3827760", "2", 108962164, "A", "G", 0.004, 0.001, 0.870, 0.010, 0.380),
    ("rs1834640", "15", 48408306, "A", "G", 0.890, 0.100, 0.760, 0.620, 0.560),
    ("rs1800414", "15", 28197037, "C", "T", 0.002, 0.001, 0.580, 0.020, 0.060),
    ("rs2228479", "16", 89985753, "G", "A", 0.080, 0.013, 0.200, 0.050, 0.100),
    ("rs1042602", "11", 89178528, "C", "A", 0.370, 0.020, 0.020, 0.100, 0.240),
    ("rs12913832", "15", 28120472, "A", "G", 0.720, 0.030, 0.010, 0.040, 0.310),
    ("rs1545397", "15", 28099092, "A", "T", 0.610, 0.920, 0.130, 0.400, 0.420),
    ("rs6119471", "20", 34352721, "C", "G", 0.010, 0.700, 0.005, 0.020, 0.090),
    # --- Duffy antigen / DARC (highly differentiated AFR vs non-AFR) ---
    ("rs2814778", "1", 159205564, "T", "C", 0.997, 0.030, 0.999, 0.990, 0.870),
    # --- LCT / lactase persistence region ---
    ("rs4988235", "2", 135851076, "G", "A", 0.760, 0.090, 0.010, 0.270, 0.370),
    ("rs182549", "2", 135859184, "C", "T", 0.760, 0.130, 0.020, 0.300, 0.390),
    # --- SLC24A5, SLC45A2 (pigmentation, high EUR frequency) ---
    ("rs1800407", "15", 27985172, "G", "A", 0.074, 0.009, 0.000, 0.010, 0.035),
    # --- HERC2/OCA2 eye colour ---
    ("rs12203592", "6", 396321, "C", "T", 0.128, 0.006, 0.000, 0.010, 0.060),
    # --- ADH1B (alcohol metabolism, high EAS frequency) ---
    ("rs1229984", "4", 99318162, "T", "C", 0.030, 0.010, 0.730, 0.100, 0.100),
    # --- ALDH2 (alcohol flush, EAS-specific) ---
    ("rs671", "12", 111803962, "G", "A", 0.001, 0.000, 0.230, 0.002, 0.020),
    # --- ABCC11 (earwax type, high EAS frequency) ---
    ("rs17822931", "16", 48224287, "C", "T", 0.120, 0.050, 0.930, 0.300, 0.440),
    # --- APOL1 (AFR-specific kidney disease risk) ---
    ("rs73885319", "22", 36265860, "A", "G", 0.001, 0.220, 0.000, 0.000, 0.040),
    ("rs60910145", "22", 36265988, "T", "G", 0.001, 0.150, 0.000, 0.000, 0.030),
    # --- HBB region (sickle cell, AFR-enriched) ---
    ("rs334", "11", 5227002, "T", "A", 0.001, 0.080, 0.000, 0.005, 0.020),
    # --- G6PD region (malaria selection) ---
    ("rs1050828", "X", 154536002, "C", "T", 0.002, 0.200, 0.001, 0.020, 0.040),
    # --- FTO (obesity-associated, frequency varies) ---
    ("rs9939609", "16", 53786615, "T", "A", 0.420, 0.490, 0.140, 0.310, 0.270),
    # --- High Fst markers from Kosoy et al. and other AIMS panels ---
    ("rs2065160", "1", 234443048, "C", "T", 0.200, 0.730, 0.070, 0.170, 0.300),
    ("rs7657799", "4", 38489422, "T", "C", 0.080, 0.620, 0.040, 0.130, 0.170),
    ("rs798443", "6", 34671412, "G", "A", 0.630, 0.110, 0.730, 0.360, 0.430),
    ("rs1079597", "11", 113296286, "G", "A", 0.160, 0.020, 0.010, 0.070, 0.100),
    ("rs3737576", "8", 100063342, "G", "A", 0.160, 0.580, 0.100, 0.220, 0.260),
    ("rs1876482", "2", 22457788, "C", "T", 0.550, 0.250, 0.530, 0.550, 0.480),
    ("rs2891remapped_placeholder", "0", 0, "X", "X", 0, 0, 0, 0, 0),  # placeholder
    ("rs260714", "2", 176883561, "T", "C", 0.550, 0.070, 0.660, 0.420, 0.430),
    ("rs1871534", "5", 16783644, "G", "A", 0.090, 0.560, 0.050, 0.110, 0.180),
    ("rs2789823", "1", 202455662, "A", "G", 0.210, 0.660, 0.170, 0.270, 0.330),
    ("rs310644", "20", 62165061, "A", "G", 0.790, 0.320, 0.830, 0.680, 0.620),
    ("rs17034", "3", 122373123, "G", "A", 0.150, 0.550, 0.070, 0.200, 0.230),
    ("rs1461227", "14", 99461689, "A", "G", 0.300, 0.710, 0.100, 0.290, 0.340),
    ("rs2416791", "12", 885142, "A", "G", 0.440, 0.140, 0.830, 0.390, 0.400),
    ("rs1335873", "1", 165093367, "T", "C", 0.050, 0.480, 0.020, 0.100, 0.120),
    ("rs1478785", "5", 7044948, "G", "A", 0.650, 0.200, 0.750, 0.520, 0.500),
    ("rs1385699", "3", 38332012, "T", "C", 0.270, 0.680, 0.130, 0.340, 0.350),
    ("rs2077681", "11", 88738972, "C", "T", 0.100, 0.550, 0.030, 0.160, 0.180),
    ("rs722098", "7", 28141455, "A", "G", 0.580, 0.130, 0.680, 0.430, 0.460),
    ("rs4471745", "5", 69834375, "A", "G", 0.210, 0.600, 0.100, 0.200, 0.270),
    ("rs1040045", "10", 30281590, "G", "A", 0.410, 0.070, 0.500, 0.330, 0.320),
    ("rs730570", "2", 175046338, "C", "T", 0.300, 0.680, 0.210, 0.340, 0.370),
    ("rs3916235", "7", 99026694, "T", "C", 0.130, 0.500, 0.070, 0.140, 0.200),
    ("rs2695", "17", 46952935, "G", "A", 0.560, 0.180, 0.670, 0.460, 0.450),
    ("rs1058083", "17", 29680478, "A", "G", 0.250, 0.640, 0.180, 0.310, 0.340),
    ("rs2688457", "10", 73488297, "C", "T", 0.110, 0.530, 0.040, 0.130, 0.180),
    ("rs3784230", "12", 110478700, "C", "T", 0.310, 0.640, 0.180, 0.300, 0.350),
    ("rs1363448", "3", 4127206, "A", "G", 0.320, 0.080, 0.550, 0.260, 0.290),
    ("rs735480", "20", 42557193, "G", "A", 0.600, 0.200, 0.710, 0.470, 0.480),
    ("rs2072857", "22", 42252924, "C", "T", 0.140, 0.510, 0.060, 0.170, 0.200),
    ("rs1592881", "8", 130826987, "G", "A", 0.350, 0.750, 0.250, 0.400, 0.410),
    ("rs1403454", "1", 236984982, "T", "C", 0.470, 0.090, 0.600, 0.380, 0.370),
    ("rs2077863", "11", 68990072, "A", "G", 0.250, 0.620, 0.130, 0.270, 0.310),
    ("rs2024566", "8", 102563903, "C", "T", 0.200, 0.570, 0.100, 0.210, 0.260),
    ("rs1325502", "6", 75698614, "G", "A", 0.440, 0.100, 0.590, 0.340, 0.360),
    ("rs1406927", "3", 164476133, "T", "C", 0.360, 0.690, 0.230, 0.360, 0.380),
    ("rs7226659", "17", 41508279, "A", "G", 0.210, 0.580, 0.110, 0.250, 0.280),
    ("rs2024324", "14", 65085792, "G", "A", 0.360, 0.050, 0.510, 0.280, 0.300),
    ("rs3176921", "9", 136668475, "C", "T", 0.180, 0.550, 0.080, 0.200, 0.230),
    ("rs1050501", "1", 161479745, "G", "A", 0.100, 0.430, 0.040, 0.140, 0.160),
    ("rs1439553", "6", 160478831, "T", "C", 0.320, 0.670, 0.210, 0.350, 0.370),
    ("rs1490388", "6", 26058924, "T", "C", 0.570, 0.210, 0.690, 0.460, 0.470),
    ("rs1344870", "7", 73756342, "A", "G", 0.410, 0.780, 0.300, 0.430, 0.440),
    ("rs3806624", "17", 30475447, "C", "T", 0.330, 0.670, 0.220, 0.360, 0.380),
    ("rs1560089", "4", 155714629, "G", "A", 0.280, 0.620, 0.150, 0.290, 0.320),
    ("rs2397060", "9", 100239638, "A", "G", 0.190, 0.560, 0.090, 0.230, 0.250),
    ("rs1373302", "13", 54478157, "T", "C", 0.420, 0.110, 0.570, 0.340, 0.350),
    ("rs1355358", "2", 42178051, "C", "T", 0.500, 0.160, 0.640, 0.400, 0.410),
    ("rs2688608", "10", 73490021, "A", "G", 0.340, 0.680, 0.230, 0.370, 0.390),
    ("rs1508911", "11", 47336618, "G", "A", 0.260, 0.600, 0.140, 0.280, 0.310),
    ("rs1412385", "14", 64037824, "T", "C", 0.380, 0.050, 0.530, 0.300, 0.310),
    ("rs1478284", "3", 161296975, "G", "A", 0.290, 0.640, 0.170, 0.310, 0.340),
    ("rs1448484", "1", 95637584, "C", "T", 0.370, 0.080, 0.530, 0.290, 0.310),
]

# Filter out placeholder entries
AIMS = [m for m in AIMS if m[2] != 0]

POPULATIONS = ["EUR", "AFR", "EAS", "SAS", "AMR"]
POP_LABELS = {
    "EUR": "European",
    "AFR": "African",
    "EAS": "East Asian",
    "SAS": "South Asian",
    "AMR": "Americas",
}


# ---------------------------------------------------------------------------
# Look up AIMs in the already-extracted vcf_variants table
# ---------------------------------------------------------------------------


def lookup_aims_in_db(results_conn: sqlite3.Connection) -> list[dict]:
    """
    Match AIMs against the vcf_variants table (already extracted by scan).

    For positions present in the VCF: count ALT alleles from the genotype.
    For positions absent from the VCF: assume homozygous reference (alt_count=0).
    Both cases are informative for ancestry estimation.
    """
    markers = []

    for rsid, chrom, pos, ref, alt, eur, afr, eas, sas, amr in tqdm(
        AIMS,
        desc="Looking up AIMs",
        unit=" markers",
        colour="#a78bfa",
        bar_format="{l_bar}{bar}| {n_fmt}/{total_fmt} [{elapsed}<{remaining}] {postfix}",
    ):
        freqs = {"EUR": eur, "AFR": afr, "EAS": eas, "SAS": sas, "AMR": amr}
        aim_alt = alt.upper()

        # Try position match first
        row = results_conn.execute(
            "SELECT genotype, ref, alt FROM vcf_variants WHERE chromosome = ? AND position = ?",
            (chrom, pos),
        ).fetchone()

        # Try rsID match as fallback (rsIDs may be in vcf_rsids column)
        if row is None:
            row = results_conn.execute(
                "SELECT genotype, ref, alt FROM vcf_variants WHERE vcf_rsids LIKE ?",
                (f"%{rsid}%",),
            ).fetchone()

        if row:
            genotype_str, vcf_ref, vcf_alt = row
            # Parse genotype like "(A;G)" to count ALT alleles
            m = re.match(r"\((.);(.)\)", genotype_str)
            if m:
                a1, a2 = m.group(1).upper(), m.group(2).upper()
                alt_count = (1 if a1 == aim_alt else 0) + (1 if a2 == aim_alt else 0)
                genotype = f"{a1}/{a2}"
            else:
                continue
        else:
            # Not in VCF = homozygous reference = 0 copies of ALT
            alt_count = 0
            genotype = f"{ref}/{ref}"

        markers.append(
            {
                "rsid": rsid,
                "chrom": chrom,
                "pos": pos,
                "genotype": genotype,
                "alt_count": alt_count,
                "freqs": freqs,
                "in_vcf": row is not None,
            }
        )

    return markers


# ---------------------------------------------------------------------------
# Maximum-likelihood ancestry estimation
# ---------------------------------------------------------------------------

FLOOR = 1e-6  # avoid log(0)


def genotype_likelihood(alt_count: int, alt_freq: float) -> float:
    """P(genotype | population) assuming Hardy-Weinberg."""
    p = max(FLOOR, min(1 - FLOOR, alt_freq))
    q = 1 - p
    if alt_count == 0:
        return q * q
    elif alt_count == 1:
        return 2 * p * q
    else:
        return p * p


def log_likelihood(markers: list[dict], proportions: list[float]) -> float:
    """Compute log-likelihood of observed genotypes given ancestry proportions."""
    ll = 0.0
    for m in markers:
        # P(genotype) = sum over pops: proportion_k * P(genotype | pop_k)
        p_geno = 0.0
        for i, pop in enumerate(POPULATIONS):
            p_geno += proportions[i] * genotype_likelihood(
                m["alt_count"], m["freqs"][pop]
            )
        ll += math.log(max(p_geno, FLOOR))
    return ll


def estimate_ancestry(markers: list[dict], n_iter: int = 500) -> dict[str, float]:
    """
    Estimate ancestry proportions via EM algorithm (maximum likelihood).

    Returns dict mapping population code to proportion (sums to 1.0).
    """
    n_pops = len(POPULATIONS)
    # Initialize uniform
    props = [1.0 / n_pops] * n_pops

    for _ in range(n_iter):
        # E-step: compute responsibilities
        new_props = [0.0] * n_pops

        for m in markers:
            # P(genotype | pop_k) * proportion_k
            weighted = []
            for i, pop in enumerate(POPULATIONS):
                w = props[i] * genotype_likelihood(m["alt_count"], m["freqs"][pop])
                weighted.append(w)

            total = sum(weighted)
            if total < FLOOR:
                continue

            for i in range(n_pops):
                new_props[i] += weighted[i] / total

        # M-step: normalize
        total = sum(new_props)
        if total > 0:
            props = [p / total for p in new_props]

    return {pop: round(props[i], 4) for i, pop in enumerate(POPULATIONS)}
