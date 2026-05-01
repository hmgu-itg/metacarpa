#!/usr/bin/env python3
"""Generate synthetic GWAS summary statistics for metacarpa performance/regression testing.

Produces study1.txt and study2.txt in the current directory.
Format: space-delimited, 12 columns matching the existing test convention:
  ID CHR POS RPOS AE AO EAF BETA SE Z P NEFF

Extreme-p variants at chr1:100-500 exercise the arbitrary-precision fallback.
Regular variants span chromosomes 1-22, sorted by chr:pos.

Usage: python3 generate_gwas.py [n_variants]
"""

import math
import random
import sys

N = int(sys.argv[1]) if len(sys.argv) > 1 else 50000
N1 = 50000
N2 = 30000
SEED = 42

# Hardcoded extreme p-values: representative points along the precision tiers.
# p=1e-50   → z≈15.2,  within double range, fast long-double path
# p=1e-300  → z≈37.2,  near double underflow, long double handles it
# p=1e-1000 → z≈67.9,  exceeds double, within long double; output p triggers asymptotic
# p=1e-4000 → z≈135,   deep long double territory
# p=1e-4900 → z≈150,   near long double limit (~3.4e-4932 min normal)
# Note: 1e-5000 is below the minimum representable long double subnormal (~3.6e-4951)
#       and causes stold() to throw std::out_of_range in metacarpa's first pass.
EXTREME = [
    (1, 100, "A", "G", 0.30,  1.0, 0.02, "1e-50"),
    (1, 200, "C", "T", 0.40, -0.8, 0.02, "1e-300"),
    (1, 300, "G", "A", 0.50,  1.2, 0.02, "1e-1000"),
    (1, 400, "T", "C", 0.20,  0.9, 0.02, "1e-4000"),
    (1, 500, "A", "C", 0.60, -1.5, 0.02, "1e-4900"),
]


def two_sided_p(z):
    p = math.erfc(abs(z) / math.sqrt(2))
    return max(p, 5e-324)


def write_study(fname, shared_variants, neff, beta_seed):
    rng = random.Random(beta_seed)
    rows = []

    for chrom, pos, ae, ao, eaf, beta, se, p_str in EXTREME:
        z = beta / se
        vid = f"chr{chrom}:{pos}:{ae}:{ao}"
        rows.append((chrom, pos,
                      f"{vid} {chrom} {pos} {pos} {ae} {ao} {eaf:.4f}"
                      f" {beta:.6f} {se:.6f} {z:.4f} {p_str} {neff}\n"))

    for chrom, pos, ae, ao, eaf in shared_variants:
        beta = rng.gauss(0, 0.03)
        se = 1.0 / math.sqrt(2 * eaf * (1 - eaf) * neff)
        z = beta / se
        p = two_sided_p(z)
        vid = f"chr{chrom}:{pos}:{ae}:{ao}"
        rows.append((chrom, pos,
                      f"{vid} {chrom} {pos} {pos} {ae} {ao} {eaf:.4f}"
                      f" {beta:.6f} {se:.6f} {z:.4f} {p:.6e} {neff}\n"))

    rows.sort(key=lambda r: (r[0], r[1]))

    with open(fname, "w") as f:
        f.write("ID CHR POS RPOS AE AO EAF BETA SE Z P NEFF\n")
        for _, _, line in rows:
            f.write(line)


def main():
    rng = random.Random(SEED)
    allele_pairs = [("A", "G"), ("C", "T"), ("G", "A"), ("T", "C")]
    chroms = list(range(1, 23))

    per_chrom = N // 22
    remainder = N % 22
    shared = []
    for i, chrom in enumerate(chroms):
        n_this = per_chrom + (1 if i < remainder else 0)
        positions = sorted({rng.randint(600_000, 250_000_000) for _ in range(n_this * 2)})[:n_this]
        for pos in positions:
            ae, ao = rng.choice(allele_pairs)
            eaf = round(rng.uniform(0.05, 0.95), 4)
            shared.append((chrom, pos, ae, ao, eaf))

    write_study("study1.txt", shared, N1, beta_seed=SEED + 1)
    write_study("study2.txt", shared, N2, beta_seed=SEED + 2)
    print(f"Generated {len(shared) + len(EXTREME)} variants per study "
          f"({len(EXTREME)} extreme-p, {len(shared)} regular).")


if __name__ == "__main__":
    main()
