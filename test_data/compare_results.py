#!/usr/bin/env python3
"""
Compare two metacarpa output files.

Usage:
    python3 compare_results.py A.txt B.txt [label_a] [label_b] [out_prefix]
    python3 compare_results.py A.txt B.txt A B results/compare \
        --matrix-a A.matrix.txt --matrix-b B.matrix.txt

Produces:
    {out_prefix}.png   — 4-panel figure
    {out_prefix}.txt   — text summary (also printed to stdout)
"""

import sys
import math
import re
import argparse
from collections import defaultdict


def parse_log10p(s):
    """Parse a p-value string (including extreme exponents) → -log10(p).
    Returns None for missing / single-study (-1) values."""
    s = s.strip()
    if s in ("-1", "-1.0", "0", "0.0", ""):
        return None
    m = re.match(r"([+-]?\d+\.?\d*)[eE]([+-]?\d+)$", s)
    if m:
        log10_mantissa = math.log10(abs(float(m.group(1))))
        return -(log10_mantissa + int(m.group(2)))   # -log10(p)
    try:
        v = float(s)
        if v <= 0:
            return None
        return -math.log10(v)
    except ValueError:
        return None


def read_output(path):
    rows = {}
    with open(path) as f:
        header = f.readline().rstrip("\n").split("\t")
        col = {h: i for i, h in enumerate(header)}
        for line in f:
            fields = line.rstrip("\n").split("\t")
            if not fields[0]:
                continue
            rows[fields[0]] = fields
    return rows, col


def read_matrix(path):
    """Read a metacarpa matrix file → list of floats (upper triangle, row-major)."""
    if path is None:
        return None
    values = []
    try:
        with open(path) as f:
            for line in f:
                line = line.strip()
                if line:
                    try:
                        values.append(float(line))
                    except ValueError:
                        pass
    except FileNotFoundError:
        return None
    return values


def summarise(lines):
    for line in lines:
        print(line)


def run(file_a, file_b, label_a, label_b, out_prefix, matrix_a, matrix_b):
    rows_a, col_a = read_output(file_a)
    rows_b, col_b = read_output(file_b)

    common = sorted(set(rows_a) & set(rows_b))
    multi = []   # (lp_a, lp_b, lz_a, lz_b, beta_b, z_b) for multi-study variants

    for rsid in common:
        ra, rb = rows_a[rsid], rows_b[rsid]
        lp_a = parse_log10p(ra[col_a["p_corrected"]])
        lp_b = parse_log10p(rb[col_b["p_corrected"]])
        if lp_a is None or lp_b is None:
            continue
        try:
            z_a = float(ra[col_a["z"]])
            z_b = float(rb[col_b["z"]])
            beta_b = float(rb[col_b["beta"]])
        except (ValueError, KeyError):
            continue
        multi.append((lp_a, lp_b, abs(z_a), abs(z_b), beta_b, z_b))

    n = len(multi)
    if n == 0:
        print("No matched multi-study variants found.")
        return

    deltas   = [abs(b - a) for a, b, *_ in multi]
    lp_a_all = [a for a, *_ in multi]
    lp_b_all = [b for _, b, *_ in multi]

    n_sign_wrong = sum(
        1 for *_, beta, z in multi
        if (beta >= 0 and z < 0) or (beta < 0 and z >= 0)
    )
    n_large  = sum(1 for d in deltas if d > 0.5)
    n_huge   = sum(1 for d in deltas if d > 1.0)
    mean_d   = sum(deltas) / n
    median_d = sorted(deltas)[n // 2]
    max_d    = max(deltas)

    summary = [
        "",
        f"{'─'*52}",
        f"  Comparison: {label_a}  vs  {label_b}",
        f"{'─'*52}",
        f"  Matched multi-study variants:  {n:,}",
        f"  Mean  |Δlog₁₀(p_corrected)|:  {mean_d:.4f}",
        f"  Median|Δlog₁₀(p_corrected)|:  {median_d:.4f}",
        f"  Max   |Δlog₁₀(p_corrected)|:  {max_d:.4f}",
        f"  Variants with |Δ| > 0.5:       {n_large:,}",
        f"  Variants with |Δ| > 1.0:       {n_huge:,}",
        f"  z sign wrong in {label_b}:  {n_sign_wrong:,}",
    ]

    mat_a = read_matrix(matrix_a)
    mat_b = read_matrix(matrix_b)
    if mat_a or mat_b:
        summary.append("")
        summary.append("  Tetrachoric correlation matrices:")
        if mat_a:
            summary.append(f"    {label_a}: {mat_a}")
        if mat_b:
            summary.append(f"    {label_b}: {mat_b}")
        if mat_a and mat_b and len(mat_a) == len(mat_b):
            diffs = [abs(a - b) for a, b in zip(mat_a, mat_b)]
            summary.append(f"    Max |Δρ|: {max(diffs):.6f}")

    summary.append(f"{'─'*52}")
    summarise(summary)

    # ── Plots ──────────────────────────────────────────────────────────────
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import numpy as np
    except ImportError:
        print("matplotlib not available — skipping plots")
        return

    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    fig.suptitle(f"metacarpa comparison: {label_a} vs {label_b}", fontsize=13)

    lp_a = np.array(lp_a_all)
    lp_b = np.array(lp_b_all)
    delta = lp_b - lp_a   # positive = more significant in B

    cap = 50   # cap display at -log10(p) = 50 to keep extreme variants from dominating

    # ── 1. p-value scatter ─────────────────────────────────────────────────
    ax = axes[0, 0]
    x = np.clip(lp_a, 0, cap)
    y = np.clip(lp_b, 0, cap)
    ax.scatter(x, y, s=1, alpha=0.3, color="steelblue", rasterized=True)
    lim = max(x.max(), y.max()) * 1.05
    ax.plot([0, lim], [0, lim], "r--", lw=0.8, label="identity")
    ax.set_xlabel(f"−log₁₀(p)  [{label_a}]")
    ax.set_ylabel(f"−log₁₀(p)  [{label_b}]")
    ax.set_title("p_corrected scatter (capped at 50)")
    ax.legend(fontsize=8)

    # ── 2. QQ plot ─────────────────────────────────────────────────────────
    ax = axes[0, 1]
    expected = -np.log10(np.arange(1, n + 1) / (n + 1))
    for vals, label, color in [
        (lp_a, label_a, "steelblue"),
        (lp_b, label_b, "darkorange"),
    ]:
        observed = np.sort(np.clip(vals, 0, cap))[::-1]
        ax.scatter(np.clip(expected, 0, cap), observed,
                   s=1, alpha=0.5, color=color, label=label, rasterized=True)
    top = min(cap, max(expected.max(), lp_a.max(), lp_b.max())) * 1.05
    ax.plot([0, top], [0, top], "r--", lw=0.8)
    ax.set_xlabel("Expected −log₁₀(p)")
    ax.set_ylabel("Observed −log₁₀(p)")
    ax.set_title("QQ plot (capped at 50)")
    ax.legend(fontsize=8, markerscale=5)

    # ── 3. Δlog₁₀(p) histogram ────────────────────────────────────────────
    ax = axes[1, 0]
    d_clipped = np.clip(delta, -10, 10)
    ax.hist(d_clipped, bins=100, color="steelblue", edgecolor="none")
    ax.axvline(0, color="red", lw=0.8, linestyle="--")
    ax.set_xlabel(f"Δlog₁₀(p)  [{label_b} − {label_a}]  (clipped ±10)")
    ax.set_ylabel("Variants")
    ax.set_title("Distribution of p_corrected change")
    ax.text(0.97, 0.95, f"mean={delta.mean():.3f}\nmedian={np.median(delta):.3f}",
            transform=ax.transAxes, ha="right", va="top", fontsize=8)

    # ── 4. |z| scatter ─────────────────────────────────────────────────────
    ax = axes[1, 1]
    az_a = np.array([row[2] for row in multi])   # |z_a|
    az_b = np.array([row[3] for row in multi])   # |z_b|
    z_cap = min(50, max(az_a.max(), az_b.max()) * 1.05)
    ax.scatter(np.clip(az_a, 0, z_cap), np.clip(az_b, 0, z_cap),
               s=1, alpha=0.3, color="steelblue", rasterized=True)
    ax.plot([0, z_cap], [0, z_cap], "r--", lw=0.8)
    ax.set_xlabel(f"|z|  [{label_a}]")
    ax.set_ylabel(f"|z|  [{label_b}]")
    ax.set_title("|z| scatter (capped at 50)")

    plt.tight_layout()
    out_png = out_prefix + ".png"
    plt.savefig(out_png, dpi=150)
    print(f"  Plot saved to: {out_png}")


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("file_a")
    p.add_argument("file_b")
    p.add_argument("label_a", nargs="?", default="A")
    p.add_argument("label_b", nargs="?", default="B")
    p.add_argument("out_prefix", nargs="?", default="compare")
    p.add_argument("--matrix-a", default=None)
    p.add_argument("--matrix-b", default=None)
    args = p.parse_args()
    run(args.file_a, args.file_b, args.label_a, args.label_b,
        args.out_prefix, args.matrix_a, args.matrix_b)


if __name__ == "__main__":
    main()
