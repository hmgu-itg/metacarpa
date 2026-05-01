#!/usr/bin/env bash
# Compare two metacarpa binaries on synthetic GWAS data.
# Runs both, then calls compare_results.py for metrics and plots.
#
# Usage:
#   METACARPA_A=/path/to/old  METACARPA_B=/path/to/new  bash compare.sh
#
# Required env vars:
#   METACARPA_A   path to first binary
#   METACARPA_B   path to second binary  (defaults to ../src/metacarpa)
#
# Optional env vars:
#   LABEL_A       label for first binary   (default: "A")
#   LABEL_B       label for second binary  (default: "B")
#   N_VARIANTS    synthetic variant count  (default: 50000)
#   OUT           output prefix for .png/.txt  (default: compare_out)
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"

METACARPA_A=${METACARPA_A:?Set METACARPA_A to the first binary path}
METACARPA_B=${METACARPA_B:-"$SCRIPT_DIR/../src/metacarpa"}
LABEL_A=${LABEL_A:-A}
LABEL_B=${LABEL_B:-B}
N_VARIANTS=${N_VARIANTS:-10000}
OUT=${OUT:-compare_out}

for bin in "$METACARPA_A" "$METACARPA_B"; do
  [ -x "$bin" ] || { echo "ERROR: not executable: $bin" >&2; exit 1; }
done

WORKDIR=$(mktemp -d)
trap 'rm -rf "$WORKDIR"' EXIT

echo "Generating ${N_VARIANTS} synthetic variants..." >&2
( cd "$WORKDIR" && python3 "$SCRIPT_DIR/perf_gwas/generate_gwas.py" "$N_VARIANTS" >&2 )

run_binary() {
  local binary="$1" label="$2" outfile="$3"
  echo "Running ${label}..." >&2
  ( cd "$WORKDIR" && "$binary" \
      -I study1.txt -I study2.txt -O "$outfile" \
      -t ' ' -c 2 -q 4 -u 5 -v 6 -a 7 -b 8 -s 9 -p 11 -n 12 -i 1 \
      2>/dev/null ) || true
}

run_binary "$METACARPA_A" "$LABEL_A" "out_a.txt"
run_binary "$METACARPA_B" "$LABEL_B" "out_b.txt"

[ "$(dirname "$OUT")" != "." ] && mkdir -p "$(dirname "$OUT")"

# Matrix files can be passed explicitly if wanted:
#   MATRIX_A=/path/to/a.matrix.txt  MATRIX_B=/path/to/b.matrix.txt  bash compare.sh
MATRIX_A=${MATRIX_A:-}
MATRIX_B=${MATRIX_B:-}
matrix_args=()
[ -n "$MATRIX_A" ] && matrix_args+=(--matrix-a "$MATRIX_A")
[ -n "$MATRIX_B" ] && matrix_args+=(--matrix-b "$MATRIX_B")

python3 "$SCRIPT_DIR/compare_results.py" \
  "$WORKDIR/out_a.txt" "$WORKDIR/out_b.txt" \
  "$LABEL_A" "$LABEL_B" "$OUT" \
  "${matrix_args[@]}"
