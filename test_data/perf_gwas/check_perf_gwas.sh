#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "$0")"

METACARPA=../../src/metacarpa
METACARPA_BASELINE=${METACARPA_BASELINE:-}
N_VARIANTS=${N_VARIANTS:-50000}
SLOWDOWN_THRESHOLD=${SLOWDOWN_THRESHOLD:-1.2}  # fail if >20% slower than baseline

if [ ! -x "$METACARPA" ]; then
  echo "ERROR: metacarpa binary not found at $METACARPA. Run 'make -C ../../src' first." >&2
  exit 1
fi

if [ -n "$METACARPA_BASELINE" ] && [ ! -x "$METACARPA_BASELINE" ]; then
  echo "ERROR: baseline binary not found at $METACARPA_BASELINE" >&2
  exit 1
fi

echo "Generating synthetic GWAS data (${N_VARIANTS} variants)..."
python3 generate_gwas.py "$N_VARIANTS"

run_metacarpa() {
  local binary="$1"
  local outfile="$2"
  local label="$3"
  local start end elapsed_ms

  rm -f "$outfile" .matrix.txt

  start=$(date +%s%N)
  "$binary" \
    -I study1.txt \
    -I study2.txt \
    -O "$outfile" \
    -t ' ' \
    -c 2 -q 4 -u 5 -v 6 -a 7 -b 8 -s 9 -p 11 -n 12 -i 1 \
    >/tmp/metacarpa_perf_"${label}".stdout 2>/tmp/metacarpa_perf_"${label}".stderr
  end=$(date +%s%N)

  elapsed_ms=$(( (end - start) / 1000000 ))
  echo "  $label: ${elapsed_ms}ms"
  echo "$elapsed_ms"
}

echo "Running metacarpa (current)..."
current_ms=$(run_metacarpa "$METACARPA" out_current.txt "current" | tail -1)

if [ -n "$METACARPA_BASELINE" ]; then
  echo "Running metacarpa (baseline)..."
  baseline_ms=$(run_metacarpa "$METACARPA_BASELINE" out_baseline.txt "baseline" | tail -1)

  ratio=$(awk "BEGIN { printf \"%.2f\", $current_ms / $baseline_ms }")
  echo "  Speedup ratio: ${ratio}x (current / baseline; <1 = faster)"

  ok=$(awk "BEGIN { print ($ratio <= $SLOWDOWN_THRESHOLD) ? \"yes\" : \"no\" }")
  if [ "$ok" != "yes" ]; then
    echo "FAIL: current binary is ${ratio}x slower than baseline (threshold: ${SLOWDOWN_THRESHOLD}x)" >&2
    exit 1
  fi
fi

# Correctness checks on current output
if [ ! -s out_current.txt ]; then
  echo "FAIL: output file is empty" >&2
  exit 1
fi

awk '
BEGIN {
  FS = "\t"
  ok = 1
  n_zeros = 0
  n_extreme_ok = 0
}
NR == 1 {
  for (i = 1; i <= NF; i++) col[$i] = i
  next
}
{
  p = $(col["p_corrected"])
  if (p == "0") {
    n_zeros++
    printf("UNDERFLOW: %s has p_corrected=0\n", $1) > "/dev/stderr"
    ok = 0
  }
}
$1 ~ /^chr1:(100|200|300|400|500):/ && $(col["p_corrected"]) != "-1" {
  p = $(col["p_corrected"])
  if (p != "0") n_extreme_ok++
  else {
    printf("UNDERFLOW in extreme-p variant: %s p_corrected=%s\n", $1, p) > "/dev/stderr"
    ok = 0
  }
}
END {
  if (n_extreme_ok < 5) {
    printf("FAIL: only %d/5 extreme-p variants have valid p_corrected\n", n_extreme_ok) > "/dev/stderr"
    ok = 0
  } else {
    printf("  Extreme-p variants: %d/5 passed precision check\n", n_extreme_ok)
  }
  if (n_zeros == 0) {
    printf("  Underflow check: passed (no p_corrected=0 in %d variants)\n", NR - 1)
  }
  exit(ok ? 0 : 1)
}
' out_current.txt

rm -f study1.txt study2.txt out_current.txt out_baseline.txt .matrix.txt
