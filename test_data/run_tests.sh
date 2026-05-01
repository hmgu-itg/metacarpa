#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "$0")"

PASS=0
FAIL=0

run_test() {
  local name="$1"
  local script="$2"
  if bash "$script"; then
    echo "PASS: $name"
    PASS=$((PASS + 1))
  else
    echo "FAIL: $name"
    FAIL=$((FAIL + 1))
  fi
}

run_test "extreme_p"  repro_extreme_p/check_extreme_p.sh
run_test "z_sign"     repro_z_sign/check_z_sign.sh
run_test "perf_gwas"  perf_gwas/check_perf_gwas.sh

echo ""
echo "Results: $PASS passed, $FAIL failed"
[ "$FAIL" -eq 0 ]
