#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "$0")"

rm -f T2DGGI.MGI.EUR.T2D.metacarpa.txt

set +e
timeout "${METACARPA_REPRO_TIMEOUT:-5s}" ../../src/metacarpa \
  -I T2DGGI.MGI.EUR.T2D.F.ldsc.txt \
  -I T2DGGI.MGI.EUR.T2D.M.ldsc.txt \
  -O T2DGGI.MGI.EUR.T2D.metacarpa.txt \
  -t ' ' \
  -c 2 \
  -q 4 \
  -u 5 \
  -v 6 \
  -a 7 \
  -b 8 \
  -s 9 \
  -p 11 \
  -n 12 \
  -i 1 \
  -m T2DGGI.MGI.EUR.T2D.matrix.txt
status=$?
set -e

if [[ $status -ne 0 && $status -ne 124 ]]; then
  exit "$status"
fi

if [[ ! -s T2DGGI.MGI.EUR.T2D.metacarpa.txt ]]; then
  echo "METACARPA did not produce output" >&2
  exit 1
fi

awk 'BEGIN { FS = OFS = "\t" } NR == 1 || $1 == "chr1:414350:C:T" || $1 == "chr1:701697:TCA:T" { print }' \
  T2DGGI.MGI.EUR.T2D.metacarpa.txt
