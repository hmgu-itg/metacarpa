#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "$0")"

./run_repro.sh >/tmp/metacarpa_z_sign_repro.tsv

awk '
BEGIN { FS = "\t"; ok = 1 }
NR == 1 {
  for (i = 1; i <= NF; i++) col[$i] = i
  next
}
$1 == "chr1:414350:C:T" {
  if ($(col["beta"]) <= 0 || $(col["z"]) <= 0) {
    printf("Expected chr1:414350:C:T beta and z to be positive, got beta=%s z=%s\n", $(col["beta"]), $(col["z"])) > "/dev/stderr"
    ok = 0
  }
}
$1 == "chr1:701697:TCA:T" {
  if ($(col["beta"]) >= 0 || $(col["z"]) >= 0) {
    printf("Expected chr1:701697:TCA:T beta and z to be negative, got beta=%s z=%s\n", $(col["beta"]), $(col["z"])) > "/dev/stderr"
    ok = 0
  }
}
END { exit(ok ? 0 : 1) }
' T2DGGI.MGI.EUR.T2D.metacarpa.txt
