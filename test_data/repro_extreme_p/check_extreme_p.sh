#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "$0")"

rm -f out.txt

../../src/metacarpa \
  -I study1.txt \
  -I study2.txt \
  -O out.txt \
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
  -m matrix.txt >/tmp/metacarpa_extreme_p.stdout 2>/tmp/metacarpa_extreme_p.stderr

awk '
BEGIN { FS = "\t"; ok = 0 }
NR == 1 {
  for (i = 1; i <= NF; i++) col[$i] = i
  next
}
$1 == "chr1:1000:A:G" {
  if ($(col["z"]) > 60 && $(col["p_corrected"]) != "0") ok = 1
  else {
    printf("Extreme p-value check failed: z=%s p_corrected=%s\n", $(col["z"]), $(col["p_corrected"])) > "/dev/stderr"
  }
}
END { exit(ok ? 0 : 1) }
' out.txt
