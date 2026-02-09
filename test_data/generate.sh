#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"

METACARPA=../src/metacarpa
LIB_PATH=../src/lib

if [ ! -x "$METACARPA" ]; then
  echo "ERROR: metacarpa binary not found at $METACARPA. Run 'make -C ../src' first."
  exit 1
fi

# ============================================================
# 1. Download GIANT 2022 Height EUR summary statistics
# ============================================================
GIANT_BASE="https://portals.broadinstitute.org/collaboration/giant/images"
WITH_UKB_URL="${GIANT_BASE}/2/21/GIANT_HEIGHT_YENGO_2022_GWAS_SUMMARY_STATS_EUR.gz"
NO_UKB_URL="${GIANT_BASE}/4/40/GIANT_HEIGHT_YENGO_2022_GWAS_SUMMARY_STATS_EUR_EXCLUDING_UKB.gz"

for f in height_eur_with_ukb.gz height_eur_no_ukb.gz; do
  if [ -f "$f" ]; then
    echo "SKIP: $f already exists"
  else
    case "$f" in
      *with_ukb*) url="$WITH_UKB_URL" ;;
      *no_ukb*)   url="$NO_UKB_URL" ;;
    esac
    echo "Downloading $f ..."
    wget -q -O "$f" "$url"
  fi
done

# ============================================================
# 2. Filter and subsample
#    - Keep MAF > 0.05 in with-UKB cohort
#    - Randomly subsample to 1M variants
#    - Filter no-UKB to matching SNPIDs
#    - Sort both by chr:pos
# ============================================================
if [ -f with_ukb_filtered.tsv ] && [ -f no_ukb_filtered.tsv ]; then
  echo "SKIP: filtered files already exist"
else
  echo "Filtering and subsampling..."

  # Filter with-UKB for MAF > 0.05, sample 1M, sort
  zcat height_eur_with_ukb.gz \
    | awk -F'\t' 'NR==1 || ($9 > 0.05 && $9 < 0.95)' \
    | awk 'NR==1{print; next} {print | "shuf -n 1000000"}' \
    | sort -t$'\t' -k3,3n -k4,4n \
    > with_ukb_filtered.tsv

  # Extract SNPID list
  awk -F'\t' 'NR>1{print $1}' with_ukb_filtered.tsv > snpids_keep.txt

  # Filter no-UKB to matching SNPIDs, sort
  zcat height_eur_no_ukb.gz \
    | awk -F'\t' 'NR==1{print; next} FNR==NR{ids[$1]; next} ($1 in ids)' \
        snpids_keep.txt - \
    | sort -t$'\t' -k3,3n -k4,4n \
    > no_ukb_filtered.tsv

  rm -f snpids_keep.txt
  echo "  with_ukb_filtered.tsv: $(wc -l < with_ukb_filtered.tsv) lines"
  echo "  no_ukb_filtered.tsv: $(wc -l < no_ukb_filtered.tsv) lines"
fi

# ============================================================
# 3. Run METACARPA — current method
# ============================================================
TAB=$'\t'
METACARPA_ARGS="-t '${TAB}' -c 3 -q 4 -u 5 -v 6 -a 7 -b 8 -s 9 -p 10 -n 11 -i 1"

if [ -f output_current.tsv ]; then
  echo "SKIP: output_current.tsv already exists"
else
  echo "Running METACARPA (current method)..."
  LD_LIBRARY_PATH="$LIB_PATH" "$METACARPA" \
    -I with_ukb_filtered.tsv,423955 \
    -I no_ukb_filtered.tsv,253288 \
    -O output_current.tsv \
    -t "$TAB" -c 3 -q 4 -u 5 -v 6 -a 7 -b 8 -s 9 -p 10 -n 11 -i 1
fi

# ============================================================
# 4. Run METACARPA — beta-sign method (with allele harmonization)
# ============================================================
if [ -f output_betasign_fixed.tsv ]; then
  echo "SKIP: output_betasign_fixed.tsv already exists"
else
  echo "Running METACARPA (beta-sign method)..."
  LD_LIBRARY_PATH="$LIB_PATH" "$METACARPA" \
    -I with_ukb_filtered.tsv,423955 \
    -I no_ukb_filtered.tsv,253288 \
    -O output_betasign_fixed.tsv \
    -t "$TAB" -c 3 -q 4 -u 5 -v 6 -a 7 -b 8 -s 9 -p 10 -n 11 -i 1 \
    --use-beta-sign
fi

# ============================================================
# 5. Run null simulations
# ============================================================
if [ -f null_sim_results.csv ]; then
  echo "SKIP: null_sim_results.csv already exists"
else
  echo "Running null simulations..."
  Rscript run_null_sim.R
fi

echo ""
echo "Done. Open compare_methods.ipynb to explore results."
