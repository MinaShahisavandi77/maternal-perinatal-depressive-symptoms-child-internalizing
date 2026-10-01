#!/bin/bash
## ============================================================================
## Calculate PRS with PLINK 2 for LDpred2-auto and SBayesRC weights
## ============================================================================
set -euo pipefail

## ---- 1. Paths (edit these) ------------------------------------------------
PLINK="plink2"                                   # or full path, e.g. /usr/local/bin/plink2
WEIGHTS_DIR="/path/to/weights"                   # folder with the 4 weight files
OUT_DIR="/path/to/prs_output"                    # where results go

## Genotype sets: "label|plink_prefix|format"
##   format = bfile (.bed/.bim/.fam)  or  pfile (.pgen/.pvar/.psam)
## One combined file for children, mothers and fathers:
GENO_SETS=(
  "all|/path/to/genotypes/cohort_imputed_hm3|bfile"
)
## Or separate files per family member, e.g.:
# GENO_SETS=(
#   "child|/path/to/genotypes/children|bfile"
#   "mother|/path/to/genotypes/mothers|bfile"
#   "father|/path/to/genotypes/fathers|bfile"
# )

## ---- 2. Weight files: "name|file|SNP col|A1 col|BETA col" -----------------
## LDpred2 files (no header): CHR POS SNPID A1 A2 BETA  -> cols 3 4 6
## SBayesRC file:            SNP A1 BETA SE (PIP BETAlast) -> cols 1 2 3
SCORES=(
  "PPD|PPD_adjusted_weights_hapmap3_ldpred2.txt|3|4|6"
  "MDD|MDD_adjusted_weights_hapmap3_ldpred2.txt|3|4|6"
  "MDD_female|MDD_female_adjusted_weights_hapmap3_ldpred2.txt|3|4|6"
  "PPD_unique|PPD_unique_adjusted_weights_hapmap3_sbrc.txt|1|2|3"
)

mkdir -p "$OUT_DIR"
SUMMARY="$OUT_DIR/prs_overlap_summary.txt"
printf "geno_set\tscore\tn_weights\tn_used\tpct_used\n" > "$SUMMARY"

## ---- 3. Loop over genotype sets and scores --------------------------------
for geno in "${GENO_SETS[@]}"; do
  IFS='|' read -r GLABEL GPREFIX GFORMAT <<< "$geno"

  for entry in "${SCORES[@]}"; do
    IFS='|' read -r NAME FILE SNPCOL A1COL BETACOL <<< "$entry"
    W="$WEIGHTS_DIR/$FILE"

    if [[ ! -f "$W" ]]; then
      echo "WARNING: $W not found, skipping $NAME"
      continue
    fi

    ## detect a header line: if the BETA column of line 1 is not a number
    first_val=$(head -n 1 "$W" | awk -v c="$BETACOL" '{print $c}')
    if [[ "$first_val" =~ ^-?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?$ ]]; then
      HDR=""
      n_weights=$(wc -l < "$W")
    else
      HDR="header"
      n_weights=$(( $(wc -l < "$W") - 1 ))
    fi

    OUT="$OUT_DIR/prs_${NAME}_${GLABEL}"
    echo ">>> $NAME on $GLABEL (header: ${HDR:-none})"

    $PLINK --"$GFORMAT" "$GPREFIX" \
           --score "$W" "$SNPCOL" "$A1COL" "$BETACOL" $HDR \
                   cols=+scoresums list-variants ignore-dup-ids \
           --out "$OUT"

    ## overlap check: how many weight SNPs were found in the genotypes
    if [[ -f "$OUT.sscore.vars" ]]; then
      n_used=$(wc -l < "$OUT.sscore.vars")
      pct=$(awk -v u="$n_used" -v w="$n_weights" 'BEGIN{printf "%.1f", 100*u/w}')
      printf "%s\t%s\t%s\t%s\t%s\n" "$GLABEL" "$NAME" "$n_weights" "$n_used" "$pct" >> "$SUMMARY"
    fi
  done
done

echo "Done. Scores in $OUT_DIR (*.sscore); SNP overlap in $SUMMARY"
column -t "$SUMMARY"
