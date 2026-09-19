#!/bin/bash
# Munge BRaVa exome meta-analysis gene burden results into the genebass gene-burden layout,
# then (behind --stage only) publish the served files to the daly bucket.
#
# Source: gs://daly-genetics-results/raw/brava/gene/, a copy of the consortium's requester-pays
#   gs://brava-meta-analysis/gene/. Phenotype names and sample sizes come from the
#   brava_pheno.json that scripts/brava_phenotypes.py builds, not from the result files.
#
# PHENOTYPES selects what the run covers, in the phenocode form munge_brava.py takes:
#   CODE for the cross-ancestry meta, CODE|STRATUM for one stratum, CODE|all for both.
#   Every trait in one run lands in the SAME combined BRaVa_gene_results.mlog10p_gt4.tsv.gz,
#   so a partial re-run replaces that file rather than adding to it -- run the full phenotype
#   set whenever the combined file is being refreshed.
#
# Run inside the genetics-results-munge Docker image (polars, htslib bgzip/tabix available).
# Staging to GCS is OFF by default; pass --stage explicitly only when ready to publish.
set -euo pipefail

SCRIPTS_DIR=$(dirname "$0")
SCRIPT="$SCRIPTS_DIR/munge_brava.py"

PHENOTYPES=${PHENOTYPES:-""}
STRATA=${STRATA:-""}
WORKDIR=${WORKDIR:-"$HOME/brava_munge"}
CACHE_DIR=${CACHE_DIR:-"$WORKDIR/cache"}
OUT_DIR=${OUT_DIR:-"$WORKDIR/out"}

STAGE=()
for arg in "$@"; do
  [ "$arg" = "--stage" ] && STAGE=(--stage)
done

echo "=== BRaVa gene burden -> combined + per-trait tabix TSVs ==="
# shellcheck disable=SC2086 -- PHENOTYPES and STRATA are deliberately word-split lists
python3 "$SCRIPT" \
  ${PHENOTYPES:+--phenotypes $PHENOTYPES} \
  ${STRATA:+--strata $STRATA} \
  --cache-dir "$CACHE_DIR" \
  --output-dir "$OUT_DIR" \
  --per-trait-dir "$OUT_DIR/gene_burden_per_trait" \
  "${STAGE[@]+"${STAGE[@]}"}"

if [ ${#STAGE[@]} -eq 0 ]; then
  echo ""
  echo "=== --stage not set: skipping GCS upload (produced $OUT_DIR locally) ==="
fi
