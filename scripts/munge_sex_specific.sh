#!/bin/bash
# Write the FinnGen R12 sex-specific GWAS phenotype metadata JSON from the release's own
# sample-count table, then (behind --stage only) publish it to the finngen bucket.
#
# This is all the munge there is for the sex-specific and recessive releases: results-api
# serves both sets of summary statistics straight from the green library, nothing is
# copied and nothing goes to BigQuery. The recessive release covers the core endpoint set
# and reuses the core metadata file, so it has no munge at all.
#
# Source: the release `sex_difference/` directory in the FinnGen green library. Set SOURCE
#   to point at a different release. The core phenotype metadata JSON supplies the names.
#
# Staging to GCS is OFF by default; pass --stage explicitly only when ready to publish.
set -euo pipefail

SCRIPTS_DIR=$(dirname "$0")
SCRIPT="$SCRIPTS_DIR/munge_sex_specific.py"

SOURCE=${SOURCE:-"gs://finngen-production-library-green/finngen_R12/finngen_R12_analysis_data/sex_difference"}
PHENO_METADATA=${PHENO_METADATA:-"gs://finngen-commons/results_api_data/mapping_files/finngen_r12_pheno.json"}
WORKDIR=${WORKDIR:-"$HOME/sex_specific_munge"}

# used only with --stage; the JSON goes next to the other mapping files. finngen only:
# the daly profile has no registry entry for the release
GCS_MAPPING=${GCS_MAPPING:-"gs://finngen-commons/results_api_data/mapping_files"}

STAGE=false
for arg in "$@"; do
  [ "$arg" = "--stage" ] && STAGE=true
done

SRC_DIR="$WORKDIR/source"
OUT_DIR="$WORKDIR/out"
mkdir -p "$SRC_DIR"

echo "=== Fetching source ($SOURCE) ==="
gcloud storage cp "$SOURCE/R12_sex_specific_N_cases_controls_samples.tsv" "$SRC_DIR/"
gcloud storage cp "$PHENO_METADATA" "$WORKDIR/pheno_metadata.json"

echo ""
echo "=== Phenotype metadata ==="
python3 "$SCRIPT" "$SRC_DIR/R12_sex_specific_N_cases_controls_samples.tsv" \
  --pheno-metadata "$WORKDIR/pheno_metadata.json" \
  --output "$OUT_DIR"

if [ "$STAGE" = true ]; then
  echo ""
  echo "=== Staging to GCS ==="
  gcloud storage cp "$OUT_DIR/finngen_r12_sex_specific_pheno.json" "$GCS_MAPPING/"
  echo "=== Staged ==="
else
  echo ""
  echo "=== --stage not set: skipping GCS upload (produced $OUT_DIR locally) ==="
fi
