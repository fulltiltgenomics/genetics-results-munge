#!/bin/bash
# Download and munge a gene-disease association source into the versioned TSV the results-api's
# /gene_disease/{gene} endpoint reads, then (behind --stage only) publish it to the profile buckets.
#
# Source:
#   gencc   https://search.thegencc.org/download/action/submissions-export-tsv
#   monarch https://data.monarchinitiative.org/monarch-kg/latest/  (two TSVs plus metadata.yaml)
#
# PRODUCT selects the source: `gencc` or `monarch`. Neither output carries coordinates, so neither
# is bgzipped or tabixed -- the API reads the plain TSV. See the .py header for what is dropped,
# deduplicated and renamed, and docs/gene-disease-associations.md for why.
#
# The output file name carries the source's version (the KG release for monarch, the UTC download
# date for gencc, which publishes no release identifier), so the previous file stays in place
# rather than being overwritten. Point the results-api profile at the new name to cut over.
#
# Run inside the genetics-results-munge Docker image (polars available).
# Staging to GCS is OFF by default; pass --stage explicitly only when ready to publish.
set -euo pipefail

SCRIPTS_DIR=$(dirname "$0")
SCRIPT="$SCRIPTS_DIR/munge_gene_disease.py"

PRODUCT=${PRODUCT:-gencc}
CACHE_DIR=${CACHE_DIR:-"$HOME/gene_disease_munge/cache"}
OUT_DIR=${OUT_DIR:-"$HOME/gene_disease_munge/out"}

# GCS destinations (used only with --stage): <bucket>/gene_disease/. finngen serves from
# finngen-commons/results_api_data; daly from daly-genetics-results. Override GCS_DESTS to stage to
# one bucket only -- a refresh that is not being taken up by both deployments should not write both.
GCS_DESTS=${GCS_DESTS:-"gs://finngen-commons/results_api_data/gene_disease/ gs://daly-genetics-results/gene_disease/"}

STAGE=false
for arg in "$@"; do
  [ "$arg" = "--stage" ] && STAGE=true
done

echo "=== Munging gene-disease --product $PRODUCT ==="
# the script names the output after the source version, so the path is read back rather than
# rebuilt here -- two spellings of the same version string is how they drift apart
OUTPUT=$(python3 "$SCRIPT" --product "$PRODUCT" --download --cache-dir "$CACHE_DIR" --out-dir "$OUT_DIR" \
  | sed -n 's/^wrote //p')
[ -n "$OUTPUT" ] && [ -f "$OUTPUT" ] || { echo "munge produced no output file" >&2; exit 1; }
echo "  $OUTPUT"

if [ "$STAGE" = true ]; then
  echo ""
  echo "=== Staging to GCS ==="
  for dest in $GCS_DESTS; do
    echo "  -> $dest"
    gcloud storage cp "$OUTPUT" "$dest"
  done
  echo "=== Staged ==="
  echo "Now point the results-api gene_disease profile at $(basename "$OUTPUT")."
else
  echo ""
  echo "=== --stage not set: skipping GCS upload (produced $OUTPUT locally) ==="
fi
