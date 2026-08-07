#!/bin/bash
# Munge the FinnGen R14 imputed classical HLA allele association results into the canonical
# per-phenotype tabix TSVs plus the combined BigQuery load file, then (behind --stage only)
# publish both to the finngen and daly buckets.
#
# Source: the release `hla/` directory in the FinnGen green library. Set SOURCE to point at a
#   different release. The phenotype metadata JSON decides which phenocodes survive: the source
#   ships ~96 extra _WIDE endpoints with no R14 metadata and those are dropped.
#
# Run inside the genetics-results-munge Docker image (htslib bgzip/tabix available).
# Staging to GCS is OFF by default; pass --stage explicitly only when ready to publish.
set -euo pipefail

SCRIPTS_DIR=$(dirname "$0")
SCRIPT="$SCRIPTS_DIR/munge_hla.py"

SOURCE=${SOURCE:-"gs://finngen-production-library-green/finngen_R14/finngen_R14_analysis_data/hla"}
PHENO_METADATA=${PHENO_METADATA:-"gs://finngen-commons/results_api_data/mapping_files/finngen_r14_pheno.json"}
DATASET_ID=${DATASET_ID:-"finngen_hla"}
WORKDIR=${WORKDIR:-"$HOME/hla_munge"}

# GCS destinations (used only with --stage). Both mirror the results-api hla profile config:
#   <bucket>/hla/<dataset-id>/{summary_stats/<PHENO>.tsv.gz, <dataset-id>.tsv.gz}
# finngen serves from finngen-commons/results_api_data; daly from daly-genetics-results.
GCS_FINNGEN=${GCS_FINNGEN:-"gs://finngen-commons/results_api_data/hla/$DATASET_ID"}
GCS_DALY=${GCS_DALY:-"gs://daly-genetics-results/hla/$DATASET_ID"}

STAGE=false
for arg in "$@"; do
  [ "$arg" = "--stage" ] && STAGE=true
done

SRC_DIR="$WORKDIR/source"
OUT_DIR="$WORKDIR/out"
mkdir -p "$SRC_DIR"

echo "=== Fetching source ($SOURCE) ==="
# the whole release hla/ directory is ~18 MiB, so one recursive copy beats 2.8k object reads
gcloud storage cp -r "$SOURCE/*" "$SRC_DIR/"
gcloud storage cp "$PHENO_METADATA" "$WORKDIR/pheno_metadata.json"

echo ""
echo "=== HLA source -> per-phenotype tabix TSVs + combined BQ file ==="
python3 "$SCRIPT" "$SRC_DIR" \
  --pheno-metadata "$WORKDIR/pheno_metadata.json" \
  --output "$OUT_DIR" \
  --dataset-id "$DATASET_ID"

if [ "$STAGE" = true ]; then
  echo ""
  echo "=== Staging to GCS (both buckets) ==="
  for dest in "$GCS_FINNGEN" "$GCS_DALY"; do
    echo "  -> $dest"
    gcloud storage cp "$OUT_DIR/$DATASET_ID.tsv.gz" "$dest/"
    gcloud storage rsync --recursive --delete-unmatched-destination-objects \
      "$OUT_DIR/summary_stats" "$dest/summary_stats"
  done
  echo "=== Staged ==="
else
  echo ""
  echo "=== --stage not set: skipping GCS upload (produced $OUT_DIR locally) ==="
fi
