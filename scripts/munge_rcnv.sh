#!/bin/bash
# Munge the Collins et al. 2022 rare-CNV dosage sensitivity map (Zenodo record 6347673,
# CC-BY 4.0) into the suite's rcnv tables, then (behind --stage only) publish to both
# profile buckets.
#
# Source:
#   paper  : Collins et al., Cell 2022, 185(16):3041-3055 (doi:10.1016/j.cell.2022.06.036)
#   zenodo : https://zenodo.org/records/6347673
#   scores : Collins_rCNV_2022.dosage_sensitivity_scores.tsv.gz  (18,641 genes x pHaplo/pTriplo)
#
# --product selects which Zenodo product to munge; `scores` is the only one implemented.
# The scores carry no coordinates, so the output is build-independent and needs no liftOver
# and no tabix index -- it is a BigQuery load file.
#
# Inputs are cached under CACHE_DIR and downloaded on demand (--download): the scores from
# Zenodo, the gencode gene name mapping and the HGNC complete set from the daly
# mapping_files bucket (the HGNC set falls back to genenames.org if that bucket is not
# readable).
#
# Run inside the genetics-results-munge Docker image (htslib bgzip available).
# Staging to GCS is OFF by default; pass --stage explicitly only when ready to publish.
set -euo pipefail

SCRIPTS_DIR=$(dirname "$0")
SCRIPT="$SCRIPTS_DIR/munge_rcnv.py"

PRODUCT=${PRODUCT:-scores}
CACHE_DIR=${CACHE_DIR:-"$HOME/rcnv_munge/cache"}
OUT_DIR=${OUT_DIR:-"$HOME/rcnv_munge/out"}
DATASET_ID=${DATASET_ID:-"collins_rcnv_2022"}
OUTPUT=${OUTPUT:-"$OUT_DIR/${DATASET_ID}_dosage_sensitivity.tsv.gz"}

# GCS destinations (used only with --stage): <bucket>/rcnv/<dataset-id>/
# finngen serves from finngen-commons/results_api_data; daly from daly-genetics-results.
GCS_FINNGEN=${GCS_FINNGEN:-"gs://finngen-commons/results_api_data/rcnv/$DATASET_ID/$(basename "$OUTPUT")"}
GCS_DALY=${GCS_DALY:-"gs://daly-genetics-results/rcnv/$DATASET_ID/$(basename "$OUTPUT")"}

echo "=== Munging Collins rCNV 2022 --product $PRODUCT ==="
python3 "$SCRIPT" \
    --product "$PRODUCT" \
    --download \
    --cache-dir "$CACHE_DIR" \
    --output "$OUTPUT" \
    --gcs-finngen "$GCS_FINNGEN" \
    --gcs-daly "$GCS_DALY" \
    "$@"

echo ""
echo "=== Done (add --stage above to upload to GCS) ==="
