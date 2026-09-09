#!/bin/bash
# Munge the Collins et al. 2022 rare-CNV dosage sensitivity map (Zenodo record 6347673,
# CC-BY 4.0) into the suite's rcnv tables, then (behind --stage only) publish to both
# profile buckets.
#
# Source:
#   paper  : Collins et al., Cell 2022, 185(16):3041-3055 (doi:10.1016/j.cell.2022.06.036)
#   zenodo : https://zenodo.org/records/6347673
#   scores : Collins_rCNV_2022.dosage_sensitivity_scores.tsv.gz     (18,641 genes x pHaplo/pTriplo)
#   genes  : Collins_rCNV_2022.gene_association_sumstats.tar.gz    (108 phenotype x DEL/DUP BEDs)
#   segments: Cell supplement mmc3.xlsx, sheet "Table S3"          (163 disease-associated segments)
#
# PRODUCT selects what to munge: `scores` (dosage-sensitivity probabilities), `genes`
# (gene-based CNV association sumstats, one long table) or `segments` (the supplement's 163
# segments). scores and genes carry no coordinates at all; segments carries GRCh37 ones and is
# lifted to GRCh38 here. All three are BigQuery load files with no tabix index.
#
# Inputs are cached under CACHE_DIR and downloaded on demand (--download): the scores or
# gene-association tar from Zenodo, the gencode gene name mapping and the HGNC complete set
# from the daly mapping_files bucket (the HGNC set falls back to genenames.org if that bucket
# is not readable), and the UCSC liftOver binary and hg19ToHg38 chain for `segments`.
# mmc3.xlsx has no URL -- Elsevier and PMC answer a script with a bot-check page -- so place
# it in CACHE_DIR by hand before running PRODUCT=segments.
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
case "$PRODUCT" in
    scores) OUTPUT_SUFFIX=dosage_sensitivity ;;
    genes)  OUTPUT_SUFFIX=gene_associations ;;
    segments) OUTPUT_SUFFIX=segments ;;
    *) echo "unknown PRODUCT '$PRODUCT' (expected scores, genes or segments)" >&2; exit 1 ;;
esac
OUTPUT=${OUTPUT:-"$OUT_DIR/${DATASET_ID}_${OUTPUT_SUFFIX}.tsv.gz"}

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
