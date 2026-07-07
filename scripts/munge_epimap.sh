#!/bin/bash
# Munge the EpiMap (Boix et al. 2021 Nature) ChromHMM 18-state hg38 CALLS into the canonical
# open_chromatin TSV (accessible/active regulatory states only).
#
# Public source (verify build/paths before running):
#   paper   : Boix et al., Nature 2021, 590(7845):300-307 (doi:10.1038/s41586-020-03145-z)
#   CALLS   : https://personal.broadinstitute.org/cboix/epimap/ChromHMM/observed_aux_18_hg38/CALLS/
#             per-biosample files: BSS<id>_18_CALLS_segments.bed.gz  (hg38 — do NOT use hg19 CALLS)
#   metadata: https://personal.broadinstitute.org/cboix/epimap/metadata/main_metadata_table.tsv
#             (build a biosample->tissue[->life_stage] TSV, or point --meta-*-col at this table:
#              --meta-biosample-col id --meta-tissue-col SECONDARY --meta-lifestage-col lifestage)
#
# Set DATA_DIR to where you downloaded the CALLS files + metadata, then run this inside the
# genetics-results-munge Docker image (htslib bgzip/tabix available). Staging to GCS is OFF by
# default; pass --stage explicitly (see munge_epimap.py) only when ready to publish.
#
# NOTE: this does NOT download the CALLS (833 biosamples). Fetch them separately, e.g.:
#   wget -r -np -nH --cut-dirs=5 -A '*_18_CALLS_segments.bed.gz' \
#     https://personal.broadinstitute.org/cboix/epimap/ChromHMM/observed_aux_18_hg38/CALLS/
set -euo pipefail

DATA_DIR=${DATA_DIR:-data/epimap}
SCRIPTS_DIR=$(dirname "$0")

CALLS_DIR=${CALLS_DIR:-"$DATA_DIR/CALLS"}
CALLS_GLOB=${CALLS_GLOB:-"BSS*_18_CALLS_segments.bed.gz"}
BIOSAMPLE_METADATA=${BIOSAMPLE_METADATA:-"$DATA_DIR/biosample_metadata.tsv"}
LINKS=${LINKS:-""}   # optional enhancer->gene links; leave empty to skip
OUTPUT=${OUTPUT:-"$DATA_DIR/epimap_open_chromatin.tsv.gz"}

echo "=== Local validation (synthetic, no download / no upload) ==="
python3 "$SCRIPTS_DIR/munge_epimap.py" --sample

echo ""
echo "=== Munging EpiMap ChromHMM 18-state -> open_chromatin ==="
LINKS_ARG=()
if [[ -n "$LINKS" ]]; then LINKS_ARG=(--links "$LINKS"); fi

python3 "$SCRIPTS_DIR/munge_epimap.py" \
    --calls-dir "$CALLS_DIR" \
    --calls-glob "$CALLS_GLOB" \
    --biosample-metadata "$BIOSAMPLE_METADATA" \
    "${LINKS_ARG[@]}" \
    --output "$OUTPUT" \
    "$@"

echo ""
echo "=== Done (add --stage above to upload to GCS) ==="
