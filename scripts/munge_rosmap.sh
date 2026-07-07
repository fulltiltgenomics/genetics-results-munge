#!/bin/bash
# Munge the Xiong et al. 2023 Cell ROSMAP AD-brain snATAC epigenome into the canonical
# open_chromatin TSV.
#
# Public source (verify exact file names/columns before running):
#   paper : Xiong et al., Cell 2023, 186(20):4422-4437 (doi:10.1016/j.cell.2023.08.040)
#   data  : PUBLIC processed per-cell-type / per-subtype accessible-peak BEDs (+ peak->gene links)
#           from the Broad `ad_epigenome` resource. hg38. NO controlled-access / AD Knowledge
#           Portal DUC data is used (raw fragments / genotypes are NOT needed for Product A).
#   build : hg38 -> NO liftOver.
#
# TO-VERIFY: whether the public release splits peaks by AD/control (per-file or a column) or ships
# UNION peaks per cell type. Set --condition-from-filename / --condition-idx accordingly; union
# peaks default to condition=unknown.
#
# Run inside the genetics-results-munge Docker image (htslib bgzip/tabix + polars available).
# Staging to GCS is OFF by default; pass --stage only when ready to publish.
set -euo pipefail

DATA_DIR=${DATA_DIR:-data/rosmap}
SCRIPTS_DIR=$(dirname "$0")

# TO-VERIFY: real paths / naming once confirmed against the ad_epigenome download.
PEAKS_DIR=${PEAKS_DIR:-"$DATA_DIR/peaks"}
LINKS=${LINKS:-"$DATA_DIR/peak_gene_links.tsv"}
OUTPUT=${OUTPUT:-"$DATA_DIR/rosmap_open_chromatin.tsv.gz"}

echo "=== Local validation (synthetic hg38 peaks, no download / no upload) ==="
python3 "$SCRIPTS_DIR/munge_rosmap.py" --sample

echo ""
echo "=== Munging ROSMAP (Xiong 2023) -> open_chromatin (hg38, no liftOver) ==="
# adjust --celltype-from-filename / --condition-from-filename to the real ad_epigenome naming;
# add --links "$LINKS" when the peak->gene links file is available;
# add --score-idx N (+ optional --score-type) if the BEDs carry a signal column.
python3 "$SCRIPTS_DIR/munge_rosmap.py" \
    --peaks-dir "$PEAKS_DIR" \
    --celltype-from-filename '^([^.]+)\.' \
    --output "$OUTPUT" \
    "$@"

echo ""
echo "=== Done (add --stage above to upload to GCS) ==="
