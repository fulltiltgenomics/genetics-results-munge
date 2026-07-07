#!/bin/bash
# Munge the Li et al. 2023 human brain snATAC atlas into the canonical open_chromatin TSV.
#
# Public source (verify exact file names/columns against the CATlas portal before running):
#   paper : Li et al., Science 2023, 382(6667):eadf7044 (doi:10.1126/science.adf7044)
#   portal: https://catlas.org/humanbrain/   GEO: GSE244618   build: hg38
#
# Set DATA_DIR to where you downloaded the processed cCRE files, then run this inside the
# genetics-results-munge Docker image (htslib bgzip/tabix available). Staging to GCS is OFF
# by default; pass --stage explicitly (see munge_li_brain.py) only when ready to publish.
set -euo pipefail

DATA_DIR=${DATA_DIR:-data/li_brain}
SCRIPTS_DIR=$(dirname "$0")

# TO-VERIFY: real file names from the CATlas download once confirmed.
CCRE_MATRIX=${CCRE_MATRIX:-"$DATA_DIR/human_brain_cCRE_by_celltype_cpm.tsv"}
CCRE_GENE_LINKS=${CCRE_GENE_LINKS:-"$DATA_DIR/human_brain_cCRE_gene_links.tsv"}
CELL_METADATA=${CELL_METADATA:-"$DATA_DIR/human_brain_celltype_metadata.tsv"}
OUTPUT=${OUTPUT:-"$DATA_DIR/li_brain_open_chromatin.tsv.gz"}

echo "=== Local validation (synthetic, no download / no upload) ==="
python3 "$SCRIPTS_DIR/munge_li_brain.py" --sample

echo ""
echo "=== Munging Li 2023 brain -> open_chromatin ==="
python3 "$SCRIPTS_DIR/munge_li_brain.py" \
    --ccre-matrix "$CCRE_MATRIX" \
    --ccre-gene-links "$CCRE_GENE_LINKS" \
    --cell-metadata "$CELL_METADATA" \
    --output "$OUTPUT" \
    "$@"

echo ""
echo "=== Done (add --stage above to upload to GCS) ==="
