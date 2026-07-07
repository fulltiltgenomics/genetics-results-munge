#!/bin/bash
# Munge the Zhang et al. 2021 Cell "Catlas" body-wide human snATAC atlas into the canonical
# open_chromatin TSV.
#
# Public source (verify exact file names/columns against the Catlas portal before running):
#   paper : Zhang et al., Cell 2021, 184(24):5985-6001 (doi:10.1016/j.cell.2021.10.024)
#   portal: http://catlas.org/  (adult "human tissues" resource; ~1.15M cCREs x 222 cell types)
#   mirror: https://decoder-genetics.wustl.edu/catlasv1/catlas_hub/   build: hg38
#
# Set DATA_DIR to where you downloaded the processed Catlas files, then run this inside the
# genetics-results-munge Docker image (htslib bgzip/tabix available). Staging to GCS is OFF by
# default; pass --stage explicitly (see munge_catlas.py) only when ready to publish.
set -euo pipefail

DATA_DIR=${DATA_DIR:-data/catlas}
SCRIPTS_DIR=$(dirname "$0")

# TO-VERIFY: real file names from the Catlas download once confirmed.
CCRE_MATRIX=${CCRE_MATRIX:-"$DATA_DIR/catlas_cCRE_by_celltype_cpm.tsv"}
CELLTYPE_TISSUE_MAP=${CELLTYPE_TISSUE_MAP:-"$DATA_DIR/catlas_celltype_tissue_map.tsv"}
CCRE_GENE_LINKS=${CCRE_GENE_LINKS:-"$DATA_DIR/catlas_cCRE_gene_links.tsv"}
CELL_METADATA=${CELL_METADATA:-"$DATA_DIR/catlas_celltype_metadata.tsv"}
OUTPUT=${OUTPUT:-"$DATA_DIR/catlas_open_chromatin.tsv.gz"}

echo "=== Local validation (synthetic, no download / no upload) ==="
python3 "$SCRIPTS_DIR/munge_catlas.py" --sample

echo ""
echo "=== Munging Catlas (Zhang 2021) -> open_chromatin ==="
python3 "$SCRIPTS_DIR/munge_catlas.py" \
    --ccre-matrix "$CCRE_MATRIX" \
    --celltype-tissue-map "$CELLTYPE_TISSUE_MAP" \
    --ccre-gene-links "$CCRE_GENE_LINKS" \
    --cell-metadata "$CELL_METADATA" \
    --output "$OUTPUT" \
    "$@"

echo ""
echo "=== Done (add --stage above to upload to GCS) ==="
