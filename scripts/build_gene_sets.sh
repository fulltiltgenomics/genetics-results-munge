#!/bin/bash
# Download the gene set collections (MSigDB subcollections + Reactome's own GMT) into a new
# versioned directory with a manifest. Download stage only: no harmonisation, no bundle, and no
# GCS staging -- that arrives with the bundle stage. See build_gene_sets.py and docs/gene-sets.md.
#
# The output directory is named after the UTC build date and is never overwritten; a second run
# on the same day refuses.
set -euo pipefail

SCRIPTS_DIR=$(dirname "$0")
OUT_DIR=${OUT_DIR:-"$HOME/gene_sets_build"}

python3 "$SCRIPTS_DIR/build_gene_sets.py" --out-dir "$OUT_DIR" "$@"
