#!/bin/bash
# Build the gene set bundle the sandbox's enrichment analyses read, then (behind --stage only)
# publish it to the daly bucket.
#
# Stage 1, build_gene_sets.py: download MSigDB subcollections + Reactome's own GMT into a new
# versioned gene_sets.<YYYYMMDD>/ with a manifest. Skipped when RAW_DIR names an existing one.
# Stage 2, gene_set_bundle.py: harmonise members to Ensembl via HGNC_FILE/GENCODE_FILE (the
# staged copies gene_annotations_v is built from) and write gene_sets_bundle.<release>/.
# See both .py headers and docs/gene-sets.md.
#
# Neither directory is ever overwritten. Staging is OFF by default; pass --stage only when ready
# to publish, and --dry-run with it to print the target without copying.
set -euo pipefail

SCRIPTS_DIR=$(dirname "$0")
OUT_DIR=${OUT_DIR:-"$HOME/gene_sets_build"}
RAW_DIR=${RAW_DIR:-}
HGNC_FILE=${HGNC_FILE:?set HGNC_FILE (gs://daly-genetics-results/mapping_files/hgnc_complete_set.txt, copied locally)}
HGNC_VERSION=${HGNC_VERSION:?set HGNC_VERSION to the hgnc_version gene_annotations_v carries, e.g. 2026-06-02}
GENCODE_FILE=${GENCODE_FILE:?set GENCODE_FILE (gs://daly-genetics-results/mapping_files/gencode.v49.annotation.genes.tsv, copied locally)}
GENCODE_VERSION=${GENCODE_VERSION:-49}

STAGE_ARGS=()
STAGE=0 DRY_RUN=0
for arg in "$@"; do
  case "$arg" in
    --stage) STAGE=1; STAGE_ARGS+=("$arg") ;;
    --dry-run) DRY_RUN=1; STAGE_ARGS+=("$arg") ;;
    *) echo "unknown argument $arg" >&2; exit 1 ;;
  esac
done
# a bare --dry-run would otherwise run the full download and build for real
if [ "$DRY_RUN" = 1 ] && [ "$STAGE" = 0 ]; then
  echo "--dry-run only applies to --stage; pass --stage --dry-run" >&2
  exit 1
fi

if [ -z "$RAW_DIR" ]; then
  # the script names the dir after the build date, so the path is read back rather than rebuilt
  RAW_DIR=$(python3 "$SCRIPTS_DIR/build_gene_sets.py" --out-dir "$OUT_DIR" | sed -n 's/^wrote //p')
  [ -n "$RAW_DIR" ] && [ -d "$RAW_DIR" ] || { echo "download produced no directory" >&2; exit 1; }
fi

python3 "$SCRIPTS_DIR/gene_set_bundle.py" --raw-dir "$RAW_DIR" --out-dir "$OUT_DIR" \
  --hgnc "$HGNC_FILE" --hgnc-version "$HGNC_VERSION" \
  --gencode "$GENCODE_FILE" --gencode-version "$GENCODE_VERSION" "${STAGE_ARGS[@]}"
