#!/bin/bash
# Build a protein structure release from the lab's reference bundle and, with --upload, stage it
# under the rehearsal root (see munge_protein_structure.py's header for inputs and decisions,
# docs/protein-structure.md for the contract). The build holds one chromosome or one protein's
# rows at a time but the residue stage still peaks at a few GB, and an OOM kill on a shared VM
# takes the whole session with it, so it runs in its own memory-capped scope.
#
# Extra arguments pass through: `--stages variants,manifest` resumes, `--upload` stages.
set -euo pipefail

SCRIPTS_DIR=$(cd "$(dirname "$0")" && pwd)
BASE=${BASE:-/mnt/disks/data/protein-structure}
PYTHON=${PYTHON:-$BASE/venv/bin/python}
WORKDIR=${WORKDIR:-$BASE/release-build}
TED_DIR=${TED_DIR:-$BASE/ted}
MEMORY_MAX=${MEMORY_MAX:-5G}

TED_ARGS=()
if [ -f "$TED_DIR/ted_human.domain_summary.tsv" ]; then
  TED_ARGS=(--ted "$TED_DIR/ted_human.domain_summary.tsv" --ted-meta "$TED_DIR/ted_meta.json")
fi

# -I: the bundle is untrusted data, so nothing is imported from the working directory
exec systemd-run --user --scope -p MemoryMax="$MEMORY_MAX" \
  "$PYTHON" -I "$SCRIPTS_DIR/munge_protein_structure.py" \
  --bundle "$BASE/bundle/sir-reference-data" --bundle-tar "$BASE/bundle/sir-data.tar.gz" \
  --lab-files "$BASE/lab-files" \
  --gencode "${GENCODE:-/mnt/disks/data/cellxgene-pg/ref/gencode.v49.annotation.genes.tsv}" \
  --workdir "$WORKDIR" "${TED_ARGS[@]}" "$@"
