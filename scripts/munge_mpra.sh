#!/bin/bash
# Munge the Siraj et al. 2026 MPRA per-variant functional annotation from its WIDE source into the
# canonical LONG, tabix-ready `mpra` TSV, then (behind --stage only) publish it to both GCS buckets.
#
# Source: a bgzip'd WIDE per-variant TSV. Local default: ~/mpra_annotation.tsv.gz. Canonical upstream
#   provenance: gs://finngen-multiome-xavier/batch1_5/misc/MPRA_Siraj_etal/. Pass INPUT to override.
#
# Run inside the genetics-results-munge Docker image (htslib bgzip/tabix + polars available).
# Staging to GCS is OFF by default; pass --stage explicitly only when ready to publish.
set -euo pipefail

SCRIPTS_DIR=$(dirname "$0")
SCRIPT="$SCRIPTS_DIR/munge_mpra.py"

INPUT=${INPUT:-"$HOME/mpra_annotation.tsv.gz"}
OUTPUT=${OUTPUT:-"siraj_mpra.tsv.gz"}

# GCS destinations (used only with --stage). Both mirror results-api mpra profile config:
#   <bucket>/mpra/siraj_mpra/<file>. finngen serves from finngen-commons/results_api_data;
#   daly serves from daly-genetics-results.
GCS_FINNGEN=${GCS_FINNGEN:-"gs://finngen-commons/results_api_data/mpra/siraj_mpra/"}
GCS_DALY=${GCS_DALY:-"gs://daly-genetics-results/mpra/siraj_mpra/"}

STAGE=false
for arg in "$@"; do
  [ "$arg" = "--stage" ] && STAGE=true
done

echo "=== Local validation (synthetic, no upload) ==="
python3 "$SCRIPT" --sample

echo ""
echo "=== MPRA WIDE -> LONG (POINT index) ==="
python3 "$SCRIPT" "$INPUT" --output "$OUTPUT"

if [ "$STAGE" = true ]; then
  echo ""
  echo "=== Staging to GCS (both buckets) ==="
  for dest in "$GCS_FINNGEN" "$GCS_DALY"; do
    echo "  -> $dest"
    gcloud storage cp "$OUTPUT" "$dest"
    gcloud storage cp "$OUTPUT.tbi" "$dest"
  done
  echo "=== Staged ==="
else
  echo ""
  echo "=== --stage not set: skipping GCS upload (produced $OUTPUT + $OUTPUT.tbi locally) ==="
fi
