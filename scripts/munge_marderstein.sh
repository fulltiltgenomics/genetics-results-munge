#!/bin/bash
# Munge the Marderstein/Kundaje 2026 neuro-variants/FLARE release into the canonical
# open_chromatin (Product A) and variant_effect (Product B) TSVs — THREE output files.
#
# Sources (verify exact file names / column layouts against the real downloads before running):
#   Synapse project syn64693551 "neuro_variants" (+ second release entity syn73770440); build hg38
#     peaks       -> syn64716764   => open_chromatin  (marderstein_open_chromatin)
#     predictions -> syn64713923   => variant_effect  (marderstein_chrombpnet)
#     flare       -> syn64717038   => variant_effect  (marderstein_flare)
#   ENCODE adult-heart snATAC (public, separate accession) -> extra open_chromatin contexts.
#   paper: Marderstein/Kundu/Kundaje/Montgomery 2026, Nature Genetics.
#
# Data access: synapseclient reads SYNAPSE_AUTH_TOKEN from the environment. Downloading is OFF by
# default; this wrapper assumes the input files are ALREADY downloaded under DATA_DIR. To fetch,
# add --download (per product) with a fresh SYNAPSE_AUTH_TOKEN in the env (never commit the token).
#
# Run inside the genetics-results-munge Docker image (htslib bgzip/tabix + synapseclient available).
# Staging to GCS is OFF by default; pass --stage explicitly only when ready to publish.
set -euo pipefail

DATA_DIR=${DATA_DIR:-data/marderstein}
SCRIPTS_DIR=$(dirname "$0")
SCRIPT="$SCRIPTS_DIR/munge_marderstein.py"

# TO-VERIFY: real file names + column layouts once confirmed against the Synapse/ENCODE downloads.
PEAKS=${PEAKS:-"$DATA_DIR/peaks_long.tsv"}                 # long: chrom,start,end,cell_type,score
CONTEXT_MAP=${CONTEXT_MAP:-"$DATA_DIR/context_map.tsv"}    # cell_type -> tissue,life_stage[,assay] overrides
# ChromBPNet: the REAL release is a WIDE per-variant matrix (~439 cols; 132 context triplets), streamed
# in bounded batches. Point CHROMBPNET at the downloaded wide file(s), or pass --download to fetch +
# stream + delete each file one at a time (peak disk = largest single input + the filtered output).
CHROMBPNET=${CHROMBPNET:-"$DATA_DIR/asd.all_dataset.K562_bias.annot2.txt.gz"}
FLARE=${FLARE:-"$DATA_DIR/flare_scores.tsv"}

OUT_OPEN=${OUT_OPEN:-"$DATA_DIR/marderstein_open_chromatin.tsv.gz"}
OUT_CHROMBPNET=${OUT_CHROMBPNET:-"$DATA_DIR/marderstein_chrombpnet.tsv.gz"}
OUT_FLARE=${OUT_FLARE:-"$DATA_DIR/marderstein_flare.tsv.gz"}

echo "=== Local validation (synthetic, no download / no upload) ==="
python3 "$SCRIPT" --sample

echo ""
echo "=== [1/3] Marderstein peaks -> open_chromatin (INTERVAL index) ==="
python3 "$SCRIPT" --product open_chromatin \
    --peaks "$PEAKS" \
    --context-map "$CONTEXT_MAP" \
    --output "$OUT_OPEN" \
    "$@"

echo ""
echo "=== [2/3] ChromBPNet -> variant_effect (WIDE streaming reshape; thresholded; POINT index) ==="
# tissue/life_stage come from the study suffix of each "<celltype>.<study>" context (no --context-map).
# For the full common/rare/ldsc run, replace --chrombpnet with: --download [--cbp-name-filter <regex>]
python3 "$SCRIPT" --product chrombpnet \
    --chrombpnet "$CHROMBPNET" \
    --output "$OUT_CHROMBPNET" \
    "$@"

echo ""
echo "=== [3/3] FLARE -> variant_effect (pan-context; POINT index) ==="
python3 "$SCRIPT" --product flare \
    --flare "$FLARE" \
    --output "$OUT_FLARE" \
    "$@"

echo ""
echo "=== Done (add --stage above to upload to both GCS buckets) ==="
