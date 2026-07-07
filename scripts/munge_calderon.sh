#!/bin/bash
# Munge the Calderon et al. 2019 Nat Genet stimulation-responsive immune ATAC-seq dataset into the
# canonical open_chromatin TSV.
#
# Public source (verify exact file names/columns/build before running):
#   paper : Calderon et al., Nature Genetics 2019, 51:1494-1505 (doi:10.1038/s41588-019-0505-9)
#   GEO   : GSE118189 -> GSE118189_ATAC_counts.txt.gz (peak x sample count matrix)
#   supp  : Supplementary_data_3_ATAC_stimulation_DA_peaks (stimulated-vs-resting DA peaks, optional)
#   build : hg19  -> MUST liftOver to hg38 (the distinctive step for this dataset)
#
# Sample columns encode "<donor>-<cell_type>-<U|S>": -U -> condition=resting, -S -> stimulated.
#
# liftOver requires the UCSC `liftOver` binary and the hg19ToHg38 chain:
#   https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz
# Set DATA_DIR to where you downloaded the inputs + chain, then run inside the genetics-results-munge
# Docker image (htslib bgzip/tabix available; ensure liftOver is on PATH). Staging to GCS is OFF by
# default; pass --stage explicitly (see munge_calderon.py) only when ready to publish.
set -euo pipefail

DATA_DIR=${DATA_DIR:-data/calderon}
SCRIPTS_DIR=$(dirname "$0")

# TO-VERIFY: real file names once confirmed against the GEO download / paper supplement.
COUNTS=${COUNTS:-"$DATA_DIR/GSE118189_ATAC_counts.txt.gz"}
DA_PEAKS=${DA_PEAKS:-"$DATA_DIR/Supplementary_data_3_ATAC_stimulation_DA_peaks.tsv"}
CHAIN=${CHAIN:-"$DATA_DIR/hg19ToHg38.over.chain.gz"}
OUTPUT=${OUTPUT:-"$DATA_DIR/calderon_open_chromatin.tsv.gz"}

echo "=== Local validation (synthetic hg38, --skip-liftover, no download / no upload) ==="
python3 "$SCRIPTS_DIR/munge_calderon.py" --sample

echo ""
echo "=== Munging Calderon (2019) -> open_chromatin (hg19 -> liftOver -> hg38) ==="
# add --da-peaks "$DA_PEAKS" to append the stimulation DA supplement (score_type=stim_log2fc)
python3 "$SCRIPTS_DIR/munge_calderon.py" \
    --counts "$COUNTS" \
    --chain "$CHAIN" \
    --output "$OUTPUT" \
    "$@"

echo ""
echo "=== Done (add --stage above to upload to GCS) ==="
