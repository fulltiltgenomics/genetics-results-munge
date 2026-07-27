#!/bin/bash
set -euo pipefail

DATA_DIR=/mnt/disks/data/CPGMethylation/Publish
REF=/mnt/disks/data/ref/GRCh38.fa
GNOMAD=/mnt/disks/data/gnomad/gnomad.genomes.exomes.v4.0.sites.v2.tsv.bgz
SCRIPTS_DIR=$(dirname "$0")

echo "=== Normalizing Data-S1 (CpG QTLs) ==="
python "$SCRIPTS_DIR/normalize_asmqtl.py" --input "$DATA_DIR/Data-S1.tab" --ref "$REF"

echo ""
echo "=== Normalizing Data-S3 (MDS QTLs) ==="
python "$SCRIPTS_DIR/normalize_asmqtl.py" --input "$DATA_DIR/Data-S3.tab" --ref "$REF"

echo ""
echo "=== Munging Data-S1 ==="
python "$SCRIPTS_DIR/munge_asmqtl.py" --input "$DATA_DIR/Data-S1.tab" \
    --norm-mapping "$DATA_DIR/Data-S1.norm_mapping.tsv" \
    --gnomad "$GNOMAD" \
    --gnomad-af-plot

echo ""
echo "=== Munging Data-S3 (reuses gnomAD from S1) ==="
python "$SCRIPTS_DIR/munge_asmqtl.py" --input "$DATA_DIR/Data-S3.tab" \
    --norm-mapping "$DATA_DIR/Data-S3.norm_mapping.tsv" \
    --gnomad-af-plot \
    --gnomad-filtered "$DATA_DIR/Data-S1.gnomad_filtered.tsv"

echo ""
echo "=== All done ==="
