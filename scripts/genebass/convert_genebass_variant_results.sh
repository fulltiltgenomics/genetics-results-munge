#!/bin/bash
#
# Post-process genebass variant results exported by the Hail dataproc script.
# Splits per-trait files (unfiltered, bgzipped+tabixed) and creates
# a combined mlog10p>4 filtered file.
#
# Usage:
#   ./convert_genebass_variant_results.sh <input_tsv_bgz> <output_dir>
#
# Example:
#   ./convert_genebass_variant_results.sh \
#     gs://finngen-commons/genebass/variant_results.tsv.bgz \
#     gs://finngen-commons/genebass/variant_results_v2

set -euo pipefail

if [ $# -lt 2 ]; then
    echo "Usage: $0 <input_tsv_bgz> <output_dir>"
    exit 1
fi

INPUT="$1"
OUTPUT_DIR="$2"

TMPDIR=$(mktemp -d)
trap 'rm -rf "$TMPDIR"' EXIT

LOCAL_INPUT="$TMPDIR/variant_results.tsv.bgz"
PER_TRAIT_DIR="$TMPDIR/per_trait"
COMBINED="$TMPDIR/genebass_variant_results_mlog10p4.tsv"

mkdir -p "$PER_TRAIT_DIR"

# download input if on GCS
if [[ "$INPUT" == gs://* ]]; then
    echo "Downloading $INPUT..."
    gcloud storage cp "$INPUT" "$LOCAL_INPUT"
else
    LOCAL_INPUT="$INPUT"
fi

HEADER=$(zcat "$LOCAL_INPUT" | head -1)
# columns: #dataset(1) chr(2) pos(3) ref(4) alt(5) gene(6) annotation(7)
#          mlog10p(8) beta(9) se(10) af_overall(11) af_cases(12) af_controls(13)
#          ac(14) an(15) n_cases(16) n_controls(17) trait(18) trait_original(19)
# split per-trait by trait_original (unique key)
TRAIT_COL=19
MLOG10P_COL=8

echo "Splitting per-trait files..."
zcat "$LOCAL_INPUT" | tail -n +2 | awk -F'\t' -v dir="$PER_TRAIT_DIR" -v tc="$TRAIT_COL" \
    '{print >> (dir "/" $tc ".tsv")}'

echo "Bgzipping and tabixing per-trait files..."
for f in "$PER_TRAIT_DIR"/*.tsv; do
    trait=$(basename "$f" .tsv)
    # prepend header
    { echo "$HEADER"; cat "$f"; } | bgzip -@2 > "$f.gz"
    rm "$f"
    tabix -f -s2 -b3 -e3 "$f.gz"
done

echo "Creating combined mlog10p>4 filtered file..."
{
    echo "$HEADER"
    zcat "$LOCAL_INPUT" | tail -n +2 | awk -F'\t' -v col="$MLOG10P_COL" \
        '$col != "NA" && $col+0 > 4'
} | bgzip -@4 > "$COMBINED.gz"
tabix -f -s2 -b3 -e3 "$COMBINED.gz"

# upload results
if [[ "$OUTPUT_DIR" == gs://* ]]; then
    echo "Uploading per-trait files to $OUTPUT_DIR/per_trait/..."
    gcloud storage cp "$PER_TRAIT_DIR"/*.tsv.gz "$PER_TRAIT_DIR"/*.tsv.gz.tbi "$OUTPUT_DIR/per_trait/"

    echo "Uploading combined filtered file to $OUTPUT_DIR/..."
    gcloud storage cp "$COMBINED.gz" "$COMBINED.gz.tbi" "$OUTPUT_DIR/"
else
    mkdir -p "$OUTPUT_DIR/per_trait"
    cp "$PER_TRAIT_DIR"/*.tsv.gz "$PER_TRAIT_DIR"/*.tsv.gz.tbi "$OUTPUT_DIR/per_trait/"
    cp "$COMBINED.gz" "$COMBINED.gz.tbi" "$OUTPUT_DIR/"
fi

N_TRAITS=$(ls "$PER_TRAIT_DIR"/*.tsv.gz | wc -l)
echo "Done! Created $N_TRAITS per-trait files and 1 combined filtered file."
