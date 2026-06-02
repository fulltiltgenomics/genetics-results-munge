#!/bin/bash
#
# Run genebass gene burden Hail export, then bgzip+tabix both full and filtered files.
#
# Usage:
#   ./convert_genebass_gene_results.sh [input_mt] [output_prefix] [pheno_blacklist]

set -euxo pipefail

INPUT="${1:-/mnt/disks/data/results.mt}"
PREFIX="${2:-gene_burden_results}"
BLACKLIST="${3:-}"

BLACKLIST_ARG=()
if [ -n "$BLACKLIST" ]; then
    BLACKLIST_ARG=(--pheno-blacklist "$BLACKLIST")
fi

time python3 scripts/genebass/convert_genebass_gene_results.py \
    --input "$INPUT" \
    --output "${PREFIX}.tsv" \
    "${BLACKLIST_ARG[@]}"

# bgzip and tabix the full file
bgzip -f "${PREFIX}.tsv"
tabix -f -s5 -b6 -e6 "${PREFIX}.tsv.gz"

# create mlog10p > 4 filtered file
# mlog10p_burden is column 9 (after removing skat/skato)
{
    zcat "${PREFIX}.tsv.gz" | head -1
    zcat "${PREFIX}.tsv.gz" | tail -n +2 | awk -F'\t' '$9 != "NA" && $9+0 > 4'
} | bgzip > "${PREFIX}.mlog10p_gt4.tsv.gz"
tabix -f -s5 -b6 -e6 "${PREFIX}.mlog10p_gt4.tsv.gz"

# split per-trait full (unfiltered) files, named by trait_original (column 16),
# bgzipped + tabixed. Rows are already sorted by gene_chr/gene_start_pos.
PER_TRAIT_DIR="${PREFIX}_per_trait"
TRAIT_COL=16
mkdir -p "$PER_TRAIT_DIR"

HEADER=$(zcat "${PREFIX}.tsv.gz" | head -1)
echo "Splitting per-trait files..."
zcat "${PREFIX}.tsv.gz" | tail -n +2 | awk -F'\t' -v dir="$PER_TRAIT_DIR" -v tc="$TRAIT_COL" \
    '{print >> (dir "/" $tc ".tsv")}'

echo "Bgzipping and tabixing per-trait files..."
for f in "$PER_TRAIT_DIR"/*.tsv; do
    { echo "$HEADER"; cat "$f"; } | bgzip -@2 > "$f.gz"
    rm "$f"
    tabix -f -s5 -b6 -e6 "$f.gz"
done

N_TRAITS=$(ls "$PER_TRAIT_DIR"/*.tsv.gz | wc -l)
echo "Done! Created ${PREFIX}.tsv.gz, ${PREFIX}.mlog10p_gt4.tsv.gz and $N_TRAITS per-trait files in ${PER_TRAIT_DIR}/"
