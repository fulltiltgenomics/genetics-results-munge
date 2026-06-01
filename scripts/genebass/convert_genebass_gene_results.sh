#!/bin/bash
#
# Run genebass gene burden Hail export, then bgzip+tabix both full and filtered files.
#
# Usage:
#   ./convert_genebass_gene_results.sh [input_mt] [output_prefix]

set -euxo pipefail

INPUT="${1:-/mnt/disks/data/results.mt}"
PREFIX="${2:-gene_burden_results}"

time python3 scripts/genebass/convert_genebass_gene_results.py \
    --input "$INPUT" \
    --output "${PREFIX}.tsv"

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

echo "Done! Created ${PREFIX}.tsv.gz and ${PREFIX}.mlog10p_gt4.tsv.gz"
