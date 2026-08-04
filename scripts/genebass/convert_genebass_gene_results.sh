#!/bin/bash
#
# Run the genebass gene burden Hail export, then build the two products served
# downstream:
#
#   ${PREFIX}.mlog10p_gt4.tsv.gz   combined significant hits, tabixed on the gene
#                                  locus — what /gene_based/{gene} reads
#   ${PREFIX}_per_trait/<trait>.tsv.gz
#                                  full unfiltered results for one trait, tabixed
#                                  on the gene locus — what
#                                  /gene_based_results_by_phenotype reads and what
#                                  the BigQuery gene_burden_results table is loaded from
#
# There is deliberately no combined UNFILTERED tabixed file: the unfiltered export
# is ~343M rows (76k gene x annotation rows x 4.5k phenotypes), so sorting it whole
# is a shuffle no single machine handles comfortably. Every consumer wants either
# one trait or the significant hits, and both are cheap to produce from the
# unsorted export.
#
# Usage:
#   ./convert_genebass_gene_results.sh [input_mt] [output_prefix] [pheno_blacklist] [gcs_output_dir]

set -euxo pipefail

INPUT="${1:-/mnt/disks/data/results.mt}"
PREFIX="${2:-gene_burden_results}"
BLACKLIST="${3:-}"
GCS_DIR="${4:-}"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

BLACKLIST_ARG=()
if [ -n "$BLACKLIST" ]; then
    BLACKLIST_ARG=(--pheno-blacklist "$BLACKLIST")
fi

# Hail block-gzips the export directly, so the ~65 GB plain TSV never hits disk
time python3 "${SCRIPT_DIR}/convert_genebass_gene_results.py" \
    --input "$INPUT" \
    --output "${PREFIX}.tsv.bgz" \
    "${BLACKLIST_ARG[@]}"

# combined mlog10p_burden > 4 file; small enough to sort here (mlog10p_burden is
# column 9, the gene locus is columns 5-7)
# reading only the header SIGPIPEs the decompression of a 13 GB file, which
# `pipefail` would turn into a failure — tolerate that and assert we got a header
HEADER=$(gzip -cd "${PREFIX}.tsv.bgz" 2>/dev/null | head -1) || true
[ -n "$HEADER" ] || { echo "no header in ${PREFIX}.tsv.bgz" >&2; exit 1; }

{
    printf '%s\n' "$HEADER"
    gzip -cd "${PREFIX}.tsv.bgz" | tail -n +2 \
        | awk -F'\t' '$9 != "NA" && $9+0 > 4' \
        | sort -t$'\t' -k5,5n -k6,6n -k7,7n
} | bgzip -@4 > "${PREFIX}.mlog10p_gt4.tsv.gz"
tabix -f -s5 -b6 -e6 "${PREFIX}.mlog10p_gt4.tsv.gz"

# per-trait full files, keyed by trait_original (column 16). --tmp-dir keeps the
# ~11 GB of intermediates off /tmp, next to the export instead
PER_TRAIT_DIR="${PREFIX}_per_trait"
python3 "${SCRIPT_DIR}/../split_burden_per_trait.py" \
    --input "${PREFIX}.tsv.bgz" \
    --output-dir "${PER_TRAIT_DIR}" \
    --tmp-dir "$(dirname "${PREFIX}")" \
    --trait-col 16 \
    --tabix-args "-s5 -b6 -e6"

if [ -n "$GCS_DIR" ]; then
    gcloud storage cp "${PREFIX}.mlog10p_gt4.tsv.gz" "${PREFIX}.mlog10p_gt4.tsv.gz.tbi" "${GCS_DIR%/}/"
    gcloud storage rsync "${PER_TRAIT_DIR}" "${GCS_DIR%/}/gene_burden_per_trait"
fi

N_TRAITS=$(ls "$PER_TRAIT_DIR"/*.tsv.gz | wc -l)
echo "Done! Created ${PREFIX}.mlog10p_gt4.tsv.gz and $N_TRAITS per-trait files in ${PER_TRAIT_DIR}/"
