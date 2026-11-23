#!/bin/bash
# Build the gene-indexed copy of a munged coloc file (<name>.gene.tsv.gz + .tbi) that backs
# the API's /colocalization_by_gene endpoint. See create_gene_indexed_coloc_file.py for what
# the copy holds and how a QTL trait is resolved to a gene.
#
# Needs bgzip/tabix and polars. The gencode file must be the version the API's coloc profile
# entry names in gene_gencode_version. Nothing is uploaded.
#
# Usage: scripts/coloc/create_gene_indexed_coloc_file.sh <munged_coloc.tsv.gz> \
#            <gencode.vNN.annotation.genes.tsv> [trait_gene_map.tsv]
# Output: <munged_coloc>.gene.tsv.gz next to the input, with its .tbi and an
#         .unmapped.tsv of the rows no gene was found for.

set -euo pipefail

if [ $# -lt 2 ]; then
    echo "Usage: $0 <munged_coloc.tsv.gz> <gencode.genes.tsv> [trait_gene_map.tsv]"
    exit 1
fi
coloc=$1
gencode=$2
trait_gene_map=${3:-}

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
out_tsv=${coloc%.tsv.gz}.gene.tsv
out_gz=$out_tsv.gz

PYTHONPATH="$SCRIPT_DIR/..${PYTHONPATH:+:$PYTHONPATH}" python3 \
    "$SCRIPT_DIR/create_gene_indexed_coloc_file.py" \
    --coloc "$coloc" \
    --gene-mapping "$gencode" \
    ${trait_gene_map:+--trait-gene-map "$trait_gene_map"} \
    --output "$out_tsv"

col() { head -1 "$out_tsv" | tr '\t' '\n' | grep -nx "$1" | cut -d: -f1; }
gc=$(col gene_chr); gs=$(col gene_start); ge=$(col gene_end)
c=$(col chr); rs=$(col region_start_min); re=$(col region_end_max)
t1=$(col trait1); t2=$(col trait2)

# tabix needs the rows ordered by the gene locus; the keys after it are the order the API
# merges coloc files in. LC_ALL=C sorts the string keys byte-wise, as the API compares them
(
    head -1 "$out_tsv"
    tail -n +2 "$out_tsv" \
    | LC_ALL=C sort -t $'\t' -T "$(dirname "$out_tsv")" -S 2G \
        -k$gc,${gc}n -k$gs,${gs}n -k$ge,${ge}n -k$c,${c}n -k$rs,${rs}n -k$re,${re}n -k$t1,$t1 -k$t2,$t2
) | bgzip -@4 > "$out_gz"
tabix -f -s "$gc" -b "$gs" -e "$ge" "$out_gz"
rm "$out_tsv"
echo "wrote $out_gz"
