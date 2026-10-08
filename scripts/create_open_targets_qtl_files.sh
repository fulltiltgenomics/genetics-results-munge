#!/bin/bash

set -euxo pipefail

if [ $# -ne 2 ]; then
    echo "Usage: $0 <dataset_name> <data_dir>"
    echo "  <data_dir> holds the Open Targets release: credible_set/*.parquet and study_metadata/*.parquet"
    exit 1
fi
dataset=$1
data_dir=$2
output_file=${dataset}_credible_sets.tsv.gz
gene_metadata=gene_counts_Ensembl_105_phenotype_metadata.tsv.gz

# the gene names eQTL Catalogue uses, so that both QTL resources name a gene the same way
if [ ! -f $data_dir/$gene_metadata ]; then
    curl -L -o $data_dir/$gene_metadata "https://zenodo.org/records/7808390/files/${gene_metadata}?download=1"
fi

time python3 scripts/create_open_targets_qtl_files.py $dataset $data_dir

# a full sort rather than sort -m: the per-study files are in polars' order, which need not be
# sort's
export LC_ALL=C
per_study=($data_dir/opentargets_qtl_per_study/*.SUSIE.munged.tsv)
time cat \
<(echo -n "#") \
<(head -1 ${per_study[0]}) \
<(tail -q -n +2 "${per_study[@]}" \
| sort -t $'\t' -T $data_dir -S 4G -k6,6n -k7,7n -k8,8 -k9,9 -k3,3 -k4,4) \
| bgzip -@4 > $data_dir/$output_file \
&& tabix -f -s 6 -b 7 -e 7 $data_dir/$output_file
