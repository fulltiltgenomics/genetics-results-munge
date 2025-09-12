#!/bin/bash

set -euxo pipefail

if [ $# -ne 3 ]; then
    echo "Usage: $0 <dataset_name> <data_dir> <variant_annotation_file>"
    exit 1
fi
dataset=$1
data_dir=$2
variant_annotation_file=$3
output_file=$dataset.tsv.gz

mkdir -p $data_dir/eqtlcat_per_study

time python3 scripts/create_eqtl_catalogue_files.py $data_dir $variant_annotation_file

time (
    ls -1 $data_dir/eqtlcat_per_study/QTD*.SUSIE.munged.tsv \
    | tr '\n' '\0' > $data_dir/merge_these \
    && cat \
    <(echo -n "#") \
    <(sort -m -T . --files0-from=$data_dir/merge_these -k6,6g -k7,7g -k8,8 -k9,9 -k3,3 | uniq) \
    | bgzip -@4 > $data_dir/$output_file \
    && tabix -@4 -s 6 -b 7 -e 7 $data_dir/$output_file
)
