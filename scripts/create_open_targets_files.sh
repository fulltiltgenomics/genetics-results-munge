#!/bin/bash

set -euxo pipefail

if [ $# -ne 3 ]; then
    echo "Usage: $0 <dataset_name> <data_dir> <variant_annotation_file>"
    exit 1
fi
dataset=$1
data_dir=$2
variant_annotation_file=$3
output_file=${dataset}_credible_sets.tsv.gz

mkdir -p $data_dir/opentargets_per_study

# create per-parquet file tsv files
time python3 scripts/create_open_targets_files.py $dataset $data_dir $variant_annotation_file

# create combined tsv file
# use uniq to remove duplicate header lines
time cat \
<(echo -n "#") \
<(cat $data_dir/*_cs_95.tsv \
| sort -T . -k6,6g -k7,7g -k8,8 -k9,9 | uniq) \
| bgzip -@4 > $data_dir/$output_file \
&& tabix -f -@4 -s 6 -b 7 -e 7 $data_dir/$output_file

# create per-study files
time python3 <<EOF
import polars as pl
data = pl.read_csv("$data_dir/$output_file", separator="\t", null_values=["NA"])
for trait in data["trait"].unique():
    study_data = data.filter(pl.col("trait") == trait)
    study_data.write_csv(f"$data_dir/opentargets_per_study/{trait}.SUSIE.munged.tsv", separator="\t", null_value="NA")
EOF

# create study metadata file
time python3 scripts/create_open_targets_study_file.py $data_dir/study_metadata/*.parquet $data_dir/$output_file $data_dir/study_metadata/study_metadata.json
