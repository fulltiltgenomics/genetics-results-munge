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

# create per-study files with stats
time python3 <<EOF
import sys
sys.path.insert(0, "scripts")
import polars as pl
from credible_set_stats import calculate_stats, write_stats_json, get_tsv_header, stats_to_tsv_row

data = pl.read_csv("$data_dir/$output_file", separator="\t", null_values=["NA"])
all_stats = []
for trait in data["trait"].unique():
    study_data = data.filter(pl.col("trait") == trait)
    study_data.write_csv(f"$data_dir/opentargets_per_study/{trait}.SUSIE.munged.tsv", separator="\t", null_value="NA")
    stats = calculate_stats(study_data)
    write_stats_json(stats, f"$data_dir/opentargets_per_study/{trait}.SUSIE.munged.stats.json")
    all_stats.append(stats)

with open(f"$data_dir/opentargets_per_study/credible_set_stats.tsv", "w") as f:
    f.write(get_tsv_header() + "\n")
    for s in all_stats:
        f.write(stats_to_tsv_row(s) + "\n")
print(f"Wrote aggregate stats for {len(all_stats)} traits")
EOF

# create study metadata file
time python3 scripts/create_open_targets_study_file.py $data_dir/study_metadata/*.parquet $data_dir/$output_file $data_dir/study_metadata/study_metadata.json
