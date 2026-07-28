#!/bin/bash

set -euxo pipefail

if [ $# -ne 3 ]; then
    echo "Usage: $0 <dataset_name> <data_dir> <variant_annotation_file>"
    echo "  <data_dir> holds the Open Targets release: credible_set/*.parquet and study_metadata/*.parquet"
    exit 1
fi
dataset=$1
data_dir=$2
variant_annotation_file=$3
output_file=${dataset}_credible_sets.tsv.gz
# Open_Targets_26.06 -> ot_2606_data_studies.json
version=${dataset##*_}
study_json=ot_${version//./}_data_studies.json

mkdir -p $data_dir/opentargets_per_study

# convert the release to one unsorted TSV
time python3 scripts/create_open_targets_files.py $dataset $data_dir $variant_annotation_file

# sort on chr and pos for tabix, keeping the header first
time cat \
<(echo -n "#") \
<(head -1 $data_dir/${dataset}_cs_95.tsv) \
<(tail -n +2 $data_dir/${dataset}_cs_95.tsv \
| sort -T $data_dir -S 4G -k6,6n -k7,7n -k8,8 -k9,9) \
| bgzip -@4 > $data_dir/$output_file \
&& tabix -f -s 6 -b 7 -e 7 $data_dir/$output_file

# create per-study files with stats
time python3 <<EOF
import sys
sys.path.insert(0, "scripts")
import polars as pl
from credible_set_stats import calculate_stats, write_stats_json, get_tsv_header, stats_to_tsv_row

data = pl.read_csv("$data_dir/$output_file", separator="\t", null_values=["NA"])
all_stats = []
for (trait,), study_data in data.partition_by("trait", as_dict=True).items():
    study_data.write_csv(f"$data_dir/opentargets_per_study/{trait}.SUSIE.munged.tsv", separator="\t", null_value="NA")
    stats = calculate_stats(study_data)
    write_stats_json(stats, f"$data_dir/opentargets_per_study/{trait}.SUSIE.munged.stats.json")
    all_stats.append(stats)

# kept out of opentargets_per_study/ so that directory holds only the per-study files
with open(f"$data_dir/credible_set_stats.tsv", "w") as f:
    f.write(get_tsv_header() + "\n")
    for s in all_stats:
        f.write(stats_to_tsv_row(s) + "\n")
print(f"Wrote aggregate stats for {len(all_stats)} traits")
EOF

# create study metadata file
time python3 scripts/create_open_targets_study_file.py $data_dir/study_metadata/*.parquet $data_dir/$output_file $data_dir/$study_json
