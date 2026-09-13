#!/bin/bash

set -euxo pipefail

if [ $# -ne 4 ]; then
    echo "Usage: $0 <finemapping_parquet> <lead_variants> <variant_annotation_file> <data_dir>"
    echo "  <finemapping_parquet>     UKBB_EUR_fine_mapping_with_meta_EUR_lead_variants.parquet"
    echo "                            from https://zenodo.org/records/18132538"
    echo "  <lead_variants>           EUR_all_lead_variants.tsv.gz"
    echo "                            from https://zenodo.org/records/18377015"
    echo "  <variant_annotation_file> tabix indexed FinnGen annotated variants"
    exit 1
fi
parquet=$1
lead_variants=$2
variant_annotation_file=$3
data_dir=$4

dataset=nmr_ukbb_est
output_file=${dataset}_credible_sets.tsv.gz

mkdir -p "$data_dir/individual"

time python3 scripts/munge_nmr_meta.py \
--input "$parquet" \
--lead-variants "$lead_variants" \
--annotation "$variant_annotation_file" \
--output-dir "$data_dir" \
--dataset $dataset

# sort on chr and pos for tabix, keeping the header first
time cat \
<(echo -n "#") \
<(head -1 "$data_dir/${dataset}_cs_95.tsv") \
<(tail -n +2 "$data_dir/${dataset}_cs_95.tsv" \
| sort -T "$data_dir" -S 4G -k6,6n -k7,7n -k8,8 -k9,9) \
| bgzip -@4 > "$data_dir/$output_file" \
&& tabix -f -s 6 -b 7 -e 7 "$data_dir/$output_file"

# per-trait file with stats, matching the layout the API expects for the other credible set datasets
time python3 <<EOF
import sys
sys.path.insert(0, "scripts")
import polars as pl
from credible_set_stats import calculate_stats, write_stats_json, get_tsv_header, stats_to_tsv_row

data = pl.read_csv("$data_dir/$output_file", separator="\t", null_values=["NA"])
all_stats = []
for (trait,), trait_data in data.partition_by("trait", as_dict=True).items():
    trait_data.write_csv(f"$data_dir/individual/{trait}.SUSIE.munged.tsv", separator="\t", null_value="NA")
    stats = calculate_stats(trait_data)
    write_stats_json(stats, f"$data_dir/individual/{trait}.SUSIE.munged.stats.json")
    all_stats.append(stats)

with open("$data_dir/credible_set_stats.tsv", "w") as f:
    f.write(get_tsv_header() + "\n")
    for s in all_stats:
        f.write(stats_to_tsv_row(s) + "\n")
print(f"Wrote aggregate stats for {len(all_stats)} traits")
EOF
