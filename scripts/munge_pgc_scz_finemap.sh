#!/bin/bash

set -euxo pipefail

if [ $# -ne 4 ]; then
    echo "Usage: $0 <st11a_tsv> <munged_sumstats> <variant_annotation_file> <data_dir>"
    echo "  <st11a_tsv>               ST11a_95_perc_Credible_Sets.tsv from Trubetskoy et al. 2022"
    echo "  <munged_sumstats>         daner_PGC_SCZ_w3_90_0418b.munged.tsv.gz"
    echo "  <variant_annotation_file> tabix indexed FinnGen annotated variants"
    exit 1
fi
st11a=$1
sumstats=$2
variant_annotation_file=$3
data_dir=$4

dataset=PGC_SCZ_2022
output_file=${dataset}_credible_sets.tsv.gz

mkdir -p "$data_dir/individual"

time python3 scripts/munge_pgc_scz_finemap.py \
--input "$st11a" \
--sumstats "$sumstats" \
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
    trait_data.write_csv(f"$data_dir/individual/{trait}.FINEMAP.munged.tsv", separator="\t", null_value="NA")
    stats = calculate_stats(trait_data)
    write_stats_json(stats, f"$data_dir/individual/{trait}.FINEMAP.munged.stats.json")
    all_stats.append(stats)

with open("$data_dir/credible_set_stats.tsv", "w") as f:
    f.write(get_tsv_header() + "\n")
    for s in all_stats:
        f.write(stats_to_tsv_row(s) + "\n")
print(f"Wrote aggregate stats for {len(all_stats)} traits")
EOF
