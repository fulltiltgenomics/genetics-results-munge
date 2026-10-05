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
# Open_Targets_26.06 -> ot_2606_data_studies.json
version=${dataset##*_}
study_json=ot_${version//./}_data_studies.json

mkdir -p $data_dir/opentargets_per_study

# convert the release to one unsorted TSV
time python3 scripts/create_open_targets_files.py $dataset $data_dir

# sort on chr and pos for tabix, keeping the header first
time cat \
<(echo -n "#") \
<(head -1 $data_dir/${dataset}_cs_95.tsv) \
<(tail -n +2 $data_dir/${dataset}_cs_95.tsv \
| sort -T $data_dir -S 4G -k6,6n -k7,7n -k8,8 -k9,9) \
| bgzip -@4 > $data_dir/$output_file \
&& tabix -f -s 6 -b 7 -e 7 $data_dir/$output_file

# create per-study files with stats. The stats count no coding variants here because the
# consequence columns are NA until annotate_resource.sh stamps them and regenerates the stats;
# they are written anyway because that step expects this layout to exist. The aggregate stats
# stay out of opentargets_per_study/ so that directory holds only the per-study files.
time python3 scripts/create_open_targets_per_study_files.py \
    $data_dir/$output_file $data_dir/opentargets_per_study $data_dir/credible_set_stats.tsv

# create study metadata file
time python3 scripts/create_open_targets_study_file.py $data_dir/study_metadata/*.parquet $data_dir/$output_file $data_dir/$study_json
