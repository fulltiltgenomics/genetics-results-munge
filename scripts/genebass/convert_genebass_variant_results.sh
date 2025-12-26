#!/bin/bash

# first run convert_genebass_variant_results_dataproc.py in dataproc, then this script

set -euxo pipefail

mkdir -p genebass_per_trait
python3 scripts/genebass/cleanup_genebass_variant_results.py
bgzip -@4 genebass_variant_results_mlog10p4.tsv
tabix -f -s2 -b3 -e3 genebass_variant_results_mlog10p4.tsv.gz
