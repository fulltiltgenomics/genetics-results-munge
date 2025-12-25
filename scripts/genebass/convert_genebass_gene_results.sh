#!/bin/bash

set -euxo pipefail

time python3 scripts/genebass/convert_genebass_gene_results.py \
--input /mnt/disks/data/results.mt \
--output gene_burden_results.tsv

bgzip -f gene_burden_results.tsv
tabix -f -s5 -b6 -e6 gene_burden_results.tsv.gz
