#!/bin/bash
# Build the gene-indexed copy of the Open4Gene peak-to-gene table
# (open4gene.all.results.sig.by_gene.tsv.gz) that backs the API's /gene_to_peaks endpoint.
#
# Run inside the genetics-results-munge Docker image (needs bgzip/tabix + polars).
# Uploading to GCS is off by default; pass --stage to publish to both profile buckets.
#
# Usage: scripts/create_gene_indexed_peak_gene_file.sh <data_dir> [--stage]

set -euxo pipefail

if [ $# -lt 1 ]; then
    echo "Usage: $0 <data_dir> [--stage]"
    exit 1
fi
data_dir=$1
stage=${2:-}

open4gene_gs=gs://cascade-browser/cascade_results/open4gene.all.results.sig.tsv.gz
gencode_gs=gs://finngen-commons/results_api_data/mapping_files/gencode.v32.annotation.genes.tsv

finngen_out_gs=gs://finngen-commons/results_api_data/atacseq/
daly_out_gs=gs://daly-genetics-results/atacseq/

open4gene=$data_dir/$(basename $open4gene_gs)
gencode=$data_dir/$(basename $gencode_gs)
out_tsv=$data_dir/open4gene.all.results.sig.by_gene.tsv
out_gz=$out_tsv.gz

mkdir -p $data_dir

for src in $open4gene_gs $gencode_gs; do
    dst=$data_dir/$(basename $src)
    [ -f $dst ] || gcloud storage cp $src $dst
done

time python3 scripts/create_gene_indexed_peak_gene_file.py \
    --open4gene $open4gene \
    --gene-mapping $gencode \
    --output $out_tsv

# tabix needs the rows ordered by the appended gene locus (cols 21-23); the peak keys after
# it keep each gene's block in peak order, which is the order the API merges resources in
time (
    head -1 $out_tsv
    tail -n +2 $out_tsv \
    | LC_ALL=C sort -T $data_dir -S 2G -k21,21n -k22,22n -k23,23n -k1,1 -k2,2n -k3,3n
) | bgzip -@4 > $out_gz

tabix -f -s 21 -b 22 -e 23 $out_gz

echo "rows out: $(($(wc -l < $out_tsv) - 1))"

if [ "$stage" = "--stage" ]; then
    gcloud storage cp $out_gz $out_gz.tbi $finngen_out_gs
    gcloud storage cp $out_gz $out_gz.tbi $daly_out_gs
else
    echo "not staged; to publish:"
    echo "  gcloud storage cp $out_gz $out_gz.tbi $finngen_out_gs"
    echo "  gcloud storage cp $out_gz $out_gz.tbi $daly_out_gs"
fi
