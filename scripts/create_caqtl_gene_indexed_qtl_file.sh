#!/bin/bash
# Build the gene-indexed caQTL credible set file (FinnGen_ATACseq_*_credible_sets.qtl.tsv.gz)
# that backs the API's /credible_sets_by_qtl_gene endpoint for finngen_caqtl.
#
# Run inside the genetics-results-munge Docker image (needs bgzip/tabix + polars).
# Needs ~15 GB free disk in <data_dir>: the uncompressed credible sets are ~4 GB, the
# gene-expanded output ~2 GB, plus sort temp space.
#
# Uploading to GCS is off by default; pass --stage to publish to both profile buckets.
#
# Usage: scripts/create_caqtl_gene_indexed_qtl_file.sh <data_dir> [--stage]

set -euxo pipefail

if [ $# -lt 1 ]; then
    echo "Usage: $0 <data_dir> [--stage]"
    exit 1
fi
data_dir=$1
stage=${2:-}

cs_gs=gs://finngen-commons/results_api_data/credible_sets/finngen_atacseq/20251118/FinnGen_ATACseq_202509_credible_sets.tsv.gz
open4gene_gs=gs://cascade-browser/cascade_results/open4gene.all.results.sig.tsv.gz
gencode_gs=gs://finngen-commons/results_api_data/mapping_files/gencode.v32.annotation.genes.tsv

finngen_out_gs=gs://finngen-commons/results_api_data/credible_sets/finngen_atacseq/20251118/
daly_out_gs=gs://daly-genetics-results/credible_sets/finngen_atacseq/20251118/

cs_gz=$data_dir/$(basename $cs_gs)
cs_tsv=${cs_gz%.gz}
open4gene=$data_dir/$(basename $open4gene_gs)
gencode=$data_dir/$(basename $gencode_gs)
out_tsv=$data_dir/FinnGen_ATACseq_202509_credible_sets.qtl.tsv
out_gz=$out_tsv.gz

mkdir -p $data_dir

for src in $cs_gs $open4gene_gs $gencode_gs; do
    dst=$data_dir/$(basename $src)
    [ -f $dst ] || gcloud storage cp $src $dst
done

# polars cannot stream a compressed CSV, so the credible sets are read uncompressed
[ -f $cs_tsv ] || bgzip -dk -@4 $cs_gz

time python3 scripts/create_caqtl_gene_indexed_qtl_file.py \
    --credible-sets $cs_tsv \
    --open4gene $open4gene \
    --gene-mapping $gencode \
    --output $out_tsv

# tabix needs the rows ordered by the gene locus (cols 20-22); the variant keys after it
# reproduce the order the API merges resources in (chr, pos, ref, alt, trait). LC_ALL=C
# makes the string keys sort byte-wise, which is how the API compares them.
time (
    head -1 $out_tsv
    tail -n +2 $out_tsv \
    | LC_ALL=C sort -T $data_dir -S 4G -k20,20g -k21,21g -k22,22g -k6,6g -k7,7g -k8,8 -k9,9 -k3,3
) | bgzip -@4 > $out_gz

tabix -f -s 20 -b 21 -e 21 $out_gz

echo "rows in: $(($(wc -l < $cs_tsv) - 1)), rows out: $(($(wc -l < $out_tsv) - 1))"

if [ "$stage" = "--stage" ]; then
    gcloud storage cp $out_gz $out_gz.tbi $finngen_out_gs
    gcloud storage cp $out_gz $out_gz.tbi $daly_out_gs
else
    echo "not staged; to publish:"
    echo "  gcloud storage cp $out_gz $out_gz.tbi $finngen_out_gs"
    echo "  gcloud storage cp $out_gz $out_gz.tbi $daly_out_gs"
fi
