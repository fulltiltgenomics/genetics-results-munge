#!/bin/bash

set -euxo pipefail

mkdir -p data
cd data

while read url; do
    until curl -L -o $(basename $url) $url; do
        echo "Download failed, retrying in 10 seconds..."
        sleep 10
    done
done < <(cut -f11 ../metadata/eqtl_catalogue_tabix_ftp_paths.tsv | tail -n+2)

# molecular trait metadata
urls=(
    "https://zenodo.org/records/7808390/files/Affy_Human_Gene_1_0_ST_Ensembl_96_phenotype_metadata.tsv.gz?download=1"
    "https://zenodo.org/records/7808390/files/exon_counts_Ensembl_105_phenotype_metadata.tsv.gz?download=1"
    "https://zenodo.org/records/7808390/files/gene_counts_Ensembl_105_phenotype_metadata.tsv.gz?download=1"
    "https://zenodo.org/records/7808390/files/HumanHT-12_V4_Ensembl_96_phenotype_metadata.tsv.gz?download=1"
    "https://zenodo.org/records/7808390/files/SomaLogic_Ensembl_96_phenotype_metadata.tsv.gz?download=1"
    "https://zenodo.org/records/7808390/files/transcript_usage_Ensembl_105_phenotype_metadata.tsv.gz?download=1"
    "https://zenodo.org/records/7808390/files/txrevise_Ensembl_105_phenotype_metadata.tsv.gz?download=1"
)
for url in "${urls[@]}"; do
    curl -L -o $(basename $url | sed 's/?download=1$//') $url
done

# leafcutter sQTL metadata is stored per study
while read study; do
    curl -L -o leafcutter_${study}_Ensembl_105_phenotype_metadata.tsv.gz "https://zenodo.org/records/7850746/files/leafcutter_${study}_Ensembl_105_phenotype_metadata.tsv.gz?download=1"
done < <(awk -F'\t' '$9=="leafcutter"' ../metadata/eqtl_catalogue_studies.tsv | cut -f2 | sort)
