#!/usr/bin/env python3
### this script creates a gene name mapping file across all gencode versions
### output file:
### ensg    gene_name_49    gene_name_43    gene_name_39    gene_name_32
### ENSG00000290825 DDX11L16        DDX11L2 NA      NA
### ENSG00000223972 DDX11L1 DDX11L1 DDX11L1 DDX11L1
### ...
### ENSG00000234290 NA      NA      NA      AC116366.1
###
### the version list must match `gencode_versions` in the results-api genes config,
### which reads both this file (by its version-stamped name) and one
### gencode.v{version}.annotation.genes.tsv per version. usage:
###   python3 create_gene_name_mapping_across_gencode_versions.py [gencode_dir] [output_dir]

import sys
from pathlib import Path

import polars as pl

# descending; the API takes the first as its default version
gencode_versions = ["49", "45", "43", "39", "35", "32", "19"]

gencode_dir = Path(sys.argv[1] if len(sys.argv) > 1 else "/mnt/disks/data")
output_dir = Path(sys.argv[2] if len(sys.argv) > 2 else gencode_dir)

data = None
for version in gencode_versions:
    df = (
        pl.read_csv(
            gencode_dir / f"gencode.v{version}.annotation.genes.tsv",
            separator="\t",
        )
        .with_columns(pl.col("gene_id").str.split(".").list.first().alias("ensg"))
        .select(["ensg", "gene_name"])
        .rename({"gene_name": f"gene_name_{version}"})
    )
    if data is None:
        data = df
    else:
        data = data.join(df, on="ensg", how="full", coalesce=True)

output_file = output_dir / f"gencode_gene_name_mapping_{'-'.join(gencode_versions)}.tsv"
data.write_csv(output_file, separator="\t", null_value="NA")
print(f"wrote {data.height} genes to {output_file}")
