#!/usr/bin/env python3
"""build the two munge_li_brain inputs from the real CATlas human-brain files:
  - wide cCRE x cell-type CPM matrix from the 44 subclass per-cell-type bed files
  - cCRE->gene links TSV from Table S12 (gene-cCRE correlations bedpe): anchorB = cCRE
coords keep the "chr" prefix so the munge script exercises its numeric-chrom conversion.
"""
import glob
import os
import polars as pl

BEDS = "data/li_brain/beds"
S12 = "data/li_brain/gene_cCRE_correlations.bedpe.gz"
OUT_MATRIX = "data/li_brain/human_brain_cCRE_by_celltype_cpm.tsv"
OUT_LINKS = "data/li_brain/human_brain_cCRE_gene_links.tsv"

frames = []
for path in sorted(glob.glob(f"{BEDS}/*.bed")):
    cell = os.path.basename(path)[:-4]
    df = pl.read_csv(path, separator="\t", has_header=False,
                     new_columns=["chrom", "start", "end", "cCRE_id", "score"],
                     schema_overrides={"chrom": pl.Utf8, "start": pl.Int64, "end": pl.Int64, "score": pl.Float64})
    frames.append(df.select("chrom", "start", "end", pl.lit(cell).alias("cell_type"), "score"))

long = pl.concat(frames)
print(f"long entries across {len(frames)} cell types: {long.height}")

wide = long.pivot(index=["chrom", "start", "end"], on="cell_type", values="score", aggregate_function="first")
print(f"wide matrix: {wide.height} unique cCREs x {wide.width - 3} cell types")
wide.write_csv(OUT_MATRIX, separator="\t")

links = pl.read_csv(S12, separator="\t", has_header=False,
                    columns=[3, 4, 5, 6],
                    new_columns=["chrom", "start", "end", "name"],
                    schema_overrides={"column_4": pl.Utf8, "column_7": pl.Utf8})
links = links.with_columns(pl.col("name").str.split("|").list.first().alias("gene_name")) \
             .select("chrom", "start", "end", "gene_name").unique()
print(f"links: {links.height} unique (cCRE, gene) pairs; {links['gene_name'].n_unique()} genes")
links.write_csv(OUT_LINKS, separator="\t")
print("wrote", OUT_MATRIX, "and", OUT_LINKS)
