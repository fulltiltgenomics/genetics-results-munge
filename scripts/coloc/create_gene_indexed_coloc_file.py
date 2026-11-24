"""Build the gene-indexed copy of a munged colocalization file (`*.gene.tsv`) served by the
API's `/colocalization_by_gene/{gene}` endpoint.

The munged coloc file is indexed on the colocalizing region, so a gene's colocalizations
can only be found there by already knowing where its QTL signals are, and a trans QTL sits
on another chromosome altogether. This copy carries every pair that has an eQTL, pQTL or
sQTL on one side, once per such side, with the coordinates of that QTL's gene appended, and
is sorted and tabix-indexed on them: one range read at the gene returns its cis and trans
colocalizations alike.

Pairs with no gene QTL side (GWAS-GWAS, metaboQTL, caQTL) are left out. A caQTL trait is a
peak, and linking it to genes needs the peak-to-gene table, as
`create_caqtl_gene_indexed_qtl_file.py` does for credible sets.

Gene resolution, first match wins (`resolve_gene`):
  1. an Ensembl gene id in the QTL's trait_original, or as its trait
  2. `--trait-gene-map`: the QTL's trait_original looked up in a two-column table, for
     platforms whose trait is a reagent id (SomaScan `seq.*`)
  3. the trait as a gene symbol
  4. the trait with an Olink batch suffix (`_cross_batch_normalised_<n>`) removed
  5. the part of the trait before its first `.`, for SomaScan `SYMBOL.<aptamer>` names
Rows that resolve to nothing go to `<output>.unmapped.tsv`, and the counts are printed by
dataset so a platform that silently stopped resolving shows up.

Gene coordinates MUST come from the GENCODE file the API looks query genes up in for coloc
(`gene_gencode_version` on the API's coloc profile entry): the API keeps the rows whose
gene_start and gene_end equal the query gene's, so coordinates from another version match
nothing.
"""

import argparse

import polars as pl

from sumstat_utils import chrom_code_expr

GENE_QTL_TYPES = ["eQTL", "pQTL", "sQTL"]
GENE_COLS = ["qtl_gene", "qtl_side", "gene_chr", "gene_start", "gene_end"]


def read_gene_table(path: str) -> pl.DataFrame:
    """GENCODE genes on chromosomes 1-23 with `ensg` (unversioned id), one row per id."""
    return (
        pl.read_csv(path, separator="\t", infer_schema_length=0)
        .with_columns(chrom_code_expr("chrom").cast(pl.Int8, strict=False).alias("gene_chr"))
        .filter(pl.col("gene_chr") < 24)
        .select(
            pl.col("gene_id").str.split(".").list.first().alias("ensg"),
            "gene_name",
            "gene_chr",
            pl.col("gene_start").cast(pl.Int32),
            pl.col("gene_end").cast(pl.Int32),
        )
        .sort("gene_chr", "gene_start")
        .unique(subset="ensg", keep="first", maintain_order=True)
    )


def resolve_gene(
    side: pl.DataFrame, genes: pl.DataFrame, trait_gene_map: pl.DataFrame | None
) -> pl.DataFrame:
    """Add GENE_COLS (minus qtl_side) to rows carrying `trait` and `trait_original`.

    A row that resolves to no gene keeps nulls in them.
    """
    by_ensg = genes.select(
        "ensg", pl.col("gene_name").alias("n0"), pl.col("gene_chr").alias("c0"),
        pl.col("gene_start").alias("s0"), pl.col("gene_end").alias("e0"),
    )
    # a symbol two ids share resolves to the first by position, the same row every run
    by_symbol = genes.unique(subset="gene_name", keep="first", maintain_order=True)

    candidates = {
        "k1": pl.col("mapped_symbol") if trait_gene_map is not None else pl.lit(None, pl.Utf8),
        "k2": pl.col("trait"),
        "k3": pl.col("trait").str.replace(r"_cross_batch_normalised_\d+$", ""),
        "k4": pl.col("trait").str.extract(r"^([^.]+)\."),
    }
    out = side.with_columns(
        pl.coalesce(
            pl.col("trait_original").str.extract(r"(ENSG\d+)"),
            pl.col("trait").str.extract(r"^(ENSG\d+)"),
        ).alias("ensg")
    )
    if trait_gene_map is not None:
        out = out.join(trait_gene_map, on="trait_original", how="left")
    out = out.with_columns(**candidates).join(by_ensg, on="ensg", how="left")
    for k in candidates:
        i = k[1]
        out = out.join(
            by_symbol.select(
                pl.col("gene_name").alias(k), pl.col("gene_name").alias(f"n{i}"),
                pl.col("gene_chr").alias(f"c{i}"), pl.col("gene_start").alias(f"s{i}"),
                pl.col("gene_end").alias(f"e{i}"),
            ),
            on=k,
            how="left",
        )
    order = ["0", "1", "2", "3", "4"]
    return out.with_columns(
        pl.coalesce([f"n{i}" for i in order]).alias("qtl_gene"),
        pl.coalesce([f"c{i}" for i in order]).alias("gene_chr"),
        pl.coalesce([f"s{i}" for i in order]).alias("gene_start"),
        pl.coalesce([f"e{i}" for i in order]).alias("gene_end"),
    ).select(side.columns + ["qtl_gene", "gene_chr", "gene_start", "gene_end"])


def gene_indexed(
    coloc: pl.DataFrame, genes: pl.DataFrame, trait_gene_map: pl.DataFrame | None
) -> pl.DataFrame:
    """Every coloc row once per gene QTL side, with GENE_COLS appended; unresolved rows keep
    null gene columns."""
    sides = []
    for n in (1, 2):
        rows = coloc.filter(pl.col(f"data_type{n}").is_in(GENE_QTL_TYPES)).with_columns(
            pl.col(f"trait{n}").alias("trait"),
            pl.col(f"trait{n}_original").alias("trait_original"),
        )
        resolved = resolve_gene(rows, genes, trait_gene_map).with_columns(
            pl.lit(n, pl.Int8).alias("qtl_side")
        )
        sides.append(resolved.select(coloc.columns + GENE_COLS))
    return pl.concat(sides)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--coloc", required=True, help="munged colocQC file (.tsv or .tsv.gz)")
    ap.add_argument("--gene-mapping", required=True, help="gencode.vNN.annotation.genes.tsv")
    ap.add_argument(
        "--trait-gene-map",
        help="headerless TSV of trait_original and gene symbol, for traits that name a reagent",
    )
    ap.add_argument("--output", required=True, help="unsorted output .tsv")
    args = ap.parse_args()

    # every column as a string: the rows are copied through and must not be reformatted
    coloc = pl.read_csv(args.coloc, separator="\t", infer_schema_length=0)
    genes = read_gene_table(args.gene_mapping)
    trait_gene_map = None
    if args.trait_gene_map:
        trait_gene_map = pl.read_csv(
            args.trait_gene_map, separator="\t", has_header=False, infer_schema_length=0,
            new_columns=["trait_original", "mapped_symbol"],
        ).unique(subset="trait_original", keep="first", maintain_order=True)

    out = gene_indexed(coloc, genes, trait_gene_map)
    mapped = out.filter(pl.col("gene_chr").is_not_null())
    unmapped = out.filter(pl.col("gene_chr").is_null())
    mapped.write_csv(args.output, separator="\t")
    unmapped.drop(GENE_COLS).write_csv(
        args.output.removesuffix(".tsv") + ".unmapped.tsv", separator="\t"
    )

    print(f"rows in: {coloc.height}, gene QTL sides: {out.height}, mapped: {mapped.height}")
    per_dataset = (
        out.with_columns(
            pl.when(pl.col("qtl_side") == 1).then(pl.col("#dataset1")).otherwise(pl.col("dataset2"))
            .alias("dataset"),
            pl.when(pl.col("qtl_side") == 1).then(pl.col("trait1")).otherwise(pl.col("trait2"))
            .alias("trait"),
        )
        .group_by("dataset")
        .agg(
            pl.len().alias("rows"),
            pl.col("gene_chr").is_null().sum().alias("unmapped_rows"),
            pl.col("trait").filter(pl.col("gene_chr").is_null()).n_unique().alias("unmapped_traits"),
            pl.col("trait").filter(pl.col("gene_chr").is_null()).first().alias("example"),
        )
        .filter(pl.col("unmapped_rows") > 0)
        .sort("unmapped_rows", descending=True)
    )
    print(f"datasets with unmapped rows: {per_dataset.height}; the largest:")
    with pl.Config(tbl_rows=15, tbl_width_chars=200, fmt_str_lengths=60):
        print(per_dataset.head(15))


if __name__ == "__main__":
    main()
