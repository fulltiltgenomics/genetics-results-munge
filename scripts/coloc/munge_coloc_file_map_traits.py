import polars as pl
import glob
import os
import sys

# Second pass over the colocQC munged intermediate: split trait1/trait2 into
# trait/trait_original and map eQTL Catalogue molecular-trait IDs (trait2) to gene names.
#
# R14 layout assumptions (verified against the data):
#   - side 1 is always FinnGen GWAS -> trait1 is an endpoint, never gene-mapped
#   - eQTL Catalogue datasets are on side 2 and carry a QTD... dataset id (set by
#     munge_coloc_file.py); their trait2 is a phenotype id (ge/exon/tx/txrev/
#     microarray/aptamer/leafcutter/majiq) that we map to a gene name
#   - non-catalogue side-2 datasets (FinnLiver, FinnGen QTL/GWAS) keep their raw trait
#
# Gene-name mapping uses the eQTL Catalogue phenotype-metadata files next to the
# dataset-metadata file (in <metadata_dir>/metadata/). majiq has no phenotype-metadata
# file: its gene id is the trait prefix (ENSG...:s:...), so it is resolved via the
# gene-level lookup. leafcutter clusters from studies whose per-study phenotype file is
# absent (e.g. R8 INTERVAL_RNA) stay unmapped and keep the raw cluster id.

LEADING_COLS = [
    "#dataset1",
    "dataset2",
    "data_type1",
    "data_type2",
    "trait1",
    "trait1_original",
    "trait2",
    "trait2_original",
    "cell_type1",
    "cell_type2",
]


def _build_lookups(pheno_dir):
    files = glob.glob(os.path.join(pheno_dir, "*_phenotype_metadata.tsv.gz"))
    if not files:
        raise FileNotFoundError(f"no phenotype metadata files found in {pheno_dir}")
    frames = []
    for f in files:
        frames.append(
            pl.read_csv(f, separator="\t", schema_overrides={"chromosome": pl.Utf8})
            .select("phenotype_id", "gene_id", "gene_name")
        )
    combined = pl.concat(frames)
    pheno_lookup = (
        combined.select("phenotype_id", pl.col("gene_name").alias("gene_by_pheno"))
        .unique(subset="phenotype_id")
    )
    # gene-level lookup keyed on bare gene id (strip Ensembl version) for majiq/ge/exon
    gene_lookup = (
        combined.select(
            pl.col("gene_id").str.replace(r"\..*$", "").alias("gene_key"),
            pl.col("gene_name").alias("gene_by_id"),
        )
        .drop_nulls("gene_key")
        .unique(subset="gene_key")
    )
    return pheno_lookup, gene_lookup


def main(filename, eqtl_cat_metadata_filename, output_filename):
    meta_dir = os.path.dirname(os.path.abspath(eqtl_cat_metadata_filename))
    pheno_dir = os.path.join(meta_dir, "metadata")
    pheno_lookup, gene_lookup = _build_lookups(pheno_dir)

    # QTD dataset id -> quant method, for the trait2_original suffix
    qtd_quant = dict(
        pl.read_csv(eqtl_cat_metadata_filename, separator="\t")
        .select("dataset_id", "quant_method")
        .iter_rows()
    )

    data = pl.read_csv(filename, separator="\t", infer_schema_length=100000, null_values=["NA"])

    is_cat = pl.col("dataset2").str.starts_with("QTD")
    data = (
        data.with_columns(
            pl.col("trait1").alias("trait1_original"),
            # bare gene id from an ENSG-prefixed trait (ge/exon/txrev/majiq)
            pl.col("trait2").str.extract(r"^(ENSG\d+)").alias("_gene_key"),
            pl.when(is_cat)
            .then(pl.col("dataset2").replace_strict(qtd_quant, default=None))
            .otherwise(None)
            .alias("_quant2"),
        )
        .join(pheno_lookup, left_on="trait2", right_on="phenotype_id", how="left")
        .join(gene_lookup, left_on="_gene_key", right_on="gene_key", how="left")
        .with_columns(
            pl.when(is_cat)
            .then(pl.concat_str(pl.col("trait2"), pl.col("_quant2"), separator="|"))
            .otherwise(pl.col("trait2"))
            .alias("trait2_original"),
            # mapped gene name for catalogue traits; raw trait kept if unmapped or non-catalogue
            pl.when(is_cat)
            .then(pl.coalesce(pl.col("gene_by_pheno"), pl.col("gene_by_id"), pl.col("trait2")))
            .otherwise(pl.col("trait2"))
            .alias("trait2"),
        )
    )

    rest = [c for c in data.columns if c not in LEADING_COLS and not c.startswith("_") and c not in ("gene_by_pheno", "gene_by_id")]
    data.select(LEADING_COLS + rest).write_csv(output_filename, separator="\t", null_value="NA")


if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python3 munge_coloc_file_map_traits.py <filename> <eqtl_cat_metadata_filename> <output_filename>")
        sys.exit(1)
    main(sys.argv[1], sys.argv[2], sys.argv[3])
