import polars as pl

MLOG10P_THRESHOLD = 4

data = (
    pl.scan_csv(
        "gs://finngen-commons/genebass/variant_results_mlog10p3.tsv.bgz",
        separator="\t",
        schema_overrides={
            "mlog10p": pl.Utf8,
        },
        null_values=["NA"],
    )
    .with_columns(pl.lit("genebass").alias("#dataset"))
    .with_columns(
        pl.when(pl.col("mlog10p").str.to_lowercase().str.starts_with("inf"))
        .then(324)  # se is 0 for these so can't compute real mlog10p
        .otherwise(pl.col("mlog10p").cast(pl.Float64))
        .alias("mlog10p")
    )
    .with_columns(pl.col("mlog10p").round(4).alias("mlog10p"))
    .with_columns(
        pl.col("beta")
        .map_elements(lambda x: f"{x:.3e}", return_dtype=pl.Utf8)
        .alias("beta")
    )
    .with_columns(
        pl.col("se")
        .map_elements(lambda x: f"{x:.3e}", return_dtype=pl.Utf8)
        .alias("se")
    )
    .with_columns(
        pl.col("af_overall")
        .map_elements(lambda x: f"{x:.3e}", return_dtype=pl.Utf8)
        .alias("af_overall")
    )
    .with_columns(
        pl.col("af_cases")
        .map_elements(lambda x: f"{x:.3e}", return_dtype=pl.Utf8)
        .alias("af_cases")
    )
    .with_columns(
        pl.col("af_controls")
        .map_elements(lambda x: f"{x:.3e}", return_dtype=pl.Utf8)
        .alias("af_controls")
    )
    .with_columns(
        pl.col("heritability")
        .map_elements(lambda x: f"{x:.3e}", return_dtype=pl.Utf8)
        .alias("heritability")
    )
    .select(
        [
            "#dataset",
            "chr",
            "pos",
            "ref",
            "alt",
            "gene",
            "annotation",
            "mlog10p",
            "beta",
            "se",
            "af_overall",
            "af_cases",
            "af_controls",
            "ac",
            "an",
            "heritability",
            "trait",
        ]
    )
    .filter(pl.col("mlog10p") > MLOG10P_THRESHOLD)
    .collect()
)

data.write_csv(
    f"genebass_variant_results_mlog10p{MLOG10P_THRESHOLD}.tsv",
    separator="\t",
    null_value="NA",
)

for trait in data["trait"].unique():
    study_data = data.filter(pl.col("trait") == trait)
    study_data.write_csv(
        f"genebass_per_trait/{trait}.mlog10p{MLOG10P_THRESHOLD}.tsv",
        separator="\t",
        null_value="NA",
    )
