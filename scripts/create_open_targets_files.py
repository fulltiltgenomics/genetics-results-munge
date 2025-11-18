import os
import sys
import numpy as np
import pyarrow.parquet as pq
import polars as pl


def convert_pq_to_tsv(
    dataset: str, parquet_path: str, anno: pl.DataFrame
) -> dict[str, list[tuple[str, float]]]:
    """
    Convert Open Targets credible set parquet file to 95 % credible set tsv files.
    """
    table = pq.read_table(parquet_path)
    df = table.to_pandas()

    cs_95_outfile = open(parquet_path.replace(".parquet", "_cs_95_noanno.tsv"), "wt")
    cs_95_outfile.write(
        "dataset\tdata_type\ttrait\ttrait_original\tcell_type\tchr\tpos\tref\talt\tmlog10p\tbeta\tse\tpip\tcs_id\tcs_size\tcs_min_r2\taaf\n"
    )
    # skipping 99 % credible sets as the data doesn't seem to differentiate between 95 % and 99 % credible sets
    # cs_99_outfile = open(parquet_path.replace(".parquet", "_cs_99_noanno.tsv"), "wt")
    # cs_99_outfile.write(
    #     "dataset\tdata_type\ttrait\ttrait_original\tcell_type\tchr\tpos\tref\talt\tmlog10p\tbeta\tse\tpip\tcs_id\tcs_size\tcs_min_r2\n"
    # )

    for _, row in df.iterrows():
        study_id = row["studyId"]
        data_type = row["studyType"].upper()
        method = row["finemappingMethod"]
        if (
            data_type != "GWAS"
            or "FINNGEN" in study_id
            or not "susie" in method.lower()
        ):
            continue
        cs_id = row["studyLocusId"]
        if row["locus"] is not None:
            variants_95 = []
            variants_99 = []
            for variant in row["locus"]:
                variant_id = variant.get("variantId")
                cpra = variant_id.split("_")
                # make sure no alternate contigs or other non-chromosomes
                cpra[0] = str(int(cpra[0].replace("X", "23")))
                beta = (
                    f"{variant.get('beta'):.3e}"
                    if variant.get("beta") is not None
                    else "NA"
                )
                se = (
                    f"{variant.get('standardError'):.3e}"
                    if variant.get("standardError") is not None
                    else "NA"
                )
                pip = (
                    round(variant.get("posteriorProbability"), 4)
                    if variant.get("posteriorProbability") is not None
                    else "NA"
                )
                cs_min_r2 = (
                    round(row["purityMinR2"], 4)
                    if not np.isnan(row["purityMinR2"])
                    else "NA"
                )
                mlog10p = "NA"
                af = "NA"  # not available in the data currently
                mantissa = variant.get("pValueMantissa")
                exponent = variant.get("pValueExponent")
                if mantissa is not None and exponent is not None:
                    mlog10p = round(-np.log10(mantissa) - exponent, 4)
                elif variant_id == row["variantId"]:
                    mantissa = row["pValueMantissa"]
                    exponent = row["pValueExponent"]
                    if mantissa is not None and exponent is not None:
                        mlog10p = round(-np.log10(mantissa) - exponent, 4)
                our_variant = [
                    dataset,
                    data_type,
                    study_id,
                    study_id,
                    "NA",
                    cpra[0],
                    cpra[1],
                    cpra[2],
                    cpra[3],
                    str(mlog10p),
                    str(beta),
                    str(se),
                    str(pip),
                    str(cs_id),
                ]
                if variant.get("is95CredibleSet"):
                    variants_95.append(our_variant)
                if variant.get("is99CredibleSet"):
                    variants_99.append(our_variant)
            cs_size_95 = len(variants_95)
            cs_size_99 = len(variants_99)
            for v in variants_95:
                cs_95_outfile.write(
                    "\t".join(v + [str(cs_size_95), str(cs_min_r2), str(af)]) + "\n"
                )
            # for v in variants_99:
            #     cs_99_outfile.write(
            #         "\t".join(v + [str(cs_size_99), str(cs_min_r2), str(af)]) + "\n"
            #     )
    cs_95_outfile.close()
    # cs_99_outfile.close()

    # read file and join annotations
    df = pl.read_csv(
        parquet_path.replace(".parquet", "_cs_95_noanno.tsv"),
        separator="\t",
        null_values=["NA"],
    ).with_columns(
        pl.concat_str(
            [
                pl.col("chr"),
                pl.col("pos").cast(pl.Utf8),
                pl.col("ref"),
                pl.col("alt"),
            ],
            separator=":",
        ).alias("variant_id"),
        pl.col("beta")
        .map_elements(lambda x: f"{x:.3e}", return_dtype=pl.Utf8)
        .alias("beta"),
        pl.col("se")
        .map_elements(lambda x: f"{x:.3e}", return_dtype=pl.Utf8)
        .alias("se"),
        pl.col("pip").round(4),
        pl.col("cs_min_r2").round(4),
    )
    df = (
        df.drop(
            [
                col
                for col in df.columns
                if col.endswith("most_severe") or col.startswith("aaf")
            ]
        )
        .join(anno, on="variant_id", how="left")
        .drop("variant_id")
    )
    df.write_csv(
        parquet_path.replace(".parquet", "_cs_95.tsv"),
        separator="\t",
        null_value="NA",
    )


if __name__ == "__main__":
    if len(sys.argv) != 4:
        print(
            "Usage: python create_open_targets_files.py <dataset_name> <data_dir> <variant_annotation_file>"
        )
        sys.exit(1)
    dataset = sys.argv[1]
    data_dir = sys.argv[2]
    variant_annotation_file = sys.argv[3]
    anno = (
        pl.scan_csv(variant_annotation_file, separator="\t", null_values=["NA"])
        .rename({"#variant": "variant_id"})
        .with_columns(
            pl.col("AF")
            .map_elements(lambda x: f"{x:.3e}", return_dtype=pl.Utf8)
            .alias("aaf"),
        )
        .select("variant_id", "aaf", "most_severe", "gene_most_severe")
        .collect()
    )
    files = [f for f in os.listdir(data_dir) if f.endswith(".parquet")]
    for file in files:
        convert_pq_to_tsv(dataset, f"{data_dir}/{file}", anno)
