"""
Module for calculating credible set statistics (risk/protective counts with coding/LoF breakdowns).
"""

import json
import polars as pl
from dataclasses import dataclass, asdict

CODING_VARIANTS = {
    "missense_variant",
    "synonymous_variant",
    "stop_gained",
    "stop_lost",
    "start_lost",
    "frameshift_variant",
    "inframe_insertion",
    "inframe_deletion",
    "splice_acceptor_variant",
    "splice_donor_variant",
}

LOF_VARIANTS = {
    "stop_gained",
    "stop_lost",
    "start_lost",
    "frameshift_variant",
    "splice_acceptor_variant",
    "splice_donor_variant",
}


@dataclass
class CredibleSetStats:
    trait: str
    trait_original: str
    dataset: str
    data_type: str
    n_risk_cs: int
    n_risk_cs_with_coding: int
    n_risk_cs_with_coding_pip_gt_0_05: int
    n_risk_cs_with_lof: int
    n_risk_cs_with_lof_pip_gt_0_05: int
    n_protective_cs: int
    n_protective_cs_with_coding: int
    n_protective_cs_with_coding_pip_gt_0_05: int
    n_protective_cs_with_lof: int
    n_protective_cs_with_lof_pip_gt_0_05: int


def calculate_stats(df: pl.DataFrame) -> CredibleSetStats:
    """
    Calculate credible set statistics from a per-trait DataFrame.

    Lead variant = variant with highest PIP in each credible set.
    Risk classification:
      - If AAF available: (beta > 0 AND aaf < 0.5) OR (beta < 0 AND aaf > 0.5)
      - If AAF is NA: beta > 0
    Protective classification: complement of risk.

    Always returns a CredibleSetStats object (with 0s if no credible sets).
    """
    trait = df["trait"][0] if len(df) > 0 else ""
    trait_original = (
        df["trait_original"][0]
        if len(df) > 0 and "trait_original" in df.columns
        else trait
    )
    if "dataset" in df.columns:
        dataset = df["dataset"][0] if len(df) > 0 else ""
    elif "#dataset" in df.columns:
        dataset = df["#dataset"][0] if len(df) > 0 else ""
    else:
        raise KeyError("Neither 'dataset' nor '#dataset' found in DataFrame columns.")
    data_type = df["data_type"][0] if len(df) > 0 else ""

    if len(df) == 0:
        return CredibleSetStats(
            trait=trait,
            trait_original=trait_original,
            dataset=dataset,
            data_type=data_type,
            n_risk_cs=0,
            n_risk_cs_with_coding=0,
            n_risk_cs_with_coding_pip_gt_0_05=0,
            n_risk_cs_with_lof=0,
            n_risk_cs_with_lof_pip_gt_0_05=0,
            n_protective_cs=0,
            n_protective_cs_with_coding=0,
            n_protective_cs_with_coding_pip_gt_0_05=0,
            n_protective_cs_with_lof=0,
            n_protective_cs_with_lof_pip_gt_0_05=0,
        )

    # cast string columns to float if needed
    if df["pip"].dtype == pl.Utf8:
        df = df.with_columns(pl.col("pip").cast(pl.Float64, strict=False))
    if df["beta"].dtype == pl.Utf8:
        df = df.with_columns(pl.col("beta").cast(pl.Float64, strict=False))
    if "aaf" in df.columns and df["aaf"].dtype == pl.Utf8:
        df = df.with_columns(pl.col("aaf").cast(pl.Float64, strict=False))

    # find lead variant per credible set (highest PIP)
    lead_variants = df.sort("pip", descending=True).group_by("cs_id").first()

    # classify credible sets as risk or protective based on lead variant
    # when aaf is available: (beta > 0 AND aaf < 0.5) OR (beta < 0 AND aaf > 0.5) = risk
    # when aaf is NA: beta > 0 = risk
    if "aaf" in lead_variants.columns:
        lead_variants = lead_variants.with_columns(
            pl.when(pl.col("aaf").is_not_null())
            .then(
                pl.when(
                    ((pl.col("beta") > 0) & (pl.col("aaf") < 0.5))
                    | ((pl.col("beta") < 0) & (pl.col("aaf") > 0.5))
                )
                .then(pl.lit("risk"))
                .otherwise(pl.lit("protective"))
            )
            .otherwise(
                pl.when(pl.col("beta") > 0)
                .then(pl.lit("risk"))
                .otherwise(pl.lit("protective"))
            )
            .alias("direction")
        )
    else:
        # no aaf column, use beta only
        lead_variants = lead_variants.with_columns(
            pl.when(pl.col("beta") > 0)
            .then(pl.lit("risk"))
            .otherwise(pl.lit("protective"))
            .alias("direction")
        )

    risk_cs_ids = set(
        lead_variants.filter(pl.col("direction") == "risk")["cs_id"].to_list()
    )
    protective_cs_ids = set(
        lead_variants.filter(pl.col("direction") == "protective")["cs_id"].to_list()
    )

    # get all variants in risk/protective credible sets for coding/LoF analysis
    risk_variants = df.filter(pl.col("cs_id").is_in(risk_cs_ids))
    protective_variants = df.filter(pl.col("cs_id").is_in(protective_cs_ids))

    def count_cs_with_variants(
        variants_df: pl.DataFrame, variant_set: set, pip_threshold: float | None = None
    ) -> int:
        """Count unique credible sets containing variants from variant_set."""
        if "most_severe" not in variants_df.columns:
            return 0
        filtered = variants_df.filter(pl.col("most_severe").is_in(variant_set))
        if pip_threshold is not None:
            filtered = filtered.filter(pl.col("pip") > pip_threshold)
        return filtered["cs_id"].n_unique()

    return CredibleSetStats(
        trait=trait,
        trait_original=trait_original,
        dataset=dataset,
        data_type=data_type,
        n_risk_cs=len(risk_cs_ids),
        n_risk_cs_with_coding=count_cs_with_variants(risk_variants, CODING_VARIANTS),
        n_risk_cs_with_coding_pip_gt_0_05=count_cs_with_variants(
            risk_variants, CODING_VARIANTS, 0.05
        ),
        n_risk_cs_with_lof=count_cs_with_variants(risk_variants, LOF_VARIANTS),
        n_risk_cs_with_lof_pip_gt_0_05=count_cs_with_variants(
            risk_variants, LOF_VARIANTS, 0.05
        ),
        n_protective_cs=len(protective_cs_ids),
        n_protective_cs_with_coding=count_cs_with_variants(
            protective_variants, CODING_VARIANTS
        ),
        n_protective_cs_with_coding_pip_gt_0_05=count_cs_with_variants(
            protective_variants, CODING_VARIANTS, 0.05
        ),
        n_protective_cs_with_lof=count_cs_with_variants(
            protective_variants, LOF_VARIANTS
        ),
        n_protective_cs_with_lof_pip_gt_0_05=count_cs_with_variants(
            protective_variants, LOF_VARIANTS, 0.05
        ),
    )


def write_stats_json(stats: CredibleSetStats, output_path: str) -> None:
    """Write statistics to JSON file."""
    with open(output_path, "w") as f:
        json.dump(asdict(stats), f, indent=2)


def stats_to_tsv_row(stats: CredibleSetStats) -> str:
    """Convert stats to TSV row string."""
    d = asdict(stats)
    return "\t".join(str(v) for v in d.values())


def get_tsv_header() -> str:
    """Get TSV header for stats file."""
    return "\t".join(
        [
            "trait",
            "trait_original",
            "dataset",
            "data_type",
            "n_risk_cs",
            "n_risk_cs_with_coding",
            "n_risk_cs_with_coding_pip_gt_0_05",
            "n_risk_cs_with_lof",
            "n_risk_cs_with_lof_pip_gt_0_05",
            "n_protective_cs",
            "n_protective_cs_with_coding",
            "n_protective_cs_with_coding_pip_gt_0_05",
            "n_protective_cs_with_lof",
            "n_protective_cs_with_lof_pip_gt_0_05",
        ]
    )
