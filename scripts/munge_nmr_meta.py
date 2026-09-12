"""
Convert the published EstBB-UKBB NMR metabolic trait SuSiE fine-mapping (Tambets et al. 2026,
Zenodo 10.5281/zenodo.18132538) to the credible set TSV format used by the rest of this
repository.

Input
  --input          UKBB_EUR_fine_mapping_with_meta_EUR_lead_variants.parquet, one row per
                   credible set variant, GRCh38
  --lead-variants  EUR_all_lead_variants.tsv.gz from the companion record
                   (Zenodo 10.5281/zenodo.18377015), the same cohort's published BETA, SE
                   and LOG10P at every genome-wide significant lead variant
  --annotation     FinnGen annotated variants (tabix indexed), for most_severe /
                   gene_most_severe and for deciding which allele the published MAF belongs to

Output
  <output_dir>/<dataset>_cs_95.tsv   unsorted, with a header; the shell driver sorts, bgzips
                                     and tabix indexes it

The fine-mapping is SuSiE run on the UKBB_EUR subset alone (413,897 samples) at regions defined
by the EstBB-UKBB meta_EUR lead variants, so --n-samples describes UKBB_EUR and not the 619,372
of the meta-analysis. The metabolic traits are inverse-normal transformed, which is what makes
the beta/se derivation below meaningful: the effect is in trait standard deviations.

THE FINE-MAPPING FILE'S z IS THE REFERENCE ALLELE'S, NOT THE ALTERNATIVE ALLELE'S -- and the
lead-variant file's is the other way round. Both label their alleles against GRCh38 and agree on
which is which (`ALL0` == `ref`, `ALL1` == `alt`, on all 49,815 rows they share), but their
z-scores are exact negatives of each other: over the 48,742 shared rows with a finite z,
max |Z_lead + z_fm| = 0.0046, and not one pair shares a sign. The lead file is the one that is
right: its BETA is positive for alt at HMGCR rs12916, which raises LDL-C, and its sign agrees
with published direction at APOE, PCSK9, SORT1, LPL and CETP too.

So `beta` is oriented on `alt`, which means negating the fine-mapping file's z. Taking that sign
at face value would invert every direction of effect in the dataset, silently, in a file that
loads and indexes perfectly.

beta, se and mlog10p come from the LEAD-VARIANT FILE where it has them, and are DERIVED
elsewhere. The fine-mapping file carries no effect size, standard error or p-value at all -- only
a signed z and the minor allele frequency. For an inverse-normal transformed trait the
standardised effect follows from those plus the sample size (Zhu et al. 2016):

    se   = s_trait / sqrt(2 f (1 - f) (n + z^2))
    beta = -z * se

`f` is the published UKBB_EUR MAF, and 2f(1-f) is invariant under f -> 1-f, so the unknown allele
orientation of the MAF does not enter. `s_trait` is NOT assumed to be 1: it is fitted per trait
against that trait's published lead variants, as the median of
SE_published * sqrt(2 f (1-f) (n + z^2)). Assuming 1 -- i.e. that the transformed trait has unit
residual variance after age, sex and PC covariates -- makes |beta| and |se| too large by 8 % on
the median trait and by 21 % on the worst (the fitted scale ranges 0.79 to 0.996 across the 249).
Every trait has at least 23 lead variants to fit from, the median has 192, and the fit is tight:
the within-trait spread of the scale is 1.6 % of its median, 4.9 % at worst. mlog10p needs no
scale, since z is exact -- the derived value reproduces the published LOG10P.

`z` OVERFLOWS for the strongest associations. The fine-mapping file's z is a two-sided p-to-z
conversion, and its largest finite value is 38.47 -- exactly where a double-precision p-value
underflows to zero -- so 6,077 of 3,792,183 rows carry +-inf. 1,073 of those are lead variants
and take the published numbers. The other 4,984 get a null `mlog10p`, `beta` and `se`: their true
|z| is not recoverable, and the derivation would not fail on them, it would CONVERGE -- as
|z| -> inf, |beta| -> s_trait/sqrt(2f(1-f)), a finite and entirely fictional effect size. Their
`pip` and credible set membership are unaffected.

SuSiE's own alphas do not rescue those 4,984 either, and the reason is worth recording so nobody
re-derives it: log(alpha) is very nearly linear in z^2 within a credible set (median R^2 0.989
over 4,000 sets), so an overflowed z can in principle be inverted from the set's other members --
but only 11 of the 2,585 credible sets containing an overflowed row have any member with a finite
z, covering 120 rows. 1,753 of those sets are singletons. The information is not there.

`pip` is the credible set's own alpha, not the overall PIP. SuSiE reports both: `pip` in the
source file is 1 - prod(1 - alpha_l) over all L single effects, while `alpha{cs_index}` is the
variant's posterior weight within the one signal this credible set represents. The rest of this
repository writes the latter (`cs_specific_prob` in the FinnGen fine-mapping pipeline), and it is
also what makes the per-set probabilities sum to the ~0.95 that "95 % credible set" names. The
two differ by a median of 0.0004 and by more than 0.01 on 286 rows.

`aaf` keeps the study's own frequency and takes only its ORIENTATION from the annotation. The
source publishes MAF, which has lost the allele, so `aaf` is `maf` where the annotation says alt
is the minor allele and `1 - maf` where it says alt is the major one -- the annotation decides a
bit, not a number, and its own AF never reaches the output.

Which annotation AF decides that bit is a measured choice, not a stylistic one. The annotation's
`AF` is Finnish, and this is a UK cohort; its `*_enrichment_nfe` columns are AF_fin / AF_nfe, so
dividing recovers a non-Finnish European frequency. Against the published UKBB_EUR MAF over the
280,058 annotated variants, orienting on the Finnish AF gives a median |dMAF| of 0.0285 and
disagrees by more than 0.1 on 4.7 % of them; the recovered NFE frequency gives 0.0097 and 0.03 %.
The count of variants where the decision both matters and is nearly a coin toss -- reference AF
within 0.05 of 0.5 while the published MAF is not -- falls from 11,316 to 4,276. So the NFE
frequency decides, with the Finnish AF as the fallback where neither enrichment column is usable.

A wrong orientation does not corrupt `maf` in `credible_sets_v`, which is LEAST(aaf, 1 - aaf) and
survives the flip, but it does invert the risk/protective call in `credible_set_stats.py`.
Variants absent from the annotation get a null `aaf`, and that classifier then falls back to the
sign of beta alone.

`cs_min_r2` has no source and is written as NA: the published table carries no within-set LD.
"""

import argparse
import io
import subprocess
import sys
from pathlib import Path

import numpy as np
import polars as pl
from scipy.special import log_ndtr

CHR_MAP = {"X": "23"}

ALPHA_COLUMNS = [f"alpha{i}" for i in range(1, 11)]

INPUT_COLUMNS = [
    "molecular_trait_id",
    "region",
    # the source's own chr_pos_ref_alt string, kept verbatim as the join key to the lead-variant
    # file: rebuilding it here would have to reproduce that file's chromosome spelling too
    "variant",
    "chromosome",
    "position",
    "ref",
    "alt",
    "maf",
    "cs_index",
    "pip",
    "z",
    *ALPHA_COLUMNS,
]

LEAD_COLUMNS = ["metabolite", "SNP", "ALL0", "ALL1", "BETA", "SE", "LOG10P", "Z", "MAF"]

ANNOTATION_COLUMNS = [
    "variant",
    "AF",
    "GENOME_enrichment_nfe",
    "EXOME_enrichment_nfe",
    "most_severe",
    "gene_most_severe",
]

OUTPUT_COLUMNS = [
    "dataset",
    "data_type",
    "trait",
    "trait_original",
    "cell_type",
    "chr",
    "pos",
    "ref",
    "alt",
    "mlog10p",
    "beta",
    "se",
    "pip",
    "cs_id",
    "cs_size",
    "cs_min_r2",
    "aaf",
    "most_severe",
    "gene_most_severe",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, help="Zenodo fine-mapping parquet")
    parser.add_argument(
        "--lead-variants",
        required=True,
        help="EUR_all_lead_variants.tsv.gz, the same cohort's published BETA/SE/LOG10P",
    )
    parser.add_argument("--annotation", required=True, help="Tabix indexed FinnGen annotated variants")
    parser.add_argument("--output-dir", default=".", help="Directory to write the unsorted TSV into")
    parser.add_argument("--dataset", default="UKBB_EUR_NMR_2026", help="Value of the dataset column")
    parser.add_argument("--cell-type", default="plasma", help="Value of the cell_type column")
    parser.add_argument(
        "--n-samples",
        type=int,
        default=413897,
        help="UKBB_EUR sample size the fine-mapping was run on, used to scale beta and se",
    )
    parser.add_argument(
        "--af-ambiguous",
        type=float,
        default=0.05,
        help="Report how many rows the annotation puts within this distance of AF 0.5 while the "
        "published MAF is not, i.e. where orienting aaf both matters and is least supported",
    )
    return parser.parse_args()


def _sci(expr: pl.Expr) -> pl.Expr:
    """Format a float expression as %.3e in one numpy call rather than per row, keeping nulls null."""
    formatted = expr.fill_null(0.0).map_batches(
        lambda s: pl.Series(np.char.mod("%.3e", s.to_numpy())), return_dtype=pl.String
    )
    return pl.when(expr.is_null()).then(None).otherwise(formatted)


def read_finemapping(path: str) -> pl.DataFrame:
    """Read the credible set variants, mapping X to 23 and picking each set's own alpha.

    cs_index names which of SuSiE's L single effects this credible set is ("L1".."L10"), so the
    variant's weight within it is the correspondingly numbered alpha column. Everything after
    this point uses that, not the file's overall `pip`.
    """
    df = pl.read_parquet(path, columns=INPUT_COLUMNS).rename({"variant": "variant_key"}).with_columns(
        pl.col("chromosome").replace(CHR_MAP).cast(pl.Int32, strict=False).alias("chr"),
        pl.col("position").cast(pl.Int64).alias("pos"),
        pl.col("cs_index").str.strip_prefix("L").cast(pl.Int32).alias("cs_number"),
    )

    n_bad_chr = df.filter(pl.col("chr").is_null()).height
    if n_bad_chr:
        sys.exit(f"{n_bad_chr} rows have a chromosome that is neither an autosome nor X")

    return df.with_columns(
        pl.concat_list(ALPHA_COLUMNS).list.get(pl.col("cs_number") - 1).alias("alpha_cs"),
        pl.concat_str("region", "cs_number", separator="_").alias("cs_id"),
        pl.concat_str("chr", "pos", "ref", "alt", separator=":").alias("variant_id"),
    ).drop(ALPHA_COLUMNS)


def annotation_header(path: str) -> list[str]:
    return (
        subprocess.run(["tabix", "-H", path], capture_output=True, text=True, check=True)
        .stdout.strip()
        .lstrip("#")
        .split("\t")
    )


def read_annotation(path: str, variants: pl.DataFrame) -> pl.DataFrame:
    """Fetch the annotation of the given variants with tabix.

    The regions are the credible set positions themselves, so single-base lookups through the
    index replace a full scan of a 700 MB file.
    """
    regions = "\n".join(
        f"{chrom}\t{pos - 1}\t{pos}"
        for chrom, pos in variants.select("chr", "pos").unique().sort("chr", "pos").iter_rows()
    )
    proc = subprocess.run(
        ["tabix", "-R", "-", path], input=regions, capture_output=True, text=True, check=True
    )
    anno = (
        pl.read_csv(
            io.BytesIO(proc.stdout.encode()),
            separator="\t",
            has_header=False,
            new_columns=annotation_header(path),
            null_values=["NA"],
            schema_overrides={
                "variant": pl.Utf8,
                "AF": pl.Float64,
                "GENOME_enrichment_nfe": pl.Float64,
                "EXOME_enrichment_nfe": pl.Float64,
            },
        )
        .select(ANNOTATION_COLUMNS)
        .rename({"variant": "variant_id"})
        .with_columns(
            # the enrichment columns are AF_fin / AF_nfe, so dividing recovers a non-Finnish
            # European frequency; a quotient outside (0, 1) is not one and falls back
            (pl.col("AF") / pl.coalesce("GENOME_enrichment_nfe", "EXOME_enrichment_nfe"))
            .pipe(lambda e: pl.when(e.is_between(0.0, 1.0, closed="none")).then(e).otherwise(None))
            .alias("AF_nfe")
        )
    )

    return anno.join(variants.select("variant_id"), on="variant_id", how="semi").unique(
        subset="variant_id"
    )


def read_lead_variants(path: str) -> pl.DataFrame:
    """Read the companion record's published statistics, keyed like the fine-mapping table."""
    return (
        pl.read_csv(path, separator="\t", columns=LEAD_COLUMNS, null_values=["NA"])
        .rename({"metabolite": "molecular_trait_id", "SNP": "variant_key"})
        .unique(subset=["molecular_trait_id", "variant_key"])
    )


def check_lead_agreement(df: pl.DataFrame) -> None:
    """Fail unless the two files describe the same variants with exactly opposed z-scores.

    This is the assumption the whole beta orientation rests on, so it is asserted on every run
    rather than trusted: if a future release of either file changes its convention, the negation
    below would silently invert every effect direction.
    """
    matched = df.filter(pl.col("BETA").is_not_null())
    if matched.is_empty():
        sys.exit("no fine-mapping row matched a published lead variant; check --lead-variants")

    mismatched_alleles = matched.filter(
        (pl.col("ALL0") != pl.col("ref")) | (pl.col("ALL1") != pl.col("alt"))
    ).height
    if mismatched_alleles:
        sys.exit(f"{mismatched_alleles} lead variants disagree with the fine-mapping on ref/alt")

    mismatched_maf = matched.filter((pl.col("MAF") - pl.col("maf")).abs() > 1e-6).height
    if mismatched_maf:
        sys.exit(
            f"{mismatched_maf} lead variants disagree with the fine-mapping on MAF; the two "
            "files are not describing the same cohort"
        )

    finite = matched.filter(pl.col("z").is_finite())
    worst = finite.select((pl.col("Z") + pl.col("z")).abs().max()).item()
    if worst is not None and worst > 0.01:
        sys.exit(
            f"the lead file's Z is not the negation of the fine-mapping's z "
            f"(max |Z + z| = {worst}); the effect orientation must be re-established"
        )
    print(
        f"  {matched.height} rows matched a published lead variant; "
        f"max |Z_lead + z| = {worst:.4g} over {finite.height} with a finite z"
    )


def derive_statistics(df: pl.DataFrame, n_samples: int) -> pl.DataFrame:
    """Fill mlog10p, beta and se: published where published, derived from z elsewhere.

    The derivation's per-trait scale is fitted against that trait's published standard errors,
    so `s_trait` absorbs whatever the inverse-normal transformed trait's residual variance
    really is after covariates rather than assuming 1.

    All three stay null where z overflowed and no published value exists: the source cannot say
    how large those |z| are, and the beta formula converges to a finite ceiling instead of
    failing, which would read as a measurement.
    """
    finite_z = pl.col("z").is_finite()
    # substituted before the arithmetic rather than masked after it: an infinite z turns the
    # root into inf and beta into a NaN, and numpy warns on the way there
    z = pl.when(finite_z).then(pl.col("z")).otherwise(0.0)
    root = (2.0 * pl.col("maf") * (1.0 - pl.col("maf")) * (n_samples + z**2)).sqrt()

    scale = (
        df.filter(pl.col("SE").is_not_null() & pl.col("z").is_finite())
        .with_columns((pl.col("SE") * root).alias("s"))
        .group_by("molecular_trait_id")
        .agg(pl.col("s").median().alias("s_trait"), pl.len().alias("n_lead"))
    )
    missing = df.select("molecular_trait_id").unique().join(scale, on="molecular_trait_id", how="anti")
    if missing.height:
        sys.exit(
            f"{missing.height} traits have no published lead variant to fit the effect scale "
            f"from: {sorted(missing['molecular_trait_id'])[:5]}"
        )
    print(
        f"  per-trait effect scale fitted from {scale['n_lead'].min()}-{scale['n_lead'].max()} "
        f"lead variants each, median {scale['s_trait'].median():.4f}, "
        f"range {scale['s_trait'].min():.4f}..{scale['s_trait'].max():.4f}"
    )

    df = df.join(scale.drop("n_lead"), on="molecular_trait_id", how="left")
    se = pl.col("s_trait") / root

    return df.with_columns(
        pl.coalesce(
            pl.col("LOG10P"),
            pl.when(finite_z).then(((-log_ndtr(-z.abs()) - np.log(2)) / np.log(10)).round(4)),
        ).alias("mlog10p"),
        pl.coalesce(pl.col("BETA"), pl.when(finite_z).then(-z * se)).alias("beta_num"),
        pl.coalesce(pl.col("SE"), pl.when(finite_z).then(se)).alias("se_num"),
    )


def orient_aaf(df: pl.DataFrame) -> pl.DataFrame:
    """Turn the published MAF into an alternative allele frequency.

    The annotation supplies only the answer to "is alt the minor allele here?"; the frequency
    itself stays the study's own. Variants the annotation does not carry keep a null aaf rather
    than borrowing a frequency, which is a different quantity.
    """
    reference_af = pl.coalesce("AF_nfe", "AF")
    return df.with_columns(
        pl.when(reference_af.is_null())
        .then(None)
        .when(reference_af < 0.5)
        .then(pl.col("maf"))
        .otherwise(1.0 - pl.col("maf"))
        .alias("aaf_num")
    )


def main() -> None:
    args = parse_args()

    print(f"Reading {args.input}...")
    cs = read_finemapping(args.input)
    print(
        f"  {cs.height} credible set variants in "
        f"{cs.select('molecular_trait_id', 'cs_id').n_unique()} sets "
        f"({cs['cs_id'].n_unique()} distinct cs_id values, which repeat across traits) "
        f"over {cs['molecular_trait_id'].n_unique()} traits and "
        f"{cs['variant_id'].n_unique()} distinct variants"
    )

    n_overflow = cs.filter(~pl.col("z").is_finite()).height
    print(f"  {n_overflow} rows have an overflowed z")

    print(f"Reading {args.lead_variants}...")
    lead = read_lead_variants(args.lead_variants)
    print(f"  {lead.height} published lead variants over {lead['molecular_trait_id'].n_unique()} traits")
    cs = cs.join(lead, on=["molecular_trait_id", "variant_key"], how="left")
    check_lead_agreement(cs)

    print(f"Reading annotation from {args.annotation}...")
    anno = read_annotation(args.annotation, cs.select("chr", "pos", "variant_id").unique())
    print(f"  {anno.height} of {cs['variant_id'].n_unique()} variants annotated")

    cs = cs.join(anno, on="variant_id", how="left")
    cs = derive_statistics(cs, args.n_samples)
    rescued = cs.filter(~pl.col("z").is_finite() & pl.col("beta_num").is_not_null()).height
    print(
        f"  {rescued} of the {n_overflow} overflowed rows take published statistics; "
        f"{n_overflow - rescued} keep a null mlog10p, beta and se"
    )
    cs = orient_aaf(cs)

    reference_af = pl.coalesce("AF_nfe", "AF")
    annotated = cs.filter(reference_af.is_not_null())
    disagreement = (
        annotated.select(
            (pl.min_horizontal(reference_af, 1.0 - reference_af) - pl.col("maf")).abs()
        )
        .to_series()
        .median()
    )
    print(
        f"  median |MAF(annotation) - MAF(published)| = {disagreement:.4f} over the "
        f"{annotated.height} rows the annotation orients"
    )
    n_ambiguous = annotated.filter(
        ((reference_af - 0.5).abs() < args.af_ambiguous)
        & ((pl.col("maf") - 0.5).abs() >= args.af_ambiguous)
    ).height
    print(
        f"  {n_ambiguous} of them have a reference AF within {args.af_ambiguous} of 0.5 while "
        "their published MAF is not, so orienting aaf both matters there and is least supported"
    )

    output_path = Path(args.output_dir) / f"{args.dataset}_cs_95.tsv"
    cs.with_columns(
        pl.lit(args.dataset).alias("dataset"),
        pl.lit("metaboQTL").alias("data_type"),
        pl.col("molecular_trait_id").alias("trait"),
        pl.col("molecular_trait_id").alias("trait_original"),
        pl.lit(args.cell_type).alias("cell_type"),
        _sci(pl.col("beta_num")).alias("beta"),
        _sci(pl.col("se_num")).alias("se"),
        _sci(pl.col("aaf_num")).alias("aaf"),
        pl.col("alpha_cs").round(4).alias("pip"),
        # region_cs is the house cs_id spelling and repeats across traits by design, so the
        # size has to be counted per trait as well -- over cs_id alone it merges 123,899
        # credible sets into 21,923
        pl.len().over("molecular_trait_id", "cs_id").cast(pl.Int32).alias("cs_size"),
        pl.lit(None, dtype=pl.Float64).alias("cs_min_r2"),
    ).select(OUTPUT_COLUMNS).write_csv(output_path, separator="\t", null_value="NA")
    print(f"wrote {output_path}")


if __name__ == "__main__":
    main()
