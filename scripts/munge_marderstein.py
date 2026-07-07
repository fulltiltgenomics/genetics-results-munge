#!/usr/bin/env python3
"""Munge the Marderstein/Kundu/Kundaje/Montgomery 2026 neuro-variants/FLARE release into the
canonical `open_chromatin` (Product A) and `variant_effect` (Product B) TSVs.

Source (Synapse project "neuro_variants" syn64693551; second release entity syn73770440;
       + public ENCODE adult-heart snATAC):
  - Paper: Marderstein, Kundu, Kundaje, Montgomery et al. 2026, Nature Genetics ("neuro-variants" /
    FLARE). scATAC peak atlases across 132 fetal+adult brain & heart cellular contexts, plus
    ChromBPNet and FLARE in-silico variant-effect predictions. Genome build hg38 (NO liftOver).
  - Synapse folder map (VERIFIED as folder ids; per-file layouts NOT verified — see TO-VERIFY):
      peaks            -> syn64716764   => Product A open_chromatin (peak calls per context)
      predictions      -> syn64713923   => Product B variant_effect (ChromBPNet variant scores)
      chrombpnet_models-> syn64713927   (models only; provenance, NOT loaded)
      flare            -> syn64717038   => Product B variant_effect (FLARE variant scores)
    Second release entity syn73770440 — enumerate its children at munge time (may add contexts).
  - ENCODE adult-heart snATAC is a separate PUBLIC accession (assay snATAC, not scATAC).

This ONE script emits THREE files, selected by --product:
  --product open_chromatin  -> marderstein_open_chromatin.tsv.gz   (18 open_chromatin cols; INTERVAL index)
  --product chrombpnet      -> marderstein_chrombpnet.tsv.gz       (18 variant_effect cols; POINT index)
  --product flare           -> marderstein_flare.tsv.gz            (18 variant_effect cols; POINT index)

>>> TO-VERIFY (format assumptions) <<<
No Synapse token is available in this environment, so the REAL per-file column layouts could not be
inspected. Every source column name below is a DEFAULT that the user MUST verify against the real
downloaded files and override with the corresponding --*-col flag (no code change needed for a
column-name difference). The synthetic --sample run exercises the transform end-to-end so the logic
is validated even though the real column names are unconfirmed.

INDEXING CONTRACT (from the epic design — the results-api overlap/point engines depend on it):
  - open_chromatin (Product A): INTERVAL-indexed  `tabix -p bed` (-s1 -b2 -e3, distinct start/end).
    Point-indexing would make the API's variant-overlap fast path SILENTLY MISS peaks whose
    interval contains pos.
  - variant_effect (Product B, chrombpnet + flare): POINT-indexed  `tabix -s1 -b2 -e2`.
  All files: `sort -k1,1 -k2,2n` before bgzip; chr-prefixed seqnames; exact canonical column order.
  The API prepends `resource` itself; it is NOT written into the file.

DATA ACCESS (documented; download is OFF by default and NOT run here):
  Fetching the Synapse folders requires synapseclient + a SYNAPSE_AUTH_TOKEN in the environment.
  This script assumes the input files are ALREADY downloaded (paths passed via flags). With
  --download (and only then) it would call synapseclient — see synapse_download() for the exact,
  documented commands. Never commit the token; put a fresh token in the runner env as a secret.

STAGING (off by default): with --stage the .tsv.gz + .tsv.gz.tbi are uploaded to BOTH GCS buckets
  (finngen-commons + daly-genetics-results). Without --stage nothing is uploaded.

Local validation without any real data / download / upload:
  python3 scripts/munge_marderstein.py --sample     # runs all three products on synthetic input
"""

import argparse
import re
import subprocess
import tempfile
from pathlib import Path

import polars as pl

RESOURCE = "marderstein"
VERSION = "2026"

# ---- Product A: open_chromatin -------------------------------------------------------------------
DATASET_OPEN = "marderstein_open_chromatin"
OC_COLUMNS = [
    "chrom", "start", "end", "peak_id", "dataset", "cell_type", "tissue", "life_stage",
    "condition", "assay", "score", "score_type", "n_cells", "cell_ontology_id",
    "uberon_id", "target_gene", "target_gene_id", "version",
]
OC_CONDITION = "unknown"                 # atlases default to unknown biological condition
OC_DEFAULT_ASSAY = "scATAC"              # Synapse contexts; ENCODE adult-heart overridden to snATAC

# ---- Product B: variant_effect (shared 18-col schema) --------------------------------------------
VE_COLUMNS = [
    "chrom", "pos", "ref", "alt", "variant", "rsid", "dataset", "model", "cell_type", "tissue",
    "life_stage", "score", "score_type", "mlog10p", "predicted_direction", "quantile_rank",
    "is_significant", "version",
]
DATASET_CHROMBPNET = "marderstein_chrombpnet"
MODEL_CHROMBPNET = "chrombpnet"
SCORE_TYPE_CHROMBPNET = "chrombpnet_logfc"
DATASET_FLARE = "marderstein_flare"
MODEL_FLARE = "flare"
SCORE_TYPE_FLARE = "flare_score"

# GCS destinations (documented; used only with --stage). Paths must match the results-api profile
# config files (open_chromatin/<resource>/<dataset>.tsv.gz and variant_effect/<resource>/...).
_GCS_FINNGEN = "gs://finngen-commons/results_api_data"
_GCS_DALY = "gs://daly-genetics-results"
GCS = {
    "open_chromatin": (
        f"{_GCS_FINNGEN}/open_chromatin/{RESOURCE}/{DATASET_OPEN}.tsv.gz",
        f"{_GCS_DALY}/open_chromatin/{RESOURCE}/{DATASET_OPEN}.tsv.gz",
    ),
    "chrombpnet": (
        f"{_GCS_FINNGEN}/variant_effect/{RESOURCE}/{DATASET_CHROMBPNET}.tsv.gz",
        f"{_GCS_DALY}/variant_effect/{RESOURCE}/{DATASET_CHROMBPNET}.tsv.gz",
    ),
    "flare": (
        f"{_GCS_FINNGEN}/variant_effect/{RESOURCE}/{DATASET_FLARE}.tsv.gz",
        f"{_GCS_DALY}/variant_effect/{RESOURCE}/{DATASET_FLARE}.tsv.gz",
    ),
}

# ------------------------------------------------------------------------------------------------
# Context-name parser: derive harmonized tissue {brain,heart} + life_stage {fetal,adult} from the
# 132 free-text context labels (30 fetal-brain / 24 adult-brain / 33 fetal-heart / 45 adult-heart).
# >>> TO-VERIFY: these keyword sets are a best-effort guess at the real label vocabulary. The
# authoritative path is the --context-map override (a TSV of cell_type -> tissue,life_stage[,assay]);
# build it from the release's context manifest and it takes precedence over the keyword parser. <<<
# ------------------------------------------------------------------------------------------------
_BRAIN_KEYWORDS = [
    "brain", "cortex", "cortical", "cerebr", "neuro", "neuron", "gaba", "glut", "astro",
    "oligo", "opc", "microglia", "hippocamp", "thalam", "striat", "interneuron", "npc",
    "radial glia", "rg ", "excitatory", "inhibitory",
]
_HEART_KEYWORDS = [
    "heart", "cardi", "cardiac", "cardiomyocyte", "ventric", "atri", "atrial", "myocard",
    "endocard", "pericard", "cm ", "_cm", "-cm", "smc", "pacemaker",
]
_FETAL_KEYWORDS = ["fetal", "fetus", "foetal", "prenatal", "embry", "developing", "gestation", "gw"]
_ADULT_KEYWORDS = ["adult", "postnatal", "aging", "aged", "mature"]


def _match_keyword(label: str, keywords: list[str]) -> bool:
    lo = label.lower()
    return any(kw in lo for kw in keywords)


def parse_context(label: str) -> tuple[str, str]:
    """Best-effort (tissue, life_stage) from a context label; 'unknown' when undetermined.

    tissue    in {brain, heart, unknown}
    life_stage in {fetal, adult, unknown}
    """
    if _match_keyword(label, _BRAIN_KEYWORDS):
        tissue = "brain"
    elif _match_keyword(label, _HEART_KEYWORDS):
        tissue = "heart"
    else:
        tissue = "unknown"

    if _match_keyword(label, _FETAL_KEYWORDS):
        life = "fetal"
    elif _match_keyword(label, _ADULT_KEYWORDS):
        life = "adult"
    else:
        life = "unknown"
    return tissue, life


def derive_context_map(contexts: list[str], args: argparse.Namespace) -> pl.DataFrame:
    """cell_type -> (tissue, life_stage, assay). Override map wins; else keyword parser; else default.

    assay defaults to args.assay (scATAC) and can be overridden per context by the map's assay
    column so the ENCODE adult-heart contexts are tagged snATAC.
    """
    map_tissue: dict[str, str] = {}
    map_life: dict[str, str] = {}
    map_assay: dict[str, str] = {}
    if args.context_map:
        m = pl.read_csv(args.context_map, separator="\t", infer_schema_length=10_000)
        has_assay = args.map_assay_col in m.columns
        for row in m.iter_rows(named=True):
            ct = str(row[args.map_context_col])
            if row.get(args.map_tissue_col) not in (None, ""):
                map_tissue[ct] = str(row[args.map_tissue_col])
            if row.get(args.map_life_col) not in (None, ""):
                map_life[ct] = str(row[args.map_life_col])
            if has_assay and row.get(args.map_assay_col) not in (None, ""):
                map_assay[ct] = str(row[args.map_assay_col])

    rows = []
    for ct in contexts:
        kw_tissue, kw_life = parse_context(ct)
        tissue = map_tissue.get(ct) or kw_tissue
        life = map_life.get(ct) or kw_life
        assay = map_assay.get(ct) or args.assay
        rows.append((ct, tissue, life, assay))
    return pl.DataFrame(
        rows,
        schema={"cell_type": pl.Utf8, "tissue": pl.Utf8, "life_stage": pl.Utf8, "assay": pl.Utf8},
        orient="row",
    )


# ------------------------------------------------------------------------------------------------
# Shared coordinate helpers
# ------------------------------------------------------------------------------------------------
# "chr1:1000-2000", "chr1-1000-2000", "chr1_1000_2000"; "chr" prefix optional (normalized on).
_INTERVAL_ID_RE = r"^(?:chr)?([^:_\-]+)[:_\-](\d+)[_\-](\d+)$"
# "chr1:12345:A:G" and separators - _ : | ; "chr" prefix optional.
_VARIANT_ID_RE = r"^(?:chr)?([^:_\-|]+)[:_\-|](\d+)[:_\-|]([A-Za-z]+)[:_\-|]([A-Za-z]+)$"


def _chr_prefixed(col: str) -> pl.Expr:
    c = pl.col(col).cast(pl.Utf8)
    return pl.when(c.str.starts_with("chr")).then(c).otherwise("chr" + c).alias(col)


def _coords_from_interval_id(df: pl.DataFrame, id_col: str) -> pl.DataFrame:
    g = df.select(pl.col(id_col).str.extract_groups(_INTERVAL_ID_RE).alias("_g")).unnest("_g")
    chrom = g["1"]
    chrom = pl.when(chrom.str.starts_with("chr")).then(chrom).otherwise("chr" + chrom)
    return df.with_columns(
        chrom.alias("chrom"),
        g["2"].cast(pl.Int64).alias("start"),
        g["3"].cast(pl.Int64).alias("end"),
    )


def _coords_from_variant_id(df: pl.DataFrame, var_col: str) -> pl.DataFrame:
    g = df.select(pl.col(var_col).str.extract_groups(_VARIANT_ID_RE).alias("_g")).unnest("_g")
    chrom = g["1"]
    chrom = pl.when(chrom.str.starts_with("chr")).then(chrom).otherwise("chr" + chrom)
    return df.with_columns(
        chrom.alias("chrom"),
        g["2"].cast(pl.Int64).alias("pos"),
        g["3"].str.to_uppercase().alias("ref"),
        g["4"].str.to_uppercase().alias("alt"),
    )


def _peak_key() -> pl.Expr:
    return pl.format("{}-{}-{}", pl.col("chrom"), pl.col("start"), pl.col("end"))


def _opt_expr(df: pl.DataFrame, col: str | None, dtype=pl.Utf8) -> pl.Expr:
    """Column expr if the (optional) column exists, else a typed null literal."""
    if col and col in df.columns:
        return pl.col(col).cast(dtype, strict=False)
    return pl.lit(None, dtype=dtype)


# ------------------------------------------------------------------------------------------------
# Product A: open_chromatin
# ------------------------------------------------------------------------------------------------
def load_peaks(args: argparse.Namespace) -> pl.DataFrame:
    """Read peak calls into long rows: one per (peak, context) with chrom/start/end/cell_type/score.

    Two accepted input shapes (TO-VERIFY which the release ships):
      --peaks       LONG TSV: one row per (peak, context); coords + a context column + optional score.
      --peak-matrix MATRIX TSV (catlas-style): one row per peak, one value column per context;
                    unpivoted here (kept where value > --min-score).
    """
    if args.peak_matrix:
        df = pl.read_csv(args.peak_matrix, separator="\t", infer_schema_length=10_000)
        if args.peak_id_col:
            df = _coords_from_interval_id(df, args.peak_id_col)
            coord_cols = [args.peak_id_col]
        else:
            df = df.rename({args.peak_chrom_col: "chrom", args.peak_start_col: "start", args.peak_end_col: "end"})
            df = df.with_columns(_chr_prefixed("chrom"), pl.col("start").cast(pl.Int64), pl.col("end").cast(pl.Int64))
            coord_cols = [args.peak_chrom_col, args.peak_start_col, args.peak_end_col]
        coord_cols = [c for c in coord_cols if c in df.columns]
        ctx_cols = [c for c in df.columns if c not in {"chrom", "start", "end", *coord_cols}]
        if not ctx_cols:
            raise ValueError("no per-context value columns found in --peak-matrix (check coordinate/id flags)")
        long = df.unpivot(index=["chrom", "start", "end"], on=ctx_cols, variable_name="cell_type", value_name="_value")
        long = long.with_columns(pl.col("_value").cast(pl.Float64, strict=False)).filter(pl.col("_value") > args.min_score)
        if args.score_type == "presence":
            long = long.with_columns(pl.lit(None, dtype=pl.Utf8).alias("score"))
        else:
            long = long.with_columns(pl.col("_value").alias("score"))
        return long.drop("_value")

    # LONG mode
    df = pl.read_csv(args.peaks, separator="\t", infer_schema_length=10_000)
    if args.peak_id_col:
        df = _coords_from_interval_id(df, args.peak_id_col)
    else:
        df = df.rename({args.peak_chrom_col: "chrom", args.peak_start_col: "start", args.peak_end_col: "end"})
        df = df.with_columns(_chr_prefixed("chrom"), pl.col("start").cast(pl.Int64), pl.col("end").cast(pl.Int64))
    df = df.rename({args.peak_context_col: "cell_type"})
    if args.score_type == "presence":
        df = df.with_columns(pl.lit(None, dtype=pl.Utf8).alias("score"))
    else:
        df = df.with_columns(_opt_expr(df, args.peak_score_col, pl.Float64).alias("score"))
        if args.peak_score_col in df.columns and args.min_score is not None:
            df = df.filter(pl.col("score") > args.min_score)
    return df.select("chrom", "start", "end", "cell_type", "score")


def build_open_chromatin(long: pl.DataFrame, args: argparse.Namespace) -> pl.DataFrame:
    df = long.with_columns(_peak_key().alias("peak_id"))
    ctx_map = derive_context_map(df["cell_type"].unique().to_list(), args)
    df = df.join(ctx_map, on="cell_type", how="left")
    df = df.with_columns(
        pl.lit(DATASET_OPEN).alias("dataset"),
        pl.lit(OC_CONDITION).alias("condition"),
        pl.lit(args.score_type).alias("score_type"),
        pl.lit(None, dtype=pl.Utf8).alias("n_cells"),
        pl.lit(None, dtype=pl.Utf8).alias("cell_ontology_id"),
        pl.lit(None, dtype=pl.Utf8).alias("uberon_id"),
        pl.lit(None, dtype=pl.Utf8).alias("target_gene"),
        pl.lit(None, dtype=pl.Utf8).alias("target_gene_id"),
        pl.lit(VERSION).alias("version"),
    )
    return df.select(OC_COLUMNS)


# ------------------------------------------------------------------------------------------------
# Product B: variant_effect (chrombpnet + flare)
# ------------------------------------------------------------------------------------------------
def _load_variant_coords(df: pl.DataFrame, args: argparse.Namespace) -> pl.DataFrame:
    """Normalize to chrom(chr-prefixed)/pos(Int)/ref/alt from either a single variant id or 4 cols."""
    if args.variant_col and args.variant_col in df.columns:
        return _coords_from_variant_id(df, args.variant_col)
    df = df.rename({args.ve_chrom_col: "chrom", args.ve_pos_col: "pos", args.ve_ref_col: "ref", args.ve_alt_col: "alt"})
    return df.with_columns(
        _chr_prefixed("chrom"),
        pl.col("pos").cast(pl.Int64),
        pl.col("ref").cast(pl.Utf8).str.to_uppercase(),
        pl.col("alt").cast(pl.Utf8).str.to_uppercase(),
    )


def _numeric_chrom_expr() -> pl.Expr:
    """Numeric chromosome token from the chr-prefixed chrom column: strip 'chr', then X=23, Y=24,
    M/MT=25, else the numeric chromosome as-is. Mirrors CHR_STRING_TO_INT_SQL / chrom_to_int() in
    genetics-results-db so the `variant` token matches the variant_effect table's chr INT64 encoding.
    """
    base = pl.col("chrom").cast(pl.Utf8).str.replace("(?i)^chr", "").str.to_uppercase()
    return (
        pl.when(base == "X").then(pl.lit("23"))
        .when(base == "Y").then(pl.lit("24"))
        .when(base == "M").then(pl.lit("25"))
        .when(base == "MT").then(pl.lit("25"))
        .otherwise(base)
    )


def _variant_string(numeric_chrom: bool) -> pl.Expr:
    # variant = "chr:pos:ref:alt". Canonical platform convention (default): numeric chromosome with
    # X=23/Y=24/M/MT=25 and NO "chr" prefix (e.g. "1:1000:A:G", "23:100:A:G"), matching the
    # variant_effect table's chr INT64 encoding. The file's `chrom` column stays chr-prefixed.
    chrom = _numeric_chrom_expr() if numeric_chrom else pl.col("chrom").cast(pl.Utf8)
    return pl.format("{}:{}:{}:{}", chrom, pl.col("pos"), pl.col("ref"), pl.col("alt"))


def _finalize_variant_effect(df: pl.DataFrame, args: argparse.Namespace) -> pl.DataFrame:
    """Add the `variant` column (numeric-chrom by default) and select the 18 canonical columns."""
    variant = _variant_string(numeric_chrom=not args.variant_keep_chr)
    df = df.with_columns(variant.alias("variant"))
    return df.select(VE_COLUMNS)


def _quantile_rank_over(value: pl.Expr, group: str | None) -> pl.Expr:
    """Fractional rank of abs(value) in (0,1]; per-group if group given, else global."""
    rank = value.abs().rank(method="average")
    n = pl.len()
    if group is not None:
        return (rank.over(group) / n.over(group))
    return rank / n


def build_chrombpnet(args: argparse.Namespace) -> pl.DataFrame:
    """ChromBPNet variant scores -> thresholded variant_effect rows (one per passing (variant,context)).

    Thresholding policy (task .2): emit a row only for (variant, context) pairs passing significance.
      1. If the release ships a per-context significance call (--cbp-sig-col), use it verbatim.
      2. Else is_significant = (mlog10p >= --mlog10p-thresh) OR (abs(logfc) in the top --top-quantile
         per context, i.e. quantile_rank >= 1 - top_quantile).
    is_significant AND mlog10p are ALWAYS set so the cutoff is transparent/re-filterable. A variant is
    retained if significant in ANY context (its passing contexts are emitted).
    """
    df = pl.read_csv(args.chrombpnet, separator="\t", infer_schema_length=10_000)
    df = _load_variant_coords(df, args)
    df = df.rename({args.cbp_context_col: "cell_type"})
    df = df.with_columns(pl.col(args.cbp_score_col).cast(pl.Float64, strict=False).alias("score"))

    # mlog10p: explicit column, else -log10(pval), else null
    if args.cbp_mlog10p_col and args.cbp_mlog10p_col in df.columns:
        df = df.with_columns(pl.col(args.cbp_mlog10p_col).cast(pl.Float64, strict=False).alias("mlog10p"))
    elif args.cbp_pval_col and args.cbp_pval_col in df.columns:
        df = df.with_columns((-1.0 * pl.col(args.cbp_pval_col).cast(pl.Float64, strict=False).log10()).alias("mlog10p"))
    else:
        df = df.with_columns(pl.lit(None, dtype=pl.Float64).alias("mlog10p"))

    # quantile_rank of abs(logfc) per context (from release column if provided, else computed)
    if args.cbp_quantile_col and args.cbp_quantile_col in df.columns:
        df = df.with_columns(pl.col(args.cbp_quantile_col).cast(pl.Float64, strict=False).alias("quantile_rank"))
    else:
        df = df.with_columns(_quantile_rank_over(pl.col("score"), "cell_type").alias("quantile_rank"))

    # predicted_direction: release column, else derived from the sign of the effect
    if args.cbp_direction_col and args.cbp_direction_col in df.columns:
        df = df.with_columns(pl.col(args.cbp_direction_col).cast(pl.Utf8).alias("predicted_direction"))
    else:
        df = df.with_columns(
            pl.when(pl.col("score") > 0).then(pl.lit("gain"))
            .when(pl.col("score") < 0).then(pl.lit("loss"))
            .otherwise(pl.lit(None, dtype=pl.Utf8)).alias("predicted_direction")
        )

    # is_significant: release call if present, else derived cutoff
    if args.cbp_sig_col and args.cbp_sig_col in df.columns:
        df = df.with_columns(pl.col(args.cbp_sig_col).cast(pl.Boolean, strict=False).alias("is_significant"))
    else:
        top_cut = 1.0 - args.top_quantile
        df = df.with_columns(
            (
                (pl.col("mlog10p").fill_null(float("-inf")) >= args.mlog10p_thresh)
                | (pl.col("quantile_rank") >= top_cut)
            ).alias("is_significant")
        )

    n_before = df.height
    df = df.filter(pl.col("is_significant"))
    print(f"  chrombpnet thresholding: kept {df.height}/{n_before} (variant,context) rows after significance filter")

    df = df.with_columns(
        _opt_expr(df, args.rsid_col).alias("rsid"),
        pl.lit(DATASET_CHROMBPNET).alias("dataset"),
        pl.lit(MODEL_CHROMBPNET).alias("model"),
        pl.lit(args.score_type_chrombpnet).alias("score_type"),
        pl.lit(VERSION).alias("version"),
    )
    # cell_type is verbatim; tissue/life_stage from the same context parser as Product A
    ctx_map = derive_context_map(df["cell_type"].unique().to_list(), args)
    df = df.join(ctx_map.select("cell_type", "tissue", "life_stage"), on="cell_type", how="left")
    return _finalize_variant_effect(df, args)


def build_flare(args: argparse.Namespace) -> pl.DataFrame:
    """FLARE variant scores -> pan-context variant_effect rows (cell_type/tissue/life_stage empty).

    FLARE is already a curated prioritized-variant set (far smaller): keep the FULL released set, NO
    thresholding. quantile_rank answers "how strongly": from the release column if present, else a
    global fractional rank of abs(score).
    """
    df = pl.read_csv(args.flare, separator="\t", infer_schema_length=10_000)
    df = _load_variant_coords(df, args)
    df = df.with_columns(pl.col(args.flare_score_col).cast(pl.Float64, strict=False).alias("score"))

    if args.flare_mlog10p_col and args.flare_mlog10p_col in df.columns:
        df = df.with_columns(pl.col(args.flare_mlog10p_col).cast(pl.Float64, strict=False).alias("mlog10p"))
    elif args.flare_pval_col and args.flare_pval_col in df.columns:
        df = df.with_columns((-1.0 * pl.col(args.flare_pval_col).cast(pl.Float64, strict=False).log10()).alias("mlog10p"))
    else:
        df = df.with_columns(pl.lit(None, dtype=pl.Float64).alias("mlog10p"))

    if args.flare_quantile_col and args.flare_quantile_col in df.columns:
        df = df.with_columns(pl.col(args.flare_quantile_col).cast(pl.Float64, strict=False).alias("quantile_rank"))
    else:
        df = df.with_columns(_quantile_rank_over(pl.col("score"), None).alias("quantile_rank"))

    df = df.with_columns(
        _opt_expr(df, args.flare_direction_col).alias("predicted_direction"),
        _opt_expr(df, args.flare_sig_col, pl.Boolean).alias("is_significant"),
        _opt_expr(df, args.rsid_col).alias("rsid"),
        # pan-context: no cell_type / tissue / life_stage
        pl.lit(None, dtype=pl.Utf8).alias("cell_type"),
        pl.lit(None, dtype=pl.Utf8).alias("tissue"),
        pl.lit(None, dtype=pl.Utf8).alias("life_stage"),
        pl.lit(DATASET_FLARE).alias("dataset"),
        pl.lit(MODEL_FLARE).alias("model"),
        pl.lit(SCORE_TYPE_FLARE).alias("score_type"),
        pl.lit(VERSION).alias("version"),
    )
    return _finalize_variant_effect(df, args)


# ------------------------------------------------------------------------------------------------
# Writers (distinct index modes) + staging
# ------------------------------------------------------------------------------------------------
def _shell_quote(s: str) -> str:
    return "'" + s.replace("'", "'\\''") + "'"


def _write_sorted_bgzip(df: pl.DataFrame, columns: list[str], output_path: str) -> str:
    """sort -k1,1 -k2,2n -> bgzip; header (first token '#'-prefixed) kept on top. Returns path."""
    tmpdir = tempfile.mkdtemp()
    body = Path(tmpdir) / "body.tsv"
    df.write_csv(body, separator="\t", include_header=False, null_value="")
    header = "#" + "\t".join(columns)
    pipeline = (
        f'( printf "%s\\n" {_shell_quote(header)}; '
        f'LC_ALL=C sort -k1,1 -k2,2n {_shell_quote(str(body))} ) | bgzip -c > {_shell_quote(output_path)}'
    )
    subprocess.run(pipeline, shell=True, check=True, executable="/bin/bash")
    return output_path


def write_interval(df: pl.DataFrame, output_path: str) -> None:
    """open_chromatin: INTERVAL index (tabix -p bed / -s1 -b2 -e3)."""
    _write_sorted_bgzip(df, OC_COLUMNS, output_path)
    subprocess.run(["tabix", "-f", "-p", "bed", output_path], check=True)
    print(f"  wrote {df.height} rows -> {output_path}")
    print(f"  indexed {output_path}.tbi (INTERVAL: tabix -p bed / -s1 -b2 -e3)")


def write_point(df: pl.DataFrame, output_path: str) -> None:
    """variant_effect: POINT index (tabix -s1 -b2 -e2)."""
    _write_sorted_bgzip(df, VE_COLUMNS, output_path)
    subprocess.run(["tabix", "-f", "-s", "1", "-b", "2", "-e", "2", output_path], check=True)
    print(f"  wrote {df.height} rows -> {output_path}")
    print(f"  indexed {output_path}.tbi (POINT: tabix -s1 -b2 -e2)")


def upload_to_gcs(local_path: str, gcs_path: str) -> None:
    subprocess.run(["gcloud", "storage", "cp", local_path, gcs_path], check=True)
    print(f"  uploaded {gcs_path}")
    tbi = local_path + ".tbi"
    if Path(tbi).exists():
        subprocess.run(["gcloud", "storage", "cp", tbi, gcs_path + ".tbi"], check=True)
        print(f"  uploaded {gcs_path}.tbi")


def synapse_download(args: argparse.Namespace) -> None:
    """Documented Synapse fetch — OFF by default and INTENTIONALLY not exercised here (no token).

    With --download and SYNAPSE_AUTH_TOKEN in the environment this would fetch the folders below.
    The exact folder for the current --product is the one this script actually consumes; the others
    are listed for provenance. Enumerate syn73770440's children too (second release entity).
    """
    folder = {
        "open_chromatin": "syn64716764",  # peaks
        "chrombpnet": "syn64713923",      # predictions (ChromBPNet variant scores)
        "flare": "syn64717038",           # FLARE variant scores
    }[args.product]
    print("=== Synapse download (documented) ===")
    print(f"  product {args.product} <- Synapse folder {folder}  (project syn64693551; also syn73770440)")
    print("  requires: pip install synapseclient ; export SYNAPSE_AUTH_TOKEN=<fresh token>")
    print("  equivalent CLI: synapse get -r " + folder + " --downloadLocation <DATA_DIR>")
    if not args.download:
        print("  --download NOT set: skipping actual download (this is the default).")
        return
    import os
    token = os.environ.get("SYNAPSE_AUTH_TOKEN")
    if not token:
        raise SystemExit("--download set but SYNAPSE_AUTH_TOKEN is not in the environment.")
    import synapseclient  # imported lazily so the script runs without the dep when not downloading
    syn = synapseclient.Synapse()
    syn.login(authToken=token)
    dest = args.download_dir or "."
    for child in syn.getChildren(folder):
        syn.get(child["id"], downloadLocation=dest)
    print(f"  downloaded folder {folder} -> {dest}")


# ------------------------------------------------------------------------------------------------
# CLI
# ------------------------------------------------------------------------------------------------
def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--product", choices=["open_chromatin", "chrombpnet", "flare"],
                   help="which of the three outputs to build")

    # Product A: peaks
    p.add_argument("--peaks", help="[open_chromatin] LONG peaks TSV (one row per peak x context)")
    p.add_argument("--peak-matrix", help="[open_chromatin] MATRIX peaks TSV (peak x context value columns)")
    p.add_argument("--peak-id-col", help="[open_chromatin] single interval-id column (e.g. chr1:1000-2000)")
    p.add_argument("--peak-chrom-col", default="chrom", help="[open_chromatin] chrom column (default: chrom)")
    p.add_argument("--peak-start-col", default="start", help="[open_chromatin] start column (default: start)")
    p.add_argument("--peak-end-col", default="end", help="[open_chromatin] end column (default: end)")
    p.add_argument("--peak-context-col", default="cell_type", help="[open_chromatin, LONG] context column (default: cell_type)")
    p.add_argument("--peak-score-col", default="score", help="[open_chromatin, LONG] score column (default: score)")
    p.add_argument("--score-type", default="signal", help="open_chromatin score_type token (default: signal; 'presence' emits empty score)")
    p.add_argument("--assay", default=OC_DEFAULT_ASSAY, help=f"open_chromatin default assay (default: {OC_DEFAULT_ASSAY}; ENCODE heart -> snATAC via context map)")
    p.add_argument("--min-score", type=float, default=0.0, help="open_chromatin: keep entries with score > min-score (default: 0.0)")

    # context override map (Product A + used for chrombpnet tissue/life)
    p.add_argument("--context-map", help="TSV overriding context -> tissue/life_stage[/assay]")
    p.add_argument("--map-context-col", default="cell_type", help="context map: context column (default: cell_type)")
    p.add_argument("--map-tissue-col", default="tissue", help="context map: tissue column (default: tissue)")
    p.add_argument("--map-life-col", default="life_stage", help="context map: life_stage column (default: life_stage)")
    p.add_argument("--map-assay-col", default="assay", help="context map: optional assay column (default: assay)")

    # Product B: shared variant coordinate flags
    p.add_argument("--variant-col", help="[variant_effect] single variant-id column (chr:pos:ref:alt)")
    p.add_argument("--ve-chrom-col", default="chrom", help="[variant_effect] chrom column (default: chrom)")
    p.add_argument("--ve-pos-col", default="pos", help="[variant_effect] pos column (default: pos)")
    p.add_argument("--ve-ref-col", default="ref", help="[variant_effect] ref column (default: ref)")
    p.add_argument("--ve-alt-col", default="alt", help="[variant_effect] alt column (default: alt)")
    p.add_argument("--rsid-col", default="rsid", help="[variant_effect] optional rsid column (default: rsid)")
    p.add_argument("--variant-keep-chr", action="store_true", help="escape hatch: keep the chr-prefixed `variant` token (e.g. chr1:1000:A:G) instead of the canonical numeric-chrom X=23 form")

    # Product B1: chrombpnet
    p.add_argument("--chrombpnet", help="[chrombpnet] LONG variant-effect TSV (one row per variant x context)")
    p.add_argument("--cbp-context-col", default="cell_type", help="[chrombpnet] context column (default: cell_type)")
    p.add_argument("--cbp-score-col", default="logfc", help="[chrombpnet] effect column (default: logfc)")
    p.add_argument("--cbp-mlog10p-col", default="mlog10p", help="[chrombpnet] optional mlog10p column (default: mlog10p)")
    p.add_argument("--cbp-pval-col", default="pval", help="[chrombpnet] optional p-value column (default: pval)")
    p.add_argument("--cbp-sig-col", default="is_significant", help="[chrombpnet] optional release significance call (default: is_significant)")
    p.add_argument("--cbp-direction-col", default="predicted_direction", help="[chrombpnet] optional direction column")
    p.add_argument("--cbp-quantile-col", default="quantile_rank", help="[chrombpnet] optional quantile column")
    p.add_argument("--score-type-chrombpnet", default=SCORE_TYPE_CHROMBPNET, help=f"[chrombpnet] score_type token (default: {SCORE_TYPE_CHROMBPNET})")
    p.add_argument("--mlog10p-thresh", type=float, default=3.0, help="[chrombpnet] mlog10p significance cutoff (default: 3.0 => p<=1e-3)")
    p.add_argument("--top-quantile", type=float, default=0.01, help="[chrombpnet] keep top fraction of abs(effect) per context (default: 0.01 => top 1%%)")

    # Product B2: flare
    p.add_argument("--flare", help="[flare] variant-effect TSV (pan-context)")
    p.add_argument("--flare-score-col", default="flare_score", help="[flare] score column (default: flare_score)")
    p.add_argument("--flare-mlog10p-col", default="mlog10p", help="[flare] optional mlog10p column")
    p.add_argument("--flare-pval-col", default="pval", help="[flare] optional p-value column")
    p.add_argument("--flare-sig-col", default="is_significant", help="[flare] optional significance column")
    p.add_argument("--flare-direction-col", default="predicted_direction", help="[flare] optional direction column")
    p.add_argument("--flare-quantile-col", default="quantile_rank", help="[flare] optional quantile column")

    # download (Synapse; OFF by default) + output + staging
    p.add_argument("--download", action="store_true", help="fetch the Synapse folder for --product (needs SYNAPSE_AUTH_TOKEN; OFF by default)")
    p.add_argument("--download-dir", help="download destination directory")
    p.add_argument("--output", help="output .tsv.gz path (default derived from --product)")
    p.add_argument("--stage", action="store_true", help="upload .tsv.gz + .tbi to BOTH GCS buckets (OFF by default)")
    p.add_argument("--gcs-finngen", help="override finngen GCS destination")
    p.add_argument("--gcs-daly", help="override daly GCS destination")

    p.add_argument("--sample", "--dry-run", dest="sample", action="store_true",
                   help="run all three products on synthetic input, validate both index modes; no download/upload")
    return p.parse_args()


def _default_output(product: str) -> str:
    return {
        "open_chromatin": f"./{DATASET_OPEN}.tsv.gz",
        "chrombpnet": f"./{DATASET_CHROMBPNET}.tsv.gz",
        "flare": f"./{DATASET_FLARE}.tsv.gz",
    }[product]


def run_product(args: argparse.Namespace) -> None:
    if args.download:
        synapse_download(args)

    if args.product == "open_chromatin":
        if not (args.peaks or args.peak_matrix):
            raise SystemExit("--product open_chromatin needs --peaks or --peak-matrix (or use --sample).")
        print(f"Building {DATASET_OPEN} ...")
        out = build_open_chromatin(load_peaks(args), args)
        assert out.columns == OC_COLUMNS, f"open_chromatin column order mismatch: {out.columns}"
        output_path = args.output or _default_output("open_chromatin")
        write_interval(out, output_path)
    elif args.product == "chrombpnet":
        if not args.chrombpnet:
            raise SystemExit("--product chrombpnet needs --chrombpnet (or use --sample).")
        print(f"Building {DATASET_CHROMBPNET} ...")
        out = build_chrombpnet(args)
        assert out.columns == VE_COLUMNS, f"variant_effect column order mismatch: {out.columns}"
        output_path = args.output or _default_output("chrombpnet")
        write_point(out, output_path)
    elif args.product == "flare":
        if not args.flare:
            raise SystemExit("--product flare needs --flare (or use --sample).")
        print(f"Building {DATASET_FLARE} ...")
        out = build_flare(args)
        assert out.columns == VE_COLUMNS, f"variant_effect column order mismatch: {out.columns}"
        output_path = args.output or _default_output("flare")
        write_point(out, output_path)
    else:
        raise SystemExit("--product is required (open_chromatin | chrombpnet | flare), or use --sample.")

    if args.stage:
        finngen, daly = GCS[args.product]
        print("Staging to GCS (both buckets) ...")
        upload_to_gcs(output_path, args.gcs_finngen or finngen)
        upload_to_gcs(output_path, args.gcs_daly or daly)
    else:
        print("  --stage not set: skipping GCS upload")
    print("Done.")


# ------------------------------------------------------------------------------------------------
# Sample / dry-run: synthetic inputs exercising all three products + both index modes
# ------------------------------------------------------------------------------------------------
def _sample_open_chromatin(tmpdir: Path) -> None:
    print("\n--- [1/3] open_chromatin (peaks -> INTERVAL index) ---")
    peaks = tmpdir / "peaks_long.tsv"
    # contexts span brain/heart x fetal/adult; the last is an ENCODE adult-heart context whose
    # assay (snATAC) comes from the override map, not the keyword parser.
    peaks.write_text(
        "chrom\tstart\tend\tcell_type\tscore\n"
        "chr1\t100100\t100600\tfetal_brain_GABAergic_neuron\t5.1\n"
        "chr1\t200200\t200700\tadult_brain_Astrocyte\t7.2\n"
        "chr2\t50050\t50550\tfetal_heart_Cardiomyocyte\t4.4\n"
        "chr2\t80080\t80580\tadult_heart_Ventricular_CM\t3.3\n"
        "chrX\t900900\t901400\tENCODE_adult_heart_snATAC_CM\t2.0\n"
    )
    ctx_map = tmpdir / "context_map.tsv"
    ctx_map.write_text(
        "cell_type\ttissue\tlife_stage\tassay\n"
        "ENCODE_adult_heart_snATAC_CM\theart\tadult\tsnATAC\n"
    )
    args = parse_args()
    args.product = "open_chromatin"
    args.peaks = str(peaks)
    args.peak_matrix = None
    args.peak_id_col = None
    args.peak_chrom_col, args.peak_start_col, args.peak_end_col = "chrom", "start", "end"
    args.peak_context_col, args.peak_score_col = "cell_type", "score"
    args.score_type, args.assay, args.min_score = "signal", "scATAC", 0.0
    args.context_map = str(ctx_map)
    args.map_context_col, args.map_tissue_col, args.map_life_col, args.map_assay_col = "cell_type", "tissue", "life_stage", "assay"

    out = build_open_chromatin(load_peaks(args), args)
    assert out.columns == OC_COLUMNS, out.columns
    tl = {r["cell_type"]: (r["tissue"], r["life_stage"], r["assay"]) for r in out.iter_rows(named=True)}
    assert tl["fetal_brain_GABAergic_neuron"] == ("brain", "fetal", "scATAC"), tl
    assert tl["adult_brain_Astrocyte"] == ("brain", "adult", "scATAC"), tl
    assert tl["fetal_heart_Cardiomyocyte"] == ("heart", "fetal", "scATAC"), tl
    assert tl["adult_heart_Ventricular_CM"] == ("heart", "adult", "scATAC"), tl
    assert tl["ENCODE_adult_heart_snATAC_CM"] == ("heart", "adult", "snATAC"), tl  # assay from map
    print(f"  context parsing (brain/heart x fetal/adult + ENCODE snATAC): OK")
    with pl.Config(tbl_cols=-1, tbl_width_chars=240, fmt_str_lengths=40):
        print(out.head())

    out_path = str(tmpdir / f"{DATASET_OPEN}.sample.tsv.gz")
    write_interval(out, out_path)
    # INTERVAL index: a position STRICTLY INSIDE a peak must be returned (point index would miss it)
    q = subprocess.run(["tabix", out_path, "chr1:200300-200400"], capture_output=True, text=True, check=True)
    assert q.stdout.strip(), "INTERVAL overlap query returned nothing — wrong index mode"
    print(f"  INTERVAL overlap query chr1:200300-200400 (inside chr1:200200-200700): OK")


def _sample_chrombpnet(tmpdir: Path) -> None:
    print("\n--- [2/3] chrombpnet (variant scores -> thresholded, POINT index) ---")
    cbp = tmpdir / "chrombpnet.tsv"
    # Focal variants v1 (chr1:1000 A/G) and v2 (chr1:2000 C/T), plus per-context background rows so
    # the top-quantile path has a real distribution (singleton contexts would trivially rank 1.0):
    #   v1 in fetal_brain -> sig via mlog10p; in adult_brain -> DROPPED (low mlog10p, bottom effect);
    #      in adult_heart -> sig via top-quantile abs(effect); v2 in fetal_brain -> sig via mlog10p;
    #      in adult_heart -> DROPPED (bottom effect, low mlog10p).
    rows = [
        # fetal_brain: two focal rows kept via mlog10p; background small effects
        ("chr1", 1000, "A", "G", "fetal_brain_GABAergic_neuron", 1.80, 6.0),
        ("chr1", 2000, "C", "T", "fetal_brain_GABAergic_neuron", 1.20, 4.0),
        ("chr1", 3100, "A", "G", "fetal_brain_GABAergic_neuron", 0.10, 0.5),
        ("chr1", 3200, "A", "G", "fetal_brain_GABAergic_neuron", 0.15, 0.4),
        ("chr1", 3300, "A", "G", "fetal_brain_GABAergic_neuron", 0.20, 0.3),
        ("chr1", 3400, "A", "G", "fetal_brain_GABAergic_neuron", 0.25, 0.2),
        # adult_brain: focal v1 is the bottom effect + low mlog10p -> DROPPED
        ("chr1", 1000, "A", "G", "adult_brain_Astrocyte", 0.02, 0.3),
        ("chr1", 4100, "A", "G", "adult_brain_Astrocyte", 0.80, 0.5),
        ("chr1", 4200, "A", "G", "adult_brain_Astrocyte", 0.90, 0.4),
        ("chr1", 4300, "A", "G", "adult_brain_Astrocyte", 1.00, 0.6),
        ("chr1", 4400, "A", "G", "adult_brain_Astrocyte", 1.10, 0.5),
        ("chr1", 4500, "A", "G", "adult_brain_Astrocyte", 1.20, 0.4),
        # adult_heart: focal v1 is the top effect -> KEPT via quantile; focal v2 is bottom -> DROPPED
        ("chr1", 1000, "A", "G", "adult_heart_Ventricular_CM", 2.50, 1.0),
        ("chr1", 2000, "C", "T", "adult_heart_Ventricular_CM", 0.01, 0.1),
        ("chr1", 5100, "A", "G", "adult_heart_Ventricular_CM", 0.50, 0.2),
        ("chr1", 5200, "A", "G", "adult_heart_Ventricular_CM", 0.60, 0.3),
        ("chr1", 5300, "A", "G", "adult_heart_Ventricular_CM", 0.70, 0.2),
        ("chr1", 5400, "A", "G", "adult_heart_Ventricular_CM", 0.80, 0.4),
    ]
    cbp.write_text(
        "chrom\tpos\tref\talt\tcell_type\tlogfc\tmlog10p\n"
        + "\n".join("\t".join(str(c) for c in r) for r in rows) + "\n"
    )
    args = parse_args()
    args.product = "chrombpnet"
    args.chrombpnet = str(cbp)
    args.variant_col = None
    args.ve_chrom_col, args.ve_pos_col, args.ve_ref_col, args.ve_alt_col = "chrom", "pos", "ref", "alt"
    args.rsid_col = "rsid"
    args.cbp_context_col, args.cbp_score_col = "cell_type", "logfc"
    args.cbp_mlog10p_col, args.cbp_pval_col = "mlog10p", "pval"
    args.cbp_sig_col = "is_significant"   # absent from synthetic file -> derived cutoff used
    args.cbp_direction_col, args.cbp_quantile_col = "predicted_direction", "quantile_rank"
    args.score_type_chrombpnet = SCORE_TYPE_CHROMBPNET
    args.mlog10p_thresh, args.top_quantile = 3.0, 0.2   # top 20% per context (real cutoff default: 0.01)
    args.context_map = None
    args.map_context_col, args.map_tissue_col, args.map_life_col, args.map_assay_col = "cell_type", "tissue", "life_stage", "assay"
    args.assay = OC_DEFAULT_ASSAY
    args.variant_keep_chr = False

    out = build_chrombpnet(args)
    assert out.columns == VE_COLUMNS, out.columns
    # `variant` is canonical numeric-chrom (no "chr" prefix); file's chrom column stays chr-prefixed
    assert set(out["variant"].to_list()) and all(not v.startswith("chr") for v in out["variant"].to_list()), out["variant"].to_list()
    assert all(c.startswith("chr") for c in out["chrom"].to_list()), out["chrom"].to_list()
    # sub-threshold (variant,context) rows dropped
    kept = {(r["variant"], r["cell_type"]) for r in out.iter_rows(named=True)}
    assert ("1:1000:A:G", "adult_brain_Astrocyte") not in kept, "sub-threshold row NOT dropped"
    assert ("1:2000:C:T", "adult_heart_Ventricular_CM") not in kept, "sub-threshold row NOT dropped"
    assert ("1:1000:A:G", "fetal_brain_GABAergic_neuron") in kept, kept
    assert ("1:1000:A:G", "adult_heart_Ventricular_CM") in kept, "top-effect row wrongly dropped"
    # is_significant + mlog10p always set; tissue/life parsed; direction derived
    for r in out.iter_rows(named=True):
        assert r["is_significant"] is True, r
        assert r["mlog10p"] is not None, r
        assert r["predicted_direction"] in ("gain", "loss"), r
    tissue_by_ct = {r["cell_type"]: (r["tissue"], r["life_stage"]) for r in out.iter_rows(named=True)}
    assert tissue_by_ct["fetal_brain_GABAergic_neuron"] == ("brain", "fetal"), tissue_by_ct
    assert tissue_by_ct["adult_heart_Ventricular_CM"] == ("heart", "adult"), tissue_by_ct
    print(f"  thresholding drops sub-threshold (variant,context) rows; is_significant/mlog10p set: OK")
    with pl.Config(tbl_cols=-1, tbl_width_chars=240, fmt_str_lengths=40):
        print(out.head())

    out_path = str(tmpdir / f"{DATASET_CHROMBPNET}.sample.tsv.gz")
    write_point(out, out_path)
    # POINT index: exact-position query
    q = subprocess.run(["tabix", out_path, "chr1:1000-1000"], capture_output=True, text=True, check=True)
    assert q.stdout.strip(), "POINT exact-pos query returned nothing — wrong index mode"
    assert "adult_brain_Astrocyte" not in q.stdout, "dropped row leaked into the file"
    print(f"  POINT exact-pos query chr1:1000-1000: OK")


def _sample_flare(tmpdir: Path) -> None:
    print("\n--- [3/3] flare (pan-context -> POINT index) ---")
    flare = tmpdir / "flare.tsv"
    flare.write_text(
        "chrom\tpos\tref\talt\tflare_score\trsid\n"
        "chr1\t1500\tA\tT\t0.92\trs111\n"
        "chr1\t3000\tG\tC\t0.55\trs222\n"
        "chr2\t7000\tC\tG\t0.10\trs333\n"
        "chrX\t100\tA\tT\t0.40\trs444\n"
    )
    args = parse_args()
    args.product = "flare"
    args.flare = str(flare)
    args.variant_col = None
    args.ve_chrom_col, args.ve_pos_col, args.ve_ref_col, args.ve_alt_col = "chrom", "pos", "ref", "alt"
    args.rsid_col = "rsid"
    args.flare_score_col = "flare_score"
    args.flare_mlog10p_col, args.flare_pval_col = "mlog10p", "pval"
    args.flare_sig_col = "is_significant"
    args.flare_direction_col, args.flare_quantile_col = "predicted_direction", "quantile_rank"
    args.variant_keep_chr = False

    out = build_flare(args)
    assert out.columns == VE_COLUMNS, out.columns
    # pan-context: cell_type / tissue / life_stage all empty (null)
    for r in out.iter_rows(named=True):
        assert r["cell_type"] is None and r["tissue"] is None and r["life_stage"] is None, r
        assert r["model"] == MODEL_FLARE and r["score_type"] == SCORE_TYPE_FLARE, r
    assert out["quantile_rank"].null_count() == 0, "quantile_rank should be computed globally"
    # canonical numeric-chrom `variant`: chr1 -> "1:...", chrX -> "23:..."; chrom column stays chr-prefixed
    var_by_chrom = {r["chrom"]: r["variant"] for r in out.iter_rows(named=True)}
    assert var_by_chrom["chr1"].startswith("1:"), var_by_chrom
    assert var_by_chrom["chrX"] == "23:100:A:T", var_by_chrom
    assert all(c.startswith("chr") for c in out["chrom"].to_list()), out["chrom"].to_list()
    print(f"  pan-context (cell_type/tissue/life_stage empty) + global quantile_rank + numeric-chrom variant (chrX->23): OK")
    with pl.Config(tbl_cols=-1, tbl_width_chars=240, fmt_str_lengths=40):
        print(out.head())

    out_path = str(tmpdir / f"{DATASET_FLARE}.sample.tsv.gz")
    write_point(out, out_path)
    q = subprocess.run(["tabix", out_path, "chr1:1500-1500"], capture_output=True, text=True, check=True)
    assert q.stdout.strip(), "POINT exact-pos query returned nothing — wrong index mode"
    print(f"  POINT exact-pos query chr1:1500-1500: OK")


def run_sample() -> None:
    print("=== SAMPLE / DRY-RUN: synthetic input, no Synapse download, no GCS upload ===")
    tmpdir = Path(tempfile.mkdtemp())
    _sample_open_chromatin(tmpdir)
    _sample_chrombpnet(tmpdir)
    _sample_flare(tmpdir)
    print("\n=== SAMPLE OK — open_chromatin=INTERVAL index; chrombpnet+flare=POINT index; "
          "no download / no upload performed ===")


def main() -> None:
    args = parse_args()
    if args.sample:
        run_sample()
        return
    if not args.product:
        raise SystemExit("--product is required (open_chromatin | chrombpnet | flare), or use --sample.")
    run_product(args)


if __name__ == "__main__":
    main()
