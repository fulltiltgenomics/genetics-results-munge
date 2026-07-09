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

REAL peaks layout (VERIFIED against Synapse syn64716764): 5 study SUBFOLDERS, each with per-context
ENCODE narrowPeak BEDs "<study>.<context>.overlap.peaks.bed.gz" (134 files, ~800 MB). narrowPeak is
10 headerless cols: chrom(chr) start(0-based) end name score(0-1000) strand signalValue pValue qValue
summit. Read directly with --peaks-dir (signalValue col 7 -> `score`, score_type=signal). cell_type =
"<study>.<context>" from the filename; the study prefix drives the baked-in DEFAULT_STUDY_CONTEXT map
(tissue/life_stage/assay; domcke tissue split from the filename). --context-map still overrides it.
The variant_effect (chrombpnet/flare) source column names remain DEFAULTS overridable via --*-col.

REAL chrombpnet layout (VERIFIED as a WIDE per-variant matrix, syn64713923): one row per variant, id in
`snp_id` (e.g. "chr19:58326671:C:A"), ~43 annotation columns + 132 context TRIPLETS of columns:
`abs_logfc.mean.<ctx>` (ABSOLUTE effect, >=0), `abs_logfc.mean.pval.<ctx>` (its p-value), and
`peak_overlap.<ctx>`; the context label is "<celltype>.<study>". These files are large (asd 160MB,
common 12GB, rare 16.6GB, ldsc 19GB compressed => >100GB uncompressed EACH). The --product chrombpnet
path STREAMS them in bounded row batches (--cbp-batch-rows): each batch is reshaped wide->long (melt the
132 triplets to one row per (variant,context) with abs_logfc+pval), mapped context->cell_type +
tissue/life_stage via the study suffix, thresholded single-pass on mlog10p = -log10(pval)
(>= --mlog10p-thresh, default 3), and APPENDED to an on-disk temp; peak RAM stays ~one batch, NOT the
100GB uncompress. Input files are processed ONE AT A TIME and (with --download) each raw file is deleted
before the next, so peak disk is bounded by the largest single input + the small filtered output. After
all files the temp is external-sorted (LC_ALL=C sort -k1,1 -k2,2n), bgzipped and POINT-indexed. A
legacy already-LONG single TSV is still accepted (auto-detected; falls back to the in-memory builder).

OUTPUT CONVENTIONS (applied to all three products):
  - chrom NUMERIC: the `chrom` column, `peak_id`, and the `variant` token strip any "chr" prefix and
    map X->23, Y->24, M/MT->25 (else the numeric contig). peak_id = "<numchrom>-<start>-<end>".
  - EMPTY -> "NA": every empty/missing output cell is written as the literal string "NA".

INDEXING CONTRACT (from the epic design — the results-api overlap/point engines depend on it):
  - open_chromatin (Product A): INTERVAL-indexed  `tabix -p bed` (-s1 -b2 -e3, distinct start/end).
    Point-indexing would make the API's variant-overlap fast path SILENTLY MISS peaks whose
    interval contains pos.
  - variant_effect (Product B, chrombpnet + flare): POINT-indexed  `tabix -s1 -b2 -e2`.
  All files: numeric-aware `sort -k1,1n -k2,2n` before bgzip (numeric seqnames group per tabix);
  exact canonical column order. The API prepends `resource` itself; it is NOT written into the file.

DATA ACCESS (download is OFF by default):
  synapseclient authenticates from ~/.synapseConfig (a SYNAPSE_AUTH_TOKEN env var, if set, wins).
  With --download the peaks folder is pulled RECURSIVELY (synapseutils.syncFromSynapse) to
  --download-dir, subfolders and all; --peaks-dir then reads those narrowPeak files. Without
  --download the script assumes inputs are already present. The token is never logged/printed.

STAGING (off by default): with --stage the .tsv.gz + .tsv.gz.tbi are uploaded to BOTH GCS buckets
  (finngen-commons + daly-genetics-results). Without --stage nothing is uploaded.

Local validation without any real data / download / upload:
  python3 scripts/munge_marderstein.py --sample     # runs all three products on synthetic input
"""

import argparse
import io
import itertools
import os
import re
import shutil
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
# the REAL Synapse ChromBPNet release scores `abs_logfc.mean.<ctx>` are ABSOLUTE (>=0), so the
# wide-streaming path tags them distinctly and leaves predicted_direction NA.
SCORE_TYPE_CHROMBPNET_ABS = "chrombpnet_abs_logfc"
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

# ------------------------------------------------------------------------------------------------
# DEFAULT study-level context map (VERIFIED against the real Synapse peaks layout syn64716764).
# The 5 study subfolders each carry per-context narrowPeak files named "<study>.<context>...";
# cell_type = "<study>.<context>" and the study prefix determines tissue/life_stage/assay. This is
# baked in so tissue/life_stage/assay resolve with NO manual --context-map (which still overrides).
# domcke_2020 mixes tissues, so its tissue is split from the filename (fetal_brain vs fetal_heart).
# ------------------------------------------------------------------------------------------------
DEFAULT_STUDY_CONTEXT: dict[str, tuple[str, str, str]] = {
    # study        -> (tissue, life_stage, assay)
    "ameen_2022":  ("heart", "fetal", "scATAC"),
    "corces_2020": ("brain", "adult", "scATAC"),
    "trevino_2021": ("brain", "fetal", "scATAC"),
    "encode_2024": ("heart", "adult", "snATAC"),
    # domcke_2020 handled specially below (life_stage/assay fixed, tissue from filename)
}


def _default_study_context(cell_type: str) -> tuple[str | None, str | None, str | None]:
    """(tissue, life_stage, assay) from the study prefix of a "<study>.<context>" cell_type.

    Returns None for any field that the default map cannot resolve so lower-priority resolvers
    (an explicit --context-map, then the keyword parser) can still fill it in.
    """
    study = cell_type.split(".", 1)[0].lower()
    if study == "domcke_2020":
        lo = cell_type.lower()
        tissue = "brain" if "fetal_brain" in lo else "heart" if "fetal_heart" in lo else None
        return tissue, "fetal", "scATAC"
    return DEFAULT_STUDY_CONTEXT.get(study, (None, None, None))


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
    """cell_type -> (tissue, life_stage, assay).

    Resolution precedence (highest first):
      1. explicit --context-map override (per-context TSV)
      2. baked-in DEFAULT_STUDY_CONTEXT (study prefix; domcke tissue split from filename)
      3. free-text keyword parser
      4. "unknown" (tissue/life_stage) / args.assay (assay) — never crashes
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
        d_tissue, d_life, d_assay = _default_study_context(ct)
        kw_tissue, kw_life = parse_context(ct)
        tissue = map_tissue.get(ct) or d_tissue or kw_tissue or "unknown"
        life = map_life.get(ct) or d_life or kw_life or "unknown"
        assay = map_assay.get(ct) or d_assay or args.assay
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
# narrowPeak files are named "<study>.<context>.overlap.peaks.bed.gz"; strip this to get cell_type.
_PEAKS_SUFFIX = ".overlap.peaks.bed.gz"


def load_peaks_dir(args: argparse.Namespace) -> pl.DataFrame:
    """Read the REAL Synapse peaks layout: per-context ENCODE narrowPeak BEDs read directly.

    Each "<study>.<context>.overlap.peaks.bed.gz" is a headerless 10-col narrowPeak:
      chrom(1, chr-prefixed) start(2, 0-based) end(3) name(4) score/0-1000(5) strand(6)
      signalValue(7) pValue(8) qValue(9) summit(10).
    We keep chrom/start/end, use signalValue (col 7) as `score` (score_type=signal), and derive
    cell_type = "<study>.<context>" from the filename (study = the "<study>." prefix). Files are
    found recursively so the 5 study subfolders produced by --download are all picked up.
    """
    root = Path(args.peaks_dir)
    if not root.exists():
        raise SystemExit(f"--peaks-dir {root} does not exist (run --download first, or fix the path)")
    files = sorted(root.rglob(args.peaks_glob))
    if not files:
        raise SystemExit(f"no files matching '{args.peaks_glob}' under {root}")
    print(f"  --peaks-dir: {len(files)} narrowPeak files under {root}")
    frames = []
    for fp in files:
        name = fp.name
        cell_type = name[:-len(_PEAKS_SUFFIX)] if name.endswith(_PEAKS_SUFFIX) else re.sub(r"\.bed(\.gz)?$", "", name)
        df = pl.read_csv(
            fp, separator="\t", has_header=False, infer_schema_length=10_000,
            columns=[0, 1, 2, 6], new_columns=["chrom", "start", "end", "score"],
        )
        df = df.with_columns(
            pl.col("start").cast(pl.Int64),
            pl.col("end").cast(pl.Int64),
            pl.col("score").cast(pl.Float64, strict=False),
            pl.lit(cell_type).alias("cell_type"),
        )
        if args.min_score is not None:
            df = df.filter(pl.col("score") > args.min_score)
        frames.append(df.select("chrom", "start", "end", "cell_type", "score"))
    return pl.concat(frames, how="vertical", rechunk=True)


def load_peaks(args: argparse.Namespace) -> pl.DataFrame:
    """Read peak calls into long rows: one per (peak, context) with chrom/start/end/cell_type/score.

    Three accepted input shapes:
      --peaks-dir   DIR of per-context narrowPeak "*.overlap.peaks.bed.gz" (the REAL Synapse layout).
      --peaks       LONG TSV: one row per (peak, context); coords + a context column + optional score.
      --peak-matrix MATRIX TSV (catlas-style): one row per peak, one value column per context;
                    unpivoted here (kept where value > --min-score).
    """
    if args.peaks_dir:
        return load_peaks_dir(args)
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
    # convention A: chrom + peak_id are NUMERIC (chrX=23, chrY=24, chrM/MT=25). load_peaks yields a
    # chr-prefixed chrom; normalize it first so peak_id = "<numchrom>-<start>-<end>" (e.g. 23-100-200).
    df = long.with_columns(_numeric_chrom_expr().alias("chrom"))
    # drop non-canonical hg38 contigs (alt/random/scaffold/Un) so the INT64 chr load never sees them
    df = df.filter(pl.col("chrom").is_in(CANONICAL_CHROMS))
    df = df.with_columns(_peak_key().alias("peak_id"))
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


def _numeric_chrom_expr(col: str = "chrom") -> pl.Expr:
    """Numeric chromosome token from a chrom column: strip 'chr', then X=23, Y=24, M/MT=25, else the
    numeric chromosome as-is. Mirrors CHR_STRING_TO_INT_SQL / chrom_to_int() in genetics-results-db so
    the output `chrom`/`peak_id`/`variant` tokens match the tables' chr INT64 encoding. Result is the
    numeric token as a string (kept as string so nulls / any unexpected contig survive untouched).
    """
    base = pl.col(col).cast(pl.Utf8).str.replace("(?i)^chr", "").str.to_uppercase()
    return (
        pl.when(base == "X").then(pl.lit("23"))
        .when(base == "Y").then(pl.lit("24"))
        .when(base == "M").then(pl.lit("25"))
        .when(base == "MT").then(pl.lit("25"))
        .otherwise(base)
    )


# canonical primary-assembly chromosomes as numeric strings (1..22, X=23, Y=24, M/MT=25). the
# platform is primary-assembly only and loads chrom as INT64, so non-canonical hg38 contigs
# (alt/random/scaffold/unplaced/Un_*) must be DROPPED from BOTH products or they break the chr load.
CANONICAL_CHROMS = frozenset(str(c) for c in range(1, 26))


def _variant_string(numeric_chrom: bool) -> pl.Expr:
    # variant = "chr:pos:ref:alt". Canonical platform convention (default): numeric chromosome with
    # X=23/Y=24/M/MT=25 and NO "chr" prefix (e.g. "1:1000:A:G", "23:100:A:G"), matching the
    # variant_effect table's chr INT64 encoding. The file's `chrom` column stays chr-prefixed.
    chrom = _numeric_chrom_expr() if numeric_chrom else pl.col("chrom").cast(pl.Utf8)
    return pl.format("{}:{}:{}:{}", chrom, pl.col("pos"), pl.col("ref"), pl.col("alt"))


def _finalize_variant_effect(df: pl.DataFrame, args: argparse.Namespace) -> pl.DataFrame:
    """Add the `variant` column and select the 18 canonical columns.

    convention A: the output `chrom` column is made NUMERIC (chrX=23 etc.), consistent with the
    numeric-chrom `variant` token, so tabix seqnames match the variant_effect table's chr encoding.
    """
    # drop non-canonical hg38 contigs (alt/random/scaffold/Un) regardless of --variant-keep-chr:
    # filter on the numeric mapping so the platform's INT64 chr load never sees a scaffold/alt name
    df = df.filter(_numeric_chrom_expr().is_in(CANONICAL_CHROMS))
    variant = _variant_string(numeric_chrom=not args.variant_keep_chr)
    df = df.with_columns(variant.alias("variant"))
    if not args.variant_keep_chr:
        df = df.with_columns(_numeric_chrom_expr().alias("chrom"))
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

    LEGACY long-input path: the real release is WIDE and handled by the streaming path
    (run_chrombpnet_streaming); this builder is retained for an already-LONG single TSV.
    """
    path = args.chrombpnet[0] if isinstance(args.chrombpnet, list) else args.chrombpnet
    df = pl.read_csv(path, separator="\t", infer_schema_length=10_000)
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


# ------------------------------------------------------------------------------------------------
# Product B1 (REAL data): ChromBPNet WIDE per-variant matrix -> STREAMING, memory-bounded reshape.
# The release files are >100GB uncompressed each, so we NEVER load a whole file: we read raw lines in
# bounded batches, melt each batch's 132 context triplets to long, single-pass threshold on mlog10p,
# and APPEND the passing rows to an on-disk temp. Files are processed one at a time (and deleted after
# download), then the temp is external-sorted + bgzipped + POINT-indexed.
# ------------------------------------------------------------------------------------------------
class _WideLayout:
    """Column layout of a wide ChromBPNet file: which columns are score/pval and their contexts."""

    def __init__(self, header_cols: list[str], args: argparse.Namespace):
        self.var_col = args.cbp_variant_col
        if self.var_col not in header_cols:
            raise SystemExit(
                f"--cbp-variant-col '{self.var_col}' not in the wide header (have e.g. {header_cols[:5]}...)"
            )
        sp = args.cbp_score_prefix                 # "abs_logfc.mean."
        pm = args.cbp_pval_infix.strip(".")        # ".pval." -> "pval"
        self.ctx_of_score: dict[str, str] = {}     # score column name -> context label
        self.ctx_of_pval: dict[str, str] = {}      # pval  column name -> context label
        for c in header_cols:
            if not c.startswith(sp):
                continue
            rem = c[len(sp):]
            if rem.startswith(pm + "."):
                self.ctx_of_pval[c] = rem[len(pm) + 1:]
            else:
                self.ctx_of_score[c] = rem
        score_ctx = set(self.ctx_of_score.values())
        if not score_ctx:
            raise SystemExit(
                f"no wide score columns start with --cbp-score-prefix '{sp}' (is this a wide file?)"
            )
        # only contexts with BOTH a score and a pval column are usable (need pval to threshold)
        pval_ctx = set(self.ctx_of_pval.values())
        self.contexts = sorted(score_ctx & pval_ctx)
        ctx_ok = set(self.contexts)
        self.score_cols = [c for c, ctx in self.ctx_of_score.items() if ctx in ctx_ok]
        self.pval_cols = [c for c, ctx in self.ctx_of_pval.items() if ctx in ctx_ok]
        self.needed = [self.var_col] + self.score_cols + self.pval_cols


def is_wide_chrombpnet_header(header_cols: list[str], args: argparse.Namespace) -> bool:
    return any(c.startswith(args.cbp_score_prefix) for c in header_cols)


def _line_reader(path: str):
    """Yield raw bytes lines from a possibly-gzipped file, STREAMING (bounded memory).

    gzip/bgzip is detected by magic bytes and decompressed via `gzip -dc` (bgzip output is
    gzip-compatible); plain files are read directly. Nothing is ever fully buffered.
    """
    with open(path, "rb") as fh:
        magic = fh.read(2)
    if magic == b"\x1f\x8b":
        proc = subprocess.Popen(["gzip", "-dc", str(path)], stdout=subprocess.PIPE, bufsize=1 << 20)
        stream, closer = proc.stdout, proc
    else:
        stream, closer = open(path, "rb"), None
    try:
        yield from stream
    finally:
        stream.close()
        if closer is not None:
            closer.wait()


def _read_header_cols(path: str, sep: str) -> list[str]:
    for line in _line_reader(path):
        return line.decode().rstrip("\r\n").split(sep)
    raise SystemExit(f"empty input file: {path}")


def _wide_batches(path: str, layout: _WideLayout, args: argparse.Namespace):
    """Yield polars DataFrames of ~--cbp-batch-rows variants, reading only the needed columns.

    We accumulate raw byte lines and hand each batch (header + lines) to polars, projecting to the
    needed columns with an explicit schema so peak RAM is bounded by one batch, not the whole file.
    """
    sep = args.cbp_sep
    overrides = {c: pl.Float64 for c in layout.needed if c != layout.var_col}
    overrides[layout.var_col] = pl.Utf8
    it = _line_reader(path)
    header = next(it)
    while True:
        lines = list(itertools.islice(it, args.cbp_batch_rows))
        if not lines:
            break
        buf = header + b"".join(lines)
        yield pl.read_csv(
            io.BytesIO(buf), separator=sep, columns=layout.needed, schema_overrides=overrides,
            infer_schema_length=0, null_values=["NA", "NaN", "nan", ""],
        )


def _wide_study_lookup(study: str, label: str) -> tuple[str | None, str | None]:
    """(tissue, life_stage) from the study SUFFIX of a "<celltype>.<study>" context label.

    Reuses the baked-in DEFAULT_STUDY_CONTEXT (domcke_2020 tissue split from the celltype text).
    """
    if study == "domcke_2020":
        lo = label.lower()
        tissue = "brain" if "brain" in lo else "heart" if "heart" in lo else None
        return tissue, "fetal"
    t = DEFAULT_STUDY_CONTEXT.get(study)
    return (t[0], t[1]) if t else (None, None)


def _wide_context_map(contexts: list[str]) -> pl.DataFrame:
    """context ("<celltype>.<study>") -> (cell_type verbatim, tissue, life_stage).

    Precedence: study-suffix DEFAULT_STUDY_CONTEXT, then the free-text keyword parser, then "unknown".
    """
    rows = []
    for ct in contexts:
        study = ct.rsplit(".", 1)[-1].lower() if "." in ct else ct.lower()
        d_tissue, d_life = _wide_study_lookup(study, ct)
        kw_tissue, kw_life = parse_context(ct)
        rows.append((ct, d_tissue or kw_tissue or "unknown", d_life or kw_life or "unknown"))
    return pl.DataFrame(
        rows, schema={"context": pl.Utf8, "tissue": pl.Utf8, "life_stage": pl.Utf8}, orient="row",
    )


def _reshape_wide_batch(
    df: pl.DataFrame, layout: _WideLayout, ctx_map: pl.DataFrame, args: argparse.Namespace,
) -> tuple[pl.DataFrame, int]:
    """One wide batch -> filtered 18-col variant_effect long rows. Returns (rows, total_pairs_seen).

    Streamable single-pass thresholding: melt the pval triplet, compute mlog10p = -log10(pval), keep
    only rows with mlog10p >= --mlog10p-thresh, THEN join the (much smaller) survivors to the melted
    scores — so we never materialize the full unfiltered variant x context product beyond the melt.
    """
    var_col = layout.var_col
    total_pairs = df.height * len(layout.contexts)

    plong = (
        df.select([var_col, *layout.pval_cols]).rename(layout.ctx_of_pval)
        .unpivot(index=var_col, variable_name="context", value_name="pval")
        .with_columns(
            pl.when(pl.col("pval") > 0).then(-pl.col("pval").log10())
            .when(pl.col("pval") == 0).then(pl.lit(float("inf")))
            .otherwise(pl.lit(None, dtype=pl.Float64)).alias("mlog10p")
        )
        .filter(pl.col("mlog10p") >= args.mlog10p_thresh)
    )
    if plong.height == 0:
        return plong.clear(), total_pairs

    slong = (
        df.select([var_col, *layout.score_cols]).rename(layout.ctx_of_score)
        .unpivot(index=var_col, variable_name="context", value_name="score")
    )
    long = plong.join(slong, on=[var_col, "context"], how="left")
    long = long.join(ctx_map, on="context", how="left").rename({"context": "cell_type"})
    long = _coords_from_variant_id(long, var_col)
    long = long.with_columns(
        pl.lit(None, dtype=pl.Utf8).alias("rsid"),
        pl.lit(DATASET_CHROMBPNET).alias("dataset"),
        pl.lit(MODEL_CHROMBPNET).alias("model"),
        pl.lit(SCORE_TYPE_CHROMBPNET_ABS).alias("score_type"),
        pl.lit(None, dtype=pl.Utf8).alias("predicted_direction"),   # abs_logfc has no direction
        pl.lit(None, dtype=pl.Float64).alias("quantile_rank"),      # not computable single-pass
        pl.lit(True).alias("is_significant"),
        pl.lit(VERSION).alias("version"),
    )
    return _finalize_variant_effect(long, args), total_pairs


def _append_ve_batch(df: pl.DataFrame, temp_fh) -> None:
    """Append 18-col rows to the temp (no header; empty cells coerced to null -> written as "NA")."""
    df = df.select(VE_COLUMNS).with_columns(
        pl.when(pl.col(c).cast(pl.Utf8).str.len_chars() == 0).then(None).otherwise(pl.col(c)).alias(c)
        for c in VE_COLUMNS
    )
    df.write_csv(temp_fh, separator="\t", include_header=False, null_value="NA")


def _iter_chrombpnet_inputs(args: argparse.Namespace):
    """Yield (path, delete_after) for each wide input, ONE AT A TIME.

    With --download each file is fetched from the Synapse predictions folder just before it is needed
    (and delete_after=True so the raw file is removed before the next), bounding peak disk usage.
    """
    if args.download:
        syn = _synapse_login()
        dest = args.download_dir or "."
        Path(dest).mkdir(parents=True, exist_ok=True)
        for child in _synapse_file_children("syn64713923", syn, args.cbp_name_filter):
            print(f"  downloading {child['name']} ({child['id']}) -> {dest}")
            ent = syn.get(child["id"], downloadLocation=dest)
            yield ent.path, True
    else:
        if not args.chrombpnet:
            raise SystemExit("--product chrombpnet needs --chrombpnet FILE(s) or --download (or --sample).")
        for p in args.chrombpnet:
            yield p, False


def run_chrombpnet_streaming(args: argparse.Namespace, output_path: str) -> None:
    """Full wide->long ChromBPNet build in ONE pass over bounded batches, per input file."""
    tmpdir = args.cbp_tmpdir or os.path.dirname(os.path.abspath(output_path)) or "."
    Path(tmpdir).mkdir(parents=True, exist_ok=True)
    temp_path = os.path.join(tmpdir, f"{DATASET_CHROMBPNET}.long.tmp.tsv")
    ctx_cache: dict[tuple[str, ...], pl.DataFrame] = {}
    kept = total = 0
    print(f"Building {DATASET_CHROMBPNET} (WIDE streaming; batch={args.cbp_batch_rows} variants) ...")
    with open(temp_path, "wb") as temp_fh:
        for path, delete_after in _iter_chrombpnet_inputs(args):
            layout = _WideLayout(_read_header_cols(path, args.cbp_sep), args)
            key = tuple(layout.contexts)
            if key not in ctx_cache:
                ctx_cache[key] = _wide_context_map(layout.contexts)
            ctx_map = ctx_cache[key]
            print(f"  processing {path}: {len(layout.contexts)} contexts")
            f_kept = f_total = 0
            for batch in _wide_batches(path, layout, args):
                out, pairs = _reshape_wide_batch(batch, layout, ctx_map, args)
                f_total += pairs
                f_kept += out.height
                if out.height:
                    _append_ve_batch(out, temp_fh)
            kept += f_kept
            total += f_total
            print(f"    kept {f_kept}/{f_total} (variant,context) pairs")
            if delete_after:
                os.remove(path)
                print(f"    deleted raw input {path}")
    print(f"  chrombpnet thresholding TOTAL: kept {kept}/{total} (variant,context) pairs "
          f"(mlog10p >= {args.mlog10p_thresh})")
    _finalize_point_external_sort(temp_path, VE_COLUMNS, output_path, tmpdir)
    os.remove(temp_path)


def _finalize_point_external_sort(temp_path: str, columns: list[str], output_path: str, tmpdir: str) -> None:
    """External `LC_ALL=C sort -k1,1 -k2,2n` (numeric pos, seqnames grouped) -> bgzip -> POINT index.

    Sorting on disk (not in RAM) keeps the finalize step memory-bounded regardless of the temp size.
    """
    sorted_path = temp_path + ".sorted"
    env = {**os.environ, "LC_ALL": "C"}
    with open(sorted_path, "wb") as out:
        subprocess.run(["sort", "-T", tmpdir, "-k1,1", "-k2,2n", temp_path],
                       env=env, stdout=out, check=True)
    header = ("#" + "\t".join(columns) + "\n").encode()
    with open(output_path, "wb") as out_fh:
        proc = subprocess.Popen(["bgzip", "-c"], stdin=subprocess.PIPE, stdout=out_fh)
        assert proc.stdin is not None
        proc.stdin.write(header)
        with open(sorted_path, "rb") as sf:
            shutil.copyfileobj(sf, proc.stdin, length=1 << 20)
        proc.stdin.close()
        rc = proc.wait()
    if rc != 0:
        raise subprocess.CalledProcessError(rc, "bgzip -c")
    os.remove(sorted_path)
    subprocess.run(["tabix", "-f", "-s", "1", "-b", "2", "-e", "2", output_path], check=True)
    print(f"  wrote {output_path}")
    print(f"  indexed {output_path}.tbi (POINT: tabix -s1 -b2 -e2)")


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
def _write_sorted_bgzip(df: pl.DataFrame, columns: list[str], output_path: str) -> str:
    """numeric-aware sort -> bgzip; header (first token '#'-prefixed) kept on top. Returns path.

    convention B: every empty/missing cell is written as the literal "NA" (nulls via null_value, and
    any residual empty strings coerced to null first) so no output cell is ever the empty string.

    Sorting is done IN-MEMORY by (numeric chrom, position) — chrom is the numeric token (chrX=23) so
    all records of a seqname group contiguously and positions ascend within, exactly as tabix needs.
    The sorted rows are streamed straight into `bgzip` (no multi-GB uncompressed temp file and no
    external `sort` disk-spill), which matters for large atlases on tight disks.
    """
    df = df.with_columns(
        pl.when(pl.col(c).cast(pl.Utf8).str.len_chars() == 0)
        .then(None)
        .otherwise(pl.col(c))
        .alias(c)
        for c in df.columns
    )
    # columns[1] is the position column (start for open_chromatin, pos for variant_effect).
    df = df.sort(
        by=[pl.col("chrom").cast(pl.Int64, strict=False), pl.col(columns[1]).cast(pl.Int64, strict=False)],
        nulls_last=True,
    )
    header = ("#" + "\t".join(columns) + "\n").encode()
    with open(output_path, "wb") as out_fh:
        proc = subprocess.Popen(["bgzip", "-c"], stdin=subprocess.PIPE, stdout=out_fh)
        assert proc.stdin is not None
        proc.stdin.write(header)
        df.write_csv(proc.stdin, separator="\t", include_header=False, null_value="NA")
        proc.stdin.close()
        rc = proc.wait()
    if rc != 0:
        raise subprocess.CalledProcessError(rc, "bgzip -c")
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


def _synapse_login():
    """Authenticate synapseclient. SYNAPSE_AUTH_TOKEN wins, else ~/.synapseConfig. Never prints it."""
    import synapseclient   # imported lazily so the script runs without the dep when not downloading

    syn = synapseclient.Synapse()
    token = os.environ.get("SYNAPSE_AUTH_TOKEN")
    if token:
        syn.login(authToken=token)            # token never logged/printed
    else:
        syn.login()                           # reads ~/.synapseConfig
    return syn


def _synapse_file_children(folder_id: str, syn, name_filter: str | None):
    """List FILE entities under a Synapse folder (recursing into subfolders), optional name regex.

    Used by the chrombpnet streaming path to fetch wide files one at a time; --cbp-name-filter narrows
    the set (e.g. 'asd' for the bounded real test) so we never pull all 52GB at once.
    """
    rx = re.compile(name_filter, re.I) if name_filter else None
    out = []
    for child in syn.getChildren(folder_id, includeTypes=["file", "folder"]):
        if child.get("type", "").endswith("Folder") or "Folder" in child.get("concreteType", ""):
            out.extend(_synapse_file_children(child["id"], syn, name_filter))
        elif rx is None or rx.search(child["name"]):
            out.append(child)
    return out


def synapse_download(args: argparse.Namespace) -> None:
    """Fetch the Synapse folder for --product RECURSIVELY (peaks has per-study SUBFOLDERS).

    Auth: synapseclient reads ~/.synapseConfig by default; a SYNAPSE_AUTH_TOKEN in the environment
    takes precedence if set. The token is NEVER printed. The peaks folder syn64716764 holds 5 study
    subfolders of per-context narrowPeak files, so we use synapseutils.syncFromSynapse to walk the
    whole tree (a flat getChildren() would miss the subfolders).
    """
    folder = {
        "open_chromatin": "syn64716764",  # peaks (5 study subfolders of narrowPeak BEDs)
        "chrombpnet": "syn64713923",      # predictions (ChromBPNet variant scores)
        "flare": "syn64717038",           # FLARE variant scores
    }[args.product]
    dest = args.download_dir or "."
    print("=== Synapse download (recursive) ===")
    print(f"  product {args.product} <- Synapse folder {folder}  (project syn64693551)")
    if not args.download:
        print("  --download NOT set: skipping actual download (this is the default).")
        return

    import time

    import synapseutils

    Path(dest).mkdir(parents=True, exist_ok=True)
    syn = _synapse_login()
    t0 = time.time()
    results = synapseutils.syncFromSynapse(syn, folder, path=dest)
    elapsed = time.time() - t0
    total_bytes = sum(f.stat().st_size for f in Path(dest).rglob("*") if f.is_file())
    print(f"  downloaded {len(results)} files ({total_bytes / 1e6:.1f} MB) -> {dest} in {elapsed:.0f}s")


# ------------------------------------------------------------------------------------------------
# CLI
# ------------------------------------------------------------------------------------------------
def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--product", choices=["open_chromatin", "chrombpnet", "flare"],
                   help="which of the three outputs to build")

    # Product A: peaks
    p.add_argument("--peaks-dir", help="[open_chromatin] DIR of per-context narrowPeak '*.overlap.peaks.bed.gz' (REAL Synapse layout; searched recursively)")
    p.add_argument("--peaks-glob", default="*.overlap.peaks.bed.gz", help="[open_chromatin] glob for --peaks-dir narrowPeak files (default: *.overlap.peaks.bed.gz)")
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
    p.add_argument("--chrombpnet", nargs="+", help="[chrombpnet] wide per-variant matrix FILE(s), processed one at a time (auto-detected; a legacy LONG TSV also works)")
    p.add_argument("--cbp-wide", action="store_true", help="[chrombpnet] force the WIDE streaming reshape (else auto-detected from the header)")
    p.add_argument("--cbp-variant-col", default="snp_id", help="[chrombpnet, wide] variant-id column (default: snp_id, e.g. chr19:58326671:C:A)")
    p.add_argument("--cbp-score-prefix", default="abs_logfc.mean.", help="[chrombpnet, wide] per-context score column prefix (default: 'abs_logfc.mean.'); context = the remainder")
    p.add_argument("--cbp-pval-infix", default=".pval.", help="[chrombpnet, wide] token marking a p-value column right after the score prefix (default: '.pval.')")
    p.add_argument("--cbp-batch-rows", type=int, default=100_000, help="[chrombpnet, wide] variants per streamed batch (default: 100000; bounds peak RAM)")
    p.add_argument("--cbp-sep", default="\t", help="[chrombpnet, wide] field separator of the wide file (default: tab)")
    p.add_argument("--cbp-name-filter", help="[chrombpnet, --download] regex to select which predictions files to fetch (e.g. 'asd')")
    p.add_argument("--cbp-tmpdir", help="[chrombpnet, wide] dir for the on-disk temp + external sort (default: alongside --output)")
    p.add_argument("--cbp-context-col", default="cell_type", help="[chrombpnet, LONG] context column (default: cell_type)")
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


def _use_wide_chrombpnet(args: argparse.Namespace) -> bool:
    """Wide streaming vs legacy long: forced (--cbp-wide), implied by --download (real files are wide),
    else auto-detected from the first local file's header."""
    if args.cbp_wide or args.download:
        return True
    if not args.chrombpnet:
        return False
    return is_wide_chrombpnet_header(_read_header_cols(args.chrombpnet[0], args.cbp_sep), args)


def run_product(args: argparse.Namespace) -> None:
    # chrombpnet-wide manages its own per-file Synapse download; only the other products (and the
    # legacy long path) use the whole-folder recursive sync.
    chrombpnet_wide = args.product == "chrombpnet" and _use_wide_chrombpnet(args)
    if args.download and not chrombpnet_wide:
        synapse_download(args)

    if args.product == "open_chromatin":
        if not (args.peaks_dir or args.peaks or args.peak_matrix):
            raise SystemExit("--product open_chromatin needs --peaks-dir, --peaks or --peak-matrix (or use --sample).")
        print(f"Building {DATASET_OPEN} ...")
        out = build_open_chromatin(load_peaks(args), args)
        assert out.columns == OC_COLUMNS, f"open_chromatin column order mismatch: {out.columns}"
        output_path = args.output or _default_output("open_chromatin")
        write_interval(out, output_path)
    elif args.product == "chrombpnet":
        output_path = args.output or _default_output("chrombpnet")
        if chrombpnet_wide:
            run_chrombpnet_streaming(args, output_path)   # streams wide files; POINT-indexes the output
        else:
            if not args.chrombpnet:
                raise SystemExit("--product chrombpnet needs --chrombpnet (or use --sample).")
            print(f"Building {DATASET_CHROMBPNET} (legacy LONG) ...")
            out = build_chrombpnet(args)
            assert out.columns == VE_COLUMNS, f"variant_effect column order mismatch: {out.columns}"
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
        # non-canonical hg38 contig (alt/random/scaffold) — MUST be dropped by the canonical filter
        "chr1_KI270706v1_random\t300300\t300800\tadult_brain_Astrocyte\t9.0\n"
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
    # convention A: chrom + peak_id are NUMERIC (chrX -> 23), no "chr" prefix
    assert all(not c.startswith("chr") for c in out["chrom"].to_list()), out["chrom"].to_list()
    pk = {r["chrom"]: r["peak_id"] for r in out.iter_rows(named=True)}
    assert pk["23"] == "23-900900-901400", pk  # chrX row -> numeric peak_id
    assert out.filter(pl.col("chrom") == "23").height == 1, "chrX should map to numeric chrom 23"
    print(f"  numeric chrom + peak_id (chrX->23, peak_id 23-900900-901400): OK")
    # non-canonical hg38 contigs (alt/random/scaffold/Un) are DROPPED (Fix A) or they break the
    # INT64 chr load; the "chr1_KI270706v1_random" synthetic peak must not survive
    assert not any("KI270706" in p for p in out["peak_id"].to_list()), "non-canonical contig leaked"
    assert set(out["chrom"].to_list()) == {"1", "2", "23"}, out["chrom"].unique().to_list()
    print(f"  non-canonical contig (chr1_KI270706v1_random) dropped: OK")
    with pl.Config(tbl_cols=-1, tbl_width_chars=240, fmt_str_lengths=40):
        print(out.head())

    out_path = str(tmpdir / f"{DATASET_OPEN}.sample.tsv.gz")
    write_interval(out, out_path)
    # INTERVAL index: a position STRICTLY INSIDE a peak must be returned (point index would miss it)
    q = subprocess.run(["tabix", out_path, "1:200300-200400"], capture_output=True, text=True, check=True)
    assert q.stdout.strip(), "INTERVAL overlap query returned nothing — wrong index mode"
    print(f"  INTERVAL overlap query 1:200300-200400 (inside 1:200200-200700): OK")


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
    # convention A: chrom column is numeric too, consistent with the numeric-chrom variant token
    assert all(not c.startswith("chr") for c in out["chrom"].to_list()), out["chrom"].to_list()
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
    # POINT index: exact-position query (numeric seqname)
    q = subprocess.run(["tabix", out_path, "1:1000-1000"], capture_output=True, text=True, check=True)
    assert q.stdout.strip(), "POINT exact-pos query returned nothing — wrong index mode"
    assert "adult_brain_Astrocyte" not in q.stdout, "dropped row leaked into the file"
    print(f"  POINT exact-pos query 1:1000-1000: OK")


def _sample_flare(tmpdir: Path) -> None:
    print("\n--- [3/3] flare (pan-context -> POINT index) ---")
    flare = tmpdir / "flare.tsv"
    flare.write_text(
        "chrom\tpos\tref\talt\tflare_score\trsid\n"
        "chr1\t1500\tA\tT\t0.92\trs111\n"
        "chr1\t3000\tG\tC\t0.55\trs222\n"
        "chr2\t7000\tC\tG\t0.10\trs333\n"
        "chrX\t100\tA\tT\t0.40\trs444\n"
        # non-canonical hg38 contig — MUST be dropped from the variant_effect output too (Fix A)
        "chr1_KI270706v1_random\t500\tA\tT\t0.30\trs555\n"
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
    # convention A: numeric-chrom `variant` AND numeric `chrom` column: chr1 -> "1", chrX -> "23"
    var_by_chrom = {r["chrom"]: r["variant"] for r in out.iter_rows(named=True)}
    assert var_by_chrom["1"].startswith("1:"), var_by_chrom
    assert var_by_chrom["23"] == "23:100:A:T", var_by_chrom
    assert all(not c.startswith("chr") for c in out["chrom"].to_list()), out["chrom"].to_list()
    # non-canonical hg38 contigs are DROPPED from variant_effect too (Fix A): rs555 must not survive
    assert set(out["chrom"].to_list()) == {"1", "2", "23"}, out["chrom"].unique().to_list()
    assert "rs555" not in out["rsid"].to_list(), "non-canonical contig variant leaked"
    print(f"  pan-context (cell_type/tissue/life_stage empty) + global quantile_rank + numeric chrom (chrX->23) + non-canonical dropped: OK")
    with pl.Config(tbl_cols=-1, tbl_width_chars=240, fmt_str_lengths=40):
        print(out.head())

    out_path = str(tmpdir / f"{DATASET_FLARE}.sample.tsv.gz")
    write_point(out, out_path)
    q = subprocess.run(["tabix", out_path, "1:1500-1500"], capture_output=True, text=True, check=True)
    assert q.stdout.strip(), "POINT exact-pos query returned nothing — wrong index mode"
    print(f"  POINT exact-pos query 1:1500-1500: OK")


def _sample_chrombpnet_wide(tmpdir: Path) -> None:
    print("\n--- [2b/3] chrombpnet WIDE streaming (wide->long reshape + threshold, POINT index) ---")
    import gzip

    # two contexts x 3 columns each (score / score.pval / peak_overlap); context = "<celltype>.<study>".
    # trevino_2021 -> brain/fetal, encode_2024 -> heart/adult (from the study suffix). Rows:
    #   chr1:1000:A:G significant in ctx1 (pval 1e-6), SUB-THRESHOLD in ctx2 (pval 0.5 -> dropped);
    #   chrX:200:C:T  significant in ctx1 (pval 1e-5), sub-threshold in ctx2 (pval 0.9 -> dropped).
    ctx1, ctx2 = "GABAergic_neuron.trevino_2021", "Ventricular_CM.encode_2024"
    wide = tmpdir / "chrombpnet_wide.tsv.gz"
    header = "\t".join([
        "snp_id", f"abs_logfc.mean.{ctx1}", f"abs_logfc.mean.{ctx2}",
        f"abs_logfc.mean.pval.{ctx1}", f"abs_logfc.mean.pval.{ctx2}",
        f"peak_overlap.{ctx1}", f"peak_overlap.{ctx2}",
    ])
    rows = [
        ("chr1:1000:A:G", "1.8", "0.5", "1e-6", "0.5", "1", "0"),
        ("chrX:200:C:T", "2.0", "0.1", "1e-5", "0.9", "1", "0"),
    ]
    with gzip.open(wide, "wt") as fh:   # gzipped -> exercises the streaming decompression path too
        fh.write(header + "\n")
        for r in rows:
            fh.write("\t".join(r) + "\n")

    args = parse_args()
    args.product, args.chrombpnet, args.cbp_wide, args.download = "chrombpnet", [str(wide)], True, False
    args.cbp_variant_col, args.cbp_score_prefix, args.cbp_pval_infix = "snp_id", "abs_logfc.mean.", ".pval."
    args.cbp_sep, args.cbp_batch_rows, args.cbp_tmpdir = "\t", 1, str(tmpdir)   # batch=1 forces >1 batch
    args.cbp_name_filter, args.mlog10p_thresh = None, 3.0
    args.variant_keep_chr = False

    out_path = str(tmpdir / f"{DATASET_CHROMBPNET}.wide.sample.tsv.gz")
    run_chrombpnet_streaming(args, out_path)

    # POINT index present + exact-pos queries
    assert Path(out_path + ".tbi").exists(), "POINT index .tbi not written"
    q1 = subprocess.run(["tabix", out_path, "1:1000-1000"], capture_output=True, text=True, check=True)
    assert q1.stdout.strip(), "POINT exact-pos query 1:1000 returned nothing"
    qx = subprocess.run(["tabix", out_path, "23:200-200"], capture_output=True, text=True, check=True)
    assert qx.stdout.strip(), "chrX variant not indexed under numeric seqname 23"

    body = subprocess.run(["bgzip", "-dc", out_path], capture_output=True, text=True, check=True).stdout
    lines = [l for l in body.strip().split("\n") if l]
    cols = lines[0].lstrip("#").split("\t")
    assert cols == VE_COLUMNS, f"18-col order mismatch: {cols}"
    data = [dict(zip(cols, l.split("\t"))) for l in lines[1:]]
    kept = {(d["variant"], d["cell_type"]) for d in data}
    assert kept == {("1:1000:A:G", ctx1), ("23:200:C:T", ctx1)}, kept   # ctx2 sub-threshold dropped
    for d in data:
        assert not d["chrom"].startswith("chr"), d          # numeric chrom
        assert not d["variant"].startswith("chr"), d        # numeric variant token
        assert d["is_significant"] == "true", d
        assert d["mlog10p"] not in ("", "NA") and float(d["mlog10p"]) >= 3.0, d
        assert d["score_type"] == SCORE_TYPE_CHROMBPNET_ABS, d
        assert d["predicted_direction"] == "NA" and d["rsid"] == "NA", d   # absolute + no rsid -> NA
        assert d["model"] == MODEL_CHROMBPNET, d
    chrom_of = {d["variant"]: d["chrom"] for d in data}
    assert chrom_of["23:200:C:T"] == "23", chrom_of        # chrX -> 23
    tissue_of = {d["cell_type"]: (d["tissue"], d["life_stage"]) for d in data}
    assert tissue_of[ctx1] == ("brain", "fetal"), tissue_of
    print("  wide->long reshape: sub-threshold (variant,context) dropped; is_significant/mlog10p set; "
          "numeric chrom+variant (chrX->23); NA nulls; 18-col order + POINT index: OK")
    print("  kept rows:")
    for d in data:
        print(f"    {d['chrom']}\t{d['variant']}\t{d['cell_type']}\t{d['score']}\t{d['mlog10p']}")


def run_sample() -> None:
    print("=== SAMPLE / DRY-RUN: synthetic input, no Synapse download, no GCS upload ===")
    tmpdir = Path(tempfile.mkdtemp())
    _sample_open_chromatin(tmpdir)
    _sample_chrombpnet(tmpdir)
    _sample_chrombpnet_wide(tmpdir)
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
