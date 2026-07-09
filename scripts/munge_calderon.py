#!/usr/bin/env python3
"""Munge the Calderon et al. 2019 Nat Genet stimulation-responsive immune ATAC-seq dataset into
the canonical `open_chromatin` TSV.

Source (public):
  - Paper: Calderon, Nguyen, Mezger, ..., Pritchard, Corces, Greenleaf, "Landscape of stimulation-
    responsive chromatin across diverse human immune cells", Nature Genetics 2019, 51:1494-1505.
    doi:10.1038/s41588-019-0505-9
  - Bulk ATAC-seq across ~25-32 sorted immune populations, each profiled RESTING (unstimulated)
    and STIMULATED, from multiple donors.
  - GEO: GSE118189. The single processed supplementary file is the peak x sample count matrix
    GSE118189_ATAC_counts.txt.gz (~111 MB). Row = peak (id in the row-name / first column,
    e.g. "chr1_100_200"); column = one sample named like "1001-Bulk_B-U" / "1001-Bulk_B-S"
    (donor - cell_type - stimulation) OR "Bulk_B-U" when no donor prefix is present.
  - Stimulation differential-accessibility supplement (paper supplementary table
    "Supplementary_data_3_ATAC_stimulation_DA_peaks"): per cell-type stimulated-vs-resting
    differential-accessibility peaks with a log2 fold-change / p-value. Optional secondary input.

>>> TO-VERIFY (format assumptions) <<<
  - GENOME BUILD: the Calderon peaks are hg19 (TO-VERIFY against the GEO record / methods). This
    dataset MUST be lifted over to hg38 (see LIFTOVER below) — this is the distinctive step.
  - The count-matrix peak-id and sample-name layouts (row-name peak id; "donor-cell-stim" column
    names) and the DA-supplement column names are documented below; verify against the real files
    and adjust the --*-col / --*-regex flags (no code change needed for column-name differences).

INPUT 1 — peak x sample count matrix  (--counts, REQUIRED for a real run)
  Assumed: a TSV/whitespace matrix whose first column (or --peak-col) holds the peak id
  ("chr1_100_200" / "chr1-100-200" / "chr1:100-200", hg19, BED 0-based half-open) and every other
  column is one sample. Sample column names encode cell type + stimulation state (and usually a
  donor prefix): "<donor>-<cell_type>-<U|S>" or "<cell_type>-<U|S>". Parsing (see _parse_sample):
    * split the column name on "-"; the LAST token is the stimulation state, the FIRST token is
      the donor (only when >=3 tokens), and everything in between is the cell_type (verbatim,
      e.g. "Bulk_B", "CD8pos_T", "Monocytes"). Calderon cell-type labels use underscores
      internally so "-" is a reliable separator.
    * condition: map {U, Unstim, Unstimulated, resting} -> "resting" and
      {S, Stim, Stimulated} -> "stimulated"; anything else -> "unknown".
  Values are raw ATAC counts. The atlas rows are aggregated ACROSS DONORS within each
  (peak, cell_type, condition) group (mean by default) and kept when the aggregate > --min-score.
  --score-type names the emitted score (default "mean_count"); pass "presence" to emit an empty
  score and keep any peak present (>0) in the group.

INPUT 2 — stimulation DA supplement  (--da-peaks, OPTIONAL)
  Assumed: a TSV of stimulated-vs-resting differential-accessibility peaks. Columns (configurable):
  a peak id/coords, a cell_type, a log2 fold-change (--da-log2fc-col), and an optional p-value.
  These describe the stimulation RESPONSE, so rows are emitted with condition="stimulated",
  score=log2FC and score_type="stim_log2fc" (independent of --score-type). They are APPENDED to the
  primary accessible-peak atlas (they do not replace it) and are lifted over together with it.

INPUT 3 — peak -> gene links  (--links, OPTIONAL)
  Usually absent for Calderon. If given: a TSV linking each peak to one or more target genes
  (peak id/coords + gene symbol + Ensembl id). Multiple genes per peak are aggregated into
  comma-separated lists with symbol and id paired POSITIONALLY (i-th symbol <-> i-th id), matching
  munge_li_brain.py / munge_catlas.py. Links are matched on the hg38 peak key AFTER liftOver.

Genome build: hg19 -> LIFTED OVER to hg38 (the distinctive requirement for this dataset).

OUTPUT — canonical `open_chromatin` long TSV, 18 columns IN THIS ORDER (tab-separated):
  chrom, start, end, peak_id, dataset, cell_type, tissue, life_stage, condition, assay,
  score, score_type, n_cells, cell_ontology_id, uberon_id, target_gene, target_gene_id, version
The header's first token is written as `#chrom` (comment-prefixed) following this repo's
convention (cf. munge_hpa.py `#dataset`, sumstats `#chr`): it makes the header a tabix meta line so
`tabix -p bed` skips it while the logical column name stays `chrom`.

Field rules for Calderon (see task genetics-results-suite-bzl.17):
  chrom            NUMERIC hg38, no "chr" prefix (1..22, X->23, Y->24, M/MT->25); converted from
                   the chr-prefixed liftOver output at the final write. REQUIRED by the tabix contract.
  start,end        peak interval, hg38 (AFTER liftOver), BED 0-based half-open.
  peak_id          f"{chrom}-{start}-{end}" using the numeric hg38 chrom (e.g. "23-100-200").
  dataset          "calderon_open_chromatin"  (drives resource mapping to calderon_immune).
  cell_type        verbatim immune population label parsed from the sample name (stim suffix
                   stripped into `condition`).
  tissue           "immune".
  life_stage       "adult".
  condition        {"resting","stimulated"} parsed from the sample name; "unknown" if unparseable.
  assay            "bulk_ATAC".
  score/score_type mean count + "mean_count" (or empty + "presence"); DA rows -> log2FC +
                   "stim_log2fc".
  n_cells          empty (bulk assay, not single cell).
  cell_ontology_id empty.
  uberon_id        empty by default; optional UBERON:0000178 (blood) via --blood-uberon.
  target_gene(_id) from --links when present, else empty.
  version          "2019".

LIFTOVER (hg19 -> hg38) — the distinctive step:
  Uses UCSC liftOver with the hg19ToHg38 chain. Default chain URL (documented, NOT auto-downloaded):
    https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz
  Binary and chain paths are parameterized (--liftover-bin, --chain). The step:
    1. write the UNIQUE hg19 peak intervals to a BED4 with name = hg19 peak key;
    2. run `liftOver in.bed chain out.bed unmapped.bed`;
    3. DROP unmapped intervals (the .unmapped output is discarded);
    4. DROP any peak whose name maps MORE THAN ONCE (multi-mapped / split) or maps to a DIFFERENT
       chromosome than its hg19 chrom (inconsistent);
    5. join the surviving hg38 coords back onto the long table by hg19 key, rebuild peak_id, and
       RE-SORT (-k1,1 -k2,2n) before bgzip (liftOver does not preserve coordinate order).
  --skip-liftover: for --sample / testing on ALREADY-hg38 synthetic input, bypass the liftOver call
  and pass coordinates through unchanged. THE REAL PIPELINE MUST run liftOver (do not pass
  --skip-liftover on the true hg19 data).

INDEXING CONTRACT (from the epic design, surfaced by the results-api open_chromatin review):
  open_chromatin files MUST be INTERVAL-indexed: `tabix -p bed` (i.e. -s1 -b2 -e3, distinct
  start/end columns). Do NOT point-index (-b2 -e2) or the API variant-overlap fast path would
  silently miss peaks whose interval contains pos. Pipeline: sort -k1,1 -k2,2n -> bgzip -> tabix.

STAGING (off by default): with --stage the .tsv.gz and .tsv.gz.tbi are uploaded to BOTH
  gs://finngen-commons/results_api_data/open_chromatin/calderon_immune/calderon_open_chromatin.tsv.gz
  gs://daly-genetics-results/open_chromatin/calderon_immune/calderon_open_chromatin.tsv.gz
  These paths must match the results-api profile config files. Without --stage nothing is uploaded.

Local validation without the full dataset or any upload (uses --skip-liftover on synthetic hg38):
  python3 scripts/munge_calderon.py --sample
"""

import argparse
import re
import shutil
import subprocess
import tempfile
from pathlib import Path

import polars as pl

DATASET = "calderon_open_chromatin"
RESOURCE = "calderon_immune"
VERSION = "2019"

TISSUE = "immune"
LIFE_STAGE = "adult"
ASSAY = "bulk_ATAC"

BLOOD_UBERON = "UBERON:0000178"  # optional, via --blood-uberon
DA_SCORE_TYPE = "stim_log2fc"

# documented default chain (NOT auto-downloaded; the real run must provide it via --chain)
DEFAULT_CHAIN_URL = "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz"

# canonical open_chromatin column order; first header token is comment-prefixed for tabix
OUTPUT_COLUMNS = [
    "chrom", "start", "end", "peak_id", "dataset", "cell_type", "tissue", "life_stage",
    "condition", "assay", "score", "score_type", "n_cells", "cell_ontology_id",
    "uberon_id", "target_gene", "target_gene_id", "version",
]

GCS_FINNGEN = f"gs://finngen-commons/results_api_data/open_chromatin/{RESOURCE}/{DATASET}.tsv.gz"
GCS_DALY = f"gs://daly-genetics-results/open_chromatin/{RESOURCE}/{DATASET}.tsv.gz"

# "chr1:1000-2000", "chr1-1000-2000", "chr1_1000_2000"; the "chr" prefix is optional and gets
# normalized on in _coords_from_id so ids like "1:1000-2000" are not dropped.
_ID_RE = re.compile(r"^(?:chr)?([^:_\-]+)[:_\-](\d+)[_\-](\d+)$")

# stimulation-state token (last "-"-token of a sample name) -> harmonized condition
_RESTING_TOKENS = {"u", "unstim", "unstimulated", "resting", "rest", "ctrl", "control"}
_STIM_TOKENS = {"s", "stim", "stimulated"}


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--counts", help="peak x sample count matrix TSV (row=peak, col=sample; see header)")
    p.add_argument("--peak-col", help="peak-id column in --counts (default: the first column)")
    p.add_argument("--score-type", default="mean_count",
                   help="score_type token for the count atlas (default: mean_count; use 'presence' for a binary atlas)")
    p.add_argument("--agg", choices=["mean", "sum", "max"], default="mean",
                   help="how to aggregate a sample's value across donors within a (peak,cell_type,condition) group (default: mean)")
    p.add_argument("--min-score", type=float, default=0.0,
                   help="keep (peak,cell_type,condition) groups with aggregate value > min-score (default: 0.0)")
    p.add_argument("--chunk-size", type=int, default=100_000,
                   help="peaks per batch when melting/aggregating the count matrix (default: 100000). "
                        "Each peak row is fully contained in one batch, so cross-donor aggregation stays "
                        "exact while peak RAM is bounded to ~one batch instead of the full melt.")
    p.add_argument("--blood-uberon", action="store_true",
                   help="fill uberon_id with UBERON:0000178 (blood) instead of leaving it empty")

    p.add_argument("--da-peaks", help="optional stimulation DA supplement TSV (Supplementary_data_3)")
    p.add_argument("--da-peak-col", help="DA: peak-id column to parse coords from")
    p.add_argument("--da-chrom-col", default="chrom", help="DA: chrom column (default: chrom)")
    p.add_argument("--da-start-col", default="start", help="DA: start column (default: start)")
    p.add_argument("--da-end-col", default="end", help="DA: end column (default: end)")
    p.add_argument("--da-celltype-col", default="cell_type", help="DA: cell_type column (default: cell_type)")
    p.add_argument("--da-log2fc-col", default="log2FoldChange", help="DA: log2 fold-change column (default: log2FoldChange)")

    p.add_argument("--links", help="optional peak->gene links TSV")
    p.add_argument("--links-id-col", help="links: peak id column to parse coords from")
    p.add_argument("--links-chrom-col", default="chrom", help="links: chrom column (default: chrom)")
    p.add_argument("--links-start-col", default="start", help="links: start column (default: start)")
    p.add_argument("--links-end-col", default="end", help="links: end column (default: end)")
    p.add_argument("--links-gene-col", default="gene_name", help="links: gene symbol column (default: gene_name)")
    p.add_argument("--links-geneid-col", default="gene_id", help="links: Ensembl gene id column (default: gene_id)")

    p.add_argument("--liftover-bin", default="liftOver", help="UCSC liftOver binary (default: liftOver on PATH)")
    p.add_argument("--chain", help=f"hg19ToHg38 chain(.gz) file. Default source: {DEFAULT_CHAIN_URL}")
    p.add_argument("--skip-liftover", action="store_true",
                   help="TESTING ONLY: input is already hg38; bypass liftOver. The real hg19 run must NOT set this.")

    p.add_argument("--output", help="output .tsv.gz path (default: ./calderon_open_chromatin.tsv.gz)")
    p.add_argument("--stage", action="store_true", help="upload .tsv.gz + .tbi to BOTH GCS buckets (OFF by default)")
    p.add_argument("--gcs-finngen", default=GCS_FINNGEN, help="finngen GCS destination")
    p.add_argument("--gcs-daly", default=GCS_DALY, help="daly GCS destination")

    p.add_argument("--sample", "--dry-run", dest="sample", action="store_true",
                   help="run the transform on tiny synthetic hg38 input (--skip-liftover), print rows, validate tabix; no upload")
    return p.parse_args()


# ---------------------------------------------------------------------------
# coordinate helpers (shared shape with munge_li_brain.py / munge_catlas.py)
# ---------------------------------------------------------------------------
def _coords_from_id(df: pl.DataFrame, id_col: str) -> pl.DataFrame:
    """Parse chrom/start/end from a single id column like chr1:1000-2000 (chr prefix normalized on)."""
    parsed = df.select(
        pl.col(id_col).cast(pl.Utf8).str.extract_groups(_ID_RE.pattern).alias("_g")
    ).unnest("_g")
    chrom = parsed["1"]
    chrom = pl.when(chrom.str.starts_with("chr")).then(chrom).otherwise("chr" + chrom)
    return df.with_columns(
        chrom.alias("chrom"),
        parsed["2"].cast(pl.Int64).alias("start"),
        parsed["3"].cast(pl.Int64).alias("end"),
    )


def _normalize_coords(df: pl.DataFrame, chrom_col: str, start_col: str, end_col: str) -> pl.DataFrame:
    """Rename explicit coordinate columns to chrom/start/end and enforce chr-prefix + int coords."""
    df = df.rename({chrom_col: "chrom", start_col: "start", end_col: "end"})
    return df.with_columns(
        pl.when(pl.col("chrom").cast(pl.Utf8).str.starts_with("chr"))
        .then(pl.col("chrom").cast(pl.Utf8))
        .otherwise("chr" + pl.col("chrom").cast(pl.Utf8))
        .alias("chrom"),
        pl.col("start").cast(pl.Int64),
        pl.col("end").cast(pl.Int64),
    )


def _peak_key() -> pl.Expr:
    return pl.format("{}-{}-{}", pl.col("chrom"), pl.col("start"), pl.col("end"))


def _numeric_chrom(expr: pl.Expr) -> pl.Expr:
    """Strip any 'chr' prefix and map to a numeric-string chrom (X=23, Y=24, M/MT=25).

    Applied only at the FINAL (hg38) write: the coord helpers keep chr-prefixed names because
    the liftOver chain (hg19ToHg38) is keyed on chr-prefixed seqnames.
    """
    e = expr.cast(pl.Utf8).str.replace(r"^chr", "")
    up = e.str.to_uppercase()
    return (
        pl.when(up == "X").then(pl.lit("23"))
        .when(up == "Y").then(pl.lit("24"))
        .when(up.is_in(["M", "MT"])).then(pl.lit("25"))
        .otherwise(e)
    )


# canonical primary-assembly chromosomes as numeric strings (1..22, X=23, Y=24, M/MT=25). the
# platform is primary-assembly only and loads chrom as INT64, so non-canonical hg38 contigs
# (alt/random/scaffold/unplaced/Un_*) must be DROPPED or they break the BigQuery chr load. liftOver
# can emit alt/random contigs from an hg19 primary locus, so this is applied AFTER liftOver.
CANONICAL_CHROMS = frozenset(str(c) for c in range(1, 26))


# ---------------------------------------------------------------------------
# sample-name parsing: "<donor>-<cell_type>-<U|S>" or "<cell_type>-<U|S>"
# ---------------------------------------------------------------------------
def _parse_sample(sample: str) -> tuple[str, str]:
    """Return (cell_type, condition) parsed from a Calderon sample column name.

    Split on "-": last token = stimulation state, first token = donor (only when >=3 tokens),
    the middle = cell_type verbatim. Calderon cell-type labels use "_" internally so "-" is a safe
    separator. condition maps to resting/stimulated/unknown.
    """
    parts = sample.split("-")
    if len(parts) >= 3:
        cell_type = "-".join(parts[1:-1])
        stim = parts[-1]
    elif len(parts) == 2:
        cell_type, stim = parts[0], parts[1]
    else:
        return sample, "unknown"

    tok = stim.strip().lower()
    if tok in _STIM_TOKENS:
        condition = "stimulated"
    elif tok in _RESTING_TOKENS:
        condition = "resting"
    else:
        condition = "unknown"
    return cell_type, condition


def _sample_map(sample_cols: list[str]) -> pl.DataFrame:
    """One-shot sample-column -> (cell_type, condition) lookup table for a join."""
    mapping = {s: _parse_sample(s) for s in sample_cols}
    return pl.DataFrame(
        [(s, ct, cond) for s, (ct, cond) in mapping.items()],
        schema={"_sample": pl.Utf8, "cell_type": pl.Utf8, "condition": pl.Utf8},
        orient="row",
    )


def _aggregate_batch(df: pl.DataFrame, peak_col: str, sample_cols: list[str],
                     map_df: pl.DataFrame, args: argparse.Namespace) -> pl.DataFrame:
    """Melt one peak-row batch to long, parse sample -> (cell_type, condition), aggregate across
    donors within (peak, cell_type, condition), keep aggregate > min-score. Because a peak is a
    single matrix ROW, every sample for that peak is inside this batch, so the aggregation is exact
    per batch — no cross-batch merge is ever required.
    """
    df = _coords_from_id(df, peak_col)
    long = df.unpivot(
        index=["chrom", "start", "end"], on=sample_cols,
        variable_name="_sample", value_name="_value",
    ).with_columns(pl.col("_value").cast(pl.Float64, strict=False))
    long = long.join(map_df, on="_sample", how="left").drop("_sample")

    agg_expr = {"mean": pl.mean, "sum": pl.sum, "max": pl.max}[args.agg]("_value")
    grouped = long.group_by(["chrom", "start", "end", "cell_type", "condition"]).agg(
        agg_expr.alias("_agg")
    ).filter(pl.col("_agg") > args.min_score)

    if args.score_type == "presence":
        grouped = grouped.with_columns(pl.lit(None, dtype=pl.Utf8).alias("score"))
    else:
        grouped = grouped.with_columns(pl.col("_agg").cast(pl.Utf8).alias("score"))
    return grouped.drop("_agg").with_columns(pl.lit(args.score_type).alias("score_type"))


def load_counts_chunked(args: argparse.Namespace, work_dir: Path) -> None:
    """Stream the peak x sample matrix in --chunk-size peak batches, aggregating each batch to the
    long hg19 shape (chrom, start, end, cell_type, condition, score, score_type) and writing it to
    an on-disk parquet part in `work_dir` INCREMENTALLY. Peak RAM stays bounded to ~one batch's melt
    rather than the full ~145M-row melt that OOM-killed the single-shot path.
    """
    header = pl.read_csv(args.counts, separator="\t", n_rows=0).columns
    peak_col = args.peak_col if (args.peak_col and args.peak_col in header) else header[0]
    sample_cols = [c for c in header if c != peak_col]
    if not sample_cols:
        raise ValueError("no sample value columns found in --counts (check --peak-col)")
    map_df = _sample_map(sample_cols)

    reader = pl.read_csv_batched(
        args.counts, separator="\t", batch_size=args.chunk_size, infer_schema_length=10_000,
    )
    part = n_peaks = n_rows = 0
    while True:
        batches = reader.next_batches(1)
        if not batches:
            break
        for batch in batches:
            n_peaks += batch.height
            grouped = _aggregate_batch(batch, peak_col, sample_cols, map_df, args)
            grouped.write_parquet(work_dir / f"counts_{part:05d}.parquet")
            n_rows += grouped.height
            part += 1
    print(f"  chunked melt/aggregate: {n_peaks} peaks in {part} batches "
          f"(<= {args.chunk_size}/batch) -> {n_rows} aggregated hg19 rows on disk")


def load_da_peaks(args: argparse.Namespace) -> pl.DataFrame:
    """Read the stimulation DA supplement -> long hg19 rows (condition=stimulated, score=log2FC)."""
    df = pl.read_csv(args.da_peaks, separator="\t", infer_schema_length=10_000)
    if args.da_peak_col:
        df = _coords_from_id(df, args.da_peak_col)
    else:
        df = _normalize_coords(df, args.da_chrom_col, args.da_start_col, args.da_end_col)

    cell_type = (
        pl.col(args.da_celltype_col).cast(pl.Utf8) if args.da_celltype_col in df.columns
        else pl.lit("unknown")
    )
    score = (
        pl.col(args.da_log2fc_col).cast(pl.Utf8) if args.da_log2fc_col in df.columns
        else pl.lit(None, dtype=pl.Utf8)
    )
    return df.select(
        "chrom", "start", "end",
        cell_type.alias("cell_type"),
        pl.lit("stimulated").alias("condition"),
        score.alias("score"),
        pl.lit(DA_SCORE_TYPE).alias("score_type"),
    )


def load_gene_links(args: argparse.Namespace) -> pl.DataFrame:
    """Read peak->gene links (hg38) and aggregate to one comma-joined gene list per hg38 peak key,
    pairing symbol and id positionally (identical logic to munge_li_brain.py / munge_catlas.py)."""
    df = pl.read_csv(args.links, separator="\t", infer_schema_length=10_000)
    if args.links_id_col:
        df = _coords_from_id(df, args.links_id_col)
    else:
        df = _normalize_coords(df, args.links_chrom_col, args.links_start_col, args.links_end_col)

    # links are already hg38; convert to numeric chrom so the key matches the final hg38 peak_key
    df = df.with_columns(_numeric_chrom(pl.col("chrom")).alias("chrom"))
    df = df.with_columns(_peak_key().alias("_key"))
    gene = pl.col(args.links_gene_col).cast(pl.Utf8) if args.links_gene_col in df.columns else pl.lit(None, dtype=pl.Utf8)
    geneid = pl.col(args.links_geneid_col).cast(pl.Utf8) if args.links_geneid_col in df.columns else pl.lit(None, dtype=pl.Utf8)
    # normalize nulls to empty strings so a missing symbol or id still keeps its pair and emits an
    # empty token on the missing side, keeping the two comma-joined lists aligned
    df = df.with_columns(gene.fill_null("").alias("_gene"), geneid.fill_null("").alias("_geneid"))
    return (
        df.group_by("_key")
        .agg(pl.struct("_gene", "_geneid").unique().sort().alias("_pairs"))
        .with_columns(
            pl.col("_pairs").list.eval(pl.element().struct.field("_gene")).list.join(",").alias("target_gene"),
            pl.col("_pairs").list.eval(pl.element().struct.field("_geneid")).list.join(",").alias("target_gene_id"),
        )
        .drop("_pairs")
    )


# ---------------------------------------------------------------------------
# hg19 -> hg38 liftOver
# ---------------------------------------------------------------------------
def liftover_intervals(intervals: pl.DataFrame, args: argparse.Namespace) -> pl.DataFrame:
    """Lift the UNIQUE hg19 intervals to hg38 via UCSC liftOver.

    intervals: distinct chrom/start/end (hg19). Returns a map table:
      _hg19key (chrom-start-end, hg19) -> new_chrom, new_start, new_end (hg38).
    Drops unmapped, multi-mapped, and cross-chromosome (inconsistent) intervals.
    When --skip-liftover is set the coordinates pass through unchanged (testing on hg38 input).
    """
    intervals = intervals.with_columns(_peak_key().alias("_hg19key")).unique(subset=["_hg19key"])

    if args.skip_liftover:
        print("  --skip-liftover: passing coordinates through unchanged (TESTING; real run must liftOver)")
        return intervals.select(
            pl.col("_hg19key"),
            pl.col("chrom").alias("new_chrom"),
            pl.col("start").alias("new_start"),
            pl.col("end").alias("new_end"),
        )

    if not args.chain:
        raise SystemExit(f"liftOver requires --chain (hg19ToHg38). Default source: {DEFAULT_CHAIN_URL}")
    if not Path(args.chain).exists():
        raise SystemExit(f"chain file not found: {args.chain} (download from {DEFAULT_CHAIN_URL})")
    if shutil.which(args.liftover_bin) is None and not Path(args.liftover_bin).exists():
        raise SystemExit(f"liftOver binary not found: {args.liftover_bin} (UCSC liftOver required for the real run)")

    tmp = Path(tempfile.mkdtemp())
    bed_in, bed_out, unmapped = tmp / "hg19.bed", tmp / "hg38.bed", tmp / "unmapped.bed"
    # BED4: name column carries the hg19 key so we can join the hg38 result back
    intervals.select("chrom", "start", "end", "_hg19key").write_csv(
        bed_in, separator="\t", include_header=False
    )
    subprocess.run(
        [args.liftover_bin, str(bed_in), args.chain, str(bed_out), str(unmapped)], check=True
    )

    mapped = pl.read_csv(
        bed_out, separator="\t", has_header=False,
        new_columns=["new_chrom", "new_start", "new_end", "_hg19key"],
    )
    n_lifted = mapped.height
    # drop intervals whose hg19 key maps MORE THAN ONCE (multi-mapped / split -> inconsistent)
    counts = mapped.group_by("_hg19key").len()
    unique_keys = counts.filter(pl.col("len") == 1).select("_hg19key")
    mapped = mapped.join(unique_keys, on="_hg19key", how="inner")
    n_after_multi = mapped.height
    # drop intervals that landed on a DIFFERENT chromosome than the hg19 chrom (encoded in the key)
    mapped = mapped.with_columns(
        pl.col("_hg19key").str.extract(r"^(.+)-\d+-\d+$", 1).alias("_hg19chrom")
    ).filter(pl.col("new_chrom") == pl.col("_hg19chrom")).drop("_hg19chrom")

    print(f"  liftOver: {intervals.height} hg19 intervals -> {n_lifted} lifted; "
          f"dropped {n_lifted - n_after_multi} multi-mapped, "
          f"{n_after_multi - mapped.height} cross-chromosome; kept {mapped.height}")
    return mapped.select("_hg19key", "new_chrom", "new_start", "new_end")


# ---------------------------------------------------------------------------
# assemble canonical output
# ---------------------------------------------------------------------------
def build_output(long_hg38: pl.LazyFrame, args: argparse.Namespace) -> pl.LazyFrame:
    """Assemble the 18 canonical columns from a lifted (hg38) long table.

    long_hg38 carries: chrom, start, end, cell_type, condition, score, score_type. Operates lazily
    so the full atlas is never materialized in RAM (the count-matrix melt/aggregate already happened
    on disk); the caller streams the result to CSV via sink.
    """
    # convert the final hg38 chrom to numeric (X->23, ...) AFTER liftOver, then build peak_id/keys
    df = long_hg38.with_columns(_numeric_chrom(pl.col("chrom")).alias("chrom"))
    # liftOver can map an hg19 primary locus onto an hg38 alt/random contig; drop those non-canonical
    # seqnames so the platform's INT64 chr load never sees a scaffold/alt/Un name
    df = df.filter(pl.col("chrom").is_in(CANONICAL_CHROMS))
    df = df.with_columns(_peak_key().alias("peak_id"))

    if args.links:
        links = load_gene_links(args).lazy()
        df = df.with_columns(_peak_key().alias("_key")).join(links, on="_key", how="left").drop("_key")
    else:
        df = df.with_columns(
            pl.lit(None, dtype=pl.Utf8).alias("target_gene"),
            pl.lit(None, dtype=pl.Utf8).alias("target_gene_id"),
        )

    df = df.with_columns(
        pl.lit(DATASET).alias("dataset"),
        pl.lit(TISSUE).alias("tissue"),
        pl.lit(LIFE_STAGE).alias("life_stage"),
        pl.lit(ASSAY).alias("assay"),
        pl.lit(None, dtype=pl.Utf8).alias("n_cells"),
        pl.lit(None, dtype=pl.Utf8).alias("cell_ontology_id"),
        pl.lit(BLOOD_UBERON if args.blood_uberon else None, dtype=pl.Utf8).alias("uberon_id"),
        pl.lit(VERSION).alias("version"),
    )
    return df.select(OUTPUT_COLUMNS)


def _sort_index_bgzip(body: Path, output_path: str) -> None:
    """external LC_ALL=C sort -k1,1 -k2,2n -> bgzip -> tabix -p bed (interval index) of a body TSV.

    The on-disk sort keeps memory bounded and absorbs liftOver's coordinate reordering before bgzip."""
    header = "#" + "\t".join(OUTPUT_COLUMNS)
    pipeline = (
        f'( printf "%s\\n" {shell_quote(header)}; '
        f'LC_ALL=C sort -k1,1 -k2,2n {shell_quote(str(body))} ) | bgzip -c > {shell_quote(output_path)}'
    )
    subprocess.run(pipeline, shell=True, check=True, executable="/bin/bash")
    subprocess.run(["tabix", "-f", "-p", "bed", output_path], check=True)
    print(f"  indexed {output_path}.tbi (tabix -p bed / -s1 -b2 -e3)")


def write_open_chromatin(df: pl.DataFrame, output_path: str) -> None:
    """In-memory writer (used by --sample): serialize the small frame then sort/index/bgzip."""
    tmpdir = tempfile.mkdtemp()
    body = Path(tmpdir) / "body.tsv"
    # every empty/missing cell serialized as the literal "NA"
    df.write_csv(body, separator="\t", include_header=False, null_value="NA")
    _sort_index_bgzip(body, output_path)
    print(f"  wrote {df.height} rows -> {output_path}")


def write_open_chromatin_lazy(lf: pl.LazyFrame, output_path: str) -> None:
    """Streaming writer (real run): sink the lazy atlas to an on-disk body TSV without materializing
    it in RAM, then sort/index/bgzip. Keeps the final write memory-bounded like the melt/aggregate."""
    tmpdir = tempfile.mkdtemp()
    body = Path(tmpdir) / "body.tsv"
    lf.sink_csv(body, separator="\t", include_header=False, null_value="NA",
                maintain_order=False, engine="streaming")
    _sort_index_bgzip(body, output_path)
    print(f"  streamed atlas -> {output_path}")


def shell_quote(s: str) -> str:
    return "'" + s.replace("'", "'\\''") + "'"


def upload_to_gcs(local_path: str, gcs_path: str) -> None:
    subprocess.run(["gcloud", "storage", "cp", local_path, gcs_path], check=True)
    print(f"  uploaded {gcs_path}")
    tbi = local_path + ".tbi"
    if Path(tbi).exists():
        subprocess.run(["gcloud", "storage", "cp", tbi, gcs_path + ".tbi"], check=True)
        print(f"  uploaded {gcs_path}.tbi")


def transform(args: argparse.Namespace, work_dir: Path) -> pl.LazyFrame:
    """counts (+ optional DA) -> chunked/disk-backed long hg19 -> liftOver -> hg38 -> 18 columns.

    Returns a LazyFrame scanning the on-disk aggregated parts joined to the hg38 liftOver mapping.
    Nothing except the small unique-peak set and the liftOver mapping is held in RAM; `work_dir`
    holds the parquet parts and must outlive the returned frame (until it is sunk/collected).
    """
    load_counts_chunked(args, work_dir)
    if args.da_peaks:
        # DA rows share the long hg19 schema and are lifted over together with the count atlas
        load_da_peaks(args).write_parquet(work_dir / "da_peaks.parquet")

    lf = pl.scan_parquet(str(work_dir / "*.parquet"))
    peaks = lf.select("chrom", "start", "end").unique().collect(engine="streaming")
    mapping = liftover_intervals(peaks, args)  # small: one row per unique hg19 interval

    lifted = (
        lf.with_columns(_peak_key().alias("_hg19key"))
        .join(mapping.lazy(), on="_hg19key", how="inner")
        .drop("_hg19key", "chrom", "start", "end")
        .rename({"new_chrom": "chrom", "new_start": "start", "new_end": "end"})
    )
    out = build_output(lifted, args)
    names = out.collect_schema().names()
    assert names == OUTPUT_COLUMNS, f"column order mismatch: {names}"
    return out


# ---------------------------------------------------------------------------
# sample / dry-run
# ---------------------------------------------------------------------------
def _synthetic_inputs(tmpdir: Path) -> str:
    """Write a tiny synthetic hg38 count matrix mimicking Calderon sample-name conventions.

    Exercises: donor prefixes ("1001-"/"1002-") aggregated across donors, the no-donor 2-token
    form ("Monocytes-U"), and both stimulation states (-U resting / -S stimulated).
    """
    counts = tmpdir / "atac_counts.tsv"
    counts.write_text(
        "peak\t1001-Bulk_B-U\t1001-Bulk_B-S\t1002-Bulk_B-U\t1001-CD8pos_T-S\tMonocytes-U\n"
        "chr1_100100_100600\t4\t20\t6\t0\t0\n"
        "chr1_200200_200700\t10\t0\t8\t15\t0\n"
        "chr2_50050_50550\t0\t0\t0\t0\t9\n"
        "chrX_900900_901400\t3\t3\t5\t0\t0\n"
        # non-canonical unplaced contig (liftOver can emit alt/random/Un contigs) — MUST be dropped
        # by the canonical filter that runs AFTER the numeric-chrom conversion in build_output
        "chrUn_300300_300800\t8\t8\t0\t0\t0\n"
    )
    return str(counts)


def run_sample() -> None:
    print("=== SAMPLE / DRY-RUN: synthetic hg38 input, --skip-liftover, no GCS upload ===")
    tmpdir = Path(tempfile.mkdtemp())
    counts = _synthetic_inputs(tmpdir)

    args = parse_args()
    args.counts = counts
    args.peak_col = "peak"
    args.score_type = "mean_count"
    args.agg = "mean"
    args.min_score = 0.0
    args.da_peaks = None
    args.links = None
    args.blood_uberon = False
    args.skip_liftover = True  # synthetic input is already hg38
    if not getattr(args, "chunk_size", None):
        args.chunk_size = 100_000

    # condition parsing from sample names
    assert _parse_sample("1001-Bulk_B-U") == ("Bulk_B", "resting"), _parse_sample("1001-Bulk_B-U")
    assert _parse_sample("1001-Bulk_B-S") == ("Bulk_B", "stimulated"), _parse_sample("1001-Bulk_B-S")
    assert _parse_sample("Monocytes-U") == ("Monocytes", "resting"), _parse_sample("Monocytes-U")
    assert _parse_sample("1001-CD8pos_T-S") == ("CD8pos_T", "stimulated"), _parse_sample("1001-CD8pos_T-S")
    print("  sample-name -> (cell_type, condition) parsing: OK")

    # exercise the real chunked/disk-backed path, then collect the tiny result for assertions
    out = transform(args, Path(tempfile.mkdtemp())).collect(engine="streaming")
    assert out.columns == OUTPUT_COLUMNS, f"column order mismatch: {out.columns}"
    print(f"  output has {len(out.columns)} columns in canonical order: OK")

    # conditions parsed correctly into the resting/stimulated axis
    conds = set(out["condition"].unique().to_list())
    assert conds <= {"resting", "stimulated", "unknown"}, conds
    assert conds == {"resting", "stimulated"}, conds
    print(f"  condition axis: OK  {sorted(conds)}")

    # chrom / peak_id are numeric with chrX -> 23
    chroms = set(out["chrom"].to_list())
    assert chroms <= {"1", "2", "23"}, chroms
    assert "23" in chroms, f"chrX should map to 23; got {chroms}"
    xrow = out.filter(pl.col("chrom") == "23").row(0, named=True)
    assert xrow["peak_id"] == "23-900900-901400", xrow["peak_id"]
    print(f"  numeric chrom (chrX -> 23): OK  chroms={sorted(chroms)}  peak_id={xrow['peak_id']}")

    # non-canonical contigs (alt/random/scaffold/Un that liftOver can emit) are DROPPED (Fix A)
    # AFTER the numeric-chrom conversion; the "chrUn" synthetic peak must not survive while chrX does
    assert not any("Un" in str(p) for p in out["peak_id"].to_list()), "non-canonical contig leaked"
    assert chroms == {"1", "2", "23"}, f"expected canonical chroms only; got {chroms}"
    print("  non-canonical contig (chrUn) dropped: OK")

    # Bulk_B resting aggregates donors 1001(4) + 1002(6) -> mean 5.0 at 1-100100-100600
    bulk_b_rest = out.filter(
        (pl.col("peak_id") == "1-100100-100600")
        & (pl.col("cell_type") == "Bulk_B")
        & (pl.col("condition") == "resting")
    )
    assert bulk_b_rest.height == 1, bulk_b_rest
    assert abs(float(bulk_b_rest["score"][0]) - 5.0) < 1e-9, bulk_b_rest["score"][0]
    print(f"  cross-donor aggregation (mean of 4,6 -> 5.0): OK")

    print("\nfirst output rows:")
    with pl.Config(tbl_cols=-1, tbl_width_chars=260, fmt_str_lengths=40):
        print(out.sort("chrom", "start").head(10))

    out_path = str(tmpdir / f"{DATASET}.sample.tsv.gz")
    write_open_chromatin(out, out_path)

    # empty cells serialize as the literal "NA" (n_cells/cell_ontology_id/uberon_id/target_gene*);
    # every row still has 18 columns
    import gzip
    with gzip.open(out_path, "rt") as fh:
        body_lines = [ln for ln in fh if not ln.startswith("#")]
    assert all(len(ln.rstrip("\n").split("\t")) == 18 for ln in body_lines), "not 18 columns"
    assert any("\tNA\t" in ln or ln.rstrip("\n").endswith("\tNA") for ln in body_lines), "no NA emitted"
    print(f"  empty cells rendered as NA, 18 columns per row: OK ({len(body_lines)} rows)")

    # validate the interval index answers overlap queries (numeric seqnames, incl. chrX->23)
    q = subprocess.run(["tabix", out_path, "1:200300-200400"], capture_output=True, text=True, check=True)
    print("\ntabix overlap query 1:200300-200400 ->")
    print(q.stdout.rstrip() or "  (no rows)")
    assert q.stdout.strip(), "interval query returned nothing — indexing is wrong"
    qx = subprocess.run(["tabix", out_path, "23:901000-901100"], capture_output=True, text=True, check=True)
    print("tabix overlap query 23:901000-901100 (chrX) ->")
    print(qx.stdout.rstrip() or "  (no rows)")
    assert qx.stdout.strip(), "chrX (23) interval query returned nothing — indexing is wrong"

    if shutil.which("liftOver") is not None:
        print("\n  note: liftOver binary IS present; provide --chain to smoke-test the real liftOver path.")
    else:
        print("\n  note: liftOver binary NOT present here — the real hg19->hg38 liftOver runs in the "
              "pipeline image with UCSC liftOver + the hg19ToHg38 chain (see --help).")
    print("\n=== SAMPLE OK (no upload performed) ===")


def main() -> None:
    args = parse_args()
    if args.sample:
        run_sample()
        return

    if not args.counts:
        raise SystemExit("--counts is required (or use --sample). See --help.")
    if not args.skip_liftover and not args.chain:
        raise SystemExit(
            "Calderon is hg19 and MUST be lifted over. Pass --chain <hg19ToHg38.over.chain.gz> "
            f"(source: {DEFAULT_CHAIN_URL}), or --skip-liftover for TESTING on hg38 input only."
        )

    output_path = args.output or f"./{DATASET}.tsv.gz"
    print(f"Reading counts {args.counts} (chunk-size={args.chunk_size} peaks/batch) ...")
    work_dir = Path(tempfile.mkdtemp())
    out = transform(args, work_dir)
    write_open_chromatin_lazy(out, output_path)

    if args.stage:
        print("Staging to GCS (both buckets) ...")
        upload_to_gcs(output_path, args.gcs_finngen)
        upload_to_gcs(output_path, args.gcs_daly)
    else:
        print("  --stage not set: skipping GCS upload")
    print("Done.")


if __name__ == "__main__":
    main()
