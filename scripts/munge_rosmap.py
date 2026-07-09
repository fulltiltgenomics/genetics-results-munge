#!/usr/bin/env python3
"""Munge the Xiong et al. 2023 Cell ROSMAP AD-brain snATAC epigenome into the canonical
`open_chromatin` TSV.

Source (public):
  - Paper: Xiong et al., "Epigenomic dissection of Alzheimer's disease pinpoints causal
    variants and reveals epigenome erosion", Cell 2023, 186(20):4422-4437.
    doi:10.1016/j.cell.2023.08.040
  - snATAC-seq of the aged human PREFRONTAL CORTEX: ~850k nuclei, 92 ROSMAP donors spanning
    Alzheimer's disease (AD) and non-AD controls. Genome build hg38.
  - Processed per-cell-type / per-subtype accessible-peak BEDs and peak->gene links are the
    PUBLIC derivative distributed via the Broad `ad_epigenome` resource. Only these public
    processed peaks are consumed here — NO controlled-access / AD Knowledge Portal DUC data
    (raw fragments / genotypes) is touched.

>>> TO-VERIFY (format assumptions) <<<
The exact `ad_epigenome` download URLs, file naming, and per-file column layout could not be
scraped programmatically. The inputs this script consumes and the layouts it assumes are
documented below. The user MUST verify these against the real files and adjust the --*-idx /
--*-from-filename flags (no code change needed for layout differences).

  - CELL TYPE: derived per-file from the filename (default), or from a peak-BED column
    (--celltype-idx), or a fixed --cell-type. ROSMAP cell-type labels are kept VERBATIM
    (e.g. microglia, "excitatory neurons", oligodendrocytes, astrocytes, OPCs).
  - CONDITION (AD-vs-control axis): >>> TO-VERIFY whether the public processed release splits
    peaks by AD/control or ships a single UNION peak set per cell type. If per-condition, set it
    per-file (--condition-from-filename) or per-row (--condition-idx); if union peaks, leave the
    default and every row gets condition="unknown". Synonyms are normalized to {AD, control,
    unknown} (see _normalize_condition; ROSMAP "NCI" / "non-AD" -> control).

INPUT 1 — per-cell-type accessible-peak BEDs  (--peaks-dir OR --peaks, REQUIRED for a real run)
  Assumed: one BED per cell type (optionally per condition), headerless BED, hg38, 0-based
  half-open, columns: chrom, start, end [, name [, score ...]]. Coordinate columns default to
  positions 0/1/2; an optional signal/score column is read with --score-idx (e.g. narrowPeak
  signalValue at index 6, or a plain score at index 4). Pass --bed-has-header if the files carry
  a header row (coordinates are still read positionally; the header line is skipped).
    * --peaks-dir DIR   : glob DIR for --glob (default "*.bed*") and load every match.
    * --peaks A B C     : explicit list of BED files.
  cell_type / condition are resolved per file (see below).

INPUT 2 — peak -> gene links  (--links, OPTIONAL)
  Xiong 2023 provides peak-to-gene links. Assumed: a TSV linking each peak to one or more target
  genes. Columns (configurable): a peak id/coords, a gene symbol, and an Ensembl gene id.
  Multiple genes per peak are aggregated into comma-separated lists with symbol and id paired
  POSITIONALLY (i-th symbol <-> i-th id), matching munge_li_brain.py / munge_catlas.py — the
  pairs are aggregated JOINTLY (never independently sorted).

Genome build: hg38 (per the paper and the epic design) — NO liftOver.

OUTPUT — canonical `open_chromatin` long TSV, 18 columns IN THIS ORDER (tab-separated):
  chrom, start, end, peak_id, dataset, cell_type, tissue, life_stage, condition, assay,
  score, score_type, n_cells, cell_ontology_id, uberon_id, target_gene, target_gene_id, version
The header's first token is written as `#chrom` (comment-prefixed) following this repo's
convention (cf. munge_hpa.py `#dataset`, sumstats `#chr`): it makes the header a tabix meta line
so `tabix -p bed` skips it while the logical column name stays `chrom`.

Field rules for ROSMAP / Xiong 2023 (see task genetics-results-suite-bzl.18):
  chrom            NUMERIC hg38, no "chr" prefix (1..22, X->23, Y->24, M/MT->25); REQUIRED by
                   the tabix contract.
  start,end        peak interval, hg38, BED 0-based half-open (kept verbatim; no liftOver).
  peak_id          f"{chrom}-{start}-{end}" with the numeric chrom (e.g. "23-100-200").
  dataset          "rosmap_open_chromatin"  (drives resource mapping to rosmap_brain).
  cell_type        verbatim brain cell-type / subtype label (from filename or a column).
  tissue           "brain" (prefrontal cortex; override with --tissue).
  life_stage       "adult" (aged cohort).
  condition        {"AD","control"} when the peaks are per-condition, else "unknown" (union).
  assay            "snATAC".
  score/score_type peak signal + --score-type (default "score" when --score-idx is given), else
                   empty + "presence".
  n_cells          empty.
  cell_ontology_id empty.
  uberon_id        empty by default; optional UBERON:0000451 (prefrontal cortex) via --pfc-uberon.
  target_gene(_id) from --links when present, else empty.
  version          "2023".

INDEXING CONTRACT (from the epic design, surfaced by the results-api open_chromatin review):
  open_chromatin files MUST be INTERVAL-indexed: `tabix -p bed` (i.e. -s1 -b2 -e3, distinct
  start/end columns). Do NOT point-index (-b2 -e2) or the API variant-overlap fast path would
  silently miss peaks whose interval contains pos. Pipeline: sort -k1,1 -k2,2n -> bgzip -> tabix.

STAGING (off by default): with --stage the .tsv.gz and .tsv.gz.tbi are uploaded to BOTH
  gs://finngen-commons/results_api_data/open_chromatin/rosmap_brain/rosmap_open_chromatin.tsv.gz
  gs://daly-genetics-results/open_chromatin/rosmap_brain/rosmap_open_chromatin.tsv.gz
  These paths must match the results-api profile config files. Without --stage nothing is uploaded.

Local validation without the full dataset or any upload:
  python3 scripts/munge_rosmap.py --sample
"""

import argparse
import glob as globmod
import re
import subprocess
import tempfile
from pathlib import Path

import polars as pl

DATASET = "rosmap_open_chromatin"
RESOURCE = "rosmap_brain"
VERSION = "2023"

TISSUE = "brain"
LIFE_STAGE = "adult"
ASSAY = "snATAC"

# canonical open_chromatin column order; first header token is comment-prefixed for tabix
OUTPUT_COLUMNS = [
    "chrom", "start", "end", "peak_id", "dataset", "cell_type", "tissue", "life_stage",
    "condition", "assay", "score", "score_type", "n_cells", "cell_ontology_id",
    "uberon_id", "target_gene", "target_gene_id", "version",
]

GCS_FINNGEN = f"gs://finngen-commons/results_api_data/open_chromatin/{RESOURCE}/{DATASET}.tsv.gz"
GCS_DALY = f"gs://daly-genetics-results/open_chromatin/{RESOURCE}/{DATASET}.tsv.gz"

# "chr1:1000-2000", "chr1-1000-2000", "chr1_1000_2000"; "chr" prefix optional, normalized on.
_ID_RE = re.compile(r"^(?:chr)?([^:_\-]+)[:_\-](\d+)[_\-](\d+)$")

# bed-ish suffixes stripped when a cell-type label defaults to the filename stem
_BED_SUFFIXES = (
    ".bed.gz", ".narrowPeak.gz", ".broadPeak.gz", ".txt.gz", ".tsv.gz",
    ".bed", ".narrowPeak", ".broadPeak", ".txt", ".tsv", ".gz",
)


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)

    src = p.add_argument_group("peak inputs")
    src.add_argument("--peaks-dir", help="directory of per-cell-type peak BEDs (globbed by --glob)")
    src.add_argument("--glob", default="*.bed*", help="glob within --peaks-dir (default: *.bed*)")
    src.add_argument("--peaks", nargs="+", help="explicit peak BED file(s) (alternative to --peaks-dir)")
    src.add_argument("--bed-has-header", action="store_true", help="peak BEDs carry a header row (coords still read positionally)")
    src.add_argument("--comment-prefix", default="#", help="skip lines starting with this prefix (default: #)")

    cols = p.add_argument_group("peak columns (0-based indices)")
    cols.add_argument("--chrom-idx", type=int, default=0, help="chrom column index (default: 0)")
    cols.add_argument("--start-idx", type=int, default=1, help="start column index (default: 1)")
    cols.add_argument("--end-idx", type=int, default=2, help="end column index (default: 2)")
    cols.add_argument("--score-idx", type=int, help="optional signal/score column index (e.g. 4, or narrowPeak signalValue 6)")
    cols.add_argument("--score-type", help="score_type token (default: 'score' if --score-idx set, else 'presence')")
    cols.add_argument("--min-score", type=float, default=0.0, help="keep scored peaks with value > min-score (default: 0.0; presence rows always kept)")

    ct = p.add_argument_group("cell_type derivation (priority: idx > filename regex > fixed > stem)")
    ct.add_argument("--celltype-idx", type=int, help="per-row cell_type column index in the BED")
    ct.add_argument("--celltype-from-filename", help="regex with a capture group extracting cell_type from the filename")
    ct.add_argument("--cell-type", help="fixed cell_type for all peaks (used when no idx/regex)")

    cond = p.add_argument_group("condition derivation (priority: idx > filename regex > fixed)")
    cond.add_argument("--condition-idx", type=int, help="per-row condition column index in the BED")
    cond.add_argument("--condition-from-filename", help="regex with a capture group extracting AD/control from the filename")
    cond.add_argument("--condition", default="unknown", help="fixed condition when no idx/regex (default: unknown; normalized to AD/control/unknown)")

    links = p.add_argument_group("peak -> gene links (optional; joint pairing)")
    links.add_argument("--links", help="optional peak->gene links TSV")
    links.add_argument("--links-id-col", help="links: peak id column to parse coords from (e.g. chr1:1000-2000)")
    links.add_argument("--links-chrom-col", default="chrom", help="links: chrom column (default: chrom)")
    links.add_argument("--links-start-col", default="start", help="links: start column (default: start)")
    links.add_argument("--links-end-col", default="end", help="links: end column (default: end)")
    links.add_argument("--links-gene-col", default="gene_name", help="links: gene symbol column (default: gene_name)")
    links.add_argument("--links-geneid-col", default="gene_id", help="links: Ensembl gene id column (default: gene_id)")

    meta = p.add_argument_group("constants / metadata")
    meta.add_argument("--tissue", default=TISSUE, help=f"tissue override (default: {TISSUE})")
    meta.add_argument("--pfc-uberon", action="store_true", help="set uberon_id=UBERON:0000451 (prefrontal cortex)")

    out = p.add_argument_group("output / staging")
    out.add_argument("--output", help=f"output .tsv.gz path (default: ./{DATASET}.tsv.gz)")
    out.add_argument("--stage", action="store_true", help="upload .tsv.gz + .tbi to BOTH GCS buckets (OFF by default)")
    out.add_argument("--gcs-finngen", default=GCS_FINNGEN, help="finngen GCS destination")
    out.add_argument("--gcs-daly", default=GCS_DALY, help="daly GCS destination")

    p.add_argument("--sample", "--dry-run", dest="sample", action="store_true",
                   help="run the transform on tiny synthetic input, print first rows, validate tabix; no upload")
    return p.parse_args()


def _numeric_chrom(expr: pl.Expr) -> pl.Expr:
    """Strip any 'chr' prefix and map to a numeric-string chrom (X=23, Y=24, M/MT=25)."""
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
# (alt/random/scaffold/unplaced/Un_*) must be DROPPED here or they break the BigQuery chr load.
CANONICAL_CHROMS = frozenset(str(c) for c in range(1, 26))


def _drop_noncanonical(df: pl.DataFrame) -> pl.DataFrame:
    """Keep only rows whose (already numeric) chrom is a canonical primary chromosome."""
    return df.filter(pl.col("chrom").is_in(CANONICAL_CHROMS))


def _coords_from_id(df: pl.DataFrame, id_col: str) -> pl.DataFrame:
    """Parse chrom/start/end from a single id column like chr1:1000-2000 (chrom -> numeric)."""
    parsed = df.select(
        pl.col(id_col).str.extract_groups(_ID_RE.pattern).alias("_g")
    ).unnest("_g")
    return _drop_noncanonical(df.with_columns(
        _numeric_chrom(parsed["1"]).alias("chrom"),
        parsed["2"].cast(pl.Int64).alias("start"),
        parsed["3"].cast(pl.Int64).alias("end"),
    ))


def _normalize_coords(df: pl.DataFrame, chrom_col: str, start_col: str, end_col: str) -> pl.DataFrame:
    """Rename explicit coordinate columns to chrom/start/end and enforce numeric chrom + int coords."""
    df = df.rename({chrom_col: "chrom", start_col: "start", end_col: "end"})
    return _drop_noncanonical(df.with_columns(
        _numeric_chrom(pl.col("chrom")).alias("chrom"),
        pl.col("start").cast(pl.Int64),
        pl.col("end").cast(pl.Int64),
    ))


def _peak_key() -> pl.Expr:
    return pl.format("{}-{}-{}", pl.col("chrom"), pl.col("start"), pl.col("end"))


def _normalize_condition(col: pl.Expr) -> pl.Expr:
    """Map free-text condition labels to the canonical {AD, control, unknown} axis."""
    c = col.cast(pl.Utf8).str.to_lowercase().str.strip_chars()
    ad = ["ad", "alzheimer", "alzheimers", "alzheimer's", "case", "disease", "dementia", "l[oad]", "load"]
    control = ["control", "ctrl", "normal", "healthy", "non-ad", "nonad", "no-ad",
               "cognitively normal", "nci", "non_ad", "ctl"]
    return (
        pl.when(c.is_in(ad)).then(pl.lit("AD"))
        .when(c.is_in(control)).then(pl.lit("control"))
        .otherwise(pl.lit("unknown"))
    )


def _from_filename(path: str, regex: str | None) -> str | None:
    """Extract a label from a filename: regex capture group if given, else the stripped stem."""
    name = Path(path).name
    if regex:
        m = re.search(regex, name)
        if not m:
            return None
        return m.group(1) if m.groups() else m.group(0)
    stem = name
    for suf in _BED_SUFFIXES:
        if stem.endswith(suf):
            stem = stem[: -len(suf)]
            break
    return stem


def _resolve_peak_files(args: argparse.Namespace) -> list[str]:
    if args.peaks_dir:
        files = sorted(globmod.glob(str(Path(args.peaks_dir) / args.glob)))
        if not files:
            raise SystemExit(f"no files matched {args.glob!r} under {args.peaks_dir}")
        return files
    if args.peaks:
        return list(args.peaks)
    raise SystemExit("provide --peaks-dir or --peaks (or use --sample). See --help.")


def load_peak_file(path: str, args: argparse.Namespace) -> pl.DataFrame:
    """Read one peak BED -> long rows: chrom, start, end, cell_type, _condition_raw, score."""
    df = pl.read_csv(
        path, separator="\t", has_header=args.bed_has_header,
        comment_prefix=args.comment_prefix or None, infer_schema_length=10_000,
        truncate_ragged_lines=True,
    )
    cols = df.columns
    # capture source column names by index BEFORE renaming coordinates
    chrom_c, start_c, end_c = cols[args.chrom_idx], cols[args.start_idx], cols[args.end_idx]
    score_c = cols[args.score_idx] if args.score_idx is not None else None
    celltype_c = cols[args.celltype_idx] if args.celltype_idx is not None else None
    condition_c = cols[args.condition_idx] if args.condition_idx is not None else None

    df = _normalize_coords(df, chrom_c, start_c, end_c)

    if score_c is not None:
        df = df.with_columns(pl.col(score_c).cast(pl.Float64, strict=False).alias("score"))
    else:
        df = df.with_columns(pl.lit(None, dtype=pl.Float64).alias("score"))

    # cell_type: per-row column > filename regex > fixed > filename stem
    if celltype_c is not None:
        df = df.with_columns(pl.col(celltype_c).cast(pl.Utf8).alias("cell_type"))
    else:
        ct = args.cell_type if args.cell_type else _from_filename(path, args.celltype_from_filename)
        if ct is None:
            raise SystemExit(f"could not derive cell_type for {path} (check --celltype-from-filename/--cell-type)")
        df = df.with_columns(pl.lit(ct).alias("cell_type"))

    # condition: per-row column > filename regex > fixed
    if condition_c is not None:
        df = df.with_columns(pl.col(condition_c).cast(pl.Utf8).alias("_condition_raw"))
    else:
        cond_raw = _from_filename(path, args.condition_from_filename) if args.condition_from_filename else args.condition
        df = df.with_columns(pl.lit(cond_raw).alias("_condition_raw"))

    return df.select("chrom", "start", "end", "cell_type", "_condition_raw", "score")


def load_peaks(args: argparse.Namespace) -> pl.DataFrame:
    """Load and concatenate every peak BED into one long table, applying the score filter."""
    files = _resolve_peak_files(args)
    frames = [load_peak_file(f, args) for f in files]
    long = pl.concat(frames, how="vertical_relaxed")
    # keep presence rows (null score) and scored rows above the threshold
    long = long.filter(pl.col("score").is_null() | (pl.col("score") > args.min_score))
    return long


def load_gene_links(args: argparse.Namespace) -> pl.DataFrame:
    """Read peak->gene links and aggregate to one comma-joined gene list per peak key.

    Symbol and id are paired JOINTLY (never independently sorted): a (symbol, id) struct is built
    per link, deduped and sorted as a unit, so the i-th target_gene token always corresponds to
    the i-th target_gene_id token.
    """
    df = pl.read_csv(args.links, separator="\t", infer_schema_length=10_000)
    if args.links_id_col:
        df = _coords_from_id(df, args.links_id_col)
    else:
        df = _normalize_coords(df, args.links_chrom_col, args.links_start_col, args.links_end_col)

    df = df.with_columns(_peak_key().alias("_key"))
    gene = pl.col(args.links_gene_col).cast(pl.Utf8) if args.links_gene_col in df.columns else pl.lit(None, dtype=pl.Utf8)
    geneid = pl.col(args.links_geneid_col).cast(pl.Utf8) if args.links_geneid_col in df.columns else pl.lit(None, dtype=pl.Utf8)
    # normalize nulls to empty strings so a missing symbol or id still keeps its pair and emits an
    # empty token on the missing side, keeping the two comma-joined lists aligned
    df = df.with_columns(
        gene.fill_null("").alias("_gene"),
        geneid.fill_null("").alias("_geneid"),
    )
    agg = (
        df.group_by("_key")
        .agg(pl.struct("_gene", "_geneid").unique().sort().alias("_pairs"))
        .with_columns(
            pl.col("_pairs").list.eval(pl.element().struct.field("_gene")).list.join(",").alias("target_gene"),
            pl.col("_pairs").list.eval(pl.element().struct.field("_geneid")).list.join(",").alias("target_gene_id"),
        )
        .drop("_pairs")
    )
    # empty-cell -> NA convention: a joined list with no real token (only commas/whitespace, e.g.
    # when the source provides gene symbols but no Ensembl ids) becomes NA, not "" / "," / ",,".
    # applied symmetrically; rows that DO have real tokens keep their positional pairing intact.
    return agg.with_columns(
        *[
            pl.when(pl.col(c).str.replace_all(r"[,\s]", "").str.len_chars() == 0)
            .then(pl.lit(None, dtype=pl.Utf8))
            .otherwise(pl.col(c))
            .alias(c)
            for c in ("target_gene", "target_gene_id")
        ]
    )


def build_output(long: pl.DataFrame, args: argparse.Namespace) -> pl.DataFrame:
    """Assemble the 18 canonical columns in order from the long (peak, cell_type) table."""
    score_type = args.score_type or ("score" if args.score_idx is not None else "presence")

    df = long.with_columns(
        _peak_key().alias("peak_id"),
        _normalize_condition(pl.col("_condition_raw")).alias("condition"),
    )
    # presence semantics: never emit a numeric score under a presence score_type
    if score_type == "presence":
        df = df.with_columns(pl.lit(None, dtype=pl.Float64).alias("score"))

    if args.links:
        links = load_gene_links(args)
        df = df.with_columns(_peak_key().alias("_key")).join(links, on="_key", how="left").drop("_key")
    else:
        df = df.with_columns(
            pl.lit(None, dtype=pl.Utf8).alias("target_gene"),
            pl.lit(None, dtype=pl.Utf8).alias("target_gene_id"),
        )

    df = df.with_columns(
        pl.lit(DATASET).alias("dataset"),
        pl.lit(args.tissue).alias("tissue"),
        pl.lit(LIFE_STAGE).alias("life_stage"),
        pl.lit(ASSAY).alias("assay"),
        pl.lit(score_type).alias("score_type"),
        pl.lit(None, dtype=pl.Utf8).alias("n_cells"),
        pl.lit(None, dtype=pl.Utf8).alias("cell_ontology_id"),
        pl.lit("UBERON:0000451" if args.pfc_uberon else None, dtype=pl.Utf8).alias("uberon_id"),
        pl.lit(VERSION).alias("version"),
    )
    return df.select(OUTPUT_COLUMNS)


def write_open_chromatin(df: pl.DataFrame, output_path: str) -> None:
    """sort -k1,1 -k2,2n -> bgzip -> tabix -p bed (interval index), missing values as "NA"."""
    tmpdir = tempfile.mkdtemp()
    body = Path(tmpdir) / "body.tsv"
    # every empty/missing cell serialized as the literal "NA"
    df.write_csv(body, separator="\t", include_header=False, null_value="NA")

    header = "#" + "\t".join(OUTPUT_COLUMNS)
    # prepend header AFTER sorting the body (header must stay on top, not be sorted)
    pipeline = (
        f'( printf "%s\\n" {shell_quote(header)}; '
        f'LC_ALL=C sort -k1,1 -k2,2n {shell_quote(str(body))} ) | bgzip -c > {shell_quote(output_path)}'
    )
    subprocess.run(pipeline, shell=True, check=True, executable="/bin/bash")
    subprocess.run(["tabix", "-f", "-p", "bed", output_path], check=True)
    print(f"  wrote {df.height} rows -> {output_path}")
    print(f"  indexed {output_path}.tbi (tabix -p bed / -s1 -b2 -e3)")


def shell_quote(s: str) -> str:
    return "'" + s.replace("'", "'\\''") + "'"


def upload_to_gcs(local_path: str, gcs_path: str) -> None:
    subprocess.run(["gcloud", "storage", "cp", local_path, gcs_path], check=True)
    print(f"  uploaded {gcs_path}")
    tbi = local_path + ".tbi"
    if Path(tbi).exists():
        subprocess.run(["gcloud", "storage", "cp", tbi, gcs_path + ".tbi"], check=True)
        print(f"  uploaded {gcs_path}.tbi")


def _synthetic_inputs(tmpdir: Path) -> tuple[str, str]:
    """Write tiny synthetic per-cell-type/condition BEDs + a multi-gene peak->gene links file.

    Filenames encode "<cell_type>.<condition>.bed" so the sample exercises both
    --celltype-from-filename and --condition-from-filename. BEDs are headerless
    chrom/start/end/name/score (score at index 4).
    """
    peaks = tmpdir / "peaks"
    peaks.mkdir()
    (peaks / "microglia.AD.bed").write_text(
        "chr1\t100100\t100600\tp1\t50\n"
        "chr1\t200200\t200700\tp2\t80\n"
    )
    (peaks / "microglia.control.bed").write_text(
        "chr1\t200200\t200700\tp2\t65\n"
        "chr2\t50050\t50550\tp3\t40\n"
    )
    # "OPCs.union.bed": union peaks (no AD/control split) -> condition normalizes to "unknown".
    # the second row is a non-canonical hg38 contig (alt/random) that MUST be dropped by Fix A.
    (peaks / "OPCs.union.bed").write_text(
        "chrX\t900900\t901400\tp4\t12\n"
        "chr1_KI270706v1_random\t300300\t300800\tp5\t30\n"
    )
    links = tmpdir / "peak_gene_links.tsv"
    links.write_text(
        "chrom\tstart\tend\tgene_name\tgene_id\n"
        "chr1\t200200\t200700\tSAMD11\tENSG00000187634\n"
        "chr1\t200200\t200700\tNOC2L\tENSG00000188976\n"
        "chr2\t50050\t50550\tSOX11\tENSG00000176887\n"
        # symbols but NO Ensembl ids -> target_gene keeps the symbols, target_gene_id -> NA
        "chr1\t100100\t100600\tFOXA1\t\n"
        "chr1\t100100\t100600\tGATA3\t\n"
    )
    return str(peaks), str(links)


def run_sample() -> None:
    print("=== SAMPLE / DRY-RUN: synthetic input, no GCS upload ===")
    tmpdir = Path(tempfile.mkdtemp())
    peaks_dir, links = _synthetic_inputs(tmpdir)

    args = parse_args()
    args.peaks_dir = peaks_dir
    args.peaks = None
    args.glob = "*.bed"
    args.bed_has_header = False
    args.comment_prefix = "#"
    args.chrom_idx, args.start_idx, args.end_idx = 0, 1, 2
    args.score_idx = 4
    args.score_type = None  # -> resolves to "score" since score_idx is set
    args.min_score = 0.0
    args.celltype_idx = None
    args.celltype_from_filename = r"^([^.]+)\."
    args.cell_type = None
    args.condition_idx = None
    # capture AD/control from the filename; "union" files fall through to "unknown"
    args.condition_from_filename = r"\.(AD|control|union)\."
    args.condition = "unknown"
    args.links = links
    args.links_id_col = None
    args.links_chrom_col, args.links_start_col, args.links_end_col = "chrom", "start", "end"
    args.links_gene_col, args.links_geneid_col = "gene_name", "gene_id"
    args.tissue = TISSUE
    args.pfc_uberon = False

    long = load_peaks(args)
    out = build_output(long, args)

    assert out.columns == OUTPUT_COLUMNS, f"column order mismatch: {out.columns}"
    assert len(out.columns) == 18, f"expected 18 columns, got {len(out.columns)}"
    print(f"  output has {len(out.columns)} columns in canonical order: OK")

    conds = set(out["condition"].unique().to_list())
    assert conds <= {"AD", "control", "unknown"}, f"unexpected condition values: {conds}"
    assert {"AD", "control", "unknown"} <= conds, f"sample should exercise AD/control/unknown; got {conds}"
    print(f"  condition axis derived: {sorted(conds)} (AD/control from filename, union -> unknown): OK")

    cts = sorted(set(out["cell_type"].unique().to_list()))
    print(f"  cell_type labels derived from filename: {cts}")

    # chrom / peak_id are numeric with chrX -> 23
    chroms = set(out["chrom"].to_list())
    assert chroms <= {"1", "2", "23"}, chroms
    assert "23" in chroms, f"chrX should map to 23; got {chroms}"
    xrow = out.filter(pl.col("chrom") == "23").row(0, named=True)
    assert xrow["peak_id"] == "23-900900-901400", xrow["peak_id"]
    print(f"  numeric chrom (chrX -> 23): OK  chroms={sorted(chroms)}  peak_id={xrow['peak_id']}")

    # non-canonical hg38 contigs (alt/random/scaffold/Un) are DROPPED (Fix A) or they break the
    # INT64 chr load; the "chr1_KI270706v1_random" synthetic peak must not survive while chrX does
    assert not any("KI270706" in p for p in out["peak_id"].to_list()), "non-canonical contig leaked"
    assert chroms == {"1", "2", "23"}, f"expected canonical chroms only; got {chroms}"
    print("  non-canonical contig (chr1_KI270706v1_random) dropped: OK")

    # verify joint (positional) gene/id pairing on the multi-gene peak
    row = out.filter((pl.col("peak_id") == "1-200200-200700") & (pl.col("condition") == "AD")).row(0, named=True)
    genes = row["target_gene"].split(",")
    ids = row["target_gene_id"].split(",")
    assert len(genes) == len(ids) == 2, f"expected 2 paired genes, got {genes} / {ids}"
    pairs = dict(zip(genes, ids))
    assert pairs["SAMD11"] == "ENSG00000187634" and pairs["NOC2L"] == "ENSG00000188976", f"mispaired: {pairs}"
    print(f"  peak->gene joint pairing verified: {list(pairs.items())}")

    # symbols-but-no-ids link: target_gene keeps the symbols, target_gene_id renders as NA (fix)
    sym_only = out.filter(pl.col("peak_id") == "1-100100-100600").row(0, named=True)
    assert sym_only["target_gene"] == "FOXA1,GATA3", sym_only["target_gene"]
    assert sym_only["target_gene_id"] is None, sym_only["target_gene_id"]
    print(f"  symbols-without-ids -> target_gene={sym_only['target_gene']!r}, target_gene_id=NA: OK")

    print("\nfirst output rows:")
    with pl.Config(tbl_cols=-1, tbl_width_chars=260, fmt_str_lengths=60):
        print(out.sort("chrom", "start", "cell_type", "condition").head(10))

    out_path = str(tmpdir / f"{DATASET}.sample.tsv.gz")
    write_open_chromatin(out, out_path)

    # empty cells serialize as the literal "NA" (n_cells/cell_ontology_id/uberon_id and, on the
    # union peak, target_gene*); every row still has 18 columns
    import gzip
    with gzip.open(out_path, "rt") as fh:
        body_lines = [ln for ln in fh if not ln.startswith("#")]
    assert all(len(ln.rstrip("\n").split("\t")) == 18 for ln in body_lines), "not 18 columns"
    assert any("\tNA\t" in ln or ln.rstrip("\n").endswith("\tNA") for ln in body_lines), "no NA emitted"
    print(f"  empty cells rendered as NA, 18 columns per row: OK ({len(body_lines)} rows)")

    # the symbols-without-ids peak serializes target_gene_id as the literal NA (not "" or ",")
    na_fields = next(ln.rstrip("\n").split("\t") for ln in body_lines if "1-100100-100600" in ln)
    assert na_fields[15] == "FOXA1,GATA3" and na_fields[16] == "NA", na_fields[15:17]
    print(f"  serialized symbols/NA pair: target_gene={na_fields[15]!r}, target_gene_id={na_fields[16]!r}: OK")

    q = subprocess.run(["tabix", out_path, "1:200300-200400"], capture_output=True, text=True, check=True)
    print("\ntabix overlap query 1:200300-200400 ->")
    print(q.stdout.rstrip() or "  (no rows)")
    assert q.stdout.strip(), "interval query returned nothing — indexing is wrong"
    # peak 200200-200700 must be returned for both AD and control (overlap, not point)
    assert q.stdout.count("1-200200-200700") >= 2, "interval index missed a peak containing pos"
    qx = subprocess.run(["tabix", out_path, "23:901000-901100"], capture_output=True, text=True, check=True)
    print("tabix overlap query 23:901000-901100 (chrX) ->")
    print(qx.stdout.rstrip() or "  (no rows)")
    assert qx.stdout.strip(), "chrX (23) interval query returned nothing — indexing is wrong"
    print("\n=== SAMPLE OK (no upload performed) ===")


def main() -> None:
    args = parse_args()
    if args.sample:
        run_sample()
        return

    output_path = args.output or f"./{DATASET}.tsv.gz"
    print("Loading ROSMAP peak BEDs ...")
    long = load_peaks(args)
    print(f"  {long.height} (peak, cell_type) accessible entries")

    out = build_output(long, args)
    assert out.columns == OUTPUT_COLUMNS, f"column order mismatch: {out.columns}"

    write_open_chromatin(out, output_path)

    if args.stage:
        print("Staging to GCS (both buckets) ...")
        upload_to_gcs(output_path, args.gcs_finngen)
        upload_to_gcs(output_path, args.gcs_daly)
    else:
        print("  --stage not set: skipping GCS upload")
    print("Done.")


if __name__ == "__main__":
    main()
