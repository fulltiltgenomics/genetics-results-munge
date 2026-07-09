#!/usr/bin/env python3
"""Munge Li et al. 2023 Science human brain snATAC atlas into the canonical `open_chromatin` TSV.

Source (public):
  - Paper: Li et al., "A comparative atlas of single-cell chromatin accessibility in the
    human brain", Science 2023, 382(6667):eadf7044. doi:10.1126/science.adf7044
  - 1.1M cells, 42 brain regions, 107 cell types, 544,735 cCREs, genome build hg38.
  - GEO: GSE244618 (raw). Processed data (cCRE BEDs, per-cell-type accessibility, cCRE->gene
    links, bigwig pileups) on the CATlas web portal: https://catlas.org/humanbrain/
  - NEMO archive holds demultiplexed data.

>>> TO-VERIFY (format assumptions) <<<
The CATlas portal is a JavaScript app; the exact download URLs and per-file column layouts
could not be scraped programmatically. The two inputs this script consumes and the column
layouts it assumes are documented below. The user MUST verify these against the real files
and adjust the --*-col flags (no code change needed for column-name differences).

INPUT 1 — per-cell-type cCRE accessibility matrix  (--ccre-matrix, REQUIRED)
  Assumed: a TSV where each row is one cCRE and there is one column per brain cell type
  holding an accessibility value (CPM = counts/signal per million, as the portal states
  pileups are "normalized signal per million reads"). cCRE coordinates identified either by
  three explicit columns (default names chrom/start/end, hg38, BED 0-based half-open) OR by a
  single id column like "chr1:1000-2000" / "chr1-1000-2000" / "chr1_1000_2000" (--id-col).
  Any column that is not a coordinate/id column is treated as a cell-type value column.
  ALTERNATIVE the user may hit instead: a binary presence matrix (0/1) or a union cCRE BED
  with a separate membership table. For a binary/presence matrix pass --score-type presence
  (score is emitted empty and rows with value>0 are kept).

INPUT 2 — cCRE -> gene links  (--ccre-gene-links, OPTIONAL)
  Assumed: a TSV linking each cCRE to one or more target genes. Columns (configurable):
  a cCRE id/coords, a gene symbol, and an Ensembl gene id. Multiple genes per cCRE are
  aggregated into a comma-separated list (both symbol and ENSG), documented so the user can
  switch to row-explosion if the downstream schema prefers one gene per row.

INPUT 3 — cell-type metadata  (--cell-metadata, OPTIONAL)
  Assumed: a TSV mapping cell_type -> supporting cell count, used to fill `n_cells`.
  Absent -> n_cells left empty.

Genome build: hg38 (per the paper and the epic design) — NO liftOver.

OUTPUT — canonical `open_chromatin` long TSV, 18 columns IN THIS ORDER (tab-separated):
  chrom, start, end, peak_id, dataset, cell_type, tissue, life_stage, condition, assay,
  score, score_type, n_cells, cell_ontology_id, uberon_id, target_gene, target_gene_id, version
The header's first token is written as `#chrom` (comment-prefixed) following this repo's
convention (cf. munge_hpa.py `#dataset`, sumstats `#chr`): it makes the header a tabix meta
line so `tabix -p bed` skips it while the logical column name stays `chrom`.

Field rules for Li 2023 (see task genetics-results-suite-bzl.14):
  chrom            NUMERIC, no "chr" prefix (1..22, X->23, Y->24, M/MT->25); REQUIRED by the
                   tabix indexing contract.
  start,end        cCRE interval, hg38, BED 0-based half-open (kept verbatim from source).
  peak_id          f"{chrom}-{start}-{end}" with the numeric chrom (e.g. "23-100-200").
  dataset          "li_brain_open_chromatin"  (drives resource mapping to li_brain_atac).
  cell_type        verbatim source brain cell-type label (the matrix column name).
  tissue           "brain".
  life_stage       "adult"   (the atlas is from three adult donors).
  condition        "unknown" (default for atlases).
  assay            "snATAC".
  score/score_type CPM accessibility value + "cpm"  (or empty + "presence" for a binary matrix).
  n_cells          from --cell-metadata if given, else empty.
  cell_ontology_id empty (not trivially mappable here).
  uberon_id        empty.
  target_gene(_id) from --ccre-gene-links when present, else empty.
  version          "2023".

INDEXING CONTRACT (from the epic design, surfaced by the results-api open_chromatin review):
  open_chromatin files MUST be INTERVAL-indexed: `tabix -p bed` (i.e. -s1 -b2 -e3, distinct
  start/end columns). Do NOT point-index (-b2 -e2) or the API variant-overlap fast path would
  silently miss peaks whose interval contains pos. Pipeline: sort -k1,1 -k2,2n -> bgzip -> tabix.

STAGING (off by default): with --stage the .tsv.gz and .tsv.gz.tbi are uploaded to BOTH
  gs://finngen-commons/results_api_data/open_chromatin/li_brain_atac/li_brain_open_chromatin.tsv.gz
  gs://daly-genetics-results/open_chromatin/li_brain_atac/li_brain_open_chromatin.tsv.gz
  These paths must match the results-api profile config files. Without --stage nothing is uploaded.

Local validation without the full dataset or any upload:
  python3 scripts/munge_li_brain.py --sample
"""

import argparse
import re
import subprocess
import tempfile
from pathlib import Path

import polars as pl

DATASET = "li_brain_open_chromatin"
RESOURCE = "li_brain_atac"
VERSION = "2023"

TISSUE = "brain"
LIFE_STAGE = "adult"
CONDITION = "unknown"
ASSAY = "snATAC"

# canonical open_chromatin column order; first header token is comment-prefixed for tabix
OUTPUT_COLUMNS = [
    "chrom", "start", "end", "peak_id", "dataset", "cell_type", "tissue", "life_stage",
    "condition", "assay", "score", "score_type", "n_cells", "cell_ontology_id",
    "uberon_id", "target_gene", "target_gene_id", "version",
]

GCS_FINNGEN = f"gs://finngen-commons/results_api_data/open_chromatin/{RESOURCE}/{DATASET}.tsv.gz"
GCS_DALY = f"gs://daly-genetics-results/open_chromatin/{RESOURCE}/{DATASET}.tsv.gz"

# "chr1:1000-2000", "chr1-1000-2000", "chr1_1000_2000"; the "chr" prefix is optional and
# gets stripped in _coords_from_id so ids like "1:1000-2000" are not dropped.
_ID_RE = re.compile(r"^(?:chr)?([^:_\-]+)[:_\-](\d+)[_\-](\d+)$")


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


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--ccre-matrix", help="cCRE x cell-type accessibility matrix TSV (see header)")
    p.add_argument("--id-col", help="single cCRE id column to parse coords from (e.g. chr1:1000-2000)")
    p.add_argument("--chrom-col", default="chrom", help="coordinate column: chrom (default: chrom)")
    p.add_argument("--start-col", default="start", help="coordinate column: start (default: start)")
    p.add_argument("--end-col", default="end", help="coordinate column: end (default: end)")
    p.add_argument("--score-type", default="cpm", help="score_type token (default: cpm; use 'presence' for a binary matrix)")
    p.add_argument("--min-score", type=float, default=0.0, help="keep cell-type entries with value > min-score (default: 0.0)")

    p.add_argument("--ccre-gene-links", help="optional cCRE->gene links TSV")
    p.add_argument("--links-id-col", help="links: cCRE id column to parse coords from")
    p.add_argument("--links-chrom-col", default="chrom", help="links: chrom column (default: chrom)")
    p.add_argument("--links-start-col", default="start", help="links: start column (default: start)")
    p.add_argument("--links-end-col", default="end", help="links: end column (default: end)")
    p.add_argument("--links-gene-col", default="gene_name", help="links: gene symbol column (default: gene_name)")
    p.add_argument("--links-geneid-col", default="gene_id", help="links: Ensembl gene id column (default: gene_id)")

    p.add_argument("--cell-metadata", help="optional cell_type -> n_cells TSV")
    p.add_argument("--meta-celltype-col", default="cell_type", help="metadata: cell_type column (default: cell_type)")
    p.add_argument("--meta-ncells-col", default="n_cells", help="metadata: n_cells column (default: n_cells)")

    p.add_argument("--output", help="output .tsv.gz path (default: ./li_brain_open_chromatin.tsv.gz)")
    p.add_argument("--stage", action="store_true", help="upload .tsv.gz + .tbi to BOTH GCS buckets (OFF by default)")
    p.add_argument("--gcs-finngen", default=GCS_FINNGEN, help="finngen GCS destination")
    p.add_argument("--gcs-daly", default=GCS_DALY, help="daly GCS destination")

    p.add_argument("--sample", "--dry-run", dest="sample", action="store_true",
                   help="run the transform on tiny synthetic input, print first rows, validate tabix; no upload")
    return p.parse_args()


def _coords_from_id(df: pl.DataFrame, id_col: str) -> pl.DataFrame:
    """Parse chrom/start/end from a single id column like chr1:1000-2000 (chrom -> numeric)."""
    parsed = df.select(
        pl.col(id_col).str.extract_groups(_ID_RE.pattern).alias("_g")
    ).unnest("_g")
    return df.with_columns(
        _numeric_chrom(parsed["1"]).alias("chrom"),
        parsed["2"].cast(pl.Int64).alias("start"),
        parsed["3"].cast(pl.Int64).alias("end"),
    )


def _normalize_coords(df: pl.DataFrame, chrom_col: str, start_col: str, end_col: str) -> pl.DataFrame:
    """Rename explicit coordinate columns to chrom/start/end and enforce numeric chrom + int coords."""
    df = df.rename({chrom_col: "chrom", start_col: "start", end_col: "end"})
    return df.with_columns(
        _numeric_chrom(pl.col("chrom")).alias("chrom"),
        pl.col("start").cast(pl.Int64),
        pl.col("end").cast(pl.Int64),
    )


def load_matrix(args: argparse.Namespace) -> pl.DataFrame:
    """Read the cCRE x cell-type matrix and unpivot to long: one row per (cCRE, cell_type)."""
    df = pl.read_csv(args.ccre_matrix, separator="\t", infer_schema_length=10_000)

    if args.id_col:
        df = _coords_from_id(df, args.id_col)
        coord_source_cols = [args.id_col]
    else:
        df = _normalize_coords(df, args.chrom_col, args.start_col, args.end_col)
        coord_source_cols = [args.chrom_col, args.start_col, args.end_col]

    coord_source_cols = [c for c in coord_source_cols if c in df.columns]
    cell_cols = [c for c in df.columns if c not in {"chrom", "start", "end", *coord_source_cols}]
    if not cell_cols:
        raise ValueError("no cell-type value columns found in --ccre-matrix (check coordinate/id flags)")

    long = df.unpivot(
        index=["chrom", "start", "end"], on=cell_cols,
        variable_name="cell_type", value_name="_value",
    )
    long = long.with_columns(pl.col("_value").cast(pl.Float64, strict=False))
    long = long.filter(pl.col("_value") > args.min_score)

    if args.score_type == "presence":
        long = long.with_columns(pl.lit(None, dtype=pl.Utf8).alias("score"))
    else:
        long = long.with_columns(pl.col("_value").alias("score"))
    return long.drop("_value")


def load_gene_links(args: argparse.Namespace) -> pl.DataFrame:
    """Read cCRE->gene links and aggregate to one comma-joined gene list per cCRE key."""
    df = pl.read_csv(args.ccre_gene_links, separator="\t", infer_schema_length=10_000)
    if args.links_id_col:
        df = _coords_from_id(df, args.links_id_col)
    else:
        df = _normalize_coords(df, args.links_chrom_col, args.links_start_col, args.links_end_col)

    df = df.with_columns(_peak_key().alias("_key"))
    gene = pl.col(args.links_gene_col).cast(pl.Utf8) if args.links_gene_col in df.columns else pl.lit(None, dtype=pl.Utf8)
    geneid = pl.col(args.links_geneid_col).cast(pl.Utf8) if args.links_geneid_col in df.columns else pl.lit(None, dtype=pl.Utf8)
    # normalize nulls to empty strings so a missing symbol or id still keeps its pair and
    # emits an empty token on the missing side, keeping the two comma-joined lists aligned
    df = df.with_columns(
        gene.fill_null("").alias("_gene"),
        geneid.fill_null("").alias("_geneid"),
    )

    # aggregate the (symbol, id) pairs JOINTLY: build a struct per link, dedup the pairs, and
    # sort the structs (lexicographic by _gene first, then _geneid) so the i-th token of
    # target_gene always corresponds to the i-th token of target_gene_id
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


def _peak_key() -> pl.Expr:
    return pl.format("{}-{}-{}", pl.col("chrom"), pl.col("start"), pl.col("end"))


def build_output(long: pl.DataFrame, args: argparse.Namespace) -> pl.DataFrame:
    """Assemble the 18 canonical columns in order from the long (cCRE, cell_type) table."""
    df = long.with_columns(_peak_key().alias("peak_id"))

    if args.ccre_gene_links:
        links = load_gene_links(args)
        df = df.with_columns(_peak_key().alias("_key")).join(links, on="_key", how="left").drop("_key")
    else:
        df = df.with_columns(
            pl.lit(None, dtype=pl.Utf8).alias("target_gene"),
            pl.lit(None, dtype=pl.Utf8).alias("target_gene_id"),
        )

    if args.cell_metadata:
        meta = pl.read_csv(args.cell_metadata, separator="\t", infer_schema_length=10_000)
        meta = meta.select(
            pl.col(args.meta_celltype_col).cast(pl.Utf8).alias("cell_type"),
            pl.col(args.meta_ncells_col).cast(pl.Utf8).alias("n_cells"),
        )
        df = df.join(meta, on="cell_type", how="left")
    else:
        df = df.with_columns(pl.lit(None, dtype=pl.Utf8).alias("n_cells"))

    df = df.with_columns(
        pl.lit(DATASET).alias("dataset"),
        pl.lit(TISSUE).alias("tissue"),
        pl.lit(LIFE_STAGE).alias("life_stage"),
        pl.lit(CONDITION).alias("condition"),
        pl.lit(ASSAY).alias("assay"),
        pl.lit(args.score_type).alias("score_type"),
        pl.lit(None, dtype=pl.Utf8).alias("cell_ontology_id"),
        pl.lit(None, dtype=pl.Utf8).alias("uberon_id"),
        pl.lit(VERSION).alias("version"),
    )
    # note: load_matrix always emits a "score" column (null for presence matrices), so no fallback needed

    return df.select(OUTPUT_COLUMNS)


def write_open_chromatin(df: pl.DataFrame, output_path: str) -> None:
    """sort -k1,1 -k2,2n -> bgzip -> tabix -p bed (interval index), missing values as "NA"."""
    tmpdir = tempfile.mkdtemp()
    body = Path(tmpdir) / "body.tsv"
    # write body without header; every empty/missing cell serialized as the literal "NA"
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


def _synthetic_inputs(tmpdir: Path) -> tuple[str, str, str]:
    """Write a tiny synthetic matrix + links + cell metadata mimicking the assumed layouts."""
    matrix = tmpdir / "ccre_matrix.tsv"
    matrix.write_text(
        "chrom\tstart\tend\tOligodendrocyte\tMicroglia\tAstrocyte\n"
        "chr1\t100100\t100600\t0.0\t3.5\t0.0\n"
        "chr1\t200200\t200700\t7.2\t0.0\t1.1\n"
        "chr2\t50050\t50550\t0.0\t0.0\t4.4\n"
        "chrX\t900900\t901400\t2.0\t0.0\t0.0\n"
    )
    links = tmpdir / "ccre_gene_links.tsv"
    links.write_text(
        "chrom\tstart\tend\tgene_name\tgene_id\n"
        "chr1\t200200\t200700\tSAMD11\tENSG00000187634\n"
        "chr1\t200200\t200700\tNOC2L\tENSG00000188976\n"
        "chr2\t50050\t50550\tSOX11\tENSG00000176887\n"
        # symbols but NO Ensembl ids -> target_gene keeps the symbols, target_gene_id -> NA
        "chr1\t100100\t100600\tFOXA1\t\n"
        "chr1\t100100\t100600\tGATA3\t\n"
    )
    meta = tmpdir / "cell_metadata.tsv"
    meta.write_text(
        "cell_type\tn_cells\n"
        "Oligodendrocyte\t120000\n"
        "Microglia\t35000\n"
        "Astrocyte\t80000\n"
    )
    return str(matrix), str(links), str(meta)


def run_sample() -> None:
    print("=== SAMPLE / DRY-RUN: synthetic input, no GCS upload ===")
    tmpdir = Path(tempfile.mkdtemp())
    matrix, links, meta = _synthetic_inputs(tmpdir)

    args = parse_args()
    args.ccre_matrix = matrix
    args.id_col = None
    args.chrom_col, args.start_col, args.end_col = "chrom", "start", "end"
    args.ccre_gene_links = links
    args.links_id_col = None
    args.links_chrom_col, args.links_start_col, args.links_end_col = "chrom", "start", "end"
    args.links_gene_col, args.links_geneid_col = "gene_name", "gene_id"
    args.cell_metadata = meta
    args.meta_celltype_col, args.meta_ncells_col = "cell_type", "n_cells"
    args.score_type = "cpm"
    args.min_score = 0.0

    long = load_matrix(args)
    out = build_output(long, args)

    assert out.columns == OUTPUT_COLUMNS, f"column order mismatch: {out.columns}"
    print(f"  output has {len(out.columns)} columns in canonical order: OK")

    # chrom / peak_id are numeric with chrX -> 23
    chroms = set(out["chrom"].to_list())
    assert chroms <= {"1", "2", "23"}, chroms
    assert "23" in chroms, f"chrX should map to 23; got {chroms}"
    xrow = out.filter(pl.col("chrom") == "23").row(0, named=True)
    assert xrow["peak_id"] == "23-900900-901400", xrow["peak_id"]
    print(f"  numeric chrom (chrX -> 23): OK  chroms={sorted(chroms)}  peak_id={xrow['peak_id']}")

    # symbols-but-no-ids link: target_gene keeps the symbols, target_gene_id renders as NA (fix)
    sym_only = out.filter(pl.col("peak_id") == "1-100100-100600").row(0, named=True)
    assert sym_only["target_gene"] == "FOXA1,GATA3", sym_only["target_gene"]
    assert sym_only["target_gene_id"] is None, sym_only["target_gene_id"]
    print(f"  symbols-without-ids -> target_gene={sym_only['target_gene']!r}, target_gene_id=NA: OK")

    print("\nfirst output rows:")
    with pl.Config(tbl_cols=-1, tbl_width_chars=240, fmt_str_lengths=60):
        print(out.head())

    out_path = str(tmpdir / f"{DATASET}.sample.tsv.gz")
    write_open_chromatin(out, out_path)

    # empty cells serialize as the literal "NA" (n_cells is present here, but cell_ontology_id/
    # uberon_id are empty -> must render NA)
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

    # validate the interval index answers overlap queries (numeric seqnames, incl. chrX->23)
    q = subprocess.run(["tabix", out_path, "1:200300-200400"], capture_output=True, text=True, check=True)
    print("\ntabix overlap query 1:200300-200400 ->")
    print(q.stdout.rstrip() or "  (no rows)")
    assert q.stdout.strip(), "interval query returned nothing — indexing is wrong"
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

    if not args.ccre_matrix:
        raise SystemExit("--ccre-matrix is required (or use --sample). See --help.")

    output_path = args.output or f"./{DATASET}.tsv.gz"
    print(f"Reading matrix {args.ccre_matrix} ...")
    long = load_matrix(args)
    print(f"  {long.height} (cCRE, cell_type) accessible entries")

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
