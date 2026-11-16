"""Shared output contract for the peak / open-chromatin munging scripts.

What lives here is the half of a peak munge that is *not* about the source paper: the numeric
chromosome encoding the BigQuery loaders and the tabix seqnames both expect, and the
sort -> bgzip -> tabix write. The per-source `load_*` / `transform` / `select()` stay in each
script — a paper's column meanings are not shared facts, and the nearest script is meant to be
readable as that dataset's specification.

Counterpart to `sumstat_utils.py`, which does the same job for the sumstat family. Where a peak
munge stages at all, the upload is byte-for-byte the sumstat one, so `upload_to_gcs` is re-exported
from there rather than restated — but staging is not universal here, and `--stage` is the script's
own flag to have or not.
"""

import argparse
import gzip
import os
import re
import subprocess
import tempfile
from pathlib import Path
from typing import NamedTuple

import polars as pl

from sumstat_utils import CHROM_CODES, chrom_code_expr
from sumstat_utils import upload_to_gcs  # noqa: F401  re-exported: peak munges stage the same way


# the autosomes plus whatever CHROM_CODES codes; 25 is absent because nothing maps to it
CANONICAL_CHROMS = frozenset(str(c) for c in range(1, 23)) | frozenset(CHROM_CODES.values())

# the open_chromatin table's columns, in the order every product of it is written
OPEN_CHROMATIN_COLUMNS = [
    "chrom", "start", "end", "peak_id", "dataset", "cell_type", "tissue", "life_stage",
    "condition", "assay", "score", "score_type", "n_cells", "cell_ontology_id",
    "uberon_id", "target_gene", "target_gene_id", "version",
]

# "chr1:1000-2000", "chr1-1000-2000", "chr1_1000_2000"; the "chr" prefix is optional, so ids like
# "1:1000-2000" are not dropped
_INTERVAL_ID_RE = re.compile(r"^(?:chr)?([^:_\-]+)[:_\-](\d+)[_\-](\d+)$")


def open_chromatin_gcs_paths(resource: str, dataset: str) -> tuple[str, str]:
    """(finngen, daly) staging destinations of an open_chromatin product."""
    return (
        f"gs://finngen-commons/results_api_data/open_chromatin/{resource}/{dataset}.tsv.gz",
        f"gs://daly-genetics-results/open_chromatin/{resource}/{dataset}.tsv.gz",
    )


class TabixIndex(NamedTuple):
    """A tabix invocation plus the wording used to report it, so the two cannot drift apart."""

    args: list[str]
    label: str


# the two output shapes. open_chromatin rows are intervals (start, end); variant-level products
# index one position twice, and the API's point lookups depend on which of these was written.
INTERVAL_INDEX = TabixIndex(["-p", "bed"], "INTERVAL: tabix -p bed / -s1 -b2 -e3")
POINT_INDEX = TabixIndex(["-s", "1", "-b", "2", "-e", "2"], "POINT: tabix -s1 -b2 -e2")


def numeric_chrom_expr(col: str | pl.Expr | pl.Series = "chrom") -> pl.Expr:
    """`sumstat_utils.chrom_code_expr` under the peak family's column default: strip any `chr`,
    then X=23, Y=24, M/MT=26, else the contig as-is, so the tabix seqnames match the tables'
    chr INT64 encoding.

    Result stays a string so nulls and non-canonical contigs survive untouched to `filter_canonical`,
    which is the only place they are dropped.
    """
    return chrom_code_expr(col)


def filter_canonical(df, chrom: str | pl.Expr = "chrom"):
    """Keep only rows whose (already numeric) chrom is a canonical primary chromosome.

    alt/random/scaffold/Un contigs would fail the BigQuery chr INT64 load, and a tabix file
    carrying them indexes fine while the loaded table silently lacks those rows.
    """
    expr = pl.col(chrom) if isinstance(chrom, str) else chrom
    return df.filter(expr.is_in(CANONICAL_CHROMS))


def coords_from_id(df: pl.DataFrame, id_col: str) -> pl.DataFrame:
    """Parse chrom/start/end from a single id column like chr1:1000-2000 (chrom -> numeric)."""
    parsed = df.select(
        pl.col(id_col).str.extract_groups(_INTERVAL_ID_RE.pattern).alias("_g")
    ).unnest("_g")
    return filter_canonical(df.with_columns(
        numeric_chrom_expr(parsed["1"]).alias("chrom"),
        parsed["2"].cast(pl.Int64).alias("start"),
        parsed["3"].cast(pl.Int64).alias("end"),
    ))


def normalize_coords(df: pl.DataFrame, chrom_col: str, start_col: str, end_col: str) -> pl.DataFrame:
    """Rename explicit coordinate columns to chrom/start/end and enforce numeric chrom + int coords."""
    df = df.rename({chrom_col: "chrom", start_col: "start", end_col: "end"})
    return filter_canonical(df.with_columns(
        numeric_chrom_expr().alias("chrom"),
        pl.col("start").cast(pl.Int64),
        pl.col("end").cast(pl.Int64),
    ))


def peak_key() -> pl.Expr:
    return pl.format("{}-{}-{}", pl.col("chrom"), pl.col("start"), pl.col("end"))


def load_gene_links(path: str, args: argparse.Namespace) -> pl.DataFrame:
    """Read peak->gene links and aggregate to one comma-joined gene list per peak key (`_key`).

    The link coordinates come from `args.links_id_col`, or from the three
    `args.links_{chrom,start,end}_col` when that is unset; the genes from `args.links_gene_col`
    and `args.links_geneid_col`, either of which may be absent from the file.
    """
    df = pl.read_csv(path, separator="\t", infer_schema_length=10_000)
    if args.links_id_col:
        df = coords_from_id(df, args.links_id_col)
    else:
        df = normalize_coords(df, args.links_chrom_col, args.links_start_col, args.links_end_col)

    df = df.with_columns(peak_key().alias("_key"))
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


def build_matrix_output(long: pl.DataFrame, args: argparse.Namespace, **constants) -> pl.DataFrame:
    """Assemble the open_chromatin columns of a cCRE-by-cell-type matrix munge from its long
    (cCRE, cell_type) table: `peak_id`, the gene links of `args.ccre_gene_links`, `n_cells` from
    `args.cell_metadata`, then one column per entry of `constants` (None for an all-null column).

    Every other output column, `score` included, has to be on `long` already.
    """
    df = long.with_columns(peak_key().alias("peak_id"))

    if args.ccre_gene_links:
        links = load_gene_links(args.ccre_gene_links, args)
        df = df.with_columns(peak_key().alias("_key")).join(links, on="_key", how="left").drop("_key")
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
        (pl.lit(None, dtype=pl.Utf8) if value is None else pl.lit(value)).alias(name)
        for name, value in constants.items()
    )
    return df.select(OPEN_CHROMATIN_COLUMNS)


def blank_to_null(df: pl.DataFrame, columns: list[str] | None = None) -> pl.DataFrame:
    """Coerce empty-string cells to null so `null_value="NA"` reaches them too.

    The output convention is that no cell is ever the empty string: a blank cell and a missing one
    are the same fact, and only "NA" survives the BigQuery load as NULL.
    """
    columns = columns if columns is not None else df.columns
    return df.with_columns(
        pl.when(pl.col(c).cast(pl.Utf8).str.len_chars() == 0).then(None).otherwise(pl.col(c)).alias(c)
        for c in columns
    )


def sort_numeric(df: pl.DataFrame, columns: list[str]) -> pl.DataFrame:
    """In-memory sort by (numeric chrom, position) — the order tabix requires.

    columns[1] is the position column (start for interval products, pos for point products).
    """
    return df.sort(
        by=[pl.col("chrom").cast(pl.Int64, strict=False), pl.col(columns[1]).cast(pl.Int64, strict=False)],
        nulls_last=True,
    )


def write_body(df: pl.DataFrame, body_path: str | Path, columns: list[str]) -> None:
    """Serialize rows to a headerless on-disk body TSV in `columns` order, missing cells as "NA".

    `columns` is what the header line is later built from, so binding the body to it here is what
    keeps the two in the same order.

    Unlike `append_body` / `write_bgzip_index` this does NOT run `blank_to_null` first: an
    empty-string cell would serialize as an empty field, not "NA". The peak transforms already
    resolve blanks to null themselves (each has its own "empty-cell -> NA" step, because only they
    know which joined-token results are logically empty), and a second full-frame pass over a
    whole-atlas frame here would cost a copy to catch nothing.
    """
    df.select(columns).write_csv(body_path, separator="\t", include_header=False, null_value="NA")


def append_body(df: pl.DataFrame, columns: list[str], body_fh) -> None:
    """Append headerless rows to an open body handle, blanks coerced to null so they serialize "NA"."""
    blank_to_null(df.select(columns), columns).write_csv(
        body_fh, separator="\t", include_header=False, null_value="NA"
    )


def _bgzip_with_header(columns: list[str], output_path: str, write_rows) -> None:
    """bgzip `write_rows` under a '#'-prefixed header line into output_path."""
    header = ("#" + "\t".join(columns) + "\n").encode()
    with open(output_path, "wb") as out_fh:
        proc = subprocess.Popen(["bgzip", "-c"], stdin=subprocess.PIPE, stdout=out_fh)
        assert proc.stdin is not None
        try:
            proc.stdin.write(header)
            write_rows(proc.stdin)
        finally:
            # a raising write_rows must not leave bgzip running on a half-written output
            try:
                if not proc.stdin.closed:
                    proc.stdin.close()
            except BrokenPipeError:
                pass
            rc = proc.wait()
    if rc != 0:
        raise subprocess.CalledProcessError(rc, "bgzip -c")


def _sort_into(stdin, body_path: str, sort_tmp: str) -> None:
    """`LC_ALL=C sort -k1,1 -k2,2n body_path` straight down an already-open pipe."""
    stdin.flush()  # the header sits in Python's buffer; sort writes to the fd underneath it
    proc = subprocess.Popen(
        ["sort", "-T", sort_tmp, "-k1,1", "-k2,2n", body_path],
        env={**os.environ, "LC_ALL": "C"}, stdout=stdin,
    )
    rc = proc.wait()
    if rc != 0:
        raise subprocess.CalledProcessError(rc, f"sort -k1,1 -k2,2n {body_path}")


def external_sort_bgzip_index(
    body_path: str | Path,
    output_path: str,
    columns: list[str],
    index: TabixIndex = INTERVAL_INDEX,
    sort_tmp: str | Path | None = None,
) -> None:
    """External `LC_ALL=C sort -k1,1 -k2,2n` of a body TSV, streamed under the header into bgzip,
    then tabix.

    What it bounds: RAM held by this process, which is a pipe buffer — the row count is sort's
    problem, not the frame's. What it costs: the body has to exist on disk first. sort's output is
    never materialized, so peak disk is the body, plus sort's own spill (nothing while it fits in
    sort's buffer, about one more body once it merges externally), plus the compressed output.
    Landing the sorted stream in a file before compressing it would add a second full uncompressed
    copy — on a tens-of-GB atlas that is the difference between fitting and ENOSPC, and these are
    exactly the callers with no room to spare.

    Both exit codes are checked. The equivalent shell pipeline reported only bgzip's, so a sort that
    died mid-merge — ENOSPC, on the disks this path is for — produced a header-only file that tabix
    then indexed cleanly and the run exited 0.

    `sort_tmp` moves sort's scratch off a small root disk; it defaults to the body's own directory,
    which the caller already chose.
    """
    body_path = str(body_path)
    sort_tmp = str(sort_tmp) if sort_tmp is not None else str(Path(body_path).parent)
    _bgzip_with_header(columns, output_path, lambda stdin: _sort_into(stdin, body_path, sort_tmp))
    subprocess.run(["tabix", "-f"] + index.args + [output_path], check=True)
    print(f"  indexed {output_path}.tbi ({index.label})")


def write_bgzip_index(
    df: pl.DataFrame,
    output_path: str,
    columns: list[str],
    index: TabixIndex = INTERVAL_INDEX,
    sort: bool = True,
) -> None:
    """In-memory sort -> bgzip -> tabix, the frame streamed straight into bgzip.

    What it bounds: disk — nothing uncompressed is ever written, so peak disk is the output alone.
    What it costs: the whole frame must be in RAM. Products that already fit there take this;
    whole-atlas builds take `external_sort_bgzip_index`, which trades a body on disk for the bound.

    `sort=False` is for callers whose transform already emitted rows in (numeric chrom, position)
    order; re-sorting would only risk reordering ties.
    """
    df = blank_to_null(df.select(columns), columns)
    if sort:
        df = sort_numeric(df, columns)
    _bgzip_with_header(
        columns, output_path,
        lambda stdin: df.write_csv(stdin, separator="\t", include_header=False, null_value="NA"),
    )
    subprocess.run(["tabix", "-f"] + index.args + [output_path], check=True)
    print(f"  wrote {df.height} rows -> {output_path}")
    print(f"  indexed {output_path}.tbi ({index.label})")


# both writers below put the body under `mkdtemp()`, and `external_sort_bgzip_index` then defaults
# sort's scratch to that same directory. So one TMPDIR in the environment moves the body AND the
# spill together — which is the knob a caller on a small root disk actually wants, and is why
# neither takes a sort_tmp of its own. Only the callers that build their own body (catlas, epimap,
# marderstein) pass sort_tmp, and they pass it straight to `external_sort_bgzip_index`.
def write_open_chromatin(df: pl.DataFrame, output_path: str, columns: list[str]) -> None:
    """The open_chromatin write: body TSV -> external sort -> bgzip -> interval index."""
    body = Path(tempfile.mkdtemp()) / "body.tsv"
    write_body(df, body, columns)
    external_sort_bgzip_index(body, output_path, columns, INTERVAL_INDEX)
    print(f"  wrote {df.height} rows -> {output_path}")


def write_open_chromatin_lazy(lf: pl.LazyFrame, output_path: str, columns: list[str]) -> None:
    """Streaming counterpart: sink the lazy frame to an on-disk body without materializing it in
    RAM, then take the same sort/bgzip/index path."""
    body = Path(tempfile.mkdtemp()) / "body.tsv"
    lf.select(columns).sink_csv(body, separator="\t", include_header=False, null_value="NA",
                                maintain_order=False, engine="streaming")
    external_sort_bgzip_index(body, output_path, columns, INTERVAL_INDEX)
    print(f"  streamed atlas -> {output_path}")


def stage_or_skip(output_path: str, args: argparse.Namespace) -> None:
    """Upload to both buckets (`args.gcs_finngen`, `args.gcs_daly`) when `args.stage` is set."""
    if args.stage:
        print("Staging to GCS (both buckets) ...")
        upload_to_gcs(output_path, args.gcs_finngen)
        upload_to_gcs(output_path, args.gcs_daly)
    else:
        print("  --stage not set: skipping GCS upload")


# the three checks below are the part of a cCRE-matrix munge's --sample run that does not depend
# on its cell types. They assert on the synthetic cCREs those scripts' `_synthetic_inputs` write
# (chr1:100100, chr1:200200, chrX:900900 and one chr1_KI270706v1_random), so a script whose
# fixture has other peaks cannot use them.
def check_sample_chroms(out: pl.DataFrame) -> None:
    """chrom / peak_id are numeric with chrX -> 23, and the non-canonical contig is gone."""
    chroms = set(out["chrom"].to_list())
    assert chroms <= {"1", "2", "23"}, chroms
    assert "23" in chroms, f"chrX should map to 23; got {chroms}"
    xrow = out.filter(pl.col("chrom") == "23").row(0, named=True)
    assert xrow["peak_id"] == "23-900900-901400", xrow["peak_id"]
    print(f"  numeric chrom (chrX -> 23): OK  chroms={sorted(chroms)}  peak_id={xrow['peak_id']}")

    # a non-canonical hg38 contig breaks the INT64 chr load, so it must not survive while chrX does
    assert not any("KI270706" in p for p in out["peak_id"].to_list()), "non-canonical contig leaked"
    assert chroms == {"1", "2", "23"}, f"expected canonical chroms only; got {chroms}"
    print("  non-canonical contig (chr1_KI270706v1_random) dropped: OK")


def check_sample_symbols_only(out: pl.DataFrame) -> None:
    """A link with symbols but no ids keeps the symbols in target_gene and nulls target_gene_id."""
    sym_only = out.filter(pl.col("peak_id") == "1-100100-100600").row(0, named=True)
    assert sym_only["target_gene"] == "FOXA1,GATA3", sym_only["target_gene"]
    assert sym_only["target_gene_id"] is None, sym_only["target_gene_id"]
    print(f"  symbols-without-ids -> target_gene={sym_only['target_gene']!r}, target_gene_id=NA: OK")


def write_and_check_sample(out: pl.DataFrame, out_path: str, tbl_width_chars: int) -> None:
    """Write the sample output and check what only the file can show: "NA" for empty cells, the
    column count, and that the interval index answers overlap queries."""
    print("\nfirst output rows:")
    with pl.Config(tbl_cols=-1, tbl_width_chars=tbl_width_chars, fmt_str_lengths=60):
        print(out.head())

    write_open_chromatin(out, out_path, OPEN_CHROMATIN_COLUMNS)

    with gzip.open(out_path, "rt") as fh:
        body_lines = [ln for ln in fh if not ln.startswith("#")]
    assert all(len(ln.rstrip("\n").split("\t")) == 18 for ln in body_lines), "not 18 columns"
    assert any("\tNA\t" in ln or ln.rstrip("\n").endswith("\tNA") for ln in body_lines), "no NA emitted"
    print(f"  empty cells rendered as NA, 18 columns per row: OK ({len(body_lines)} rows)")

    # the symbols-without-ids cCRE serializes target_gene_id as the literal NA (not "" or ",")
    na_fields = next(ln.rstrip("\n").split("\t") for ln in body_lines if "1-100100-100600" in ln)
    assert na_fields[15] == "FOXA1,GATA3" and na_fields[16] == "NA", na_fields[15:17]
    print(f"  serialized symbols/NA pair: target_gene={na_fields[15]!r}, target_gene_id={na_fields[16]!r}: OK")

    # numeric seqnames, chrX -> 23 included
    q = subprocess.run(["tabix", out_path, "1:200300-200400"], capture_output=True, text=True, check=True)
    print("\ntabix overlap query 1:200300-200400 ->")
    print(q.stdout.rstrip() or "  (no rows)")
    assert q.stdout.strip(), "interval query returned nothing — indexing is wrong"
    qx = subprocess.run(["tabix", out_path, "23:901000-901100"], capture_output=True, text=True, check=True)
    print("tabix overlap query 23:901000-901100 (chrX) ->")
    print(qx.stdout.rstrip() or "  (no rows)")
    assert qx.stdout.strip(), "chrX (23) interval query returned nothing — indexing is wrong"
    print("\n=== SAMPLE OK (no upload performed) ===")
