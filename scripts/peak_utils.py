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

import os
import subprocess
import tempfile
from pathlib import Path
from typing import NamedTuple

import polars as pl

from sumstat_utils import upload_to_gcs  # noqa: F401  re-exported: peak munges stage the same way


CANONICAL_CHROMS = frozenset(str(c) for c in range(1, 26))  # 1..22, X=23, Y=24, M/MT=25


class TabixIndex(NamedTuple):
    """A tabix invocation plus the wording used to report it, so the two cannot drift apart."""

    args: list[str]
    label: str


# the two output shapes. open_chromatin rows are intervals (start, end); variant-level products
# index one position twice, and the API's point lookups depend on which of these was written.
INTERVAL_INDEX = TabixIndex(["-p", "bed"], "INTERVAL: tabix -p bed / -s1 -b2 -e3")
POINT_INDEX = TabixIndex(["-s", "1", "-b", "2", "-e", "2"], "POINT: tabix -s1 -b2 -e2")


def numeric_chrom_expr(col: str | pl.Expr | pl.Series = "chrom") -> pl.Expr:
    """Numeric chromosome token from a seqname: strip any 'chr', then X=23, Y=24, M/MT=25, else the
    contig as-is. Mirrors CHR_STRING_TO_INT_SQL / chrom_to_int() in genetics-results-db so the tabix
    seqnames match the tables' chr INT64 encoding.

    Result stays a string so nulls and non-canonical contigs survive untouched to `filter_canonical`,
    which is the only place they are dropped.
    """
    expr = pl.col(col) if isinstance(col, str) else col
    # the strip is case-insensitive to match CHR_STRING_TO_INT_SQL, so a seqname kept here is one the
    # BigQuery load also accepts. a case-sensitive strip drops an unbounded class instead: ANY
    # case-variant `chr` prefix on an otherwise canonical chromosome falls through to
    # filter_canonical and is discarded as if it were a scaffold. wrong only if a source used case
    # to distinguish two different contigs, which no reference assembly does.
    base = expr.cast(pl.Utf8).str.replace("(?i)^chr", "").str.to_uppercase()
    return (
        pl.when(base == "X").then(pl.lit("23"))
        .when(base == "Y").then(pl.lit("24"))
        .when(base == "M").then(pl.lit("25"))
        .when(base == "MT").then(pl.lit("25"))
        .otherwise(base)
    )


def filter_canonical(df, chrom: str | pl.Expr = "chrom"):
    """Keep only rows whose (already numeric) chrom is a canonical primary chromosome.

    alt/random/scaffold/Un contigs would fail the BigQuery chr INT64 load, and a tabix file
    carrying them indexes fine while the loaded table silently lacks those rows.
    """
    expr = pl.col(chrom) if isinstance(chrom, str) else chrom
    return df.filter(expr.is_in(CANONICAL_CHROMS))


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
