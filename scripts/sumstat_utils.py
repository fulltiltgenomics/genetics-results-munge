"""Shared utilities for sumstat munging scripts."""

import subprocess
import tempfile
from pathlib import Path

import polars as pl


GNOMAD_DEFAULT = "gs://finngen-commons/results_api_data/gnomad/gnomad.genomes.exomes.v4.0.sites.v2.tsv.bgz"

GNOMAD_AF_COLS = ["AF", "AF_afr", "AF_amr", "AF_asj", "AF_eas", "AF_fin", "AF_mid", "AF_nfe", "AF_remaining", "AF_sas"]

GNOMAD_KEEP_COLS = ["#chr", "pos", "ref", "alt", "rsids", "genome_or_exome"] + GNOMAD_AF_COLS


def upload_to_gcs(local_path: str, gcs_path: str) -> None:
    """Copy a local file (and its .tbi if present) to GCS."""
    subprocess.run(["gcloud", "storage", "cp", local_path, gcs_path], check=True)
    print(f"  uploaded {gcs_path}")
    tbi_local = local_path + ".tbi"
    if Path(tbi_local).exists():
        subprocess.run(["gcloud", "storage", "cp", tbi_local, gcs_path + ".tbi"], check=True)
        print(f"  uploaded {gcs_path}.tbi")


def write_bgzip(df: pl.DataFrame, local_path: str) -> None:
    """Write DataFrame as bgzipped TSV with tabix index."""
    with subprocess.Popen(
        ["bgzip", "-c"],
        stdin=subprocess.PIPE,
        stdout=open(local_path, "wb"),
    ) as proc:
        df.write_csv(proc.stdin, separator="\t", null_value="NA")
    subprocess.run(["tabix", "-f", "-s1", "-b2", "-e2", local_path], check=True)


def write_sumstat_output(df: pl.DataFrame, output_path: str) -> None:
    """Write full + mlog10p>4 filtered sumstat files, uploading to GCS if needed.
    Expects df already prepared with final column names/order."""
    is_gcs = output_path.startswith("gs://")

    if is_gcs:
        tmpdir = tempfile.mkdtemp()
        local_full = f"{tmpdir}/full.tsv.gz"
        local_filt = f"{tmpdir}/filtered.tsv.gz"
    else:
        local_full = output_path
        local_filt = output_path.replace(".tsv.gz", ".mlog10p_gt4.tsv.gz")

    write_bgzip(df, local_full)
    print(f"  wrote {df.height} variants")

    filtered = df.filter(pl.col("mlog10p") > 4)
    write_bgzip(filtered, local_filt)
    print(f"  wrote {filtered.height} variants with mlog10p > 4")

    if is_gcs:
        upload_to_gcs(local_full, output_path)
        filtered_gcs = output_path.replace(".tsv.gz", ".mlog10p_gt4.tsv.gz")
        upload_to_gcs(local_filt, filtered_gcs)
    else:
        print(f"  output: {output_path}")
        print(f"  filtered: {local_filt}")


def write_exome_output(
    df: pl.DataFrame,
    output_path: str,
    tabix_args: list[str],
    mlog10p_col: str = "mlog10p",
    mlog10p_threshold: float = 4,
) -> None:
    """Write full + mlog10p-filtered exome result files, uploading to GCS if needed.

    Args:
        df: DataFrame with final columns
        output_path: local or gs:// path for full results (.tsv.gz)
        tabix_args: tabix -s/-b/-e args, e.g. ["-s2", "-b3", "-e3"] for variants
        mlog10p_col: column name to filter on
        mlog10p_threshold: threshold for filtered file
    """
    is_gcs = output_path.startswith("gs://")
    filtered_suffix = f".mlog10p_gt{int(mlog10p_threshold)}.tsv.gz"

    if is_gcs:
        tmpdir = tempfile.mkdtemp()
        local_full = f"{tmpdir}/full.tsv.gz"
        local_filt = f"{tmpdir}/filtered.tsv.gz"
    else:
        local_full = output_path
        local_filt = output_path.replace(".tsv.gz", filtered_suffix)

    _write_bgzip_tabix(df, local_full, tabix_args)
    print(f"  wrote {df.height} rows to {output_path}")

    filtered = df.filter(pl.col(mlog10p_col) > mlog10p_threshold)
    _write_bgzip_tabix(filtered, local_filt, tabix_args)
    filtered_gcs = output_path.replace(".tsv.gz", filtered_suffix)
    print(f"  wrote {filtered.height} rows with {mlog10p_col} > {mlog10p_threshold} to {filtered_gcs}")

    if is_gcs:
        upload_to_gcs(local_full, output_path)
        upload_to_gcs(local_filt, filtered_gcs)


def _write_bgzip_tabix(df: pl.DataFrame, local_path: str, tabix_args: list[str]) -> None:
    """Write DataFrame as bgzipped TSV with tabix index using given tabix args."""
    with subprocess.Popen(
        ["bgzip", "-c"],
        stdin=subprocess.PIPE,
        stdout=open(local_path, "wb"),
    ) as proc:
        df.write_csv(proc.stdin, separator="\t", null_value="NA")
    subprocess.run(["tabix", "-f"] + tabix_args + [local_path], check=True)


def read_gnomad_filtered(filepath: str, columns: list[str] | None = None) -> pl.DataFrame:
    """Read a previously saved filtered gnomAD TSV (plain or gzipped)."""
    schema_overrides = {
        "#chr": pl.Int32, "pos": pl.Int32,
        **{c: pl.Float64 for c in GNOMAD_AF_COLS},
    }
    if filepath.endswith(".gz"):
        proc = subprocess.Popen(["zcat", filepath], stdout=subprocess.PIPE)
        df = pl.read_csv(
            proc.stdout, separator="\t", null_values=["NA"],
            columns=columns, schema_overrides=schema_overrides,
            ignore_errors=True,
        )
        proc.wait()
    else:
        df = pl.read_csv(
            filepath, separator="\t", null_values=["NA"],
            columns=columns, schema_overrides=schema_overrides,
            ignore_errors=True,
        )
    return df.rename({"#chr": "chr"})


def build_rsid_set(df: pl.DataFrame) -> set[str]:
    """Extract set of rsids for fast gnomAD filtering."""
    return set(df["rsid"].drop_nulls().to_list())


def stream_gnomad_by_rsid(
    gnomad_path: str,
    rsids: set[str],
    save_path: str,
) -> pl.DataFrame:
    """Stream gnomAD from GCS, keeping only rows with matching rsids.
    Saves filtered rows to save_path as a TSV for reuse."""
    proc = subprocess.Popen(
        f"gsutil cat {gnomad_path} | zcat",
        shell=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        bufsize=4_000_000,
    )

    header_line = proc.stdout.readline().strip()
    header = header_line.split("\t")
    col_idx = {}
    for col in GNOMAD_KEEP_COLS:
        if col in header:
            col_idx[col] = header.index(col)

    rsids_idx = col_idx["rsids"]

    out_header = "\t".join(col_idx.keys())
    outf = open(save_path, "w")
    outf.write(out_header + "\n")

    n_lines = 0
    n_matches = 0
    for line in proc.stdout:
        n_lines += 1
        if n_lines % 10_000_000 == 0:
            print(f"  gnomAD: {n_lines / 1e6:.0f}M lines scanned, {n_matches} matches", flush=True)

        fields = line.split("\t")
        if len(fields) <= rsids_idx:
            continue

        rsid_field = fields[rsids_idx]
        if rsid_field and rsids.intersection(rsid_field.split(",")):
            out_fields = [fields[idx].rstrip("\n") for idx in col_idx.values()]
            outf.write("\t".join(out_fields) + "\n")
            n_matches += 1

    outf.close()
    proc.wait()
    print(f"  gnomAD: done. {n_lines / 1e6:.1f}M lines scanned, {n_matches} matches", flush=True)
    print(f"  saved filtered gnomAD to {save_path}")

    return read_gnomad_filtered(save_path)
