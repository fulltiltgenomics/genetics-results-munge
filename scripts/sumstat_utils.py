"""Shared utilities for sumstat munging scripts, and the IO helpers every family shares.

`upload_to_gcs` and `fetch` are not sumstat-specific: peak_utils re-exports the upload, and
`fetch` is shared by the coordinate-free products, whose inputs are public downloads. They
live here because a second module for two functions is how the copies start.
"""

import argparse
import hashlib
import shutil
import subprocess
import sys
import tempfile
import urllib.error
import urllib.request
from pathlib import Path

import numpy as np
import polars as pl


# the served gnomAD sites file, as build_gnomad_annotation.py writes it: one row per variant,
# sex chromosomes coded 23/24, consequence from the canonical transcript
GNOMAD_FILE = "gnomad.genomes.exomes.v4.1.1.sites.vep115.canonical.tsv.bgz"
GNOMAD_DEFAULT = f"gs://daly-genetics-results/gnomad/{GNOMAD_FILE}"
GNOMAD_LOCAL_DEFAULT = f"/mnt/disks/data/gnomad/4.1.1/out/{GNOMAD_FILE}"

GNOMAD_SOURCES = ["any", "g", "e"]
GNOMAD_SOURCE_HELP = (
    "gnomAD rows to use: any (default) takes whichever row a variant has, genomes first where the "
    "file carries both; g / e restrict to genomes / exomes"
)

GNOMAD_AF_COLS = ["AF", "AF_afr", "AF_amr", "AF_asj", "AF_eas", "AF_fin", "AF_mid", "AF_nfe", "AF_remaining", "AF_sas"]

GNOMAD_KEEP_COLS = ["#chr", "pos", "ref", "alt", "rsids", "genome_or_exome"] + GNOMAD_AF_COLS

# the suite's one numeric coding of the non-autosomes, as the `chr` INT64 columns in BigQuery and
# the tabix seqnames carry it. 25 is not a code.
CHROM_CODES = {"X": "23", "Y": "24", "M": "26", "MT": "26"}

# code -> name, taking each code's first spelling above, so the mitochondrion is M: the reverse
# direction is used to address `chr`-prefixed references (gnomAD, GENCODE, UCSC), which say chrM
CHROM_NAMES: dict[str, str] = {}
for _name, _code in CHROM_CODES.items():
    CHROM_NAMES.setdefault(_code, _name)

# CHROM_CODES as awk statements filling an array `code`, for the scripts whose row loop is awk
AWK_CHROM_CODES = " ".join(f'code["{_name}"] = {_code};' for _name, _code in CHROM_CODES.items())


def chrom_code(name: str) -> str:
    """`chrom_code_expr` for one value, for the scripts that read a file line by line."""
    base = name[3:] if name[:3].lower() == "chr" else name
    return CHROM_CODES.get(base.upper(), base)


def chrom_code_expr(col: str | pl.Expr | pl.Series = "chr") -> pl.Expr:
    """Chromosome token from a seqname: strip a leading `chr`, then X=23, Y=24, M/MT=26
    (`CHROM_CODES`), anything else as it was.

    The result stays a string, so the caller chooses the cast and what a contig that will not
    cast means: a strict cast makes a scaffold an error, a non-strict one a null to drop or count.
    """
    expr = pl.col(col) if isinstance(col, str) else col
    # case-insensitive like CHR_STRING_TO_INT_SQL in genetics-results-db, so a seqname coded here
    # is one the BigQuery load also accepts. a case-sensitive match drops an unbounded class
    # instead: ANY case-variant `chr` prefix on an otherwise canonical chromosome casts to null or
    # fails the canonical filter as if it were a scaffold. wrong only if a source used case to
    # distinguish two different contigs, which no reference assembly does.
    base = expr.cast(pl.Utf8).str.replace("(?i)^chr", "")
    upper = base.str.to_uppercase()
    return pl.when(upper.is_in(list(CHROM_CODES))).then(upper.replace(CHROM_CODES)).otherwise(base)


def upload_to_gcs(local_path: str, gcs_path: str) -> None:
    """Copy a local file (and its .tbi if present) to GCS."""
    subprocess.run(["gcloud", "storage", "cp", local_path, gcs_path], check=True)
    print(f"  uploaded {gcs_path}")
    tbi_local = local_path + ".tbi"
    if Path(tbi_local).exists():
        subprocess.run(["gcloud", "storage", "cp", tbi_local, gcs_path + ".tbi"], check=True)
        print(f"  uploaded {gcs_path}.tbi")


def fetch(url: str, dest: Path, *, download: bool = True, timeout: int = 120) -> Path:
    """Download `url` to `dest` unless it is already cached there.

    Written to a `.part` sibling and renamed, so an interrupted run leaves no half file that
    the next one would treat as cached. `download=False` refuses to reach the network at all,
    for a caller whose own --download flag gates the whole run rather than each call site.
    """
    if dest.exists():
        print(f"  cached: {dest}", file=sys.stderr)
        return dest
    if not download:
        raise SystemExit(f"{dest} is not cached and --download was not given")
    dest.parent.mkdir(parents=True, exist_ok=True)
    print(f"  downloading {url}", file=sys.stderr)
    tmp = dest.with_suffix(dest.suffix + ".part")
    request = urllib.request.Request(url, headers={"User-Agent": "genetics-results-munge"})
    try:
        with urllib.request.urlopen(request, timeout=timeout) as response, tmp.open("wb") as out:
            shutil.copyfileobj(response, out)
    except (urllib.error.URLError, TimeoutError) as err:
        tmp.unlink(missing_ok=True)
        raise SystemExit(
            f"could not fetch {url}: {err}\n"
            f"download it by hand and put it at {dest}, then rerun without --download"
        ) from err
    tmp.rename(dest)
    return dest


def fetch_gs(uri: str, cache_dir: Path, *, billing_project: str | None = None) -> Path:
    """Local copy of a `gs://` object, downloaded once; a local path is passed through.

    `fetch` above is the cache for http sources and cannot speak `gs://`. Written to a `.part`
    sibling and renamed so an interrupted copy is not mistaken for a cached file.
    """
    if not uri.startswith("gs://"):
        return Path(uri)
    cache_dir.mkdir(parents=True, exist_ok=True)
    local = cache_dir / uri.rsplit("/", 1)[-1]
    if local.exists():
        print(f"  cached: {local}", file=sys.stderr)
        return local
    part = local.with_name(local.name + ".part")
    command = ["gcloud", "storage", "cp", uri, str(part)]
    if billing_project:
        command.append(f"--billing-project={billing_project}")
    subprocess.run(command, check=True)
    part.rename(local)
    return local


def sha256(path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def tabix_header(path: str) -> list[str]:
    """Column names from the header line of a tabix-indexed file."""
    return (
        subprocess.run(["tabix", "-H", path], capture_output=True, text=True, check=True)
        .stdout.strip()
        .lstrip("#")
        .split("\t")
    )


def sci_expr(expr: pl.Expr) -> pl.Expr:
    """Format a float expression as %.3e in one numpy call rather than per row, keeping nulls null."""
    formatted = expr.fill_null(0.0).map_batches(
        lambda s: pl.Series(np.char.mod("%.3e", s.to_numpy())), return_dtype=pl.String
    )
    return pl.when(expr.is_null()).then(None).otherwise(formatted)


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
    mlog10p_col: str | None = "mlog10p",
    mlog10p_threshold: float = 4,
    per_trait_dir: str | None = None,
    trait_col: str = "trait_original",
) -> None:
    """Write full + mlog10p-filtered exome result files, uploading to GCS if needed.

    Args:
        df: DataFrame with final columns
        output_path: local or gs:// path for full results (.tsv.gz)
        tabix_args: tabix -s/-b/-e args, e.g. ["-s2", "-b3", "-e3"] for variants
        mlog10p_col: column name to filter on; None writes no filtered companion, for a
            product that carries no p-value at all (ASC counts) rather than one whose
            filter would be empty
        mlog10p_threshold: threshold for filtered file
        per_trait_dir: if given, also write one unfiltered `<trait>.tsv.gz` per
            distinct trait there, so the API can serve the dataset by trait the
            same way it serves genebass (see split_burden_per_trait.py, which
            does this for files too large to hold in a DataFrame)
        trait_col: column holding the trait key that names the per-trait files
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

    if mlog10p_col is not None:
        filtered = df.filter(pl.col(mlog10p_col) > mlog10p_threshold)
        _write_bgzip_tabix(filtered, local_filt, tabix_args)
        filtered_gcs = output_path.replace(".tsv.gz", filtered_suffix)
        print(f"  wrote {filtered.height} rows with {mlog10p_col} > {mlog10p_threshold} to {filtered_gcs}")

    if is_gcs:
        upload_to_gcs(local_full, output_path)
        if mlog10p_col is not None:
            upload_to_gcs(local_filt, filtered_gcs)

    if per_trait_dir:
        _write_per_trait(df, per_trait_dir, tabix_args, trait_col)


def _write_per_trait(
    df: pl.DataFrame, per_trait_dir: str, tabix_args: list[str], trait_col: str
) -> None:
    """Write one unfiltered `<trait>.tsv.gz` (+ .tbi) per distinct trait."""
    is_gcs = per_trait_dir.startswith("gs://")
    local_dir = Path(tempfile.mkdtemp()) if is_gcs else Path(per_trait_dir)
    local_dir.mkdir(parents=True, exist_ok=True)

    for (trait,), part in df.partition_by(trait_col, as_dict=True).items():
        # the trait becomes a file name and, in the API, a path component
        if "/" in trait or trait in (".", ".."):
            raise ValueError(f"trait {trait!r} is not safe to use as a file name")
        local_path = str(local_dir / f"{trait}.tsv.gz")
        _write_bgzip_tabix(part, local_path, tabix_args)
        print(f"  wrote {part.height} rows for trait {trait}")
        if is_gcs:
            upload_to_gcs(local_path, f"{per_trait_dir.rstrip('/')}/{trait}.tsv.gz")


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


def save_gnomad_by_variant_key(
    gnomad_path: str,
    variant_keys: set[str],
    save_path: str,
    keep_cols: list[str] = GNOMAD_KEEP_COLS,
) -> None:
    """Stream a local gnomAD file, saving the `keep_cols` of rows whose chr:pos:ref:alt is in
    variant_keys to save_path as a TSV for reuse."""
    proc = subprocess.Popen(
        f"zcat '{gnomad_path}'",
        shell=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        bufsize=4_000_000,
    )

    header_line = proc.stdout.readline().strip()
    header = header_line.split("\t")
    col_idx = {}
    for col in keep_cols:
        if col in header:
            col_idx[col] = header.index(col)

    chr_idx = col_idx["#chr"]
    pos_idx = col_idx["pos"]
    ref_idx = col_idx["ref"]
    alt_idx = col_idx["alt"]

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
        if len(fields) <= alt_idx:
            continue

        key = f"{fields[chr_idx]}:{fields[pos_idx]}:{fields[ref_idx]}:{fields[alt_idx]}"
        if key in variant_keys:
            out_fields = [fields[idx].rstrip("\n") for idx in col_idx.values()]
            outf.write("\t".join(out_fields) + "\n")
            n_matches += 1

    outf.close()
    proc.wait()
    print(f"  gnomAD: done. {n_lines / 1e6:.1f}M lines scanned, {n_matches} matches", flush=True)
    print(f"  saved filtered gnomAD to {save_path}")


def gnomad_plot_parser(
    description: str | None,
    input_help: str,
    *,
    gnomad_default: str = GNOMAD_DEFAULT,
    gnomad_help: str = "GCS path to gnomAD bgz file",
    af_plot_flag: bool = True,
) -> argparse.ArgumentParser:
    """The flags of a munge that aligns one sumstat file against gnomAD and plots AF against AF.

    `af_plot_flag=False` is for a munge that cannot write its output without the gnomAD join, so
    has nothing to make optional.
    """
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("--input", required=True, help=input_help)
    parser.add_argument("--gnomad", default=gnomad_default, help=gnomad_help)
    parser.add_argument("--af-col", default="AF", choices=GNOMAD_AF_COLS, help="gnomAD AF column for plot (default: AF)")
    parser.add_argument("--output", help="Output munged TSV path (default: derived from input)")
    if af_plot_flag:
        parser.add_argument("--gnomad-af-plot", action="store_true", help="Stream gnomAD and create AF-AF plot")
    parser.add_argument("--gnomad-filtered", help="Path to previously saved filtered gnomAD TSV (skip streaming)")
    parser.add_argument("--gnomad-source", default="any", choices=GNOMAD_SOURCES, help=GNOMAD_SOURCE_HELP)
    parser.add_argument("--plot", help="Output AF-AF plot PNG path (default: derived from input)")
    return parser


def pick_gnomad_rows(
    gnomad: pl.DataFrame, source: str = "any", key: tuple[str, ...] = ("chr", "pos", "ref", "alt"),
) -> pl.DataFrame:
    """The gnomAD rows a munge joins against, from rows that carry `genome_or_exome`.

    `any` keeps one row per `key`: the only row where the file has one per variant, the genomes
    row where it carries a genomes and an exomes row. `g` / `e` keep that source's rows only. On
    a file already cut to one row per variant that drops every variant whose surviving row came
    from the other source, with nothing to show for it but a lower match count, which is why it
    is a restriction to ask for and not the default.
    """
    if source == "any":
        # "g" sorts after "e", so descending puts the genomes row first
        out = gnomad.sort("genome_or_exome", descending=True, maintain_order=True).unique(
            subset=list(key), keep="first", maintain_order=True,
        )
        print(f"  gnomAD, one row per variant (genomes before exomes): {out.height} rows (of {gnomad.height})")
    else:
        label = "genomes" if source == "g" else "exomes"
        out = gnomad.filter(pl.col("genome_or_exome") == source)
        print(f"  gnomAD {label} only: {out.height} rows (of {gnomad.height})")
    return out


def _pyplot():
    # imported on use: matplotlib is not installed everywhere this module is (the cellxgene image)
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    return plt


def _af_af_axis(ax, df: pl.DataFrame, gnomad_af_col: str, xlabel: str, ylabel: str, title: str) -> None:
    sample = df.sample(n=min(200_000, df.height), seed=42)
    ax.scatter(
        sample["af"].to_numpy(),
        sample[gnomad_af_col].to_numpy(),
        alpha=0.3, s=1, c="blue",
        rasterized=True,
    )
    ax.plot([0, 1], [0, 1], "k--", alpha=0.3, linewidth=0.5)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_aspect("equal")


def create_af_af_plot(df: pl.DataFrame, af_col: str, output_path: str, input_name: str, xlabel: str) -> None:
    """AF-AF scatter: the study's alt AF (`af`) against gnomAD's `af_col`."""
    plt = _pyplot()
    fig, ax = plt.subplots(figsize=(8, 8))
    _af_af_axis(ax, df, af_col, xlabel, f"gnomAD {af_col}", f"AF-AF: {input_name} (n={df.height:,})")
    fig.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  plot saved to {output_path}")


def create_af_af_plot_by_type(
    df: pl.DataFrame,
    af_col: str,
    output_path: str,
    input_name: str,
    xlabel: str,
    *,
    gnomad_af_col: str | None = None,
    indel_label: str = "Indels",
) -> None:
    """AF-AF scatter with separate panels for SNPs and indels, split on the boolean `is_snp` column.

    `gnomad_af_col` names the column holding the gnomAD frequency when the join renamed it away
    from `af_col`, which still labels the axis.
    """
    plt = _pyplot()
    fig, axes = plt.subplots(1, 2, figsize=(16, 8))

    for ax, subset, label in [
        (axes[0], df.filter(pl.col("is_snp")), "SNPs"),
        (axes[1], df.filter(~pl.col("is_snp")), indel_label),
    ]:
        if subset.height == 0:
            ax.set_title(f"{label}: no data")
            continue
        _af_af_axis(ax, subset, gnomad_af_col or af_col, xlabel, f"gnomAD {af_col}", f"{label} (n={subset.height:,})")

    fig.suptitle(f"AF-AF: {input_name}", fontsize=14)
    fig.tight_layout()
    fig.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  plot saved to {output_path}")


def join_gnomad_by_rsid(
    sumstats: pl.DataFrame, gnomad: pl.DataFrame, freq_cols: list[str], gnomad_source: str = "any",
) -> pl.DataFrame:
    """Join daner-style sumstats (A1 effect allele, A2 other) with gnomAD on rsid, classify allele
    orientation, and flip `beta`, `af` and every column of `freq_cols` to the alt allele."""
    gnomad_dedup = pick_gnomad_rows(gnomad, gnomad_source)

    # gnomAD rsids are comma-separated: one row per rsid
    gnomad_exploded = gnomad_dedup.with_columns(
        pl.col("rsids").str.split(",").alias("rsid_list"),
    ).explode("rsid_list").rename({"rsid_list": "rsid"})

    joined = sumstats.join(
        gnomad_exploded,
        on="rsid",
        how="inner",
        suffix="_gnomad",
    )

    joined = joined.with_columns(
        pl.when(
            (pl.col("A1").str.to_uppercase() == pl.col("ref").str.to_uppercase())
            & (pl.col("A2").str.to_uppercase() == pl.col("alt").str.to_uppercase())
        ).then(pl.lit("A1_is_ref"))
        .when(
            (pl.col("A1").str.to_uppercase() == pl.col("alt").str.to_uppercase())
            & (pl.col("A2").str.to_uppercase() == pl.col("ref").str.to_uppercase())
        ).then(pl.lit("A1_is_alt"))
        .otherwise(pl.lit("mismatch"))
        .alias("allele_match")
    )

    n_a1_ref = joined.filter(pl.col("allele_match") == "A1_is_ref").height
    n_a1_alt = joined.filter(pl.col("allele_match") == "A1_is_alt").height
    n_mismatch = joined.filter(pl.col("allele_match") == "mismatch").height
    matched_rsids = joined["rsid"].unique().len()
    n_unmatched = sumstats["rsid"].unique().len() - matched_rsids
    print(f"  A1=ref (flip): {n_a1_ref}, A1=alt (no flip): {n_a1_alt}, mismatch: {n_mismatch}, no gnomAD match: {n_unmatched}")

    result = joined.filter(pl.col("allele_match") != "mismatch")

    # use gnomAD chr/pos (build 38) for all rows
    result = result.with_columns(
        pl.col("chr_gnomad").alias("chr"),
        pl.col("pos_gnomad").alias("pos"),
    )

    # A1_is_ref: A1 is ref, A2 is alt — flip beta and AFs so effect is on alt
    # A1_is_alt: A1 is alt, A2 is ref — beta/AFs already represent alt effect
    is_flip = pl.col("allele_match") == "A1_is_ref"
    return result.with_columns(
        pl.when(is_flip).then(-pl.col("beta")).otherwise(pl.col("beta")).alias("beta"),
        *[
            pl.when(is_flip).then((1.0 - pl.col(c)).round(6)).otherwise(pl.col(c)).alias(c)
            for c in ["af", *freq_cols]
        ],
    )


def prepare_rsid_output(df: pl.DataFrame, freq_cols: list[str]) -> pl.DataFrame:
    """Format beta/SE and select the output columns, `freq_cols` between `af` and `INFO`."""
    return df.with_columns(
        pl.col("beta").map_elements(lambda x: f"{x:.3e}" if x is not None else None, return_dtype=pl.Utf8).alias("beta"),
        pl.col("SE").map_elements(lambda x: f"{x:.3e}" if x is not None else None, return_dtype=pl.Utf8).alias("se"),
    ).select(
        "chr", "pos", "ref", "alt", "rsid",
        "beta", "se", "mlog10p",
        "af", *freq_cols, "INFO",
        "neff", "n_cases", "n_controls",
    ).rename({"chr": "#chr"}).sort(
        "#chr", "pos", "ref", "alt",
    )


def run_rsid_munge(description: str | None, input_help: str, read_sumstats, freq_cols: list[str], plot_xlabel: str) -> None:
    """The whole command line of a daner-style sumstat munge: parse the `gnomad_plot_parser` flags,
    read --input with `read_sumstats`, optionally align to gnomAD by rsid and plot, write the output."""
    args = gnomad_plot_parser(description, input_help).parse_args()
    input_stem = Path(args.input).name.replace(".gz", "")
    output_dir = Path(args.input).parent
    output_path = args.output or str(output_dir / f"{input_stem}.munged.tsv.gz")
    plot_path = args.plot or str(output_dir / f"{input_stem}.af_af.png")
    gnomad_save_path = str(output_dir / f"{input_stem}.gnomad_filtered.tsv")

    print(f"Reading {args.input}...")
    sumstats = read_sumstats(args.input)
    print(f"  {sumstats.height} variants loaded")

    if args.gnomad_af_plot:
        if args.gnomad_filtered:
            print(f"Reading filtered gnomAD from {args.gnomad_filtered}...")
            gnomad = read_gnomad_filtered(args.gnomad_filtered)
        else:
            print("Building rsid set...")
            rsids = build_rsid_set(sumstats)
            print(f"  {len(rsids)} unique rsids")

            print(f"Streaming gnomAD from {args.gnomad}...")
            gnomad = stream_gnomad_by_rsid(args.gnomad, rsids, gnomad_save_path)

        print(f"  {gnomad.height} gnomAD rows loaded")

        print("Joining with gnomAD...")
        sumstats = join_gnomad_by_rsid(sumstats, gnomad, freq_cols, args.gnomad_source)

        print("Creating AF-AF plot...")
        create_af_af_plot(sumstats, args.af_col, plot_path, input_stem, plot_xlabel)

    print("Writing output...")
    write_sumstat_output(prepare_rsid_output(sumstats, freq_cols), output_path)

    print("Done.")
