#!/usr/bin/env python3
"""Munge the AIH autoimmune meta-analysis summary statistics (aih, aitt1, aitt2).

Source files
------------
`/mnt/disks/data/{aih,aitt1,aitt2}_meta.txt.gz`, one per phenotype, tab separated with header

    SNP CHR BP A1 A0 A1FREQ BETA SE P N Direction HetISq HetChiSq HetDf HetPVal

Format assumptions (checked against the delivered files)
--------------------------------------------------------
- positions are already GRCh38 and chrX is already coded as 23 (`X` is still mapped
  defensively);
- `SNP` is always a dbSNP rs id;
- `A1` is the effect allele — `A1FREQ` and `BETA` refer to it — and `A0` is the other
  allele. Neither is guaranteed to be the reference allele, so the orientation is
  resolved against gnomAD;
- indels are left-aligned, VCF style (`CCT`/`C`), i.e. the same representation gnomAD
  uses, so they can be matched by exact allele string;
- `P` never underflows to 0 in these files, but the beta/se fallback is applied anyway;
- `N` is the per-variant total sample size; the files carry no case/control split;
- ~24k variants per file appear twice, under two different rs ids and meta-analysing
  disjoint cohort subsets (`Direction` `??-+` vs `++??`). Only the row with the larger
  `N` is kept — see `dedup_variants`.

Allele orientation
------------------
Each variant is matched to gnomAD on `(chr, pos)` plus an exact match of the allele
*pair* in either orientation. `ref`/`alt` are then taken from gnomAD, and `beta` and
`af` are flipped whenever `A1` turned out to be the reference allele. Indels are
flipped just like SNPs — unlike `munge_gp2.py`, which cannot, because here both alleles
must match gnomAD's strings exactly for the variant to be kept at all. No strand
flipping is attempted, so a variant reported on the opposite strand simply fails to
match and is dropped rather than being silently mis-oriented.

About 1% of variants — all of them indels — match gnomAD in *both* orientations and are
dropped. gnomAD legitimately carries a deletion (`TA` → `T`) and an insertion
(`T` → `TA`) at the same position, so a study row reporting alleles `T` and `TA` there
fits either record, once as the alt allele and once as the ref, which is the difference
between `beta` and `-beta`. The two records are indistinguishable on the evidence the
file provides: the reference genome satisfies both (`ref` is a prefix either way), and
dbSNP gives both the same rs id, so the study's `SNP` column does not break the tie
either (it matched both records for 1620 of the 1743 ambiguous variants on chr21).
Guessing from allele frequency would silently mis-sign whichever ones it got wrong, so
they are dropped and counted instead.

Because the three files share almost all of their variants, gnomAD is streamed once for
the union of their positions (`--gnomad-filter-only`) and the saved filtered file is
then reused for every phenotype with `--gnomad-filtered`.

Flags
-----
--input               one or more `*_meta.txt.gz` files (all are munged in one run)
--output-dir          directory for `<stem>.munged.tsv.gz`, may be a gs:// path
--output              explicit output path, only allowed with a single --input
--gnomad              gnomAD sites bgz to stream (default: the local copy)
--gnomad-filtered     previously saved filtered gnomAD TSV, skips the streaming pass
--gnomad-filter-only  only stream gnomAD for the union of the inputs' positions and exit
--jobs                chromosomes scanned in parallel by --gnomad-filter-only (default 2)
--gnomad-source       g / e / any (default any: prefer genomes, fall back to exomes)
--af-col              gnomAD AF column used on the AF-AF plot (default AF)
--gnomad-af-plot      also draw the AF-AF plot against gnomAD
"""

import argparse
import shutil
import subprocess
import tempfile
import warnings
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import polars as pl
from scipy.special import log_ndtr

from sumstat_utils import (
    GNOMAD_AF_COLS, read_gnomad_filtered, write_sumstat_output,
)


GNOMAD_PATH = "/mnt/disks/data/gnomad/gnomad.genomes.exomes.v4.0.sites.v2.tsv.bgz"

# the subset of gnomAD columns the filtered file keeps, in output order
GNOMAD_FILTER_COLS = ["#chr", "pos", "ref", "alt", "rsids", "genome_or_exome"] + GNOMAD_AF_COLS

OUTPUT_COLS = [
    "#chr", "pos", "ref", "alt", "rsid",
    "beta", "se", "mlog10p", "af", "n",
    "Direction", "HetISq", "HetChiSq", "HetDf", "HetPVal",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--input", required=True, nargs="+", help="One or more *_meta.txt.gz files")
    parser.add_argument("--output-dir", help="Output directory (local or gs://, default: directory of the first input)")
    parser.add_argument("--output", help="Explicit output path, only valid with a single --input")
    parser.add_argument("--gnomad", default=GNOMAD_PATH, help=f"gnomAD sites bgz (default: {GNOMAD_PATH})")
    parser.add_argument("--gnomad-filtered", help="Path to previously saved filtered gnomAD TSV (skip streaming)")
    parser.add_argument("--gnomad-filter-only", action="store_true", help="Stream gnomAD for the union of the inputs' positions, save it and exit")
    parser.add_argument("--jobs", type=int, default=2, help="Chromosomes scanned in parallel by --gnomad-filter-only (default: 2)")
    parser.add_argument("--gnomad-source", default="any", choices=["g", "e", "any"], help="gnomAD genomes (g), exomes (e), or prefer genomes and fall back to exomes (any, default)")
    parser.add_argument("--af-col", default="AF", choices=GNOMAD_AF_COLS, help="gnomAD AF column for the plot (default: AF)")
    parser.add_argument("--gnomad-af-plot", action="store_true", help="Draw the AF-AF plot against gnomAD")
    parser.add_argument("--plot-dir", help="Directory for AF-AF plots (default: directory of the first input)")
    return parser.parse_args()


def read_sumstats(filepath: str) -> pl.DataFrame:
    """Read one meta file, map chr, compute mlog10p. Keeps BETA/A1FREQ unflipped."""
    df = pl.read_csv(
        filepath,
        separator="\t",
        null_values=["NA"],
        truncate_ragged_lines=True,
        schema_overrides={
            "SNP": pl.Utf8,
            "CHR": pl.Utf8,
            "BP": pl.Int32,
            "A1": pl.Utf8,
            "A0": pl.Utf8,
            "A1FREQ": pl.Float64,
            "BETA": pl.Float64,
            "SE": pl.Float64,
            "P": pl.Float64,
            "N": pl.Int32,
            "HetISq": pl.Float64,
            "HetChiSq": pl.Float64,
            "HetDf": pl.Int32,
            "HetPVal": pl.Float64,
        },
        ignore_errors=True,
    ).rename({"SNP": "rsid", "BP": "pos", "N": "n"})

    # chrX is already 23 in these files; map it anyway in case a delivery differs
    df = df.with_columns(
        pl.col("CHR").str.replace("X", "23").cast(pl.Int32, strict=False).alias("chr"),
    ).drop("CHR")
    n_bad_chr = df.filter(pl.col("chr").is_null()).height
    if n_bad_chr:
        print(f"  dropped {n_bad_chr} rows with unparseable chromosome")
    df = df.filter(pl.col("chr").is_not_null())

    with warnings.catch_warnings():
        warnings.filterwarnings("ignore")
        df = df.with_columns(
            # mlog10p with underflow handling
            pl.when(pl.col("P").eq(0) | pl.col("P").is_null())
            .then(
                ((-log_ndtr(-(pl.col("BETA") / pl.col("SE")).abs()) - np.log(2)) / np.log(10)).round(4)
            )
            .otherwise((-np.log10(pl.col("P"))).round(4))
            .alias("mlog10p"),
        )

    n_underflow = df.filter(pl.col("P").eq(0) | pl.col("P").is_null()).height
    print(f"  {df.height} variants, {n_underflow} with mlog10p derived from beta/se")
    return df


def write_position_union(inputs: list[str], positions_path: str) -> int:
    """Write the sorted union of (chr, pos) over all inputs, X mapped to 23."""
    cat = " ".join(f"zcat {shell_quote(f)};" for f in inputs)
    cmd = (
        f"({cat}) | mawk 'NR>1 && $1!=\"SNP\" {{ c=($2==\"X\")?23:$2; print c\"\\t\"$3 }}' "
        f"| sort -u -S 2G -k1,1n -k2,2n > {shell_quote(positions_path)}"
    )
    subprocess.run(["bash", "-c", cmd], check=True)
    n = int(subprocess.run(["wc", "-l", positions_path], capture_output=True, text=True, check=True).stdout.split()[0])
    print(f"  {n} unique positions across {len(inputs)} file(s) -> {positions_path}")
    return n


def shell_quote(s: str) -> str:
    return "'" + s.replace("'", "'\\''") + "'"


def stream_gnomad(gnomad_path: str, positions_path: str, save_path: str, jobs: int = 2) -> None:
    """Keep the gnomAD rows whose (chr, pos) is in positions_path.

    gnomAD v4 is ~500 GB of text, and a single decompress-and-filter pipeline runs at
    the speed of one awk process, leaving most of the machine idle. So the scan is done
    one chromosome at a time through `tabix`, `jobs` of them at once. That also shrinks
    each filter's position set to a single chromosome, which matters more than it looks:
    the whole-genome set is ~11M keys and thrashes the cache on every lookup.
    """
    tmpdir = tempfile.mkdtemp(prefix="gnomad_filter_", dir=str(Path(save_path).parent))
    chroms = subprocess.run(
        ["tabix", "-l", gnomad_path], capture_output=True, text=True, check=True,
    ).stdout.split()

    # one position file per chromosome, named by the gnomAD contig name (X, not 23)
    subprocess.run(
        ["bash", "-c",
         f"mawk -v d={shell_quote(tmpdir)} '{{ c = ($1 == 23) ? \"X\" : $1; "
         f"print $2 > (d \"/pos_\" c) }}' {shell_quote(positions_path)}"],
        check=True,
    )
    chroms = [c for c in chroms if Path(tmpdir, f"pos_{c}").exists()]
    print(f"  scanning {gnomad_path} for {len(chroms)} chromosomes, {jobs} at a time...", flush=True)

    # gnomAD column numbers, 1-based, in GNOMAD_FILTER_COLS order
    fields = "c,$2,$3,$4,$5,$21,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17"
    # keys stay strings on both sides: awk would render a numeric subscript through
    # CONVFMT, and positions are plain decimal integers, so the raw text already matches
    awk = (
        'BEGIN { FS=OFS="\\t"; while ((getline p < posfile) > 0) keep[p] = 1; close(posfile) } '
        '$2 in keep { print ' + fields + ' }'
    )
    worker = (
        f"tabix {shell_quote(gnomad_path)} \"$1\" "
        f"| mawk -v posfile={shell_quote(tmpdir)}\"/pos_$1\" "
        f"-v c=\"$(test \"$1\" = X && echo 23 || echo \"$1\")\" {shell_quote(awk)} "
        f"| bgzip -c > {shell_quote(tmpdir)}\"/out_$1.gz\""
    )
    subprocess.run(
        ["bash", "-c",
         f"printf '%s\\n' {' '.join(shell_quote(c) for c in chroms)} "
         f"| xargs -P {jobs} -I{{}} bash -c {shell_quote(worker)} _ {{}}"],
        check=True,
    )

    header = "\t".join(GNOMAD_FILTER_COLS)
    parts = " ".join(shell_quote(f"{tmpdir}/out_{c}.gz") for c in chroms)
    subprocess.run(
        ["bash", "-c",
         f"( printf '%s\\n' {shell_quote(header)}; zcat {parts} ) "
         f"| bgzip -c -@2 > {shell_quote(save_path)}"],
        check=True,
    )
    shutil.rmtree(tmpdir)
    print(f"  saved filtered gnomAD to {save_path}")


def dedup_gnomad(gnomad: pl.DataFrame, gnomad_source: str) -> pl.DataFrame:
    """One row per (chr, pos, ref, alt). With 'any', prefer genomes over exomes."""
    if gnomad_source in ("g", "e"):
        label = "genomes" if gnomad_source == "g" else "exomes"
        out = gnomad.filter(pl.col("genome_or_exome") == gnomad_source)
        print(f"  gnomAD {label}: {out.height} rows (of {gnomad.height})")
    else:
        # "g" sorts after "e", so descending puts the genomes row first
        out = gnomad.sort("genome_or_exome", descending=True).unique(
            subset=["chr", "pos", "ref", "alt"], keep="first", maintain_order=True,
        )
        print(f"  gnomAD genomes preferred, exomes as fallback: {out.height} rows (of {gnomad.height})")
    return out.drop("genome_or_exome")


def join_with_gnomad(df: pl.DataFrame, gnomad: pl.DataFrame, af_col: str) -> pl.DataFrame:
    """Match each variant to gnomAD on (chr, pos) and an exact allele pair in either
    orientation, take ref/alt from gnomAD and flip beta/af when A1 is the ref allele.

    A variant matching in *both* orientations is dropped as ambiguous — see the
    module docstring."""
    df = df.with_row_index("_row").with_columns(
        pl.col("A1").str.to_uppercase().alias("a1_upper"),
        pl.col("A0").str.to_uppercase().alias("a0_upper"),
    )
    gnomad = gnomad.with_columns(
        pl.col("ref").str.to_uppercase().alias("ref_upper"),
        pl.col("alt").str.to_uppercase().alias("alt_upper"),
    )

    keep = ["ref", "alt", af_col]

    # A1 is the gnomAD alt allele — beta and A1FREQ already refer to alt
    matched = df.join(
        gnomad.select(["chr", "pos", "ref_upper", "alt_upper"] + keep),
        left_on=["chr", "pos", "a1_upper", "a0_upper"],
        right_on=["chr", "pos", "alt_upper", "ref_upper"],
        how="inner",
    ).with_columns(pl.lit(False).alias("flipped"))

    # A1 is the gnomAD ref allele — beta and A1FREQ have to be flipped onto alt
    swapped = df.join(
        gnomad.select(["chr", "pos", "ref_upper", "alt_upper"] + keep),
        left_on=["chr", "pos", "a1_upper", "a0_upper"],
        right_on=["chr", "pos", "ref_upper", "alt_upper"],
        how="inner",
    ).with_columns(pl.lit(True).alias("flipped"))

    result = pl.concat([matched, swapped], how="vertical")

    ambiguous = result.filter(pl.col("_row").is_duplicated())
    result = result.filter(~pl.col("_row").is_duplicated()).drop("_row")
    n_unmatched = df.height - result.height - ambiguous["_row"].n_unique()
    n_indel = result.filter(
        (pl.col("ref").str.len_chars() > 1) | (pl.col("alt").str.len_chars() > 1)
    ).height
    print(f"  A1=alt (no flip): {result.filter(~pl.col('flipped')).height}, "
          f"A1=ref (flip): {result.filter(pl.col('flipped')).height}, "
          f"no gnomAD allele match: {n_unmatched}, "
          f"ambiguous in both orientations (dropped): {ambiguous['_row'].n_unique()}")
    print(f"  matched indels: {n_indel}")

    result = result.with_columns(
        pl.when(pl.col("flipped")).then(-pl.col("BETA")).otherwise(pl.col("BETA")).alias("BETA"),
        pl.when(pl.col("flipped")).then((1.0 - pl.col("A1FREQ")).round(6)).otherwise(pl.col("A1FREQ")).alias("af"),
    )
    return result.rename({af_col: "af_gnomad"})


def dedup_variants(df: pl.DataFrame) -> pl.DataFrame:
    """Keep one row per variant, the one built from the most cohorts.

    ~24k variants per file are listed twice under two different rs ids, and `Direction`
    shows the two rows meta-analyse disjoint cohort subsets. They cannot be combined
    after the fact, so the larger `n` — the more informative partial meta-analysis —
    wins; mlog10p and rsid only break ties, to keep the choice deterministic."""
    before = df.height
    out = df.sort(
        ["n", "mlog10p", "rsid"], descending=[True, True, False], nulls_last=True,
    ).unique(subset=["chr", "pos", "ref", "alt"], keep="first", maintain_order=True)
    if out.height < before:
        print(f"  dropped {before - out.height} duplicate variant rows, keeping the largest n")
    return out


def create_af_af_plot(df: pl.DataFrame, output_path: str, input_name: str, af_col: str) -> None:
    """AF-AF scatter of the flipped study alt AF against gnomAD, SNPs and indels apart."""
    df = df.with_columns(
        ((pl.col("ref").str.len_chars() == 1) & (pl.col("alt").str.len_chars() == 1)).alias("is_snp"),
    )
    fig, axes = plt.subplots(1, 2, figsize=(16, 8))

    for ax, subset, label in [
        (axes[0], df.filter(pl.col("is_snp")), "SNPs"),
        (axes[1], df.filter(~pl.col("is_snp")), "indels"),
    ]:
        sample_n = min(200_000, subset.height)
        if sample_n == 0:
            ax.set_title(f"{label}: no data")
            continue
        sample = subset.sample(n=sample_n, seed=42)
        ax.scatter(
            sample["af"].to_numpy(), sample["af_gnomad"].to_numpy(),
            alpha=0.3, s=1, c="blue", rasterized=True,
        )
        ax.plot([0, 1], [0, 1], "k--", alpha=0.3, linewidth=0.5)
        ax.set_xlabel("study alt AF (after flip)")
        ax.set_ylabel(f"gnomAD {af_col}")
        ax.set_title(f"{label} (n={subset.height:,})")
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        ax.set_aspect("equal")

    fig.suptitle(f"AF-AF: {input_name}", fontsize=14)
    fig.tight_layout()
    fig.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  plot saved to {output_path}")


def prepare_output(df: pl.DataFrame) -> pl.DataFrame:
    return df.with_columns(
        pl.col("BETA").map_elements(lambda x: f"{x:.3e}" if x is not None else None, return_dtype=pl.Utf8).alias("beta"),
        pl.col("SE").map_elements(lambda x: f"{x:.3e}" if x is not None else None, return_dtype=pl.Utf8).alias("se"),
    ).select(
        pl.col("chr").alias("#chr"),
        *[pl.col(c) for c in OUTPUT_COLS if c != "#chr"],
    ).sort("#chr", "pos", "ref", "alt")


def main():
    args = parse_args()

    if args.output and len(args.input) > 1:
        raise SystemExit("--output can only be used with a single --input; use --output-dir")

    first_dir = Path(args.input[0]).parent
    output_dir = (args.output_dir or str(first_dir)).rstrip("/")
    plot_dir = (args.plot_dir or str(first_dir)).rstrip("/")

    if args.gnomad_filter_only:
        positions_path = str(first_dir / "ai_meta.positions.tsv")
        save_path = str(first_dir / "ai_meta.gnomad_filtered.tsv.gz")
        print("Building the union position set...")
        write_position_union(args.input, positions_path)
        print("Streaming gnomAD...")
        stream_gnomad(args.gnomad, positions_path, save_path, jobs=args.jobs)
        print("Done.")
        return

    if not args.gnomad_filtered:
        raise SystemExit(
            "--gnomad-filtered is required; run once with --gnomad-filter-only over all "
            "inputs first, then pass the saved file here"
        )

    print(f"Reading filtered gnomAD from {args.gnomad_filtered}...")
    gnomad = read_gnomad_filtered(
        args.gnomad_filtered,
        columns=["#chr", "pos", "ref", "alt", "genome_or_exome", args.af_col],
    )
    print(f"  {gnomad.height} gnomAD rows loaded")
    gnomad = dedup_gnomad(gnomad, args.gnomad_source)

    for filepath in args.input:
        stem = Path(filepath).name.replace(".txt.gz", "").replace(".gz", "")
        output_path = args.output or f"{output_dir}/{stem}.munged.tsv.gz"
        print(f"\n=== {stem} ===")
        print(f"Reading {filepath}...")
        df = read_sumstats(filepath)

        print("Joining with gnomAD...")
        df = join_with_gnomad(df, gnomad, args.af_col)
        df = dedup_variants(df)

        if args.gnomad_af_plot:
            print("Creating AF-AF plot...")
            create_af_af_plot(df, f"{plot_dir}/{stem}.af_af.png", stem, args.af_col)

        print("Writing output...")
        write_sumstat_output(prepare_output(df), output_path)
        del df

    print("\nDone.")


if __name__ == "__main__":
    main()
