#!/usr/bin/env python3
"""Munge the Siraj et al. 2026 MPRA per-variant functional annotation from its WIDE per-variant
source into the canonical LONG, tabix-ready `mpra` TSV (one row per variant x cell_line).

Source (input):
  A bgzip'd WIDE TSV, one row per variant. Local default: ~/mpra_annotation.tsv.gz. The canonical
  upstream provenance is gs://finngen-multiome-xavier/batch1_5/misc/MPRA_Siraj_etal/ (built by the
  upstream munge that joins core_mpra.hg38.tsv.bgz to mpra_meta.txt.gz). The input path is a
  positional argument so this same script can later read straight from GCS.

  Layout (see ~/mpra_annotation_columns.md for the full dictionary):
    identity (7): #chrom, pos, variant_id, variant_hg38, ref, alt, cohort
    meta block (10): emVar_meta, active_meta, log2Skew_meta, log2Skew_meta_SE, log2Skew_meta_nlogp,
                     log2Skew_meta_padj, log2FC_meta, log2FC_meta_SE, log2FC_meta_nlogp,
                     log2FC_meta_padj
    5 cell lines x 8 dotted fields: <LINE>.{emVar,active,log2Skew,Skew_logPadj,log2FC,logPadj_BF,
                     mean_RNA_ref,mean_RNA_alt} for LINE in K562,HEPG2,SKNSH,HCT116,A549
  Missing = the literal "NA".

Output (LONG): one row per variant x cell_line where cell_line in {meta,K562,HEPG2,SKNSH,HCT116,A549}.
  The meta and per-line blocks carry mismatched fields, harmonized into one common column set (see
  OUT_COLUMNS). A per-line row is skipped when the variant was not tested in that line (all 8 source
  fields NA); the meta row is emitted whenever any meta field is present.

Output conventions (mirroring munge_marderstein.py so the results-api tabix engine + the BQ
  CHR_STRING_TABLES load treat this file identically to variant_effect):
  - `chrom` is NUMERIC (strip "chr"; X->23, Y->24, M/MT->25) and is column 1 so `tabix -s1 -b2 -e2`
    indexes on (chrom, pos). Non-canonical hg38 contigs (alt/random/scaffold) are DROPPED — the
    platform is primary-assembly only and loads chr as INT64.
  - `variant` is the colon form with numeric chrom and NO "chr" prefix (e.g. "19:662327:G:C"),
    matching the variant_effect `variant` token.
  - numeric-sorted by (chrom, pos), bgzip-compressed, first line a "#"-prefixed tabix header marker,
    every empty/missing cell written as the literal "NA".

Staging (off by default): the .sh wrapper handles GCS upload behind --stage; this script never uploads.

Local validation without any real data:
  python3 scripts/munge_mpra.py --sample
"""

import argparse
import subprocess
import tempfile
from pathlib import Path

import polars as pl

RESOURCE = "siraj_mpra"
DATASET = "siraj_mpra"

DEFAULT_INPUT = str(Path.home() / "mpra_annotation.tsv.gz")

CELL_LINES = ["K562", "HEPG2", "SKNSH", "HCT116", "A549"]
META = "meta"

# canonical LONG column order. `chrom` first + `pos` second so `tabix -s1 -b2 -e2` works. This exact
# set is mirrored by the downstream BQ base table and the results-api mpra vertical. `dataset`/
# `resource` are NOT written here: resource is prepended by the API, and dataset is a per-file
# constant injected at BQ load (load_data.py --const-column), exactly as for the other products.
OUT_COLUMNS = [
    "chrom", "pos", "variant", "ref", "alt", "cohort", "cell_line",
    "emVar", "active",
    "log2Skew", "log2Skew_se", "log2Skew_mlog10p",
    "log2FC", "log2FC_mlog10p",
    "mean_RNA_ref", "mean_RNA_alt",
]

# meta significance uses the raw -log10 p (…_nlogp), NOT the adjusted p (…_padj): padj can be "Inf"
# (adjusted p underflow) and the platform stores a single mlog10p per effect. The two …_padj columns
# are intentionally DROPPED. log2Skew_se comes from the meta SE; per-line rows have no SE (-> NA).
META_SOURCE = {
    "emVar": "emVar_meta",
    "active": "active_meta",
    "log2Skew": "log2Skew_meta",
    "log2Skew_se": "log2Skew_meta_SE",
    "log2Skew_mlog10p": "log2Skew_meta_nlogp",
    "log2FC": "log2FC_meta",
    "log2FC_mlog10p": "log2FC_meta_nlogp",
}

# per-line source suffixes -> harmonized names. Skew_logPadj/logPadj_BF are the per-line significances
# for skew/activity; mean_RNA_ref/alt carry the per-line reporter RNA levels (meta has none).
LINE_SOURCE = {
    "emVar": "emVar",
    "active": "active",
    "log2Skew": "log2Skew",
    "log2Skew_mlog10p": "Skew_logPadj",
    "log2FC": "log2FC",
    "log2FC_mlog10p": "logPadj_BF",
    "mean_RNA_ref": "mean_RNA_ref",
    "mean_RNA_alt": "mean_RNA_alt",
}

CANONICAL_CHROMS = frozenset(str(c) for c in range(1, 26))  # 1..22, X=23, Y=24, M/MT=25


def _numeric_chrom(col: str) -> pl.Expr:
    """Numeric chromosome token from a chr-prefixed seqname: strip 'chr', X=23, Y=24, M/MT=25, else
    the numeric contig. Kept as a string so nulls / non-canonical contigs survive to the drop filter.
    Mirrors CHR_STRING_TO_INT_SQL in genetics-results-db so the tabix seqnames match the chr INT64
    encoding.
    """
    base = pl.col(col).cast(pl.Utf8).str.replace("(?i)^chr", "").str.to_uppercase()
    return (
        pl.when(base == "X").then(pl.lit("23"))
        .when(base == "Y").then(pl.lit("24"))
        .when(base == "M").then(pl.lit("25"))
        .when(base == "MT").then(pl.lit("25"))
        .otherwise(base)
    )


def _block(df: pl.DataFrame, cell_line: str, source: dict[str, str], prefix: str | None) -> pl.DataFrame:
    """One cell_line's slice as OUT_COLUMNS rows.

    `source` maps harmonized name -> source column (per-line names are prefixed "<LINE>."). Any
    OUT_COLUMNS measurement absent from `source` is a null literal (meta has no mean_RNA; per-line has
    no log2Skew_se). Rows where every mapped source value is null are dropped: for a per-line block
    that means the variant was not tested in that line; the meta block is effectively always present.
    """
    src_cols = {name: (f"{prefix}.{col}" if prefix else col) for name, col in source.items()}
    measures = []
    for name in OUT_COLUMNS[7:]:  # measurement columns after the identity + cell_line block
        if name in src_cols:
            measures.append(pl.col(src_cols[name]).alias(name))
        else:
            measures.append(pl.lit(None, dtype=pl.Utf8).alias(name))
    out = df.select(
        pl.col("chrom"), pl.col("pos"), pl.col("variant"), pl.col("ref"), pl.col("alt"),
        pl.col("cohort"), pl.lit(cell_line).alias("cell_line"), *measures,
    )
    present = pl.any_horizontal(pl.col(name).is_not_null() for name in src_cols)
    return out.filter(present).select(OUT_COLUMNS)


def build_long(df: pl.DataFrame) -> pl.DataFrame:
    """WIDE per-variant -> LONG per variant x cell_line, canonical column order, non-canonical dropped."""
    df = df.rename({"#chrom": "chrom"})
    df = df.with_columns(_numeric_chrom("chrom").alias("chrom"))
    df = df.filter(pl.col("chrom").is_in(CANONICAL_CHROMS))
    df = df.with_columns(
        pl.format("{}:{}:{}:{}", pl.col("chrom"), pl.col("pos"), pl.col("ref"), pl.col("alt")).alias("variant")
    )
    frames = [_block(df, META, META_SOURCE, None)]
    frames += [_block(df, line, LINE_SOURCE, line) for line in CELL_LINES]
    long = pl.concat(frames, how="vertical", rechunk=True)
    return long.sort(
        by=[pl.col("chrom").cast(pl.Int64, strict=False), pl.col("pos").cast(pl.Int64, strict=False)],
        nulls_last=True,
    )


def read_wide(path: str) -> pl.DataFrame:
    """Read the WIDE source as all-Utf8 (infer_schema_length=0) so numeric string formats (e.g.
    "6.7138e-01") and TRUE/FALSE booleans pass through verbatim; "NA" becomes null.
    """
    return pl.read_csv(path, separator="\t", infer_schema_length=0, null_values=["NA"])


def write_point(df: pl.DataFrame, output_path: str) -> None:
    """Numeric-sort is already applied; bgzip with a '#'-prefixed header then POINT-index (-s1 -b2 -e2).

    Every empty/missing cell is written as the literal "NA" (residual empty strings coerced to null
    first) so no output cell is ever the empty string.
    """
    df = df.select(OUT_COLUMNS).with_columns(
        pl.when(pl.col(c).cast(pl.Utf8).str.len_chars() == 0).then(None).otherwise(pl.col(c)).alias(c)
        for c in OUT_COLUMNS
    )
    header = ("#" + "\t".join(OUT_COLUMNS) + "\n").encode()
    with open(output_path, "wb") as out_fh:
        proc = subprocess.Popen(["bgzip", "-c"], stdin=subprocess.PIPE, stdout=out_fh)
        assert proc.stdin is not None
        proc.stdin.write(header)
        df.write_csv(proc.stdin, separator="\t", include_header=False, null_value="NA")
        proc.stdin.close()
        rc = proc.wait()
    if rc != 0:
        raise subprocess.CalledProcessError(rc, "bgzip -c")
    subprocess.run(["tabix", "-f", "-s", "1", "-b", "2", "-e", "2", output_path], check=True)
    print(f"  wrote {df.height} rows -> {output_path}")
    print(f"  indexed {output_path}.tbi (POINT: tabix -s1 -b2 -e2)")


def run(args: argparse.Namespace) -> None:
    print(f"Building {DATASET} LONG from {args.input} ...")
    long = build_long(read_wide(args.input))
    assert long.columns == OUT_COLUMNS, f"column order mismatch: {long.columns}"
    counts = long.group_by("cell_line").len().sort("cell_line")
    print(f"  {long.height} LONG rows across {long['cell_line'].n_unique()} cell_line values:")
    for row in counts.iter_rows(named=True):
        print(f"    {row['cell_line']}: {row['len']}")
    write_point(long, args.output)


# ------------------------------------------------------------------------------------------------
# Sample / dry-run: a synthetic WIDE input exercising the pivot, the tested/untested per-line skip,
# the meta<->per-line field harmonization, non-canonical drop, and the POINT index.
# ------------------------------------------------------------------------------------------------
def run_sample() -> None:
    print("=== SAMPLE / DRY-RUN: synthetic WIDE input, no GCS upload ===")
    tmpdir = Path(tempfile.mkdtemp())
    header = ["#chrom", "pos", "variant_id", "variant_hg38", "ref", "alt", "cohort"]
    header += ["emVar_meta", "active_meta", "log2Skew_meta", "log2Skew_meta_SE", "log2Skew_meta_nlogp",
               "log2Skew_meta_padj", "log2FC_meta", "log2FC_meta_SE", "log2FC_meta_nlogp", "log2FC_meta_padj"]
    for line in CELL_LINES:
        header += [f"{line}.emVar", f"{line}.active", f"{line}.log2Skew", f"{line}.Skew_logPadj",
                   f"{line}.log2FC", f"{line}.logPadj_BF", f"{line}.mean_RNA_ref", f"{line}.mean_RNA_alt"]

    def line_block(tested: bool, emv: str) -> list[str]:
        if not tested:
            return ["NA"] * 8
        return [emv, "TRUE", "-0.2", "3.5", "1.1", "5.0", "150.0", "120.0"]

    # v1 (chr1) tested only in K562 (emVar) + HEPG2; v2 (chrX) tested only in A549; v3 on a
    # non-canonical alt contig -> dropped entirely.
    rows = []
    r1 = ["chr1", "63697", "chr1_63697_T_C", "chr1:63697:T:C", "T", "C", "GTEx",
          "TRUE", "TRUE", "-0.18", "0.14", "2.31", "2.01", "1.09", "0.09", "78.7", "73.2"]
    r1 += line_block(True, "TRUE") + line_block(True, "FALSE") + line_block(False, "") + line_block(False, "") + line_block(False, "")
    rows.append(r1)
    r2 = ["chrX", "1000", "chrX_1000_A_G", "chrX:1000:A:G", "A", "G", "UKBB",
          "FALSE", "TRUE", "0.05", "0.10", "0.40", "0.20", "0.90", "0.10", "9.9", "1.0"]
    r2 += line_block(False, "") + line_block(False, "") + line_block(False, "") + line_block(False, "") + line_block(True, "FALSE")
    rows.append(r2)
    r3 = ["chr8_KI270821v1_alt", "500", "x", "x", "A", "T", "control",
          "FALSE", "FALSE", "0.0", "0.1", "0.1", "0.1", "0.1", "0.1", "1.0", "1.0"]
    r3 += line_block(False, "") * 5
    rows.append(r3)

    src = tmpdir / "mpra_wide.tsv"
    src.write_text("\t".join(header) + "\n" + "\n".join("\t".join(r) for r in rows) + "\n")

    long = build_long(read_wide(str(src)))
    assert long.columns == OUT_COLUMNS, long.columns
    recs = long.to_dicts()
    keyed = {(r["variant"], r["cell_line"]): r for r in recs}

    # non-canonical alt contig dropped entirely
    assert not any(r["chrom"] not in CANONICAL_CHROMS for r in recs), recs
    assert not any("KI270821" in (r["variant"] or "") for r in recs), "non-canonical contig leaked"
    # v1: meta + only the two tested per-line rows (K562, HEPG2)
    assert ("1:63697:T:C", "meta") in keyed
    assert ("1:63697:T:C", "K562") in keyed
    assert ("1:63697:T:C", "HEPG2") in keyed
    assert ("1:63697:T:C", "SKNSH") not in keyed, "untested per-line row should be skipped"
    # meta harmonization: se from meta SE, mlog10p from nlogp (NOT padj), mean_RNA null (-> "NA" on write)
    m = keyed[("1:63697:T:C", "meta")]
    assert m["log2Skew_se"] == "0.14" and m["log2Skew_mlog10p"] == "2.31", m
    assert m["log2FC_mlog10p"] == "78.7" and m["mean_RNA_ref"] is None, m
    # per-line harmonization: se null (-> "NA" on write), mlog10p from Skew_logPadj, mean_RNA present
    k = keyed[("1:63697:T:C", "K562")]
    assert k["log2Skew_se"] is None and k["log2Skew_mlog10p"] == "3.5", k
    assert k["mean_RNA_ref"] == "150.0" and k["log2FC_mlog10p"] == "5.0", k
    # chrX -> numeric 23 in both chrom and variant token
    assert ("23:1000:A:G", "A549") in keyed and keyed[("23:1000:A:G", "A549")]["chrom"] == "23", keyed
    print("  pivot + per-line skip + meta/per-line harmonization + chrX->23 + non-canonical drop: OK")
    with pl.Config(tbl_cols=-1, tbl_width_chars=240, fmt_str_lengths=30):
        print(long)

    out_path = str(tmpdir / f"{DATASET}.long.sample.tsv.gz")
    write_point(long, out_path)
    q = subprocess.run(["tabix", out_path, "1:63697-63697"], capture_output=True, text=True, check=True)
    assert q.stdout.strip(), "POINT exact-pos query returned nothing — wrong index mode"
    assert q.stdout.count("\n") >= 3, "expected meta + K562 + HEPG2 rows at 1:63697"
    qx = subprocess.run(["tabix", out_path, "23:1000-1000"], capture_output=True, text=True, check=True)
    assert qx.stdout.strip(), "chrX variant not indexed under numeric seqname 23"
    print("  POINT queries 1:63697-63697 (meta+K562+HEPG2) and 23:1000-1000 (chrX->23): OK")
    print("\n=== SAMPLE OK — no GCS upload performed ===")


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("input", nargs="?", default=DEFAULT_INPUT,
                   help=f"WIDE source TSV (bgzip ok; a gs:// URI works if gcsfuse/gcloud mounts it). Default: {DEFAULT_INPUT}")
    p.add_argument("--output", default=f"{DATASET}.tsv.gz", help="output LONG .tsv.gz path")
    p.add_argument("--sample", "--dry-run", dest="sample", action="store_true",
                   help="run on synthetic input, validate the pivot + POINT index; no upload")
    return p.parse_args()


def main() -> None:
    args = parse_args()
    if args.sample:
        run_sample()
        return
    run(args)


if __name__ == "__main__":
    main()
