#!/usr/bin/env python3
"""
Measure how many credible set rows of the non-FinnGen resources would get a consequence
annotation from gnomAD, compared to what they carry today.

Presence of the variant in a gnomAD sites file is used as the measure, so the result is an
upper bound on annotation coverage: a present variant can still lack a consequence.

Usage: measure_gnomad_coverage.py <workdir> <gnomad_sites_tsv_gz>...

The sites files need #CHROM/POS/REF/ALT as their first four columns, chromosomes spelled
chr1..chr22/chrX, and may be gs:// paths. They are streamed once each and never held in
memory; only the credible set variant keys are.
"""

import os
import subprocess
import sys

import polars as pl
from google.cloud import bigquery
from google.oauth2.credentials import Credentials

PROJECT = "phewas-development"
TABLE = f"{PROJECT}.genetics_results.credible_sets"

# the resources whose consequence columns are to be replaced. UKB-PPP is reported but kept
# out of the pooled figure: its annotation is passed through from its source files rather
# than joined from the FinnGen file, so it is not what the pooled gain is a statement about
POOLED = ["eqtl_catalogue", "open_targets", "nmr_ukbb_est", "pgc_scz_finemap"]
RESOURCE_CASE = """
    CASE
        WHEN STARTS_WITH(resource, 'qtd') THEN 'eqtl_catalogue'
        WHEN resource = 'open_targets' THEN 'open_targets'
        WHEN resource = 'nmr_ukbb_est' THEN 'nmr_ukbb_est'
        WHEN dataset = 'PGC_SCZ_2022' THEN 'pgc_scz_finemap'
        WHEN resource = 'ukbb' THEN 'ukb_ppp'
    END
"""

# the key set is the only thing held in memory while gnomAD streams past; refuse to start
# rather than let awk grow without bound if the table turns out far larger than expected
MAX_KEYS = 20_000_000
AWK_VMEM_KB = 8 * 1024 * 1024

AWK_PROGRAM = r"""
NR == FNR { keys[$1] = 1; next }
FNR == 1 { for (i = 1; i <= NF; i++) if ($i ~ /^AN_/) an[i] = 1; next }
{
    key = $1 ":" $2 ":" $3 ":" $4
    if (key in keys) {
        total = 0
        for (i in an) total += $i
        print key "\t" total
    }
}
"""


def export_variants(path: str) -> None:
    token = subprocess.run(
        ["gcloud", "auth", "print-access-token"], check=True, capture_output=True, text=True
    ).stdout.strip()
    client = bigquery.Client(project=PROJECT, credentials=Credentials(token))
    sql = f"""
        SELECT {RESOURCE_CASE} AS resource, chr, pos, ref, alt,
               COUNT(*) AS n_rows,
               COUNTIF(most_severe IS NOT NULL) AS n_annotated
        FROM `{TABLE}`
        GROUP BY 1, 2, 3, 4, 5
        HAVING resource IS NOT NULL
    """
    pl.from_arrow(client.query(sql).to_arrow()).write_parquet(path)


def stream_matches(sites: str, keys_path: str, out_path: str) -> subprocess.Popen:
    cat = f"gcloud storage cat {sites}" if sites.startswith("gs://") else f"cat {sites}"
    command = (
        f"set -o pipefail; ulimit -v {AWK_VMEM_KB}; "
        f"{cat} | zcat | mawk -F'\\t' '{AWK_PROGRAM}' {keys_path} - > {out_path}"
    )
    return subprocess.Popen(["bash", "-c", command])


def main(workdir: str, sites_files: list[str]) -> None:
    os.makedirs(workdir, exist_ok=True)
    variants_path = os.path.join(workdir, "cs_variants.parquet")
    if not os.path.exists(variants_path):
        export_variants(variants_path)
    variants = pl.read_parquet(variants_path).with_columns(
        pl.concat_str(
            pl.lit("chr"),
            pl.col("chr").cast(pl.String).replace("23", "X"),
            pl.lit(":"),
            pl.col("pos").cast(pl.String),
            pl.lit(":"),
            pl.col("ref"),
            pl.lit(":"),
            pl.col("alt"),
        ).alias("key")
    )

    keys = variants.select("key").unique()
    if keys.height > MAX_KEYS:
        sys.exit(f"{keys.height} distinct variants exceeds the cap of {MAX_KEYS}")
    keys_path = os.path.join(workdir, "keys.txt")
    keys.write_csv(keys_path, include_header=False)
    print(f"{keys.height} distinct variants over {variants['n_rows'].sum()} rows", flush=True)

    match_paths = [os.path.join(workdir, f"matches_{i}.tsv") for i in range(len(sites_files))]
    procs = [stream_matches(s, keys_path, p) for s, p in zip(sites_files, match_paths)]
    for sites, proc in zip(sites_files, procs):
        if proc.wait() != 0:
            sys.exit(f"streaming {sites} failed with exit code {proc.returncode}")

    matches = pl.concat(
        [
            pl.read_csv(p, separator="\t", has_header=False, new_columns=["key", "an"])
            for p in match_paths
        ]
    )
    present = matches.group_by("key").agg(pl.col("an").max().alias("an"))
    joined = variants.join(present, on="key", how="left").with_columns(
        pl.col("an").is_not_null().alias("present"),
        (pl.col("an").fill_null(0) > 0).alias("present_called"),
    )

    def summarise(frame: pl.DataFrame, label: str) -> dict:
        rows = frame["n_rows"].sum()
        annotated = frame["n_annotated"].sum()
        in_gnomad = frame.filter("present")["n_rows"].sum()
        called = frame.filter("present_called")["n_rows"].sum()
        lost = frame.filter(~pl.col("present"))["n_annotated"].sum()
        recovered = frame.filter("present").select(
            (pl.col("n_rows") - pl.col("n_annotated")).sum()
        ).item()
        return {
            "resource": label,
            "rows": rows,
            "variants": frame["key"].n_unique(),
            "pct_annotated_today": round(100 * annotated / rows, 2),
            "pct_in_gnomad": round(100 * in_gnomad / rows, 2),
            "gain_points": round(100 * (in_gnomad - annotated) / rows, 2),
            "pct_in_gnomad_an_gt_0": round(100 * called / rows, 2),
            "pct_rows_annotated_to_na": round(100 * lost / rows, 2),
            "pct_rows_na_to_present": round(100 * recovered / rows, 2),
        }

    report = [summarise(f, name) for (name,), f in joined.partition_by("resource", as_dict=True).items()]
    report.append(summarise(joined.filter(pl.col("resource").is_in(POOLED)), "POOLED"))
    out = pl.DataFrame(report).sort("rows", descending=True)
    out.write_csv(os.path.join(workdir, "coverage.tsv"), separator="\t")
    with pl.Config(tbl_cols=-1, tbl_rows=-1, tbl_width_chars=250):
        print(out)


if __name__ == "__main__":
    if len(sys.argv) < 3:
        sys.exit(__doc__)
    main(sys.argv[1], sys.argv[2:])
