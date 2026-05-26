version 1.0

workflow create_pseudo_credible_sets {

    input {
        String docker
        File input_files_list
        String flags
        String output_suffix
        Float r2_high_ld_thres = 0.95
        Boolean filter_hla = true
        String output_file
        File? phenotype_json
        Map[String,String] column_aliases = {}
    }

    meta {
        description: "Create pseudo credible sets from GWAS report files, one task per input file"
    }

    parameter_meta {
        docker: {
            help: "Docker image to use"
        }
        input_files_list: {
            help: "File of filenames (one GCS path per line) to process"
        }
        flags: {
            help: "Flags passed to the script (e.g. '--dataset FinnGen_R13 --no-proximity-filter --min-r2 0.05 --mlog10p-diff 2 --r2-to-lead-thres 0.6')"
        }
        output_suffix: {
            help: "Appended to basename of input file to form output filename (e.g. '.pseudo_cs.tsv')"
        }
        r2_high_ld_thres: {
            help: "r2 threshold for unconditional inclusion in pseudo CS (high-LD exception), default 0.95"
        }
        filter_hla: {
            help: "Keep only the most significant pseudo CS in the HLA region (chr6:25-34Mb), default true"
        }
        output_file: {
            help: "Output filename for the combined bgzipped TSV file (e.g. 'FinnGen_R13_pseudo_credible_sets.tsv.gz')"
        }
        phenotype_json: {
            help: "Optional JSON file with phenotype metadata (array of objects with phenocode and phenostring fields) for mapping trait names"
        }
        column_aliases: {
            help: "Optional map of alternate input column name -> canonical name; merged on top of script defaults. Use when inputs name columns like se/af instead of all_inv_var_meta_sebeta/fg_af_alt."
        }
    }

    Array[String] input_files = read_lines(input_files_list)

    scatter (input_file in input_files) {
        call create_pseudo_cs {
            input:
            docker = docker,
            input_file = input_file,
            flags = flags,
            output_suffix = output_suffix,
            r2_high_ld_thres = r2_high_ld_thres,
            filter_hla = filter_hla,
            phenotype_json = phenotype_json,
            column_aliases = column_aliases
        }
    }

    call collect_results {
        input:
        docker = docker,
        pseudo_cs_files = create_pseudo_cs.output_file,
        output_file = output_file
    }

    output {
        Array[File] pseudo_cs_files = create_pseudo_cs.output_file
        File collected_file = collect_results.collected_file
        File collected_file_tbi = collect_results.collected_file_tbi
    }
}

task create_pseudo_cs {

    input {
        String docker
        File input_file
        String flags
        String output_suffix
        Float r2_high_ld_thres
        Boolean filter_hla
        File? phenotype_json
        Map[String,String] column_aliases
    }

    String output_filename = basename(input_file) + output_suffix

    meta {
        description: "Create pseudo credible sets from a single GWAS report file"
    }

    output {
        File output_file = output_filename
    }

    runtime {
        docker: docker
        memory: "4 GB"
        cpu: 1
        disks: "local-disk 20 HDD"
        preemptible: 2
        noAddress: true
    }

    command <<<

        set -euxo pipefail

        cat <<'PYEOF' > create_pseudo_credible_sets.py
import argparse
import json
import math
import os
import subprocess
import sys
import time
import urllib.request

import polars as pl


NEEDED_COLUMNS = [
    "#CHR", "POS", "REF", "ALT",
    "#variant", "locus_id", "r2_to_lead",
    "all_inv_var_meta_mlogp", "all_inv_var_meta_beta", "all_inv_var_meta_sebeta",
    "fg_af_alt", "most_severe_gene", "most_severe_consequence",
]

# default alternate column names found in some input files -> canonical names
# extra aliases can be merged in per-dataset via --column-aliases
DEFAULT_COLUMN_ALIASES = {
    "#chrom": "#CHR",
    "pos": "POS",
    "ref": "REF",
    "alt": "ALT",
    "mlogp": "all_inv_var_meta_mlogp",
    "beta": "all_inv_var_meta_beta",
    "sebeta": "all_inv_var_meta_sebeta",
    "af_alt": "fg_af_alt",
}

# canonical column -> dtype; aliased source columns inherit the target dtype at read time
CANONICAL_DTYPES = {
    "#CHR": pl.Utf8,
    "POS": pl.Int64,
    "REF": pl.Utf8,
    "ALT": pl.Utf8,
    "#variant": pl.Utf8,
    "locus_id": pl.Utf8,
    "r2_to_lead": pl.Float64,
    "all_inv_var_meta_mlogp": pl.Float64,
    "all_inv_var_meta_beta": pl.Float64,
    "all_inv_var_meta_sebeta": pl.Float64,
    "fg_af_alt": pl.Float64,
    "most_severe_gene": pl.Utf8,
    "most_severe_consequence": pl.Utf8,
}

# extended MHC region (GRCh38) — keep only the most significant pseudo CS here
HLA_REGION_CHR = "6"
HLA_REGION_START = 25_000_000
HLA_REGION_END = 34_000_000


def parse_args():
    parser = argparse.ArgumentParser(
        description="Create pseudo credible sets from GWAS report files"
    )
    parser.add_argument("input", help="Input report file path")
    parser.add_argument("output", help="Output TSV file path")
    parser.add_argument("--dataset", default="FinnGen_R13")
    parser.add_argument("--data-type", default="GWAS")
    parser.add_argument("--cell-type", default="NA")
    parser.add_argument("--mlog10p-diff", type=float, default=3.0)
    parser.add_argument("--r2-to-lead-thres", type=float, default=0.6)
    parser.add_argument("--total-pip", type=float, default=0.99)
    parser.add_argument("--minimum-pip", type=float, default=0.01)
    parser.add_argument("--min-mlog10p", type=float, default=None,
                        help="minimum lead mlog10p to form a pseudo CS (e.g. 7.3 = p<5e-8)")
    parser.add_argument(
        "--no-proximity-filter", action="store_true",
        help="disable 1Mb proximity suppression of secondary leads near very strong signals",
    )
    parser.add_argument("--r2-high-ld-thres", type=float, default=0.95,
                        help="r2 threshold for unconditional inclusion in pseudo CS (high-LD exception)")
    parser.add_argument("--min-r2", type=float, default=None,
                        help="minimum pairwise LD r2 between all variants in a pseudo CS")
    parser.add_argument("--ld-panel",
                        default="gs://r12-data/ld_matrix/finngen_r12_chr{chr}_ld.tsv.gz",
                        help="LD matrix path pattern with {chr} placeholder")
    parser.add_argument("--filter-hla", action=argparse.BooleanOptionalAction, default=True,
                        help="keep only the most significant pseudo CS in the HLA region (chr6:25-34Mb)")
    parser.add_argument("--phenotype-json", default=None,
                        help="JSON file with phenotype metadata for mapping phenocode to phenostring")
    parser.add_argument("--column-aliases", default=None,
                        help="JSON file mapping alt column name -> canonical name; merged on top of defaults")
    return parser.parse_args()


def extract_phenotype(filepath: str) -> str:
    basename = os.path.basename(filepath)
    return basename.split(".")[0]


def _get_gce_access_token() -> str:
    """get access token from GCE metadata server (service account)"""
    req = urllib.request.Request(
        "http://metadata.google.internal/computeMetadata/v1/instance/service-accounts/default/token",
        headers={"Metadata-Flavor": "Google"},
    )
    with urllib.request.urlopen(req) as resp:
        return json.loads(resp.read())["access_token"]


def read_report(filepath: str, extra_aliases: dict | None = None) -> pl.DataFrame:
    aliases = {**DEFAULT_COLUMN_ALIASES, **(extra_aliases or {})}
    # pin both canonical and aliased source columns to canonical dtypes so polars
    # inference doesn't pick e.g. Int64 for a #chr column that lacks "X" rows
    schema_overrides = dict(CANONICAL_DTYPES)
    for src, canon in aliases.items():
        if canon in CANONICAL_DTYPES:
            schema_overrides[src] = CANONICAL_DTYPES[canon]

    # in WDL context, input files are localized by Cromwell
    df = pl.read_csv(
        filepath,
        separator="\t",
        schema_overrides=schema_overrides,
        null_values=["NA"],
        infer_schema_length=100000,
    )

    # normalize alternate column names to canonical names
    renames = {alt: canon for alt, canon in aliases.items() if alt in df.columns and canon not in df.columns}
    if renames:
        df = df.rename(renames)

    missing = [c for c in NEEDED_COLUMNS if c not in df.columns]
    if missing:
        print(f"error: missing columns in {filepath}: {missing}", file=sys.stderr)
        sys.exit(1)

    return df.select(NEEDED_COLUMNS).unique(subset=["#variant", "locus_id"])


def calculate_pip(
    mlog10p_values: list[float],
    total_pip: float,
    min_pip: float,
) -> list[float]:
    n = len(mlog10p_values)
    if n == 0:
        return []
    if n == 1:
        return [total_pip]

    if n * min_pip > total_pip:
        return [total_pip / n] * n

    # weight = 10^mlog10p = 1/p-value
    # use log-space arithmetic to avoid overflow for large mlog10p values:
    # subtract max mlog10p so the largest weight is 1.0, then exponentiate
    max_mlp = max(mlog10p_values)
    log_weights = [(mlp - max_mlp) * math.log(10) for mlp in mlog10p_values]
    weights = [math.exp(lw) for lw in log_weights]

    clamped = [False] * n
    pips = [0.0] * n

    while True:
        unclamped_weight_sum = sum(w for w, c in zip(weights, clamped) if not c)
        clamped_pip_sum = sum(p for p, c in zip(pips, clamped) if c)
        remaining_pip = total_pip - clamped_pip_sum

        if unclamped_weight_sum == 0:
            break

        newly_clamped = False
        for i in range(n):
            if clamped[i]:
                continue
            pips[i] = (weights[i] / unclamped_weight_sum) * remaining_pip
            if pips[i] < min_pip:
                pips[i] = min_pip
                clamped[i] = True
                newly_clamped = True

        if not newly_clamped:
            break

    return pips


def filter_lead_variants(
    df: pl.DataFrame,
    t_threshold: float = 0.02,
    proximity_bp: int = 1_000_000,
) -> pl.DataFrame:
    """Suppress secondary lead variants within proximity_bp of very strong signals.

    For each lead where T = min(0.1, 5/chi2) < t_threshold, remove any less-significant
    lead on the same chromosome within proximity_bp.
    """
    leads = (
        df.filter(pl.col("#variant") == pl.col("locus_id"))
        .unique(subset=["locus_id"])
    )

    beta = leads["all_inv_var_meta_beta"].to_list()
    se = leads["all_inv_var_meta_sebeta"].to_list()
    mlp = leads["all_inv_var_meta_mlogp"].to_list()
    chrs = leads["#CHR"].to_list()
    positions = leads["POS"].to_list()
    locus_ids = leads["locus_id"].to_list()

    # compute T for each lead
    t_values = []
    for b, s in zip(beta, se):
        if b is not None and s is not None and s != 0:
            chi2 = (b / s) ** 2
            t_values.append(min(0.1, 5.0 / chi2) if chi2 > 0 else 0.1)
        else:
            t_values.append(None)

    # sort by mlog10p descending
    order = sorted(range(len(mlp)), key=lambda i: mlp[i] if mlp[i] is not None else -1, reverse=True)

    suppressed = set()
    for idx in order:
        lid = locus_ids[idx]
        if lid in suppressed:
            continue
        t = t_values[idx]
        if t is None or t >= t_threshold:
            continue
        # suppress less-significant leads within proximity_bp on same chr
        for other_idx in order:
            other_lid = locus_ids[other_idx]
            if other_lid == lid or other_lid in suppressed:
                continue
            if chrs[other_idx] != chrs[idx]:
                continue
            if abs(positions[other_idx] - positions[idx]) > proximity_bp:
                continue
            if (mlp[other_idx] or 0) < (mlp[idx] or 0):
                suppressed.add(other_lid)
                print(f"suppressed {other_lid} because it is within {proximity_bp}bp of {lid} and has mlog10p {mlp[other_idx]} < {mlp[idx]} and T {t} < {t_threshold}", file=sys.stderr)

    if suppressed:
        print(
            f"proximity filter: suppressed {len(suppressed)} secondary loci near strong signals",
            file=sys.stderr,
        )
        df = df.filter(~pl.col("locus_id").is_in(list(suppressed)))

    return df


def _query_ld_for_variants(
    ld_panel: str, chrom: str, min_pos: int, max_pos: int, token: str,
    variant_ids: set[str],
    max_retries: int = 3, retry_delay: float = 5.0,
) -> dict:
    """Query tabix for LD between specific variants, streaming to limit memory use."""
    region = f"{chrom}:{min_pos}-{max_pos}"
    path = ld_panel.format(chr=chrom)
    env = {**os.environ, "GCS_OAUTH_TOKEN": token}

    for attempt in range(1, max_retries + 1):
        proc = subprocess.Popen(
            ["tabix", path, region],
            stdout=subprocess.PIPE, stderr=subprocess.PIPE,
            text=True, env=env,
        )
        ld = {}
        for line in proc.stdout:
            fields = line.split("\t")
            v1, v2 = fields[2], fields[3]
            if v1 in variant_ids and v2 in variant_ids:
                r2 = float(fields[5])
                ld[(v1, v2)] = r2
                ld[(v2, v1)] = r2
        stderr = proc.stderr.read()
        if proc.wait() == 0:
            return ld
        print(
            f"warning: tabix failed for {region} (attempt {attempt}/{max_retries}): "
            f"{stderr.strip()}",
            file=sys.stderr,
        )
        if attempt < max_retries:
            time.sleep(retry_delay * attempt)
    return {}


def compute_and_filter_cs_pairwise_r2(
    df: pl.DataFrame,
    ld_panel: str,
    min_r2: float | None = None,
) -> pl.DataFrame:
    """Compute true pairwise min r2 for each CS and update cs_min_r2.

    If min_r2 is given, also drop CSes where any variant pair has r2 < min_r2.
    """
    from itertools import combinations

    token = _get_gce_access_token()

    # group CSes by chromosome; only check multi-variant CSes
    cs_info = []
    for (cs_id,), group in df.group_by("cs_id"):
        variants = group["#variant"].to_list()
        if len(variants) <= 1:
            continue
        chrom = group["#CHR"][0]
        positions = group["POS"].to_list()
        cs_info.append({
            "cs_id": cs_id, "chr": chrom,
            "min_pos": min(positions), "max_pos": max(positions),
            "variants": set(variants),
        })

    if not cs_info:
        return df

    # query per CS — each CS spans a small region, keeping memory use low
    cs_min_r2_map = {}
    failed_cs_ids = set()
    for cs in cs_info:
        ld = _query_ld_for_variants(
            ld_panel, cs["chr"], cs["min_pos"], cs["max_pos"], token,
            variant_ids=cs["variants"],
        )
        pair_min = None
        for v1, v2 in combinations(cs["variants"], 2):
            r2 = ld.get((v1, v2))
            if r2 is None:
                pair_min = None
                break
            if pair_min is None or r2 < pair_min:
                pair_min = r2
        cs_min_r2_map[cs["cs_id"]] = pair_min
        if min_r2 is not None and (pair_min is None or pair_min < min_r2):
            failed_cs_ids.add(cs["cs_id"])

    # update cs_min_r2 column with true pairwise values
    df = df.with_columns(
        pl.struct("cs_id", "cs_min_r2").map_elements(
            lambda row: cs_min_r2_map.get(row["cs_id"], row["cs_min_r2"]),
            return_dtype=pl.Float64,
        ).alias("cs_min_r2")
    )

    if failed_cs_ids:
        print(
            f"LD r2 filter: removed {len(failed_cs_ids)} pseudo CSes with min pairwise r2 < {min_r2}",
            file=sys.stderr,
        )
        df = df.filter(~pl.col("cs_id").is_in(list(failed_cs_ids)))

    return df


def build_pseudo_credible_sets(
    df: pl.DataFrame,
    mlog10p_diff: float,
    r2_thres: float,
    total_pip: float,
    min_pip: float,
    min_mlog10p: float | None = None,
    r2_high_ld_thres: float = 0.95,
) -> pl.DataFrame:
    results = []

    for locus_id, group in df.group_by("locus_id"):
        locus_id = locus_id[0]
        lead = group.filter(pl.col("#variant") == locus_id)
        if len(lead) == 0:
            print(f"warning: no lead variant found for locus {locus_id}, skipping", file=sys.stderr)
            continue

        lead_mlog10p = lead["all_inv_var_meta_mlogp"][0]
        if lead_mlog10p is None:
            print(f"warning: lead variant {locus_id} has null mlog10p, skipping", file=sys.stderr)
            continue

        if min_mlog10p is not None and lead_mlog10p < min_mlog10p:
            print(f"skipping {locus_id} because mlog10p {lead_mlog10p} < {min_mlog10p}", file=sys.stderr)
            continue

        # filter to pseudo CS members
        mlog10p_col = pl.col("all_inv_var_meta_mlogp")
        r2_col = pl.col("r2_to_lead")

        cs_members = group.filter(
            pl.col("all_inv_var_meta_mlogp").is_not_null()
            & pl.col("r2_to_lead").is_not_null()
            & (
                (pl.col("#variant") == locus_id)
                | (r2_col > r2_high_ld_thres)
                | (
                    ((lead_mlog10p - mlog10p_col).abs() < mlog10p_diff)
                    & (r2_col > r2_thres)
                )
            )
        )

        if len(cs_members) == 0:
            continue

        mlog10p_values = cs_members["all_inv_var_meta_mlogp"].to_list()
        pips = calculate_pip(mlog10p_values, total_pip, min_pip)

        cs_size = len(cs_members)

        cs_members = cs_members.with_columns(
            pl.Series("pip", pips),
            pl.lit(locus_id).alias("cs_id"),
            pl.lit(cs_size).alias("cs_size").cast(pl.Int32),
            pl.lit(1.0 if cs_size == 1 else None).cast(pl.Float64).alias("cs_min_r2"),
        )

        results.append(cs_members)

    if not results:
        print("warning: no pseudo credible sets created", file=sys.stderr)
        return pl.DataFrame(schema={
            "#dataset": pl.Utf8, "data_type": pl.Utf8, "trait": pl.Utf8,
            "trait_original": pl.Utf8, "cell_type": pl.Utf8,
            "chr": pl.Int32, "pos": pl.Int64, "ref": pl.Utf8, "alt": pl.Utf8,
            "mlog10p": pl.Float64, "beta": pl.Utf8, "se": pl.Utf8,
            "pip": pl.Float64, "cs_id": pl.Utf8, "cs_size": pl.Int32,
            "cs_min_r2": pl.Float64, "aaf": pl.Float64,
            "most_severe": pl.Utf8, "gene_most_severe": pl.Utf8,
        })

    return pl.concat(results)


def filter_hla_credible_sets(df: pl.DataFrame) -> pl.DataFrame:
    """Keep only the most significant pseudo credible set in the HLA region."""
    if len(df) == 0:
        return df

    leads = df.filter(pl.col("#variant") == pl.col("cs_id"))
    hla_leads = leads.filter(
        (pl.col("#CHR") == HLA_REGION_CHR)
        & (pl.col("POS") >= HLA_REGION_START)
        & (pl.col("POS") <= HLA_REGION_END)
    )

    if len(hla_leads) <= 1:
        return df

    best_cs_id = hla_leads.sort("all_inv_var_meta_mlogp", descending=True)["cs_id"][0]
    remove_cs_ids = [
        cs_id for cs_id in hla_leads["cs_id"].to_list() if cs_id != best_cs_id
    ]

    print(
        f"HLA filter: keeping {best_cs_id}, removing {len(remove_cs_ids)} "
        f"other credible sets in chr{HLA_REGION_CHR}:{HLA_REGION_START}-{HLA_REGION_END}",
        file=sys.stderr,
    )
    return df.filter(~pl.col("cs_id").is_in(remove_cs_ids))


def format_output(
    df: pl.DataFrame,
    dataset: str,
    data_type: str,
    cell_type: str,
    phenotype: str,
    phenotype_mapping: dict | None = None,
) -> pl.DataFrame:
    trait = (phenotype_mapping.get(phenotype, phenotype) if phenotype_mapping else phenotype).replace(" ", "_")
    return (
        df.with_columns(
            pl.lit(dataset).alias("#dataset"),
            pl.lit(data_type).alias("data_type"),
            pl.lit(trait).alias("trait"),
            pl.lit(phenotype).alias("trait_original"),
            pl.lit(cell_type).alias("cell_type"),
            pl.col("#CHR").str.replace("X", "23").cast(pl.Int32).alias("chr"),
            pl.col("POS").alias("pos"),
            pl.col("REF").alias("ref"),
            pl.col("ALT").alias("alt"),
            pl.col("all_inv_var_meta_mlogp").round(4).alias("mlog10p"),
            pl.col("all_inv_var_meta_beta")
                .map_elements(lambda x: f"{x:.3e}", return_dtype=pl.Utf8)
                .alias("beta"),
            pl.col("all_inv_var_meta_sebeta")
                .map_elements(lambda x: f"{x:.3e}", return_dtype=pl.Utf8)
                .alias("se"),
            pl.col("pip").round(4),
            pl.col("cs_min_r2").round(4),
            pl.col("fg_af_alt").alias("aaf"),
            pl.col("most_severe_consequence").alias("most_severe"),
            pl.col("most_severe_gene").alias("gene_most_severe"),
        )
        .select(
            "#dataset", "data_type", "trait", "trait_original", "cell_type",
            "chr", "pos", "ref", "alt", "mlog10p", "beta", "se",
            "pip", "cs_id", "cs_size", "cs_min_r2",
            "aaf", "most_severe", "gene_most_severe",
        )
        .sort("chr", "pos", "ref", "alt")
    )


def main():
    args = parse_args()
    phenotype = extract_phenotype(args.input)

    extra_aliases = None
    if args.column_aliases:
        with open(args.column_aliases) as f:
            extra_aliases = json.load(f)

    df = read_report(args.input, extra_aliases=extra_aliases)

    phenotype_mapping = None
    if args.phenotype_json:
        with open(args.phenotype_json) as f:
            phenotype_mapping = {
                entry["phenocode"]: entry["phenostring"]
                for entry in json.load(f)
            }

    if len(df) == 0:
        print(f"warning: {args.input} has no data rows, writing header only", file=sys.stderr)
        result = build_pseudo_credible_sets(df, mlog10p_diff=0, r2_thres=0, total_pip=0, min_pip=0)
    else:
        if not args.no_proximity_filter:
            df = filter_lead_variants(df)

        df = build_pseudo_credible_sets(
            df,
            mlog10p_diff=args.mlog10p_diff,
            r2_thres=args.r2_to_lead_thres,
            total_pip=args.total_pip,
            min_pip=args.minimum_pip,
            min_mlog10p=args.min_mlog10p,
            r2_high_ld_thres=args.r2_high_ld_thres,
        )

        if args.filter_hla:
            df = filter_hla_credible_sets(df)

        df = compute_and_filter_cs_pairwise_r2(df, args.ld_panel, min_r2=args.min_r2)

        result = format_output(df, args.dataset, args.data_type, args.cell_type, phenotype, phenotype_mapping)

    result.write_csv(args.output, separator="\t", null_value="NA")
    print(
        f"wrote {len(result)} variants in {result['cs_id'].n_unique()} "
        f"pseudo credible sets to {args.output}"
    )


if __name__ == "__main__":
    main()
PYEOF

        python3 create_pseudo_credible_sets.py ~{input_file} ~{output_filename} ~{flags} --r2-high-ld-thres ~{r2_high_ld_thres} ~{if filter_hla then "--filter-hla" else "--no-filter-hla"} ~{"--phenotype-json " + phenotype_json} --column-aliases ~{write_json(column_aliases)}

    >>>
}

task collect_results {

    input {
        String docker
        Array[File] pseudo_cs_files
        String output_file
    }

    meta {
        description: "Combine per-trait pseudo credible set files into one sorted, bgzipped, tabix-indexed file"
    }

    output {
        File collected_file = output_file
        File collected_file_tbi = output_file + ".tbi"
    }

    runtime {
        docker: docker
        memory: "4 GB"
        cpu: 1
        disks: "local-disk 50 HDD"
        preemptible: 2
    }

    command <<<

        set -euxo pipefail

        # list input files for sort
        for f in ~{sep=' ' pseudo_cs_files}; do
            echo "$f"
        done | tr '\n' '\0' > merge_these

        # merge-sort by chr(6), pos(7), ref(8), alt(9), trait(3)
        # header lines sort to top via -g treating "chr" as 0; uniq removes duplicates
        sort -m -T . --files0-from=merge_these -k6,6g -k7,7g -k8,8 -k9,9 -k3,3 \
            | uniq \
            | bgzip > ~{output_file}
        tabix -s6 -b7 -e7 ~{output_file}

        # sanity check row counts
        n_lines_orig=0
        for f in ~{sep=' ' pseudo_cs_files}; do
            n_lines_orig=$((n_lines_orig + $(tail -n+2 "$f" | wc -l)))
        done
        n_lines_output=$(zcat ~{output_file} | tail -n+2 | wc -l)
        if [ "$n_lines_orig" -ne "$n_lines_output" ]; then
            echo "ERROR: row count mismatch: expected $n_lines_orig, got $n_lines_output" >&2
            exit 1
        fi

    >>>
}
