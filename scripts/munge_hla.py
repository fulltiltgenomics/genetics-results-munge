#!/usr/bin/env python3
"""Munge FinnGen R14 imputed classical HLA allele association results.

Source (FinnGen green library, R14 analysis data, `hla/`):
  - ``summary_stats/<PHENO>.gz``          per-phenotype REGENIE results, one row per HLA allele
  - ``R14_imputed_hla.snpstats``          qctool per-allele QC (AF, HWE, imputation INFO)

The association files are already bgzip+tabix TSVs with the FinnGen sumstats column
layout, but they model an HLA allele as a variant: ``ref`` is the literal string
``<absent>`` and ``alt`` carries the allele name (``A*02:01``). That is faithful to the
dosage test actually run (presence of the allele vs its absence) but it reads as a
nucleotide variant everywhere downstream, and it leaves the HLA gene implicit in the
allele string. This munge therefore rewrites the pair into explicit ``gene``/``allele``
columns and joins the per-allele imputation INFO onto every row, so a returned
association is self-describing about how well its allele was imputed.

Two artifacts come out, because the data has two query axes:
  1. per-phenotype tabix TSVs -- the results-api ``hla`` reads (all alleles for a trait)
  2. one combined TSV        -- the BigQuery ``hla_associations`` load (all traits for an
     allele, which no per-phenotype file can answer)

Phenotype scope: only phenocodes present in the R14 core phenotype metadata are kept.
The source ships ~96 extra ``_WIDE`` endpoints that have no entry in
``finngen_r14_pheno.json``; without metadata they would surface in the API as traits
with no name, case count or category, so they are dropped (see --pheno-metadata).
"""

import argparse
import csv
import gzip
import json
import shutil
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

# Column order of both outputs. The per-phenotype files omit `phenotype` (it is the
# file identity and is prepended by results-api from the request); the combined file
# carries it so BigQuery rows are self-contained.
PER_PHENO_COLUMNS = [
    "#chrom", "pos", "gene", "allele", "pval", "mlogp", "beta", "sebeta",
    "af_alt", "af_alt_cases", "af_alt_controls", "info",
]
COMBINED_COLUMNS = [
    "chrom", "pos", "gene", "allele", "phenotype", "pval", "mlogp", "beta", "sebeta",
    "af_alt", "af_alt_cases", "af_alt_controls", "info",
]

# columns read straight through from the source association file. The case/control
# allele frequencies are absent for quantitative endpoints (BMI, height, Kanta labs),
# exactly as in the core FinnGen sumstats, so they are filled with NA rather than required.
_PASSTHROUGH = [
    "pval", "mlogp", "beta", "sebeta", "af_alt", "af_alt_cases", "af_alt_controls",
]
_REQUIRED_PASSTHROUGH = ["pval", "mlogp", "beta", "sebeta", "af_alt"]

NA = "NA"


def allele_to_gene(allele: str) -> str:
    """``A*02:01`` -> ``HLA-A``. The locus is the allele string up to the ``*``."""
    locus = allele.split("*", 1)[0].strip()
    if not locus:
        raise ValueError(f"cannot derive an HLA gene from allele {allele!r}")
    return f"HLA-{locus}"


def read_snpstats(path: Path) -> dict[str, dict[str, str]]:
    """Parse the qctool ``.snpstats`` sidecar into {allele: {position, info, ...}}.

    qctool writes a block of ``#``-prefixed analysis provenance before the header, so
    the first non-comment line is the header rather than line 1. The allele name is in
    ``rsid`` (``alleleB`` repeats it); ``impute_info`` duplicates ``info``.
    """
    stats: dict[str, dict[str, str]] = {}
    with path.open() as fh:
        header = None
        for line in fh:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if header is None:
                header = fields
                for required in ("rsid", "position", "info"):
                    if required not in header:
                        raise ValueError(
                            f"{path} is missing the {required!r} column; got {header}"
                        )
                continue
            row = dict(zip(header, fields))
            stats[row["rsid"]] = row
    if not stats:
        raise ValueError(f"{path} contained no allele rows")
    return stats


def load_core_phenotypes(path: Path) -> set[str]:
    """Phenocodes from the R14 core phenotype metadata JSON (a list of pheno dicts)."""
    with path.open() as fh:
        entries = json.load(fh)
    codes = {e["phenocode"] for e in entries if e.get("phenocode")}
    if not codes:
        raise ValueError(f"{path} contained no phenocodes")
    return codes


def _fmt_info(raw: str | None) -> str:
    """qctool writes INFO in scientific notation; keep it short but lossless enough."""
    if raw is None or raw == "":
        return NA
    try:
        return f"{float(raw):.5g}"
    except ValueError:
        return NA


def read_associations(
    path: Path, snpstats: dict[str, dict[str, str]], phenotype: str
) -> list[list[str]]:
    """Read one per-phenotype source file into rows of PER_PHENO_COLUMNS order.

    Rows are returned sorted by (pos, allele) so the per-phenotype output is a valid
    tabix input and the combined output has a deterministic within-position order.
    """
    rows: list[list[str]] = []
    with gzip.open(path, "rt") as fh:
        reader = csv.reader(fh, delimiter="\t")
        header = next(reader, None)
        if header is None:
            raise ValueError(f"{path} is empty")
        header[0] = header[0].lstrip("#")
        idx = {name: i for i, name in enumerate(header)}
        for required in ("chrom", "pos", "alt", *_REQUIRED_PASSTHROUGH):
            if required not in idx:
                raise ValueError(
                    f"{path} is missing the {required!r} column; got {header}"
                )
        for fields in reader:
            if not fields or fields[0].startswith("#"):
                continue
            allele = fields[idx["alt"]]
            stat = snpstats.get(allele)
            if stat is None:
                raise ValueError(
                    f"{path}: allele {allele!r} has no row in the snpstats sidecar, so "
                    f"its imputation INFO is unknown -- the sidecar and the association "
                    f"files are out of sync"
                )
            rows.append(
                [
                    fields[idx["chrom"]],
                    fields[idx["pos"]],
                    allele_to_gene(allele),
                    allele,
                    *(
                        fields[idx[c]] if c in idx and idx[c] < len(fields) else NA
                        for c in _PASSTHROUGH
                    ),
                    _fmt_info(stat.get("info")),
                ]
            )
    if not rows:
        raise ValueError(f"{path} contained no association rows for {phenotype}")
    rows.sort(key=lambda r: (int(r[1]), r[3]))
    return rows


def write_tabix_tsv(rows: list[list[str]], columns: list[str], out_path: Path) -> None:
    """Write rows as a bgzipped TSV and tabix-index it (point index, -s1 -b2 -e2)."""
    plain = out_path.with_suffix("")  # strip .gz; bgzip writes <plain>.gz
    with plain.open("w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
        writer.writerow(columns)
        writer.writerows(rows)
    subprocess.run(["bgzip", "-f", str(plain)], check=True)
    subprocess.run(
        ["tabix", "-f", "-s", "1", "-b", "2", "-e", "2", "-c", "#", str(out_path)],
        check=True,
    )


def write_combined(rows: list[list[str]], out_path: Path) -> None:
    """Write the BigQuery load file: gzipped TSV, no index (BQ does not read one)."""
    with gzip.open(out_path, "wt", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
        writer.writerow(COMBINED_COLUMNS)
        writer.writerows(rows)


def munge(
    source_dir: Path,
    snpstats_path: Path,
    pheno_metadata: Path,
    out_dir: Path,
    dataset_id: str,
) -> None:
    snpstats = read_snpstats(snpstats_path)
    core = load_core_phenotypes(pheno_metadata)
    print(f"snpstats alleles: {len(snpstats)}", file=sys.stderr)
    print(f"core R14 phenocodes: {len(core)}", file=sys.stderr)

    sumstats_dir = source_dir / "summary_stats"
    source_files = sorted(
        p for p in sumstats_dir.glob("*.gz") if not p.name.endswith(".tbi")
    )
    if not source_files:
        raise SystemExit(f"no per-phenotype files found under {sumstats_dir}")

    per_pheno_dir = out_dir / "summary_stats"
    per_pheno_dir.mkdir(parents=True, exist_ok=True)

    combined: list[list[str]] = []
    genes: set[str] = set()
    alleles_by_gene: dict[str, set[str]] = defaultdict(set)
    kept = 0
    skipped_no_metadata: list[str] = []

    for src in source_files:
        phenotype = src.name[: -len(".gz")]
        if phenotype not in core:
            skipped_no_metadata.append(phenotype)
            continue
        rows = read_associations(src, snpstats, phenotype)
        write_tabix_tsv(rows, PER_PHENO_COLUMNS, per_pheno_dir / f"{phenotype}.tsv.gz")
        for r in rows:
            genes.add(r[2])
            alleles_by_gene[r[2]].add(r[3])
            # combined layout inserts `phenotype` after `allele`
            combined.append(r[:4] + [phenotype] + r[4:])
        kept += 1
        if kept % 250 == 0:
            print(f"  {kept} phenotypes written", file=sys.stderr)

    # chrom is constant (6) across the MHC, so (pos, allele, phenotype) fully orders it
    combined.sort(key=lambda r: (int(r[1]), r[3], r[4]))
    write_combined(combined, out_dir / f"{dataset_id}.tsv.gz")

    missing = sorted(core - {p.name[: -len(".gz")] for p in source_files})
    print("", file=sys.stderr)
    print(f"phenotypes written:            {kept}", file=sys.stderr)
    print(f"dropped (no R14 metadata):     {len(skipped_no_metadata)}", file=sys.stderr)
    print(f"core phenos with no HLA run:   {len(missing)} {missing}", file=sys.stderr)
    print(f"combined rows:                 {len(combined)}", file=sys.stderr)
    print(f"HLA genes:                     {len(genes)}", file=sys.stderr)
    for gene in sorted(genes):
        print(f"  {gene}: {len(alleles_by_gene[gene])} alleles", file=sys.stderr)
    print(f"\nwrote {out_dir}", file=sys.stderr)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "source_dir",
        type=Path,
        help="Local copy of the release hla/ directory (summary_stats/ + .snpstats)",
    )
    parser.add_argument(
        "--snpstats",
        type=Path,
        help="Path to the .snpstats sidecar (default: the single *.snpstats in source_dir)",
    )
    parser.add_argument(
        "--pheno-metadata",
        type=Path,
        required=True,
        help="R14 core phenotype metadata JSON (finngen_r14_pheno.json); phenocodes "
        "absent from it are dropped",
    )
    parser.add_argument(
        "--output", type=Path, default=Path("hla_out"), help="Output directory"
    )
    parser.add_argument(
        "--dataset-id",
        default="finngen_hla",
        help="Dataset id; names the combined BigQuery load file",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    for tool in ("bgzip", "tabix"):
        if shutil.which(tool) is None:
            raise SystemExit(f"{tool} not found on PATH (run inside the munge image)")

    snpstats_path = args.snpstats
    if snpstats_path is None:
        candidates = sorted(args.source_dir.glob("*.snpstats"))
        if len(candidates) != 1:
            raise SystemExit(
                f"expected exactly one *.snpstats in {args.source_dir}, found "
                f"{[c.name for c in candidates]}; pass --snpstats"
            )
        snpstats_path = candidates[0]

    args.output.mkdir(parents=True, exist_ok=True)
    munge(
        args.source_dir,
        snpstats_path,
        args.pheno_metadata,
        args.output,
        args.dataset_id,
    )


if __name__ == "__main__":
    main()
