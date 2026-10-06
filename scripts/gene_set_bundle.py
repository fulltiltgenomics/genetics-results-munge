#!/usr/bin/env python3
"""Harmonise a downloaded gene set release (build_gene_sets.py's `gene_sets.<YYYYMMDD>/`) to
unversioned Ensembl gene ids and write the parquet bundle the sandbox reads; optionally stage it.

Second stage of build_gene_sets.py, kept in its own module because its inputs are local files
(the raw release dir plus the HGNC and GENCODE files below) and it never downloads.

Inputs
  --raw-dir   a gene_sets.<YYYYMMDD>/ from build_gene_sets.py. Every file's sha256 is checked
              against its manifest and every GMT is re-run through build_gene_sets.check_gmt
              (the KEGG/BIOCARTA guard) before anything is read.
  --hgnc      hgnc_complete_set.txt, the same staged copy gene_annotations_v is built from
              (gs://daly-genetics-results/mapping_files/hgnc_complete_set.txt, loaded by
              genetics-results-db scripts/load_gene_annotations.sh). Tab-separated, every row
              status "Approved"; pipe-delimited prev/alias fields are sometimes double-quoted,
              so quoting is off and quotes are stripped. The file carries no version, hence the
              required --hgnc-version (the `hgnc_version` gene_annotations_v rows carry).
  --gencode   gencode.v<N>.annotation.genes.tsv from the same mapping_files/ folder. `chrom` is
              already the integer encoding (23=X, 24=Y, 26=MT) and passes through unchanged, as
              in gene_annotations_v; gene_id is versioned and stripped.

Mapping to Ensembl: a member resolves to one HGNC gene, then to that gene's Ensembl id.
  MSigDB (entrez GMTs)  entrez_id == the member.
  Reactome (symbols)    tiers approved symbol, prev_symbol, alias_symbol, in that order. The
                        member resolves at the FIRST tier in which any HGNC gene lists it and
                        never falls through: if two genes list it in that tier, or the one gene
                        has no Ensembl id (or an ambiguous one), it is unmapped. Every index is
                        built over all HGNC rows, Ensembl id or not, so a symbol naming a gene
                        without an Ensembl id cannot fall through to another gene's alias.
                        Matching is case-sensitive: HGNC symbols are upper case, and Reactome's
                        lower-case members are viral or bacterial genes (gag, vif, fliC).
A member that does not map is dropped and counted in its set's n_unmapped. Two members mapping
to one Ensembl id collapse to one row, so n_mapped counts distinct Ensembl ids.
A collection whose membership-level mapping rate (mapped (set, member) pairs over all pairs) is
below MIN_MAPPING_RATE for its source refuses the whole bundle. An Ensembl id two HGNC genes
claim is ambiguous: it leaves the universe and its members count as unmapped.

Set ids: MSigDB standard names (HALLMARK_*, GOBP_*, ...) and Reactome stable ids (R-HSA-<n>).
The two cannot collide (MSigDB names are [A-Z0-9_]+), and uniqueness across the bundle is
asserted anyway.

Output: <out-dir>/gene_sets_bundle.<raw release>.hgnc<YYYYMMDD>.gencode<N>/ with sets.parquet,
members.parquet, genes.parquet and manifest.json; no raw GMTs (they stay in --raw-dir, named
by sha256 in the manifest). Assembled under `.partial` and renamed once complete; never
overwritten. genes.parquet is every HGNC gene with an Ensembl id, column names as in
gene_annotations_v (chr, gene_start, gene_end, strand; null where GENCODE has no such gene).

--stage copies the bundle to gs://daly-genetics-results/gene_sets/<bundle release>/ (daly only:
c5.hpo is LicenseRef-HPO, which the finngen profile does not allow), refusing a non-empty
remote path. Each object is copied on its own with --if-generation-match=0 (create only), the
parquet files first and manifest.json last, so a remote path without manifest.json is a partial
upload: consumers must require manifest.json and check every file's sha256 against it.
--dry-run prints the target and the copies without running them.
"""

import argparse
import json
import shutil
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

import polars as pl

from build_gene_sets import check_gmt, sha256_of, write_manifest

# Membership-level floor per source. Measured on MSigDB 2026.1.Hs / Reactome 97 against HGNC
# 2026-06-02: every MSigDB collection >= 0.9994, Reactome 0.976. Reactome's unmapped members are
# mostly not human genes (viral and bacterial proteins, rRNAs, "SARS coronavirus, complete
# genome", generic IGHV), plus 166 pairs whose HGNC gene has no Ensembl id. The floors leave
# room for HGNC drift; a wrong id column or a mis-parsed file maps near 0 and cannot pass.
MIN_MAPPING_RATE = {"msigdb": 0.99, "reactome": 0.95}

STAGE_BUCKET = "gs://daly-genetics-results/gene_sets"
BUNDLE_FILES = ("sets.parquet", "members.parquet", "genes.parquet")


def _clean(value: str | None) -> str | None:
    if value is None:
        return None
    value = value.replace('"', "").strip()
    return value or None


def read_hgnc(path: Path) -> pl.DataFrame:
    hgnc = pl.read_csv(path, separator="\t", infer_schema_length=0, null_values=[""],
                       quote_char=None)
    cols = ["hgnc_id", "symbol", "locus_group", "locus_type", "status", "entrez_id",
            "ensembl_gene_id", "prev_symbol", "alias_symbol"]
    missing = set(cols) - set(hgnc.columns)
    if missing:
        raise SystemExit(f"{path}: missing HGNC columns {sorted(missing)}")
    hgnc = hgnc.select(cols).with_columns(pl.col(c).map_elements(_clean, return_dtype=pl.String)
                                          for c in cols)
    statuses = set(hgnc["status"].unique())
    if statuses != {"Approved"}:
        raise SystemExit(f"{path}: expected only Approved rows, found statuses {statuses}")
    return hgnc


def read_gencode(path: Path) -> pl.DataFrame:
    gencode = pl.read_csv(path, separator="\t", infer_schema_length=0, null_values=[""])
    gencode = gencode.select(
        pl.col("gene_id").str.replace(r"\..*$", "").alias("ensembl_gene_id"),
        pl.col("chrom").cast(pl.Int64).alias("chr"),
        pl.col("gene_start").cast(pl.Int64),
        pl.col("gene_end").cast(pl.Int64),
        pl.col("gene_strand").alias("strand"),
        pl.col("gene_type").alias("gencode_gene_type"),
    )
    dup = gencode.filter(pl.col("ensembl_gene_id").is_duplicated())
    if dup.height:
        raise SystemExit(f"{path}: {dup.height} rows share a version-stripped gene id, "
                         f"e.g. {dup['ensembl_gene_id'][0]}")
    return gencode


def gene_universe(hgnc: pl.DataFrame, gencode: pl.DataFrame) -> pl.DataFrame:
    # an Ensembl id two HGNC genes claim (HGNC 2026-06-02: three pseudogene/lncRNA pairs) is
    # ambiguous, so it leaves the universe and its members count as unmapped
    genes = hgnc.filter(pl.col("ensembl_gene_id").is_not_null()
                        & ~pl.col("ensembl_gene_id").is_duplicated())
    return (genes.join(gencode, on="ensembl_gene_id", how="left")
            .select("ensembl_gene_id", "hgnc_id", "symbol", "entrez_id", "locus_group",
                    "locus_type", "gencode_gene_type", "chr", "gene_start", "gene_end",
                    "strand", "prev_symbol", "alias_symbol")
            .sort("ensembl_gene_id"))


def _index(pairs) -> dict[str, set[str]]:
    seen: dict[str, set[str]] = defaultdict(set)
    for key, value in pairs:
        if key:
            seen[key].add(value)
    return seen


def build_resolvers(hgnc: pl.DataFrame, genes: pl.DataFrame) -> dict[str, dict[str, str]]:
    """Per id type, member -> Ensembl id. Tiers are indexed over every HGNC row (`hgnc`) and a
    member is decided at the first tier that names it; only then is its one gene looked up in
    the Ensembl-bearing universe (`genes`)."""
    ensembl_of = dict(zip(genes["hgnc_id"], genes["ensembl_gene_id"]))
    rows = hgnc.select("hgnc_id", "symbol", "entrez_id", "prev_symbol", "alias_symbol").to_dicts()

    def split(field):
        return ((s.strip(), r["hgnc_id"]) for r in rows
                for s in (r[field] or "").split("|") if s.strip())

    def resolve(*tiers):
        out: dict[str, str] = {}
        decided: set[str] = set()
        for tier in tiers:
            for key, ids in tier.items():
                if key in decided:
                    continue
                decided.add(key)
                if len(ids) == 1 and (hid := next(iter(ids))) in ensembl_of:
                    out[key] = ensembl_of[hid]
        return out

    entrez = resolve(_index((r["entrez_id"], r["hgnc_id"]) for r in rows))
    symbol = resolve(_index((r["symbol"], r["hgnc_id"]) for r in rows),
                     _index(split("prev_symbol")), _index(split("alias_symbol")))
    return {"entrez": entrez, "symbol": symbol}


def parse_gmt(path: Path):
    for line in path.read_text(encoding="utf-8").split("\n"):
        if line:
            name, link, *members = line.split("\t")
            yield name, link, list(dict.fromkeys(m.strip() for m in members if m.strip()))


def harmonise(collection: str, gmt: Path, resolver: dict[str, str], symbol_of: dict[str, str],
              *, reactome: bool, min_rate: float):
    sets, members = [], []
    pairs = mapped_pairs = 0
    for name, link, given in parse_gmt(gmt):
        ensg = [resolver[m] for m in given if m in resolver]
        distinct = list(dict.fromkeys(ensg))
        set_id = link if reactome else name
        sets.append({"set_id": set_id, "collection": collection, "name": name,
                     "url_or_stable_id": link, "n_given": len(given),
                     "n_mapped": len(distinct), "n_unmapped": len(given) - len(ensg)})
        members.extend({"set_id": set_id, "ensembl_gene_id": e, "symbol": symbol_of[e]}
                       for e in distinct)
        pairs += len(given)
        mapped_pairs += len(ensg)
    rate = mapped_pairs / pairs if pairs else 0.0
    stats = {"n_sets": len(sets), "member_pairs": pairs, "mapped_pairs": mapped_pairs,
             "mapping_rate": round(rate, 6), "min_mapping_rate": min_rate}
    return sets, members, stats


def verified_gmts(raw_dir: Path, manifest: dict) -> list[tuple[dict, Path]]:
    """Each collection's GMT, after checking every recorded checksum and re-running the guard."""
    out = []
    for entry in manifest["collections"]:
        checks = [(entry["file"], entry["sha256"])]
        if "gmt_file" in entry:
            checks.append((entry["gmt_file"], entry["gmt_sha256"]))
        for name, digest in checks:
            if sha256_of(raw_dir / name) != digest:
                raise SystemExit(f"{raw_dir / name}: sha256 differs from the raw manifest")
        gmt = raw_dir / entry.get("gmt_file", entry["file"])
        check_gmt(gmt, entry["collection"])
        out.append((entry, gmt))
    return out


def bundle_release_id(raw_release: str, hgnc_version: str, gencode_version: str) -> str:
    return f"{raw_release}.hgnc{hgnc_version.replace('-', '')}.gencode{gencode_version}"


def build_bundle(raw_dir: Path, hgnc_path: Path, gencode_path: Path, out_dir: Path, *,
                 hgnc_version: str, gencode_version: str) -> Path:
    raw_manifest = json.loads((raw_dir / "manifest.json").read_text(encoding="utf-8"))
    release = bundle_release_id(raw_manifest["release_id"], hgnc_version, gencode_version)
    final = out_dir / f"gene_sets_bundle.{release}"
    if final.exists():
        raise SystemExit(f"{final} already exists; refusing to overwrite (versions are immutable)")
    work = final.with_name(final.name + ".partial")
    if work.exists():
        raise SystemExit(f"{work} exists from an interrupted run; remove it by hand")

    gmts = verified_gmts(raw_dir, raw_manifest)
    hgnc = read_hgnc(hgnc_path)
    genes = gene_universe(hgnc, read_gencode(gencode_path))
    resolvers = build_resolvers(hgnc, genes)
    symbol_of = dict(zip(genes["ensembl_gene_id"], genes["symbol"]))

    all_sets, all_members, mapping, refused = [], [], {}, []
    for entry, gmt in gmts:
        reactome = entry["source"] == "reactome"
        sets, members, stats = harmonise(
            entry["collection"], gmt, resolvers["symbol" if reactome else "entrez"], symbol_of,
            reactome=reactome, min_rate=MIN_MAPPING_RATE[entry["source"]])
        stats["member_id_type"] = "HGNC symbol" if reactome else "NCBI Entrez gene id"
        mapping[entry["collection"]] = stats
        # compared on counts, not on the rounded rate written to the manifest
        if stats["mapped_pairs"] < stats["min_mapping_rate"] * stats["member_pairs"]:
            refused.append(f"{entry['collection']} {stats['mapping_rate']:.4f} "
                           f"< {stats['min_mapping_rate']}")
        all_sets += sets
        all_members += members
    for name, stats in mapping.items():
        print(f"  {name:20s} {stats['n_sets']:6d} sets  mapping rate {stats['mapping_rate']:.4f}",
              file=sys.stderr)
    if refused:
        raise SystemExit(f"mapping rate below its floor: {', '.join(refused)}")

    sets_df = pl.DataFrame(all_sets)
    dup = sets_df.filter(pl.col("set_id").is_duplicated())
    if dup.height:
        raise SystemExit(f"set ids not unique across collections: {dup['set_id'].unique()[:5].to_list()}")
    members_df = pl.DataFrame(all_members, schema={"set_id": pl.String,
                                                   "ensembl_gene_id": pl.String,
                                                   "symbol": pl.String})
    tables = {
        "sets.parquet": sets_df.sort("set_id"),
        "members.parquet": members_df.sort("set_id", "ensembl_gene_id"),
        "genes.parquet": genes.drop("prev_symbol", "alias_symbol"),
    }

    work.mkdir(parents=True)
    try:
        files = {}
        for name, df in tables.items():
            df.write_parquet(work / name, compression="zstd")
            files[name] = {"sha256": sha256_of(work / name), "bytes": (work / name).stat().st_size,
                           "rows": df.height}
        manifest = dict(raw_manifest)
        manifest["bundle"] = {
            "release_id": release,
            "raw_release_id": raw_manifest["release_id"],
            "raw_manifest_sha256": sha256_of(raw_dir / "manifest.json"),
            "hgnc": {"file": hgnc_path.name, "version": hgnc_version,
                     "sha256": sha256_of(hgnc_path)},
            "gencode": {"file": gencode_path.name, "version": gencode_version,
                        "sha256": sha256_of(gencode_path)},
            "gene_id": "unversioned Ensembl gene id via HGNC",
            "set_id_scheme": "MSigDB standard name; Reactome stable id (R-HSA-<n>)",
            "mapping": mapping,
            "files": files,
        }
        write_manifest(work / "manifest.json", manifest)
    except BaseException:
        shutil.rmtree(work, ignore_errors=True)
        raise
    work.rename(final)
    return final


def verify_bundle(bundle: Path) -> dict:
    manifest = json.loads((bundle / "manifest.json").read_text(encoding="utf-8"))
    files = manifest["bundle"]["files"]
    if sorted(files) != sorted(BUNDLE_FILES):
        raise SystemExit(f"{bundle}: manifest lists {sorted(files)}, expected {sorted(BUNDLE_FILES)}")
    for name, meta in files.items():
        if sha256_of(bundle / name) != meta["sha256"]:
            raise SystemExit(f"{bundle / name}: sha256 differs from the bundle manifest")
    return manifest


def stage(bundle: Path, *, dry_run: bool, run=subprocess.run) -> str:
    manifest = verify_bundle(bundle)
    target = f"{STAGE_BUCKET}/{manifest['bundle']['release_id']}/"
    # one object per copy: the generation-0 precondition makes each a create-only write that
    # fails rather than skips, and manifest.json goes last so its presence marks completion
    copies = [["gcloud", "storage", "cp", "--if-generation-match=0", str(bundle / name),
               target + name] for name in (*BUNDLE_FILES, "manifest.json")]
    if dry_run:
        print(f"would stage to {target}")
        for copy in copies:
            print(f"  {' '.join(copy)}")
        return target
    listing = run(["gcloud", "storage", "ls", target], capture_output=True, text=True)
    if listing.returncode == 0 and listing.stdout.strip():
        if target + "manifest.json" not in listing.stdout.split():
            raise SystemExit(f"{target} has objects but no manifest.json: a partial upload. "
                             f"Delete it by hand (gcloud storage rm -r {target}) and re-run")
        raise SystemExit(f"{target} already exists; refusing to overwrite (versions are immutable)")
    if listing.returncode != 0 and "matched no objects" not in listing.stderr:
        raise SystemExit(f"could not list {target}: {listing.stderr.strip()}")
    for copy in copies:
        if run(copy).returncode != 0:
            raise SystemExit(f"copy failed: {' '.join(copy)}. Without manifest.json, {target} "
                             f"is a partial upload: delete it by hand (gcloud storage rm -r "
                             f"{target}) and re-run")
    print(f"staged {target}")
    return target


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--raw-dir", type=Path, help="gene_sets.<YYYYMMDD>/ to bundle")
    parser.add_argument("--hgnc", type=Path)
    parser.add_argument("--hgnc-version", help="HGNC release date, YYYY-MM-DD")
    parser.add_argument("--gencode", type=Path)
    parser.add_argument("--gencode-version", default="49")
    parser.add_argument("--out-dir", type=Path)
    parser.add_argument("--bundle", type=Path, help="stage this existing bundle instead of building")
    parser.add_argument("--stage", action="store_true")
    parser.add_argument("--dry-run", action="store_true", help="with --stage: print, do not copy")
    args = parser.parse_args()
    bundle = args.bundle
    if bundle is None:
        needed = ("raw_dir", "hgnc", "hgnc_version", "gencode", "out_dir")
        if any(getattr(args, n) is None for n in needed):
            parser.error("building needs --raw-dir --hgnc --hgnc-version --gencode --out-dir")
        args.out_dir.mkdir(parents=True, exist_ok=True)
        bundle = build_bundle(args.raw_dir, args.hgnc, args.gencode, args.out_dir,
                              hgnc_version=args.hgnc_version, gencode_version=args.gencode_version)
        print(f"wrote {bundle}")
    if args.stage:
        stage(bundle, dry_run=args.dry_run)


if __name__ == "__main__":
    main()
