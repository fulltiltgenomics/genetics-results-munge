#!/usr/bin/env python3
"""Download and munge the two gene-disease association sources the results-api serves from its
`/gene_disease/{gene}` endpoint: GenCC submissions and the Monarch Initiative KG.

Neither source carries coordinates, so neither output is bgzipped or tabixed -- the API reads the
plain TSV straight out of GCS at startup and holds it in memory. That is also why the outputs stay
uncompressed: `app/services/gene_disease_data.py` reads them through
`pl.read_csv(BytesIO(blob.download_as_bytes()))`, which is given no decompression hint.

Sources (input):
  gencc   https://search.thegencc.org/download/action/submissions-export-tsv
          One row per submission, 30 columns. The export carries no release identifier of its own,
          so the output is versioned by the UTC download date. Fields are quoted only where they
          need to be (an older export quoted every field); both parse the same.
  monarch https://data.monarchinitiative.org/monarch-kg/latest/tsv/
          TWO files, concatenated here into one:
            all_associations/causal_gene_to_disease_association.all.tsv.gz
              biolink:causes + biolink:associated_with_increased_likelihood_of, from OMIM and
              ClinGen.
            gene_associations/gene_disease.noncausal.tsv.gz
              biolink:gene_associated_with_condition (Orphanet) + biolink:contributes_to (OMIM).
          The causal file alone is what the suite served before, and the current release of it
          covers ~26% fewer gene-disease pairs than the release it replaces; the non-causal file
          is taken as well so the endpoint's coverage grows rather than shrinks, which is why
          `predicate` has to reach the API -- causal and non-causal rows are no longer
          distinguishable once they are in the same table. `latest/metadata.yaml` names the KG
          release, and that string versions the output.

Output (one TSV per product, no index):
  gencc   <out-dir>/gencc-submissions-export.<download date>.tsv
          The source's 30 columns unchanged -- byte-identical to the upstream export, as of the
          release this was written against. The parse/serialize round trip is the point: it fails
          here rather than in the API if an export ever stops parsing.
  monarch <out-dir>/monarch-gene_to_disease.<kg release>.tsv
          The sources' 15 columns unchanged, except that `predicate` loses its `biolink:` prefix.
          Rows are dropped when:
            - `subject_category` is not `biolink:Gene`. The non-causal file carries a handful of
              disease-to-disease rows (an OMIM subtype `contributes_to` its parent); left in, they
              reach the endpoint as a gene named "major depressive disorder 1".
            - `subject_taxon` is not human. Both files are human-only today; the filter is here
              because the KG's other gene association files are not.
            - `negated` is set. Nothing is negated today either.
          Exact duplicate rows are then collapsed. Both files repeat rows verbatim, and the API's
          Monarch `uuid` is a composite of subject, object, source and predicate, so a repeat
          reaches the endpoint as two rows claiming the same identifier. The run prints how many
          it collapsed; the same composite is asserted unique before the file is written.

Staging (off by default): the .sh wrapper uploads behind --stage; this script never uploads.
"""

import argparse
import gzip
import re
import shutil
import sys
import urllib.error
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

import polars as pl

GENCC_URL = "https://search.thegencc.org/download/action/submissions-export-tsv"

MONARCH_BASE = "https://data.monarchinitiative.org/monarch-kg/latest"
MONARCH_METADATA_URL = f"{MONARCH_BASE}/metadata.yaml"
MONARCH_FILES = {
    "causal": f"{MONARCH_BASE}/tsv/all_associations/causal_gene_to_disease_association.all.tsv.gz",
    "noncausal": f"{MONARCH_BASE}/tsv/gene_associations/gene_disease.noncausal.tsv.gz",
}

# the API selects from these by name, so a source that drops one fails loudly here instead of
# serving a column of nulls
GENCC_COLUMNS = [
    "uuid", "gene_curie", "gene_symbol", "disease_curie", "disease_title",
    "disease_original_curie", "disease_original_title",
    "classification_curie", "classification_title", "moi_curie", "moi_title",
    "submitter_curie", "submitter_title",
    "submitted_as_hgnc_id", "submitted_as_hgnc_symbol",
    "submitted_as_disease_id", "submitted_as_disease_name",
    "submitted_as_moi_id", "submitted_as_moi_name",
    "submitted_as_submitter_id", "submitted_as_submitter_name",
    "submitted_as_classification_id", "submitted_as_classification_name",
    "submitted_as_date", "submitted_as_public_report_url", "submitted_as_notes",
    "submitted_as_pmids", "submitted_as_assertion_criteria_url",
    "submitted_as_submission_id", "submitted_run_date",
]

MONARCH_COLUMNS = [
    "subject", "subject_label", "subject_category", "subject_taxon", "subject_taxon_label",
    "negated", "predicate", "object", "object_label", "object_category",
    "qualifiers", "publications", "has_evidence",
    "primary_knowledge_source", "aggregator_knowledge_source",
]

HUMAN_TAXON = "NCBITaxon:9606"
GENE_CATEGORY = "biolink:Gene"


def fetch(url: str, dest: Path, download: bool) -> Path:
    """Download `url` to `dest` unless it is already cached there."""
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
        with urllib.request.urlopen(request, timeout=300) as response, tmp.open("wb") as out:
            shutil.copyfileobj(response, out)
    except (urllib.error.URLError, TimeoutError) as err:
        tmp.unlink(missing_ok=True)
        raise SystemExit(
            f"could not fetch {url}: {err}\n"
            f"download it by hand and put it at {dest}, then rerun without --download"
        ) from err
    tmp.rename(dest)
    return dest


def read_tsv(path: Path, expected: list[str]) -> pl.DataFrame:
    """Read a (possibly gzipped) TSV as all-strings and assert its column set."""
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rb") as handle:
        df = pl.read_csv(handle.read(), separator="\t", infer_schema_length=0, quote_char='"')
    if df.columns != expected:
        missing = [c for c in expected if c not in df.columns]
        extra = [c for c in df.columns if c not in expected]
        raise SystemExit(
            f"{path}: unexpected columns\n  missing: {missing}\n  extra: {extra}\n"
            f"  the source layout changed -- update this script and the results-api column map"
        )
    return df


def monarch_release(cache_dir: Path, download: bool) -> str:
    """Read the KG release string out of latest/metadata.yaml.

    Scanned with a regex rather than a YAML parser: this is the only YAML the repo reads and it
    is not worth a dependency for one top-level scalar.
    """
    path = fetch(MONARCH_METADATA_URL, cache_dir / "monarch-metadata.yaml", download)
    text = path.read_text()
    match = re.search(r"^version:\s*'?([0-9]{4}-[0-9]{2}-[0-9]{2})'?\s*$", text, re.MULTILINE)
    if not match:
        raise SystemExit(f"{path}: no top-level `version:` date found")
    return match.group(1)


def munge_gencc(cache_dir: Path, out_dir: Path, download: bool) -> Path:
    source = fetch(GENCC_URL, cache_dir / "gencc-submissions-export.tsv", download)
    df = read_tsv(source, GENCC_COLUMNS)

    version = datetime.now(timezone.utc).strftime("%Y-%m-%d")
    output = out_dir / f"gencc-submissions-export.{version}.tsv"
    write(df, output)

    print(f"  submissions      : {df.height}", file=sys.stderr)
    print(f"  genes            : {df['gene_symbol'].n_unique()}", file=sys.stderr)
    print(f"  diseases         : {df['disease_curie'].n_unique()}", file=sys.stderr)
    print(f"  submitters       : {df['submitter_title'].n_unique()}", file=sys.stderr)
    report_counts(df, "classification_title")
    return output


def munge_monarch(cache_dir: Path, out_dir: Path, download: bool) -> Path:
    version = monarch_release(cache_dir, download)

    frames = []
    for name, url in MONARCH_FILES.items():
        source = fetch(url, cache_dir / f"monarch-{name}.tsv.gz", download)
        frame = read_tsv(source, MONARCH_COLUMNS)
        print(f"  {name:<10}     : {frame.height} rows", file=sys.stderr)
        frames.append(frame)
    df = pl.concat(frames)

    kept = df.filter(
        (pl.col("subject_category") == GENE_CATEGORY)
        & (pl.col("subject_taxon") == HUMAN_TAXON)
        & (pl.col("negated").is_null() | (pl.col("negated").str.to_lowercase() != "true"))
    )
    print(f"  non-gene/non-human/negated dropped: {df.height - kept.height}", file=sys.stderr)

    kept = kept.with_columns(pl.col("predicate").str.strip_prefix("biolink:"))
    deduped = kept.unique(maintain_order=False).sort(
        ["subject_label", "object", "predicate", "primary_knowledge_source"]
    )
    print(f"  duplicate rows collapsed         : {kept.height - deduped.height}", file=sys.stderr)

    output = out_dir / f"monarch-gene_to_disease.{version}.tsv"
    write(deduped, output)

    print(f"  associations     : {deduped.height}", file=sys.stderr)
    print(f"  genes            : {deduped['subject_label'].n_unique()}", file=sys.stderr)
    print(f"  diseases         : {deduped['object'].n_unique()}", file=sys.stderr)
    report_counts(deduped, "predicate")
    report_counts(deduped, "primary_knowledge_source")

    # the API's Monarch uuid is subject|object|source|predicate; a collision there would make two
    # different assertions share an identifier rather than merely repeat
    key = deduped.select(["subject", "object", "primary_knowledge_source", "predicate"])
    if key.height != key.unique().height:
        raise SystemExit(
            f"{key.height - key.unique().height} rows share subject/object/source/predicate "
            f"but differ elsewhere -- the API's uuid would collide"
        )
    return output


def report_counts(df: pl.DataFrame, column: str) -> None:
    counts = df.group_by(column).len().sort("len", descending=True)
    print(f"  by {column}:", file=sys.stderr)
    for row in counts.iter_rows():
        print(f"    {row[0]}: {row[1]}", file=sys.stderr)


def write(df: pl.DataFrame, output: Path) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    df.write_csv(output, separator="\t", null_value="")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--product", choices=["gencc", "monarch"], required=True)
    parser.add_argument("--cache-dir", default=str(Path.home() / "gene_disease_munge" / "cache"))
    parser.add_argument("--out-dir", default=str(Path.home() / "gene_disease_munge" / "out"))
    parser.add_argument("--download", action="store_true", help="fetch anything missing from the cache")
    args = parser.parse_args()

    cache_dir = Path(args.cache_dir)
    out_dir = Path(args.out_dir)
    munge = munge_gencc if args.product == "gencc" else munge_monarch
    output = munge(cache_dir, out_dir, args.download)

    # the wrapper reads this line to learn the versioned file name it has to stage
    print(f"wrote {output}")


if __name__ == "__main__":
    main()
