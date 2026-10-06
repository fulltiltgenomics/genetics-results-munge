#!/usr/bin/env python3
"""Download the gene set collections the sandbox's enrichment analyses are built from, guard
them against unlicensable sets, and write them with a manifest into a new versioned directory.

This is the download stage only. Gene ID harmonisation and the parquet bundle are a later stage
(docs/gene-sets.md).

Sources (input), and nothing else -- the list is `MSIGDB_COLLECTIONS` plus Reactome:
  MSigDB  https://data.broadinstitute.org/gsea-msigdb/msigdb/release/<release>/
          h.all, c2.cp.wikipathways, c2.cp.pid, c5.go.bp, c5.go.cc, c5.go.mf, c5.hpo, as
          `<collection>.v<release>.entrez.gmt`. Per-subcollection files are taken on purpose:
          c2.cp, c2.all and the msigdb.* bundles contain KEGG_ and BIOCARTA_ sets, whose terms
          do not allow redistribution to every deployment.
  Reactome  https://download.reactome.org/<n>/ReactomePathways.gmt.zip
          The only Reactome source (MSigDB's c2.cp.reactome is a derived copy of it). HGNC
          symbols, not Entrez ids. <n> is read from the ContentService at download time and
          pinned in the URL, so the manifest names a release that stays fetchable.

Output: <out-dir>/gene_sets.<YYYYMMDD>/ holding the raw files and manifest.json (per collection:
source URL, upstream release, download date, sha256, SPDX licence, attribution). The directory
is named after the UTC build date and is never overwritten; it is assembled under a `.partial`
name and renamed only once every guard has passed, so a refused build leaves no versioned dir.

Licences come from `LICENCES`, committed here per (source, upstream release) and never read from
the downloaded data; a release with no entry there refuses to build. c5.hpo carries the HPO's own
terms as `LicenseRef-HPO`, not MSigDB's CC BY 4.0.

Staging (off by default): this script never uploads.
"""

import argparse
import hashlib
import json
import re
import shutil
import sys
import urllib.error
import urllib.request
import zipfile
from datetime import datetime, timezone
from pathlib import Path

MSIGDB_RELEASE = "2026.1.Hs"
MSIGDB_BASE = "https://data.broadinstitute.org/gsea-msigdb/msigdb/release"
REACTOME_VERSION_URL = "https://reactome.org/ContentService/data/database/version"
REACTOME_GMT_URL = "https://download.reactome.org/{release}/ReactomePathways.gmt.zip"
REACTOME_GMT_MEMBER = "ReactomePathways.gmt"

# collection -> the prefix every set name in it must carry; a name outside it is refused
MSIGDB_COLLECTIONS = {
    "h.all": "HALLMARK_",
    "c2.cp.wikipathways": "WP_",
    "c2.cp.pid": "PID_",
    "c5.go.bp": "GOBP_",
    "c5.go.cc": "GOCC_",
    "c5.go.mf": "GOMF_",
    "c5.hpo": "HP_",
}
MSIGDB_SET_NAME = re.compile(r"[A-Z0-9_]+")
REACTOME_STABLE_ID = re.compile(r"R-HSA-[0-9]+")

# sources we may not redistribute, matched against a name with every non-alphanumeric removed so
# " KEGG_X", "KEGG-MEDICUS_X" and "KEGG X" all reduce to the same refused prefix
REFUSED_SOURCES = ("KEGG", "BIOCARTA")

MSIGDB_ATTRIBUTION = (
    "Gene sets from the Molecular Signatures Database (MSigDB), Broad Institute and "
    "UC San Diego, licensed CC BY 4.0. Cite Liberzon et al. 2015 Cell Systems "
    "(doi:10.1016/j.cels.2015.12.004) and Subramanian et al. 2005 PNAS "
    "(doi:10.1073/pnas.0506580102)."
)
HPO_LICENCE_URL = "https://human-phenotype-ontology.github.io/license.html"

# Keyed by (source, upstream release) and committed, never inferred from the downloaded files: a
# release with no entry here refuses to build, so a new upstream release reaches a bundle only
# through a reviewed edit of this table. The GO and HPO releases behind an MSigDB release are
# stated in that release's notes, which is what `release_notes` points at.
LICENCES = {
    ("msigdb", "2026.1.Hs"): {
        "spdx": "CC-BY-4.0",
        "attribution": MSIGDB_ATTRIBUTION,
        "release_notes": "https://docs.gsea-msigdb.org/MSigDB/Release_Notes/MSigDB_2026.1.Hs/",
        "go_release": "GO-basic obo file released on 2025-10-10 and NCBI gene2go annotations "
                      "downloaded on 2025-12-14",
        "hpo_release": "2026-01-08",
    },
    ("reactome", "97"): {
        "spdx": "CC0-1.0",
        "attribution": "Reactome pathway data, https://reactome.org/license (CC0 1.0). "
                       "Citation is requested, not required.",
    },
}


def sha256_of(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as fh:
        for block in iter(lambda: fh.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def licence_entry(source: str, release: str) -> dict:
    entry = LICENCES.get((source, release))
    if entry is None:
        known = sorted(r for s, r in LICENCES if s == source)
        raise SystemExit(
            f"no committed licence entry for {source} release {release!r} (have {known}); "
            f"review that release's terms and add it to LICENCES before building it"
        )
    return entry


def collection_licence(source: str, release: str, collection: str) -> tuple[str, str, dict]:
    """(SPDX id, attribution, extra manifest fields) for one collection of one release."""
    entry = licence_entry(source, release)
    if source != "msigdb":
        return entry["spdx"], entry["attribution"], {}
    notes = entry["release_notes"]
    if collection.startswith("c5.go."):
        attribution = (
            f"{entry['attribution']} Contains Gene Ontology data (Gene Ontology Consortium, "
            f"CC BY 4.0; Ashburner et al. 2000, doi:10.1038/75556; The Gene Ontology Consortium "
            f"2023, doi:10.1093/genetics/iyad031), GO release: {entry['go_release']}, as stated "
            f"in the MSigDB {release} release notes ({notes}), which give no DOI for that "
            f"GO release."
        )
        return entry["spdx"], attribution, {"go_release": entry["go_release"],
                                            "go_release_source": notes}
    if collection == "c5.hpo":
        attribution = (
            f"{entry['attribution']} Contains data from the Human Phenotype Ontology "
            f"(https://hpo.jax.org), HPO release {entry['hpo_release']} as stated in the MSigDB "
            f"{release} release notes ({notes}). The Human Phenotype Ontology Consortium is "
            f"acknowledged; cite Köhler et al. 2014 Nucleic Acids Res. 42(D1):D966-D974 "
            f"(doi:10.1093/nar/gkt1026). Under the HPO licence ({HPO_LICENCE_URL}) the HPO "
            f"must be cited, the HPO version displayed, and its content and logical relationships "
            f"not altered."
        )
        return "LicenseRef-HPO", attribution, {"hpo_release": entry["hpo_release"],
                                               "hpo_release_source": notes,
                                               "licence_url": HPO_LICENCE_URL}
    return entry["spdx"], entry["attribution"], {}


def refused_source(text: str) -> bool:
    return re.sub(r"[^A-Z0-9]", "", text.upper()).startswith(REFUSED_SOURCES)


def check_gmt(path: Path, collection: str) -> int:
    """Return the number of sets in a GMT, refusing any set from a source we may not ship.

    The collection list already avoids the bundles that contain these; this is the second guard,
    for the day an upstream changes a file's contents or a collection is added carelessly. The
    deny check alone can be dodged by a name spelt differently, so MSigDB names must also match
    their collection's exact naming scheme, and Reactome's free-text names must sit next to a
    Reactome stable id.
    """
    is_reactome = collection == "reactome"
    prefix = None if is_reactome else MSIGDB_COLLECTIONS[collection]
    n_sets = 0
    for lineno, line in enumerate(path.read_text(encoding="utf-8").split("\n"), 1):
        if line == "":
            continue
        fields = line.split("\t")
        where = f"{path}:{lineno}"
        if len(fields) < 3:
            raise SystemExit(f"{where}: not a GMT line (expected >= 3 tab-separated fields)")
        name, link = fields[0], fields[1]
        if not name.strip():
            raise SystemExit(f"{where}: empty set name")
        if refused_source(name) or (not is_reactome and refused_source_in_url(link)):
            raise SystemExit(
                f"{where}: set {name!r} ({link!r}) is from a refused source "
                f"{REFUSED_SOURCES}; its source cannot be redistributed"
            )
        if is_reactome:
            if not REACTOME_STABLE_ID.fullmatch(link):
                raise SystemExit(f"{where}: set {name!r} has {link!r}, not a Reactome stable id")
        elif not (MSIGDB_SET_NAME.fullmatch(name) and name.startswith(prefix)):
            raise SystemExit(f"{where}: set name {name!r} is not a {prefix}* {collection} name")
        n_sets += 1
    if n_sets == 0:
        raise SystemExit(f"{path}: no gene sets")
    return n_sets


def refused_source_in_url(url: str) -> bool:
    upper = url.upper()
    return any(source in upper for source in REFUSED_SOURCES)


def download(url: str, dest: Path, *, timeout: int = 300) -> Path:
    """Download `url` to `dest`, refusing a body shorter or longer than its Content-Length.

    Not sumstat_utils.fetch: http.client returns a short body without raising when the
    connection drops early, and that file would then be checksummed into the manifest as good.
    `dest` is inside the build's `.partial` directory, which a failure removes whole.
    """
    print(f"  downloading {url}", file=sys.stderr)
    request = urllib.request.Request(url, headers={"User-Agent": "genetics-results-munge"})
    try:
        with urllib.request.urlopen(request, timeout=timeout) as response, dest.open("wb") as out:
            expected = response.headers.get("Content-Length")
            shutil.copyfileobj(response, out)
            written = out.tell()
    except (urllib.error.URLError, TimeoutError) as err:
        raise SystemExit(f"could not fetch {url}: {err}") from err
    if expected is not None:
        if not expected.strip().isdigit():
            raise SystemExit(f"{url}: unparseable Content-Length {expected!r}")
        if written != int(expected):
            raise SystemExit(f"{url}: got {written} bytes, Content-Length says {expected}; "
                             f"truncated download")
    return dest


def require_trailing_newline(path: Path) -> None:
    """MSigDB GMTs end with a newline; one that does not was cut short (Reactome's may not)."""
    with path.open("rb") as fh:
        fh.seek(0, 2)
        if fh.tell() == 0 or (fh.seek(-1, 2), fh.read(1))[1] != b"\n":
            raise SystemExit(f"{path}: does not end with a newline; truncated download?")


def reactome_release() -> str:
    request = urllib.request.Request(REACTOME_VERSION_URL, headers={"User-Agent": "genetics-results-munge"})
    with urllib.request.urlopen(request, timeout=60) as response:
        text = response.read().decode().strip()
    if not re.fullmatch(r"[0-9]{1,4}", text):
        raise SystemExit(f"{REACTOME_VERSION_URL}: expected a release number, got {text[:80]!r}")
    return text


def extract_member(zip_path: Path, member: str, dest: Path, max_bytes: int = 100 << 20) -> None:
    """Copy the one expected member out of an untrusted zip; any other content is refused."""
    with zipfile.ZipFile(zip_path) as archive:
        names = archive.namelist()
        if names != [member]:
            raise SystemExit(f"{zip_path}: expected exactly [{member!r}], found {names[:10]}")
        info = archive.getinfo(member)
        if info.file_size > max_bytes:
            raise SystemExit(f"{zip_path}:{member} is {info.file_size} bytes, over the {max_bytes} cap")
        with archive.open(info) as src, dest.open("wb") as out:
            shutil.copyfileobj(src, out)


def collection_entry(*, name, source, url, release, downloaded, path, n_sets, extra=None) -> dict:
    spdx, attribution, licence_extra = collection_licence(source, release, name)
    return {
        "collection": name,
        "source": source,
        "url": url,
        "upstream_release": release,
        "download_date": downloaded,
        "file": path.name,
        "bytes": path.stat().st_size,
        "sha256": sha256_of(path),
        "n_sets": n_sets,
        "spdx_license": spdx,
        "attribution": attribution,
        **licence_extra,
        **(extra or {}),
    }


def write_manifest(path: Path, manifest: dict) -> None:
    path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def build(out_dir: Path, *, msigdb_release: str, reactome: str | None, today: datetime) -> Path:
    release_id = today.strftime("%Y%m%d")
    final = out_dir / f"gene_sets.{release_id}"
    if final.exists():
        raise SystemExit(f"{final} already exists; refusing to overwrite (versions are immutable)")
    work = out_dir / f"gene_sets.{release_id}.partial"
    if work.exists():
        raise SystemExit(f"{work} exists from an interrupted run; remove it by hand")
    # both releases are checked against the committed table before anything is downloaded
    licence_entry("msigdb", msigdb_release)
    reactome = reactome or reactome_release()
    licence_entry("reactome", reactome)
    work.mkdir(parents=True)
    downloaded = today.strftime("%Y-%m-%d")
    entries = []
    try:
        for name in MSIGDB_COLLECTIONS:
            filename = f"{name}.v{msigdb_release}.entrez.gmt"
            url = f"{MSIGDB_BASE}/{msigdb_release}/{filename}"
            path = download(url, work / filename)
            require_trailing_newline(path)
            entries.append(collection_entry(
                name=name, source="msigdb", url=url, release=msigdb_release,
                downloaded=downloaded, path=path, n_sets=check_gmt(path, name)))

        url = REACTOME_GMT_URL.format(release=reactome)
        zip_path = download(url, work / "ReactomePathways.gmt.zip")
        gmt = work / REACTOME_GMT_MEMBER
        extract_member(zip_path, REACTOME_GMT_MEMBER, gmt)
        entries.append(collection_entry(
            name="reactome", source="reactome", url=url, release=reactome,
            downloaded=downloaded, path=zip_path, n_sets=check_gmt(gmt, "reactome"),
            extra={"gmt_file": gmt.name, "gmt_sha256": sha256_of(gmt),
                   "gene_id_type": "HGNC symbol"}))

        write_manifest(work / "manifest.json", {
            "release_id": release_id, "build_date": downloaded, "collections": entries})
    except BaseException:
        # a refused build must not leave a directory that looks like a release
        shutil.rmtree(work, ignore_errors=True)
        raise
    work.rename(final)
    return final


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--msigdb-release", default=MSIGDB_RELEASE)
    parser.add_argument("--reactome-release", help="default: the current release, resolved online")
    args = parser.parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)
    final = build(args.out_dir, msigdb_release=args.msigdb_release, reactome=args.reactome_release,
                  today=datetime.now(timezone.utc))
    for entry in json.loads((final / "manifest.json").read_text())["collections"]:
        print(f"  {entry['collection']:20s} {entry['n_sets']:6d} sets {entry['bytes']:10d} bytes",
              file=sys.stderr)
    print(f"wrote {final}")


if __name__ == "__main__":
    main()
