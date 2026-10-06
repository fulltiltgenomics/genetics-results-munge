# Gene sets

`scripts/build_gene_sets.sh` builds the gene set bundle the sandbox's enrichment analyses read,
in two stages: `build_gene_sets.py` downloads and guards the raw collections,
`gene_set_bundle.py` harmonises them to Ensembl and writes the parquet bundle. Each script's
header docstring is the reference for its inputs and rules; this page says why.

## What is downloaded

Only the collections in `MSIGDB_COLLECTIONS` plus Reactome's own GMT
(`ReactomePathways.gmt.zip` from a release-pinned `download.reactome.org/<n>/` URL; the release is
read from the ContentService at run time and recorded). Per-subcollection MSigDB files are taken
because `c2.cp`, `c2.all` and the `msigdb.*` bundles contain KEGG and BioCarta sets, which cannot
be redistributed to every deployment. A second guard fails the build on any set whose name, with
punctuation and whitespace removed, starts with `KEGG` or `BIOCARTA` (or, for MSigDB, whose URL
names either); MSigDB names must also match their collection's own prefix, and Reactome names must
sit next to an `R-HSA-` stable id. A download shorter than its `Content-Length`, or an MSigDB GMT
without its final newline, is refused as truncated.

`<out-dir>/gene_sets.<YYYYMMDD>/` holds the raw files and `manifest.json`. The directory is never
overwritten and appears only once every guard has passed. Licences are the committed `LICENCES`
table, keyed by source and upstream release and not read from the data: a release with no entry
refuses to build, so adopting a new MSigDB or Reactome release means reviewing its terms first.
The GO and HPO releases behind an MSigDB release come from that release's notes and go into the
GO and HPO collections' attribution. `c5.hpo` is `LicenseRef-HPO` (the HPO licence: cite, show
the version, do not alter), which only deployments listing it in `gene_set_licences` may carry.

## The bundle

Members are mapped to unversioned Ensembl gene ids through the same staged
`hgnc_complete_set.txt` and GENCODE genes file that `gene_annotations_v` is built from, so a gene
set and a gene annotation agree on what a gene is. MSigDB maps by Entrez id, Reactome by symbol:
approved, then previous, then alias symbols, decided at the first tier that names the symbol at
all (over every HGNC gene, with or without an Ensembl id) and never falling through to a later
one. A fallback that skipped a tier because its gene lacked an Ensembl id would hand the symbol
to an unrelated gene's alias. What does not map is dropped and counted
per set (`n_unmapped`); a collection below its source's `MIN_MAPPING_RATE` refuses the bundle.
Reactome's floor is lower because its GMT carries viral, bacterial and rRNA members that are not
human genes at all.

`gene_sets_bundle.<raw release>.hgnc<YYYYMMDD>.gencode<N>/` holds `sets.parquet`, `members.parquet`,
`genes.parquet` and `manifest.json` (the raw manifest's per-collection licences unchanged, plus a
`bundle` block with the HGNC/GENCODE versions and checksums, per-collection mapping rates and a
sha256 per file). The tables are flat and long so they load into BigQuery as they are. The raw
GMTs are not in the bundle; the sandbox needs only the parquet and the manifest.

`--stage` copies the bundle to `gs://daly-genetics-results/gene_sets/<bundle release>/`, refusing
a remote path that already has objects. Each object is a separate create-only copy
(`--if-generation-match=0`), the parquet files first and `manifest.json` last, so an interrupted
stage leaves a path without `manifest.json`; the next stage refuses it as a partial upload to be
deleted by hand. Consumers must therefore treat a path without `manifest.json` as absent and
check each file's sha256 against it before use. Daly only: the bundle includes `c5.hpo`, which the finngen
profile does not allow.

Run (inputs copied locally from `gs://daly-genetics-results/mapping_files/`):

    OUT_DIR=/mnt/disks/data/gene_sets_build HGNC_FILE=.../hgnc_complete_set.txt \
      HGNC_VERSION=2026-06-02 GENCODE_FILE=.../gencode.v49.annotation.genes.tsv \
      scripts/build_gene_sets.sh [--stage [--dry-run]]

`RAW_DIR=<existing gene_sets.<YYYYMMDD>/>` skips the download.
