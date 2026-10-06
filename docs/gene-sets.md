# Gene sets

`scripts/build_gene_sets.{py,sh}` is the download stage of the gene set bundle the sandbox's
enrichment analyses read. Harmonisation of gene ids and the parquet bundle are a later stage that
extends this document.

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

## Output

`<out-dir>/gene_sets.<YYYYMMDD>/` with the raw files and `manifest.json`. The directory is never
overwritten and appears only once every guard has passed. Licences are the committed `LICENCES`
table, keyed by source and upstream release and not read from the data: a release with no entry
refuses to build, so adopting a new MSigDB or Reactome release means reviewing its terms first.
The GO and HPO releases behind an MSigDB release come from that release's notes and go into the
GO and HPO collections' attribution. `c5.hpo` is `LicenseRef-HPO` (the HPO licence: cite, show
the version, do not alter), which only deployments listing it in `gene_set_licences` may carry.

Run: `OUT_DIR=/mnt/disks/data/gene_sets_build scripts/build_gene_sets.sh`
