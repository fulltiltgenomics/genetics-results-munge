# Gene-disease associations (GenCC, Monarch Initiative)

`scripts/munge_gene_disease.{py,sh}` downloads the two curated gene-disease sources the
results-api serves from `/gene_disease/{gene}` and writes one versioned TSV per source.
This note records the decisions that are not obvious from the code.

## Sources

| product | URL | what it is |
|---|---|---|
| `gencc` | `https://search.thegencc.org/download/action/submissions-export-tsv` | one row per curation submission, 30 columns; ClinGen, Genomics England PanelApp, Orphanet, Labcorp Genetics, G2P and others, each with its own validity classification and mode of inheritance |
| `monarch` | `https://data.monarchinitiative.org/monarch-kg/latest/tsv/all_associations/causal_gene_to_disease_association.all.tsv.gz` | Biolink `causes` and `associated_with_increased_likelihood_of`, from OMIM and ClinGen |
| `monarch` | `https://data.monarchinitiative.org/monarch-kg/latest/tsv/gene_associations/gene_disease.noncausal.tsv.gz` | Biolink `gene_associated_with_condition` (Orphanet) and `contributes_to` (OMIM) |

Both Monarch files share the KG's 15-column association layout and are concatenated into
one output. `latest/metadata.yaml` names the release and versions that output.

Neither source carries genomic coordinates, so neither output is bgzipped or tabixed: the
API reads the plain TSV out of GCS once at startup and answers from memory.

## Versioning, and why the old file is never overwritten

The output name carries the source's version — the KG release for Monarch, the UTC
download date for GenCC, whose export publishes no release identifier of its own. Both
deployments read their file path from a results-api profile
(`app/config/profiles/<profile>/gene_disease.py`), so publishing under a new name and
editing that path are two separate steps, and the previous file stays readable until the
second one lands. Overwriting in place would cut both deployments over at upload time,
with nothing to roll back to.

## What the munge changes, and why

**GenCC is passed through.** The only transformation is a polars parse/serialize round
trip — byte-identical to the upstream export as of the release this was written against.
Its value is that it fails here rather than in the API when an export stops parsing, and
it asserts the 30 column names the API selects by. The export's quoting has already
changed once without notice (every field quoted, then only the fields that need it); both
parse the same, and a round trip that stops being byte-identical is not by itself a
problem.

**Monarch takes the non-causal file as well as the causal one.** The suite served the
causal file alone until 2026-09. Monarch's causal export contracted sharply between the
release the suite was serving and the current one — about a quarter of the gene-disease
pairs it used to carry are in no current gene-disease export, causal or not. Taking both
files grows the endpoint's coverage rather than shrinking it, and adds Orphanet, which
the dataset description had always claimed and the file had never actually contained.

The cost is that "Monarch says so" no longer means one thing: a `causes` row and a
`gene_associated_with_condition` row are different claims. That is why `predicate`
survives the munge and reaches the API, which maps it onto the `classification` column
that carries GenCC's validity terms — the endpoint's columns are a harmonization of
per-source vocabularies, not a shared one, exactly as `submitter` already holds both
`Ambry Genetics` and `infores:omim`.

`predicate` loses its `biolink:` prefix. Every value carries it, so it distinguishes
nothing, and it is the one column of the fifteen whose values are read by a person.

**Rows are dropped in three cases**, each of which is empty or nearly empty today and is
filtered anyway because the KG's shape is not the suite's to fix:

- `subject_category` is not `biolink:Gene`. The non-causal file carries a handful of
  disease-to-disease rows — an OMIM subtype `contributes_to` its parent. Left in, they
  reach the endpoint as a gene named `major depressive disorder 1`.
- `subject_taxon` is not human. Both files are human-only; the KG's other gene
  association files are not.
- `negated` is set.

**Exact duplicate rows are collapsed.** Both files repeat rows verbatim. The API builds
its Monarch `uuid` as `subject|object|primary_knowledge_source|predicate`, so a repeat
does not merely show twice — it shows twice under one identifier. The script asserts that
composite is unique before writing, which is the check that would catch Monarch starting
to distinguish two rows by a field the uuid does not include.

## Running it

```
# produce locally, no upload
PRODUCT=gencc   scripts/munge_gene_disease.sh
PRODUCT=monarch scripts/munge_gene_disease.sh

# produce and publish to both profile buckets
PRODUCT=monarch scripts/munge_gene_disease.sh --stage

# publish to one bucket only, when a refresh is being taken up by one deployment
PRODUCT=monarch GCS_DESTS="gs://daly-genetics-results/gene_disease/" \
  scripts/munge_gene_disease.sh --stage
```

Inputs are cached under `CACHE_DIR` and re-downloaded only when missing, so a rerun after
a schema surprise does not re-fetch 26 MB. Staging does not cut anything over: point the
results-api profile at the new file name, and update the dataset's `version` and
`description` in the suite's `configs/datasets.yaml`.
