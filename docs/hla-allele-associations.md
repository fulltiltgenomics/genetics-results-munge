# Classical HLA allele associations (FinnGen)

`scripts/munge_hla.{py,sh}` turns a FinnGen release's imputed classical HLA association
results into the two artifacts the suite serves. This note records the decisions that are
not obvious from the code.

## Source

The release `hla/` directory in the FinnGen green library, e.g.

```
gs://finngen-production-library-green/finngen_R14/finngen_R14_analysis_data/hla/
  summary_stats/<PHENO>.gz(.tbi)   per-phenotype REGENIE results, one row per HLA allele
  R14_imputed_hla.snpstats         qctool per-allele QC (AF, HWE, imputation INFO)
  r14_hla_regenie_summary.tsv.gz   top allele per phenotype (not used — derivable)
```

The whole directory is ~18 MiB, so the wrapper pulls it in one recursive copy rather than
reading 2,800 objects individually.

R10 through R14 all ship this directory. Only R14 is munged and served: R12/R13 HLA would
add a second copy of the same 507k rows for a release the suite no longer serves for
anything else.

## What the munge changes, and why

**1. An allele is not a variant.** The source encodes the dosage test as a biallelic
variant: `ref` is the literal string `<absent>` and `alt` carries the allele name
(`A*02:01`). That is faithful to the test actually run — presence of the allele versus its
absence — but downstream it reads as a nucleotide variant everywhere, and it leaves the
HLA gene implicit inside the allele string. The munge drops `ref`/`alt` and writes
explicit `gene` (`HLA-A`, derived from the locus before the `*`) and `allele` columns.

The consequence is deliberate and worth stating: these rows have no `chr:pos:ref:alt`
identity, so they do not join to credible sets or variant annotations, and no by-variant
API path can ever match one.

**2. Imputation INFO is joined onto every row.** The `.snpstats` sidecar carries each
allele's imputation quality, constant across phenotypes. Without it a returned association
is not interpretable: rare alleles imputed below ~0.5 produce enormous unstable betas that
look like spectacular findings. `B*07:01` (INFO 0.016, AF 0.0004) reaches |beta| > 40 in
several endpoints. One joined column makes every row self-describing about that, at the
cost of duplicating 187 values across 507k rows.

qctool writes a block of `#`-prefixed provenance lines before the header, so the sidecar
parser skips to the first non-comment line rather than reading line 1 as the header.

**3. Phenotypes without release metadata are dropped.** The source ships more endpoints
than the release phenotype metadata JSON describes — 2,808 versus 2,716 for R14. The 96
extras are all `_WIDE` variants of endpoints that are themselves present. Serving them
would surface traits in the API with no name, case count or category, so they are
filtered out. 2,712 phenotypes survive; 4 core endpoints (`AB1_EBV`,
`AB1_PNEUMOCYSTOSIS`, `AB1_SYPHILIS`, `C3_ASTROCYTOMA_WIDE`) have no HLA run at all and
are reported by the script.

**4. Case/control allele frequencies are optional.** Quantitative endpoints (BMI, height,
Kanta labs) have no `af_alt_cases`/`af_alt_controls` columns, exactly as in the core
FinnGen sumstats. They are filled with `NA` rather than treated as a schema violation.

## Outputs

Two artifacts, because the data has two query axes and no single file serves both.

| Artifact | Layout | Consumer |
|---|---|---|
| `summary_stats/<PHENO>.tsv.gz` (+ `.tbi`) | one file per phenotype, ~187 rows | results-api `/hla` — all alleles for a trait |
| `<dataset-id>.tsv.gz` | one combined file, ~507k rows, carries `phenotype` | BigQuery `hla_associations` — all traits for an allele |

Per-phenotype columns:

```
#chrom  pos  gene  allele  pval  mlogp  beta  sebeta  af_alt  af_alt_cases  af_alt_controls  info
```

The combined file inserts `phenotype` after `allele` and drops the leading `#`. Neither
carries `dataset`; the BigQuery loader injects it with `--const-column`.

Both are tabix-compatible point indexes (`-s1 -b2 -e2`) on the gene anchor position. Only
the per-phenotype files are actually indexed — BigQuery does not read an index.

**Positions are gene anchors, not allele locations.** Every allele of a gene is written at
one position, which is what the tabix index is built on. HLA-DRB3/DRB4/DRB5 share the
placeholder `32500000` (the pipeline gives the secondary DRB loci one synthetic anchor),
so they cannot be separated positionally — only by the `gene` column.

**`pval` underflows to 0** for the strongest signals (coeliac `DQB1*02:01` is mlogp 1596).
`mlogp` is the field to rank and threshold on, everywhere.

## Running it

```bash
# produce locally, no upload
scripts/munge_hla.sh

# produce and publish to both profile buckets
scripts/munge_hla.sh --stage
```

`SOURCE`, `PHENO_METADATA`, `DATASET_ID` and `WORKDIR` override the defaults; `--stage` is
never implied. Staged layout mirrors the results-api profile config:

```
gs://finngen-commons/results_api_data/hla/finngen_hla/    # finngen profile
gs://daly-genetics-results/hla/finngen_hla/               # daly profile
  finngen_hla.tsv.gz
  summary_stats/<PHENO>.tsv.gz(.tbi)
```

After staging, load BigQuery from `genetics-results-db`:

```bash
PROJECT_ID=<project> GCS_BUCKET=finngen-commons GCS_PREFIX=results_api_data/ scripts/load_hla.sh
```
