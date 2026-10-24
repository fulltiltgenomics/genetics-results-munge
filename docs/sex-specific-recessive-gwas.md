# Sex-specific and recessive GWAS (FinnGen R12)

FinnGen ran two model-specific GWAS over the whole endpoint set once, in release 12:
every endpoint in females and males separately with a female-male meta-analysis, and
every endpoint under a recessive genotype model. R13 and R14 ship neither, so these are
R12 products served beside the R14 core. This note records what the suite does with them
and the decisions that are not obvious from the code. `scripts/munge_sex_specific.{py,sh}`
is the only munge, and it writes one phenotype metadata file; the recessive release needs
none.

## Source

```
gs://finngen-production-library-green/finngen_R12/finngen_R12_analysis_data/
  sex_difference/
    meta_analysis/meta_summary_stats/<PHENO>_meta_out_filtered.tsv.gz(.tbi)   2,207 endpoints
    meta_analysis/sex_specific_meta_gwsig_in_either_sex_all_endpoints.tsv.gz   the <= 5e-8 tail
    summary_stats/<PHENO>.{females,males}.gz(.tbi)                             per-sex REGENIE
    R12_sex_specific_N_cases_controls_samples.tsv                               per-sex counts
    regenie_single_sex_drop_list_plus_reason_262_endpoints.txt                  the exclusions
  summary_stats_recessive/
    results/<PHENO>.gz(.tbi)                                                    2,469 endpoints
    data/<PHENO>.add_rec_combined_full.tsv.gz                                   superseded
```

The release notes (`finngen_R12_analysis_documentation/finngen_R12_sex_specific_meta.md`
and `finngen_R12_recessive_gwas.md`) describe the methods. Two things they get wrong about
their own delivery: the per-endpoint meta files are named `_meta_out_filtered.tsv.gz`, not
`.regenie.sex_diff.gz`; and the per-sex p-values do underflow to 0 in the aggregate.

## What is served, and from where

**Per-endpoint summary statistics are read in place, on request.** results-api serves the
sex-specific meta files and the recessive `results/` files straight from the green
library, by variant and by region, as two on-request resources of their own,
`finngen_sex_specific` (`data_type` `gwas_sex_specific`) and `finngen_recessive`
(`gwas_recessive`), in the finngen deployment only. On request means what it means for
the sandbox custom GWAS: absent from the default dataset catalogue and phenotype search,
reached when a request names the resource, and the agent is told to use them only when
the user asks for a sex-specific or recessive result — a disease or variant question is
answered from the core release. Nothing is copied: together the two products are ~7 TiB,
and a tabix index here is ~2 MiB, smaller than an R14 core index. The per-sex REGENIE
files are not registered — the meta file carries both sexes' statistics and the
heterogeneity test, so they add nothing a query needs.

**Nothing goes to BigQuery.** The sex-specific release ships one cross-endpoint
artifact, the genome-wide significant tail of every endpoint (`meta_analysis/
sex_specific_meta_gwsig_in_either_sex_all_endpoints.tsv.gz`, 2.12M rows over 1,388
endpoints), and it is not loaded: SQL has no resource gate, so a table of it would answer a
locus question from the sandbox with nobody having asked for the sex-specific model, and
the on-request rule could reach it only as prose. The dataset has no `datasets_v` or
`phenotypes_v` row either. A cross-endpoint sex-difference question is answered from the
per-endpoint files.

**The recessive release's lead tables are not carried.** The clumped recessive-beats-additive
tables exist only in the R12 working bucket (`gs://r12-data-green/regenie/add_rec_difference/`),
not in the library, so by the decision to load only what the library holds they are left
out. A user who wants recessive hits across endpoints has the per-endpoint files and
`mlog10p_diff`.

## Decisions in the munge

**1. The combined counts are what the phenotype metadata records.** `num_cases` and
`num_controls` sum the sexes, which is the meta-analysis sample; the three quantitative
endpoints carry their sample as `num_cases` with `num_controls` 0, the shape the R12 core
file uses. Only the 2,207 analysed endpoints are listed, so the 262 single-sex endpoints
cannot surface as traits with no data. The recessive dataset reuses the R12 core metadata
unchanged: its 2,469 endpoints are the core set, and the core file's 33 extra codes have no
R12 core summary statistics either.

**2. The per-sex p-values are left as the release wrote them.** In the per-endpoint files
`female_pval`/`male_pval` underflow to 0 on the strongest signals (APOE for `AD_LO_EXMORE`
among them, where the combined mlog10p reads 665); results-api serves the file's own
columns, so a reader ranking within a sex should use the combined `mlog10p` or recompute
from beta/se. Nothing is rewritten because nothing is copied.

## Staging

```
scripts/munge_sex_specific.sh --stage
```

stages `gs://finngen-commons/results_api_data/mapping_files/finngen_r12_sex_specific_pheno.json`,
which results-api reads for `/resource_metadata/finngen_sex_specific` and the named search.
finngen only: the daly profile has no registry entry, and nothing is staged in the daly
bucket.
