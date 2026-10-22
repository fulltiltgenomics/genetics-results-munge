# FinnGen R14 multi-cohort meta (FinnGen + All of Us + UKBB + MVP)

## Shape of the release

`gs://finngen-production-library-green/finngen_R14/finngen_R14_analysis_data/meta_analysis/`
holds one tabix-indexed `summary_stats/<phenocode>_meta_out.tsv.gz` per phenotype, a
phenotype mapping TSV, and the release's own autoreporting outputs under `autoreporting/`.
Each sumstats file carries the inverse-variance meta (`all_inv_var_meta_*`) over whichever
cohorts contributed that phenotype with its sample size `all_meta_N`, a
`<cohort>_{beta,sebeta,pval,af_alt}` block per contributing ancestry-split cohort (an MVP
block also carries `_r2`, its imputation quality), and `leave_cohort_<family>_*`
leave-one-out blocks. The contributing cohorts therefore vary per phenotype; results-api's
`app/config/profiles/r14_meta_columns.py` is the column mapping that serves them.

## How it is served

- **One resource, one dataset**: resource `finngen_mvp_ukbb_aou`, a single dataset over every
  phenocode in the mapping TSV (1188 on 2026-10-12). WIDE endpoints and Labs are ordinary
  phenocodes in that dataset — there is no separate data type for them.
- **Sumstats** are served straight from the green `summary_stats/` prefix; nothing is munged
  or copied.
- **Pseudo credible sets** come from the green `autoreporting/` reports through step 2 of
  [the pseudo-CS pipeline](pseudo-credible-sets.md) only. Those reports contain FinnGen
  variants only, so the R14 meta pseudo CS carry the same restriction the R13 meta pseudo
  CS do: a locus driven by a variant FinnGen lacks has no pseudo CS.

## Phenotype metadata

[`scripts/r14_meta_phenotypes.py`](../scripts/r14_meta_phenotypes.py) turns the mapping TSV
into the pheweb-shaped phenotype JSON the API serves as the dataset's metadata file, staged
at `gs://finngen-commons/results_api_data/mapping_files/finngen_r14_mvp_ukbb_aou_meta.json`.
Its docstring lists the input columns and output keys. Two choices are not obvious from
the code:

- **`num_controls` is omitted, not written as `NA`, for quantitative traits** (Labs and
  Quantitative endpoints: 85 rows on 2026-10-12), and so is a cohort's `n_controls` inside
  `n_by_cohort`; `num_cases` is then the sample size. This is the contract the JSON is written to for one reader: a dedicated `finngen_meta`
  harmonizer, to be added to genetics-results-api's `metadata_harmonizer.py` and mirrored in
  genetics-results-db's `build_phenotypes.py` (neither exists yet). It treats an absent
  `num_controls` as quantitative, reads the label from `name`, and passes
  `contributing_cohorts`, `cohort_families` and `n_by_cohort` through. No existing harmonizer
  may be pointed at this file: `finngen_r13` types every row binary (and defaults
  `num_controls` to 0), and `pheweb` reads only `phenostring` and types `num_cases`-only rows
  binary.
- **The per-phenotype field is `contributing_cohorts`, not `cohorts`.** The registry already
  has a dataset-level `cohorts` field holding collection tokens (finngen/ukb/mvp/aou) with a
  different meaning; `contributing_cohorts` is the ancestry-split list (`fg`, `aou_eur`, …,
  `mvp_amr`) for one phenotype. `cohort_families` is the fg/aou/mvp/ukbb roll-up of it.

## Cohort presence comes from the sumstats headers

The mapping TSV is not authoritative for which cohorts contributed: the script reads every
phenotype's sumstats header and keeps a cohort only if its `<cohort>_beta` column is there,
taking names, categories and n from the TSV. Measured 2026-10-12 over all 1188 files:
**8 disagree**, all Labs, all one direction — the TSV lists an MVP non-European cohort with
non-NA n that the sumstats do not carry; no header has a cohort the TSV omits.

| phenocode | in the TSV, not in the sumstats (dropped) |
|---|---|
| median_3005785_age_adjusted_IRN | mvp_amr |
| median_3011960_age_adjusted_IRN | mvp_amr |
| median_3013707_age_adjusted_IRN | mvp_amr |
| median_3018251_age_adjusted_IRN | mvp_afr |
| median_3020460_age_adjusted_IRN | mvp_amr |
| median_3021337_age_adjusted_IRN | mvp_amr |
| median_3021614_age_adjusted_IRN | mvp_afr, mvp_amr |
| median_3029187_age_adjusted_IRN | mvp_amr |

The TSV's `num_cases` equals the sum of its per-cohort `<cohort>_n_cases` on all 1188 rows,
and `num_controls` the sum of `<cohort>_n_controls` on all 1103 binary rows, including these 8.
So for a phenotype that loses a cohort the script recomputes the totals from the cohorts kept,
rather than serving a total that counts samples the meta did not use.
