# PGC schizophrenia published fine-mapping

## What this is

[Trubetskoy et al. 2022](https://www.nature.com/articles/s41586-022-04434-5) published FINEMAP
95 % credible sets for the PGC schizophrenia GWAS as supplementary table **ST11a**. This is
genuine fine-mapping, unlike the [pseudo credible sets](pseudo-credible-sets.md) we derive from
the wave 3 summary statistics — and both are served, side by side, under the same `pgc` resource:

| | dataset (API) | `dataset` column | source |
|---|---|---|---|
| pseudo | `pgc_scz` | `PGC` | LD clumping of the wave 3 daner file |
| fine-mapped | `pgc_scz_finemap` | `PGC_SCZ_2022` | ST11a |

They share the trait code `SCZ`, so a query that wants only genuine fine-mapping has to filter on
`dataset`, not on `resource`.

## Running it

```
scripts/munge_pgc_scz_finemap.sh \
  ST11a_95_perc_Credible_Sets.tsv \
  daner_PGC_SCZ_w3_90_0418b.munged.tsv.gz \
  R14_annotated_variants_v0.small.gz \
  data/pgc_scz_finemap
```

The three inputs are the supplementary table, the munged wave 3 summary statistics
(`gs://<bucket>/results_api_data/sumstats/PGC/`, see [`run_sumstats.sh`](../run_sumstats.sh)) and
a tabix-indexed FinnGen variant annotation
(`gs://<bucket>/results_api_data/variant_annotations/`). It writes the merged bgzipped and
tabix-indexed file, the per-trait file under `individual/` and `credible_set_stats.tsv`, in the
layout the API expects.

## How ST11a is mapped onto the credible set schema

| output column | source |
|---|---|
| `chr`, `pos`, `ref`, `alt` | munged wave 3 summary statistics, matched on `rsid` |
| `mlog10p` | `-log10(pval)` from ST11a |
| `beta` | `log(or)` from ST11a, negated where the effect allele is the reference allele |
| `se` | ST11a `se` |
| `pip` | ST11a `finemap_posterior_probability` |
| `cs_id` | `chr<chr>_<pos>_<ref>_<alt>` of the set's index SNP |
| `cs_size` | members surviving the rsid match |
| `cs_min_r2` | **not available**, always `NA` |
| `aaf` | wave 3 alt allele frequency; `NA` on chromosome X |
| `most_severe`, `gene_most_severe` | FinnGen variant annotation, matched on the variant id |

ST11a is on **GRCh37** and carries no reference/alternative alleles, only an effect and an other
allele. Rather than a liftover chain, the GRCh38 locus comes from the munged wave 3 summary
statistics by rsid — the same variant space the pseudo credible sets live in, which keeps the two
directly comparable, and which brings the alleles and the allele frequency with it. The match is
clean: no rsid lands on a different chromosome than ST11a reports, and no rsid has an
allele-incompatible locus.

Effect sizes, standard errors and p-values are kept from ST11a rather than taken from the summary
statistics, because most sets were fine-mapped on the *extended* GWAS and only ST11a's statistics
correspond to the reported posterior probabilities.

## Caveats

**FINEMAP was run with several causal variants allowed per locus.** ST11a carries an
`expected_causals_k` per set (up to 5), and posterior probabilities within one set sum to roughly
*k* rather than to 1. One `cs_id` can therefore span several independent signals, and sets are
much larger than SuSiE per-signal sets — 1 to 1742 members, median 33. Neither `expected_causals_k`
nor the `extended_gwas` flag has a column in the credible set schema, so both are lost here.

**The published odds ratios have two decimals.** `beta` is quantised in steps of roughly 0.005,
and rounds to exactly 0 for the weakest members. Against the wave 3 summary statistics the betas
still correlate at 0.9985 and the p-values at 0.9999.

**Chromosome X has no allele frequency.** The wave 3 daner file is autosomes only, so the five
X-linked sets are located through the FinnGen variant annotation instead. Those rows keep a null
`aaf`: the Finnish allele frequency is a different quantity from the case/control weighted one the
rest of the column holds, and mixing them silently in one column would be worse than a null. It
also means `maf` is null for them in `credible_sets_v`, and `credible_set_stats.py` classifies
them risk/protective by the sign of beta alone.

**`cs_min_r2` is unavailable.** ST11a reports no within-set LD and the script does not query an LD
panel, so the column is `NA` for every row of this dataset.

## What is dropped

Of the 20766 ST11a rows, **20508 (98.8 %) survive** in all **255** sets:

- 20381 are matched in the wave 3 summary statistics;
- 127 more — chromosome X, plus a handful of autosomal ones — are recovered from the variant
  annotation;
- 258 have an rsid in neither and are dropped, which trims 36 sets (worst: `rs2532240`, 190 of its
  1742 members);
- 20205 of the surviving 20508 get a `most_severe` / `gene_most_severe` annotation; the rest are
  outside the FinnGen imputation panel.

`cs_size` counts the members that survive, not the published set size.

For 11 of the 255 sets the index SNP is not itself a member of its own credible set — ST11a allows
this — so those sets are named after their highest-PIP member instead.
