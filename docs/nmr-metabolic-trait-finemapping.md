# EstBB-UKBB NMR metabolic trait fine-mapping

## What this is

[Tambets et al. 2026](https://www.nature.com/articles/s41586-026-10532-5) meta-analysed 249
circulating Nightingale NMR metabolic traits in the Estonian Biobank and UK Biobank (619,372
individuals) and published the SuSiE credible sets on
[Zenodo](https://zenodo.org/records/18132538). Those credible sets are what this munge converts.

The distinction that matters for every number downstream: the **GWAS** is the meta-analysis, the
**fine-mapping** is not. To keep the LD in-sample, SuSiE was run on the UKBB_EUR subset alone —
413,897 individuals — at 3 Mb windows around the meta-analysis European lead variants with
MAF > 0.1 %, MHC excluded. So the dataset's sample size is 413,897 and not 619,372, and the
effect sizes describe UKBB Europeans.

| | dataset (API) | `dataset` column | resource |
|---|---|---|---|
| | `nmr_ukbb_est` | `nmr_ukbb_est` | `nmr_ukbb_est` |

Registry key, `dataset` column and resource carry the same string, so the name an agent sees in a
query result is the name the dataset catalogue lists. The resource is **not** `ukbb`, which holds
UKB-PPP and UKB Finucane, and the value deliberately does not start with `UKB`: the
`dataset_to_resource_rules` entry for it is an exact match, and a `UKB`-prefixed value would be
claimed by the `UKB%` rule.

## Running it

```
scripts/munge_nmr_meta.sh \
  UKBB_EUR_fine_mapping_with_meta_EUR_lead_variants.parquet \
  EUR_all_lead_variants.tsv.gz \
  R14_annotated_variants_v0.small.gz \
  data/nmr_meta
```

The three inputs are the Zenodo fine-mapping parquet, the companion record's published
lead-variant statistics ([10.5281/zenodo.18377015](https://zenodo.org/records/18377015) —
`EUR`, not `meta_EUR`, because that is the cohort that was fine-mapped) and a tabix-indexed
FinnGen variant annotation (`gs://<bucket>/variant_annotations/`). It writes the merged
bgzipped and tabix-indexed file, the 249 per-trait files under `individual/` and
`credible_set_stats.tsv`, in the layout the API expects.

The phenotype metadata is a separate script, because it is derived from a different source:

```
scripts/nmr_meta_phenotypes.py \
  --input UKBB_EUR_fine_mapping_with_meta_EUR_lead_variants.parquet \
  --output configs/nmr_ukbb_est_pheno.json
```

## How the Zenodo table is mapped onto the credible set schema

| output column | source |
|---|---|
| `chr`, `pos`, `ref`, `alt` | as published, GRCh38, X mapped to 23 |
| `trait`, `trait_original` | `molecular_trait_id`, the Nightingale biomarker code |
| `cell_type` | `plasma` — the NMR panel is measured in EDTA plasma |
| `pip` | `alpha{cs_index}`, the credible set's own SuSiE weight |
| `cs_id` | `<region>_<cs number>`, the same spelling the FinnGen pipeline writes |
| `cs_size` | members of the set, counted per (trait, cs_id) |
| `mlog10p` | published `LOG10P`, else from \|z\| via `log_ndtr`; **null where neither exists** |
| `beta` | published `BETA`, else `-z * se`; **null where neither exists** |
| `se` | published `SE`, else `s_trait / sqrt(2 f (1-f) (n + z²))`; **null where neither exists** |
| `cs_min_r2` | **not available**, always `NA` |
| `aaf` | published MAF, oriented by the annotation's non-Finnish European AF |
| `most_severe`, `gene_most_severe` | FinnGen variant annotation, matched on the variant id |

## The four things the source does not say, and how each was settled

### The fine-mapping file's z belongs to the reference allele

The fine-mapping file publishes a signed `z` and no effect allele. `ref` and `alt` are correctly
labelled against GRCh38 — but the sign is the **reference** allele's.

The companion lead-variant file settles it, because the two overlap on 49,815 (trait, variant)
pairs. They agree on the alleles (`ALL0` == `ref` and `ALL1` == `alt`, on every one) and on the
frequency (`MAF` identical to the last digit), and their z-scores are **exact negatives**: over
the 48,742 pairs with a finite z, `max |Z_lead + z_fm| = 0.0046`, and not a single pair shares a
sign. The lead file is the correct one — its `BETA` is positive for the alt allele at HMGCR
rs12916, which raises LDL-C, and agrees with the published direction at APOE ε4 and ε2, PCSK9
R46L, SORT1, LPL S447X and CETP as well.

So `beta` is `-z · se`. `munge_nmr_meta.py` re-asserts all three checks on every run rather than
trusting this paragraph, because it is the failure that would be invisible: the file loads,
indexes and queries perfectly with every direction of effect reversed.

### beta, se and mlog10p are published where they exist, derived where they do not

The fine-mapping file carries no effect size, standard error or p-value — only `z` and the MAF.
The companion record does, at every genome-wide significant lead variant, in the same cohort:
49,815 of the 3,792,183 rows (1.3 %) therefore carry the authors' own `BETA`, `SE` and `LOG10P`.

Everywhere else they are derived. For an inverse-normal transformed trait the standardised
effect follows from the z-score, the MAF and the sample size (Zhu et al. 2016):
`se = s_trait/sqrt(2f(1-f)(n + z²))`, `beta = -z·se`. `2f(1-f)` is symmetric in `f`, so the MAF's
unknown orientation does not enter.

`s_trait` is **fitted, not assumed**. Setting it to 1 — asserting the transformed trait has unit
residual variance after age, sex and 20 PCs — makes `|beta|` and `|se|` 8 % too large on the
median trait and 21 % too large on the worst. Fitting it per trait against that trait's published
standard errors removes that bias:

| | value |
|---|---|
| lead variants available per trait | 23 minimum, 192 median, 356 maximum |
| fitted scale across the 249 traits | 0.792 to 0.996, median 0.918 |
| within-trait spread of the scale | 1.6 % of its median, 4.9 % at worst |

`mlog10p` needs no scale: `z` is exact, so the derived value reproduces the published `LOG10P`
(146.537 for HMGCR rs12916, from both routes).

### 5,004 rows carry no statistics at all, and the alphas cannot rescue them

The fine-mapping file's `z` is a two-sided p-to-z conversion, and its largest finite value is
38.47 — precisely where a double-precision p-value underflows to zero. 6,077 of 3,792,183 rows
are therefore `±inf`. 1,073 of those are lead variants and take the published numbers; the
remaining **5,004 get a null `mlog10p`, `beta` and `se`**.

Nothing weaker would do. The beta formula **converges** as `|z| → ∞`, to
`s_trait/sqrt(2f(1-f))`, so extrapolating would have produced a finite, plausible, entirely
fictional effect size rather than failing.

SuSiE's own alphas look like a way out and are not, which is worth recording so nobody
re-derives it. `log(alpha)` is very nearly linear in `z²` within a credible set — median R² of
0.989 over 4,000 sets — so an overflowed `z` can in principle be inverted from the set's other
members. But of the 2,585 credible sets containing an overflowed row, only **11** have any member
with a finite `z` to anchor against, covering 120 rows; 1,753 of them are singletons. The
information is not in the file.

Those rows keep their `pip` and their credible set membership, and they include strong lipid
signals such as APOE. Ranking a metabolic locus by `mlog10p` therefore drops some of the top
hits; rank by `pip`.

### aaf takes its orientation, not its value, from the annotation

The source publishes MAF, which has lost the allele. The suite's column is the alternative
allele frequency, so something has to say which allele is the minor one. The FinnGen annotation
answers that — but its own `AF` is Finnish, and this is a UK cohort, so it supplies a bit and
not a number: `aaf` is `maf` where alt is the minor allele and `1 - maf` where it is not.

Which of the annotation's frequencies decides the bit was measured rather than assumed. Its
`*_enrichment_nfe` columns are AF_fin / AF_nfe, so dividing recovers a non-Finnish European
frequency. Against the published UKBB_EUR MAF over the 280,058 annotated variants:

| orienting on | median \|ΔMAF\| | disagrees > 0.1 | near-coin-toss variants |
|---|---|---|---|
| Finnish `AF` | 0.0285 | 4.7 % | 11,316 |
| recovered NFE AF | 0.0097 | 0.03 % | 4,276 |

so the NFE frequency decides, with the Finnish AF as the fallback where neither enrichment
column is usable. A wrong orientation leaves `maf` in `credible_sets_v` correct — it is
`LEAST(aaf, 1-aaf)` — but flips the risk/protective call in `credible_set_stats.py`.

## Phenotype names are matched, not typed

`scripts/nmr_meta_phenotypes.py` builds the 249-trait metadata JSON. It does **not** contain a
hand-written list of biomarker descriptions: the GWAS Catalog holds 2,241 study records for this
publication (249 traits × 9 ancestry groups), and those supply the strings the authors
registered. The mapping onto the codes is derived — a grammar for the compositional lipoprotein
codes (`[size_]class_measure[_pct]`), an explicit table for the rest — and then checked to be a
bijection onto the catalogue's 249 names. A wrong entry either steals a name another code needs
or leaves one unclaimed, and the script exits non-zero either way.

## What is dropped, and what is not

Nothing is dropped. All 3,792,183 credible set variants, all 123,899 credible sets (over 21,923
distinct `cs_id` values, which repeat across traits by design) and all 249 traits reach the
output. What is *missing* is per-column:

- 5,004 rows (0.13 %) have no `mlog10p`, `beta` or `se`;
- 14,957 of the 295,015 distinct variants (5.1 %) are outside the FinnGen imputation panel, so
  they have no `aaf`, `most_severe` or `gene_most_severe`;
- `cs_min_r2` is null everywhere.

244 rows repeat a variant within one credible set, which the source also does; `credible_sets_v`
documents that shape, so they are kept rather than deduplicated.
