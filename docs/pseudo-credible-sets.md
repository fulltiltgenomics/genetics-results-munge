# Pseudo credible sets

## What they are and why

Statistical fine-mapping (e.g. SuSiE) produces *credible sets* — small groups of
variants that, with high probability, contain the causal variant at a locus, each
member carrying a posterior inclusion probability (PIP). Many external GWAS we ingest
only ship summary statistics, with no fine-mapping output. For those datasets we build
**pseudo credible sets**: LD- and significance-based approximations of real credible
sets, with a heuristic PIP, so they can live in the same `credible_sets` table and be
browsed alongside genuinely fine-mapped results.

A pseudo CS is **not** a Bayesian fine-mapping result. It is an LD clump around a
genome-wide-significant lead variant, trimmed to the variants that plausibly tag the
same signal, with PIP assigned from the relative association strength of the members.
The `pip` and `cs_min_r2` columns should be read with that caveat in mind.

## The two-step flow

Pseudo credible sets are produced from munged summary statistics in two pipeline runs:

```
munged sumstats (.munged.tsv.gz)
        │
        ▼
[1] autoreporting WDL   (~/autoreporting, FinnGen autoreporting tool)
        │   LD clumping → per-variant report assigning each variant to a locus
        ▼
  <pheno>.report.out    (one row per variant, with locus_id and r2_to_lead)
        │
        ▼
[2] create_pseudo_credible_sets.wdl   (this repo, wdl/)
        │   trim each locus to pseudo-CS members, assign PIP, filter
        ▼
  <dataset>_pseudo_credible_sets.mlog10p_2.r2_0.6.tsv.gz  (+ .tbi)
```

The naming suffix `mlog10p_2.r2_0.6` records the two main step-2 thresholds
(`--mlog10p-diff 2`, `--r2-to-lead-thres 0.6`).

---

## Step 1 — LD clumping in autoreporting

Step 1 runs the FinnGen autoreporting tool (`~/autoreporting`). The example input used
for the external-sumstats pseudo-CS case is [`wdl/autoreporting_external.json`](../wdl/autoreporting_external.json).
Relevant settings from that input:

| autoreporting setting | value | meaning |
|---|---|---|
| `primary_grouping_method` | `ld` | group variants by LD clumping (not by SuSiE credible sets) |
| `sign_treshold` | `5e-8` | primary (p1) threshold — a variant must reach this to be a clump **lead** |
| `alt_sign_treshold` | `1e-2` | secondary (p2) threshold — a variant must reach this to be an LD **partner** |
| `locus_width_kb` | `2000` | half-window: only variants within ±2 Mb of the lead are considered |
| `ld_opts` | `--dynamic-r2-chisq 5.0` | dynamic per-lead r² threshold (see below) |
| `pval_is_mlog10p` | `true` | the `mlog10p` column holds −log10(p); thresholds are converted to that scale (5e-8 → 7.30, 1e-2 → 2.0) |
| `ld_panel` / `ld_api` | FinnGen R12 LD matrix, `tabix` | source of pairwise r² |
| `overlap` | `false` | each variant is assigned to at most one clump |

### Dynamic r² threshold

With `--dynamic-r2-chisq 5.0`, the r² cutoff for LD partners is **per-lead**, not fixed:

```
r2_threshold(lead) = min( 5.0 / chi2.isf(p_lead, df=1), 1.0 )
```

(`Scripts/grouping.py::ld_threshold`). The stronger the lead's signal (larger χ²), the
*lower* the r² needed to be pulled into the clump, so strong signals gather broader
clumps. A lead right at p = 5e-8 (χ² ≈ 29.7) needs partners at r² ≳ 0.17.

### What the threshold means — expected χ² of a partner

The LD reader keeps a partner only when **`r2 > ld_threshold`** (strict;
`data_access/linkage.py`). Substituting the dynamic threshold, a partner is included iff:

```
r2 > 5.0 / chi2_lead   ⟺   chi2_lead · r2 > 5.0   ⟺   E(χ²_partner) > 5
```

`chi2_lead · r2` is the chi-square a variant would be **expected** to show if its
association were entirely due to tagging the lead (under LD, χ² scales with r² to the
lead). So `--dynamic-r2-chisq 5.0` means: *group a nearby variant with the lead when its
expected χ² — explained purely by the lead — exceeds 5.* Consequences:

- A genome-wide-significant variant with `E(χ²) > 5` relative to a stronger nearby lead is
  **absorbed into that lead's clump**: it appears in the report as an LD-partner row
  (`locus_id` = the stronger lead, with `r2_to_lead` set), and `ld_grouping` removes it
  from the candidate-lead pool so it never seeds its own locus. It is therefore **not**
  reported as an independent secondary signal — but it is *not* dropped from the report.
- A variant with `E(χ²) ≤ 5` (r² too low to clear the threshold) is **left out of the
  clump**; if it independently reaches p1 = 5e-8 it survives as a candidate lead and forms
  its **own** locus — a genuine secondary signal.
- The cut is strict, so a variant at exactly `E(χ²) = 5` is **not** absorbed.

This per-peak rule is what distinguishes a true independent secondary signal from an LD
shadow of a stronger one, and it is the only place clump membership is decided — upstream
of, and separate from, the optional 1 Mb proximity filter in step 2.

### Clumping algorithm (`Scripts/grouping.py::ld_grouping`)

Per chromosome:

1. Load every variant passing **p2** (1e-2) into the partner pool; those also passing
   **p1** (5e-8) become candidate leads.
2. Pop the most significant remaining candidate lead.
3. Query the LD panel for variants within ±2 Mb at r² ≥ the lead's dynamic threshold;
   keep those that are also in the partner pool.
4. Remove the lead and all its partners from both pools (because `overlap = false`), so
   they cannot seed or join any later clump.
5. Repeat until no candidate leads remain.

Each clump becomes a *locus*. The per-variant report (`<pheno>.report.out`) then carries,
for every variant, the columns step 2 relies on:

- `#variant` — the variant id (`chr:pos:ref:alt`)
- `locus_id` — the **lead** variant id of the clump the variant belongs to
  (a variant is its own clump's lead when `#variant == locus_id`)
- `r2_to_lead` — pairwise r² between the variant and its clump lead
- `all_inv_var_meta_mlogp`, `all_inv_var_meta_beta`, `all_inv_var_meta_sebeta`,
  `fg_af_alt`, `most_severe_gene`, `most_severe_consequence`

(For non-FinnGen inputs whose columns are named differently, step 2 maps them via
`column_aliases` — e.g. `#chr`→`#CHR`, `mlog10p`→`all_inv_var_meta_mlogp`,
`se`→`all_inv_var_meta_sebeta`, `af`→`fg_af_alt`.)

---

## Step 2 — building pseudo credible sets

Step 2 is [`wdl/create_pseudo_credible_sets.wdl`](../wdl/create_pseudo_credible_sets.wdl)
(embedded `create_pseudo_credible_sets.py`). It scatters one task per `.report.out` file,
then merges, sorts, bgzips and tabix-indexes the results.

Each input locus (grouped by `locus_id`) becomes at most one pseudo credible set.

### Membership rule

For a locus with lead −log10(p) = `L`, a variant of that locus is **kept** as a pseudo-CS
member if it satisfies any of:

1. it is the lead (`#variant == locus_id`), **or**
2. `r2_to_lead > r2_high_ld_thres` (default **0.95**) — the high-LD exception: a variant
   in near-perfect LD with the lead is always kept regardless of its p-value, **or**
3. `|L − mlog10p| < mlog10p_diff` (**2** in every production config; the script's own default
   is 3) **and** `r2_to_lead > r2_to_lead_thres` (default **0.6**) — i.e. comparable association
   strength *and* sufficient LD with the lead.

(`build_pseudo_credible_sets`.) Variants with null `mlog10p` or null `r2_to_lead` are
never kept.

### PIP assignment (`calculate_pip`)

PIP is a heuristic, derived purely from association strength within the set, not from
fine-mapping:

- weight of each member = 10^mlog10p (= 1/p), computed in log-space to avoid overflow;
- weights are normalised so PIPs sum to `total_pip` (**0.99**);
- each PIP is clamped to a floor of `minimum_pip` (**0.01**), redistributing the remainder
  iteratively; a singleton set gets PIP = `total_pip`.

### Per-set annotations and filters

- `cs_size` — number of member variants.
- `cs_min_r2` — minimum pairwise r² among all members, computed by
  `compute_and_filter_cs_pairwise_r2`, which queries the LD panel per set (singletons get
  1.0). If `--min-r2` is supplied, sets whose minimum pairwise r² falls below it are
  dropped; the production configs do **not** set `--min-r2`, so this only annotates.
- HLA region: with `--filter-hla` (default on), among all pseudo CS whose lead falls in
  chr6:25–34 Mb, only the single most significant one is kept; the rest are dropped.

### Output columns

`format_output` emits the standard credible-set schema (see the [outputs section of the
README](../README.md#outputs)): `#dataset, data_type, trait, trait_original, cell_type,
chr, pos, ref, alt, mlog10p, beta, se, pip, cs_id, cs_size, cs_min_r2, aaf, most_severe,
gene_most_severe`. `cs_id` is the lead variant id; chr `X` is mapped to `23`. The
`collect_results` task merge-sorts the per-trait files, de-duplicates, bgzips, and
indexes with `tabix -s6 -b7 -e7`.

---

## What variants / sets are excluded

**Excluded already in step 1 (never reach the report):**
- variants not passing the secondary threshold p2 = 1e-2 (mlog10p ≤ 2);
- variants farther than ±2 Mb from any lead;
- variants assigned to another clump (one-clump-per-variant, `overlap = false`);
- variants below the lead's dynamic r² threshold relative to every lead.

**Excluded in step 2, within a kept locus:**
- variants with null `mlog10p` or null `r2_to_lead`;
- variants with `r2_to_lead ≤ 0.6` — **unless** `r2_to_lead > 0.95` (high-LD exception);
- variants with `0.6 < r2_to_lead ≤ 0.95` whose `|lead_mlog10p − mlog10p| ≥ 2`
  (associated too weakly relative to the lead).

**Whole loci / sets excluded in step 2:**
- loci whose lead row is missing or has null `mlog10p`;
- loci whose lead `mlog10p < --min-mlog10p`, when that flag is set (not set in the
  production configs, so the effective lead-significance floor is autoreporting's
  p1 = 5e-8);
- all but the most significant pseudo CS within the HLA region chr6:25–34 Mb
  (when `--filter-hla`, the default);
- sets with minimum pairwise r² `< --min-r2`, when that flag is set (not set in
  production).

### Optional proximity filter

`filter_lead_variants` runs **before** set construction (only when *not* given
`--no-proximity-filter`) and drops weaker secondary lead loci sitting next to a very
strong signal. All production configs pass `--no-proximity-filter`, so this step is
currently **disabled** for every dataset; it is documented here for completeness and in
case it is ever re-enabled.

It operates on the locus leads only (rows where `#variant == locus_id`, deduped). For
each lead it computes a "strength" statistic from the lead's own association:

```
chi2 = (beta / se)**2
T    = min(0.1, 5.0 / chi2)      # None if beta/se missing or se == 0
```

`T` is small for strong signals (large χ²). It is the same `5/χ²` quantity autoreporting
uses for its dynamic r² threshold, capped at 0.1 instead of 1.0, used here purely as a
strength gate. A lead is treated as a strong **anchor** when `T < t_threshold`. Two
parameters are **hardcoded** in the function (not exposed as CLI flags):

- `t_threshold = 0.02` — `T < 0.02` means `5/χ² < 0.02`, i.e. **χ² > 250**
  (|z| ≳ 15.8, roughly mlog10p ≳ 55), so only genuinely huge signals act as anchors;
- `proximity_bp = 1_000_000` — the suppression radius (±1 Mb).

The rule: **walking leads from most to least significant, each anchor (`T < 0.02`)
suppresses every other lead on the same chromosome within ±1 Mb that is strictly less
significant than it.** Suppressed loci are removed wholesale (all their variants), and a
suppressed lead can no longer act as an anchor itself.

Properties that follow from the ordering:

- only anchors suppress — an ordinary lead (`T ≥ 0.02`) removes nothing;
- an anchor never removes a *more*-significant lead (strict `mlog10p_other < mlog10p_anchor`
  guard), and already-suppressed leads are skipped;
- when two anchors are within 1 Mb, the stronger is processed first and suppresses the
  weaker; an anchor more than 1 Mb away survives and can still prune its own weaker neighbors;
- suppression is purely **distance + significance** based — it ignores LD between the two
  leads and whether they tag the same haplotype.

Rationale: near an extremely strong association the step-1 LD clumping can leave residual
sub-peaks that are not truly independent; this prunes them so they do not become spurious
pseudo credible sets. The caveat — and likely why production disables it — is that the cut
is LD-agnostic, so a genuinely independent secondary signal within 1 Mb of a blockbuster
lead would also be removed.

---

## Per-dataset invocations

Each dataset has an input JSON in [`wdl/`](../wdl/) named
`create_pseudo_credible_sets.<dataset>.json`. All current datasets share the same core
flags — `--no-proximity-filter --mlog10p-diff 2 --r2-to-lead-thres 0.6`,
`r2_high_ld_thres = 0.95`, `filter_hla = true` — and differ only in the input
file-of-filenames, the phenotype-name JSON, and (for non-FinnGen inputs) the
`column_aliases` map. `--dataset` is not in the input JSON: the fofn has two columns per
line, dataset name and report file path, so one run can mix several datasets. Outputs land under
`gs://finngen-commons/results_api_data/credible_sets/.../` per the matching
`*.cromwell_options.*.json`, and are loaded into BigQuery by
`genetics-results-db/scripts/load_pseudo.sh`.

| dataset(s) | notes |
|---|---|
| `FinnGen_R13`, `FinnGen_R13_UKBB(_labs)`, `FinnGen_R13_MVP_UKBB(_labs)` | FinnGen and meta-analysis pseudo CS; phenotype names from `finngen_r13_pheno_202509.json` |
| `COVID19_HGI` / `PGC` / `GP2` (external) | bundled external file (`ext`); inputs use `column_aliases` to map munged-sumstat column names to the canonical report columns |
