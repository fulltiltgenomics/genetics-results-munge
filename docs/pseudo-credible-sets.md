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
for the external-sumstats pseudo-CS case is [`wdl/autoreporting_external.json`](../wdl/autoreporting_external.json);
[`wdl/autoreporting_aih.json`](../wdl/autoreporting_aih.json) is the same input with the AIH
file-of-filenames and phenotype info, and with `extra_columns` cut down to the columns
`munge_aih.py` actually writes (`se af rsid n`) — an extra column that is missing from the
sumstat is only warned about, not fatal, but listing only real ones keeps the report honest.
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
`se`→`all_inv_var_meta_sebeta`, `af`→`fg_af_alt`. `fg_af_alt` is the one column step 2
tolerates missing: a sumstat with no allele frequency, such as deCODE, yields a report
without `af`, and the output's `aaf` is then NA.)

A phenotype with no locus at all gets an empty (0-byte) `.report.out` from the WDL; step 2
writes a header-only file for it, so such phenotypes may stay in the file-of-filenames.

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
  dropped; the production configs do **not** set `--min-r2`, so this only annotates. The
  query maps chromosome `23` back to `X`, the name the panel's file and contig carry; before
  that mapping every chrX set silently got NA here (44 rows of `EXT_20260610`).
- HLA region: with `--filter-hla` (default on), among all pseudo CS whose lead falls in
  chr6:25–34 Mb, only the single most significant one is kept; the rest are dropped.

### Output columns

`format_output` emits the standard credible-set schema (see the [outputs section of the
README](../README.md#outputs)): `#dataset, data_type, trait, trait_original, cell_type,
chr, pos, ref, alt, mlog10p, beta, se, pip, cs_id, cs_size, cs_min_r2, aaf, most_severe,
gene_most_severe`. `cs_id` is the lead variant id; chr `X` is mapped to `23`.
`trait_original` is the report's file name with `.report.out` removed — autoreporting
names the report after the phenotype id it was given, and that id may itself contain dots
(SomaScan aptamers, `seq.10000.28`) — and `trait` is that id's `phenostring` in the
`phenotype_json`, or the id itself without one. `data_type` and `cell_type` are `GWAS` and
`NA` unless `--data-type`/`--cell-type` are passed in `flags`. The `collect_results` task
merge-sorts the per-trait files, de-duplicates, bgzips, and indexes with
`tabix -s6 -b7 -e7`.

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
strong signal. Every GWAS config passes `--no-proximity-filter`, so this step is
**disabled** for those datasets; the deCODE pQTL config is the one that leaves it on (see
the deCODE section for why).

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
flags — `--mlog10p-diff 2 --r2-to-lead-thres 0.6`, `r2_high_ld_thres = 0.95`,
`filter_hla = true` — and differ in the input file-of-filenames, the phenotype-name JSON,
(for non-FinnGen inputs) the `column_aliases` map, and whether `--no-proximity-filter` is
passed: every GWAS config passes it, the deCODE pQTL config does not. `--dataset` is not in the input JSON: the fofn has two columns per
line, dataset name and report file path, so one run can mix several datasets. Outputs land under
`gs://finngen-commons/results_api_data/credible_sets/.../` per the matching
`*.cromwell_options.*.json`, and are loaded into BigQuery by
`genetics-results-db/scripts/load_pseudo.sh`.

| dataset(s) | notes |
|---|---|
| `FinnGen_R13`, `FinnGen_R13_UKBB(_labs)`, `FinnGen_R13_MVP_UKBB(_labs)` | FinnGen and meta-analysis pseudo CS; phenotype names from `finngen_r13_pheno_202509.json` |
| `COVID19_HGI` / `PGC` / `GP2` (external) | bundled external file (`ext`); inputs use `column_aliases` to map munged-sumstat column names to the canonical report columns |
| `AIH` | the three autoimmune hypothyroidism meta-analysis phenotypes (`AIH`, `AITT1`, `AITT2`) munged by `scripts/munge_aih.py`; run separately from the `ext` bundle because its sumstats carry a different column set, but with the same `column_aliases` — autoreporting names the report columns after the input columns, which are the same `#chr`/`mlog10p`/`se`/`af` |
| `deCODE_pQTL_2021` | 4,907 SomaScan aptamers, one phenotype each, from `scripts/munge_decode_pqtl.py`; `--data-type pQTL --cell-type plasma`, `trait` = gene symbol via `configs/decode_pqtl_pheno.json`. See [the deCODE pQTL run](#the-decode-pqtl-run) |

---

## The deCODE pQTL run

The deCODE 2021 plasma pQTLs (Ferkingstad et al., 35,559 Icelanders) are the one pseudo-CS
input that is a QTL study rather than a GWAS, and the one that ships as a single file:
`gs://finngen-commons/decode/deCODE_pQTLs_NatGen2021_aligned_p0.005.tsv.gz`, every
aptamer's p < 0.005 rows together, already GRCh38 and aligned to gnomAD `ref`/`alt`. The
pipeline above is per phenotype, so the run has a step 0.

**Step 0 — one sumstat per aptamer.** `scripts/munge_decode_pqtl.py` streams the file once
into `<aptamer>.munged.tsv.gz` (`#chr pos ref alt mlog10p beta se`, X as 23), keeping one
row per variant where the alignment had folded two indel representations onto the same
gnomAD variant, and `--stage`s the directory to
`gs://finngen-commons/results_api_data/sumstats/deCODE_pQTL_2021/`. Its `--input-array`
writes the autoreporting input array with only the aptamers that have at least one row at
`sign_treshold` (p ≤ 5e-8): an aptamer with none cannot seed a locus, so a shard for it
would only produce an empty report. `scripts/decode_pqtl_phenotypes.py` writes the
phenotype-info TSV and the phenotype JSON from the aptamer → gene table the FinnGen
SomaScan credible sets use, so an aptamer carries the same `trait` in both datasets (the
seven aptamers that table maps to `NA` keep their id).

**Step 1** is [`wdl/autoreporting_decode.json`](../wdl/autoreporting_decode.json), the AIH
input with these differences:

| setting | value | why |
|---|---|---|
| `extra_columns` | `se` | the delivery has no allele frequency (its source carries only `maf`, which cannot be oriented to `alt`), no rsid and no per-variant `n` worth carrying |
| `post_process_top_reports.af_col`, `in_fg_col` | `lead_af_alt` | `meta_filter_top.py` raises if its AF column is absent from the top report and the WDL runs it on every shard that had results; `lead_af_alt` is one of the columns the top report always carries (it is NA here), so the top-report post-processing runs and filters nothing. Step 2 reads only `.report.out`, never the top report |
| `input_array_file`, `phenotype_info_file` | `external_sumstats_input.decode.tsv`, `external_pheno_info.decode.tsv` | from step 0 |
| `ld_assume_variant1_indexed` | `true` | autoreporting fetches LD for **every** candidate lead (every p ≤ 5e-8 variant) before clumping, and by default reads the whole ±2 Mb window of the panel for each one. A GWAS has tens to hundreds of candidates per phenotype; a cis-pQTL has thousands (median 929 per aptamer, p99 21,565, 11.5M over all aptamers), and a wide fetch costs 30 s to 5 min each, so the extreme aptamers would run for days. The FinnGen R12 panel is indexed by `variant1` position and lists every partner under it, so the 1 bp fetch this flag enables returns the same rows in well under a second (checked on the chr21:46003475 lead: 1,945 partner rows both ways). The flag is an input the WDL did not have; see the autoreporting change |

Everything else — the FinnGen R12 LD panel, `finngen_variants_only`, the thresholds — is
as for the other external datasets, with the same consequence: an Icelandic study is
clumped with Finnish LD, so leads whose Icelandic LD partners are absent or unlinked in
FinnGen become singleton sets (see the LD-panel discussion for the external GWAS). The
reports do not stay in the Cromwell bucket: step 1 took several submissions (below), so
`/mnt/disks/data/decode/stage_decode_reports.py <run id> ...` gathers them. For every
aptamer with a `.report.out` in any run listed (a later run wins), it copies the report and
its `.top.out` to
`gs://finngen-commons/results_api_data/sumstats/autoreporting/deCODE_pQTL_2021/`, writes
the step-2 fofn `gs://finngen-commons/results_api_metadata/autoreporting.decode.fofn`
from that location, and stages `external_sumstats_input.decode.rerun3.tsv`, the array of
the aptamers still without a report, which
[`wdl/autoreporting_decode.rerun3.json`](../wdl/autoreporting_decode.rerun3.json) submits.
A 0-byte report is a shard that ran and found no locus, and is kept.

**What the first submission taught.** Run `2f2fafa8` (2026-09-21) launched all 4,844
shards at once and lost most of them to one mechanism: within the first ten minutes,
about a hundred shards failed opening a `gs://` file from pysam — the LD panel, gnomAD or
the annotation — with `Invalid argument`, which htslib returns for an HTTP 4xx other than
401/403/404 (a wrong token gives `Operation not permitted`; checked against the image), so
almost certainly 429 from thousands of shards opening the same 23 panel files at the same
moment. `main.py` retries a failing *fetch* without bound but opens the panel exactly
once, so the shard exits 1. Cromwell's default failure mode is `NoNewCalls`: from the
first exit 1 it issued no more attempts, so every shard preempted after that point stayed
dead, and 4,210 of the 4,844 ended without a return code. Only 524 produced a report. Two
things are set differently for the rerun:

- [`wdl/autoreporting.cromwell_options.decode.json`](../wdl/autoreporting.cromwell_options.decode.json)
  sets `workflow_failure_mode` to `ContinueWhilePossible`, and the `report` task runtime in
  the autoreporting WDL carries `maxRetries: 2`, so a 4xx at startup costs one retry rather
  than the run. Both are needed together: the second run had `maxRetries` but the default
  failure mode, and only 5 of its 546 failed shards were retried — the ones that failed
  before the workflow entered its failing state;
- the remaining aptamers are resubmitted from a new input array that excludes every
  aptamer with a report in any earlier run (the staging script above builds it).

**What the second submission taught.** Run `965f9ac8` reran the 4,320 leftovers with the
narrow fetch. Shards now finished in minutes rather than hours, and the startup 4xx
recurred on only 23 shards, but 546 of the first 2,122 to finish were OOM-killed at 4 GB
during the LD prefetch. That image's `ld_grouping` fetched the partners of **every**
candidate lead at threshold 0 and held all of them until the greedy loop had consumed the
pile, so memory grew with the candidate count: `monitoring.log` showed about 2.3 GB before
the prefetch and roughly 0.9 GB more per 1,000 candidates, and no shard above ~1,300
candidates survived at 4 GB. The `external_sumstats_jk` branch of autoreporting (image
`autorep:20260610.1`) fetches LD lazily, 500 leads at a time in significance order,
drops each entry once its lead is processed or consumed, and never fetches a lead that a
stronger one absorbed as a partner; the narrow fetch is its default, so the WDL input for
it is gone. Rehearsed locally under a 4 GB cap: `seq.16828.8` (9,707 candidates) peaked at
1.37 GB with a report identical to the prefetching image's; the heaviest aptamer
(`seq.2730.58`, 75,660 candidates, 27k of them in the MHC) peaked at 3.4 GB after three
hours on four shared cores, where the prefetch would have needed on the order of 70 GB,
and `seq.17692.2` (55,566) at 1.7 GB. A batch of 500 MHC leads is what sets the peak, so
the `report` task gets 6 GB rather than 4. The config now names that image;
the leftover array comes from the staging script above. The branch still opened the
panel once without a retry at that point, so `maxRetries` and the options file stayed.

**What the third submission taught.** Run `f851d0ec` (image `20260610.1`, the full array,
4 GB, `maxRetries` and the options file) finished 4,421 of 4,838 shards, and 409 of the
417 failures were again OOM kills during grouping: 97 of them before the first LD batch
had returned, on any chromosome, on aptamers with as few as 421 candidates. The lazy
fetch was not the cause. The report command loads the local GWAS catalog (488 MB TSV)
into pandas as strings **before** grouping, because `main.py` builds the catalog
annotation up front; measured in the image that frame is 0.94 GB in the parent, and each
of the four forked LD workers ends up with a private copy of about 0.9 GB, because the
garbage collector walks the inherited objects and copies their pages. Four workers plus
the parent is about 4.5 GB on a VM with 3.8 GB. The local rehearsals never saw it because
they left the catalog out (its allele VCF is unreadable to the VM service account).
Two changes in the `external_sumstats_jk` working tree, image `autorep:20260922.1`
(`20260610.1` with `Scripts/` laid over it, recipe in
`/mnt/disks/data/decode/test/derived_image/`):

- `LocalDB` keeps the catalog path and reads the frame on first query, in the annotation
  stage after the worker pool is gone;
- `LD_FETCH_BATCH` 500 → 100: a batch's partner rows exist in the workers, in flight and
  in the parent at once, and at 500 leads that transient alone reached 2.2 GB on a
  1,047-lead aptamer.

Measured on `seq.12536.46` (1,047 candidates, killed six times in the cloud) with the
catalog in the command and a 3 GB cap: the unpatched image was killed at 3.0 GB during
grouping; the patched image grouped at a 1.0 GB peak and finished at 1.8 GB, the catalog
stage in a single process. The heaviest aptamer, `seq.2730.58`, finished under the same
cap at 2.7 GB with the same 1,109 loci as before, so the `report` task keeps 6 GB: the
4 GB VM has about 2.8 GB left after the OS and the Batch agent. Test suite unchanged: 33
pass, the same 8 pre-existing failures as the branch tip.

**What the fourth submission taught.** Run `aefa2c90` (image `20260922.1`, the full array,
6 GB) had no OOM kill at all; its only failures, 72 of the first 2,937 shards, were the
startup 4xx opens again, all in the first quarter hour, and none was retried. The panel
open in `TabixLD.__init__` and the annotation and sumstat opens in `load_tabix.py` now go
through `data_access.db.open_tabix`, which retries a `gs://` open with exponential backoff
(1 to 32 s, six tries) and fails a local path on the first try; image `autorep:20260922.2`.
`maxRetries` and the options file stay as a second line.

The submitted WDL also lacked the `ld_assume_variant1_indexed` input (the JSON carried it;
Cromwell ignores an input the WDL does not declare), so the 524 reports came from wide
fetches: a median of 65 minutes per shard for a median of 120 candidate leads, and two
shards with 2,000+ candidates were OOM-killed at 4 GB. The shards that never finished have a
median of 1,200 candidates. The WDL in the `external_sumstats_only` checkout of
`~/autoreporting` declares the input.

**Step 2** is [`wdl/create_pseudo_credible_sets.decode.json`](../wdl/create_pseudo_credible_sets.decode.json):
the `ext` flags plus `--data-type pQTL --cell-type plasma`, the phenotype JSON, and output
`deCODE_pQTL_2021_pseudo_credible_sets.mlog10p_2.r2_0.6.tsv.gz` under
`credible_sets/decode_pseudo/`.

Two properties of pQTL sumstats matter for how the thresholds behave:

- the delivery is filtered at p < 0.005 while autoreporting's partner threshold p2 is 0.01,
  so partners with 0.005 ≤ p < 0.01 are never seen. None of them could have become a
  pseudo-CS member: membership needs `mlog10p` within 2 of a lead that is at least 7.3, or
  r² > 0.95 to it, and a variant at r² > 0.95 to a genome-wide-significant lead has an
  expected χ² above 28;
- cis-pQTL leads reach `mlog10p` in the tens of thousands (21,686 for `seq.16828.8`,
  COL6A1, at chr21:46.0 Mb). The dynamic r² threshold `5/χ²` is then far below the
  panel's 0.01 floor, so the lead's clump takes every variant the panel lists as a partner
  at all — and **everything it does not list becomes a locus of its own**. Measured on the
  local rehearsal of `seq.16828.8`: 9,707 candidate leads, 3,948 loci out of clumping, 482
  after `--finngen-variants-only`, of which 471 on chr21 spread over more than 4 Mb around
  COL6A1 (228 within 1 Mb of the top lead, 129 at 1–2 Mb, 114 beyond); step 2 turns them
  into 482 pseudo credible sets, 320 of them singletons, and the second- and third-
  strongest "independent" leads (`mlog10p` 3,610 and 1,148, 20–50 kb from the top lead)
  are simply variants the Finnish panel does not list at r² ≥ 0.01 with it. An expected
  χ² of 5 is reached at r² ≈ 0.00005 to a lead this strong, so the LD shadow of a cis-pQTL
  extends far below anything an LD panel records, and neither step distinguishes shadow
  from signal there. Re-enabling step 2's proximity filter (anchors at χ² > 250, ±1 Mb)
  halves this — 255 sets, 176 singletons, 244 still on chr21 — because leads more than
  1 Mb from the anchor survive and are not anchors themselves. **Decision (2026-09-22):**
  the deCODE step-2 config leaves the proximity filter on — the one production config
  that does — and a pseudo CS at a cis-pQTL is to be read as "the top signal". Secondary
  sets that survive inside a cis window are limited by the LD panel, not evidence of
  independent signals; deCODE's own conditional analysis is the source for those. A
  significance-relative suppression radius and a one-set-per-window rule were considered
  and rejected: both would invent a rule the panel cannot verify, since the panel is why
  the shadow exists.
- the modest aptamers behave like a GWAS: `seq.4876.32` (F9) gave 7 sets, three of them
  cis at chrX:139.5 Mb with real `cs_min_r2`; `seq.7085.81` (CDSN) 58 sets, 23 singletons.
