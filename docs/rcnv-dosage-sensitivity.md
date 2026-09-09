# rCNV dosage sensitivity, gene associations and segments (Collins et al. 2022)

`scripts/munge_rcnv.{py,sh}` turns three published products of
[Collins et al. 2022](https://doi.org/10.1016/j.cell.2022.06.036) into the suite's `rcnv`
tables: the dosage-sensitivity scores (`--product scores`), the gene-based CNV association
summary statistics (`--product genes`) and the 163 disease-associated large segments
(`--product segments`). This note records the decisions that are not obvious from the code.

**Citation.** Collins RL, Glessner JT, Porcu E, et al. *A cross-disorder dosage
sensitivity map of the human genome.* Cell 2022;185(16):3041-3055.e25.
doi:[10.1016/j.cell.2022.06.036](https://doi.org/10.1016/j.cell.2022.06.036).
Data from Zenodo record [6347673](https://zenodo.org/records/6347673) (v0.2, 2022-03-11),
released under **CC-BY 4.0** — attribution is required wherever these values are served.

## Source

```
https://zenodo.org/records/6347673
  Collins_rCNV_2022.dosage_sensitivity_scores.tsv.gz   --product scores  (386 kB)
  Collins_rCNV_2022.gene_association_sumstats.tar.gz   --product genes   (108 tabixed BEDs)
  Collins_rCNV_2022.sliding_window_sumstats.tar.gz     not munged here
  Collins_rCNV_2022.gene_features_matrix.tar.gz        not loaded

Cell supplement (doi:10.1016/j.cell.2022.06.036), NOT on Zenodo
  mmc3.xlsx  sheet "Table S3"                          --product segments (163 segments)
```

`--product` exists because the record ships four products: the sliding-window sumstats are
separate work and the gene-features matrix is not loaded. `scores`, `genes` and `segments`
are implemented; nothing is stubbed for the others.

`segments` is the one product whose input is **not** on Zenodo and **not** CC-BY: it comes
from the paper's own supplement, under Elsevier's terms. The xlsx is therefore never
committed to this repo, and the script has no download URL for it — Cell and PMC answer a
scripted request with a bot-check page rather than the file (measured on the PMC
`articles/instance/9742861/bin/` path, which returns a 1.8 kB placeholder). Fetch it in a
browser and put it in `--cache-dir`.

The scores file is `#gene pHaplo pTriplo`, 18,641 rows, one per autosomal protein-coding
gene of Gencode v19, no duplicate symbols and no missing-value token. pHaplo and pTriplo are
probabilities in [0,1] at full double precision, and they are written through **verbatim**:
rounding them would let a consumer that recomputes a threshold from the rounded number
disagree with the boolean columns in the same row.

## The thresholds are in the data, not in a doc

`pHaplo >= 0.86` selects 2,987 genes and `pTriplo >= 0.94` selects 1,559 — exactly the counts
in the paper's abstract, which is how the thresholds were established from the file rather
than assumed. Both are written as the boolean columns `haploinsufficient` and
`triplosensitive`, so a query never has to carry a magic number, and all three counts are
asserted by the script.

**Those three assertions cover the input file, not the mapping.** Row count and the two
threshold counts are computed from `pHaplo`/`pTriplo` alone, so every symbol below could land
on the wrong locus and all three would still pass; what they catch is a truncated download or
a re-released Zenodo file. The mapping is asserted separately, on the properties a mapping
regression actually breaks:

- `ensembl_gene_id` is **never NA** — every v19 symbol is in the mapping file's
  `gene_name_19` column, so a missing one means a mapping file that no longer covers v19;
- `ensembl_gene_id` is **unique** — it is the primary key, so two v19 symbols resolving to
  one ENSG is a `pick_gencode` regression;
- **at most 10 symbols carry more than one row** (`MAX_DUPLICATE_SYMBOLS`), which is the
  number Gencode's own merges produce. Anything above it means the HGNC fallback renamed a
  gene onto a symbol another one already holds. Raising the number is a claim that Gencode
  merged two more names, and the failure prints the whole list to check against this doc.

The dataset has no coordinates to spot-query and no second source to cross-check against, so
these are the whole of its automated verification.

## No coordinates, deliberately

The source file has no positions at all. The scores are a per-gene property, so adding
GRCh37 coordinates (or lifting them to GRCh38) would introduce a build-dependent column
where the data has none, and a build mismatch there is the kind of error that joins cleanly
and answers wrongly. The output is keyed by gene symbol and ENSG and is build-independent by
construction; it joins to the suite's GRCh38 gene views on those keys. There is no tabix
index because there is nothing to index — the file is a BigQuery load file only.

## Symbol resolution

This section describes `--product scores`. `--product genes` (below) resolves symbols
through the exact same functions, run separately on its own 17,263-symbol set rather than
sharing this product's 18,641-symbol resolution -- see "Gene associations" for why that
matters and for that run's own counts.

The paper's symbols are Gencode v19 (GRCh37-era); the script prints how many of them are not
what the gene is called now (`symbol changed from v19` in the run summary) rather than this
doc carrying a count that would drift out from under it. Resolution is symbol -> ENSG ->
current symbol:

1. **Gencode.** `gene_name_19 -> ensg` in
   `gencode_gene_name_mapping_49-45-43-39-35-32-19.tsv` (built by
   `scripts/create_gene_name_mapping_across_gencode_versions.py`, staged at
   `gs://daly-genetics-results/mapping_files/`), then the newest non-NA `gene_name_*` for
   that ENSG. All 18,641 v19 symbols are in that column, so `ensembl_gene_id` is **never
   missing** — what misses is only the current *name*. Gencode writes the bare ENSG id as
   `gene_name` for genes it carries but does not name; those are not symbols and are skipped
   when picking the current name.
2. **Ambiguity.** 89 v19 symbols map to two ENSGs, and for 60 of them both candidates carry
   a current name, so the choice is real. Candidates are ordered: the locus whose current
   symbol **is** the v19 symbol, then one Gencode renamed, then a clone-style placeholder
   (matching `^(RP\d+-|AC\d|AL\d|AP\d|CT[ABD]-|Z\d)`), then one Gencode leaves unnamed. Only
   what is still tied breaks on the ENSG id, so a rerun produces the same file. Ordering on
   the ENSG id first — as this did — hands a gene the clone id that happens to sit on its
   symbol with the smaller number: `TUBB3` became `AC092143.1`, and `TLR9`, `SNURF`,
   `NDUFA7`, `LCN6` and `ZNF709` the same way.
3. **HGNC fallback** for genes Gencode no longer names: the v19 symbol is looked up
   **case-insensitively** — Gencode v19 writes `C2ORF15` where HGNC records `C2orf15`, and
   an exact-case lookup misses that record rather than resolving it — among
   HGNC approved symbols, then `prev_symbol`, then `alias_symbol`
   (`hgnc_complete_set.txt`, also staged in `mapping_files/`; `--hgnc-url` fetches it from
   HGNC's public download URL where that bucket is not readable). No other munge in this
   repo reads HGNC, so this is the first script here that needs it.
4. **Unresolved symbols keep their v19 spelling**, with `ensembl_gene_id` still filled.

| route | genes |
|---|---|
| Gencode `gene_name_19` -> current `gene_name` | 18,105 |
| HGNC — v19 symbol still approved | 181 |
| HGNC — `prev_symbol` | 60 |
| HGNC — `alias_symbol` | 6 |
| unresolved, kept as the v19 symbol | 289 |

(The script prints this table on every run; it is the numbers to trust if they ever differ
from the ones above.)

**Most of the unresolved are clone-based placeholder names** — `RP11-*`, `AC0*.1`,
`AL*.1`, `CTD-*`, `CTA-*`, `Z98049.1`. None ever had an approved HGNC symbol, so there is
nothing to update them to; they were never findable by name and remain so, but their ENSG is
exact. The rest are 39 cDNA-era identifiers (`FLJ*`, `FKSG*`, `DKFZP*`, `HUG1`)
and — this is the `C*orf*`/`KIAA*` class the epic expected — 14 `C*orf*` and 3 `KIAA*`
symbols. Those 17 are unresolved *by choice*, not for want of an HGNC record: see the merge
guard below.

**HGNC records a merge as a `prev_symbol`, and following it blindly corrupts the join key.**
`C2orf48`, `C16orf47` and `C17orf47` are listed as previous symbols of `RRM2`, `ZFHX3` and
`SEPTIN4` — but these are distinct ENSGs with their own, very different scores
(`C2orf48` pHaplo 0.41 against `RRM2` 0.97). Renaming them would put two rows with different
pHaplo under one symbol. The fallback therefore refuses any HGNC symbol already held by
another ENSG and keeps the v19 spelling instead. This takes the symbols carrying more than
one row from 44 to **10**.

**The guard has to see every symbol already settled, not only the Gencode-resolved ones.**
Symbols are settled in two passes before any renaming happens — Gencode's current name, then
the v19 symbol HGNC still approves — and the renaming routes then claim as they go. Guarding
against the Gencode picks alone left seven duplicates that the fallback itself manufactured:
`MRC1L1` renamed onto `MRC1`, which the next row was keeping under HGNC approval, and the
same for `ANXA8L1`, `FAM27E2`, `NBPF11`, `NBPF15` and `PPIAL4A`; `GATSL1` and `GATSL2` both
took `CASTOR2` because neither rename was recorded. Which of two competing renames wins is
decided in sorted v19 order, so it does not depend on how Zenodo sorted the file. The winner
is therefore alphabetical rather than principled, and the loser keeps its v19 symbol —
`GATSL2` in the `CASTOR2` case.

The 10 symbols that still carry two rows are irreducible, but not because Gencode gave two
ENSGs the same current name. In each pair one ENSG is a Gencode rename onto that symbol
(`CTAGE5` -> `MIA2` for `ENSG00000150527`) and the other is an ENSG Gencode no longer names at
all (`gene_name_49` is `NA`) that falls through to the unmapped route and keeps its v19
spelling — the same string the renamed row now carries (`ENSG00000150526` stays `MIA2`). The
unmapped route does not consult `claimed`, because inventing a name for a gene Gencode itself
declines to name would be worse than a duplicate. This is the script's own printed list:

```
symbols carrying >1 row:       10
  ['AC010327.2', 'AC020922.1', 'AGAP9', 'AL357673.1', 'AP002884.3', 'FAM25G', 'MIA2',
   'NBPF20', 'PRAMEF15', 'RPS17']
```

**`symbol` is therefore not a primary key; `ensembl_gene_id` is** — it is unique across all
18,641 rows, and the script asserts both.

`symbol_gencode_v19` is kept alongside `symbol` because the sibling rCNV products (gene
association sumstats, sliding windows) are keyed on the v19 symbol, and because a published
result must stay traceable to the identifier the paper actually used.

## Output

One bgzipped TSV, `collins_rcnv_2022_dosage_sensitivity.tsv.gz`, sorted by
`(symbol, symbol_gencode_v19)`:

```
symbol  symbol_gencode_v19  ensembl_gene_id  phaplo  ptriplo  haploinsufficient  triplosensitive
```

`haploinsufficient`/`triplosensitive` are written as `true`/`false` for the BigQuery BOOL
load. The file carries no `dataset` column; the loader injects it with `--const-column`, as
for the HLA combined file.

## Gene associations (`--product genes`)

`Collins_rCNV_2022.gene_association_sumstats.tar.gz` unpacks to 108 tabixed BEDs -- one
`<phenotype>.rCNV.<DEL|DUP>.gene_association.meta_analysis.stats.bed.gz` per phenotype x
CNV-type combination, 54 phenotypes (HP-code or `UNKNOWN`) times `DEL`/`DUP`. Each BED has
17,263 rows, one per autosomal protein-coding Gencode v19 gene, 21 columns per the tarball's
own `README`. Phenotype and CNV type are read from the **file name**, not a column.

**The gene set is the same across all 108 files; row order is not.** Two genes tied on
GRCh37 `(chr, start)` come out in a different relative order between a phenotype's DEL and
DUP file -- confirmed on `HP0000118`, where `APITD1`, `PMF1-BGLAP`, `CHMP3`, `LY75` and
`URGCP-MRPS24` each sit one row off between the two. This output drops chr/start/end
entirely, so the reorder is invisible downstream; what the script actually checks is the
gene **set**, sorted, of every file against the first, plus each file's own row count
against the expected 17,263, backed by the total row count (108 x 17,263 = 1,864,404) for
the whole product.

### Column mapping

GRCh37 `chr`/`start`/`end` are dropped -- coordinates come from `gene_annotations_v` at
query time, the same as `--product scores`. `gene` becomes `symbol_gencode_v19`; `symbol`
and `ensembl_gene_id` are added by the same resolution described above.

| source column (README) | output column | how |
|---|---|---|
| `gene` | `symbol_gencode_v19` | verbatim |
| -- | `symbol` | resolved (see below) |
| -- | `ensembl_gene_id` | resolved (see below) |
| -- | `dataset` | constant `Collins_rCNV_2022` (the Zenodo file-name spelling) |
| file name | `phenotype` | HP-code exactly as in the file name, or `UNKNOWN` |
| file name | `cnv_type` | `DEL` or `DUP` exactly as in the file name |
| `n_nominal_cohorts` | `n_nominal_cohorts` | verbatim |
| `top_cohort` | `top_cohort` | verbatim |
| `cohorts_excluded_from_meta` | `cohorts_excluded` | verbatim (`;`-list, or `NA` if none excluded) |
| `case_freq` | `case_freq` | verbatim |
| `control_freq` | `control_freq` | verbatim |
| `meta_lnOR` | `beta` | `:.3e` (repo's beta-formatting invariant) |
| `meta_lnOR_lower` | `beta_lower` | `:.3e` |
| `meta_lnOR_upper` | `beta_upper` | `:.3e` |
| `meta_z` | `z` | verbatim -- no house rule names `z` |
| `meta_neg_log10_p` | `mlog10p` | rounded to 4 decimals (repo's mlog10p invariant) |
| `meta_neg_log10_fdr_q` | `mlog10_fdr_q` | rounded to 4 decimals |
| `meta_lnOR_secondary` | `beta_secondary` | `:.3e` |
| `meta_lnOR_lower_secondary` | `beta_lower_secondary` | `:.3e` |
| `meta_lnOR_upper_secondary` | `beta_upper_secondary` | `:.3e` |
| `meta_z_secondary` | `z_secondary` | verbatim |
| `meta_neg_log10_p_secondary` | `mlog10p_secondary` | rounded to 4 decimals |
| `meta_neg_log10_fdr_q_secondary` | `mlog10_fdr_q_secondary` | rounded to 4 decimals |

`_secondary` columns repeat the same statistic after excluding the cohort named in
`top_cohort` from the meta-analysis (the README's own definition).

### NA rows are kept, not dropped

`meta_lnOR` onward is `NA` wherever the meta-analysis produced no estimate for that
gene/phenotype/CNV-type row; this correlates with no CNV ever observed in the cohort
(`case_freq = control_freq = 0` -- a few thousand of the NA rows have `case_freq` or
`control_freq` themselves `NA`), but it is **not** implied by `n_nominal_cohorts`:
449,103 rows have `n_nominal_cohorts = 0` with a non-`NA` `beta`, and 17,412 `NA` rows have
`n_nominal_cohorts >= 1`. Filtering on `n_nominal_cohorts` does **not** select the
analysable rows -- filter on `beta IS NOT NULL` (or `mlog10p`) instead. This munge keeps
the NA rows, all-`NA` past `control_freq`, rather than dropping them: "tested, no
meta-analysis" (the row is present, stats are `NA`) has to stay distinguishable from "gene
absent from this file" (the row doesn't exist at all), and only keeping the row can carry
that distinction. The **NA rate is not constant across phenotypes**: 65.2% overall
(1,214,820 of 1,864,404 rows), from 591 rows (3.4%) on the largest phenotype checked
(`HP0000118.DUP`) to 16,570 rows (96.0%) on the smallest (`HP0012447.DEL`); there is no
single "~N NA rows per file" to quote.

### Symbol resolution, run on this product's own gene set

`--product genes` calls the shared `resolve_symbols()` on the 17,263 v19 symbols found in
the gene-association files -- a strict subset of the scores product's 18,641 -- not on the
scores' own resolution. The functions are identical (`pick_gencode`, `resolve_fallback`);
running them on a smaller set can settle an ambiguous symbol differently, because `claimed`
only sees the competitors that are actually in the set being resolved. `ensembl_gene_id` is
asserted never `NA`, same as `--product scores`.

This product is meant to be joined on `ensembl_gene_id`, not `symbol`: a few symbols carry
two genes within a single phenotype/CNV-type file (the same collision pattern as `--product
scores` -- one ENSG Gencode renamed onto the symbol and one ENSG Gencode no longer names at
all) -- check the script's own printed symbol-resolution summary for the current run's
list rather than trusting a hand-typed one here.

Counts from the reference run (script prints these on every run; trust them over this table
if they differ):

| route | genes |
|---|---|
| Gencode `gene_name_19` -> current `gene_name` | 16,888 |
| HGNC -- v19 symbol still approved | 110 |
| HGNC -- `prev_symbol` | 41 |
| HGNC -- `alias_symbol` | 5 |
| unresolved, kept as the v19 symbol | 219 |

### Output

One bgzipped TSV, `collins_rcnv_2022_gene_associations.tsv.gz`, long format (one row per
gene x phenotype x CNV type), grouped by file in the order the 108 BEDs were read (not
sorted -- there is no natural single sort key across dataset/phenotype/cnv_type/symbol that
the API needs, unlike the scores' `(symbol, symbol_gencode_v19)`):

```
dataset  phenotype  cnv_type  symbol  symbol_gencode_v19  ensembl_gene_id
n_nominal_cohorts  top_cohort  cohorts_excluded  case_freq  control_freq
beta  beta_lower  beta_upper  z  mlog10p  mlog10_fdr_q
beta_secondary  beta_lower_secondary  beta_upper_secondary  z_secondary
mlog10p_secondary  mlog10_fdr_q_secondary
```

Reference-run totals: 1,864,404 rows (108 files x 17,263 genes), 54 distinct phenotypes,
`cnv_type` exactly `{DEL, DUP}`, 1,214,820 NA rows (65.16%), `ensembl_gene_id` never `NA`.

## Disease-associated segments (`--product segments`)

`mmc3.xlsx` sheet `Table S3` is the paper's locus-level result table: 163 large rCNV
segments reaching genome-wide significance or FDR in the cross-disorder meta-analysis, one
row each, with pooled frequencies and effect sizes, the associated HPO terms, the 95%
credible interval(s) for the association and the Gencode v19 genes inside the segment.

**Table S4 is deliberately not loaded.** It is the 178-segment *consensus* set: the same 163
rows plus 15 genomic disorders taken from the literature, with extra columns that are
annotations derived from other datasets (Size, Best P-Value, Best ln(OR), Discovery Sig.,
Known GD, gnomAD Constrained Genes, min(LOEUF), min(MisOEUF), ClinGen) rather than results of
this study. Loading it would duplicate 163 rows of results to gain 15 literature loci and a
set of derived columns the suite can compute or does not want.

### Input assertions

163 rows; 69 `DEL` and 94 `DUP`; 88 `Genome-wide` and 75 `FDR`. Both breakdowns are checked,
not just the row count — the row count alone passes on the wrong sheet of a re-released
supplement. Beyond that, each of the three source counter columns is checked against the list
it counts, per row: `# HPOs` vs `associated_hpos`, `# CredInts` vs `credints_grch37`,
`# Genes` vs `genes_gencode_v19`. That is the check that would catch a cell truncated by the
xlsx reader or a stray delimiter, and it is why `split_list` treats an empty cell as **zero**
tokens rather than one: 12 of the 163 segments contain no genes at all (`# Genes` = 0 with an
empty `Genes` cell), and `str.split(';')` on those returns `['']`.

### Column mapping

| source column | output column | note |
|---|---|---|
| — | `dataset` | constant `Collins_rCNV_2022`, as in `--product genes` |
| `Segment ID` | `segment_id` | e.g. `merged_DEL_segment_22q11.21` |
| `CNV Type` | `cnv_type` | `DEL` / `DUP` |
| `Chrom` | `chr` | bare integer 1-22, no `chr` prefix — `datasets.yaml` types every `chr` in the suite `INT64`. All 163 segments are autosomal, so the X-is-23 convention never arises here |
| lifted `Start` | `start` | GRCh38; `NA` where the lift failed |
| lifted `End` | `end` | GRCh38; `NA` where the lift failed |
| `Start` | `start_grch37` | always present |
| `End` | `end_grch37` | always present |
| `Cytoband` | `cytoband` | |
| `Best Significance` | `best_significance` | `Genome-wide` / `FDR` |
| `Pooled Control Freq.` | `control_freq` | |
| `Pooled Case Freq.` | `case_freq` | |
| `Pooled ln(OR)` | `beta` | |
| `Pooled ln(OR) Lower` / `Upper` | `beta_lower` / `beta_upper` | pooled 95% CI |
| `Min. ln(OR)` / `Max. ln(OR)` | `beta_min` / `beta_max` | across the segment's associated phenotypes |
| `# HPOs` | `n_hpos` | |
| `Associated HPOs` | `associated_hpos` | `;`-joined, colon stripped (below) |
| `# CredInts` | `n_credints` | |
| lifted `CredInts` | `credints` | `;`-joined GRCh38 `chr:start-end` |
| `CredInts` | `credints_grch37` | `;`-joined, as published |
| `CredInt Size` | `credint_size` | the source's own **GRCh37** total, in bp, of the credible intervals; it is not recomputed from the lifted ones |
| `# Genes` | `n_genes` | |
| resolved `Genes` | `genes` | current symbols |
| `Genes` | `genes_gencode_v19` | as published |
| resolved `Genes` | `gene_ensembl_ids` | ENSG per gene, same order |

The `beta*` columns are formatted `:.3e` per the repo's statistics invariant, the same as
`--product genes`. `control_freq` and `case_freq` are written at full round-trip precision:
openpyxl hands back Python floats, so unlike the text-file products there is no source
string to pass through verbatim, and `repr` is the spelling that reads back as the same
double.

### The `;` delimiter is a contract with the loader

`associated_hpos`, `credints`, `credints_grch37`, `genes`, `genes_gencode_v19` and
`gene_ensembl_ids` are `;`-joined strings, and the BigQuery loader splits them into
`ARRAY<STRING>`. `;` is the source's own delimiter and it stays `;` — changing it would
silently change what the loader produces. Nothing has to trust that no value contains one:
the per-row count checks above fail the run if a list ever gains or loses an element.

`gene_ensembl_ids` is positionally aligned with `genes` and `genes_gencode_v19`, so
`genes[i]`, `genes_gencode_v19[i]` and `gene_ensembl_ids[i]` are the same gene. The column
carries `NA` for a gene the mapping cannot place, but `munge_segments` asserts it never does
for this input, naming the offending symbols and raising `SystemExit` if it ever did: all
1,711 distinct v19 symbols across the 163 segments are in the mapping file's `gene_name_19`
column.

### HPO ids are spelled without the colon

The supplement writes `HP:0012759`; the output writes `HP0012759`. That is the spelling the
Zenodo file names use, hence the `phenotype` column of `--product genes` and the `phenocode`
of `configs/rcnv_pheno.json` in the suite. Writing the paper's spelling here would mean every
join from a segment to a phenotype or to a gene association carried a `REPLACE(hpo, ':', '')`.
`UNKNOWN` (the unaffected-phenotype-unknown group, present on 33 of the 163 segments) passes
through as itself.

### Symbol resolution

Identical to the other two products — `resolve_symbols()`, run on this product's own set of
v19 symbols. Reference run: 1,711 distinct v19 symbols over 2,200 gene mentions; 1,632
resolved through Gencode, 51 through HGNC, 28 left with their v19 spelling (clone-style
placeholder names Gencode no longer carries); 204 distinct symbols changed, 254 gene
mentions.

### liftOver GRCh37 -> GRCh38

The procedure is the one `scripts/rcnv_liftover_windows.py` established and measured for the
sliding windows, imported from it rather than re-implemented: whole-interval BED4 with a key
in the name column, one UCSC `liftOver` run at defaults (minMatch 0.95, no `-multiple`), then
drop anything that mapped more than once, landed on another chromosome, or changed length by
more than the tolerance. The tolerance is **±10% of the interval's own GRCh37 length** — the
windows measurement's 180-220 kb around a fixed 200 kb window, expressed as a fraction,
because segments run from 200 kb to 10.3 Mb and credible intervals are smaller still.

The segment span and every credible interval go through one `liftOver` call. That script's
second, endpoint-only pass is not repeated here: it never rescues an interval, it exists to
attribute a measurement's failure to one end, and liftOver's own reason string is enough for
a 163-row table. The four filters decide the dropped set, so both products drop the same
intervals for the same reasons.

The binary and the chain are not vendored; `--download` fetches them into `--cache-dir`
(`--liftover-bin` / `--chain` override):

```
https://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/liftOver
https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz
```

**A failure keeps the row.** The segment or credible interval that did not lift gets `NA` in
the GRCh38 column and keeps its GRCh37 coordinates; nothing is dropped. Inside `credints` an
`NA` holds the failed interval's place, so the GRCh38 list stays element-for-element aligned
with `credints_grch37` and both still carry `n_credints` entries.

Reference run: **153 of 163 segment spans and 214 of 225 credible intervals lifted.** The
10 segments left with NULL GRCh38 coordinates:

| segment | GRCh37 | liftOver verdict | credible intervals also lost |
|---|---|---|---|
| `merged_DUP_segment_1p36.32-p36.33` | 1:1,690,000-5,380,000 | Partially deleted in new | 1 of 5 |
| `merged_DUP_segment_1p11.2-p12` | 1:119,360,000-120,990,000 | Split in new | 1 of 1 |
| `merged_DUP_segment_1q21.1-q21.2` | 1:145,290,000-147,820,000 | Split in new | 0 of 2 |
| `merged_DUP_segment_8q24.3` | 8:145,540,000-146,360,000 | Split in new | 1 of 2 |
| `merged_DUP_segment_10q26.3` | 10:135,260,000-135,530,000 | lifted 340,926 bp vs 270,000 bp | 1 of 1 |
| `merged_DEL_segment_13q34` | 13:114,490,000-115,160,000 | Partially deleted in new | 1 of 1 |
| `merged_DUP_segment_15q11.2-q13.3` | 15:22,740,000-32,530,000 | Split in new | 0 of 4 |
| `merged_DEL_segment_17p13.3` | 17:130,000-1,900,000 | Split in new | 1 of 2 |
| `merged_DUP_segment_22q11.21` | 22:18,560,000-21,540,000 | Split in new | 1 of 1 |
| `merged_DEL_segment_22q11.21` | 22:18,820,000-21,540,000 | Split in new | 1 of 1 |

Three further segments lifted themselves but lost one credible interval each:
`merged_DEL_segment_1q43-q44` (1 of 2), `merged_DEL_segment_9q34.3` (1 of 3),
`merged_DUP_segment_16p13.3_A` (1 of 2).

**6.1% of segments fail, against 1.8% of the sliding windows, and that is not a regression in
the chain.** The ten failed spans — 1p36.32-p36.33 DUP, 1p11.2-p12, 1q21.1-q21.2, 8q24.3,
10q26.3, 13q34, 15q11.2-q13.3, 17p13.3, 22q11.21 DEL, 22q11.21 DUP — are each a recurrent
genomic disorder *because* they are flanked by segmental duplications, which is exactly the
sequence hg19 and hg38 rearranged; "Split in new" is liftOver saying the interval no longer
has one image. Three further segments lift their span but lose a credible interval — a
credible-interval-only failure — 16p13.3_A, 9q34.3 and 1q43-q44. The failure rate is higher
here than for the windows because these intervals are 10-50x longer and are, by construction,
the rearranged ones. A measurement recorded during this task and **not** acted on: lifting
the two endpoints independently and rebuilding the span would recover 4 of the 10 segments
(1p36.32-p36.33 DUP, 15q11.2-q13.3, 17p13.3 DEL, 22q11.21 DUP) and 3 of the 11 credible
intervals. Of the remaining 6 segments, 2 are rejected by the ±10% filter
(`merged_DEL_segment_22q11.21` at -13.5%, `merged_DUP_segment_1q21.1-q21.2` at -13.1%), which
is the tolerance doing its job rather than a case to widen it, and 4 have an endpoint where
liftOver reports "Deleted in new" (1p11.2-p12, 8q24.3, 10q26.3, 13q34). Changing the
procedure is a decision for the epic, not for this munge: it would make this product drop a
different set of intervals than the sliding-window product does from the same chain.

## Running it

```bash
# download inputs to $HOME/rcnv_munge/cache, produce locally, no upload
scripts/munge_rcnv.sh

# gene associations instead of the dosage-sensitivity scores
PRODUCT=genes scripts/munge_rcnv.sh

# the 163 segments; mmc3.xlsx must already be in $HOME/rcnv_munge/cache
PRODUCT=segments scripts/munge_rcnv.sh

# produce and publish to both profile buckets
scripts/munge_rcnv.sh --stage
```

`--download` for `--product genes` fetches and untars
`Collins_rCNV_2022.gene_association_sumstats.tar.gz` into `<cache-dir>/`; this path has not
been exercised against a live Zenodo download from this host (see the script's
`fetch_gene_assoc` docstring). `--gene-assoc-dir` points the script at an already-unpacked
copy of the 108 BEDs directly, bypassing both the download and the untar.

`--stage` attempts **both** destinations regardless of whether the first succeeds, reports
every failure at the end and exits non-zero if any failed — the two buckets are in different
projects and a host commonly has credentials for only one of them, so stopping at the first
error meant the second copy could never be written from such a host.

`PRODUCT`, `CACHE_DIR`, `OUT_DIR`, `DATASET_ID`, `OUTPUT`, `GCS_FINNGEN` and `GCS_DALY`
override the defaults; `--stage` is never implied. Staged layout:

```
gs://finngen-commons/results_api_data/rcnv/collins_rcnv_2022/    # finngen profile
gs://daly-genetics-results/rcnv/collins_rcnv_2022/               # daly profile
  collins_rcnv_2022_dosage_sensitivity.tsv.gz
  collins_rcnv_2022_gene_associations.tsv.gz
  collins_rcnv_2022_segments.tsv.gz
```

`--product segments` is the one product with an input the script cannot fetch: put
`mmc3.xlsx` in `--cache-dir` (or pass `--segments-xlsx`) first. `--download` still fetches
the gencode mapping, the HGNC set, and the liftOver binary and chain.

The BigQuery table and view are defined in `genetics-results-db`, and the dataset is
declared in `genetics-results-suite`'s `configs/datasets.yaml`.
