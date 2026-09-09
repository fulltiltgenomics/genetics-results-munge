# rCNV dosage sensitivity (Collins et al. 2022)

`scripts/munge_rcnv.{py,sh}` turns the published dosage-sensitivity scores of
[Collins et al. 2022](https://doi.org/10.1016/j.cell.2022.06.036) into the suite's
`rcnv` table. This note records the decisions that are not obvious from the code.

**Citation.** Collins RL, Glessner JT, Porcu E, et al. *A cross-disorder dosage
sensitivity map of the human genome.* Cell 2022;185(16):3041-3055.e25.
doi:[10.1016/j.cell.2022.06.036](https://doi.org/10.1016/j.cell.2022.06.036).
Data from Zenodo record [6347673](https://zenodo.org/records/6347673) (v0.2, 2022-03-11),
released under **CC-BY 4.0** — attribution is required wherever these values are served.

## Source

```
https://zenodo.org/records/6347673
  Collins_rCNV_2022.dosage_sensitivity_scores.tsv.gz   --product scores  (386 kB)
  Collins_rCNV_2022.gene_association_sumstats.tar.gz   not munged here
  Collins_rCNV_2022.sliding_window_sumstats.tar.gz     not munged here
  Collins_rCNV_2022.gene_features_matrix.tar.gz        not loaded
```

`--product` exists because the record ships four products and the later ones are separate
work. Only `scores` is implemented; nothing is stubbed for the others.

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

## Running it

```bash
# download inputs to $HOME/rcnv_munge/cache, produce locally, no upload
scripts/munge_rcnv.sh

# produce and publish to both profile buckets
scripts/munge_rcnv.sh --stage
```

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
```

The BigQuery table and view are defined in `genetics-results-db`, and the dataset is
declared in `genetics-results-suite`'s `configs/datasets.yaml`.
