# Running the coloc munged-file update (FinnGen R14)

How to regenerate the two coloc munged files from the R14 colocalization results:

- `gs://finngen-commons/results_api_data/coloc/coloc.credsets.munged.tsv.gz`
- `gs://finngen-commons/results_api_data/coloc/colocQC.munged.tsv.gz`

The scripts in `scripts/coloc/` (`munge_coloc_credset_file.{py,sh}`,
`munge_coloc_file.{py,sh}`, `munge_coloc_file_map_traits.py`) are written for the R14
input layout (explicit `dataset`/`tissue`/`quant` columns, precomputed `mlogp`,
explicit hit beta/se/p/mlogp). They do not support the older R13 layout.

## Inputs

| file | path |
|---|---|
| credsets | `gs://finngen-production-library-green/finngen_R14/finngen_R14_analysis_data/colocalization/coloc.credsets.tsv.gz` |
| QC / H4 | `gs://finngen-production-library-green/finngen_R14/finngen_R14_analysis_data/colocalization/colocQC.tsv.gz` |

`coloc.H4_tables.tsv.gz` and `unfiltered_summaries/` in the same folder are not used.

## eQTL Catalogue metadata (R7 + R8 pilot)

R14 colocalized against the **eQTL Catalogue R8 pre-release ("r8_beta", Jan 2026)** on
top of R7 (new studies `INTERVAL_RNA`, `INTERVAL_RNA_WGS`, the `*_reannotated`
single-cell studies, and the `majiq` quant). The merged metadata and the reusable
phenotype-metadata files are already staged under
`/mnt/disks/data/eqtl_catalogue_r8_beta/`:

- `dataset_metadata_r7_r8_merged.tsv` — **pass this as the metadata argument**
  (r7 datasets + the R8 dataset_ids not in r7; 9-column form).
- `dataset_metadata_r8_beta.tsv` — raw R8 pilot table (kept for reference).
- `metadata/` — symlinks to the r7 phenotype-metadata files (gene_counts, exon, tx,
  txrev, microarray, aptamer, per-study leafcutter). `munge_coloc_file_map_traits.py`
  reads every `*_phenotype_metadata.tsv.gz` here to map molecular-trait IDs to genes.

To refresh the metadata in a future release, re-fetch
`data_tables/dataset_metadata_r8_beta.tsv` (and the matching `_r7` table) from the
[eQTL-Catalogue-resources repo](https://github.com/eQTL-Catalogue/eQTL-Catalogue-resources/tree/master/data_tables)
and rebuild the merged table and the `metadata/` symlink directory.

### INTERVAL_RNA sQTL metadata (resolved)

The per-study splicing phenotype metadata for the INTERVAL studies — leafcutter
(`QTD001001`, `QTD001004`) and majiq (`QTD001002`, `QTD001005`) — was initially
unpublished. It is now available at <https://zenodo.org/records/17956387> and staged
in `/mnt/disks/data/eqtl_catalogue_r8_beta/metadata/` as
`leafcutter_QTD001001_Ensembl_105_phenotype_metadata.tsv.gz`,
`leafcutter_QTD001004_…`, `majiq_QTD001002_…`, `majiq_QTD001005_…` (the `*_phenotype_metadata.tsv.gz`
naming is required so `munge_coloc_file_map_traits.py`'s glob picks them up). With
these present, INTERVAL leafcutter `trait2` maps to gene names (was 3,351 unmapped
cluster ids), and all eQTL Catalogue `trait2` values are mapped (0 unmapped). To
refresh in a future release, re-download from that Zenodo record into the metadata
dir under the same naming and rerun the QC munge.

## Separate FinnGen QTL coloc files — leave as-is

The output bucket also has `FinnGen-R12.eQTL.*`, `FinnGen-KANTA.eQTL.*`,
`FinnGen-R12.caQTL.*`, `FinnGen-KANTA.caQTL.*` munged files. These come from separate
R12/R13-era QTL coloc runs and have no R14 equivalent (single-cell / FinnLiver are
folded into the combined R14 files). **Do not regenerate, delete, or touch their API
profile entries.** Only the two combined files above are updated.

## Memory and disk

- The credible-set input is large (~315 MB gzipped, ~11–12 M rows). Run from a
  working directory on a roomy disk — use `/mnt/disks/data` (≈100 GB free), **not**
  the root filesystem (only a few GB free). The shell wrappers sort with `-T .`, so
  temp files land in the working directory.
- `munge_coloc_credset_file.py` and `munge_coloc_file.py` stream
  (`scan_csv` → `sink_csv`); peak RAM stays in the low-GB range.
- `munge_coloc_file_map_traits.py` reads the (already filtered, much smaller) QC
  intermediate eagerly and loads all phenotype-metadata files into a lookup — a few
  GB peak. The machine has ~28 GB free, which is comfortable. If memory is ever
  constrained, run the credset and QC jobs sequentially (not in parallel) so only one
  large job is resident at a time.

## Run

```bash
RUN=/mnt/disks/data/coloc_r14_run        # roomy scratch dir
mkdir -p "$RUN" && cd "$RUN"
META=/mnt/disks/data/eqtl_catalogue_r8_beta/dataset_metadata_r7_r8_merged.tsv
SCRIPTS=/home/jkarjala/suite/genetics-results-munge/scripts/coloc

# fetch inputs locally
gcloud storage cp gs://finngen-production-library-green/finngen_R14/finngen_R14_analysis_data/colocalization/coloc.credsets.tsv.gz .
gcloud storage cp gs://finngen-production-library-green/finngen_R14/finngen_R14_analysis_data/colocalization/colocQC.tsv.gz .

# credible sets  -> coloc.credsets.munged.tsv.gz (+ .tbi)
zcat coloc.credsets.tsv.gz > coloc.credsets.tsv
"$SCRIPTS/munge_coloc_credset_file.sh" coloc.credsets.tsv "$META" coloc.credsets.munged.tsv

# QC / H4       -> colocQC.munged.tsv.gz (+ .tbi)
zcat colocQC.tsv.gz > colocQC.tsv
"$SCRIPTS/munge_coloc_file.sh" colocQC.tsv "$META" colocQC.munged.tsv
```

The wrappers run the python munge, an awk NA-completeness check on required columns,
then sort + bgzip + tabix. The QC wrapper also runs the trait gene-name mapping pass
and writes the final `colocQC.munged.tsv.gz`.

## Upload (only when validated)

Do **not** upload until the verification below passes.

```bash
gcloud storage cp coloc.credsets.munged.tsv.gz{,.tbi} gs://finngen-commons/results_api_data/coloc/
gcloud storage cp colocQC.munged.tsv.gz{,.tbi}        gs://finngen-commons/results_api_data/coloc/
```

## Verify

- Both `.sh` wrappers exit 0 (the awk NA check fails the run if a required column is
  `NA`).
- Output headers match the current munged files exactly:
  - credsets: `#dataset, data_type, trait, trait_original, cell_type, chr, pos, ref,
    alt, mlog10p, beta, se, pip, cs_id`
  - QC: `#dataset1, dataset2, data_type1, data_type2, trait1, trait1_original, trait2,
    trait2_original, cell_type1, cell_type2, cs1_id, cs2_id, hit1, hit2, hit1_beta,
    hit1_mlog10p, hit2_beta, hit2_mlog10p, chr, region_start_min, region_end_max,
    PP.H0–H4.abf, nsnps, nsnps1, nsnps2, cs1_log10bf, cs2_log10bf, clpp, clpa,
    cs1_size, cs2_size, cs_overlap, topInOverlap`
- `data_type` populated for all rows; `cell_type` populated for QTL rows and `NA`/
  empty for GWAS rows.
- QTL `trait`/`trait2` are gene names (not raw `ENSG…`/`ENST…`/cluster IDs); with the
  INTERVAL sQTL metadata staged (see above) there are 0 unmapped QTL trait2 values.
- `tabix` indexes built (`.tbi` present for both files).
