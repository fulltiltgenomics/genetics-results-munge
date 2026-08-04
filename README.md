# genetics-results-munge

This repository contains a WDL pipeline and scripts to harmonize and process human genetics GWAS and QTL results into a unified TSV format. The output files can be used directly or via APIs (see [genetics-results-api](https://github.com/fulltiltgenomics/genetics-results-api)).

Fine-mapping results from FinnGen, Open Targets and eQTL Catalogue are processed into credible set files. The repository also holds the munging of the other result types served alongside them: colocalization results, gene burden and exome variant results, external GWAS summary statistics and the pseudo credible sets derived from them, allele-specific methylation QTLs, and open chromatin, variant effect and MPRA datasets (see [other datasets](#other-datasets)).

## Table of Contents

- [Docker image](#docker-image)
- [WDL pipeline](#wdl-pipeline)
- [scripts](#scripts)
  - [Open Targets](#open-targets)
  - [eQTL Catalogue](#eqtl-catalogue)
  - [PGC schizophrenia fine-mapping](#pgc-schizophrenia-fine-mapping)
  - [caQTL gene-indexed credible sets](#caqtl-gene-indexed-credible-sets)
  - [gene-indexed peak-to-gene table](#gene-indexed-peak-to-gene-table)
  - [other datasets](#other-datasets)
  - [per-trait burden files](#per-trait-burden-files)
- [outputs](#outputs)
  - [pseudo credible sets](#pseudo-credible-sets)
- [variant annotation](#variant-annotation)

## Docker image

To run the WDL pipeline or scripts, first clone this repository and create a Docker image `genetics-results-munge`:

```
git clone https://github.com/fulltiltgenomics/genetics-results-munge
cd genetics-results-munge
docker build --network host -t genetics-results-munge .
```

## WDL pipeline

The [WDL](https://github.com/openwdl/wdl) munging pipeline [wdl/munge_finngen_finemapping_results.wdl](wdl/munge_finngen_finemapping_results.wdl) is used to process SuSiE fine-mapping results from the [FinnGen fine-mapping pipeline](https://github.com/FINNGEN/finemapping-pipeline). See [FinnGen documentation](https://docs.finngen.fi/finngen-data-specifics/green-library-data-aggregate-data/core-analysis-results-files/finemapping-results-format) for details on fine-mapping results format. The munging pipeline is run per resource: e.g. FinnGen core GWAS and lab value GWAS are processed in separate runs of the pipeline. JSON inputs to the pipeline (FinnGen core GWAS, Kanta lab value GWAS, drug purchase GWAS, Olink and SomaScan pQTLs, UKB-PPP pQTLs, snRNA-seq eQTLs, ATAC-seq caQTLs) are included in [wdl/](/wdl).

[Cromwell](https://cromwell.readthedocs.io/en/latest/) can be used to run the WDL pipeline.

Running the pipeline with the provided JSON inputs requires FinnGen green data and cloud access. However you can run the pipeline on publicly available data locally. 

First make sure you have Java 17+ installed and get Cromwell jar:

```
curl -LO https://github.com/broadinstitute/cromwell/releases/download/90/cromwell-90.jar
```

Download publicly available FinnGen R12 credible sets to `data/`:

```
mkdir -p data
gcloud storage cp gs://finngen-public-data-r12/finemap/summary/finngen_R12_*.SUSIE.snp.filter.tsv data/
gcloud storage cp gs://finngen-public-data-r12/finemap/summary/finngen_R12_*.SUSIE.cred.summary.tsv data/
```

Run the WDL pipeline (the input JSON points to the Docker image you created above, and the metadata file in the input JSON assumes you downloaded the files as above):

```
java -jar cromwell-90.jar run \
wdl/munge_finngen_finemapping_results.wdl \
-i wdl/munge_finngen_finemapping_results.local.json
```

Output files are written under `cromwell-executions/munge_finngen_finemapping_results`.

## scripts

There are scripts to munge data from Open Targets and eQTL Catalogue to the same format as the WDL pipeline does. These scripts are separate because of differences in the input data.

To run the scripts, make sure you have git, Docker and Google Cloud SDK installed, and a [Docker image created](#docker-image).

### Open Targets

Get Open Targets data files (requester pays bucket - see [Open Targets website](https://platform.opentargets.org/downloads/credible_set/access) for other download options) and publicly available FinnGen variant annotations. The credible set and study parquet files go in their own subdirectories of the data directory:

```
# replace [your_google_project_name] with your project name
BILLING_PROJECT=[your_google_project_name]
mkdir -p data/credible_set data/study_metadata
gcloud storage --billing-project $BILLING_PROJECT cp gs://open-targets-data-releases/26.06/output/credible_set/*.parquet data/credible_set/
gcloud storage --billing-project $BILLING_PROJECT cp gs://open-targets-data-releases/26.06/output/study/*.parquet data/study_metadata/
gcloud storage cp gs://finngen-public-data-r13/annotations/finngen_R13_annotated_variants_v0.gz data/
```

Run the Docker container you built, mounting the current directory (genetics-results-munge, root of this repository) in it:

```
docker run -v $(pwd):/munge -it genetics-results-munge /bin/bash
```

Inside the container, cut the annotation down to the four columns the script reads — `#variant`, `AF`, `most_severe` and `gene_most_severe` — and run the script. The field numbers differ between FinnGen releases, so look them up first:

```
cd /munge

zcat data/finngen_R13_annotated_variants_v0.gz | head -1 | tr '\t' '\n' \
| grep -n -x -E '#variant|AF|most_severe|gene_most_severe'

# for R13 those are fields 1, 293, 1009 and 1010
zcat data/finngen_R13_annotated_variants_v0.gz | cut -f1,293,1009,1010 | bgzip \
> data/finngen_R13_annotated_variants_v0.small.gz

scripts/create_open_targets_files.sh \
Open_Targets_26.06 \
data \
data/finngen_R13_annotated_variants_v0.small.gz
```

Output files are written under `data`. Only non-FinnGen GWAS traits fine-mapped with SuSiE are included in the output files. The 26.06 release needs about 10 GB of RAM and takes roughly 15 minutes; the FinnGen R14 annotation used for the released files is not public, so the public R13 annotation above gives slightly lower `aaf` / `most_severe` coverage.

### eQTL Catalogue

Get eQTL Catalogue data and trait metadata, and FinnGen variant annotations if you didn't download them already:

```
scripts/download_eqtl_catalogue_data_and_trait_metadata.sh
gcloud storage cp gs://finngen-public-data-r12/annotations/finnge_R12_annotated_variants_v1.gz data/
```

The datasets to munge are listed in [metadata/eqtl_catalogue_files.tsv](metadata/eqtl_catalogue_files.tsv) (headerless: dataset id, path to its credible set file), and their study, tissue and quantification metadata is read from [metadata/eqtl_catalogue_studies.tsv](metadata/eqtl_catalogue_studies.tsv). The committed file list is the R8 pilot (82 datasets, credible sets as parquet) with the paths of the machine it was run on, so edit them to point to where you downloaded the files. `download_eqtl_catalogue_data_and_trait_metadata.sh` fetches the R7 credible sets from the EBI FTP site (`.tsv.gz`, also accepted) and the phenotype metadata files that the gene name mapping needs.

Run the Docker container you built above, mounting the current directory in it:

```
docker run -v $(pwd):/munge -it genetics-results-munge /bin/bash
```

Inside the container, select only necessary columns from the variant annotation file to reduce memory use and run the script:

```
cd /munge

zcat data/finnge_R12_annotated_variants_v1.gz | cut -f1,1000,1001 | bgzip \
> data/finnge_R12_annotated_variants_v1.small.gz

scripts/create_eqtl_catalogue_files.sh \
eQTL_Catalogue_R8 \
data \
data/finnge_R12_annotated_variants_v1.small.gz
```

### PGC schizophrenia fine-mapping

`munge_pgc_scz_finemap.{py,sh}` converts the published FINEMAP 95 % credible sets of
[Trubetskoy et al. 2022](https://www.nature.com/articles/s41586-022-04434-5) (supplementary table
ST11a) into the credible set format. These are genuine fine-mapping results and are served
alongside the PGC schizophrenia [pseudo credible sets](#pseudo-credible-sets) under the same
resource, as dataset `PGC_SCZ_2022`.

ST11a is on GRCh37 and has no ref/alt alleles, so the GRCh38 locus, alleles and allele frequency
come from the munged wave 3 summary statistics matched on rsid, and the consequence annotation
from the FinnGen variant annotation:

```
scripts/munge_pgc_scz_finemap.sh \
ST11a_95_perc_Credible_Sets.tsv \
data/daner_PGC_SCZ_w3_90_0418b.munged.tsv.gz \
data/R14_annotated_variants_v0.small.gz \
data/pgc_scz_finemap
```

The output has `NA` in `cs_min_r2` throughout and in `aaf` on chromosome X, and its credible sets
can hold several independent signals because FINEMAP was run with more than one causal variant
allowed per locus. See [docs/pgc-scz-finemapping.md](docs/pgc-scz-finemapping.md) for the full
column mapping, the caveats and what is dropped.

### caQTL gene-indexed credible sets

The API's `/credible_sets_by_qtl_gene/{gene}` endpoint reads a gene-indexed copy of a QTL credible set file (`*.qtl.tsv.gz`) that carries the trait's gene coordinates in `trait_chr`/`trait_start`/`trait_end`. For eQTL/pQTL the trait is itself a gene (`create_gene_indexed_qtl_file.py`), but the FinnGen caQTL trait is a chromatin peak, so the gene link comes from the Open4Gene peak-to-gene table: each credible set row is joined to the genes its peak is linked to **in the same cell type** and emitted once per linked gene, with `trait` set to the linked gene symbol and `trait_original` keeping the peak id.

Run inside the Docker container built above (needs ~15 GB free disk in the data dir):

```
docker run -v $(pwd):/munge -it genetics-results-munge /bin/bash

cd /munge
scripts/create_caqtl_gene_indexed_qtl_file.sh data/caqtl
```

Inputs (credible sets, Open4Gene results, GENCODE v32 gene coordinates) are downloaded from GCS if not already in the data dir. Gene coordinates MUST come from the GENCODE version configured for the dataset in the API (`gencode_version: 32` for `finngen_caqtl`), because the API filters returned rows by an exact match on the trait start/end positions. Add `--stage` to upload the result and its tabix index to both profile buckets.

### gene-indexed peak-to-gene table

The Open4Gene peak-to-gene table is tabix-indexed on peak coordinates, which serves the API's `/peak_to_genes/{peak_id}` endpoint but has no inverse. This script writes a second copy with the linked gene's locus appended as three trailing columns and sorted on them, which backs `/gene_to_peaks/{gene}`:

```
docker run -v $(pwd):/munge -it genetics-results-munge /bin/bash

cd /munge
scripts/create_gene_indexed_peak_gene_file.sh data/atacseq
```

The appended columns exist only for the index — the API drops them and re-derives gene coordinates from the GENCODE version the request asks for, so both endpoints return the same columns. Add `--stage` to upload to both profile buckets.

### other datasets

Credible sets are not the only result type munged here. The scripts below write their own schemas, not the credible set columns described in [outputs](#outputs), and each one documents its source files, format assumptions and command line flags in its header comment:

- colocalization (`scripts/coloc/`): FinnGen colocalization credible set and QC files. [scripts/coloc/R14_UPDATE.md](scripts/coloc/R14_UPDATE.md) is the runbook, including where the eQTL Catalogue metadata for trait gene name mapping comes from.
- gene burden and exome variant results (`scripts/genebass/`): Genebass results read from a Hail MatrixTable. Hail is not in `requirements.txt` and not in the Docker image, so the export step runs on Dataproc or on a machine with Hail installed; the shell wrappers do the bgzip/tabix, per-trait split and `mlog10p > 4` filtering. See [per-trait burden files](#per-trait-burden-files) for what the gene burden run produces and why there is no combined unfiltered file.
- external GWAS summary statistics: `munge_pgc.py` (PGC schizophrenia wave 3 daner file), `munge_bip2024.py` (BIP 2024 multi-ancestry, per-ancestry HRC frequencies kept), `munge_gp2.py` (GP2 Parkinson's, already build 38), `munge_ibd.py` (IBD/CD/UC meta-analysis, one input file per chromosome) and `munge_covid.py` (COVID-19 HGI freeze 7). All of them harmonize to the sumstat schema described in [CLAUDE.md](CLAUDE.md) via `scripts/sumstat_utils.py` and most also draw an AF-AF plot against gnomAD with `--gnomad-af-plot`. `filter_gnomad_by_rsid.py` is a standalone helper that pre-filters gnomAD by the rsids of a daner file so a rerun can skip the streaming step.
- non-Genebass exome results: `munge_schema.py` / `munge_schema_variants.py` (SCHEMA), `munge_schema2.py` / `munge_schema2_variants.py` (SCHEMA2 — allele counts only, so log-odds and p-values are derived), `munge_bipex.py` (BipEx2 gene burden), `munge_ibd_exome.py` (IBD/CD/UC gene burden and variant results from one input each) and `munge_ibd_supp_burden.py` / `munge_ibd_supp_variants.py` (the 2026 IBD supplementary tables, whose case and control counts are passed on the command line). These reproduce the Genebass gene burden and exome variant column layouts. `munge_als.py` munges the ALS exome supplementary table into a variant file with its own slightly different column set. The gene burden scripts take `--per-trait-dir` to also emit the [per-trait burden files](#per-trait-burden-files).
- allele-specific methylation QTLs: `normalize_asmqtl.py` left-aligns and trims the indel alleles with `bcftools norm` into a mapping TSV, which `munge_asmqtl.py` then applies while munging either of the two supplementary tables (CpG or MDS methylation QTLs, auto-detected from the header). `munge_asmqtl.sh` runs both steps for both tables.
- open chromatin: `munge_calderon.{py,sh}` (Calderon 2019 immune ATAC-seq, hg19 and therefore the only one of these needing a liftOver to hg38), `munge_catlas.{py,sh}` (Zhang 2021 body-wide snATAC), `munge_epimap.{py,sh}` (EpiMap ChromHMM 18-state calls, active states only), `munge_li_brain.{py,sh}` (Li 2023 brain snATAC), `munge_rosmap.{py,sh}` (Xiong 2023 ROSMAP AD brain snATAC) and `munge_marderstein.{py,sh} --product open_chromatin` (Marderstein 2026 scATAC peaks). `build_li_brain_inputs.py` folds the 44 per-cell-type bed files and the gene-cCRE correlation bedpe into the two files `munge_li_brain.py` expects.
- variant effect: `munge_marderstein.{py,sh}` with `--product chrombpnet` or `--product flare`. Its `--download` path reads Synapse and needs `SYNAPSE_AUTH_TOKEN` in the environment.
- MPRA: `munge_mpra.{py,sh}` reshapes the Siraj 2026 per-variant MPRA annotation from one row per variant to one row per variant and cell line.
- expression: `munge_gtex.py` (GTEx v10 median TPM, written both wide and one row per gene and tissue) and `munge_hpa.py` (HPA immunohistochemistry).
- gene-indexed QTL credible sets: `create_gene_indexed_qtl_file.py` for datasets whose QTL trait is a gene (the caQTL and peak-to-gene variants have their own sections above).
- gene and trait metadata helpers: `gencode_to_gene_pos_tsv.py` (GENCODE GFF3 to a gene position TSV), `create_gene_name_mapping_across_gencode_versions.py` and `kanta_metadata_to_json.py`.

The open chromatin, variant effect and MPRA wrappers write locally by default and only upload to the two profile buckets when given `--stage`. The expression, gene mapping and gene-indexed QTL scripts have their input paths hardcoded at the top of the script.

The sumstat and exome scripts read their input from `--input` (`--input-dir` for the per-chromosome IBD meta-analysis, `--gene-input`/`--variant-input` for the IBD exome) and write next to it unless given an `--output` or `--output-dir`, which may be a `gs://` path — then the file and its tabix index are uploaded there. [run_sumstats.sh](run_sumstats.sh) and [scripts/run_exome.sh](scripts/run_exome.sh) record the invocations for the datasets they cover; their input paths are absolute paths on the machine the munging was run on, so edit them to point at your own copies.

### per-trait burden files

Gene burden results are served two ways, and the munging produces a file for each:

| file | contents | consumer |
|---|---|---|
| `<dataset>_gene_results.munged.tsv.gz` (Genebass: `gene_burden_results.mlog10p_gt4.tsv.gz`) | every trait of the dataset, tabixed on the gene locus (`-s5 -b6 -e6`) | the API's `/gene_based/{gene}` |
| `gene_burden_per_trait/<trait_original>.tsv.gz` | one trait, unfiltered, same tabix index | the API's `/gene_based_results_by_phenotype/{resource}/{trait}` and the BigQuery `gene_burden_results` load |

`trait_original` names the per-trait file, so the IBD burden files are
`inflammatory_bowel_disease` / `ulcerative_colitis` / `crohns_disease`, not the
`IBD`/`UC`/`CD` codes the IBD exome *variant* files use.

Genebass is the one dataset with no combined **unfiltered** file. Its unfiltered
export is one row per gene x annotation x phenotype — 75,767 x 4,501 = ~343M rows —
so `convert_genebass_gene_results.py` exports it unsorted (sorting it in Hail is a
~65 GB shuffle) straight to `.tsv.bgz`, and the two products are built from there:
the `mlog10p_burden > 4` file is small enough to sort with `sort`, and
[scripts/split_burden_per_trait.py](scripts/split_burden_per_trait.py) streams the
export into one gzip temp file per trait and then sorts each one on its own
(~76k rows). The other burden datasets are single-trait and small, so
`write_exome_output()` writes their per-trait copy directly when given
`--per-trait-dir`.

`split_burden_per_trait.py` also works on any other combined result file — pass
`--trait-col 19 --tabix-args "-s2 -b3 -e3"` for the exome variant layout.

## outputs

Both the WDL pipeline and the credible set scripts give output text files: 1) an uncompressed file for each trait or study, and 2) a bgzip-compressed tabixed file including all traits or studies. They also write per-trait credible set statistics (`*.stats.json` plus an aggregate `credible_set_stats.tsv`, see [scripts/credible_set_stats.py](scripts/credible_set_stats.py)). With `create_qtl_file = true` the WDL pipeline additionally writes a gene-indexed copy of the merged file (`*.qtl.tsv.gz`, indexed on the trait's gene locus) for the API's QTL gene lookups.

Columns in all output files:

```
dataset             e.g. "FinnGen_R13" (FinnGen core GWAS) or "QTD000570" (eQTL Catalogue)
data_type           GWAS/eQTL/pQTL/sQTL
trait               see below
trait_original      see below
cell_type           "NA" for GWAS, name of cell type for QTLs
chr                 variant chromosome (a number between 1 and 23)
pos                 variant chromosome position
ref                 variant reference allele
alt                 variant alternative allele
mlog10p             -log10(p-value)
beta                effect size beta for the alternative allele
se                  standard error of effect size
pip                 posterior inclusion probability
cs_id               credible set id
cs_size             credible set size
cs_min_r2           minimum LD r2 between variants in the credible set
aaf                 alternative allele frequency, joined from the variant annotation file
                    (the fine-mapping results only have MAF, so the WDL pipeline leaves MAF
                    here when run without a variant annotation file)
most_severe         most severe variant consequence (VEP)
gene_most_severe    gene of most severe consequence
```

For GWAS results, the `trait` and `trait_original` columns are the same and contain the phenotype code of the trait. For QTL results, `trait` contains the gene name while `trait_original` contains the original QTL trait name depending on the dataset, e.g. ENSG gene id.

For eQTL Catalogue, the `trait_original` column contains the QTL trait name and quantification method separated by `|`, e.g. `ENSG00000272211|ge`. Similarly, for eQTL Catalogue, the `cell_type` column contains the name of the cell or tissue and condition separated by `|`, e.g. `plasmacytoid_dendritic_cell|naive`. See [eQTL Catalogue metadata](https://github.com/eQTL-Catalogue/eQTL-Catalogue-resources/blob/master/data_tables/dataset_metadata.tsv) for metadata on the studies in eQTL Catalogue.

For Open Targets, `mlog10p` is set for about half of the variants: the release carries a per-variant p-value for some studies, and for the rest only the lead variant gets one, from the credible set level p-value. Every credible set has at least one variant with `mlog10p`. No Open Targets variants have an `se` value.

There are no spaces in the output files and missing values are represented with `NA`.

### pseudo credible sets

For external GWAS that ship only summary statistics (no fine-mapping), credible sets are
approximated as *pseudo credible sets*: LD clumps around genome-wide-significant leads,
trimmed by LD and significance and given a heuristic PIP. These are produced in two steps —
the FinnGen autoreporting tool (LD clumping) followed by
[`wdl/create_pseudo_credible_sets.wdl`](wdl/create_pseudo_credible_sets.wdl). See
[docs/pseudo-credible-sets.md](docs/pseudo-credible-sets.md) for how they are defined and
which variants are excluded.

## variant annotation

Currently `most_severe` and `gene_most_severe` annotations come from FinnGen variant annotation. This means that variants not present in FinnGen imputation panel have `NA` in the `most_severe` and `gene_most_severe` columns. We're working on more comprehensive variant annotation.
