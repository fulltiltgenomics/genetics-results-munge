# genetics-results-munge

This repository contains a WDL pipeline and scripts to harmonize and process human genetics GWAS and QTL results into a unified TSV format. The output files can be used directly or via APIs (see [genetics-results-api](https://github.com/fulltiltgenomics/genetics-results-api)).

Currently, processing of fine-mapping results from FinnGen, Open Targets and eQTL Catalogue has been implemented.

## Table of Contents

- [Docker image](#docker-image)
- [WDL pipeline](#wdl-pipeline)
- [scripts](#scripts)
  - [Open Targets](#open-targets)
  - [eQTL Catalogue](#eqtl-catalogue)
- [outputs](#outputs)
- [variant annotation](#variant-annotation)

## Docker image

To run the WDL pipeline or scripts, first clone this repository and create a Docker image `genetics-results-munge`:

```
git clone https://github.com/fulltiltgenomics/genetics-results-munge
cd genetics-results-munge
docker build --network host -t genetics-results-munge .
```

## WDL pipeline

The [WDL](https://github.com/openwdl/wdl) munging pipeline [wdl/munge_finngen_finemapping_results.wdl](wdl/munge_finngen_finemapping_results.wdl) is used to process SuSiE fine-mapping results from the [FinnGen fine-mapping pipeline](https://github.com/FINNGEN/finemapping-pipeline). See [FinnGen documentation](https://docs.finngen.fi/finngen-data-specifics/green-library-data-aggregate-data/core-analysis-results-files/finemapping-results-format) for details on fine-mapping results format. The munging pipeline is run per resource: e.g. FinnGen core GWAS and lab value GWAS are processed in separate runs of the pipeline. JSON inputs to the pipeline (FinnGen core GWAS, lab value GWAS, drug purchase GWAS, Olink pQTLs, snRNA-seq eQTLs) are included in [wdl/](/wdl).

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

Get Open Targets data files (requester pays bucket - see [Open Targets website](https://platform.opentargets.org/downloads/credible_set/access) for other download options) and publicly available FinnGen variant annotations:

```
# replace [your_google_project_name] with your project name
BILLING_PROJECT=[your_google_project_name]
mkdir -p data
gcloud storage --billing-project $BILLING_PROJECT cp gs://open-targets-data-releases/25.09/output/credible_set/*.parquet data/
gcloud storage cp gs://finngen-public-data-r12/annotations/finnge_R12_annotated_variants_v1.gz data/
```

Run the Docker container you built, mounting the current directory (genetics-results-munge, root of this repository) in it:

```
docker run -v $(pwd):/munge -it genetics-results-munge /bin/bash
```

Inside the container, select only necessary columns from the variant annotation file to reduce memory use and run the script (this may take some half an hour and it's good to have 16GB of RAM):

```
cd /munge

zcat data/finnge_R12_annotated_variants_v1.gz | cut -f1,1000,1001 | bgzip \
> data/finnge_R12_annotated_variants_v1.small.gz

scripts/create_open_targets_files.sh \
Open_Targets_25.09 \
data \
data/finnge_R12_annotated_variants_v1.small.gz
```

Output files are written under `data`. Only non-FinnGen GWAS traits fine-mapped with SuSiE are included in the output files.

### eQTL Catalogue

Get eQTL Catalogue data and trait metadata, and FinnGen variant annotations if you didn't download them already:

```
scripts/download_eqtl_catalogue_data_and_trait_metadata.sh
gcloud storage cp gs://finngen-public-data-r12/annotations/finnge_R12_annotated_variants_v1.gz data/
```

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
eQTL_Catalogue_R7 \
data \
data/finnge_R12_annotated_variants_v1.small.gz
```

## outputs

Both the WDL pipeline and the scripts give output text files: 1) an uncompressed file for each trait or study, and 2) a bgzip-compressed tabixed file including all traits or studies.

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
most_severe         most severe variant consequence (VEP)
gene_most_severe    gene of most severe consequence
```

For GWAS results, the `trait` and `trait_original` columns are the same and contain the phenotype code of the trait. For QTL results, `trait` contains the gene name while `trait_original` contains the original QTL trait name depending on the dataset, e.g. ENSG gene id.

For eQTL Catalogue, the `trait_original` column contains the QTL trait name and quantification method separated by `|`, e.g. `ENSG00000272211|ge`. Similarly, for eQTL Catalogue, the `cell_type` column contains the name of the cell or tissue and condition separated by `|`, e.g. `plasmacytoid_dendritic_cell|naive`. See [eQTL Catalogue metadata](https://github.com/eQTL-Catalogue/eQTL-Catalogue-resources/blob/master/data_tables/dataset_metadata.tsv) for metadata on the studies in eQTL Catalogue.

For Open Targets, only the lead variant in each credible set has an `mlog10p` value. Also, no Open Targets variants have an `se` value.

There are no spaces in the output files and missing values are represented with `NA`.

## variant annotation

Currently `most_severe` and `gene_most_severe` annotations come from FinnGen variant annotation. This means that variants not present in FinnGen imputation panel have `NA` in the `most_severe` and `gene_most_severe` columns. We're working on more comprehensive variant annotation.
