# genetics-results-munge

This repository contains a WDL pipeline and scripts to harmonize and process human genetics GWAS and QTL results into a unified TSV format. The output files can be used directly or via APIs (see [genetics-results-api](https://github.com/fulltiltgenomics/genetics-results-api)).

Fine-mapping results from FinnGen, Open Targets and eQTL Catalogue are processed into credible set files. The repository also holds the munging of the other result types served alongside them: colocalization results, gene burden and exome variant results, external GWAS summary statistics and the pseudo credible sets derived from them, allele-specific methylation QTLs, rare-CNV dosage sensitivity, curated gene-disease associations, and open chromatin, variant effect and MPRA datasets (see [other datasets](#other-datasets)).

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

Get the credible set and study parquet files of an Open Targets release into their own subdirectories of the data directory. They are on the EBI FTP site and in a requester pays bucket (see the [Open Targets website](https://platform.opentargets.org/downloads/credible_set/access) for other download options):

```
mkdir -p data/credible_set data/study_metadata
base=https://ftp.ebi.ac.uk/pub/databases/opentargets/platform/26.09/output
for pair in credible_set:credible_set study:study_metadata; do
    src=${pair%%:*}; dst=${pair##*:}
    curl -s $base/$src/ | grep -oE 'href="[^"]+\.parquet"' | cut -d'"' -f2 \
    | sed "s#^#$base/$src/#" | xargs -P 4 -n 1 wget -q -c -P data/$dst
done

# or, replacing [your_google_project_name] with your project name
BILLING_PROJECT=[your_google_project_name]
gcloud storage --billing-project $BILLING_PROJECT cp gs://open-targets-data-releases/26.09/output/credible_set/*.parquet data/credible_set/
gcloud storage --billing-project $BILLING_PROJECT cp gs://open-targets-data-releases/26.09/output/study/*.parquet data/study_metadata/
```

Run the Docker container you built, mounting the current directory (genetics-results-munge, root of this repository) in it:

```
docker run -v $(pwd):/munge -it genetics-results-munge /bin/bash
```

Inside the container, run the script:

```
cd /munge

scripts/create_open_targets_files.sh \
Open_Targets_26.09 \
data
```

Output files are written under `data`. Only non-FinnGen GWAS traits fine-mapped with SuSiE are included in the output files. `aaf`, `most_severe` and `gene_most_severe` are `NA` at this point: the release has no allele frequency, and consequence is stamped onto the munged files afterwards by `scripts/annotate_resource.sh`, which also regenerates the credible set stats that depend on it. The per-study files are named by study accession (`<accession>.SUSIE.munged.tsv`), which is the `trait_original` column.

### eQTL Catalogue

Get eQTL Catalogue data and trait metadata, and the publicly available FinnGen variant annotations, which `aaf` is read from:

```
scripts/download_eqtl_catalogue_data_and_trait_metadata.sh
gcloud storage cp gs://finngen-public-data-r12/annotations/finnge_R12_annotated_variants_v1.gz data/
```

The datasets to munge are listed in [metadata/eqtl_catalogue_files.tsv](metadata/eqtl_catalogue_files.tsv) (headerless: dataset id, path to its credible set file), and their study, tissue and quantification metadata is read from [metadata/eqtl_catalogue_studies.tsv](metadata/eqtl_catalogue_studies.tsv). The committed file list is the R8 pilot (82 datasets, credible sets as parquet) with the paths of the machine it was run on, so edit them to point to where you downloaded the files. `download_eqtl_catalogue_data_and_trait_metadata.sh` fetches the R7 credible sets from the EBI FTP site (`.tsv.gz`, also accepted) and the phenotype metadata files that the gene name mapping needs.

Run the Docker container you built above, mounting the current directory in it:

```
docker run -v $(pwd):/munge -it genetics-results-munge /bin/bash
```

Inside the container, cut the annotation down to the two columns the script reads, `#variant` and `AF`, and run the script. The field numbers differ between FinnGen releases, so look them up first:

```
cd /munge

zcat data/finnge_R12_annotated_variants_v1.gz | head -1 | tr '\t' '\n' \
| grep -n -x -E '#variant|AF'

# for R12 those are fields 1 and 291
zcat data/finnge_R12_annotated_variants_v1.gz | cut -f1,291 | bgzip \
> data/finnge_R12_annotated_variants_v1.small.gz

scripts/create_eqtl_catalogue_files.sh \
eQTL_Catalogue_R8 \
data \
data/finnge_R12_annotated_variants_v1.small.gz
```

`most_severe` and `gene_most_severe` are `NA` in the output: consequence is stamped onto the munged files afterwards by `scripts/annotate_resource.sh`.

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
- external GWAS summary statistics: `munge_pgc.py` (PGC schizophrenia wave 3 daner file), `munge_bip2024.py` (BIP 2024 multi-ancestry, per-ancestry HRC frequencies kept), `munge_gp2.py` (GP2 Parkinson's, already build 38), `munge_ibd.py` (IBD/CD/UC meta-analysis, one input file per chromosome), `munge_covid.py` (COVID-19 HGI freeze 7) and `munge_aih.py` (AIH autoimmune hypothyroidism meta-analysis, one file per phenotype, already build 38 with chrX coded 23; takes several `--input` files at once and resolves `A1`/`A0` against gnomAD). All of them harmonize to the sumstat schema described in [CLAUDE.md](CLAUDE.md) via `scripts/sumstat_utils.py` and most also draw an AF-AF plot against gnomAD with `--gnomad-af-plot`. `filter_gnomad_by_rsid.py` is a standalone helper that pre-filters gnomAD by the rsids of a daner file so a rerun can skip the streaming step.
- external pQTL summary statistics: `munge_decode_pqtl.py` splits the deCODE 2021 plasma pQTL delivery (one p < 0.005 file for all 4,907 SomaScan aptamers, already build 38 and aligned to gnomAD, carrying no allele frequency) into one munged sumstat per aptamer, keeping one row per variant where the alignment folded two indel representations together; `decode_pqtl_phenotypes.py` writes the aptamer-to-gene metadata (the autoreporting phenotype-info TSV and the pseudo-CS phenotype JSON) from the same mapping file the FinnGen SomaScan credible sets are named with. These exist only to feed the [pseudo credible sets](#pseudo-credible-sets); the API does not serve them.
- non-Genebass exome results: `munge_schema.py` / `munge_schema_variants.py` (SCHEMA), `munge_schema2.py` / `munge_schema2_variants.py` (SCHEMA2 — allele counts only, so log-odds and p-values are derived), `munge_bipex.py` (BipEx2 gene burden), `munge_ibd_exome.py` (IBD/CD/UC gene burden and variant results from one input each) and `munge_ibd_supp_burden.py` / `munge_ibd_supp_variants.py` (the 2026 IBD supplementary tables, whose case and control counts are passed on the command line). `munge_brava.{py,sh}` munges the BRaVa cross-ancestry exome meta-analysis gene burden results, one source file per (phenotype, meta-analysis stratum): only the `class=Burden` / `type=Inverse variance weighted` rows are kept, because they are the only ones carrying both an effect size and its standard error, and the mask and `max_MAF` are pasted into the `annotation` string as the source spells them. Trait names and sample sizes come from the `brava_pheno.json` that `brava_phenotypes.py` builds, not from the result files, and `--phenotypes`/`--strata` choose what a run covers — every trait of one run shares one combined file, so a partial re-run replaces it rather than adding to it. `--stage` publishes the combined filtered file and the per-trait files to `gs://daly-genetics-results/exome_results/brava/`. These reproduce the Genebass gene burden and exome variant column layouts. `munge_als.py` munges the ALS exome supplementary table into a variant file with its own slightly different column set. The gene burden scripts take `--per-trait-dir` to also emit the [per-trait burden files](#per-trait-burden-files).
- count-based exome results: `munge_asc.py` munges the Autism Sequencing Consortium 2026 release (gene and variant files, medRxiv 10.64898/2026.08.24.26360398) into three BigQuery-only products — a long per-gene table of counts by variant class and inheritance mode, a per-gene Bayes factor / FDR table, and a per-variant allele count table. The release carries no p-value, effect size or allele frequency and the script derives none: every statistic and count is passed through as the source spells it (the statistics are read as strings for that reason), so these files do not reproduce the Genebass column layouts and are not served by the results API. Gene symbols and coordinates come from GENCODE v29, which the release's VEP 95 gene ids match completely. `--output-dir` may be a `gs://` path; the loader is `genetics-results-db/scripts/load_asc.sh`.
- allele-specific methylation QTLs: `normalize_asmqtl.py` left-aligns and trims the indel alleles with `bcftools norm` into a mapping TSV, which `munge_asmqtl.py` then applies while munging either of the two supplementary tables (CpG or MDS methylation QTLs, auto-detected from the header). `munge_asmqtl.sh` runs both steps for both tables.
- open chromatin: `munge_calderon.{py,sh}` (Calderon 2019 immune ATAC-seq, hg19 and therefore the only one of these needing a liftOver to hg38), `munge_catlas.{py,sh}` (Zhang 2021 body-wide snATAC), `munge_epimap.{py,sh}` (EpiMap ChromHMM 18-state calls, active states only), `munge_li_brain.{py,sh}` (Li 2023 brain snATAC), `munge_rosmap.{py,sh}` (Xiong 2023 ROSMAP AD brain snATAC) and `munge_marderstein.{py,sh} --product open_chromatin` (Marderstein 2026 scATAC peaks). `build_li_brain_inputs.py` folds the 44 per-cell-type bed files and the gene-cCRE correlation bedpe into the two files `munge_li_brain.py` expects. All of these, and the variant effect and MPRA munges below, take their chromosome encoding and their sort/bgzip/tabix write from `scripts/peak_utils.py` — the peak-family counterpart to `scripts/sumstat_utils.py`.
- variant effect: `munge_marderstein.{py,sh}` with `--product chrombpnet` or `--product flare`. Its `--download` path reads Synapse and needs `SYNAPSE_AUTH_TOKEN` in the environment.
- MPRA: `munge_mpra.{py,sh}` reshapes the Siraj 2026 per-variant MPRA annotation from one row per variant to one row per variant and cell line.
- classical HLA allele associations: `munge_hla.{py,sh}` rewrites a FinnGen release's imputed-HLA results (which model an allele as a variant, `ref='<absent>'` / `alt='A*02:01'`) into explicit `gene`/`allele` columns, joins each allele's imputation INFO from the `.snpstats` sidecar onto every row, and drops phenotypes that have no entry in the release phenotype metadata. It emits both artifacts the suite needs: per-phenotype tabix files for the results-api `/hla` endpoints and one combined TSV for the BigQuery `hla_associations` load, which is the only place the cross-phenotype question is answerable. See [docs/hla-allele-associations.md](docs/hla-allele-associations.md).
- rare-CNV dosage sensitivity, gene associations, segments and sliding windows: `munge_rcnv.{py,sh} --product scores` maps the Collins 2022 pHaplo/pTriplo scores from their Gencode v19 symbols to current symbols via the cross-version gene name mapping, falling back to HGNC for the genes Gencode no longer names, and writes the paper's published thresholds into the file as boolean columns. `--product genes` reshapes the same paper's 108 gene-based CNV association BEDs (54 phenotypes x DEL/DUP, 17,263 genes each) into one long table, symbol-resolved through the same functions, keeping the ~65% of rows that carry no meta-analysis as NA rows rather than dropping them. `--product segments` reads the 163 disease-associated segments from the Cell supplement's `mmc3.xlsx` (`Table S3`, the one input with no download URL — place it in the cache dir by hand) and lifts the segment span and every 95% credible interval to GRCh38 with the sliding-window measurement's own liftOver procedure, falling back for an interval that does not lift whole to the two published 200 kb windows that begin and end on its boundaries -- every borrowed window is checked for membership in the published window set, read from the sliding-window BEDs that are therefore a second input to this product (the 10 kb grid is necessary for membership and not sufficient), so the composed coordinate is the one the windows table carries for that window, while lifting the 1 bp endpoints instead would displace a boundary into a paralogous copy of the segmental duplication that flanks it. The GRCh37 pair is kept and GRCh38 is left NULL where a boundary window is itself split, partially deleted, or lifts to a length the shared ±10% filter rejects (6 of the 163 segments on the current chain); its HPO, credible-interval and gene lists stay `;`-joined for the BigQuery loader to split into arrays. `--product windows` streams the same paper's 108 sliding-window BEDs (54 phenotypes x DEL/DUP, the same 267,237 GRCh37 200 kb windows in each) into one long table, lifting the window set once with that same procedure and asserting it reproduces the measurement's 262,357 lifted / 4,880 dropped windows exactly; unlike the gene product it drops the rows whose meta-analysis produced no estimate, because the window set is fixed and a missing row already means "no estimate". The lifted windows are not a regular grid -- adjacent pairs can reorder in GRCh38 and not all keep a 200 kb width -- so `window_start_grch37`/`window_end_grch37` ride along and no consumer should derive a step or a width from the GRCh38 pair; see [docs/rcnv-sliding-windows.md](docs/rcnv-sliding-windows.md) for the liftOver measurement itself and the width distribution. Neither `scores` nor `genes` adds coordinates: the source has none, so those outputs are build-independent; none of the four needs a tabix index. `--product` is there because the same Zenodo record ships one more product (the gene feature matrix); that one is not implemented. See [docs/rcnv-dosage-sensitivity.md](docs/rcnv-dosage-sensitivity.md).
- curated gene-disease associations: `munge_gene_disease.{py,sh}` with `--product gencc` or `--product monarch` downloads the source and writes a TSV named after its version, so publishing never overwrites the file the deployments are still reading. GenCC is passed through unchanged — the round trip exists to fail here rather than in the API when the export stops parsing. Monarch concatenates the KG's causal and non-causal gene-disease exports, keeps `predicate` so the two stay distinguishable downstream, drops the non-gene, non-human and negated rows, and collapses the verbatim duplicates the sources carry. Neither output has coordinates, so neither is bgzipped or tabixed: the results-api reads the plain TSV once at startup. See [docs/gene-disease-associations.md](docs/gene-disease-associations.md).
- rsID lookup: `build_gnomad_rsid_index.py` turns the gnomAD v4 sites file into the `rs`-pseudo-contig tabix file that results-api's `/rsid/variants` (the `lookup_variants_by_rsid` tool) queries by rs number. It emits one row per (rsID, allele) — one rsID can name several alt alleles, and the previous build's one-row-per-rsID layout resolved a multi-allelic rsID to whichever alt sorted first, often an allele with no carriers — and drops an allele filtered `AC0` only when the same rsID has an observed one. Rows must be sorted by rs number for the index, and the input is ~800M rows, so it partitions into rs-number buckets on disk and merges them; peak disk is about the gzip size of the exploded rows plus the output. `--upload` copies the result to a bucket directory, under a new name: results-api treats the served file as immutable and caches its index per path.
- gnomAD annotation: `build_gnomad_annotation.py` streams the position-sorted gnomAD genomes+exomes sites file (or one genomes and one exomes file, merged by position) into two bgzip + tabix outputs: the sites file deduplicated to one row per variant, and a `chr pos ref alt most_severe gene_most_severe` consequence file for stream-merging against credible-set files. The header docstring states which row of a genome/exome pair survives, the chromosome coding and the inputs it refuses. Tests: `cd scripts && python -m pytest test_build_gnomad_annotation.py`.
- expression: `munge_gtex.py` (GTEx v10 median TPM, written both wide and one row per gene and tissue) and `munge_hpa.py` (HPA immunohistochemistry).
- gene-indexed QTL credible sets: `create_gene_indexed_qtl_file.py` for datasets whose QTL trait is a gene (the caQTL and peak-to-gene variants have their own sections above).

- consequence stamping: `annotate_consequence.py` overwrites `most_severe` and `gene_most_severe` of an already served credible set file from the consequence file that `build_gnomad_annotation.py` writes, leaving every other byte of every row as it was, so an annotation refresh needs no re-munge. `--mode merge` is for the variant-sorted combined and per-trait files, `--mode lookup` for the gene-indexed `*.qtl.tsv.gz` copies, `--clear` writes `NA` to both columns, and `--verify ORIGINAL STAMPED` checks that nothing else changed. A bgzip input gets a bgzip output indexed with the settings read from the input's own index. Its header docstring has the guards and how the consequence file is read.
- stamping a whole resource: `annotate_resource.sh` runs `annotate_consequence.py` over every object under a served resource's prefix (combined file, gene-indexed QTL copies, per-trait files), regenerates the credible set stats from the stamped rows, and writes the result under a new prefix, never over a served path. It refuses an object it cannot classify, and checks before uploading that the original stats are reproduced from the original rows and that nothing outside the two annotation columns changed. `--dry-run` prints the classification; the header docstring of `annotate_resource.py` is the reference.
- BRaVa phenotype metadata: `brava_phenotypes.py` builds the pheweb-shaped JSON (`phenocode`, `phenostring`, `category`, `num_cases`/`num_controls` or `num_samples`) for the BRaVa exome-wide rare variant meta-analysis, one entry per (phenotype, stratum) that has a gene burden result file — the listing of `gs://daly-genetics-results/raw/brava/gene/` decides which, because the strata a phenotype has depend on its case counts. Descriptions, sex and per-biobank counts come from the preprint's supplementary workbook, downloaded by default. Every count is summed from the per-biobank tables so a phenotype's strata add up to its meta-analysis; Table S4 reproducing Table S6 exactly is the gate that proves the summing right, and Table S5's repeated rows are dropped even though that puts the quantitative sample sizes up to 30% below the published Table S7 totals — the deduplicated Height EUR sum is 710,271 against a largest per-variant NS of 710,270, and the raw sum is 920,670. `--check` prints a sample of the output beside those per-variant NS/NC maxima and never fails on it; `--stage` uploads to `gs://daly-genetics-results/mapping_files/brava_pheno.json`.
- gene and trait metadata helpers: `gencode_to_gene_pos_tsv.py` (GENCODE GFF3 to a gene position TSV, given a release URL), `gencode_to_exon_tsv_gff.py` (the same GFF3 to one row per exon — exon and coding bounds, plus the Ensembl-canonical and MANE Select flags — which the API filters to canonical to draw gene models; the chromosome encoding matches the gene position TSV so the two join on `chrom`), `create_gene_name_mapping_across_gencode_versions.py` (the cross-version `ensg -> name` table the API reads; its version list must match `gencode_versions` in the API's `genes.py`, and the output file name carries those versions) and `kanta_metadata_to_json.py`.

The open chromatin, variant effect, MPRA, HLA, rCNV and gene-disease wrappers write locally by default and only upload to the two profile buckets when given `--stage` (gene-disease takes `GCS_DESTS` to narrow that to one bucket). The expression, gene mapping and gene-indexed QTL scripts have their input paths hardcoded at the top of the script.

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

For BRaVa `trait_original` is the phenocode, which carries the meta-analysis stratum:
`AFib.tsv.gz` is the cross-ancestry meta-analysis and `AFib|EUR.tsv.gz` one ancestry
stratum of the same phenotype, so the file name — and the trait the API is asked for —
contains a `|`.

Both files are indexed on a POINT — `-s5 -b6 -e6` takes begin and end from
`gene_start_pos` — so a gene lookup in the API only hits when the coordinates in the
file come from exactly the GENCODE version configured for that dataset. Pass
`--gencode` the version the API declares (genebass v35, SCHEMA2/BipEx2/BRaVa v39,
IBD exome v45); a mismatch returns nothing rather than erroring.

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
aaf                 alternative allele frequency, joined from the variant annotation file;
                    NA for Open Targets, whose release has no per-variant frequency
                    (the fine-mapping results only have MAF, so the WDL pipeline leaves MAF
                    here when run without a variant annotation file)
most_severe         most severe variant consequence (VEP)
gene_most_severe    gene of most severe consequence
```

`trait` is the name a trait is shown and selected by, `trait_original` the identifier it has at its source. What that means depends on the resource:

- FinnGen-format fine-mapping (the WDL pipeline): `trait_original` is the phenotype code. With a phenotype metadata file `trait` is the phenostring with spaces replaced by underscores, falling back to the code for a phenotype the file does not name; without one both columns are the code.
- Open Targets: `trait_original` is the study accession and `trait` is the study's `traitFromSource` followed by `_(<accession>)`, e.g. `Type_2_diabetes_(GCST004602)`. Many accessions share one trait name, so the suffix is what makes `trait` identify a study. Whitespace in the name becomes underscores and double quotes become single quotes; other punctuation and non-ASCII letters are kept as the source has them.
- QTL results: `trait` contains the gene name while `trait_original` contains the original QTL trait name depending on the dataset, e.g. ENSG gene id.

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
which variants are excluded, and for the deCODE pQTL run, the one dataset whose input is a
QTL study rather than a GWAS.

## variant annotation

For FinnGen data, `most_severe` and `gene_most_severe` come from the FinnGen variant annotation joined in the munge, so a variant outside the FinnGen imputation panel has `NA` in both. For the resources that are not FinnGen data the munge writes `NA`, and `scripts/annotate_resource.sh` stamps both columns afterwards from the gnomAD consequence file that `build_gnomad_annotation.py` writes; a variant gnomAD does not hold stays `NA`.
