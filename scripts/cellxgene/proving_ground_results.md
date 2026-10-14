# CELLxGENE kill-criterion measurement (i8om.2) — results, 2026-10-13

Throwaway measurement, not the production worker. Code: `proving_ground.py` (download /
process / validate / gencode) and `cost_model.py` in this directory. Raw outputs from the run
are under `/mnt/disks/data/cellxgene-pg/` (`out/<dataset_id>/`, `logs/`, `validate*.txt`,
`marker_panel.*.tsv`, `gencode.txt`, `cost.txt`); downloaded H5ADs were deleted after each run.

## Premise numbers against the epic thresholds

| premise | threshold | measured | verdict |
|---|---|---|---|
| VALIDITY, mean log CP10k | pass >= 90 %, descope < 80 % | **33 / 37 = 89.2 %** of covered markers | grey zone (just under pass) |
| VALIDITY, donor pseudobulk CP10k | same | **25 / 37 = 67.6 %** | below descope line |
| VALIDITY, post-hoc sensitivity (myofibroblast rolled into fibroblast) | — | mean log 36/37 = 97.3 %, pseudobulk 29/37 = 78.4 % | informational only; a roll-up decision, not a measurement |
| COST, largest single dataset on one worker | <= 8 h, < 16 GB | slowest: MSSM_Cohort (36 GB, 15.6 B nnz) **33.5 min** incl. download, **5.4 GB** peak RSS; largest file (SEA-AD 53 GB) **21.3 min**, **1.9 GB** | pass by ~15x |
| SERVING | <= 1 GB scan, < 1 s tabix | not measured here | — |

## Setup

- Host: this dev VM, e2-standard-4 (AMD EPYC 7B12), one Python process per dataset (CPU-bound at ~97 % of one core).
  Every process run was under `systemd-run --user --scope -p MemoryMax=14G -p MemorySwapMax=0`.
- Download: 16 parallel 64 MB byte-range GETs from `datasets.cellxgene.cziscience.com` into `/mnt/disks/data/cellxgene-pg/dl`, one file on disk at a time.
- Reading: h5py only; `raw/X` if present else `X`; CSR rows walked in ~25 M-nnz chunks via `indptr`; chunks stay sparse
  (group sums are sparse indicator-matrix products). Non-integer counts refuse the dataset (none did). `is_primary_data` cells only.
- Groups: (tissue, tissue_type, cell_type, disease, assay, suspension_type). `tissue_type` was added to the five
  requested keys so that organoid and cell-culture cells never pool with tissue under the same UBERON id.
- Per group: n_cells, n_expressing, mean log1p(CP10k) over all cells (zeros included); per (group, donor):
  summed counts and library size, giving donor pseudobulk CP10k averaged over donors. The (group, donor) x gene
  accumulator is kept dense in float32 when it fits in 3 GB (it did for all eight); otherwise the script does a
  first pass over `data` for library sizes and accumulates the pseudobulk pre-weighted (implemented, not exercised).

## Per-dataset measurements

| tag | dataset_id | H5AD GB | download s | MB/s | obs/var s | main pass s | write s | peak RSS MB | nnz (B) | M nnz/s | cells | primary kept | groups | group x donor | genes | matrix |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| largest file (SEA-AD DLPFC) | 6f7fd0f1 | 53.19 | 245 | 217 | 9.6 | 1022.3 | 1.9 | 1854 | 7.924 | 7.8 | 1,395,601 | 1,395,601 | 72 | 1,869 | 35,483 | raw/X |
| heaviest nnz (MSSM_Cohort) | 37a17b78 | 36.09 | 140 | 258 | 22.2 | 1827.9 | 17.3 | 5355 | 15.599 | 8.5 | 4,140,453 | 4,140,453 | 1,119 | 21,859 | 34,176 | X (raw counts) |
| median size (Helmsley gut, colon) | 944dedde | 0.33 | 1 | 386 | 0.1 | 5.7 | 0.2 | 1386 | 0.064 | 11.2 | 25,641 | 25,641 | 13 | 128 | 32,344 | raw/X |
| blood (NYGC PBMC) | ed5d841d | 2.61 | 8 | 335 | 0.7 | 30.3 | 0.4 | 1480 | 0.353 | 11.6 | 161,764 | 161,764 | 30 | 240 | 20,264 | raw/X |
| lung (Krasnow HLCA 10x) | 8c42cfd0 | 0.60 | 3 | 214 | 0.3 | 10.4 | 0.4 | 1593 | 0.126 | 12.1 | 65,662 | 65,662 | 65 | 142 | 24,769 | raw/X |
| liver | e84f2780 | 1.03 | 4 | 296 | 0.7 | 19.2 | 0.4 | 1650 | 0.216 | 11.2 | 167,598 | 167,598 | 38 | 181 | 30,172 | raw/X |
| pancreas | 31f657dc | 0.96 | 62 | 15.6 | 0.4 | 8.8 | 0.2 | 1506 | 0.159 | 18.0 | 121,916 | 79,367 | 33 | 151 | 32,356 | raw/X |
| Visium (breast) | 60d5175f | 0.01 | 1 | 8 | 0.1 | 0.2 | 0.1 | 230 | 0.001 | 4.7 | 4,992 | 4,992 | 7 | 7 | 35,477 | X |

nnz/s for the pancreas row is over all cells; only 65 % were primary and the rest skip the arithmetic.
The listing's `cell_count x mean_genes_per_cell` matched measured nnz closely (liver 0.22 B vs 0.216 B; MSSM 15.6 B vs 15.599 B),
so it is a usable planner estimate. Output per dataset (parquet): 0.2-10 MB for the small ones, 30 MB SEA-AD, 276 MB MSSM
(24.2 M group x gene rows).

## VALIDITY — marker panel

**Roll-up level.** Each observed cell type is mapped to its most specific ancestor (is_a closure, CL release 2026-06-08,
`cl-basic.obo`) within a fixed target set: the panel's 24 expected terms plus five broad competitors (dendritic cell,
epithelial cell, granulocyte, pericyte, glial cell). Types with no such ancestor stay as themselves. A marker passes when
the top-ranked rolled type is the expected term or a descendant of it. Rolled types need >= 100 normal cells.

**Pool.** Normal disease (PATO:0000461), `tissue_type` = tissue, cell type not `unknown`, all eight datasets:
2,346,374 cells, 34 rolled types. A gene absent from a dataset's `var` counts as *unmeasured* there, not as zero
(the pancreas dataset has no PRSS1 or PRSS2 in its var at all). Mean log CP10k is pooled cell-weighted over the
groups of datasets that measured the gene; the pseudobulk sums counts and library per (rolled type, dataset, donor),
converts to CP10k and averages over donors.

**Panel (40 markers, expected term):** CD3E, CD3D, CD8A -> T cell CL:0000084; MS4A1, CD79A, CD19 -> B cell CL:0000236;
JCHAIN -> plasma cell CL:0000786; CD14, LYZ -> monocyte CL:0000576; NKG7, GNLY -> NK cell CL:0000623; C1QA -> macrophage
CL:0000235; CLEC4C -> pDC CL:0000784; TPSAB1 -> mast cell CL:0000097; HBB -> erythrocyte CL:0000232; PPBP -> platelet
CL:0000233; ALB, APOA1 -> hepatocyte CL:0000182; SFTPC -> AT2 CL:0002063; AGER -> AT1 CL:0002062; INS -> type B pancreatic
cell CL:0000169; GCG -> pancreatic A cell CL:0000171; PRSS1, CPA1 -> pancreatic acinar cell CL:0002064; MBP, PLP1, MOG ->
oligodendrocyte CL:0000128; GFAP, AQP4 -> astrocyte CL:0000127; SNAP25, SYT1 -> neuron CL:0000540; PECAM1, VWF, CDH5 ->
endothelial cell CL:0000115; COL1A1, DCN, LUM -> fibroblast CL:0000057; KRT5 -> basal cell CL:0000646; MYH11 -> smooth
muscle cell CL:0000192; MUC2 -> goblet cell CL:0000160.

**Coverage: 37 / 40.** Not testable with these datasets: TPSAB1 (no cells labelled mast cell; the lung dataset's mast
cells are labelled `myeloid leukocyte`), HBB (no erythrocytes in any normal tissue sample), PRSS1 (not in the only pancreas
dataset's var).

| | mean log CP10k | donor pseudobulk |
|---|---|---|
| pass / covered | 33 / 37 (89.2 %) | 25 / 37 (67.6 %) |
| failures | AQP4 (top: AT1 cell); COL1A1, DCN, LUM (top: myofibroblast) | the four mean-log failures plus CD3E, CD3D, NKG7, GNLY (top: `lymphocyte`, 105 cells / 2 donors); CD8A (top: hepatocyte, 7 donors); LYZ (dendritic cell); MYH11 (myofibroblast); MUC2 (epithelial cell) |

What the failures are:
- **Myofibroblast** (CL:0000186) is not an is_a descendant of fibroblast in CL, so it competes with fibroblasts and wins
  three fibroblast markers. Pinning it under fibroblast (a roll-up-table decision for i8om.4) lifts mean log to 36/37
  (97.3 %). Shown only as a sensitivity; it was decided after seeing the result.
- **AQP4** genuinely is high in lung AT1 cells; this is a panel choice more than a pipeline error.
- **Pseudobulk** is dominated by tiny, coarsely labelled types: averaging over donors gives a 2-donor `lymphocyte` group
  (labels coarser than T/NK from one dataset) equal weight to 39 T-cell donors, and ambient RNA in hepatocytes (CD8A) shows
  up in summed counts but not in the per-cell log mean. A minimum-donor rule or a coarse-label exclusion would change this;
  neither was applied here.

## GENCODE v48 -> v49

Union of Ensembl ids in the eight files' var (version stripped): **44,416; 29 have no match in GENCODE v49
(0.07 %)**. Per file 1-28. 28 of the 29 are detected in at least one group and 18 reach > 1 % expressing in some group,
including an id labelled MATR3 (expressed in every cell of some groups; MATR3 had two gene ids) and MAFIP. The rest are
unnamed `ENSG000002[5-7]...` loci. All var entries were ENSG ids (no non-Ensembl features).

## COST ESTIMATE — initial run

Work list from the curation-API listing (cached, 2026-10-13): human, `primary_cell_count > 0`, snmC-Seq2/methylation
excluded, spatial included: **1,047 datasets, 1.63 TB of H5AD, 147 M cells**. (The epic's "1,633 datasets / 2.6 TB" counts all
human datasets, including ones with no primary cells.)

Model per dataset: `boot + filesize / download_rate + cells x 4.2 us + nnz x rate + write`, with
- rate 12.1 M nnz/s for files <= 2 GB (median measured) and **7.8 M nnz/s for > 2 GB** (slowest measured; large files ran ~30 % slower per nnz);
- download 275 MB/s median (also run at 217 MB/s, the slowest multi-GB download: +0.5 worker-h);
- boot **150 s assumed** (Batch VM provisioning and image pull are not measured here; i8om.1 measures them);
- nnz from the listing (`cell_count x mean_genes_per_cell`); pd-balanced disk sized `max(30, 1.1 x file + 20)` GiB for the job's duration.

Prices read from the Cloud Billing Catalog API (`cloudbilling.googleapis.com/v1/services/*/skus`, USD, us-central1):

| SKU | description | price |
|---|---|---|
| CF4E-A0C7-E3BF / F449-33EC-A5EF | E2 core / RAM, on-demand | $0.02181159 /vCPU-h, $0.00292353 /GiB-h |
| F179-E1EA-D97A / 9B1F-1E62-4061 | E2 core / RAM, Spot | $0.01309 /vCPU-h, $0.001754 /GiB-h |
| BB77-5FDA-69D9 / 5B01-D157-A097 | N2 core / RAM, on-demand | $0.031611 /vCPU-h, $0.004237 /GiB-h |
| 2856-91B1-406A / 6276-ABCE-E764 | N2 core / RAM, Spot | $0.01896 /vCPU-h, $0.002542 /GiB-h |
| 6AE1-525F-8B80 | Balanced PD capacity | $0.10 /GiB-month |
| D716-3E59-667B | Coldline Storage Iowa | $0.004 /GiB-month |
| E5F0-6A5D-7BAD | Standard Storage US Regional | $0.02 /GiB-month (first 5 GiB free) |
| 0018-A5A0-9D6D | BigQuery active logical storage (us-central1) | $0.023 /GiB-month (first 10 GiB free) |

| class (by H5AD size) | n | GB | worker-h | Spot $ | on-demand $ | pd-balanced $ |
|---|---|---|---|---|---|---|
| e2-standard-2 (<= 2 GB) | 892 | 351 | 43.3 | $1.74 | $2.90 | $0.18 |
| e2-highmem-4 (2-20 GB) | 141 | 804 | 13.5 | $1.47 | $2.44 | $0.06 |
| n2-highmem-8 (> 20 GB) | 14 | 480 | 4.8 | $1.52 | $2.54 | $0.04 |
| **total** | **1,047** | **1,634** | **61.6** | **$4.73** | **$7.88** | **$0.27** |

(At the slowest measured download rate the total is 62.1 worker-h, $4.80 Spot / $8.00 on-demand; full output in
`/mnt/disks/data/cellxgene-pg/cost.txt`.)

- **Critical path** of the fan-out (slowest single dataset): MSSM_Cohort, **measured 33.5 min** (download 2.3 + processing
  31.1) plus boot; modelled 0.65 h. Next: a 27.7 GB / 12.5 B-nnz dataset at ~0.53 h modelled. Nothing approaches 8 h.
- **Download/egress:** ingress into GCP is free; workers have external IPs, so there is no Cloud NAT charge.
- **Storage per month:** Coldline raw H5ADs 1,522 GiB -> **$6.09**; aggregates (parquet, extrapolated per cell) ~6 GiB
  Standard -> **$0.03**; BigQuery **85 GiB -> $1.74** (3.1 group x gene rows per cell, median measured, x 147 M cells = 0.46 B
  rows at ~200 logical bytes/row) up to **1.4 TiB -> $33** (group count from the listing's cardinality product). The low
  end is the better estimate; the high end assumes every listed cell type occurs in every listed tissue x disease x assay.
  BigQuery batch load jobs are not billed (not a catalog read).
- **Weekly refresh:** the last 13 weeks changed 58 datasets / 181 GB -> **~4.5 datasets/week, ~$0.04 Spot ($0.06
  on-demand) compute per week**, plus Coldline growth for new versions (~14 GB/week, ~$0.06/month added per week).

**Total initial run: about $5 Spot (about $8 on-demand) in compute and disk, plus ~$8-40/month storage.**

## What surprised me

1. **The largest file is not the slowest.** Time follows nnz, not bytes: the 53 GB SEA-AD file has 7.9 B nnz and took 21 min;
   MSSM (36 GB, 15.6 B nnz) took 33.5 min. The 8 h threshold is missed by ~15x even single-threaded.
2. **Memory is tiny.** Peak RSS 5.4 GB on the heaviest dataset (3 GB of it the dense group x donor accumulator), under 2 GB
   everywhere else. Processing is single-threaded, so the planned e2-highmem-4 and n2-highmem-8 classes buy idle cores
   and unused RAM; every dataset measured fits an e2-standard-2 (8 GB) or e2-highmem-2.
3. **Gene panels differ between datasets** (the pancreas dataset lacks PRSS1 and PRSS2 entirely). A cross-dataset summary
   that treats a missing gene as zero is wrong; the publish step needs a per-dataset measured-gene mask.
4. **Donor pseudobulk is the fragile metric** (67.6 %) because of equal donor weighting over coarse or tiny label groups.
   Mean log CP10k is robust (89.2 %, 97.3 % with one roll-up pin).
5. One download ran at 15.6 MB/s (pancreas, 1 GB) while all others ran at 200-385 MB/s. A worker needs retry/timeout
   handling, but with 16 ranges it never failed.
6. MSSM ships raw integer counts in `X` with no `raw/X`; the integer check is what lets it through, so the check is load-bearing.
