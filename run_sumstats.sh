BUCKET=gs://finngen-commons/results_api_data/sumstats
IBD_INPUT=/mnt/disks/data/IBD/lustre/scratch124/humgen/projects_v2/ibdgwas/IIBDGC/post_imputation/2022/analysis/metaanalysis

# IBD meta-analysis (3 phenotypes)
for pheno in ibd cd uc; do
  python3 scripts/munge_ibd.py \
    --phenotype ${pheno} \
    --input-dir ${IBD_INPUT}/${pheno} \
    --output ${BUCKET}/IBD/${pheno}_meta.munged.tsv.gz
done

# PGC schizophrenia
python3 scripts/munge_pgc.py \
  --input /mnt/disks/data/PGC/daner_PGC_SCZ_w3_90_0418b.gz \
  --gnomad-af-plot \
  --gnomad-filtered /mnt/disks/data/PGC/daner_PGC_SCZ_w3_90_0418b.gnomad_filtered.tsv \
  --output ${BUCKET}/PGC/daner_PGC_SCZ_w3_90_0418b.munged.tsv.gz

# BIP 2024
python3 scripts/munge_bip2024.py \
  --input /mnt/disks/data/BIP/bip2024_multianc_no23andMe.gz \
  --gnomad-af-plot \
  --gnomad-filtered /mnt/disks/data/BIP/bip2024_multianc_no23andMe.gnomad_filtered.tsv \
  --output ${BUCKET}/BIP/bip2024_multianc_no23andMe.munged.tsv.gz

# GP2 Parkinson's
python3 scripts/munge_gp2.py \
  --input /mnt/disks/data/Parkinsons/GP2_euro_ancestry_meta_analysis_2024/GP2_ALL_EUR_ALL_DATASET_HG38_12162024.txt.gz \
  --gnomad-filtered /mnt/disks/data/Parkinsons/GP2_euro_ancestry_meta_analysis_2024/GP2_ALL_EUR_ALL_DATASET_HG38_12162024.txt.gnomad_filtered.tsv \
  --output ${BUCKET}/GP2/GP2_ALL_EUR_ALL_DATASET_HG38_12162024.txt.munged.tsv.gz

# AIH autoimmune hypothyroidism meta-analysis (3 phenotypes, one file each).
# gnomAD is streamed once for the union of the three files' positions, then reused.
AIH_INPUT="/mnt/disks/data/aih_meta.txt.gz /mnt/disks/data/aitt1_meta.txt.gz /mnt/disks/data/aitt2_meta.txt.gz"
python3 scripts/munge_aih.py --input ${AIH_INPUT} --gnomad-filter-only --jobs 4
python3 scripts/munge_aih.py \
  --input ${AIH_INPUT} \
  --gnomad-filtered /mnt/disks/data/ai_meta.gnomad_filtered.tsv.gz \
  --gnomad-af-plot \
  --plot-dir /mnt/disks/data \
  --output-dir ${BUCKET}/AIH
