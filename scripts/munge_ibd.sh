time for pheno in cd ibd uc; do
    /usr/bin/time -v python3 scripts/munge_ibd.py \
		  --phenotype ${pheno} \
		  --input-dir /mnt/disks/data/IBD/lustre/scratch124/humgen/projects_v2/ibdgwas/IIBDGC/post_imputation/2022/analysis/metaanalysis/${pheno} \
		  --gnomad-af-plot \
		  --gnomad-filtered /mnt/disks/data/IBD/lustre/scratch124/humgen/projects_v2/ibdgwas/IIBDGC/post_imputation/2022/analysis/metaanalysis/ibd/ibd_meta.munged.gnomad_filtered.tsv
done
