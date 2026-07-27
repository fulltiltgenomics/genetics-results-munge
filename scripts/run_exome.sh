GCS=gs://finngen-commons/results_api_data/exome_results5

# --- SCHEMA gene burden ---
python3 munge_schema.py \
	--input /mnt/disks/data/SCHEMA_gene_results.tsv.bgz \
	--output-dir $GCS/schema

# --- SCHEMA variant results ---
python3 munge_schema_variants.py \
	--input /mnt/disks/data/SCHEMA_variant_results.tsv.bgz \
	--output $GCS/schema/SCHEMA_variant_results.munged.tsv.gz

# --- SCHEMA2 gene burden ---
python3 munge_schema2.py \
	--input /mnt/disks/data/SCHEMA/SCHEMA2_gene_results.tsv.bgz \
	--output-dir $GCS/schema

# --- SCHEMA2 variant results ---
python3 munge_schema2_variants.py \
	--input /mnt/disks/data/SCHEMA/SCHEMA2_variant_results.tsv.bgz \
	--output $GCS/schema/SCHEMA2_variant_results.munged.tsv.gz

# --- BipEx gene burden ---
python3 munge_bipex.py \
	--input /mnt/disks/data/BibEx/BipEx2_gene_results.tsv.bgz \
	--output-dir $GCS/bipex

# --- IBD exome (old dataset) gene + variant ---
python3 munge_ibd_exome.py \
	--gene-input /mnt/disks/data/IBD/IBD_gene_results.tsv.bgz \
	--variant-input /mnt/disks/data/IBD/IBD_variant_results.tsv.bgz \
	--output-dir $GCS/ibd

# --- IBD exome 2026 (supp tables) gene burden ---
GCS=gs://finngen-commons/results_api_data/exome_results3
python3 munge_ibd_supp_burden.py \
	--input /mnt/disks/data/IBD/Supplementary_Tables_-_ST7_Burden_test.tsv \
	--n-cases 86213 44131 32748 \
	--n-controls 478363 \
	--output-dir $GCS/ibd

# --- IBD exome 2026 (supp tables) variant results ---
python3 munge_ibd_supp_variants.py \
	--input /mnt/disks/data/IBD/Supplementary_Tables_-_ST4.83_Tier1_8_LD_unresolved.tsv \
	--input-extra /mnt/disks/data/IBD/Supplementary_Tables_-_ST5.17nsny_p3e7_5e6.tsv \
	--n-cases 86213 44131 32748 \
	--n-controls 478363 \
	--output-dir $GCS/ibd
