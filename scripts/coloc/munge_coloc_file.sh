#!/bin/bash

set -euxo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

if [ $# -ne 3 ]; then
    echo "Usage: $0 <filename> <eqtl_cat_metadata_filename> <output_filename>"
    exit 1
fi
filename=$1
eqtl_cat_metadata_filename=$2
output_file=$3

intermediate="${output_file%.tsv}.intermediate.tsv"

# pass 1: filter + reshape (dataset/data_type/cell_type/hit stats), traits stay raw
time python3 "${SCRIPT_DIR}/munge_coloc_file.py" "$filename" "$eqtl_cat_metadata_filename" "$intermediate"

# check that there is no missing data in required columns of the intermediate
time awk '
BEGIN {
    non_na_cols="#dataset1,dataset2,data_type1,data_type2,trait1,trait2,cs1_id,cs2_id,hit1,hit2,chr,region_start_min,region_end_max,PP.H0.abf,PP.H1.abf,PP.H2.abf,PP.H3.abf,PP.H4.abf,nsnps,nsnps1,nsnps2,cs1_log10bf,cs2_log10bf,cs1_size,cs2_size,cs_overlap,topInOverlap"
    n=split(non_na_cols, a, ",")
    has_error=0
    num_errors=0
}
NR==1 {
    for (i=1; i<=NF; i++) {
        col_map[$i] = i
    }
}
{
    for (i=1; i<=n; i++) {
        col_num = col_map[a[i]]
        if ($col_num == "NA") {
            if (num_errors < 10) {
                print "column "a[i]" has NA on row "NR > "/dev/stderr"
            } else if (num_errors == 10) {
                print "..." > "/dev/stderr"
            }
            has_error=1
            num_errors++
        }
    }
}
END {
    exit has_error
}' "$intermediate"

# pass 2: split trait into trait/trait_original and map eQTL Catalogue trait2 to gene names
time python3 "${SCRIPT_DIR}/munge_coloc_file_map_traits.py" "$intermediate" "$eqtl_cat_metadata_filename" "$output_file"

# final schema column positions: 19=chr, 20=region_start_min, 21=region_end_max, 5=trait1, 7=trait2
time sort -T . -k19,19g -k20,20g -k21,21g -k5,5 -k7,7 "$output_file" | bgzip -@4 > "$output_file.gz"
time tabix -f -s 19 -b 20 -e 21 "$output_file.gz"

rm -f "$intermediate"
