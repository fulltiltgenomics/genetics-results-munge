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

time python3 "${SCRIPT_DIR}/munge_coloc_file.py" $filename $eqtl_cat_metadata_filename $output_file

# check that there is no missing data
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
}' $output_file

time sort -T . -k17,17g -k18,18g -k19,19g -k5,5 -k6,6 $output_file | bgzip -@4 > $output_file.gz
time tabix -f -s 17 -b 18 -e 19 $output_file.gz

# map trait (gene) names of QTL traits
bname=$(basename $output_file .tsv)
time python3 "${SCRIPT_DIR}/munge_coloc_file_map_traits.py" $output_file.gz $eqtl_cat_metadata_filename $bname
# combine non eQTL cat data with eQTL cat data
time sort -m -T . -k19,19g -k20,20g -k21,21g -k5,5 -k7,7 $bname.*.tsv | uniq | bgzip -@4 > $bname.tsv.gz && tabix -s 19 -b 20 -e 21 $bname.tsv.gz
