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

time python3 "${SCRIPT_DIR}/munge_coloc_credset_file.py" $filename $eqtl_cat_metadata_filename $output_file

# check that there is no missing data
time awk '
BEGIN {
    cell_type_col="cell_type"
    data_type_col="data_type"
    non_na_cols="#dataset,data_type,trait,trait_original,chr,pos,ref,alt,pip,cs_id"
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
    for (i=1;i<=n;i++) {
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
    if ($col_map[cell_type_col] == "NA" && $col_map[data_type_col] != "GWAS") {
        if (num_errors < 10) {
            print "missing cell type on row "NR > "/dev/stderr"
        } else if (num_errors == 10) {
            print "..." > "/dev/stderr"
        }
        has_error=1
        num_errors++
    }
}
END {
    exit has_error
}' $output_file

time sort -T . -k6,6g -k7,7g -k8,8 -k9,9 -k3,3 $output_file | bgzip -@4 > $output_file.gz
time tabix -s 6 -b 7 -e 7 $output_file.gz
