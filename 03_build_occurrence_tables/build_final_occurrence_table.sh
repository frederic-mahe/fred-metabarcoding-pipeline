#!/bin/bash
cd "${PWD}"

export LC_ALL=C

## ------------------------------------------------------------ define variables
declare -r OTU_TABLE="${1}"
declare -r SRC="${HOME}/save/src"
declare -r OTU_TABLE_UPDATER="OTU_table_updater.py"


## ------------------------------------------------------------------- functions

update_taxonomic_assignments() {
    python3 \
        "${SRC}/${OTU_TABLE_UPDATER}" \
        --old_otu_table "${OTU_TABLE}" \
        --new_taxonomy "${OTU_TABLE%.*}.results" \
        --new_otu_table "${NEW_TABLE}"
}

fix_cluster_sorting() {
    # sort by decreasing abundance value
    (head -n 1 "${NEW_TABLE}"
     tail -n +2 "${NEW_TABLE}" | \
         sort -k2,2nr | \
         nl -n'ln' -w1 | \
         cut --complement -f 2
    ) > "${OTU_TABLE}2"
}

clean_up() {
    chmod go+r,g-w "${OTU_TABLE}2"
    rm "${NEW_TABLE}"
}

## ------------------------------------------------------------------------ main
echo "build final OTU table..."
NEW_TABLE=$(mktemp)

update_taxonomic_assignments
fix_cluster_sorting
clean_up

exit 0
