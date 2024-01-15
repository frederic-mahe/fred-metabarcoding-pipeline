#!/bin/bash
cd "${PWD}"

export LC_ALL=C

## ------------------------------------------------------------ define variables
declare -r PROJECT="${1}"
declare -r DATA_FOLDER="${2}"  # can be a list of folders (space-separated)
declare -r MARKER="${3}"
declare -ri THREADS="${4:-4}"
declare -r SRC="${HOME}/save/src"
declare -r SWARM="${SRC}/swarm/bin/swarm"  # swarm 3.0 or more recent
declare -r VSEARCH="${SRC}/vsearch/bin/vsearch"  # vsearch 2.21 or more recent
declare -r MUMU="${SRC}/mumu/mumu"  # mumu 0.0.1 or more recent
declare -ri RESOLUTION=1
declare -ri FILTER=2
declare -r OTU_CLEAVER="OTU_cleaver.py"
declare -r OTU_TABLE_BUILDER="OTU_contingency_table_filtered.py"
declare -r MERGE_SUBSTRINGS="merge_sub_superstring_OTUs_with_larger_OTUs.py"
declare -r REBUILD_TABLE_AFTER_MUMU="rebuild_table_after_mumu.py"
declare -r STAMPA="stampa.sh"
declare -r OTU_TABLE_UPDATER="OTU_table_updater.py"

## variables and file names
N_SAMPLES=$(find ${DATA_FOLDER} -name "*.fas" \
                 -type f ! -empty -print0 | tr -d -c '\0' | wc -m)
FINAL_FASTA="${PROJECT}_${N_SAMPLES}_samples.fas"
QUALITY_FILE="${FINAL_FASTA%.*}.qual"
DISTRIBUTION_FILE="${FINAL_FASTA%.*}.distr"
POTENTIAL_SUB_SEEDS="${FINAL_FASTA%.*}_per_sample_OTUs.stats"
LOG="${FINAL_FASTA%.*}.log"
OUTPUT_SWARMS="${FINAL_FASTA%.*}_${RESOLUTION}f.swarms"
OUTPUT_LOG="${FINAL_FASTA%.*}_${RESOLUTION}f.log"
OUTPUT_STATS="${FINAL_FASTA%.*}_${RESOLUTION}f.stats"
OUTPUT_STRUCT="${FINAL_FASTA%.*}_${RESOLUTION}f.struct"
OUTPUT_REPRESENTATIVES="${FINAL_FASTA%.*}_${RESOLUTION}f_representatives.fas"
TAXONOMIC_ASSIGNMENTS="${OUTPUT_REPRESENTATIVES%.*}.results"
UCHIME_RESULTS="${OUTPUT_REPRESENTATIVES%.*}.uchime"
UCHIME_LOG="${OUTPUT_REPRESENTATIVES%.*}.log"
OTU_TABLE="${FINAL_FASTA%.*}.OTU.filtered.cleaved.table"
OUTPUT_TABLE="${OTU_TABLE%.*}.nosubstringOTUs.table"

## ------------------------------------------------------------------- functions

# add a function to check MARKER values (grep stampa.sh, check if
# marker is valid)

build_expected_error_file() {
    find ${DATA_FOLDER} -name "*.qual" \
         -type f ! -empty -print0 | \
        sort -k3,3n -k1,1d -k2,2n --merge --files0-from=- | \
        uniq --check-chars=40 > "${QUALITY_FILE}" &
}

build_distribution_file() {
    ## sequence <-> sample relations
    find ${DATA_FOLDER} -name "*.fas" \
     -type f ! -empty -execdir grep -H "^>" '{}' \; | \
    sed 's/.*\/// ; s/\.fas:>/\t/ ; s/;size=/\t/ ; s/;$//' | \
    awk 'BEGIN {FS = OFS = "\t"} {print $2, $1, $3}' > "${DISTRIBUTION_FILE}" &
}

list_all_cluster_seeds_of_size_greater_than_2() {
    find ${DATA_FOLDER} -name "*.stats" \
     -type f ! -empty -execdir grep -H "" '{}' \; | \
    sed 's/^\.\/// ; s/\.stats:/\t/' > "${POTENTIAL_SUB_SEEDS}" &
}

global_dereplication() {
    find ${DATA_FOLDER} -name "*.fas" \
         -type f ! -empty -execdir cat '{}' + | \
        "${VSEARCH}" \
            --derep_fulllength - \
            --sizein \
            --sizeout \
            --log "${LOG}" \
            --fasta_width 0 \
            --output "${FINAL_FASTA}"
}

clustering() {
    # swarm 3 or more recent
    "${SWARM}" \
        --differences "${RESOLUTION}" \
        --fastidious \
        --usearch-abundance \
        --threads "${THREADS}" \
        --internal-structure "${OUTPUT_STRUCT}" \
        --output-file "${OUTPUT_SWARMS}" \
        --statistics-file "${OUTPUT_STATS}" \
        --seeds "${OUTPUT_REPRESENTATIVES}" \
        "${FINAL_FASTA}" 2> "${OUTPUT_LOG}"
}

fake_taxonomic_assignment() {
    grep "^>" "${OUTPUT_REPRESENTATIVES}" | \
        sed -r 's/^>//
            s/;size=/\t/
            s/;?$/\t0.0\tNA\tNA/' > "${TAXONOMIC_ASSIGNMENTS}"
}

chimera_detection() {
    ## discard sequences with an abundance lower than FILTER
    ## and search for chimeras
    "${VSEARCH}" \
        --fastx_filter "${OUTPUT_REPRESENTATIVES}"  \
        --minsize "${FILTER}" \
        --fastaout - | \
        "${VSEARCH}" \
            --uchime_denovo - \
            --uchimeout "${UCHIME_RESULTS}" \
            2> "${UCHIME_LOG}"
}

cleaving() {
    python3 \
        "${SRC}/${OTU_CLEAVER}" \
        --global_stats "${OUTPUT_STATS}" \
        --per_sample_stats "${POTENTIAL_SUB_SEEDS}" \
        --struct "${OUTPUT_STRUCT}" \
        --swarms "${OUTPUT_SWARMS}" \
        --fasta "${FINAL_FASTA}"
}

fake_taxonomic_assignment2() {
    # refactoring: duplicated function
    grep "^>" "${OUTPUT_REPRESENTATIVES}2" | \
        sed -r 's/^>//
            s/;size=/\t/
            s/;?$/\t0.0\tNA\tNA/' > "${OUTPUT_REPRESENTATIVES%.*}.results2"
}

chimera_detection2() {
    ## chimera detection (only down to the smallest newly cleaved OTU)
    LOWEST_ABUNDANCE=$(sed -rn \
                           '/^>/ s/.*;size=([0-9]+);?/\1/p' \
                           "${OUTPUT_REPRESENTATIVES}2" | \
                           sort -n | \
                           head -n 1)

    # sort and filter by abundance (default to an abundance of 1),
    # search for chimeras
    cat "${OUTPUT_REPRESENTATIVES}" "${OUTPUT_REPRESENTATIVES}2" | \
        "${VSEARCH}" \
            --sortbysize - \
            --sizein \
            --minsize ${LOWEST_ABUNDANCE:-1} \
            --sizeout \
            --output - | \
        "${VSEARCH}" \
            --uchime_denovo - \
            --uchimeout "${UCHIME_RESULTS}2" \
            2> "${OUTPUT_REPRESENTATIVES%.*}.log2"

    unset LOWEST_ABUNDANCE
}

build_occurrence_table() {
    python3 \
        "${SRC}/${OTU_TABLE_BUILDER}" \
        --representatives <(cat "${OUTPUT_REPRESENTATIVES}"{,2}) \
        --stats <(cat "${OUTPUT_STATS}"{,2}) \
        --swarms <(cat "${OUTPUT_SWARMS}"{,2}) \
        --chimera <(cat "${UCHIME_RESULTS}"{,2}) \
        --quality "${QUALITY_FILE}" \
        --assignments <(cat "${TAXONOMIC_ASSIGNMENTS}"{2,}) \
        --distribution "${DISTRIBUTION_FILE}" > "${OTU_TABLE}"
}

extract_fasta_and_search_for_identical_sequences() {
    # excluding terminal gaps
    awk 'NR > 1 {printf ">"$1"\n"$10"\n"}' "${OTU_TABLE}" | \
        "${VSEARCH}" \
            --threads "${THREADS}" \
            --cluster_smallmem - \
            --id 1.0 \
            --qmask none \
            --usersort \
            --uc - | \
        grep "^H"
}

merge_superstrings() {
    # merge clusters that are sub-strings or super-strings of more
    # abundant clusters
    # refactoring: use long option names
    python3 \
        "${SRC}/${MERGE_SUBSTRINGS}" \
        -t "${OTU_TABLE}" \
        -m "${TMP_UC}" \
        -o /dev/stdout
}

check_that_number_of_reads_did_not_change() {
    local -ri BEFORE=$(awk 'NR > 1 {total += $2} END {print total}' "${OTU_TABLE}")
    local -ri AFTER=$(awk 'NR > 1 {total += $2} END {print total}' "${OUTPUT_TABLE}")
    if (( ${BEFORE} != ${AFTER} )) ; then
        echo "sub/superstring mergig went wrong" 1>&2
        exit 1
    fi
}

sort_occurrence_table() {
    (head -n 1 "${TMP_TABLE}"
     tail -n +2 "${TMP_TABLE}" | \
         sort -k1,1n) > "${OUTPUT_TABLE}"
}


extract_fasta_sequences_from_occurrence_table(){
    awk 'NR > 1 {printf ">"$4";size="$2";\n"$10"\n"}' "${OUTPUT_TABLE}" \
        > "${OUTPUT_TABLE/.table/.fas}"
}

trim_metadata_for_mumu() {
    cut -f 4,14- "${OUTPUT_TABLE}" > "${OUTPUT_TABLE%.*}_reduced.table"
}

find_similar_sequences() {
    # (recommended parameters for lulu), discard
    # abundance values
    "${VSEARCH}" \
        --usearch_global "${OUTPUT_TABLE/.table/.fas}" \
        --db "${OUTPUT_TABLE/.table/.fas}" \
        --self  \
        --threads "${THREADS}" \
        --id 0.84 \
        --iddef 1 \
        --userfields query+target+id \
        --maxaccepts 0 \
        --query_cov 0.9 \
        --maxhits 10 \
        --userout - | \
        sed -r 's/;size=[0-9]+;//g' > "${OUTPUT_TABLE%.*}.match_list"
}

run_mumu() {
    ${MUMU} \
        --otu_table "${OUTPUT_TABLE%.*}_reduced.table" \
        --match_list "${OUTPUT_TABLE%.*}.match_list" \
        --new_otu_table "${OUTPUT_TABLE%.*}_raw_mumu.table" \
        --log "${OUTPUT_TABLE%.*}.mumu.log"
}

rebuild_occurrence_table_after_mumu() {
    python3 \
        "${SRC}/${REBUILD_TABLE_AFTER_MUMU}" \
        --mumu_table "${OUTPUT_TABLE%.*}_raw_mumu.table" \
        --old_table "${OUTPUT_TABLE}" \
        > "${OUTPUT_TABLE%.*}.mumu.table"
}

extract_fasta_sequences_from_occurrence_table2(){
    # refactoring: duplicated function
    awk 'NR > 1 {printf ">"$4";size="$2";\n"$10"\n"}' "${OUTPUT_TABLE%.*}.mumu.table" \
    > "${OUTPUT_TABLE%.*}.mumu.fas"

}

taxonomic_assignment() {
    # method 'stampa' 
    (
        QUERY="$(readlink -f "${OUTPUT_TABLE%.*}.mumu.fas")"
        cd ${SRC}/
        bash "./${STAMPA}" "${QUERY}" "${MARKER}"
    )
}


## ------------------------------------------------------------------------ main

## ---------------------------------------------------------- global clustering
echo "run global clustering and chimera detection..."

build_expected_error_file
build_distribution_file
list_all_cluster_seeds_of_size_greater_than_2
global_dereplication
clustering
fake_taxonomic_assignment
chimera_detection
wait  # ...for background jobs to catch-up

## ------------------------------------------------------------------- cleaving
echo "run cleaving..."
cleaving
fake_taxonomic_assignment2
chimera_detection2


## ------------------------------------------------------------ first OTU table
echo "build first OTU table..."

build_occurrence_table


## ------------------------------------------------- merge sub- or superstrings
echo "merge sub- and super-strings..."
TMP_TABLE=$(mktemp)
TMP_UC=$(mktemp)

extract_fasta_and_search_for_identical_sequences > "${TMP_UC}"
merge_superstrings > "${TMP_TABLE}"

sort_occurrence_table
check_that_number_of_reads_did_not_change

# count deleted OTUs
wc -l "${OUTPUT_TABLE}" "${OTU_TABLE}"

# clean
rm "${TMP_TABLE}" "${TMP_UC}"
unset TMP_TABLE TMP_UC


## ------------------------------------------------------------- mumu (ex-lulu)
echo "run mumu..."

extract_fasta_sequences_from_occurrence_table
trim_metadata_for_mumu
find_similar_sequences
run_mumu

rm -f "${OUTPUT_TABLE%.*}_reduced.table"

rebuild_occurrence_table_after_mumu
# refactoring: fix cluster sorting after mumu
extract_fasta_sequences_from_occurrence_table2


## ------------------------------------------------------- taxonomic assignment

taxonomic_assignment


exit 0

# stop here: it is not trivial to wait for the end of the taxonomic
# assignment jobs

## ------------------------------------------------------ build final OTU table
echo "build final OTU table..."
NEW_TABLE=$(mktemp)
OTU_TABLE="${OUTPUT_TABLE%.*}.mumu.table"

python3 \
    "${SRC}/${OTU_TABLE_UPDATER}" \
    --old_otu_table "${OTU_TABLE}" \
    --new_taxonomy "${OUTPUT_TABLE%.*}.mumu.results" \
    --new_otu_table "${NEW_TABLE}"

# fix OTU sorting
(head -n 1 "${NEW_TABLE}"
 tail -n +2 "${NEW_TABLE}" | \
     sort -k2,2nr | \
     nl -n'ln' -w1 | \
     cut --complement -f 2
) > "${OTU_TABLE}2"

# clean up
chmod go+r,g-w "${OTU_TABLE}2"
rm "${NEW_TABLE}"

exit 0
