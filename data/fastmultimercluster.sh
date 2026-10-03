#!/bin/sh -e
fail() {
    echo "Error: $1"
    exit 1
}

notExists() {
	[ ! -f "$1" ]
}

abspath() {
    if [ -d "$1" ]; then
        (cd "$1"; pwd)
    elif [ -f "$1" ]; then
        if [ -z "${1##*/*}" ]; then
            echo "$(cd "${1%/*}"; pwd)/${1##*/}"
        else
            echo "$(pwd)/$1"
        fi
    elif [ -d "$(dirname "$1")" ]; then
        echo "$(cd "$(dirname "$1")"; pwd)/$(basename "$1")"
    fi
}

fake_pref() {
    QDB="$1"
    TDB="$2"
    RES="$3"
    # create link to data file which contains a list of all targets that should be aligned
    ln -s "$(abspath "${TDB}.index")" "${RES}"
    # create new index repeatedly pointing to same entry
    INDEX_SIZE="$(wc -c < "${TDB}.index")"
    awk -v size="$INDEX_SIZE" '{ print $1"\t0\t"size; }' "${QDB}.index" > "${RES}.index"
    # create dbtype (7)
    awk 'BEGIN { printf("%c%c%c%c",7,0,0,0); exit; }' > "${RES}.dbtype"
}

# Chain level prefilter. Its hits are a superset of what an alignment would keep,
# which is all the candidate filter below needs, so the chain alignment the regular
# multimer workflow runs here is skipped entirely.
if notExists "${TMP_PATH}/pref.dbtype"; then
    if [ "$PREFMODE" = "EXHAUSTIVE" ]; then
        fake_pref "${INPUT}_ss" "${INPUT}_ss" "${TMP_PATH}/pref"
    elif [ "$PREFMODE" = "UNGAPPED" ]; then
        # shellcheck disable=SC2086
        $RUNNER "$MMSEQS" ungappedprefilter "${INPUT}_ss" "${INPUT}_ss${INDEXEXT}" "${TMP_PATH}/pref" ${UNGAPPEDPREFILTER_PAR} \
            || fail "Ungapped prefilter matching step died"
    else
        # shellcheck disable=SC2086
        $RUNNER "$MMSEQS" prefilter "${INPUT}_ss" "${INPUT}_ss${INDEXEXT}" "${TMP_PATH}/pref" ${PREFILTER_PAR} \
            || fail "Kmer matching step died"
    fi
fi

# Candidate multimer pairs. Replaces expandmultimer: multimer pairs whose chains
# cannot be matched one to one are dropped here instead of after being aligned and
# scored, and only the chain pairs that actually hit are emitted rather than the
# full cross product.
if notExists "${TMP_PATH}/multimer_pref.dbtype"; then
    # shellcheck disable=SC2086
    "$MMSEQS" multimerprefilter "${INPUT}" "${INPUT}" "${TMP_PATH}/pref" "${TMP_PATH}/multimer_pref" ${MULTIMERPREFILTER_PAR} \
        || fail "multimerprefilter died"
fi

if notExists "${TMP_PATH}/multimer_aln.dbtype"; then
    # shellcheck disable=SC2086
    "$MMSEQS" structurealign "${INPUT}" "${INPUT}" "${TMP_PATH}/multimer_pref" "${TMP_PATH}/multimer_aln" ${MULTIMERALIGN_PAR} \
        || fail "structurealign died"
fi

if notExists "${TMP_PATH}/multimer_result.dbtype"; then
    # shellcheck disable=SC2086
    "$MMSEQS" scoremultimer "${INPUT}" "${INPUT}" "${TMP_PATH}/multimer_aln" "${TMP_PATH}/multimer_result" ${SCOREMULTIMER_PAR} \
        || fail "scoremultimer died"
fi

if notExists "${RESULT}.dbtype"; then
    # shellcheck disable=SC2086
    "$MMSEQS" clust "${INPUT}" "${TMP_PATH}/multimer_result" "${RESULT}" ${CLUSTER_PAR} \
        || fail "Clustering died"
fi

if [ -n "${REMOVE_TMP}" ]; then
    # shellcheck disable=SC2086
    "$MMSEQS" rmdb "${TMP_PATH}/multimer_result" ${VERBOSITY_PAR}
    # shellcheck disable=SC2086
    "$MMSEQS" rmdb "${TMP_PATH}/multimer_aln" ${VERBOSITY_PAR}
    # shellcheck disable=SC2086
    "$MMSEQS" rmdb "${TMP_PATH}/multimer_pref" ${VERBOSITY_PAR}
    # shellcheck disable=SC2086
    "$MMSEQS" rmdb "${TMP_PATH}/pref" ${VERBOSITY_PAR}
    rm -f -- "${TMP_PATH}/fastmultimercluster.sh"
fi
