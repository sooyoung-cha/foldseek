#!/bin/sh -e
fail() {
    echo "Error: $1"
    exit 1
}

notExists() {
	[ ! -f "$1" ]
}

# Chain level k-mer prefilter. Its hits are a superset of what an alignment would
# keep, which is all the candidate filter below needs, so the chain alignment the
# regular multimer workflow runs here is skipped entirely.
if notExists "${TMP_PATH}/pref.dbtype"; then
    # shellcheck disable=SC2086
    $RUNNER "$MMSEQS" prefilter "${INPUT}_ss" "${INPUT}_ss" "${TMP_PATH}/pref" ${PREFILTER_PAR} \
        || fail "Prefilter died"
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
