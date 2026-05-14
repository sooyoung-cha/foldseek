#!/bin/sh -e
fail() {
    echo "Error: $1"
    exit 1
}

notExists() {
    [ ! -f "$1" ]
}

if [ -n "${USE_ORIGINAL_MULTIMERCLUSTER}" ]; then
    if notExists "${RESULT}.dbtype"; then
        # shellcheck disable=SC2086
        "$MMSEQS" multimercluster "${INPUT}" "${RESULT}" "${TMP_PATH}" ${MULTIMERCLUSTER_PAR} \
            || fail "multimercluster died"
    fi
else
    CHAIN_CLUSTER_TM_THRESHOLD=$(awk -v thr="${CHAIN_TM_THRESHOLD}" 'BEGIN { printf "%.6f", thr - 0.15 }')

    if notExists "${TMP_PATH}/chain_clu.dbtype"; then
        # shellcheck disable=SC2086
        "$MMSEQS" cluster "${INPUT}" "${TMP_PATH}/chain_clu" "${TMP_PATH}/chaincluster_tmp" --remove-tmp-files 0 --tmscore-threshold "${CHAIN_CLUSTER_TM_THRESHOLD}" -c "${CHAIN_CLUSTER_C}" --cov-mode "${CHAIN_CLUSTER_COV_MODE}" ${CHAINCLUSTER_PAR} \
            || fail "chain cluster died"
    fi

    if notExists "${TMP_PATH}/candidate_pref.dbtype"; then
        # shellcheck disable=SC2086
        "$MMSEQS" chainmultimerprefilter "${INPUT}" "${INPUT}" "${TMP_PATH}/chain_clu" "${TMP_PATH}/candidate_pref" ${CHAINMULTIMERPREFILTER_PAR} \
            || fail "chainmultimerprefilter died"
    fi

    if notExists "${TMP_PATH}/candidate_aln.dbtype"; then
        # shellcheck disable=SC2086
        "$MMSEQS" structurealign "${INPUT}" "${INPUT}" "${TMP_PATH}/candidate_pref" "${TMP_PATH}/candidate_aln" ${STRUCTUREALIGN_PAR} \
            || fail "structurealign died"
    fi

    if notExists "${TMP_PATH}/multimer_result.dbtype"; then
        # shellcheck disable=SC2086
        "$MMSEQS" scoremultimerfast "${INPUT}" "${INPUT}" "${TMP_PATH}/candidate_aln" "${TMP_PATH}/multimer_result" ${SCOREMULTIMER_PAR} \
            || fail "scoremultimerfast died"
    fi

    if notExists "${RESULT}.dbtype"; then
        # shellcheck disable=SC2086
        "$MMSEQS" clust "${INPUT}" "${TMP_PATH}/multimer_result" "${RESULT}" ${CLUSTER_PAR} \
            || fail "Clustering died"
    fi
fi

if [ -n "${REMOVE_TMP}" ]; then
    if [ -z "${USE_ORIGINAL_MULTIMERCLUSTER}" ]; then
        # shellcheck disable=SC2086
        "$MMSEQS" rmdb "${TMP_PATH}/multimer_result" ${VERBOSITY_PAR} || fail "rmdb died"
        # shellcheck disable=SC2086
        "$MMSEQS" rmdb "${TMP_PATH}/candidate_pref" ${VERBOSITY_PAR} || fail "rmdb died"
        # shellcheck disable=SC2086
        "$MMSEQS" rmdb "${TMP_PATH}/candidate_aln" ${VERBOSITY_PAR} || fail "rmdb died"
        # shellcheck disable=SC2086
        "$MMSEQS" rmdb "${TMP_PATH}/chain_clu" ${VERBOSITY_PAR} || fail "rmdb died"
        rm -rf -- "${TMP_PATH}/chaincluster_tmp"
    fi
    rm -f -- "${TMP_PATH}/multimercluster_fast.sh"
fi
