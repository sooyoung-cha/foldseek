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

# Candidate multimer pairs, a constant number per multimer.
#
# Keeping every multimer pair within the threshold is quadratic: the number of multimers
# that close to a given one grows with the database (neighbours per multimer measured at
# 25.6, 59.7, 121.6 and 247.0 for 5.7k, 14.4k, 28.7k and 57.4k dimers). So the first round
# keeps only the best scoring --multimer-prefilter-topk candidates of each multimer, which
# fixes the work per multimer, and the expansion rounds below win back what that drops.
if notExists "${TMP_PATH}/multimer_pref_0.dbtype"; then
    # shellcheck disable=SC2086
    "$MMSEQS" multimerprefilter "${INPUT}" "${INPUT}" "${TMP_PATH}/pref" "${TMP_PATH}/multimer_pref_0" \
        ${MULTIMERPREFILTER_PAR} \
        || fail "multimerprefilter died"
fi

round () {
    R="$1"
    if notExists "${TMP_PATH}/multimer_aln_${R}.dbtype"; then
        # shellcheck disable=SC2086
        "$MMSEQS" structurealign "${INPUT}" "${INPUT}" "${TMP_PATH}/multimer_pref_${R}" \
            "${TMP_PATH}/multimer_aln_${R}" ${MULTIMERALIGN_PAR} || fail "structurealign died"
    fi
    if notExists "${TMP_PATH}/multimer_res_${R}.dbtype"; then
        # shellcheck disable=SC2086
        "$MMSEQS" scoremultimer "${INPUT}" "${INPUT}" "${TMP_PATH}/multimer_aln_${R}" \
            "${TMP_PATH}/multimer_res_${R}" ${SCOREMULTIMER_PAR} || fail "scoremultimer died"
    fi
    if [ "${R}" = "0" ]; then
        if notExists "${TMP_PATH}/multimer_result.dbtype"; then
            # shellcheck disable=SC2086
            "$MMSEQS" mvdb "${TMP_PATH}/multimer_res_0" "${TMP_PATH}/multimer_result" ${VERBOSITY_PAR}
        fi
    else
        if notExists "${TMP_PATH}/multimer_merged_${R}.dbtype"; then
            # shellcheck disable=SC2086
            "$MMSEQS" mergedbs "${INPUT}" "${TMP_PATH}/multimer_merged_${R}" \
                "${TMP_PATH}/multimer_result" "${TMP_PATH}/multimer_res_${R}" ${VERBOSITY_PAR} \
                || fail "mergedbs died"
            # shellcheck disable=SC2086
            "$MMSEQS" rmdb "${TMP_PATH}/multimer_result" ${VERBOSITY_PAR}
            # shellcheck disable=SC2086
            "$MMSEQS" mvdb "${TMP_PATH}/multimer_merged_${R}" "${TMP_PATH}/multimer_result" ${VERBOSITY_PAR}
        fi
    fi
    if notExists "${TMP_PATH}/clu_${R}.dbtype"; then
        # shellcheck disable=SC2086
        "$MMSEQS" clust "${INPUT}" "${TMP_PATH}/multimer_result" "${TMP_PATH}/clu_${R}" ${CLUSTER_PAR} \
            || fail "Clustering died"
    fi
    # everything handed to the alignment so far, so a pair is never aligned twice
    if [ "${R}" != "0" ] && notExists "${TMP_PATH}/done_${R}.dbtype"; then
        # shellcheck disable=SC2086
        "$MMSEQS" mergedbs "${INPUT}" "${TMP_PATH}/done_${R}" "${TMP_PATH}/done_$((R-1))" \
            "${TMP_PATH}/multimer_pref_${R}" ${VERBOSITY_PAR} || fail "mergedbs died"
    fi
}

# the expansion rounds read this instead of the chain prefilter
if notExists "${TMP_PATH}/nbr.dbtype"; then
    # shellcheck disable=SC2086
    "$MMSEQS" multimerneighbors "${INPUT}" "${INPUT}" "${TMP_PATH}/pref" "${TMP_PATH}/nbr" \
        ${MULTIMERNEIGHBORS_PAR} || fail "multimerneighbors died"
fi

round 0
if notExists "${TMP_PATH}/done_0.dbtype"; then
    # shellcheck disable=SC2086
    "$MMSEQS" lndb "${TMP_PATH}/multimer_pref_0" "${TMP_PATH}/done_0" ${VERBOSITY_PAR}
fi

# Expansion rounds. A top-k graph caps how large a cluster can grow, because set-cover can
# only absorb what the centre is joined to: on 57,448 dimers the largest cluster fell from
# 788 to 184 at k=50. Mode 0 puts every multimer in front of the centres of the clusters
# its own neighbours landed in, which is the far members of a big cluster meeting their
# centre, and mode 1 puts the centres in front of each other, which is what moves the
# representatives that could still have merged. Both emit a constant number of pairs per
# multimer and skip whatever earlier rounds already aligned.
R=0
for MODE in ${EXPAND_SCHEDULE}; do
    PREV="${R}"
    R="$((R+1))"
    if [ "${MODE}" = "0" ]; then
        EXPAND_PAR="${MULTIMEREXPAND_PAR0}"
    else
        EXPAND_PAR="${MULTIMEREXPAND_PAR1}"
    fi
    if notExists "${TMP_PATH}/multimer_pref_${R}.dbtype"; then
        # shellcheck disable=SC2086
        "$MMSEQS" multimerexpand "${INPUT}" "${INPUT}" "${TMP_PATH}/nbr" "${TMP_PATH}/clu_${PREV}" \
            "${TMP_PATH}/done_${PREV}" "${TMP_PATH}/multimer_pref_${R}" \
            ${EXPAND_PAR} || fail "multimerexpand died"
    fi
    round "${R}"
done

if notExists "${RESULT}.dbtype"; then
    # shellcheck disable=SC2086
    "$MMSEQS" mvdb "${TMP_PATH}/clu_${R}" "${RESULT}" ${VERBOSITY_PAR}
fi

if [ -n "${REMOVE_TMP}" ]; then
    I=0
    while [ "${I}" -le "${R}" ]; do
        # shellcheck disable=SC2086
        "$MMSEQS" rmdb "${TMP_PATH}/multimer_pref_${I}" ${VERBOSITY_PAR}
        # shellcheck disable=SC2086
        "$MMSEQS" rmdb "${TMP_PATH}/multimer_aln_${I}" ${VERBOSITY_PAR}
        # shellcheck disable=SC2086
        "$MMSEQS" rmdb "${TMP_PATH}/done_${I}" ${VERBOSITY_PAR}
        # shellcheck disable=SC2086
        "$MMSEQS" rmdb "${TMP_PATH}/clu_${I}" ${VERBOSITY_PAR}
        I="$((I+1))"
    done
    # shellcheck disable=SC2086
    "$MMSEQS" rmdb "${TMP_PATH}/multimer_result" ${VERBOSITY_PAR}
    # shellcheck disable=SC2086
    "$MMSEQS" rmdb "${TMP_PATH}/nbr" ${VERBOSITY_PAR}
    # shellcheck disable=SC2086
    "$MMSEQS" rmdb "${TMP_PATH}/pref" ${VERBOSITY_PAR}
    rm -f -- "${TMP_PATH}/linmultimercluster.sh"
fi
