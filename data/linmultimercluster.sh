#!/bin/sh -e
fail() {
    echo "Error: $1"
    exit 1
}

notExists() {
	[ ! -f "$1" ]
}

exists() {
	[ -f "$1" ]
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

# Chain groups, so the quadratic stage below only has to run over one chain per group.
#
# Near identical chains have near identical hit lists. linclust finds them for a linear
# cost, the prefilter runs over the representatives alone, and multimerprefilter and
# multimerneighbors read a member's hits off its representative and a hit on a
# representative as a hit on its whole group. The expanded list is never written: at
# 10M multimers it would be tens of billions of records.
PREF_IN="${INPUT}_ss"
PREF_TARGET="${INPUT}_ss${INDEXEXT}"
if [ -n "${GROUP_PREFILTER}" ] && [ "$PREFMODE" != "EXHAUSTIVE" ]; then
    if notExists "${TMP_PATH}/cgroups.dbtype"; then
        mkdir -p "${TMP_PATH}/cgroups_tmp"
        # shellcheck disable=SC2086
        "$MMSEQS" linclust "${INPUT}_ss" "${TMP_PATH}/cgroups" "${TMP_PATH}/cgroups_tmp" \
            ${LINCLUST_PAR} || fail "linclust of the chains died"
    fi
    if notExists "${TMP_PATH}/grp_ss.dbtype"; then
        awk '{ print $1 }' "${TMP_PATH}/cgroups.index" | sort -n > "${TMP_PATH}/cgroup_rep.keys"
        # shellcheck disable=SC2086
        "$MMSEQS" createsubdb "${TMP_PATH}/cgroup_rep.keys" "${INPUT}_ss" "${TMP_PATH}/grp_ss" \
            ${CREATESUBDB_PAR} || fail "createsubdb of the group representatives died"
    fi
    PREF_IN="${TMP_PATH}/grp_ss"
    PREF_TARGET="${TMP_PATH}/grp_ss"
fi

# Chain level prefilter. Its hits are a superset of what an alignment would keep,
# which is all the candidate filter below needs, so the chain alignment the regular
# multimer workflow runs here is skipped entirely.
if notExists "${TMP_PATH}/pref.dbtype"; then
    if [ "$PREFMODE" = "EXHAUSTIVE" ]; then
        fake_pref "${INPUT}_ss" "${INPUT}_ss" "${TMP_PATH}/pref"
    elif [ "$PREFMODE" = "UNGAPPED" ]; then
        # shellcheck disable=SC2086
        $RUNNER "$MMSEQS" ungappedprefilter "${PREF_IN}" "${PREF_TARGET}" "${TMP_PATH}/pref" ${UNGAPPEDPREFILTER_PAR} \
            || fail "Ungapped prefilter matching step died"
    else
        # shellcheck disable=SC2086
        $RUNNER "$MMSEQS" prefilter "${PREF_IN}" "${PREF_TARGET}" "${TMP_PATH}/pref" ${PREFILTER0_PAR} \
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

round 0
if notExists "${TMP_PATH}/done_0.dbtype"; then
    # shellcheck disable=SC2086
    "$MMSEQS" lndb "${TMP_PATH}/multimer_pref_0" "${TMP_PATH}/done_0" ${VERBOSITY_PAR}
fi

R=0

# Second prefilter, over the representatives of round 0 only.
#
# Round 0 is allowed to be cheap and therefore lossy. What it loses is not lost quietly:
# a multimer whose true neighbour was never proposed joined no cluster and came out of
# clust as its own representative. So the representatives are exactly the multimers that
# still have something to find, and running the full prefilter on them is the second half
# of the cascade. It is quadratic like any prefilter, but on reps/all of the database, so
# it costs (reps/all)^2 of what the same call would have cost in round 0.
#
# The hits are keyed by the original chain keys, so the round below runs against the full
# input unchanged: only representative chains have an entry, so only representative
# multimers emit candidates, and the ones round 0 already clustered are skipped for free.
if [ "$PREFMODE" != "EXHAUSTIVE" ]; then
    if notExists "${TMP_PATH}/rep_ss.dbtype"; then
        awk '{ print $1 }' "${TMP_PATH}/clu_0.index" | sort -n > "${TMP_PATH}/rep_mult.keys"
        awk 'NR == FNR { r[$1]; next } ($3 in r) { print $1 }' \
            "${TMP_PATH}/rep_mult.keys" "${INPUT}.lookup" > "${TMP_PATH}/rep_chain.keys"
        # shellcheck disable=SC2086
        "$MMSEQS" createsubdb "${TMP_PATH}/rep_chain.keys" "${INPUT}_ss" "${TMP_PATH}/rep_ss" \
            ${CREATESUBDB_PAR} || fail "createsubdb of the representative chains died"
    fi
    if notExists "${TMP_PATH}/pref2.dbtype"; then
        if [ "$PREFMODE" = "UNGAPPED" ]; then
            # shellcheck disable=SC2086
            $RUNNER "$MMSEQS" ungappedprefilter "${TMP_PATH}/rep_ss" "${TMP_PATH}/rep_ss" \
                "${TMP_PATH}/pref2" ${UNGAPPEDPREFILTER_PAR} \
                || fail "Ungapped prefilter over the representatives died"
        else
            # shellcheck disable=SC2086
            $RUNNER "$MMSEQS" prefilter "${TMP_PATH}/rep_ss" "${TMP_PATH}/rep_ss" "${TMP_PATH}/pref2" \
                ${PREFILTER_PAR} || fail "Prefilter over the representatives died"
        fi
    fi
    if notExists "${TMP_PATH}/multimer_pref_1.dbtype"; then
        # shellcheck disable=SC2086
        "$MMSEQS" multimerprefilter "${INPUT}" "${INPUT}" "${TMP_PATH}/pref2" \
            "${TMP_PATH}/multimer_pref_1" ${MULTIMERPREFILTER_PAR2} \
            || fail "multimerprefilter died"
    fi
    round 1
    R=1
fi

# The expansion rounds read this instead of the chain prefilter. It has to be folded from
# every chain prefilter the workflow ran, not just the first one: the rounds walk the
# neighbourhood of each multimer, and a cheap round 0 leaves that neighbourhood too thin
# for them to move anything. Both prefilters are keyed by chain, and multimerneighbors
# keeps the best score per target multimer, so merging them first costs one pass and
# cannot double count.
if notExists "${TMP_PATH}/nbr.dbtype"; then
    if exists "${TMP_PATH}/cgroups.dbtype"; then
        # The two prefilters live in different key spaces now: round 0's is filed under
        # the group representatives, the pass over the representatives under the real
        # chains. So each is folded to multimer level on its own terms and the two are
        # concatenated; multimerexpand keeps the best score per target multimer, so a
        # neighbour found by both costs nothing.
        if notExists "${TMP_PATH}/nbr_0.dbtype"; then
            # shellcheck disable=SC2086
            "$MMSEQS" multimerneighbors "${INPUT}" "${INPUT}" "${TMP_PATH}/pref" "${TMP_PATH}/nbr_0" \
                ${MULTIMERNEIGHBORS_PAR} || fail "multimerneighbors died"
        fi
        if exists "${TMP_PATH}/pref2.dbtype"; then
            if notExists "${TMP_PATH}/nbr_1.dbtype"; then
                # shellcheck disable=SC2086
                "$MMSEQS" multimerneighbors "${INPUT}" "${INPUT}" "${TMP_PATH}/pref2" "${TMP_PATH}/nbr_1" \
                    ${MULTIMERNEIGHBORS_PAR2} || fail "multimerneighbors died"
            fi
            # shellcheck disable=SC2086
            "$MMSEQS" mergedbs "${INPUT}_ss" "${TMP_PATH}/nbr" "${TMP_PATH}/nbr_0" \
                "${TMP_PATH}/nbr_1" ${VERBOSITY_PAR} || fail "mergedbs died"
        else
            # shellcheck disable=SC2086
            "$MMSEQS" mvdb "${TMP_PATH}/nbr_0" "${TMP_PATH}/nbr" ${VERBOSITY_PAR}
        fi
    else
        NBR_IN="${TMP_PATH}/pref"
        if exists "${TMP_PATH}/pref2.dbtype"; then
            if notExists "${TMP_PATH}/pref_all.dbtype"; then
                # shellcheck disable=SC2086
                "$MMSEQS" mergedbs "${INPUT}_ss" "${TMP_PATH}/pref_all" "${TMP_PATH}/pref" \
                    "${TMP_PATH}/pref2" ${VERBOSITY_PAR} || fail "mergedbs died"
            fi
            NBR_IN="${TMP_PATH}/pref_all"
        fi
        # shellcheck disable=SC2086
        "$MMSEQS" multimerneighbors "${INPUT}" "${INPUT}" "${NBR_IN}" "${TMP_PATH}/nbr" \
            ${MULTIMERNEIGHBORS_PAR} || fail "multimerneighbors died"
    fi
fi

# Expansion rounds. A top-k graph caps how large a cluster can grow, because set-cover can
# only absorb what the centre is joined to: on 57,448 dimers the largest cluster fell from
# 788 to 184 at k=50. Mode 0 puts every multimer in front of the centres of the clusters
# its own neighbours landed in, which is the far members of a big cluster meeting their
# centre, and mode 1 puts the centres in front of each other, which is what moves the
# representatives that could still have merged. Both emit a constant number of pairs per
# multimer and skip whatever earlier rounds already aligned.
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
    if exists "${TMP_PATH}/cgroups.dbtype"; then
        # shellcheck disable=SC2086
        "$MMSEQS" rmdb "${TMP_PATH}/cgroups" ${VERBOSITY_PAR}
        # shellcheck disable=SC2086
        "$MMSEQS" rmdb "${TMP_PATH}/grp_ss" ${VERBOSITY_PAR}
        # shellcheck disable=SC2086
        "$MMSEQS" rmdb "${TMP_PATH}/nbr_0" ${VERBOSITY_PAR}
        # shellcheck disable=SC2086
        "$MMSEQS" rmdb "${TMP_PATH}/nbr_1" ${VERBOSITY_PAR}
        rm -f -- "${TMP_PATH}/cgroup_rep.keys"
        rm -rf -- "${TMP_PATH}/cgroups_tmp"
    fi
    if exists "${TMP_PATH}/pref2.dbtype"; then
        # shellcheck disable=SC2086
        "$MMSEQS" rmdb "${TMP_PATH}/pref2" ${VERBOSITY_PAR}
        # shellcheck disable=SC2086
        "$MMSEQS" rmdb "${TMP_PATH}/rep_ss" ${VERBOSITY_PAR}
        # shellcheck disable=SC2086
        "$MMSEQS" rmdb "${TMP_PATH}/pref_all" ${VERBOSITY_PAR}
        rm -f -- "${TMP_PATH}/rep_mult.keys" "${TMP_PATH}/rep_chain.keys"
    fi
    rm -f -- "${TMP_PATH}/linmultimercluster.sh"
fi
