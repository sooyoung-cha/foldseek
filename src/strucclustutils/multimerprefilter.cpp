#include "DBReader.h"
#include "DBWriter.h"
#include "Debug.h"
#include "Util.h"
#include "Matcher.h"
#include "LocalParameters.h"
#include "MemoryMapped.h"
#include "MultimerUtil.h"
#include <map>
#include <vector>
#include <algorithm>
#ifdef OPENMP
#include <omp.h>
#endif

// Drop in replacement for expandmultimer.
//
// expandmultimer takes every target multimer that was hit by any single chain and
// hands the full |queryChains| x |targetChains| cross product to the alignment,
// which then has to align everything and scoremultimer has to score it. Most of
// that work is dead on arrival: scoremultimer only keeps an assignment when every
// chain of the smaller multimer finds a partner, so a multimer pair whose chains
// cannot be matched one to one is thrown away at the very end anyway.
//
// So decide it here instead, on the chain alignments that the search already
// produced. Two chains may be partners when the search aligned them; two multimers
// survive when that relation admits a matching saturating the smaller multimer.
// Only the chain pairs taking part in the relation are emitted, not the cross
// product.
//
// This never removes a multimer pair that scoremultimer would have accepted: the
// chain alignments used here are the same ones the old path expanded from, and
// they are produced far more permissively (-e 10) than the chain filter
// scoremultimer applies later.

static const int MULTIMER_PREFILTER_MATCHED_PAIRS = 1;

struct ChainHit {
    size_t qChainIdx;
    size_t dbChainIdx;
    int score;
};

// A multimer pair that survived, waiting to be emitted. Kept so the top-k rule can rank
// them before anything is written.
struct Survivor {
    unsigned int dbComplexId;
    int score;
    std::vector<int> matchedToQ;
};

static bool byScoreDesc(const Survivor &a, const Survivor &b) {
    if (a.score != b.score) {
        return a.score > b.score;
    }
    return a.dbComplexId < b.dbComplexId;   // deterministic
}

// Kuhn's algorithm, the graphs here have a handful of nodes per side
static bool tryAugment(
    size_t qIdx,
    const std::vector<std::vector<uint8_t> > &compatible,
    std::vector<char> &visited,
    std::vector<int> &matchedToQ
) {
    for (size_t dbIdx = 0; dbIdx < matchedToQ.size(); dbIdx++) {
        if (compatible[qIdx][dbIdx] == 0 || visited[dbIdx]) {
            continue;
        }
        visited[dbIdx] = 1;
        if (matchedToQ[dbIdx] == -1 || tryAugment(matchedToQ[dbIdx], compatible, visited, matchedToQ)) {
            matchedToQ[dbIdx] = static_cast<int>(qIdx);
            return true;
        }
    }
    return false;
}

static size_t getMaximumMatching(
    const std::vector<std::vector<uint8_t> > &compatible,
    size_t dbChainNum,
    std::vector<int> *matchedToQOut = NULL
) {
    std::vector<int> matchedToQ(dbChainNum, -1);
    std::vector<char> visited(dbChainNum);
    size_t matching = 0;
    for (size_t qIdx = 0; qIdx < compatible.size(); qIdx++) {
        std::fill(visited.begin(), visited.end(), 0);
        if (tryAugment(qIdx, compatible, visited, matchedToQ)) {
            matching++;
        }
    }
    if (matchedToQOut != NULL) {
        matchedToQOut->swap(matchedToQ);
    }
    return matching;
}

// The matching above has to hold on what the *chain* prefilter reported, and there it
// squares the per chain sensitivity: with p the chance that one chain finds its true
// partner chain, keeping a multimer pair on a single edge succeeds with 1-(1-p)^2 but
// demanding two independent edges succeeds with only p^2. p is small because max-seqs
// caps every hit list, so the filter quietly shreds the large clusters on a real run
// while looking perfect on any small benchmark.
//
// Measured on 57,448 CATH sampled dimers against what multimercluster accepts: of the
// 66.1M multimer pairs whose chains hit at all, only 13.3M admit a matching, and 95.9%
// of the failures have exactly one of the four chain pairs, its hit sitting at median
// rank 616 of a list that is 99.9% saturated at the cap. 40.0% of the accepted multimer
// pairs are lost this way.
//
// A two hop rescue over the top-30 hit lists recovered 0.0% of them: the thirty nearest
// neighbours of two chains in a large family are disjoint (99.8% share nothing). What
// does separate them is the score of the one hit that did survive, which the prefilter
// already wrote out and which costs nothing to read:
//
//     score of the surviving chain hit     p25   median   p75
//       pairs multimercluster accepts      248     296    384
//       pairs multimercluster rejects       88     124    191
//
// So a multimer pair with no matching is kept when its best chain hit scores at least
// --multimer-prefilter-min-score. Candidates against recall of multimercluster's edges:
//
// A chain pair is also thrown out when the two chains cannot be of the same fold at
// all: min(qTM, tTM) is normalised by the longer chain, so a pair can only reach the
// chain TM threshold t when min(len) / max(len) >= t. That bound is free, it costs a
// division, and it removed 16% of the candidates here without touching recall.
//
// Neither rule bounds the work: the number of multimers within the threshold of a given
// one grows with the database (neighbours per multimer measured at 25.6, 59.7, 121.6 and
// 247.0 for 5.7k, 14.4k, 28.7k and 57.4k dimers), so anything that enumerates them is
// quadratic. --multimer-prefilter-topk keeps only the N best scoring candidates of each
// query instead, which fixes the work per multimer and makes the stage linear. What it
// gives up is multimer pairs that would have merged two clusters, and the workflow wins
// those back with the cheap expansion rounds rather than by looking at everything.
//
//     off (matching only)   13.3M   61.6%
//     >= 250                22.3M   90.0%
//     >= 200                27.3M   97.1%
//     >= 150                35.0M   99.7%
//     expandmultimer        47.8M  100.0%

int multimerprefilter(int argc, const char **argv, const Command &command) {
    LocalParameters &par = LocalParameters::getLocalInstance();
    par.parseParameters(argc, argv, command, true, 0, MMseqsParameter::COMMAND_ALIGN);

    DBReader<unsigned int> alnDbr(par.db3.c_str(), par.db3Index.c_str(), par.threads, DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA);
    alnDbr.open(DBReader<unsigned int>::LINEAR_ACCCESS);

    int dbType = Parameters::DBTYPE_CLUSTER_RES;
    uint16_t extended = DBReader<unsigned int>::getExtendedDbtype(alnDbr.getDbtype());
    bool needSrc = false;
    if (extended & Parameters::DBTYPE_EXTENDED_INDEX_NEED_SRC) {
        needSrc = true;
        dbType = DBReader<unsigned int>::setExtendedDbtype(dbType, Parameters::DBTYPE_EXTENDED_INDEX_NEED_SRC);
    }
    DBWriter resultWriter(par.db4.c_str(), par.db4Index.c_str(), static_cast<unsigned int>(par.threads), par.compressed, dbType);
    resultWriter.open();

    const bool touch = par.preloadMode != Parameters::PRELOAD_MODE_MMAP;
    IndexReader tDbr(
        par.db2, par.threads,
        needSrc ? IndexReader::SRC_SEQUENCES : IndexReader::SEQUENCES,
        touch ? IndexReader::PRELOAD_INDEX : 0, DBReader<unsigned int>::USE_INDEX
    );
    IndexReader qDbr(
        par.db1, par.threads,
        needSrc ? IndexReader::SRC_SEQUENCES : IndexReader::SEQUENCES,
        touch ? IndexReader::PRELOAD_INDEX : 0, DBReader<unsigned int>::USE_INDEX
    );

    std::vector<unsigned int> qComplexIndices;
    std::vector<unsigned int> dbComplexIndices;
    chainKeyToComplexId_t qChainKeyToComplexIdMap, dbChainKeyToComplexIdMap;
    complexIdToChainKeys_t qComplexIdToChainKeysMap, dbComplexIdToChainKeysMap;
    std::string qLookupFile = par.db1 + ".lookup";
    std::string dbLookupFile = needSrc ? par.db2 + "_seq.lookup" : par.db2 + ".lookup";
    getKeyToIdMapIdToKeysMapIdVec(qDbr, qLookupFile, qChainKeyToComplexIdMap, qComplexIdToChainKeysMap, qComplexIndices);
    getKeyToIdMapIdToKeysMapIdVec(tDbr, dbLookupFile, dbChainKeyToComplexIdMap, dbComplexIdToChainKeysMap, dbComplexIndices);
    dbComplexIndices.clear();
    qChainKeyToComplexIdMap.clear();

    const int minScore = par.multimerPrefilterMinScore;
    const size_t topK = (size_t) par.multimerPrefilterTopk;
    // min(qTM, tTM) >= chainTmThr is impossible below this length ratio
    const float minLenRatio = par.filtChainTmThr;
    if (minScore > 0) {
        Debug(Debug::INFO) << "Keeping unmatched multimer pairs at chain hit score "
                           << minScore << " or better\n";
    }
    if (minLenRatio > 0.0f) {
        Debug(Debug::INFO) << "Dropping chain pairs below length ratio " << minLenRatio << "\n";
    }
    if (topK > 0) {
        Debug(Debug::INFO) << "Keeping the " << topK << " best scoring candidates per multimer\n";
    }

    Debug::Progress progress(qComplexIndices.size());
#pragma omp parallel
    {
        unsigned int thread_idx = 0;
#ifdef OPENMP
        thread_idx = static_cast<unsigned int>(omp_get_thread_num());
#endif
        resultToWrite_t result;
        // target complex -> hit chain pairs, as indices into the two chain lists
        // (query chain, target chain, hit score)
        std::map<unsigned int, std::vector<ChainHit> > dbComplexToHits;
        std::map<unsigned int, size_t> dbChainKeyToIdx;
        std::vector<std::vector<uint8_t> > compatible;
        std::vector<std::vector<uint8_t> > lenOk;
        std::vector<std::vector<unsigned int> > candidateChainKeys;
        std::vector<Survivor> survivors;
#pragma omp for schedule(dynamic, 1)
        for (size_t qCompIdx = 0; qCompIdx < qComplexIndices.size(); qCompIdx++) {
            const unsigned int qComplexId = qComplexIndices[qCompIdx];
            const std::vector<unsigned int> &qChainKeys = qComplexIdToChainKeysMap.at(qComplexId);
            const size_t qChainNum = qChainKeys.size();

            dbComplexToHits.clear();
            for (size_t qChainIdx = 0; qChainIdx < qChainNum; qChainIdx++) {
                const unsigned int qId = alnDbr.getId(qChainKeys[qChainIdx]);
                if (qId == NOT_AVAILABLE_CHAIN_KEY) {
                    continue;
                }
                char *data = alnDbr.getData(qId, thread_idx);
                while (*data != '\0') {
                    char dbKeyBuffer[255 + 1];
                    Util::parseKey(data, dbKeyBuffer);
                    const unsigned int dbChainKey = static_cast<unsigned int>(strtoul(dbKeyBuffer, NULL, 10));
                    chainKeyToComplexId_t::const_iterator cIt = dbChainKeyToComplexIdMap.find(dbChainKey);
                    if (cIt == dbChainKeyToComplexIdMap.end()) {
                        data = Util::skipLine(data);
                        continue;
                    }
                    const unsigned int dbComplexId = cIt->second;
                    const std::vector<unsigned int> &dbChainKeys = dbComplexIdToChainKeysMap.at(dbComplexId);
                    // position of this target chain inside its own multimer
                    size_t dbChainIdx = dbChainKeys.size();
                    for (size_t j = 0; j < dbChainKeys.size(); j++) {
                        if (dbChainKeys[j] == dbChainKey) {
                            dbChainIdx = j;
                            break;
                        }
                    }
                    if (dbChainIdx < dbChainKeys.size()) {
                        int score = 0;
                        if (minScore > 0 || topK > 0) {
                            // second column of a prefilter or alignment record
                            const char *p = data;
                            while (*p != '\t' && *p != '\n' && *p != '\0') {
                                p++;
                            }
                            if (*p == '\t') {
                                score = (int) strtol(p + 1, NULL, 10);
                            }
                        }
                        ChainHit hit = { qChainIdx, dbChainIdx, score };
                        dbComplexToHits[dbComplexId].push_back(hit);
                    }
                    data = Util::skipLine(data);
                }
            }

            candidateChainKeys.assign(qChainNum, std::vector<unsigned int>());
            survivors.clear();
            for (std::map<unsigned int, std::vector<ChainHit> >::const_iterator it = dbComplexToHits.begin();
                 it != dbComplexToHits.end(); ++it) {
                const std::vector<unsigned int> &dbChainKeys = dbComplexIdToChainKeysMap.at(it->first);
                const size_t dbChainNum = dbChainKeys.size();
                // scoremultimer needs every chain of both multimers aligned in cov-mode 0
                if (par.covMode == Parameters::COV_MODE_BIDIRECTIONAL && dbChainNum != qChainNum) {
                    continue;
                }

                // a pairing the chain TM threshold could still accept
                lenOk.assign(qChainNum, std::vector<uint8_t>(dbChainNum, 1));
                if (minLenRatio > 0.0f) {
                    for (size_t i = 0; i < qChainNum; i++) {
                        const unsigned int qId = qDbr.sequenceReader->getId(qChainKeys[i]);
                        if (qId == NOT_AVAILABLE_CHAIN_KEY) {
                            continue;
                        }
                        const float qLen = (float) qDbr.sequenceReader->getSeqLen(qId);
                        for (size_t j = 0; j < dbChainNum; j++) {
                            const unsigned int tId = tDbr.sequenceReader->getId(dbChainKeys[j]);
                            if (tId == NOT_AVAILABLE_CHAIN_KEY) {
                                lenOk[i][j] = 0;
                                continue;
                            }
                            const float tLen = (float) tDbr.sequenceReader->getSeqLen(tId);
                            const float r = qLen < tLen ? qLen / tLen : tLen / qLen;
                            lenOk[i][j] = (r >= minLenRatio) ? 1 : 0;
                        }
                    }
                    if (getMaximumMatching(lenOk, dbChainNum) < std::min(qChainNum, dbChainNum)) {
                        continue;
                    }
                }

                compatible.assign(qChainNum, std::vector<uint8_t>(dbChainNum, 0));
                int bestScore = 0;
                for (size_t h = 0; h < it->second.size(); h++) {
                    const ChainHit &hit = it->second[h];
                    // a hit the length bound already rules out is no evidence at all
                    if (lenOk[hit.qChainIdx][hit.dbChainIdx] == 0) {
                        continue;
                    }
                    compatible[hit.qChainIdx][hit.dbChainIdx] = 1;
                    if (hit.score > bestScore) {
                        bestScore = hit.score;
                    }
                }

                // every chain of the smaller multimer has to find its own partner
                const size_t needed = std::min(qChainNum, dbChainNum);
                std::vector<int> matchedToQ;
                bool byMatching = getMaximumMatching(compatible, dbChainNum, &matchedToQ) >= needed;
                if (byMatching == false) {
                    // no one to one matching: fall back on the one hit that did survive
                    if (topK == 0 && minScore <= 0) {
                        continue;
                    }
                    if (topK == 0 && bestScore < minScore) {
                        continue;
                    }
                    // extend the partial matching the hits gave us, without ever
                    // displacing a pair that actually hit: the chain pair that carried
                    // this multimer pair here has to be the one the alignment sees
                    std::vector<char> usedQ(qChainNum, 0);
                    for (size_t j = 0; j < dbChainNum; j++) {
                        if (matchedToQ[j] >= 0) {
                            usedQ[matchedToQ[j]] = 1;
                        }
                    }
                    for (size_t i = 0; i < qChainNum; i++) {
                        if (usedQ[i]) {
                            continue;
                        }
                        for (size_t j = 0; j < dbChainNum; j++) {
                            if (matchedToQ[j] < 0 && lenOk[i][j]) {
                                matchedToQ[j] = (int) i;
                                usedQ[i] = 1;
                                break;
                            }
                        }
                    }
                }
                Survivor sv;
                sv.dbComplexId = it->first;
                sv.score = bestScore;
                sv.matchedToQ.swap(matchedToQ);
                survivors.push_back(sv);
            }

            // a constant number of neighbours per multimer, the best scoring ones
            if (topK > 0 && survivors.size() > topK) {
                std::partial_sort(survivors.begin(), survivors.begin() + topK, survivors.end(), byScoreDesc);
                survivors.resize(topK);
            }

            for (size_t si = 0; si < survivors.size(); si++) {
                const std::vector<unsigned int> &dbChainKeys = dbComplexIdToChainKeysMap.at(survivors[si].dbComplexId);
                const size_t dbChainNum = dbChainKeys.size();
                const std::vector<int> &matchedToQ = survivors[si].matchedToQ;
                if (par.multimerPrefilterMode == MULTIMER_PREFILTER_MATCHED_PAIRS) {
                    // only the pairs that carried the multimer pair here. scoremultimer
                    // then has exactly one assignment to score instead of the square.
                    for (size_t j = 0; j < dbChainNum; j++) {
                        const int i = matchedToQ[j];
                        if (i < 0) {
                            continue;
                        }
                        if (tDbr.sequenceReader->getId(dbChainKeys[j]) == NOT_AVAILABLE_CHAIN_KEY) {
                            continue;
                        }
                        candidateChainKeys[i].emplace_back(dbChainKeys[j]);
                    }
                } else {
                    for (size_t i = 0; i < qChainNum; i++) {
                        for (size_t j = 0; j < dbChainNum; j++) {
                            if (tDbr.sequenceReader->getId(dbChainKeys[j]) == NOT_AVAILABLE_CHAIN_KEY) {
                                continue;
                            }
                            candidateChainKeys[i].emplace_back(dbChainKeys[j]);
                        }
                    }
                }
            }

            for (size_t qChainIdx = 0; qChainIdx < qChainNum; qChainIdx++) {
                std::vector<unsigned int> &keys = candidateChainKeys[qChainIdx];
                SORT_SERIAL(keys.begin(), keys.end());
                keys.erase(std::unique(keys.begin(), keys.end()), keys.end());
                result.clear();
                for (size_t k = 0; k < keys.size(); k++) {
                    result.append(SSTR(keys[k]));
                    result.push_back('\n');
                }
                // an entry is written even when empty, downstream modules expect every
                // query chain to be present
                resultWriter.writeData(result.c_str(), result.length(), qChainKeys[qChainIdx], thread_idx);
            }
            progress.updateProgress();
        }
    }

    alnDbr.close();
    resultWriter.close(false);
    return EXIT_SUCCESS;
}
