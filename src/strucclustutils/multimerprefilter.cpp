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

static size_t getMaximumMatching(const std::vector<std::vector<uint8_t> > &compatible, size_t dbChainNum) {
    std::vector<int> matchedToQ(dbChainNum, -1);
    std::vector<char> visited(dbChainNum);
    size_t matching = 0;
    for (size_t qIdx = 0; qIdx < compatible.size(); qIdx++) {
        std::fill(visited.begin(), visited.end(), 0);
        if (tryAugment(qIdx, compatible, visited, matchedToQ)) {
            matching++;
        }
    }
    return matching;
}

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

    Debug::Progress progress(qComplexIndices.size());
#pragma omp parallel
    {
        unsigned int thread_idx = 0;
#ifdef OPENMP
        thread_idx = static_cast<unsigned int>(omp_get_thread_num());
#endif
        resultToWrite_t result;
        // target complex -> hit chain pairs, as indices into the two chain lists
        std::map<unsigned int, std::vector<std::pair<size_t, size_t> > > dbComplexToHits;
        std::map<unsigned int, size_t> dbChainKeyToIdx;
        std::vector<std::vector<uint8_t> > compatible;
        std::vector<std::vector<unsigned int> > candidateChainKeys;
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
                        dbComplexToHits[dbComplexId].emplace_back(qChainIdx, dbChainIdx);
                    }
                    data = Util::skipLine(data);
                }
            }

            candidateChainKeys.assign(qChainNum, std::vector<unsigned int>());
            for (std::map<unsigned int, std::vector<std::pair<size_t, size_t> > >::const_iterator it = dbComplexToHits.begin();
                 it != dbComplexToHits.end(); ++it) {
                const std::vector<unsigned int> &dbChainKeys = dbComplexIdToChainKeysMap.at(it->first);
                const size_t dbChainNum = dbChainKeys.size();
                // scoremultimer needs every chain of both multimers aligned in cov-mode 0
                if (par.covMode == Parameters::COV_MODE_BIDIRECTIONAL && dbChainNum != qChainNum) {
                    continue;
                }

                compatible.assign(qChainNum, std::vector<uint8_t>(dbChainNum, 0));
                for (size_t h = 0; h < it->second.size(); h++) {
                    compatible[it->second[h].first][it->second[h].second] = 1;
                }

                // every chain of the smaller multimer has to find its own partner
                const size_t needed = std::min(qChainNum, dbChainNum);
                if (getMaximumMatching(compatible, dbChainNum) < needed) {
                    continue;
                }
                for (size_t i = 0; i < qChainNum; i++) {
                    for (size_t j = 0; j < dbChainNum; j++) {
                        // mode 0 keeps the whole cross product of a surviving multimer
                        // pair, so the second alignment can still rescue chain pairs
                        // that the search reported below its own E-value
                        if (par.multimerPrefilterMode != 0 && compatible[i][j] == 0) {
                            continue;
                        }
                        if (tDbr.sequenceReader->getId(dbChainKeys[j]) == NOT_AVAILABLE_CHAIN_KEY) {
                            continue;
                        }
                        candidateChainKeys[i].emplace_back(dbChainKeys[j]);
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
