#include "DBReader.h"
#include "DBWriter.h"
#include "Debug.h"
#include "LocalParameters.h"
#include "MemoryMapped.h"
#include "MultimerUtil.h"
#include "Util.h"

#include <algorithm>
#include <cmath>
#include <vector>

#ifdef OPENMP
#include <omp.h>
#endif

struct ChainGroup {
    unsigned int complexId;
    std::vector<unsigned int> chains;
};

struct QueryToTargetChains {
    unsigned int queryChainKey;
    std::vector<unsigned int> targetChainKeys;
};

struct ComplexCandidate {
    unsigned int targetComplexId;
    std::vector<unsigned int> matchedQueryChains;
    std::vector<unsigned int> matchedTargetChains;
    std::vector<QueryToTargetChains> edges;
};

static void sortUnique(std::vector<unsigned int> &values) {
    SORT_SERIAL(values.begin(), values.end());
    values.erase(std::unique(values.begin(), values.end()), values.end());
}

template <typename ReaderType>
static void loadComplexLookup(
    ReaderType &dbr,
    const std::string &file,
    std::vector<unsigned int> &chainToComplex,
    std::vector<std::vector<unsigned int> > &complexToChains,
    std::vector<unsigned int> &complexIds
) {
    if (file.empty()) {
        return;
    }

    MemoryMapped lookupDB(file, MemoryMapped::WholeFile, MemoryMapped::SequentialScan);
    char *data = (char *) lookupDB.getData();
    char *end = data + lookupDB.mappedSize();
    const char *entry[255];
    unsigned int maxComplexId = 0;
    unsigned int maxChainKey = 0;

    while (data < end && *data != '\0') {
        const size_t columns = Util::getWordsOfLine(data, entry, 255);
        if (columns >= 3) {
            maxChainKey = std::max(maxChainKey, Util::fast_atoi<unsigned int>(entry[0]));
            maxComplexId = std::max(maxComplexId, Util::fast_atoi<unsigned int>(entry[2]));
        }
        data = Util::skipLine(data);
    }

    chainToComplex.assign(maxChainKey + 1, NOT_AVAILABLE_CHAIN_KEY);
    complexToChains.assign(maxComplexId + 1, std::vector<unsigned int>());
    std::vector<char> seenComplex(maxComplexId + 1, 0);

    data = (char *) lookupDB.getData();
    end = data + lookupDB.mappedSize();
    while (data < end && *data != '\0') {
        const size_t columns = Util::getWordsOfLine(data, entry, 255);
        if (columns >= 3) {
            const unsigned int chainKey = Util::fast_atoi<unsigned int>(entry[0]);
            const unsigned int complexId = Util::fast_atoi<unsigned int>(entry[2]);
            if (getChainId(dbr, chainKey) != NOT_AVAILABLE_CHAIN_KEY) {
                chainToComplex[chainKey] = complexId;
                complexToChains[complexId].push_back(chainKey);
                if (seenComplex[complexId] == 0) {
                    complexIds.push_back(complexId);
                    seenComplex[complexId] = 1;
                }
            }
        }
        data = Util::skipLine(data);
    }
    lookupDB.close();
}

static ChainGroup &findOrAddChainGroup(std::vector<ChainGroup> &groups, unsigned int complexId) {
    for (size_t i = 0; i < groups.size(); ++i) {
        if (groups[i].complexId == complexId) {
            return groups[i];
        }
    }
    groups.push_back(ChainGroup());
    groups.back().complexId = complexId;
    return groups.back();
}

static ComplexCandidate &findOrAddCandidate(std::vector<ComplexCandidate> &candidates, unsigned int targetComplexId) {
    for (size_t i = 0; i < candidates.size(); ++i) {
        if (candidates[i].targetComplexId == targetComplexId) {
            return candidates[i];
        }
    }
    candidates.push_back(ComplexCandidate());
    candidates.back().targetComplexId = targetComplexId;
    return candidates.back();
}

static QueryToTargetChains &findOrAddEdge(std::vector<QueryToTargetChains> &edges, unsigned int queryChainKey) {
    for (size_t i = 0; i < edges.size(); ++i) {
        if (edges[i].queryChainKey == queryChainKey) {
            return edges[i];
        }
    }
    edges.push_back(QueryToTargetChains());
    edges.back().queryChainKey = queryChainKey;
    return edges.back();
}

int chainmultimerprefilter(int argc, const char **argv, const Command &command) {
    LocalParameters &par = LocalParameters::getLocalInstance();
    par.parseParameters(argc, argv, command, true, 0, MMseqsParameter::COMMAND_ALIGN);

    DBReader<unsigned int> clusterDbr(par.db3.c_str(), par.db3Index.c_str(), par.threads,
                                      DBReader<unsigned int>::USE_INDEX | DBReader<unsigned int>::USE_DATA);
    clusterDbr.open(DBReader<unsigned int>::NOSORT);

    const bool touch = par.preloadMode != Parameters::PRELOAD_MODE_MMAP;
    IndexReader queryDbr(par.db1, par.threads, IndexReader::SEQUENCES,
                         touch ? IndexReader::PRELOAD_INDEX : 0,
                         DBReader<unsigned int>::USE_INDEX);
    IndexReader targetDbr(par.db2, par.threads, IndexReader::SEQUENCES,
                          touch ? IndexReader::PRELOAD_INDEX : 0,
                          DBReader<unsigned int>::USE_INDEX);

    std::vector<unsigned int> queryComplexIds;
    std::vector<unsigned int> targetComplexIds;
    std::vector<unsigned int> queryChainToComplex;
    std::vector<unsigned int> targetChainToComplex;
    std::vector<std::vector<unsigned int> > queryComplexToChains;
    std::vector<std::vector<unsigned int> > targetComplexToChains;
    loadComplexLookup(queryDbr, par.db1 + ".lookup", queryChainToComplex, queryComplexToChains, queryComplexIds);
    loadComplexLookup(targetDbr, par.db2 + ".lookup", targetChainToComplex, targetComplexToChains, targetComplexIds);

    std::vector<std::vector<ComplexCandidate> > candidatesByQueryComplex(queryComplexToChains.size());
    Debug::Progress progress(clusterDbr.getSize());
    for (size_t entryId = 0; entryId < clusterDbr.getSize(); ++entryId) {
        std::vector<unsigned int> clusterChainKeys;
        std::vector<ChainGroup> queryGroups;
        std::vector<ChainGroup> targetGroups;

        clusterChainKeys.push_back(clusterDbr.getDbKey(entryId));
        char *data = clusterDbr.getData(entryId, 0);
        while (*data != '\0') {
            clusterChainKeys.push_back(Util::fast_atoi<unsigned int>(data));
            data = Util::skipLine(data);
        }
        sortUnique(clusterChainKeys);

        for (size_t idx = 0; idx < clusterChainKeys.size(); ++idx) {
            const unsigned int chainKey = clusterChainKeys[idx];
            if (chainKey < queryChainToComplex.size() && queryChainToComplex[chainKey] != NOT_AVAILABLE_CHAIN_KEY) {
                findOrAddChainGroup(queryGroups, queryChainToComplex[chainKey]).chains.push_back(chainKey);
            }
            if (chainKey < targetChainToComplex.size() && targetChainToComplex[chainKey] != NOT_AVAILABLE_CHAIN_KEY) {
                findOrAddChainGroup(targetGroups, targetChainToComplex[chainKey]).chains.push_back(chainKey);
            }
        }

        for (size_t qIdx = 0; qIdx < queryGroups.size(); ++qIdx) {
            for (size_t tIdx = 0; tIdx < targetGroups.size(); ++tIdx) {
                ComplexCandidate &candidate = findOrAddCandidate(candidatesByQueryComplex[queryGroups[qIdx].complexId],
                                                                 targetGroups[tIdx].complexId);
                candidate.matchedQueryChains.insert(candidate.matchedQueryChains.end(),
                                                    queryGroups[qIdx].chains.begin(), queryGroups[qIdx].chains.end());
                candidate.matchedTargetChains.insert(candidate.matchedTargetChains.end(),
                                                     targetGroups[tIdx].chains.begin(), targetGroups[tIdx].chains.end());
                for (size_t chainIdx = 0; chainIdx < queryGroups[qIdx].chains.size(); ++chainIdx) {
                    QueryToTargetChains &edge = findOrAddEdge(candidate.edges, queryGroups[qIdx].chains[chainIdx]);
                    edge.targetChainKeys.insert(edge.targetChainKeys.end(),
                                                targetGroups[tIdx].chains.begin(), targetGroups[tIdx].chains.end());
                }
            }
        }
        progress.updateProgress();
    }

    DBWriter resultWriter(par.db4.c_str(), par.db4Index.c_str(), static_cast<unsigned int>(par.threads),
                          par.compressed, Parameters::DBTYPE_PREFILTER_RES);
    resultWriter.open();

    size_t allPossibleChainPairs = 0;
    for (size_t qIdx = 0; qIdx < queryComplexIds.size(); ++qIdx) {
        const size_t queryChainCount = queryComplexToChains[queryComplexIds[qIdx]].size();
        for (size_t tIdx = 0; tIdx < targetComplexIds.size(); ++tIdx) {
            allPossibleChainPairs += queryChainCount * targetComplexToChains[targetComplexIds[tIdx]].size();
        }
    }

    size_t keptChainPairs = 0;
    Debug::Progress outputProgress(queryComplexIds.size());
#pragma omp parallel reduction(+:keptChainPairs)
    {
        unsigned int thread_idx = 0;
#ifdef OPENMP
        thread_idx = static_cast<unsigned int>(omp_get_thread_num());
#endif
#pragma omp for schedule(dynamic, 1)
        for (size_t qIdx = 0; qIdx < queryComplexIds.size(); ++qIdx) {
            const unsigned int queryComplexId = queryComplexIds[qIdx];
            const std::vector<unsigned int> &queryChains = queryComplexToChains[queryComplexId];
            if (par.monomerIncludeMode == SKIP_MONOMERS && queryChains.size() < MULTIPLE_CHAINED_COMPLEX) {
                for (size_t chainIdx = 0; chainIdx < queryChains.size(); ++chainIdx) {
                    resultWriter.writeData("", 0, queryChains[chainIdx], thread_idx);
                }
                outputProgress.updateProgress();
                continue;
            }

            std::vector<QueryToTargetChains> outputEdges;
            std::vector<ComplexCandidate> &complexCandidates = candidatesByQueryComplex[queryComplexId];
            for (size_t candIdx = 0; candIdx < complexCandidates.size(); ++candIdx) {
                ComplexCandidate &candidate = complexCandidates[candIdx];
                if (candidate.targetComplexId >= targetComplexToChains.size()) {
                    continue;
                }
                const std::vector<unsigned int> &targetChains = targetComplexToChains[candidate.targetComplexId];
                if (par.monomerIncludeMode == SKIP_MONOMERS && targetChains.size() < MULTIPLE_CHAINED_COMPLEX) {
                    continue;
                }

                for (size_t edgeIdx = 0; edgeIdx < candidate.edges.size(); ++edgeIdx) {
                    QueryToTargetChains &outEdge = findOrAddEdge(outputEdges, candidate.edges[edgeIdx].queryChainKey);
                    outEdge.targetChainKeys.insert(outEdge.targetChainKeys.end(),
                                                   candidate.edges[edgeIdx].targetChainKeys.begin(),
                                                   candidate.edges[edgeIdx].targetChainKeys.end());
                }
            }

            for (size_t chainIdx = 0; chainIdx < queryChains.size(); ++chainIdx) {
                const unsigned int queryChainKey = queryChains[chainIdx];
                std::string result;
                QueryToTargetChains &edge = findOrAddEdge(outputEdges, queryChainKey);
                sortUnique(edge.targetChainKeys);
                for (size_t targetIdx = 0; targetIdx < edge.targetChainKeys.size(); ++targetIdx) {
                    result.append(SSTR(edge.targetChainKeys[targetIdx]));
                    result.push_back('\n');
                    keptChainPairs++;
                }
                resultWriter.writeData(result.c_str(), result.size(), queryChainKey, thread_idx);
            }
            outputProgress.updateProgress();
        }
    }

    resultWriter.close(false);
    clusterDbr.close();

    const double reduction = (allPossibleChainPairs == 0)
                                 ? 0.0
                                 : 100.0 * (1.0 - static_cast<double>(keptChainPairs) / static_cast<double>(allPossibleChainPairs));
    Debug(Debug::INFO) << "Candidate chain-pair reduction: " << reduction << "% (" << keptChainPairs
                       << "/" << allPossibleChainPairs << " kept)\n";
    return EXIT_SUCCESS;
}
