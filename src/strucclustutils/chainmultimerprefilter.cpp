#include "DBReader.h"
#include "DBWriter.h"
#include "Debug.h"
#include "LocalParameters.h"
#include "MemoryMapped.h"
#include "MultimerUtil.h"
#include "Util.h"

#include <algorithm>
#include <cstdlib>
#include <string>
#include <vector>

#ifdef OPENMP
#include <omp.h>
#endif

struct ClusterChains {
    std::vector<unsigned int> queryChains;
    std::vector<unsigned int> targetChains;
    std::string targetResult;
};

static void sortUnique(std::vector<unsigned int> &values) {
    std::sort(values.begin(), values.end());
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

static void appendChainResult(std::string &result, unsigned int chainKey) {
    result.append(SSTR(chainKey));
    result.push_back('\n');
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

    std::vector<char> queryComplexAllowed(queryComplexToChains.size(), 1);
    for (size_t i = 0; i < queryComplexIds.size(); ++i) {
        const unsigned int complexId = queryComplexIds[i];
        if (par.monomerIncludeMode == SKIP_MONOMERS &&
            queryComplexToChains[complexId].size() < MULTIPLE_CHAINED_COMPLEX) {
            queryComplexAllowed[complexId] = 0;
        }
    }
    std::vector<char> targetComplexAllowed(targetComplexToChains.size(), 1);
    for (size_t i = 0; i < targetComplexIds.size(); ++i) {
        const unsigned int complexId = targetComplexIds[i];
        if (par.monomerIncludeMode == SKIP_MONOMERS &&
            targetComplexToChains[complexId].size() < MULTIPLE_CHAINED_COMPLEX) {
            targetComplexAllowed[complexId] = 0;
        }
    }

    DBWriter resultWriter(par.db4.c_str(), par.db4Index.c_str(), static_cast<unsigned int>(par.threads),
                          par.compressed, Parameters::DBTYPE_PREFILTER_RES);
    resultWriter.open();

    unsigned int maxQueryChainKey = 0;
    for (size_t qIdx = 0; qIdx < queryComplexIds.size(); ++qIdx) {
        const std::vector<unsigned int> &queryChains = queryComplexToChains[queryComplexIds[qIdx]];
        for (size_t chainIdx = 0; chainIdx < queryChains.size(); ++chainIdx) {
            maxQueryChainKey = std::max(maxQueryChainKey, queryChains[chainIdx]);
        }
    }
    std::vector<char> seenQueryChains(maxQueryChainKey + 1, 0);

    Debug(Debug::INFO) << "Pass 1/1: scanning chain clusters and writing prefilter DB\n";
    Debug::Progress scanProgress(clusterDbr.getSize());
    size_t keptChainPairs = 0;
#pragma omp parallel
    {
        unsigned int threadIdx = 0;
#ifdef OPENMP
        threadIdx = static_cast<unsigned int>(omp_get_thread_num());
#endif
        size_t localKeptChainPairs = 0;

#pragma omp for schedule(dynamic, 1)
        for (size_t entryId = 0; entryId < clusterDbr.getSize(); ++entryId) {
            std::vector<unsigned int> clusterChainKeys;
            ClusterChains cluster;

            clusterChainKeys.push_back(clusterDbr.getDbKey(entryId));
            char *data = clusterDbr.getData(entryId, threadIdx);
            while (*data != '\0') {
                clusterChainKeys.push_back(Util::fast_atoi<unsigned int>(data));
                data = Util::skipLine(data);
            }
            sortUnique(clusterChainKeys);

            for (size_t idx = 0; idx < clusterChainKeys.size(); ++idx) {
                const unsigned int chainKey = clusterChainKeys[idx];
                if (chainKey < queryChainToComplex.size() && queryChainToComplex[chainKey] != NOT_AVAILABLE_CHAIN_KEY) {
                    const unsigned int queryComplexId = queryChainToComplex[chainKey];
                    if (queryComplexAllowed[queryComplexId] != 0) {
                        cluster.queryChains.push_back(chainKey);
                    }
                }
                if (chainKey < targetChainToComplex.size() && targetChainToComplex[chainKey] != NOT_AVAILABLE_CHAIN_KEY) {
                    const unsigned int targetComplexId = targetChainToComplex[chainKey];
                    if (targetComplexAllowed[targetComplexId] != 0) {
                        cluster.targetChains.push_back(chainKey);
                    }
                }
            }

            if (!cluster.queryChains.empty() && !cluster.targetChains.empty()) {
                cluster.targetResult.reserve(cluster.targetChains.size() * 12);
                for (size_t targetIdx = 0; targetIdx < cluster.targetChains.size(); ++targetIdx) {
                    appendChainResult(cluster.targetResult, cluster.targetChains[targetIdx]);
                }
                localKeptChainPairs += cluster.queryChains.size() * cluster.targetChains.size();
                for (size_t queryIdx = 0; queryIdx < cluster.queryChains.size(); ++queryIdx) {
                    const unsigned int queryChainKey = cluster.queryChains[queryIdx];
                    resultWriter.writeData(cluster.targetResult.c_str(), cluster.targetResult.size(), queryChainKey, threadIdx);
                    if (queryChainKey < seenQueryChains.size()) {
                        seenQueryChains[queryChainKey] = 1;
                    }
                }
            }
            scanProgress.updateProgress();
        }
#pragma omp atomic
        keptChainPairs += localKeptChainPairs;
    }
    clusterDbr.close();

    size_t allPossibleChainPairs = 0;
    for (size_t qIdx = 0; qIdx < queryComplexIds.size(); ++qIdx) {
        const unsigned int queryComplexId = queryComplexIds[qIdx];
        if (queryComplexAllowed[queryComplexId] == 0) {
            continue;
        }
        const size_t queryChainCount = queryComplexToChains[queryComplexId].size();
        for (size_t tIdx = 0; tIdx < targetComplexIds.size(); ++tIdx) {
            const unsigned int targetComplexId = targetComplexIds[tIdx];
            if (targetComplexAllowed[targetComplexId] == 0) {
                continue;
            }
            allPossibleChainPairs += queryChainCount * targetComplexToChains[targetComplexId].size();
        }
    }

    for (size_t qIdx = 0; qIdx < queryComplexIds.size(); ++qIdx) {
        const unsigned int queryComplexId = queryComplexIds[qIdx];
        const std::vector<unsigned int> &queryChains = queryComplexToChains[queryComplexId];
        for (size_t chainIdx = 0; chainIdx < queryChains.size(); ++chainIdx) {
            const unsigned int queryChainKey = queryChains[chainIdx];
            if (queryChainKey >= seenQueryChains.size() || seenQueryChains[queryChainKey] == 0) {
                resultWriter.writeData("", 0, queryChainKey, 0);
            }
        }
    }

    resultWriter.close(false);

    const double reduction = (allPossibleChainPairs == 0)
                                 ? 0.0
                                 : 100.0 * (1.0 - static_cast<double>(keptChainPairs) / static_cast<double>(allPossibleChainPairs));
    Debug(Debug::INFO) << "Candidate chain-pair reduction: " << reduction << "% (" << keptChainPairs
                       << "/" << allPossibleChainPairs << " kept)\n";
    return EXIT_SUCCESS;
}
