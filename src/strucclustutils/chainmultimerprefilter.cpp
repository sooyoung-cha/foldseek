#include "DBReader.h"
#include "DBWriter.h"
#include "Debug.h"
#include "LocalParameters.h"
#include "MemoryMapped.h"
#include "MultimerUtil.h"
#include "Util.h"

#include <algorithm>
#include <string>
#include <vector>

#ifdef OPENMP
#include <omp.h>
#endif

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
    for (size_t idx = 0; idx < queryComplexIds.size(); ++idx) {
        const unsigned int complexId = queryComplexIds[idx];
        if (par.monomerIncludeMode == SKIP_MONOMERS &&
            queryComplexToChains[complexId].size() < MULTIPLE_CHAINED_COMPLEX) {
            queryComplexAllowed[complexId] = 0;
        }
    }
    std::vector<char> targetComplexAllowed(targetComplexToChains.size(), 1);
    for (size_t idx = 0; idx < targetComplexIds.size(); ++idx) {
        const unsigned int complexId = targetComplexIds[idx];
        if (par.monomerIncludeMode == SKIP_MONOMERS &&
            targetComplexToChains[complexId].size() < MULTIPLE_CHAINED_COMPLEX) {
            targetComplexAllowed[complexId] = 0;
        }
    }

    std::vector<std::vector<unsigned int> > queryToTargets(queryChainToComplex.size());

    Debug::Progress progress(clusterDbr.getSize());
    for (size_t entryId = 0; entryId < clusterDbr.getSize(); ++entryId) {
        std::vector<unsigned int> clusterChainKeys;
        std::vector<unsigned int> queryChains;
        std::vector<unsigned int> targetChains;

        clusterChainKeys.push_back(clusterDbr.getDbKey(entryId));
        char *data = clusterDbr.getData(entryId, 0);
        while (*data != '\0') {
            clusterChainKeys.push_back(Util::fast_atoi<unsigned int>(data));
            data = Util::skipLine(data);
        }
        sortUnique(clusterChainKeys);

        for (size_t idx = 0; idx < clusterChainKeys.size(); ++idx) {
            const unsigned int chainKey = clusterChainKeys[idx];
            if (chainKey < queryChainToComplex.size() &&
                queryChainToComplex[chainKey] != NOT_AVAILABLE_CHAIN_KEY) {
                queryChains.push_back(chainKey);
            }
            if (chainKey < targetChainToComplex.size() &&
                targetChainToComplex[chainKey] != NOT_AVAILABLE_CHAIN_KEY &&
                targetComplexAllowed[targetChainToComplex[chainKey]] != 0) {
                targetChains.push_back(chainKey);
            }
        }

        if (!queryChains.empty() && !targetChains.empty()) {
            for (size_t qIdx = 0; qIdx < queryChains.size(); ++qIdx) {
                std::vector<unsigned int> &targets = queryToTargets[queryChains[qIdx]];
                targets.insert(targets.end(), targetChains.begin(), targetChains.end());
            }
        }
        progress.updateProgress();
    }
    clusterDbr.close();

    size_t allPossibleChainPairs = 0;
    for (size_t qIdx = 0; qIdx < queryComplexIds.size(); ++qIdx) {
        const size_t queryChainCount = queryComplexToChains[queryComplexIds[qIdx]].size();
        for (size_t tIdx = 0; tIdx < targetComplexIds.size(); ++tIdx) {
            if (targetComplexAllowed[targetComplexIds[tIdx]] == 0) {
                continue;
            }
            allPossibleChainPairs += queryChainCount * targetComplexToChains[targetComplexIds[tIdx]].size();
        }
    }

    DBWriter resultWriter(par.db4.c_str(), par.db4Index.c_str(), static_cast<unsigned int>(par.threads),
                          par.compressed, Parameters::DBTYPE_PREFILTER_RES);
    resultWriter.open();

    size_t keptChainPairs = 0;
    Debug::Progress outputProgress(queryComplexIds.size());
#pragma omp parallel for schedule(dynamic, 1) reduction(+:keptChainPairs)
    for (size_t qIdx = 0; qIdx < queryComplexIds.size(); ++qIdx) {
        unsigned int threadIdx = 0;
#ifdef OPENMP
        threadIdx = static_cast<unsigned int>(omp_get_thread_num());
#endif
        const unsigned int queryComplexId = queryComplexIds[qIdx];
        const std::vector<unsigned int> &queryChains = queryComplexToChains[queryComplexId];

        if (queryComplexAllowed[queryComplexId] == 0) {
            for (size_t chainIdx = 0; chainIdx < queryChains.size(); ++chainIdx) {
                resultWriter.writeData("", 0, queryChains[chainIdx], threadIdx);
            }
            outputProgress.updateProgress();
            continue;
        }

        for (size_t chainIdx = 0; chainIdx < queryChains.size(); ++chainIdx) {
            const unsigned int queryChainKey = queryChains[chainIdx];
            std::vector<unsigned int> targets;
            if (queryChainKey < queryToTargets.size()) {
                targets = queryToTargets[queryChainKey];
            }
            sortUnique(targets);

            std::string result;
            for (size_t targetIdx = 0; targetIdx < targets.size(); ++targetIdx) {
                result.append(SSTR(targets[targetIdx]));
                result.push_back('\n');
                keptChainPairs++;
            }
            resultWriter.writeData(result.c_str(), result.size(), queryChainKey, threadIdx);
        }
        outputProgress.updateProgress();
    }

    resultWriter.close(false);

    const double reduction = (allPossibleChainPairs == 0)
                                 ? 0.0
                                 : 100.0 * (1.0 - static_cast<double>(keptChainPairs) / static_cast<double>(allPossibleChainPairs));
    Debug(Debug::INFO) << "Candidate chain-pair reduction: " << reduction << "% (" << keptChainPairs
                       << "/" << allPossibleChainPairs << " kept)\n";
    return EXIT_SUCCESS;
}
