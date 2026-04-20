#include "DBReader.h"
#include "DBWriter.h"
#include "Debug.h"
#include "LocalParameters.h"
#include "MemoryMapped.h"
#include "MultimerUtil.h"
#include "Util.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <sys/stat.h>
#include <sys/types.h>
#include <string>
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

struct EdgePair {
    unsigned int queryChainKey;
    unsigned int targetChainKey;
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

static bool edgePairLess(const EdgePair &lhs, const EdgePair &rhs) {
    if (lhs.queryChainKey != rhs.queryChainKey) {
        return lhs.queryChainKey < rhs.queryChainKey;
    }
    return lhs.targetChainKey < rhs.targetChainKey;
}

static void flushChunk(std::vector<EdgePair> &buffer,
                       const std::string &tempDir,
                       unsigned int bucketId,
                       unsigned int threadIdx,
                       unsigned int chunkId,
                       std::vector<std::string> &chunkFiles) {
    if (buffer.empty()) {
        return;
    }

    std::sort(buffer.begin(), buffer.end(), edgePairLess);
    const std::string path = tempDir + "/bucket" + SSTR(bucketId) + ".thread" + SSTR(threadIdx) + ".chunk" + SSTR(chunkId);
    FILE *handle = std::fopen(path.c_str(), "wb");
    if (handle == NULL) {
        Debug(Debug::ERROR) << "Could not open chunk file " << path << " for writing\n";
        EXIT(EXIT_FAILURE);
    }
    const size_t written = std::fwrite(&buffer[0], sizeof(EdgePair), buffer.size(), handle);
    std::fclose(handle);
    if (written != buffer.size()) {
        Debug(Debug::ERROR) << "Could not fully write chunk file " << path << "\n";
        EXIT(EXIT_FAILURE);
    }
    chunkFiles.push_back(path);
    buffer.clear();
}

struct ChunkCursor {
    FILE *handle;
    EdgePair value;
    bool hasValue;
    std::string path;
};

static bool readNext(ChunkCursor &cursor) {
    cursor.hasValue = (std::fread(&cursor.value, sizeof(EdgePair), 1, cursor.handle) == 1);
    return cursor.hasValue;
}

static void ensureDirExists(const std::string &path) {
    if (mkdir(path.c_str(), 0775) != 0 && errno != EEXIST) {
        Debug(Debug::ERROR) << "Could not create temp directory " << path << "\n";
        EXIT(EXIT_FAILURE);
    }
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
        if (par.monomerIncludeMode == SKIP_MONOMERS &&
            queryComplexToChains[queryComplexIds[i]].size() < MULTIPLE_CHAINED_COMPLEX) {
            queryComplexAllowed[queryComplexIds[i]] = 0;
        }
    }
    std::vector<char> targetComplexAllowed(targetComplexToChains.size(), 1);
    for (size_t i = 0; i < targetComplexIds.size(); ++i) {
        if (par.monomerIncludeMode == SKIP_MONOMERS &&
            targetComplexToChains[targetComplexIds[i]].size() < MULTIPLE_CHAINED_COMPLEX) {
            targetComplexAllowed[targetComplexIds[i]] = 0;
        }
    }

    const unsigned int bucketCount = static_cast<unsigned int>(std::max(1, par.threads));
    const std::string tempDir = par.db4 + ".tmp";
    ensureDirExists(tempDir);
    const size_t maxBufferedPairs = 1000000;
    std::vector<std::vector<std::string> > chunkFilesByBucket(bucketCount);

    Debug::Progress progress(clusterDbr.getSize());
#pragma omp parallel
    {
        unsigned int threadIdx = 0;
#ifdef OPENMP
        threadIdx = static_cast<unsigned int>(omp_get_thread_num());
#endif
        std::vector<std::vector<EdgePair> > edgeBuffers(bucketCount);
        std::vector<unsigned int> chunkIds(bucketCount, 0);
        for (unsigned int bucketId = 0; bucketId < bucketCount; ++bucketId) {
            edgeBuffers[bucketId].reserve(std::max<size_t>(1024, maxBufferedPairs / bucketCount));
        }
        std::vector<std::vector<std::string> > localChunkFiles(bucketCount);

#pragma omp for schedule(dynamic, 1)
        for (size_t entryId = 0; entryId < clusterDbr.getSize(); ++entryId) {
            std::vector<unsigned int> clusterChainKeys;
            std::vector<ChainGroup> queryGroups;
            std::vector<ChainGroup> targetGroups;

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
                    findOrAddChainGroup(queryGroups, queryChainToComplex[chainKey]).chains.push_back(chainKey);
                }
                if (chainKey < targetChainToComplex.size() && targetChainToComplex[chainKey] != NOT_AVAILABLE_CHAIN_KEY) {
                    findOrAddChainGroup(targetGroups, targetChainToComplex[chainKey]).chains.push_back(chainKey);
                }
            }

            for (size_t qIdx = 0; qIdx < queryGroups.size(); ++qIdx) {
                const unsigned int queryComplexId = queryGroups[qIdx].complexId;
                if (queryComplexId >= queryComplexAllowed.size() || queryComplexAllowed[queryComplexId] == 0) {
                    continue;
                }
                for (size_t tIdx = 0; tIdx < targetGroups.size(); ++tIdx) {
                    const unsigned int targetComplexId = targetGroups[tIdx].complexId;
                    if (targetComplexId >= targetComplexAllowed.size() || targetComplexAllowed[targetComplexId] == 0) {
                        continue;
                    }

                    for (size_t qChainIdx = 0; qChainIdx < queryGroups[qIdx].chains.size(); ++qChainIdx) {
                        const unsigned int queryChainKey = queryGroups[qIdx].chains[qChainIdx];
                        const unsigned int bucketId = queryChainKey % bucketCount;
                        for (size_t tChainIdx = 0; tChainIdx < targetGroups[tIdx].chains.size(); ++tChainIdx) {
                            EdgePair pair;
                            pair.queryChainKey = queryChainKey;
                            pair.targetChainKey = targetGroups[tIdx].chains[tChainIdx];
                            edgeBuffers[bucketId].push_back(pair);
                        }
                        if (edgeBuffers[bucketId].size() >= maxBufferedPairs) {
                            flushChunk(edgeBuffers[bucketId], tempDir, bucketId, threadIdx, chunkIds[bucketId]++, localChunkFiles[bucketId]);
                        }
                    }
                }
            }
            progress.updateProgress();
        }

        for (unsigned int bucketId = 0; bucketId < bucketCount; ++bucketId) {
            flushChunk(edgeBuffers[bucketId], tempDir, bucketId, threadIdx, chunkIds[bucketId]++, localChunkFiles[bucketId]);
        }

#pragma omp critical
        {
            for (unsigned int bucketId = 0; bucketId < bucketCount; ++bucketId) {
                chunkFilesByBucket[bucketId].insert(chunkFilesByBucket[bucketId].end(),
                                                    localChunkFiles[bucketId].begin(),
                                                    localChunkFiles[bucketId].end());
            }
        }
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
    unsigned int maxQueryChainKey = 0;
    for (size_t qIdx = 0; qIdx < queryComplexIds.size(); ++qIdx) {
        const std::vector<unsigned int> &queryChains = queryComplexToChains[queryComplexIds[qIdx]];
        for (size_t chainIdx = 0; chainIdx < queryChains.size(); ++chainIdx) {
            maxQueryChainKey = std::max(maxQueryChainKey, queryChains[chainIdx]);
        }
    }
    std::vector<char> seenQueryChains(maxQueryChainKey + 1, 0);

    std::vector<size_t> keptChainPairsByBucket(bucketCount, 0);
#pragma omp parallel for schedule(dynamic, 1)
    for (int bucketInt = 0; bucketInt < static_cast<int>(bucketCount); ++bucketInt) {
        const unsigned int bucketId = static_cast<unsigned int>(bucketInt);
        std::vector<ChunkCursor> cursors;
        cursors.reserve(chunkFilesByBucket[bucketId].size());
        for (size_t chunkIdx = 0; chunkIdx < chunkFilesByBucket[bucketId].size(); ++chunkIdx) {
            ChunkCursor cursor;
            cursor.path = chunkFilesByBucket[bucketId][chunkIdx];
            cursor.handle = std::fopen(cursor.path.c_str(), "rb");
            if (cursor.handle == NULL) {
                Debug(Debug::ERROR) << "Could not open chunk file " << cursor.path << " for reading\n";
                EXIT(EXIT_FAILURE);
            }
            if (readNext(cursor)) {
                cursors.push_back(cursor);
            } else {
                std::fclose(cursor.handle);
                std::remove(cursor.path.c_str());
            }
        }

        unsigned int currentQueryChain = NOT_AVAILABLE_CHAIN_KEY;
        std::vector<unsigned int> currentTargets;
        EdgePair lastWrittenPair;
        bool haveLastWrittenPair = false;
        while (!cursors.empty()) {
            size_t bestIdx = 0;
            for (size_t idx = 1; idx < cursors.size(); ++idx) {
                if (edgePairLess(cursors[idx].value, cursors[bestIdx].value)) {
                    bestIdx = idx;
                }
            }
            const EdgePair pair = cursors[bestIdx].value;
            if (haveLastWrittenPair == false ||
                lastWrittenPair.queryChainKey != pair.queryChainKey ||
                lastWrittenPair.targetChainKey != pair.targetChainKey) {
                if (currentQueryChain != pair.queryChainKey) {
                    if (currentQueryChain != NOT_AVAILABLE_CHAIN_KEY) {
                        std::string result;
                        for (size_t targetIdx = 0; targetIdx < currentTargets.size(); ++targetIdx) {
                            result.append(SSTR(currentTargets[targetIdx]));
                            result.push_back('\n');
                            keptChainPairsByBucket[bucketId]++;
                        }
                        resultWriter.writeData(result.c_str(), result.size(), currentQueryChain, bucketId % static_cast<unsigned int>(par.threads));
                        if (currentQueryChain < seenQueryChains.size()) {
                            seenQueryChains[currentQueryChain] = 1;
                        }
                    }
                    currentQueryChain = pair.queryChainKey;
                    currentTargets.clear();
                }
                currentTargets.push_back(pair.targetChainKey);
                lastWrittenPair = pair;
                haveLastWrittenPair = true;
            }

            if (!readNext(cursors[bestIdx])) {
                std::fclose(cursors[bestIdx].handle);
                std::remove(cursors[bestIdx].path.c_str());
                cursors[bestIdx] = cursors.back();
                cursors.pop_back();
            }
        }

        if (currentQueryChain != NOT_AVAILABLE_CHAIN_KEY) {
            std::string result;
            for (size_t targetIdx = 0; targetIdx < currentTargets.size(); ++targetIdx) {
                result.append(SSTR(currentTargets[targetIdx]));
                result.push_back('\n');
                keptChainPairsByBucket[bucketId]++;
            }
            resultWriter.writeData(result.c_str(), result.size(), currentQueryChain, bucketId % static_cast<unsigned int>(par.threads));
            if (currentQueryChain < seenQueryChains.size()) {
                seenQueryChains[currentQueryChain] = 1;
            }
        }
    }

    for (unsigned int bucketId = 0; bucketId < bucketCount; ++bucketId) {
        keptChainPairs += keptChainPairsByBucket[bucketId];
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
    rmdir(tempDir.c_str());
    clusterDbr.close();

    const double reduction = (allPossibleChainPairs == 0)
                                 ? 0.0
                                 : 100.0 * (1.0 - static_cast<double>(keptChainPairs) / static_cast<double>(allPossibleChainPairs));
    Debug(Debug::INFO) << "Candidate chain-pair reduction: " << reduction << "% (" << keptChainPairs
                       << "/" << allPossibleChainPairs << " kept)\n";
    return EXIT_SUCCESS;
}
