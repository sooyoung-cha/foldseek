#include "DBReader.h"
#include "DBWriter.h"
#include "Debug.h"
#include "Util.h"
#include "LocalParameters.h"
#include "MultimerUtil.h"
#include <map>
#include <vector>
#include <algorithm>
#ifdef OPENMP
#include <omp.h>
#endif

// Fold the chain level prefilter into a multimer level one, once.
//
// Every expansion round walks the neighbourhood of every multimer, and walking it on the
// chain prefilter means re-parsing one line per chain hit and looking the chain up in the
// key to multimer map each time: 439 M lines per round at 219 k dimers, which cost 706 s
// of a 6,527 s run over six rounds. The same information folded to multimer level is a
// few hundred entries per multimer with the best chain hit score attached, so the rounds
// read a fraction of the data and do no lookups at all.
//
// A neighbour is written as the first chain key of its multimer, which keeps the output a
// normal prefilter result that the expansion can consume unchanged.

int multimerneighbors(int argc, const char **argv, const Command &command) {
    LocalParameters &par = LocalParameters::getLocalInstance();
    par.parseParameters(argc, argv, command, true, 0, MMseqsParameter::COMMAND_ALIGN);

    DBReader<unsigned int> alnDbr(par.db3.c_str(), par.db3Index.c_str(), par.threads,
                                  DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA);
    alnDbr.open(DBReader<unsigned int>::LINEAR_ACCCESS);
    DBWriter resultWriter(par.db4.c_str(), par.db4Index.c_str(),
                          static_cast<unsigned int>(par.threads), par.compressed,
                          Parameters::DBTYPE_PREFILTER_RES);
    resultWriter.open();

    const bool touch = par.preloadMode != Parameters::PRELOAD_MODE_MMAP;
    IndexReader qDbr(par.db1, par.threads, IndexReader::SEQUENCES,
                     touch ? IndexReader::PRELOAD_INDEX : 0, DBReader<unsigned int>::USE_INDEX);
    std::vector<unsigned int> qComplexIndices;
    chainKeyToComplexId_t chainKeyToComplexIdMap;
    complexIdToChainKeys_t complexIdToChainKeysMap;
    getKeyToIdMapIdToKeysMapIdVec(qDbr, par.db1 + ".lookup", chainKeyToComplexIdMap,
                                  complexIdToChainKeysMap, qComplexIndices);

    const size_t cap = (size_t) par.maxResListLen;
    Debug(Debug::INFO) << "Folding the chain prefilter to multimer level, at most "
                       << cap << " neighbours each\n";
    Debug::Progress progress(qComplexIndices.size());
#pragma omp parallel
    {
        unsigned int thread_idx = 0;
#ifdef OPENMP
        thread_idx = static_cast<unsigned int>(omp_get_thread_num());
#endif
        std::map<unsigned int, int> best;
        std::vector<std::pair<int, unsigned int> > ranked;
        resultToWrite_t result;
#pragma omp for schedule(dynamic, 1)
        for (size_t qCompIdx = 0; qCompIdx < qComplexIndices.size(); qCompIdx++) {
            const unsigned int qComplexId = qComplexIndices[qCompIdx];
            const std::vector<unsigned int> &qChainKeys = complexIdToChainKeysMap.at(qComplexId);
            best.clear();
            for (size_t i = 0; i < qChainKeys.size(); i++) {
                const unsigned int qId = alnDbr.getId(qChainKeys[i]);
                if (qId == NOT_AVAILABLE_CHAIN_KEY) {
                    continue;
                }
                char *data = alnDbr.getData(qId, thread_idx);
                while (*data != '\0') {
                    char keyBuffer[255 + 1];
                    Util::parseKey(data, keyBuffer);
                    const unsigned int tChainKey = (unsigned int) strtoul(keyBuffer, NULL, 10);
                    chainKeyToComplexId_t::const_iterator c = chainKeyToComplexIdMap.find(tChainKey);
                    if (c == chainKeyToComplexIdMap.end() || c->second == qComplexId) {
                        data = Util::skipLine(data);
                        continue;
                    }
                    const char *p = data;
                    while (*p != '\t' && *p != '\n' && *p != '\0') {
                        p++;
                    }
                    const int score = (*p == '\t') ? (int) strtol(p + 1, NULL, 10) : 0;
                    std::map<unsigned int, int>::iterator f = best.find(c->second);
                    if (f == best.end()) {
                        best[c->second] = score;
                    } else if (score > f->second) {
                        f->second = score;
                    }
                    data = Util::skipLine(data);
                }
            }
            ranked.clear();
            for (std::map<unsigned int, int>::const_iterator it = best.begin(); it != best.end(); ++it) {
                ranked.emplace_back(-it->second, it->first);
            }
            if (cap > 0 && ranked.size() > cap) {
                std::partial_sort(ranked.begin(), ranked.begin() + cap, ranked.end());
                ranked.resize(cap);
            }
            std::sort(ranked.begin(), ranked.end());
            result.clear();
            for (size_t r = 0; r < ranked.size(); r++) {
                const std::vector<unsigned int> &tk = complexIdToChainKeysMap.at(ranked[r].second);
                result.append(SSTR(tk[0]));
                result.push_back('\t');
                result.append(SSTR(-ranked[r].first));
                result.append("\t0\n");
            }
            for (size_t i = 0; i < qChainKeys.size(); i++) {
                resultWriter.writeData(result.c_str(), result.length(), qChainKeys[i], thread_idx);
            }
            progress.updateProgress();
        }
    }
    resultWriter.close(false);
    alnDbr.close();
    return EXIT_SUCCESS;
}
