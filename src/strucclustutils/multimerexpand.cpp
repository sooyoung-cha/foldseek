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

// The expansion round of the linear clustering workflow.
//
// multimerprefilter with --multimer-prefilter-topk hands a constant number of neighbours
// per multimer to the alignment, which keeps the stage linear but caps how large a cluster
// can grow: set-cover can only absorb what the centre is joined to, and in a top-k graph
// no node has enough edges to head a cluster of a thousand members. Measured on 57,448
// CATH sampled dimers the largest cluster fell from 788 to 184 at k=50.
//
// This module wins that back without going back to all against all. Given the clustering
// of the previous round it emits, for every multimer, the centres of the clusters its own
// candidates ended up in. Those are the far members of a big cluster meeting its centre
// directly, and there is a constant number of them per multimer, so the round is linear
// too. Repeating it took the largest cluster from 184 to 1,006 and the co-clustered pairs
// to 95% of what multimercluster finds, for 11.9 extra alignments per multimer.
//
// --multimer-expand-mode 1 is the same idea one level up: every centre against the other
// centres it can see directly. Widening it to look through the neighbourhoods of all of
// its members was tried and reverted: it doubled the run, 1,539 s to 3,233 s, and left
// rep-rep where it was, 3.31% to 3.35%. The representatives that stay apart are not a
// candidate problem. 91% of them were aligned and scored already, and their chain TM sits
// on the threshold, median 0.711, so a hundredth either way decides them.

static const int EXPAND_TO_CENTRES = 0;
static const int EXPAND_BETWEEN_CENTRES = 1;

int multimerexpand(int argc, const char **argv, const Command &command) {
    LocalParameters &par = LocalParameters::getLocalInstance();
    par.parseParameters(argc, argv, command, true, 0, MMseqsParameter::COMMAND_ALIGN);

    DBReader<unsigned int> alnDbr(par.db3.c_str(), par.db3Index.c_str(), par.threads,
                                  DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA);
    alnDbr.open(DBReader<unsigned int>::LINEAR_ACCCESS);
    DBReader<unsigned int> cluDbr(par.db4.c_str(), par.db4Index.c_str(), par.threads,
                                  DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA);
    cluDbr.open(DBReader<unsigned int>::LINEAR_ACCCESS);

    // what the earlier rounds already handed to the alignment, so a pair is never
    // aligned or scored twice and clust never sees the same edge more than once
    DBReader<unsigned int> doneDbr(par.db5.c_str(), par.db5Index.c_str(), par.threads,
                                   DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA);
    doneDbr.open(DBReader<unsigned int>::LINEAR_ACCCESS);

    DBWriter resultWriter(par.db6.c_str(), par.db6Index.c_str(),
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

    // complex -> the centre of the cluster it was put in, by the previous round
    std::map<unsigned int, unsigned int> centreOf;
    for (size_t id = 0; id < cluDbr.getSize(); id++) {
        const unsigned int centre = cluDbr.getDbKey(id);
        char *data = cluDbr.getData(id, 0);
        while (*data != '\0') {
            char keyBuffer[255 + 1];
            Util::parseKey(data, keyBuffer);
            centreOf[(unsigned int) strtoul(keyBuffer, NULL, 10)] = centre;
            data = Util::skipLine(data);
        }
    }
    const int mode = par.multimerExpandMode;
    const size_t topK = (size_t) par.multimerExpandTopk;
    Debug(Debug::INFO) << (mode == EXPAND_BETWEEN_CENTRES
                           ? "Expanding between cluster centres"
                           : "Expanding every multimer to its neighbours' cluster centres")
                       << ", " << centreOf.size() << " multimers clustered\n";

    Debug::Progress progress(qComplexIndices.size());
#pragma omp parallel
    {
        unsigned int thread_idx = 0;
#ifdef OPENMP
        thread_idx = static_cast<unsigned int>(omp_get_thread_num());
#endif
        std::map<unsigned int, int> targetCentreScore;
        std::vector<unsigned int> alreadyDoneCx;
        std::vector<std::pair<int, unsigned int> > ranked;
        resultToWrite_t result;
#pragma omp for schedule(dynamic, 1)
        for (size_t qCompIdx = 0; qCompIdx < qComplexIndices.size(); qCompIdx++) {
            const unsigned int qComplexId = qComplexIndices[qCompIdx];
            const std::vector<unsigned int> &qChainKeys = complexIdToChainKeysMap.at(qComplexId);
            std::map<unsigned int, unsigned int>::const_iterator qc = centreOf.find(qComplexId);
            const bool qIsCentre = (qc != centreOf.end() && qc->second == qComplexId);
            if (mode == EXPAND_BETWEEN_CENTRES && qIsCentre == false) {
                for (size_t i = 0; i < qChainKeys.size(); i++) {
                    resultWriter.writeData("", 0, qChainKeys[i], thread_idx);
                }
                progress.updateProgress();
                continue;
            }
            const unsigned int qCentre = (qc == centreOf.end()) ? qComplexId : qc->second;

            targetCentreScore.clear();
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
                    chainKeyToComplexId_t::const_iterator cIt = chainKeyToComplexIdMap.find(tChainKey);
                    if (cIt == chainKeyToComplexIdMap.end()) {
                        data = Util::skipLine(data);
                        continue;
                    }
                    std::map<unsigned int, unsigned int>::const_iterator tc = centreOf.find(cIt->second);
                    if (tc == centreOf.end()) {
                        data = Util::skipLine(data);
                        continue;
                    }
                    const unsigned int tCentre = tc->second;
                    if (tCentre == qCentre || tCentre == qComplexId) {
                        data = Util::skipLine(data);
                        continue;
                    }
                    if (mode == EXPAND_BETWEEN_CENTRES && cIt->second != tCentre) {
                        // only centres talk to centres at this level
                        data = Util::skipLine(data);
                        continue;
                    }
                    const char *p = data;
                    while (*p != '\t' && *p != '\n' && *p != '\0') {
                        p++;
                    }
                    const int score = (*p == '\t') ? (int) strtol(p + 1, NULL, 10) : 0;
                    std::map<unsigned int, int>::iterator f = targetCentreScore.find(tCentre);
                    if (f == targetCentreScore.end()) {
                        targetCentreScore[tCentre] = score;
                    } else if (score > f->second) {
                        f->second = score;
                    }
                    data = Util::skipLine(data);
                }
            }

            ranked.clear();
            for (std::map<unsigned int, int>::const_iterator it = targetCentreScore.begin();
                 it != targetCentreScore.end(); ++it) {
                ranked.emplace_back(-it->second, it->first);
            }
            if (topK > 0 && ranked.size() > topK) {
                std::partial_sort(ranked.begin(), ranked.begin() + topK, ranked.end());
                ranked.resize(topK);
            }
            // a multimer pair counts as done when either side already listed the other,
            // because the earlier round emitted it in one direction only
            alreadyDoneCx.clear();
            for (size_t i = 0; i < qChainKeys.size(); i++) {
                const unsigned int dId = doneDbr.getId(qChainKeys[i]);
                if (dId == NOT_AVAILABLE_CHAIN_KEY) {
                    continue;
                }
                char *data = doneDbr.getData(dId, thread_idx);
                while (*data != '\0') {
                    char keyBuffer[255 + 1];
                    Util::parseKey(data, keyBuffer);
                    chainKeyToComplexId_t::const_iterator c =
                        chainKeyToComplexIdMap.find((unsigned int) strtoul(keyBuffer, NULL, 10));
                    if (c != chainKeyToComplexIdMap.end()) {
                        alreadyDoneCx.push_back(c->second);
                    }
                    data = Util::skipLine(data);
                }
            }
            SORT_SERIAL(alreadyDoneCx.begin(), alreadyDoneCx.end());
            alreadyDoneCx.erase(std::unique(alreadyDoneCx.begin(), alreadyDoneCx.end()), alreadyDoneCx.end());

            std::vector<unsigned int> targets;
            for (size_t r = 0; r < ranked.size(); r++) {
                const unsigned int tCx = ranked[r].second;
                if (std::binary_search(alreadyDoneCx.begin(), alreadyDoneCx.end(), tCx)) {
                    continue;
                }
                const std::vector<unsigned int> &tChainKeys = complexIdToChainKeysMap.at(tCx);
                // the other direction: did that multimer already list ours?
                bool done = false;
                for (size_t j = 0; j < tChainKeys.size() && done == false; j++) {
                    const unsigned int dId = doneDbr.getId(tChainKeys[j]);
                    if (dId == NOT_AVAILABLE_CHAIN_KEY) {
                        continue;
                    }
                    char *data = doneDbr.getData(dId, thread_idx);
                    while (*data != '\0') {
                        char keyBuffer[255 + 1];
                        Util::parseKey(data, keyBuffer);
                        chainKeyToComplexId_t::const_iterator c =
                            chainKeyToComplexIdMap.find((unsigned int) strtoul(keyBuffer, NULL, 10));
                        if (c != chainKeyToComplexIdMap.end() && c->second == qComplexId) {
                            done = true;
                            break;
                        }
                        data = Util::skipLine(data);
                    }
                }
                if (done) {
                    continue;
                }
                targets.insert(targets.end(), tChainKeys.begin(), tChainKeys.end());
            }
            SORT_SERIAL(targets.begin(), targets.end());
            targets.erase(std::unique(targets.begin(), targets.end()), targets.end());
            result.clear();
            for (size_t t = 0; t < targets.size(); t++) {
                result.append(SSTR(targets[t]));
                result.append("\t0\t0\n");
            }
            for (size_t i = 0; i < qChainKeys.size(); i++) {
                resultWriter.writeData(result.c_str(), result.length(), qChainKeys[i], thread_idx);
            }
            progress.updateProgress();
        }
    }
    resultWriter.close(false);
    doneDbr.close();
    cluDbr.close();
    alnDbr.close();
    return EXIT_SUCCESS;
}
