#include "DBReader.h"
#include "DBWriter.h"
#include "Debug.h"
#include "Util.h"
#include "LocalParameters.h"
#include "StructureUtil.h"
#include "Coordinate16.h"
#include "MultimerUtil.h"
#include "LDDT.h"
#include "tmalign/basic_fun.h"
#include "tmalign/TMalign.h"

#include <limits>
#include <map>
#include <set>
#include <unordered_set>
#include <vector>

#ifdef OPENMP
#include <omp.h>
#endif

#define INTERFACE_THRESHOLD 8
#define MIN_SEED_PROXIMITY_THRESHOLD 0.75
#define SEED_PROXIMITY_SCALE 1.5

struct PairAlignment {
    unsigned int qChainKey;
    unsigned int dbChainKey;
    unsigned int dbComplexId;
    Matcher::result_t res;
    TMaligner::TMscoreResult tmResult;
    std::string line;
    std::vector<float> qCaXVec;
    std::vector<float> qCaYVec;
    std::vector<float> qCaZVec;
    std::vector<float> dbCaXVec;
    std::vector<float> dbCaYVec;
    std::vector<float> dbCaZVec;
};

struct BestSeedResult {
    BestSeedResult() : valid(false),
                       rankingScore(-std::numeric_limits<double>::infinity()),
                       qTmScore(0.0),
                       tTmScore(0.0),
                       qCov(0.0),
                       tCov(0.0),
                       interfaceLddt(0.0),
                       assId(0) {}

    bool valid;
    double rankingScore;
    double qTmScore;
    double tTmScore;
    double qCov;
    double tCov;
    double interfaceLddt;
    unsigned int assId;
    std::string uString;
    std::string tString;
    std::string qChainTms;
    std::string tChainTms;
    std::string filterResult;
    std::vector<size_t> lineIndices;
};

static void appendPointWithTransform(
    const TMaligner::TMscoreResult &tm,
    float x, float y, float z,
    std::vector<float> &outX,
    std::vector<float> &outY,
    std::vector<float> &outZ
) {
    outX.emplace_back(tm.t[0] + x * tm.u[0][0] + y * tm.u[0][1] + z * tm.u[0][2]);
    outY.emplace_back(tm.t[1] + x * tm.u[1][0] + y * tm.u[1][1] + z * tm.u[1][2]);
    outZ.emplace_back(tm.t[2] + x * tm.u[2][0] + y * tm.u[2][1] + z * tm.u[2][2]);
}

static void appendTMString(std::string &result, const std::vector<double> &values) {
    for (size_t i = 0; i < values.size(); ++i) {
        result.append(SSTR(values[i]));
        if (i + 1 < values.size()) {
            result.push_back(',');
        }
    }
}

static void setTransformStrings(const TMaligner::TMscoreResult &tm, std::string &uString, std::string &tString) {
    const char sep = ',';
    uString.clear();
    tString.clear();
    tString.append(std::to_string(tm.t[0]) + sep + std::to_string(tm.t[1]) + sep + std::to_string(tm.t[2]));
    uString.append(std::to_string(tm.u[0][0]) + sep + std::to_string(tm.u[0][1]) + sep + std::to_string(tm.u[0][2]) + sep);
    uString.append(std::to_string(tm.u[1][0]) + sep + std::to_string(tm.u[1][1]) + sep + std::to_string(tm.u[1][2]) + sep);
    uString.append(std::to_string(tm.u[2][0]) + sep + std::to_string(tm.u[2][1]) + sep + std::to_string(tm.u[2][2]));
}

static double getTransformDistance(const TMaligner::TMscoreResult &lhs, const TMaligner::TMscoreResult &rhs) {
    double dist = 0.0;
    for (size_t i = 0; i < 3; ++i) {
        for (size_t j = 0; j < 3; ++j) {
            const double diff = static_cast<double>(lhs.u[i][j]) - static_cast<double>(rhs.u[i][j]);
            dist += diff * diff;
        }
    }
    for (size_t i = 0; i < 3; ++i) {
        const double diff = static_cast<double>(lhs.t[i]) - static_cast<double>(rhs.t[i]);
        dist += diff * diff;
    }
    return std::sqrt(dist);
}

static double computeSeedProximityThreshold(const std::vector<PairAlignment> &pairAlignments) {
    if (pairAlignments.size() <= 2) {
        return 0.0;
    }

    std::vector<double> nearestDistances;
    nearestDistances.reserve(pairAlignments.size());
    for (size_t i = 0; i < pairAlignments.size(); ++i) {
        double nearest = std::numeric_limits<double>::infinity();
        for (size_t j = 0; j < pairAlignments.size(); ++j) {
            if (i == j) {
                continue;
            }
            nearest = std::min(nearest, getTransformDistance(pairAlignments[i].tmResult, pairAlignments[j].tmResult));
        }
        if (std::isfinite(nearest)) {
            nearestDistances.emplace_back(nearest);
        }
    }

    if (nearestDistances.empty()) {
        return 0.0;
    }
    const size_t medianIdx = nearestDistances.size() / 2;
    std::nth_element(nearestDistances.begin(), nearestDistances.begin() + medianIdx, nearestDistances.end());
    return std::max(MIN_SEED_PROXIMITY_THRESHOLD, nearestDistances[medianIdx] * SEED_PROXIMITY_SCALE);
}

static std::vector<size_t> selectRepresentativeSeeds(const std::vector<PairAlignment> &pairAlignments) {
    std::vector<size_t> representatives;
    if (pairAlignments.empty()) {
        return representatives;
    }
    if (pairAlignments.size() <= 2) {
        representatives.reserve(pairAlignments.size());
        for (size_t i = 0; i < pairAlignments.size(); ++i) {
            representatives.emplace_back(i);
        }
        return representatives;
    }

    const double threshold = computeSeedProximityThreshold(pairAlignments);
    if (threshold <= 0.0) {
        representatives.reserve(pairAlignments.size());
        for (size_t i = 0; i < pairAlignments.size(); ++i) {
            representatives.emplace_back(i);
        }
        return representatives;
    }

    const size_t seedCount = pairAlignments.size();
    std::vector<std::vector<size_t> > neighbors(seedCount);
    for (size_t i = 0; i < seedCount; ++i) {
        neighbors[i].emplace_back(i);
        for (size_t j = i + 1; j < seedCount; ++j) {
            const double dist = getTransformDistance(pairAlignments[i].tmResult, pairAlignments[j].tmResult);
            if (dist <= threshold) {
                neighbors[i].emplace_back(j);
                neighbors[j].emplace_back(i);
            }
        }
    }

    std::vector<char> visited(seedCount, 0);
    std::vector<size_t> component;
    std::vector<size_t> stack;
    for (size_t start = 0; start < seedCount; ++start) {
        if (visited[start] != 0) {
            continue;
        }
        component.clear();
        stack.clear();
        stack.emplace_back(start);
        visited[start] = 1;
        while (stack.empty() == false) {
            const size_t curr = stack.back();
            stack.pop_back();
            component.emplace_back(curr);
            for (size_t nIdx = 0; nIdx < neighbors[curr].size(); ++nIdx) {
                const size_t next = neighbors[curr][nIdx];
                if (visited[next] == 0) {
                    visited[next] = 1;
                    stack.emplace_back(next);
                }
            }
        }

        size_t medoid = component[0];
        double bestSum = std::numeric_limits<double>::infinity();
        double bestSeedScore = -std::numeric_limits<double>::infinity();
        for (size_t cIdx = 0; cIdx < component.size(); ++cIdx) {
            const size_t candidate = component[cIdx];
            double distSum = 0.0;
            for (size_t otherIdx = 0; otherIdx < component.size(); ++otherIdx) {
                if (candidate == component[otherIdx]) {
                    continue;
                }
                distSum += getTransformDistance(pairAlignments[candidate].tmResult, pairAlignments[component[otherIdx]].tmResult);
            }
            const double seedScore = pairAlignments[candidate].tmResult.tmscore;
            if (distSum < bestSum || (distSum == bestSum && seedScore > bestSeedScore)) {
                bestSum = distSum;
                bestSeedScore = seedScore;
                medoid = candidate;
            }
        }
        representatives.emplace_back(medoid);
    }

    return representatives;
}

static float getTMd0(float normLen) {
    float D0_MIN = 0.0f;
    float Lnorm = 0.0f;
    float d0 = 0.0f;
    float d0Search = 0.0f;
    parameter_set4final(normLen, D0_MIN, Lnorm, d0, d0Search);
    return d0;
}

static double computeFixedTMscore(
    const std::vector<float> &qX,
    const std::vector<float> &qY,
    const std::vector<float> &qZ,
    const std::vector<float> &dbX,
    const std::vector<float> &dbY,
    const std::vector<float> &dbZ,
    float d0
) {
    if (qX.empty()) {
        return 0.0;
    }
    const double d02 = static_cast<double>(d0) * static_cast<double>(d0);
    double tmScore = 0.0;
    for (size_t i = 0; i < qX.size(); ++i) {
        const double dist = BasicFunction::dist(qX[i], qY[i], qZ[i], dbX[i], dbY[i], dbZ[i]);
        tmScore += 1.0 / (1.0 + dist / d02);
    }
    return tmScore;
}

static double getRankingScore(const LocalParameters &par, double minQChainTm, double minTChainTm) {
    if (par.covMode == Parameters::COV_MODE_TARGET) {
        return minTChainTm;
    }
    if (par.covMode == Parameters::COV_MODE_QUERY) {
        return minQChainTm;
    }
    return std::min(minQChainTm, minTChainTm);
}

static bool passesAssignedChainThreshold(
    const LocalParameters &par,
    size_t alignedQChains,
    size_t totalQChains,
    size_t alignedTChains,
    size_t totalTChains
) {
    const double qRatio = totalQChains == 0 ? 0.0 : static_cast<double>(alignedQChains) / static_cast<double>(totalQChains);
    const double tRatio = totalTChains == 0 ? 0.0 : static_cast<double>(alignedTChains) / static_cast<double>(totalTChains);
    if (par.covMode == Parameters::COV_MODE_TARGET) {
        return tRatio >= par.minAssignedChainsThreshold;
    }
    if (par.covMode == Parameters::COV_MODE_QUERY) {
        return qRatio >= par.minAssignedChainsThreshold;
    }
    return std::min(qRatio, tRatio) >= par.minAssignedChainsThreshold;
}

static bool passesMultimerTmThreshold(const LocalParameters &par, double qTmScore, double tTmScore) {
    if (par.covMode == Parameters::COV_MODE_BIDIRECTIONAL) {
        return qTmScore >= par.filtMultTmThr && tTmScore >= par.filtMultTmThr;
    }
    if (par.covMode == Parameters::COV_MODE_TARGET) {
        return tTmScore >= par.filtMultTmThr;
    }
    return qTmScore >= par.filtMultTmThr;
}

static bool passesChainTmThreshold(const LocalParameters &par, double minQChainTm, double minTChainTm) {
    if (par.filtChainTmThr <= 0) {
        return true;
    }
    if (par.covMode == Parameters::COV_MODE_BIDIRECTIONAL) {
        return minQChainTm >= par.filtChainTmThr && minTChainTm >= par.filtChainTmThr;
    }
    if (par.covMode == Parameters::COV_MODE_TARGET) {
        return minTChainTm >= par.filtChainTmThr;
    }
    return minQChainTm >= par.filtChainTmThr;
}

static void getlookupInfo(
    IndexReader *dbr,
    const std::string &file,
    std::map<unsigned int, unsigned int> &chainKeyToComplexIdLookup,
    std::map<unsigned int, std::vector<unsigned int> > &complexIdToChainKeysLookup,
    std::vector<unsigned int> &complexIdVec
) {
    getKeyToIdMapIdToKeysMapIdVec(*dbr, file, chainKeyToComplexIdLookup, complexIdToChainKeysLookup, complexIdVec);
}

class QueryInterfaceComputer {
public:
    QueryInterfaceComputer(
        const std::vector<unsigned int> &qChainKeys,
        IndexReader *qDbr,
        DBReader<unsigned int> *qCaDbr,
        unsigned int threadIdx
    ) : qChainKeys(qChainKeys), qDbr(qDbr), qCaDbr(qCaDbr), threadIdx(threadIdx), totalResidues(0) {
        interfaceFlags.resize(qChainKeys.size());
    }

    void compute() {
        const float d2 = INTERFACE_THRESHOLD * INTERFACE_THRESHOLD;
        Coordinate16 coords1;
        Coordinate16 coords2;

        for (size_t chainIdx = 0; chainIdx < qChainKeys.size(); ++chainIdx) {
            const unsigned int chainKey = qChainKeys[chainIdx];
            const unsigned int chainDbId = qDbr->sequenceReader->getId(chainKey);
            const size_t chainLen = qDbr->sequenceReader->getSeqLen(chainDbId);
            chainOffsets[chainKey] = totalResidues;
            chainIndexByKey[chainKey] = chainIdx;
            totalResidues += chainLen;
            interfaceFlags[chainIdx].assign(chainLen, 0);

            char *caData = qCaDbr->getData(chainDbId, threadIdx);
            size_t caLength = qCaDbr->getEntryLen(chainDbId);
            float *chainData = coords1.read(caData, chainLen, caLength);

            for (size_t chainIdx2 = 0; chainIdx2 < qChainKeys.size(); ++chainIdx2) {
                if (chainIdx == chainIdx2) {
                    continue;
                }
                const unsigned int chainKey2 = qChainKeys[chainIdx2];
                const unsigned int chainDbId2 = qDbr->sequenceReader->getId(chainKey2);
                const size_t chainLen2 = qDbr->sequenceReader->getSeqLen(chainDbId2);
                char *caData2 = qCaDbr->getData(chainDbId2, threadIdx);
                size_t caLength2 = qCaDbr->getEntryLen(chainDbId2);
                float *chainData2 = coords2.read(caData2, chainLen2, caLength2);

                for (size_t res1 = 0; res1 < chainLen; ++res1) {
                    for (size_t res2 = 0; res2 < chainLen2; ++res2) {
                        const float dist = BasicFunction::dist(
                            chainData[res1], chainData[chainLen + res1], chainData[2 * chainLen + res1],
                            chainData2[res2], chainData2[chainLen2 + res2], chainData2[2 * chainLen2 + res2]
                        );
                        if (dist < d2) {
                            interfaceFlags[chainIdx][res1] = 1;
                            break;
                        }
                    }
                }
            }
        }
    }

    const std::map<unsigned int, size_t> &getChainOffsets() const {
        return chainOffsets;
    }

    const std::map<unsigned int, size_t> &getChainIndexByKey() const {
        return chainIndexByKey;
    }

    const std::vector<std::vector<char> > &getInterfaceFlags() const {
        return interfaceFlags;
    }

    size_t getTotalResidues() const {
        return totalResidues;
    }

private:
    const std::vector<unsigned int> &qChainKeys;
    IndexReader *qDbr;
    DBReader<unsigned int> *qCaDbr;
    unsigned int threadIdx;
    std::vector<std::vector<char> > interfaceFlags;
    std::map<unsigned int, size_t> chainOffsets;
    std::map<unsigned int, size_t> chainIndexByKey;
    size_t totalResidues;
};

static bool buildPairAlignment(
    PairAlignment &pairAlignment,
    unsigned int qChainKey,
    unsigned int dbChainKey,
    unsigned int dbComplexId,
    const std::string &line,
    Coordinate16 &qCoords,
    Coordinate16 &tCoords,
    DBReader<unsigned int> *qCaDbr,
    DBReader<unsigned int> *tCaDbr,
    TMaligner &tmAligner,
    unsigned int threadIdx
) {
    pairAlignment.qChainKey = qChainKey;
    pairAlignment.dbChainKey = dbChainKey;
    pairAlignment.dbComplexId = dbComplexId;
    pairAlignment.line = line;
    pairAlignment.res = Matcher::parseAlignmentRecord(line.c_str());
    if (pairAlignment.res.backtrace.empty()) {
        return false;
    }

    const size_t qCaId = qCaDbr->getId(qChainKey);
    const size_t tCaId = tCaDbr->getId(dbChainKey);
    char *qCaData = qCaDbr->getData(qCaId, threadIdx);
    char *tCaData = tCaDbr->getData(tCaId, threadIdx);
    const size_t qCaLength = qCaDbr->getEntryLen(qCaId);
    const size_t tCaLength = tCaDbr->getEntryLen(tCaId);
    float *queryCaData = qCoords.read(qCaData, pairAlignment.res.qLen, qCaLength);
    float *targetCaData = tCoords.read(tCaData, pairAlignment.res.dbLen, tCaLength);

    tmAligner.initQuery(queryCaData, &queryCaData[pairAlignment.res.qLen], &queryCaData[pairAlignment.res.qLen * 2], NULL, pairAlignment.res.qLen);
    pairAlignment.tmResult = tmAligner.computeTMscore(
        targetCaData,
        &targetCaData[pairAlignment.res.dbLen],
        &targetCaData[pairAlignment.res.dbLen * 2],
        pairAlignment.res.dbLen,
        pairAlignment.res.qStartPos,
        pairAlignment.res.dbStartPos,
        Matcher::uncompressAlignment(pairAlignment.res.backtrace),
        pairAlignment.res.qLen
    );

    unsigned int qPos = pairAlignment.res.qStartPos;
    unsigned int dbPos = pairAlignment.res.dbStartPos;
    for (size_t btPos = 0; btPos < pairAlignment.res.backtrace.size(); ++btPos) {
        const char bt = pairAlignment.res.backtrace[btPos];
        if (bt == 'M') {
            pairAlignment.qCaXVec.emplace_back(queryCaData[qPos]);
            pairAlignment.qCaYVec.emplace_back(queryCaData[pairAlignment.res.qLen + qPos]);
            pairAlignment.qCaZVec.emplace_back(queryCaData[pairAlignment.res.qLen * 2 + qPos]);
            pairAlignment.dbCaXVec.emplace_back(targetCaData[dbPos]);
            pairAlignment.dbCaYVec.emplace_back(targetCaData[pairAlignment.res.dbLen + dbPos]);
            pairAlignment.dbCaZVec.emplace_back(targetCaData[pairAlignment.res.dbLen * 2 + dbPos]);
            ++qPos;
            ++dbPos;
        } else if (bt == 'I') {
            ++qPos;
        } else if (bt == 'D') {
            ++dbPos;
        }
    }
    return !pairAlignment.qCaXVec.empty();
}

static double evaluateInterfaceLddt(
    const std::vector<PairAlignment> &pairAlignments,
    const TMaligner::TMscoreResult &seedTm,
    const std::map<unsigned int, size_t> &qChainIndexByKey,
    const std::vector<std::vector<char> > &qInterfaceFlags
) {
    if (pairAlignments.size() <= 1) {
        return 0.0;
    }

    std::vector<float> qIntVecX;
    std::vector<float> qIntVecY;
    std::vector<float> qIntVecZ;
    std::vector<float> dbIntVecX;
    std::vector<float> dbIntVecY;
    std::vector<float> dbIntVecZ;
    size_t wholeInterfaceLen = 0;

    for (size_t alnIdx = 0; alnIdx < pairAlignments.size(); ++alnIdx) {
        const PairAlignment &pairAlignment = pairAlignments[alnIdx];
        std::map<unsigned int, size_t>::const_iterator chainIndexIt = qChainIndexByKey.find(pairAlignment.qChainKey);
        if (chainIndexIt == qChainIndexByKey.end()) {
            continue;
        }
        const std::vector<char> &interfaceFlag = qInterfaceFlags[chainIndexIt->second];
        wholeInterfaceLen += std::count(interfaceFlag.begin(), interfaceFlag.end(), 1);

        unsigned int qPos = pairAlignment.res.qStartPos;
        unsigned int dbPos = pairAlignment.res.dbStartPos;
        size_t matchIdx = 0;
        for (size_t btPos = 0; btPos < pairAlignment.res.backtrace.size(); ++btPos) {
            const char bt = pairAlignment.res.backtrace[btPos];
            if (bt == 'M') {
                if (qPos < interfaceFlag.size() && interfaceFlag[qPos] != 0) {
                    qIntVecX.emplace_back(pairAlignment.qCaXVec[matchIdx]);
                    qIntVecY.emplace_back(pairAlignment.qCaYVec[matchIdx]);
                    qIntVecZ.emplace_back(pairAlignment.qCaZVec[matchIdx]);
                    appendPointWithTransform(
                        seedTm,
                        pairAlignment.dbCaXVec[matchIdx],
                        pairAlignment.dbCaYVec[matchIdx],
                        pairAlignment.dbCaZVec[matchIdx],
                        dbIntVecX,
                        dbIntVecY,
                        dbIntVecZ
                    );
                }
                ++qPos;
                ++dbPos;
                ++matchIdx;
            } else if (bt == 'I') {
                ++qPos;
            } else if (bt == 'D') {
                ++dbPos;
            }
        }
    }

    if (dbIntVecX.empty() || wholeInterfaceLen == 0) {
        return 0.0;
    }

    std::string backtrace(dbIntVecX.size(), 'M');
    LDDTCalculator lddtCalculator(dbIntVecX.size() + 1, dbIntVecX.size() + 1);
    lddtCalculator.initQuery(dbIntVecX.size(), &qIntVecX[0], &qIntVecY[0], &qIntVecZ[0]);
    LDDTCalculator::LDDTScoreResult lddtRes = lddtCalculator.computeLDDTScore(
        dbIntVecX.size(),
        0,
        0,
        backtrace,
        &dbIntVecX[0],
        &dbIntVecY[0],
        &dbIntVecZ[0]
    );
    return lddtRes.avgLddtScore * lddtRes.scoreLength / wholeInterfaceLen;
}

static bool evaluateSeed(
    const LocalParameters &par,
    const std::vector<PairAlignment> &pairAlignments,
    const TMaligner::TMscoreResult &seedTm,
    size_t qComplexLength,
    size_t dbComplexLength,
    size_t qComplexChainCount,
    size_t dbComplexChainCount,
    const std::map<unsigned int, size_t> &qChainOffsets,
    const std::map<unsigned int, size_t> &dbChainOffsets,
    const std::map<unsigned int, size_t> &qChainIndexByKey,
    const std::vector<std::vector<char> > &qInterfaceFlags,
    BestSeedResult &result
) {
    if (pairAlignments.size() < static_cast<size_t>(par.minAlignedChains)) {
        return false;
    }

    std::vector<char> qCovered(qComplexLength, 0);
    std::vector<char> dbCovered(dbComplexLength, 0);
    std::unordered_set<unsigned int> alignedQChains;
    std::unordered_set<unsigned int> alignedDbChains;
    std::vector<double> qChainTms;
    std::vector<double> dbChainTms;
    double minQChainTm = std::numeric_limits<double>::infinity();
    double minTChainTm = std::numeric_limits<double>::infinity();
    double multimerSum = 0.0;
    const float multimerD0 = getTMd0(static_cast<float>(std::max<size_t>(1, std::min(qComplexLength, dbComplexLength))));

    for (size_t alnIdx = 0; alnIdx < pairAlignments.size(); ++alnIdx) {
        const PairAlignment &pairAlignment = pairAlignments[alnIdx];
        alignedQChains.insert(pairAlignment.qChainKey);
        alignedDbChains.insert(pairAlignment.dbChainKey);
        result.lineIndices.emplace_back(alnIdx);

        std::vector<float> transformedX;
        std::vector<float> transformedY;
        std::vector<float> transformedZ;
        transformedX.reserve(pairAlignment.dbCaXVec.size());
        transformedY.reserve(pairAlignment.dbCaYVec.size());
        transformedZ.reserve(pairAlignment.dbCaZVec.size());
        for (size_t i = 0; i < pairAlignment.dbCaXVec.size(); ++i) {
            appendPointWithTransform(
                seedTm,
                pairAlignment.dbCaXVec[i],
                pairAlignment.dbCaYVec[i],
                pairAlignment.dbCaZVec[i],
                transformedX,
                transformedY,
                transformedZ
            );
        }

        const float chainD0 = getTMd0(static_cast<float>(std::max(1u, std::min(pairAlignment.res.qLen, pairAlignment.res.dbLen))));
        const double rawChainScore = computeFixedTMscore(
            pairAlignment.qCaXVec, pairAlignment.qCaYVec, pairAlignment.qCaZVec,
            transformedX, transformedY, transformedZ,
            chainD0
        );
        const double qChainTm = rawChainScore / pairAlignment.res.qLen;
        const double dbChainTm = rawChainScore / pairAlignment.res.dbLen;
        qChainTms.emplace_back(qChainTm);
        dbChainTms.emplace_back(dbChainTm);
        minQChainTm = std::min(minQChainTm, qChainTm);
        minTChainTm = std::min(minTChainTm, dbChainTm);

        unsigned int qPos = pairAlignment.res.qStartPos;
        unsigned int dbPos = pairAlignment.res.dbStartPos;
        size_t matchIdx = 0;
        const double d02 = static_cast<double>(multimerD0) * static_cast<double>(multimerD0);
        for (size_t btPos = 0; btPos < pairAlignment.res.backtrace.size(); ++btPos) {
            const char bt = pairAlignment.res.backtrace[btPos];
            if (bt == 'M') {
                qCovered[qChainOffsets.at(pairAlignment.qChainKey) + qPos] = 1;
                dbCovered[dbChainOffsets.at(pairAlignment.dbChainKey) + dbPos] = 1;
                const double dist = BasicFunction::dist(
                    pairAlignment.qCaXVec[matchIdx], pairAlignment.qCaYVec[matchIdx], pairAlignment.qCaZVec[matchIdx],
                    transformedX[matchIdx], transformedY[matchIdx], transformedZ[matchIdx]
                );
                multimerSum += 1.0 / (1.0 + dist / d02);
                ++qPos;
                ++dbPos;
                ++matchIdx;
            } else if (bt == 'I') {
                ++qPos;
            } else if (bt == 'D') {
                ++dbPos;
            }
        }
    }

    if (alignedQChains.empty() || alignedDbChains.empty()) {
        return false;
    }
    if (!passesAssignedChainThreshold(par, alignedQChains.size(), qComplexChainCount, alignedDbChains.size(), dbComplexChainCount)) {
        return false;
    }

    result.qCov = static_cast<double>(std::count(qCovered.begin(), qCovered.end(), 1)) / static_cast<double>(qComplexLength);
    result.tCov = static_cast<double>(std::count(dbCovered.begin(), dbCovered.end(), 1)) / static_cast<double>(dbComplexLength);
    if (!Util::hasCoverage(par.covThr, par.covMode, result.qCov, result.tCov)) {
        return false;
    }

    result.qTmScore = multimerSum / static_cast<double>(qComplexLength);
    result.tTmScore = multimerSum / static_cast<double>(dbComplexLength);
    if (!passesMultimerTmThreshold(par, result.qTmScore, result.tTmScore)) {
        return false;
    }
    if (!passesChainTmThreshold(par, minQChainTm, minTChainTm)) {
        return false;
    }

    result.interfaceLddt = evaluateInterfaceLddt(pairAlignments, seedTm, qChainIndexByKey, qInterfaceFlags);
    if (pairAlignments.size() == 1) {
        if (par.filtInterfaceLddtThr > 0) {
            return false;
        }
    } else if (result.interfaceLddt < par.filtInterfaceLddtThr) {
        return false;
    }

    setTransformStrings(seedTm, result.uString, result.tString);
    appendTMString(result.qChainTms, qChainTms);
    appendTMString(result.tChainTms, dbChainTms);
    result.filterResult = SSTR(result.qCov) + "\t" + SSTR(result.tCov) + "\t" + result.qChainTms + "\t" + result.tChainTms + "\t" + SSTR(result.interfaceLddt);
    result.rankingScore = getRankingScore(par, minQChainTm, minTChainTm);
    result.valid = true;
    return true;
}

int scoremultimerfast(int argc, const char **argv, const Command &command) {
    LocalParameters &par = LocalParameters::getLocalInstance();
    par.parseParameters(argc, argv, command, true, 0, MMseqsParameter::COMMAND_ALIGN);

    DBReader<unsigned int> alnDbr(par.db3.c_str(), par.db3Index.c_str(), par.threads, DBReader<unsigned int>::USE_INDEX | DBReader<unsigned int>::USE_DATA);
    alnDbr.open(DBReader<unsigned int>::LINEAR_ACCCESS);

    int dbType = alnDbr.getDbtype();
    uint16_t extended = DBReader<unsigned int>::getExtendedDbtype(dbType);
    bool needSrc = false;
    if (extended & Parameters::DBTYPE_EXTENDED_INDEX_NEED_SRC) {
        needSrc = true;
        dbType = DBReader<unsigned int>::setExtendedDbtype(dbType, Parameters::DBTYPE_EXTENDED_INDEX_NEED_SRC);
    }
    DBWriter resultWriter(par.db4.c_str(), par.db4Index.c_str(), static_cast<unsigned int>(par.threads), par.compressed, dbType);
    resultWriter.open();

    const bool touch = (par.preloadMode != Parameters::PRELOAD_MODE_MMAP);
    std::string t3DiDbrName = StructureUtil::getIndexWithSuffix(par.db2, "_ss");
    bool is3DiIdx = Parameters::isEqualDbtype(FileUtil::parseDbType(t3DiDbrName.c_str()), Parameters::DBTYPE_INDEX_DB);
    IndexReader *t3DiDbr = new IndexReader(
        is3DiIdx ? t3DiDbrName : par.db2,
        par.threads,
        needSrc ? IndexReader::SRC_SEQUENCES : IndexReader::SEQUENCES,
        touch ? IndexReader::PRELOAD_INDEX : 0,
        DBReader<unsigned int>::USE_INDEX | DBReader<unsigned int>::USE_DATA,
        needSrc ? "_seq_ss" : "_ss"
    );
    DBReader<unsigned int> *tCaDbr = new DBReader<unsigned int>(
        needSrc ? (par.db2 + "_seq_ca").c_str() : (par.db2 + "_ca").c_str(),
        needSrc ? (par.db2 + "_seq_ca.index").c_str() : (par.db2 + "_ca.index").c_str(),
        par.threads,
        DBReader<unsigned int>::USE_INDEX | DBReader<unsigned int>::USE_DATA
    );
    tCaDbr->open(DBReader<unsigned int>::NOSORT);

    IndexReader *q3DiDbr = NULL;
    DBReader<unsigned int> *qCaDbr = NULL;
    bool sameDB = false;
    if (par.db1 == par.db2) {
        sameDB = true;
        q3DiDbr = t3DiDbr;
        qCaDbr = tCaDbr;
    } else {
        q3DiDbr = new IndexReader(
            StructureUtil::getIndexWithSuffix(par.db1, "_ss"),
            par.threads,
            IndexReader::SEQUENCES,
            touch ? IndexReader::PRELOAD_INDEX : 0,
            DBReader<unsigned int>::USE_INDEX | DBReader<unsigned int>::USE_DATA
        );
        qCaDbr = new DBReader<unsigned int>((par.db1 + "_ca").c_str(), (par.db1 + "_ca.index").c_str(), par.threads, DBReader<unsigned int>::USE_INDEX | DBReader<unsigned int>::USE_DATA);
        qCaDbr->open(DBReader<unsigned int>::NOSORT);
    }

    std::vector<unsigned int> qComplexIndices;
    std::vector<unsigned int> dbComplexIndices;
    chainKeyToComplexId_t qChainKeyToComplexIdMap, dbChainKeyToComplexIdMap;
    complexIdToChainKeys_t dbComplexIdToChainKeysMap, qComplexIdToChainKeysMap;
    getlookupInfo(q3DiDbr, par.db1 + ".lookup", qChainKeyToComplexIdMap, qComplexIdToChainKeysMap, qComplexIndices);
    if (sameDB) {
        dbChainKeyToComplexIdMap = qChainKeyToComplexIdMap;
        dbComplexIdToChainKeysMap = qComplexIdToChainKeysMap;
        dbComplexIndices = qComplexIndices;
    } else {
        getlookupInfo(t3DiDbr, par.db2 + ".lookup", dbChainKeyToComplexIdMap, dbComplexIdToChainKeysMap, dbComplexIndices);
    }

    Debug::Progress progress(qComplexIndices.size());

#pragma omp parallel
    {
        unsigned int threadIdx = 0;
#ifdef OPENMP
        threadIdx = static_cast<unsigned int>(omp_get_thread_num());
#endif
        Coordinate16 qCoords;
        Coordinate16 tCoords;
        TMaligner tmAligner(std::max(q3DiDbr->sequenceReader->getMaxSeqLen() + 1, t3DiDbr->sequenceReader->getMaxSeqLen() + 1), false, true, false);
        std::map<unsigned int, std::string> resultBuffers;

#pragma omp for schedule(dynamic, 1)
        for (size_t qCompIdx = 0; qCompIdx < qComplexIndices.size(); ++qCompIdx) {
            progress.updateProgress();
            resultBuffers.clear();

            const unsigned int qComplexId = qComplexIndices[qCompIdx];
            const std::vector<unsigned int> &qChainKeys = qComplexIdToChainKeysMap.at(qComplexId);
            if (par.monomerIncludeMode == SKIP_MONOMERS && qChainKeys.size() < MULTIPLE_CHAINED_COMPLEX) {
                for (size_t chainIdx = 0; chainIdx < qChainKeys.size(); ++chainIdx) {
                    resultWriter.writeData("", 0, qChainKeys[chainIdx], threadIdx);
                }
                continue;
            }

            QueryInterfaceComputer interfaceComputer(qChainKeys, q3DiDbr, qCaDbr, threadIdx);
            interfaceComputer.compute();

            std::map<unsigned int, std::vector<PairAlignment> > byDbComplex;
            for (size_t qChainIdx = 0; qChainIdx < qChainKeys.size(); ++qChainIdx) {
                const unsigned int qChainKey = qChainKeys[qChainIdx];
                const unsigned int qDbKey = alnDbr.getId(qChainKey);
                if (qDbKey == NOT_AVAILABLE_CHAIN_KEY) {
                    continue;
                }
                char *data = alnDbr.getData(qDbKey, threadIdx);
                const size_t dataSize = alnDbr.getDataSize();
                char lineBuffer[1024];
                while (*data != '\0') {
                    Util::getLine(data, dataSize, lineBuffer, 1024);
                    Matcher::result_t res = Matcher::parseAlignmentRecord(lineBuffer);
                    if (res.backtrace.empty()) {
                        data = Util::skipLine(data);
                        continue;
                    }
                    const unsigned int dbChainKey = res.dbKey;
                    chainKeyToComplexId_t::const_iterator dbComplexIt = dbChainKeyToComplexIdMap.find(dbChainKey);
                    if (dbComplexIt == dbChainKeyToComplexIdMap.end()) {
                        data = Util::skipLine(data);
                        continue;
                    }
                    const unsigned int dbComplexId = dbComplexIt->second;
                    const std::vector<unsigned int> &dbChainKeys = dbComplexIdToChainKeysMap.at(dbComplexId);
                    if (par.monomerIncludeMode == SKIP_MONOMERS && dbChainKeys.size() < MULTIPLE_CHAINED_COMPLEX) {
                        data = Util::skipLine(data);
                        continue;
                    }

                    PairAlignment pairAlignment;
                    if (buildPairAlignment(pairAlignment, qChainKey, dbChainKey, dbComplexId, lineBuffer, qCoords, tCoords, qCaDbr, tCaDbr, tmAligner, threadIdx)) {
                        byDbComplex[dbComplexId].emplace_back(pairAlignment);
                    }
                    data = Util::skipLine(data);
                }
            }

            unsigned int assId = 0;
            for (std::map<unsigned int, std::vector<PairAlignment> >::const_iterator dbIt = byDbComplex.begin(); dbIt != byDbComplex.end(); ++dbIt) {
                const unsigned int dbComplexId = dbIt->first;
                const std::vector<PairAlignment> &pairAlignments = dbIt->second;
                const std::vector<unsigned int> &dbChainKeys = dbComplexIdToChainKeysMap.at(dbComplexId);

                std::map<unsigned int, size_t> dbChainOffsets;
                size_t dbComplexLength = 0;
                for (size_t dbIdx = 0; dbIdx < dbChainKeys.size(); ++dbIdx) {
                    const unsigned int dbChainKey = dbChainKeys[dbIdx];
                    const size_t dbChainId = t3DiDbr->sequenceReader->getId(dbChainKey);
                    dbChainOffsets[dbChainKey] = dbComplexLength;
                    dbComplexLength += t3DiDbr->sequenceReader->getSeqLen(dbChainId);
                }

                BestSeedResult bestResult;
                const std::vector<size_t> representativeSeeds = selectRepresentativeSeeds(pairAlignments);
                for (size_t repIdx = 0; repIdx < representativeSeeds.size(); ++repIdx) {
                    const size_t seedIdx = representativeSeeds[repIdx];
                    BestSeedResult candidate;
                    if (!evaluateSeed(
                            par,
                            pairAlignments,
                            pairAlignments[seedIdx].tmResult,
                            interfaceComputer.getTotalResidues(),
                            dbComplexLength,
                            qChainKeys.size(),
                            dbChainKeys.size(),
                            interfaceComputer.getChainOffsets(),
                            dbChainOffsets,
                            interfaceComputer.getChainIndexByKey(),
                            interfaceComputer.getInterfaceFlags(),
                            candidate)) {
                        continue;
                    }
                    if (!bestResult.valid ||
                        candidate.rankingScore > bestResult.rankingScore ||
                        (candidate.rankingScore == bestResult.rankingScore && candidate.qTmScore > bestResult.qTmScore) ||
                        (candidate.rankingScore == bestResult.rankingScore && candidate.qTmScore == bestResult.qTmScore && candidate.tTmScore > bestResult.tTmScore)) {
                        bestResult = candidate;
                    }
                }

                if (!bestResult.valid) {
                    continue;
                }
                bestResult.assId = assId++;
                const std::string complexFields = "\t" + SSTR(bestResult.qTmScore) + "\t" + SSTR(bestResult.tTmScore) + "\t" + bestResult.uString + "\t" + bestResult.tString + "\t" + SSTR(bestResult.assId) + "\t" + bestResult.filterResult + "\n";
                for (size_t lineIdx = 0; lineIdx < bestResult.lineIndices.size(); ++lineIdx) {
                    const PairAlignment &pairAlignment = pairAlignments[bestResult.lineIndices[lineIdx]];
                    resultBuffers[pairAlignment.qChainKey].append(pairAlignment.line);
                    resultBuffers[pairAlignment.qChainKey].append(complexFields);
                }
            }

            for (size_t qChainIdx = 0; qChainIdx < qChainKeys.size(); ++qChainIdx) {
                const unsigned int qChainKey = qChainKeys[qChainIdx];
                const std::string &buffer = resultBuffers[qChainKey];
                resultWriter.writeData(buffer.c_str(), buffer.length(), qChainKey, threadIdx);
            }
        }
    }

    alnDbr.close();
    delete t3DiDbr;
    tCaDbr->close();
    delete tCaDbr;
    if (!sameDB) {
        delete q3DiDbr;
        qCaDbr->close();
        delete qCaDbr;
    }
    resultWriter.close(false);
    return EXIT_SUCCESS;
}
