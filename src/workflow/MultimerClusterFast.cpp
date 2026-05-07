#include <cassert>
#include <algorithm>
#include <sstream>
#include <string>
#include <vector>

#include "CommandCaller.h"
#include "Debug.h"
#include "FileUtil.h"
#include "LocalParameters.h"
#include "Util.h"
#include "multimercluster_fast.sh.h"

static std::string stripOptionAndValue(const std::string &params, const std::string &option) {
    std::istringstream input(params);
    std::vector<std::string> tokens;
    std::string token;
    while (input >> token) {
        tokens.push_back(token);
    }

    std::string result;
    for (size_t i = 0; i < tokens.size(); ++i) {
        if (tokens[i] == option) {
            ++i;
            continue;
        }
        if (result.empty() == false) {
            result.append(" ");
        }
        result.append(tokens[i]);
    }
    return result;
}

static bool shouldUseOriginalMultimerCluster(const float chainTmThreshold) {
    return chainTmThreshold < 0.5f;
}

static float mapChainClusterThreshold(const float chainTmThreshold) {
    if (chainTmThreshold < 0.7f) {
        return std::min(1.0f, chainTmThreshold + 0.1f);
    }
    return std::min(1.0f, chainTmThreshold + 0.05f);
}

void setMultimerClusterFastDefaults(LocalParameters *p) {
    p->filtMultTmThr = 0.65;
    p->filtChainTmThr = 0.001;
    p->filtInterfaceLddtThr = 0.5;
}

int multimercluster_fast(int argc, const char **argv, const Command &command) {
    LocalParameters &par = LocalParameters::getLocalInstance();
    par.PARAM_ADD_BACKTRACE.addCategory(MMseqsParameter::COMMAND_EXPERT);
    par.PARAM_MAX_SEQS.addCategory(MMseqsParameter::COMMAND_EXPERT);
    par.PARAM_MAX_REJECTED.addCategory(MMseqsParameter::COMMAND_EXPERT);
    par.PARAM_MAX_ACCEPT.addCategory(MMseqsParameter::COMMAND_EXPERT);
    par.PARAM_ZDROP.addCategory(MMseqsParameter::COMMAND_EXPERT);
    for (size_t i = 0; i < par.createdb.size(); i++) {
        par.createdb[i]->addCategory(MMseqsParameter::COMMAND_EXPERT);
    }
    par.PARAM_COMPRESSED.removeCategory(MMseqsParameter::COMMAND_EXPERT);
    par.PARAM_THREADS.removeCategory(MMseqsParameter::COMMAND_EXPERT);
    par.PARAM_V.removeCategory(MMseqsParameter::COMMAND_EXPERT);

    setMultimerClusterFastDefaults(&par);
    par.parseParameters(argc, argv, command, true, Parameters::PARSE_VARIADIC, 0);
    par.addBacktrace = true;
    par.PARAM_ADD_BACKTRACE.wasSet = true;
    par.clusteringSetMode = 1;
    const bool useOriginalMultimerCluster = shouldUseOriginalMultimerCluster(par.filtChainTmThr);
    const float chainClusterThreshold = mapChainClusterThreshold(par.filtChainTmThr);
    std::string chainTmThreshold = SSTR(chainClusterThreshold);
    std::string chainCoverageThreshold = SSTR(par.covThr);
    std::string chainCoverageMode = SSTR(par.covMode);

    if (useOriginalMultimerCluster) {
        Debug(Debug::INFO) << "chain-tm-threshold " << par.filtChainTmThr
                           << " is below 0.5, falling back to multimercluster\n";
    } else {
        Debug(Debug::INFO) << "Using multimercluster_fast with internal chain cluster TM threshold "
                           << chainClusterThreshold << " (user chain-tm-threshold " << par.filtChainTmThr << ")\n";
    }

    std::string tmpDir = par.filenames.back();
    std::string hash = SSTR(par.hashParameter(command.databases, par.filenames, *command.params));
    if (par.reuseLatest) {
        hash = FileUtil::getHashFromSymLink(tmpDir + "/latest");
    }
    tmpDir = FileUtil::createTemporaryDirectory(tmpDir, hash);
    par.filenames.pop_back();

    CommandCaller cmd;
    cmd.addVariable("TMP_PATH", tmpDir.c_str());
    cmd.addVariable("RESULT", par.filenames.back().c_str());
    par.filenames.pop_back();
    cmd.addVariable("INPUT", par.filenames.back().c_str());
    par.filenames.pop_back();

    std::string chainClusterPar = par.createParameterString(par.structureclusterworkflow, true);
    chainClusterPar = stripOptionAndValue(chainClusterPar, "-c");
    chainClusterPar = stripOptionAndValue(chainClusterPar, "--cov-mode");
    chainClusterPar = stripOptionAndValue(chainClusterPar, "--remove-tmp-files");
    const std::string multimerClusterPar = par.createParameterString(par.multimerclusterworkflow, true);
    cmd.addVariable("USE_ORIGINAL_MULTIMERCLUSTER", useOriginalMultimerCluster ? "TRUE" : NULL);
    cmd.addVariable("MULTIMERCLUSTER_PAR", multimerClusterPar.c_str());
    cmd.addVariable("CHAINCLUSTER_PAR", chainClusterPar.c_str());
    cmd.addVariable("CHAIN_TM_THRESHOLD", chainTmThreshold.c_str());
    cmd.addVariable("CHAIN_CLUSTER_C", chainCoverageThreshold.c_str());
    cmd.addVariable("CHAIN_CLUSTER_COV_MODE", chainCoverageMode.c_str());
    cmd.addVariable("CHAINMULTIMERPREFILTER_PAR", par.createParameterString(par.chainmultimerprefilter).c_str());
    cmd.addVariable("STRUCTUREALIGN_PAR", par.createParameterString(par.structurealign).c_str());
    cmd.addVariable("SCOREMULTIMER_PAR", par.createParameterString(par.scoremultimer).c_str());
    cmd.addVariable("CLUSTER_PAR", par.createParameterString(par.clust).c_str());
    cmd.addVariable("REMOVE_TMP", par.removeTmpFiles ? "TRUE" : NULL);
    cmd.addVariable("VERBOSITY_PAR", par.createParameterString(par.onlyverbosity).c_str());

    std::string program = tmpDir + "/multimercluster_fast.sh";
    FileUtil::writeFile(program, multimercluster_fast_sh, multimercluster_fast_sh_len);
    cmd.execProgram(program.c_str(), par.filenames);

    assert(false);
    return EXIT_FAILURE;
}
