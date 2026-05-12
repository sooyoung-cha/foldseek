#include <cassert>
#include <algorithm>
#include <sstream>
#include <string>
#include <vector>

#include "FileUtil.h"
#include "CommandCaller.h"
#include "Util.h"
#include "LocalParameters.h"
#include "Debug.h"
#include "multimercluster.sh.h"
#define multimercluster_fast_sh_len multimercluster_fast_sh_len_dispatch
#include "multimercluster_fast.sh.h"
#undef multimercluster_fast_sh_len

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

void setMultimerClusterDefaults(LocalParameters *p) {
    p->tmScoreThr = 0.65; // TODO
    p->filtChainTmThr = 0.001; // TODO
    p->filtInterfaceLddtThr = 0.5; // TODO
}
   

void mustsetMultimerCluster(LocalParameters *p) { 
    p->clusteringSetMode = 1;
}

static bool shouldUseFastMultimerCluster(const float chainTmThreshold) {
    return chainTmThreshold >= 0.5f;
}

int multimercluster(int argc, const char **argv, const Command &command) {
    LocalParameters &par = LocalParameters::getLocalInstance();
    par.PARAM_ADD_BACKTRACE.addCategory(MMseqsParameter::COMMAND_EXPERT); //align
    par.PARAM_MAX_SEQS.addCategory(MMseqsParameter::COMMAND_EXPERT); //prefilter
    par.PARAM_MAX_REJECTED.addCategory(MMseqsParameter::COMMAND_EXPERT); //align
    par.PARAM_MAX_ACCEPT.addCategory(MMseqsParameter::COMMAND_EXPERT);  //align
    par.PARAM_ZDROP.addCategory(MMseqsParameter::COMMAND_EXPERT); //align
    for (size_t i = 0; i < par.createdb.size(); i++){
        par.createdb[i]->addCategory(MMseqsParameter::COMMAND_EXPERT);
    }
    par.PARAM_COMPRESSED.removeCategory(MMseqsParameter::COMMAND_EXPERT);
    par.PARAM_THREADS.removeCategory(MMseqsParameter::COMMAND_EXPERT);
    par.PARAM_V.removeCategory(MMseqsParameter::COMMAND_EXPERT);

    setMultimerClusterDefaults(&par);
    par.parseParameters(argc, argv, command, true, Parameters::PARSE_VARIADIC, 0);
    mustsetMultimerCluster(&par);
    const bool useFastMultimerCluster = shouldUseFastMultimerCluster(par.filtChainTmThr);

    std::string tmpDir = par.filenames.back();
    std::string hash = SSTR(par.hashParameter(command.databases, par.filenames, *command.params));
    if (par.reuseLatest) {
        hash = FileUtil::getHashFromSymLink(tmpDir + "/latest");
    }
    tmpDir = FileUtil::createTemporaryDirectory(tmpDir, hash);
    par.filenames.pop_back();

    CommandCaller cmd;
    std::cout<<tmpDir.c_str()<<std::endl;
    cmd.addVariable("TMP_PATH", tmpDir.c_str());
    cmd.addVariable("RESULT", par.filenames.back().c_str());
    par.filenames.pop_back();
    cmd.addVariable("INPUT", par.filenames.back().c_str());
    par.filenames.pop_back();

    std::string program;
    if (useFastMultimerCluster) {
        Debug(Debug::INFO) << "chain-tm-threshold " << par.filtChainTmThr
                           << " is at least 0.5, routing multimercluster to multimercluster_fast\n";
        const std::string chainClusterThresholdString = SSTR(par.filtChainTmThr);
        const std::string chainCoverageThresholdString = SSTR(par.covThr);
        const std::string chainCoverageModeString = SSTR(par.covMode);
        std::string chainClusterPar = par.createParameterString(par.structureclusterworkflow, true);
        chainClusterPar = stripOptionAndValue(chainClusterPar, "-c");
        chainClusterPar = stripOptionAndValue(chainClusterPar, "--cov-mode");
        chainClusterPar = stripOptionAndValue(chainClusterPar, "--remove-tmp-files");
        std::string structureAlignPar = par.createParameterString(par.structurealign);
        structureAlignPar = stripOptionAndValue(structureAlignPar, "--tmscore-threshold");
        structureAlignPar.append(" --tmscore-threshold ");
        structureAlignPar.append(SSTR(par.filtChainTmThr));
        cmd.addVariable("CHAINCLUSTER_PAR", chainClusterPar.c_str());
        cmd.addVariable("CHAIN_TM_THRESHOLD", chainClusterThresholdString.c_str());
        cmd.addVariable("CHAIN_CLUSTER_C", chainCoverageThresholdString.c_str());
        cmd.addVariable("CHAIN_CLUSTER_COV_MODE", chainCoverageModeString.c_str());
        cmd.addVariable("CHAINMULTIMERPREFILTER_PAR", par.createParameterString(par.chainmultimerprefilter).c_str());
        cmd.addVariable("STRUCTUREALIGN_PAR", structureAlignPar.c_str());
        cmd.addVariable("SCOREMULTIMER_PAR", par.createParameterString(par.scoremultimer).c_str());
        cmd.addVariable("CLUSTER_PAR", par.createParameterString(par.clust).c_str());
        cmd.addVariable("REMOVE_TMP", par.removeTmpFiles ? "TRUE" : NULL);
        cmd.addVariable("VERBOSITY_PAR", par.createParameterString(par.onlyverbosity).c_str());
        program = tmpDir + "/multimercluster_fast.sh";
        FileUtil::writeFile(program, multimercluster_fast_sh, multimercluster_fast_sh_len_dispatch);
    } else {
        cmd.addVariable("MULTIMERSEARCH_PAR", par.createParameterString(par.multimersearchworkflow, true).c_str());
        cmd.addVariable("CLUSTER_PAR", par.createParameterString(par.clust).c_str());
        cmd.addVariable("REMOVE_TMP", par.removeTmpFiles ? "TRUE" : NULL);
        cmd.addVariable("VERBOSITY_PAR", par.createParameterString(par.onlyverbosity).c_str());
        program = tmpDir + "/multimercluster.sh";
        FileUtil::writeFile(program, multimercluster_sh, multimercluster_sh_len);
    }
    cmd.execProgram(program.c_str(), par.filenames);

    // Should never get here
    assert(false);
    return EXIT_FAILURE;
}
