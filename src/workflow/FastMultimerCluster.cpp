#include <cassert>

#include "FileUtil.h"
#include "CommandCaller.h"
#include "Util.h"
#include "LocalParameters.h"
#include "Debug.h"
#include "PrefilteringIndexReader.h"
#include <limits>
#include "fastmultimercluster.sh.h"

void setFastMultimerClusterDefaults(LocalParameters *p) {
    // same prefilter settings the structure search workflow uses, the chain level
    // prefilter here has to be as sensitive as the one it replaces
    p->kmerSize = 0;
    p->sensitivity = 9.5;
    p->maxResListLen = 1000;
    p->filtMultTmThr = 0.7;
    p->filtChainTmThr = 0.7;
    p->filtInterfaceLddtThr = 0.3;
    p->removeTmpFiles = true;
}

void mustsetFastMultimerCluster(LocalParameters *p) {
    p->clusteringSetMode = 1;
    if (p->filtMultTmThr + p->filtChainTmThr + p->filtInterfaceLddtThr == 0) {
        p->filtMultTmThr = 0.0001;
    }
    p->PARAM_K.wasSet = true;
    p->PARAM_S.wasSet = true;
    p->PARAM_MAX_SEQS.wasSet = true;
    p->PARAM_CLUSTER_SET_MODE.wasSet = true;
    p->PARAM_REMOVE_TMP_FILES.wasSet = true;
    p->PARAM_MULTIMER_TM_THRESHOLD.wasSet = true;
    p->PARAM_CHAIN_TM_THRESHOLD.wasSet = true;
    p->PARAM_INTERFACE_LDDT_THRESHOLD.wasSet = true;
}

int fastmultimercluster(int argc, const char **argv, const Command &command) {
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

    setFastMultimerClusterDefaults(&par);
    par.parseParameters(argc, argv, command, true, Parameters::PARSE_VARIADIC, 0);
    mustsetFastMultimerCluster(&par);

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

    // Chain level prefilter on the 3Di sequences, the same call the structure
    // search makes internally. Its hits feed the candidate filter directly.
    // --cov-mode/-c describe the multimer coverage scoremultimer applies later; the
    // chain prefilter runs uncovered, exactly as the structure search does
    float origCovThr = par.covThr;
    par.covThr = 0.0;
    par.compBiasCorrectionScale = 0.15;
    cmd.addVariable("PREFILTER_PAR", par.createParameterString(par.prefilter).c_str());
    double origPrefEvalThr = par.evalThr;
    par.evalThr = std::numeric_limits<double>::max();
    cmd.addVariable("UNGAPPEDPREFILTER_PAR", par.createParameterString(par.ungappedprefilter).c_str());
    par.evalThr = origPrefEvalThr;
    par.compBiasCorrectionScale = 1.0;
    par.covThr = origCovThr;

    const bool isIndex = PrefilteringIndexReader::searchForIndex(par.db1).empty() == false;
    cmd.addVariable("INDEXEXT", isIndex ? ".idx" : NULL);

    // GPU can only use the ungapped prefilter
    if (par.gpu == 1 && par.PARAM_PREF_MODE.wasSet == false) {
        par.prefMode = Parameters::PREF_MODE_UNGAPPED;
    }

    switch(par.prefMode){
        case LocalParameters::PREF_MODE_KMER:
            cmd.addVariable("PREFMODE", "KMER");
            break;
        case LocalParameters::PREF_MODE_UNGAPPED:
            cmd.addVariable("PREFMODE", "UNGAPPED");
            break;
        case LocalParameters::PREF_MODE_EXHAUSTIVE:
            cmd.addVariable("PREFMODE", "EXHAUSTIVE");
            break;
    }
    if(par.exhaustiveSearch){
        cmd.addVariable("PREFMODE", "EXHAUSTIVE");
    }
    cmd.addVariable("RUNNER", par.runner.c_str());

    cmd.addVariable("MULTIMERPREFILTER_PAR", par.createParameterString(par.multimerprefilter).c_str());

    // scoremultimer needs the backtrace of every candidate chain pair, and the
    // chain pairs come from the composition filter rather than from an E-value,
    // so the alignment must not throw them away again.
    double origEvalThr = par.evalThr;
    bool origAddBacktrace = par.addBacktrace;
    bool origAddBacktraceWasSet = par.PARAM_ADD_BACKTRACE.wasSet;
    float origAlnCovThr = par.covThr;
    par.evalThr = par.eValueThrExpandMultimer;
    par.addBacktrace = true;
    par.PARAM_ADD_BACKTRACE.wasSet = true;
    // the multimer coverage belongs to scoremultimer, the chain alignment must not
    // drop candidate pairs on it
    par.covThr = 0.0;
    cmd.addVariable("MULTIMERALIGN_PAR", par.createParameterString(par.structurealign).c_str());
    par.evalThr = origEvalThr;
    par.addBacktrace = origAddBacktrace;
    par.PARAM_ADD_BACKTRACE.wasSet = origAddBacktraceWasSet;
    par.covThr = origAlnCovThr;
    cmd.addVariable("SCOREMULTIMER_PAR", par.createParameterString(par.scoremultimer).c_str());
    cmd.addVariable("CLUSTER_PAR", par.createParameterString(par.clust).c_str());
    cmd.addVariable("REMOVE_TMP", par.removeTmpFiles ? "TRUE" : NULL);
    cmd.addVariable("VERBOSITY_PAR", par.createParameterString(par.onlyverbosity).c_str());

    std::string program = tmpDir + "/fastmultimercluster.sh";
    FileUtil::writeFile(program, fastmultimercluster_sh, fastmultimercluster_sh_len);
    cmd.execProgram(program.c_str(), par.filenames);

    // Should never get here
    assert(false);
    return EXIT_FAILURE;
}
