#include <cassert>

#include "FileUtil.h"
#include "CommandCaller.h"
#include "Util.h"
#include "LocalParameters.h"
#include "Debug.h"
#include "PrefilteringIndexReader.h"
#include <limits>
#include "linmultimercluster.sh.h"

void setLinMultimerClusterDefaults(LocalParameters *p) {
    // The chain prefilter does not have to be as sensitive as the one the quadratic
    // workflow needs. That one enumerates every pair within the threshold, so it cannot
    // give up reach; this one keeps the --multimer-prefilter-topk nearest neighbours of
    // each multimer and the expansion rounds take it from there, and the nearest
    // neighbours are easy to find. Measured on 57,448 dimers, the share of multimers that
    // see a true neighbour is 97.6% at -s 9.5 and 97.5% at -s 6.0, for 248 s against 26 s.
    p->kmerSize = 0;
    p->sensitivity = 6.0;
    p->maxResListLen = 1000;
    p->filtMultTmThr = 0.7;
    p->filtChainTmThr = 0.7;
    p->filtInterfaceLddtThr = 0.3;
    p->removeTmpFiles = true;
    // Neighbours kept per multimer in the first round. The work per multimer, and so the
    // whole stage, is a constant times this, whatever the database size.
    p->multimerPrefilterTopk = 50;
    p->multimerExpandTopk = 200;
    p->multimerExpandRounds = 3;
}

void mustsetLinMultimerCluster(LocalParameters *p) {
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

int linmultimercluster(int argc, const char **argv, const Command &command) {
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

    setLinMultimerClusterDefaults(&par);
    par.parseParameters(argc, argv, command, true, Parameters::PARSE_VARIADIC, 0);
    mustsetLinMultimerCluster(&par);

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

    // Rounds of mode 0, which lets a cluster grow past what a top-k graph can hold, then
    // one of mode 1. Measured round by round on 57,448 dimers the clusters settle at
    // 21,648, 19,763, 19,534, 19,484, 19,465 and then stop moving, so the mode 1 rounds
    // after the first only cost time.
    std::string schedule;
    for (int i = 0; i < par.multimerExpandRounds; i++) {
        schedule += (i == 0 ? "0" : " 0");
    }
    schedule += " 1";
    cmd.addVariable("EXPAND_SCHEDULE", schedule.c_str());

    // one string per mode, so the schedule can pick without repeating the flag
    int origMode = par.multimerExpandMode;
    par.multimerExpandMode = 0;
    cmd.addVariable("MULTIMERNEIGHBORS_PAR", par.createParameterString(par.multimerneighbors).c_str());
    cmd.addVariable("MULTIMEREXPAND_PAR0", par.createParameterString(par.multimerexpand).c_str());
    par.multimerExpandMode = 1;
    cmd.addVariable("MULTIMEREXPAND_PAR1", par.createParameterString(par.multimerexpand).c_str());
    par.multimerExpandMode = origMode;
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
    // A chain pair can only reach the chain TM threshold when the aligned region covers
    // that fraction of both chains: TM is the sum over matched residues of a term at most
    // one, divided by the chain length, so TM <= matched/length <= coverage. Under
    // cov-mode 0, where scoremultimer demands the threshold on every chain of both sides,
    // dropping anything below it here removes only pairs that could never have passed.
    // Measured: 175 of 7,069,934 accepted multimer pairs lost, 0.002%, for 10% of the time.
    par.covThr = (par.covMode == Parameters::COV_MODE_BIDIRECTIONAL) ? par.filtChainTmThr : 0.0;
    cmd.addVariable("MULTIMERALIGN_PAR", par.createParameterString(par.structurealign).c_str());
    par.evalThr = origEvalThr;
    par.addBacktrace = origAddBacktrace;
    par.PARAM_ADD_BACKTRACE.wasSet = origAddBacktraceWasSet;
    par.covThr = origAlnCovThr;
    cmd.addVariable("SCOREMULTIMER_PAR", par.createParameterString(par.scoremultimer).c_str());
    cmd.addVariable("CLUSTER_PAR", par.createParameterString(par.clust).c_str());
    cmd.addVariable("REMOVE_TMP", par.removeTmpFiles ? "TRUE" : NULL);
    cmd.addVariable("VERBOSITY_PAR", par.createParameterString(par.onlyverbosity).c_str());

    std::string program = tmpDir + "/linmultimercluster.sh";
    FileUtil::writeFile(program, linmultimercluster_sh, linmultimercluster_sh_len);
    cmd.execProgram(program.c_str(), par.filenames);

    // Should never get here
    assert(false);
    return EXIT_FAILURE;
}
