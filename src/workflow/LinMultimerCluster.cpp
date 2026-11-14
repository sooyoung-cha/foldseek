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
    // The chain prefilter runs twice, cheaply over everything and then at full -s over
    // what the first pass could not cluster, so the first pass does not have to be right
    // on its own. Measured on 51,770 family sampled dimers: one -s 6.0 prefilter costs
    // 17 s and leaves 2.22% of the representatives with a TM 0.7 neighbour still out
    // there; -s 4.0 followed by the pass over the representatives costs 5 s + 4 s and
    // leaves 2.25%. Half the candidate generation, same clustering.
    p->round0Sensitivity = 4.0;
    // The first pass is the quadratic stage, so it is also run over as few chains as
    // possible: linclust groups the near identical ones and the prefilter sees one chain
    // per group. Measured on the same set, 103,540 chains collapse to 67,097 groups and
    // the first pass goes 5 s to 2 s at the same 1.1% rep-rep. The collapse is a property
    // of the database, not of its size: the whole 20,179,006 chain teddymer database
    // collapses to 66.5%, the same ratio, and the stage it makes cheaper is the one that
    // grows quadratically, so it is worth more the bigger the input is.
    p->groupPrefilter = 1;
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
    // Two prefilter strings, because the chain prefilter runs twice: PREFILTER0_PAR over
    // the whole database at --round0-sensitivity, PREFILTER_PAR over the multimers it
    // left unclustered, at the -s the caller asked for. --round0-sensitivity 0 puts the
    // first pass back at -s as well, which is the way to pay for the cascade's second
    // pass without taking the discount on the first.
    float origSens = par.sensitivity;
    if (par.round0Sensitivity > 0.0f) {
        par.sensitivity = par.round0Sensitivity;
    }
    cmd.addVariable("PREFILTER0_PAR", par.createParameterString(par.prefilter).c_str());
    par.sensitivity = origSens;
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
    // The cheap first round leaves a multimer that found no true neighbour joined to
    // nothing, which makes it its own representative. The representatives are therefore
    // exactly the input of a second, deeper pass, and there are far fewer of them, so the
    // quadratic prefilter costs (reps/all)^2 of what it would have cost on the whole set.
    // Only the exhaustive prefilter skips it: that one already compared every pair.
    int origSubDbMode = par.subDbMode;
    par.subDbMode = Parameters::SUBDB_MODE_SOFT;
    cmd.addVariable("CREATESUBDB_PAR", par.createParameterString(par.createsubdb).c_str());
    par.subDbMode = origSubDbMode;
    // Near identical chains have near identical hit lists, so the quadratic stage is
    // computing the same list many times. linclust groups them and the chain prefilter
    // runs over one chain per group; multimerprefilter and multimerneighbors read the
    // representative's record for every member and read a hit on a representative as a
    // hit on its whole group, so the expanded list is never materialised. Measured on
    // 51,770 family sampled dimers: 103,540 chains collapse to 67,097 groups, the chain
    // prefilter goes 6 s to 2 s and rep-rep at TM 0.7 goes 1.13% to 1.07%. The grouping
    // only picks which pairs are looked at; the alignment and scoremultimer still run on
    // the real chains, so a group that is not uniform costs a rejected candidate, never
    // a wrong merge.
    //
    // A GPU run is the exception: there the prefilter reads a padded database that
    // createsubdb cannot carve the group representatives out of, so the grouping steps
    // aside and the prefilter sees every chain, as it did before.
    if (par.gpu == 1 && par.groupPrefilter == 1) {
        if (par.PARAM_GROUP_PREFILTER.wasSet) {
            Debug(Debug::WARNING) << "--group-prefilter does not work on a GPU database, switching it off\n";
        }
        par.groupPrefilter = 0;
    }
    cmd.addVariable("GROUP_PREFILTER", par.groupPrefilter == 1 ? "TRUE" : NULL);
    if (par.groupPrefilter == 1) {
        // A group has to mean "the same hit list", because its members share one.
        // That is a question about the 3Di k-mers the prefilter probes with, not about
        // structural similarity, which is why the threshold here is linclust's own and
        // is not tied to --chain-tm-threshold. 0.7 is what was measured: 103,540 chains
        // collapse to 67,097 groups on the family sampled set, the same collapse the
        // whole 20,179,006 chain database gives (66.5%), so the speedup carries over.
        //
        // Everything linclust reads has to be pinned back to linclust's own defaults,
        // because createParameterString serialises whatever the workflow has already
        // tuned for its other stages. Three of them were measured leaking in here:
        // --set-mode 1 segfaults linclust's clustering, -e 10 and --spaced-kmer-mode 1
        // each group far more loosely than the 0.7 that was measured (103,540 chains to
        // 63,210 groups instead of 67,097, and rep-rep at TM 0.7 from 1.07% to 1.44%).
        float origSeqId = par.seqIdThr, origCov = par.covThr;
        int origCovMode = par.covMode, origSetMode = par.clusteringSetMode;
        int origSpaced = par.spacedKmer;
        double origE = par.evalThr;
        par.clusteringSetMode = 0;
        par.evalThr = 0.001;
        par.spacedKmer = 0;
        par.seqIdThr = par.groupMinSeqId;
        par.covThr = 0.8;
        par.covMode = Parameters::COV_MODE_BIDIRECTIONAL;
        par.PARAM_MIN_SEQ_ID.wasSet = true;
        par.PARAM_C.wasSet = true;
        par.PARAM_COV_MODE.wasSet = true;
        par.PARAM_E.wasSet = true;
        cmd.addVariable("LINCLUST_PAR", par.createParameterString(par.linclustworkflow).c_str());
        par.seqIdThr = origSeqId; par.covThr = origCov; par.covMode = origCovMode;
        par.clusteringSetMode = origSetMode; par.evalThr = origE;
        par.spacedKmer = origSpaced;
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
    cmd.addVariable("MULTIMEREXPAND_PAR0", par.createParameterString(par.multimerexpand).c_str());
    par.multimerExpandMode = 1;
    cmd.addVariable("MULTIMEREXPAND_PAR1", par.createParameterString(par.multimerexpand).c_str());
    par.multimerExpandMode = origMode;

    // The second pass over the representatives is the safety net of the cascade, so it
    // is run exactly, on the real chain keys, and the stages that read it are given the
    // ungrouped parameter strings. Only round 0 goes through the groups.
    par.multimerExpandMode = 0;
    cmd.addVariable("MULTIMERNEIGHBORS_PAR2", par.createParameterString(par.multimerneighbors).c_str());
    par.multimerExpandMode = origMode;
    cmd.addVariable("MULTIMERPREFILTER_PAR2", par.createParameterString(par.multimerprefilter).c_str());

    // The two stages that read round 0's chain prefilter have to know the groups;
    // everything downstream of them sees ordinary chain keys again.
    if (par.groupPrefilter == 1) {
        par.chainGroups = tmpDir + "/cgroups";
        par.PARAM_CHAIN_GROUPS.wasSet = true;
    }
    par.multimerExpandMode = 0;
    cmd.addVariable("MULTIMERNEIGHBORS_PAR", par.createParameterString(par.multimerneighbors).c_str());
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
