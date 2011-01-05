/*
  Copyright (C) 2007, Lou Hafer, International Business Machines Corporation
  and others.  All Rights Reserved.

  This code is licensed under the terms of the Eclipse Public License (EPL).

  $Id$
*/
/*
  This file is part of cbc-generic.
*/

#include "CbcConfig.h"
#include "CoinPragma.hpp"

#include <cassert>

#include "CbcGenCtlBlk.hpp"

namespace {

char svnid[] = "$Id: CbcGenCtlBlk.cpp 1173 2009-06-04 09:44:10Z forrest $" ;

}

/*
  Constructor for cbc-generic control block.

  Set up defaults for the cbc-generic control block. Note that prototypes for
  cut generators and heuristics will be created on demand; see the access
  functions.

  Once this structure settles down, simple intialisation should move up to
  the standard `:' block. In the meantime, this avoids complaints about
  ordering.
*/

CbcGenCtlBlk::CbcGenCtlBlk ()

{
    version_ = CBC_GENERIC_VERSION ;
    /*
      It's unclear to me that this is a good choice for dfltDirectory. Makes
      sense for commands, but seems unnecessary for data files. Perhaps a null
      string instead?
    */
    char dirsep = CoinFindDirSeparator() ;
    dfltDirectory_ = (dirsep == '/' ? "./" : ".\\") ;
    lastMpsIn_ = "" ;
    allowImportErrors_ = false ;
    lastSolnOut_ = "stdout" ;
    printMode_ = 0 ;
    printMask_ = "" ;

    paramVec_ = 0 ;
    genParams_.first_ = 0 ;
    genParams_.last_ = 0 ;
    cbcParams_.first_ = 0 ;
    cbcParams_.last_ = 0 ;
    osiParams_.first_ = 0 ;
    osiParams_.last_ = 0 ;

    verbose_ = 0 ;
    paramsProcessed_ = 0 ;
    defaultSettings_ = true ;

    debugCreate_ = "" ;
    debugFile_ = "" ;
    debugSol_.numCols_ = -1 ;
    debugSol_.values_ = 0 ;

    printOpt_ = 0 ;

    /*
      Assigning us_en to cur_lang_ is entirely bogus, but CoinMessages::Language
      does not provide an `unspecified' code.
    */
    msgHandler_ = new CoinMessageHandler() ;
    ourMsgHandler_ = true ;
    cur_lang_ = CoinMessages::us_en ;
    msgs_ = 0 ;
    logLvl_ = 0 ;

    totalTime_ = 0.0 ;

    model_ = 0 ;
    dfltSolver_ = 0 ;
    goodModel_ = false ;
    bab_.majorStatus_ = BACInvalid ;
    bab_.minorStatus_ = BACmInvalid ;
    bab_.where_ = BACwInvalid ;
    bab_.haveAnswer_ = false ;
    bab_.answerSolver_ = 0 ;

    preProcess_ = CbcGenCtlBlk::IPPSOS ;
    cutDepth_ = -1 ;

    probing_.action_ = CbcGenCtlBlk::CGIfMove ;
    probing_.proto_ = 0 ;
    probing_.usingObjective_ = true ;
    probing_.maxPass_ = 3 ;
    probing_.maxPassRoot_ = 3 ;
    probing_.maxProbe_ = 10 ;
    probing_.maxProbeRoot_ = 50 ;
    probing_.maxLook_ = 10 ;
    probing_.maxLookRoot_ = 50 ;
    probing_.maxElements_ = 200 ;
    probing_.rowCuts_ = 3 ;

    clique_.action_ = CbcGenCtlBlk::CGIfMove ;
    clique_.proto_ = 0 ;
    clique_.starCliqueReport_ = false ;
    clique_.rowCliqueReport_ = false ;
    clique_.minViolation_ = 0.1 ;

    flow_.action_ = CbcGenCtlBlk::CGIfMove ;
    flow_.proto_ = 0 ;

    gomory_.action_ = CbcGenCtlBlk::CGIfMove ;
    gomory_.proto_ = 0 ;
    gomory_.limit_ = 50 ;
    gomory_.limitAtRoot_ = 512 ;

    knapsack_.action_ = CbcGenCtlBlk::CGIfMove ;
    knapsack_.proto_ = 0 ;

    // landp_action_ = CbcGenCtlBlk::CGOff ;
    // landp_.proto_ = 0 ;

    mir_.action_ = CbcGenCtlBlk::CGIfMove ;
    mir_.proto_ = 0 ;

    oddHole_.action_ = CbcGenCtlBlk::CGOff ;
    oddHole_.proto_ = 0 ;

    redSplit_.action_ = CbcGenCtlBlk::CGRoot ;
    redSplit_.proto_ = 0 ;

    twomir_.action_ = CbcGenCtlBlk::CGRoot ;
    twomir_.proto_ = 0 ;
    twomir_.maxElements_ = 250 ;

    fpump_.action_ = CbcGenCtlBlk::CGOn ;
    fpump_.proto_ = 0 ;

    combine_.action_ = CbcGenCtlBlk::CGOn ;
    combine_.proto_ = 0 ;
    combine_.trySwap_ = 1 ;

    greedyCover_.action_ = CbcGenCtlBlk::CGOn ;
    greedyCover_.proto_ = 0 ;
    greedyEquality_.action_ = CbcGenCtlBlk::CGOn ;
    greedyEquality_.proto_ = 0 ;

    localTree_.action_ = CbcGenCtlBlk::CGOff ;
    localTree_.proto_ = 0 ;
    localTree_.soln_ = 0 ;
    localTree_.range_ = 10 ;
    localTree_.typeCuts_ = 0 ;
    localTree_.maxDiverge_ = 0 ;
    localTree_.timeLimit_ = 10000 ;
    localTree_.nodeLimit_ = 2000 ;
    localTree_.refine_ = true ;

    rounding_.action_ = CbcGenCtlBlk::CGOn ;
    rounding_.proto_ = 0 ;

    djFix_.action_ = false ;
    djFix_.threshold_ = 1.0e100 ;

    priorityAction_ = CbcGenCtlBlk::BPOff ;
    /*
      The value for numBeforeTrust is as recommended by Achterberg. Cbc's
      implementation doesn't really have a parameter equivalent to Achterberg's
      dynamic limit on number of strong branching evaluations, so go with a fairly
      large default. As of 06.12.16, the magic number for shadow price mode meant
      `use shadow prices (penalties, I think) if there's no strong branching info'.
    */
    chooseStrong_.numBeforeTrust_ = 8 ;
    chooseStrong_.numStrong_ = 100 ;
    chooseStrong_.shadowPriceMode_ = 1 ;

    return ;
}

/*
  Note that we don't want to delete dfltSolver_ here because it's just a
  copy of the pointer held in the solvers map over in CbcGenSolvers.cpp.
*/
CbcGenCtlBlk::~CbcGenCtlBlk ()

{
    if (model_) delete model_ ;
    if (bab_.answerSolver_) delete bab_.answerSolver_ ;

    if (probing_.proto_) delete probing_.proto_ ;
    if (clique_.proto_) delete clique_.proto_ ;
    if (flow_.proto_) delete flow_.proto_ ;
    if (gomory_.proto_) delete gomory_.proto_ ;
    if (knapsack_.proto_) delete knapsack_.proto_ ;
    if (mir_.proto_) delete mir_.proto_ ;
    if (oddHole_.proto_) delete oddHole_.proto_ ;
    if (redSplit_.proto_) delete redSplit_.proto_ ;
    if (twomir_.proto_) delete twomir_.proto_ ;

    if (fpump_.proto_) delete fpump_.proto_ ;
    if (combine_.proto_) delete combine_.proto_ ;
    if (greedyCover_.proto_) delete greedyCover_.proto_ ;
    if (greedyEquality_.proto_) delete greedyEquality_.proto_ ;
    if (rounding_.proto_) delete rounding_.proto_ ;

    if (msgHandler_ && ourMsgHandler_) delete msgHandler_ ;
    if (msgs_) delete msgs_ ;

    return ;
}

/*
  Access functions for cut generators and heuristics. These support lazy
  creation --- if action_ is other than CGOff, an object is created if
  necessary and a pointer is stored in proto_. The pointer is returned as
  a generic CglCutGenerator or CbcHeuristic. The return value of the function
  is the value of action_.

  Because the model may have changed, the default for heuristics is to delete
  any existing object and create a new one. This can be suppressed if desired.
*/

CbcGenCtlBlk::CGControl CbcGenCtlBlk::getProbing (CglCutGenerator *&gen)

{
    if (probing_.action_ != CbcGenCtlBlk::CGOff && probing_.proto_ == 0) {
        probing_.proto_ = new CglProbing() ;
        probing_.proto_->setUsingObjective(probing_.usingObjective_) ;
        probing_.proto_->setMaxPass(probing_.maxPass_) ;
        probing_.proto_->setMaxPassRoot(probing_.maxPassRoot_) ;
        probing_.proto_->setMaxProbe(probing_.maxProbe_) ;
        probing_.proto_->setMaxProbeRoot(probing_.maxProbeRoot_) ;
        probing_.proto_->setMaxLook(probing_.maxLook_) ;
        probing_.proto_->setMaxLookRoot(probing_.maxLookRoot_) ;
        probing_.proto_->setMaxElements(probing_.maxElements_) ;
        probing_.proto_->setRowCuts(probing_.rowCuts_) ;
    }
    gen = dynamic_cast<CglCutGenerator *>(probing_.proto_) ;

    return (probing_.action_) ;
}


CbcGenCtlBlk::CGControl CbcGenCtlBlk::getClique (CglCutGenerator *&gen)

{
    if (clique_.action_ != CbcGenCtlBlk::CGOff && clique_.proto_ == 0) {
        clique_.proto_ = new CglClique() ;
        clique_.proto_->setStarCliqueReport(clique_.starCliqueReport_) ;
        clique_.proto_->setRowCliqueReport(clique_.rowCliqueReport_) ;
        clique_.proto_->setMinViolation(clique_.minViolation_) ;
    }
    gen = dynamic_cast<CglCutGenerator *>(clique_.proto_) ;

    return (clique_.action_) ;
}


CbcGenCtlBlk::CGControl CbcGenCtlBlk::getFlow (CglCutGenerator *&gen)

{
    if (flow_.action_ != CbcGenCtlBlk::CGOff && flow_.proto_ == 0) {
        flow_.proto_ = new CglFlowCover() ;
    }
    gen = dynamic_cast<CglCutGenerator *>(flow_.proto_) ;

    return (flow_.action_) ;
}


CbcGenCtlBlk::CGControl CbcGenCtlBlk::getGomory (CglCutGenerator *&gen)

{
    if (gomory_.action_ != CbcGenCtlBlk::CGOff && gomory_.proto_ == 0) {
        gomory_.proto_ = new CglGomory() ;
        gomory_.proto_->setLimitAtRoot(gomory_.limitAtRoot_) ;
        gomory_.proto_->setLimit(gomory_.limit_) ;
    }
    gen = dynamic_cast<CglCutGenerator *>(gomory_.proto_) ;

    return (gomory_.action_) ;
}


CbcGenCtlBlk::CGControl CbcGenCtlBlk::getKnapsack (CglCutGenerator *&gen)

{
    if (knapsack_.action_ != CbcGenCtlBlk::CGOff && knapsack_.proto_ == 0) {
        knapsack_.proto_ = new CglKnapsackCover() ;
    }
    gen = dynamic_cast<CglCutGenerator *>(knapsack_.proto_) ;

    return (knapsack_.action_) ;
}


CbcGenCtlBlk::CGControl CbcGenCtlBlk::getMir (CglCutGenerator *&gen)

{
    if (mir_.action_ != CbcGenCtlBlk::CGOff && mir_.proto_ == 0) {
        mir_.proto_ = new CglMixedIntegerRounding2() ;
    }
    gen = dynamic_cast<CglCutGenerator *>(mir_.proto_) ;

    return (mir_.action_) ;
}


CbcGenCtlBlk::CGControl CbcGenCtlBlk::getRedSplit (CglCutGenerator *&gen)

{
    if (redSplit_.action_ != CbcGenCtlBlk::CGOff && redSplit_.proto_ == 0) {
        redSplit_.proto_ = new CglRedSplit() ;
    }
    gen = dynamic_cast<CglCutGenerator *>(redSplit_.proto_) ;

    return (redSplit_.action_) ;
}


CbcGenCtlBlk::CGControl CbcGenCtlBlk::getTwomir (CglCutGenerator *&gen)

{
    if (twomir_.action_ != CbcGenCtlBlk::CGOff && twomir_.proto_ == 0) {
        twomir_.proto_ = new CglTwomir() ;
        twomir_.proto_->setMaxElements(twomir_.maxElements_) ;
    }
    gen = dynamic_cast<CglCutGenerator *>(twomir_.proto_) ;

    return (twomir_.action_) ;
}


CbcGenCtlBlk::CGControl
CbcGenCtlBlk::getFPump (CbcHeuristic *&gen, CbcModel *model,
                        bool alwaysCreate)

{
    if (fpump_.action_ != CbcGenCtlBlk::CGOff &&
            (fpump_.proto_ == 0 || alwaysCreate)) {
        if (fpump_.proto_) {
            delete fpump_.proto_ ;
        }
        fpump_.proto_ = new CbcHeuristicFPump(*model) ;
        fpump_.proto_->setMaximumPasses(fpump_.iters_) ;
    }
    gen = dynamic_cast<CbcHeuristic *>(fpump_.proto_) ;

    return (fpump_.action_) ;
}


CbcGenCtlBlk::CGControl
CbcGenCtlBlk::getCombine (CbcHeuristic *&gen, CbcModel *model,
                          bool alwaysCreate)

{
    if (combine_.action_ != CbcGenCtlBlk::CGOff &&
            (combine_.proto_ == 0 || alwaysCreate)) {
        if (combine_.proto_) {
            delete combine_.proto_ ;
        }
        combine_.proto_ = new CbcHeuristicLocal(*model) ;
        combine_.proto_->setSearchType(combine_.trySwap_) ;
    }
    gen = dynamic_cast<CbcHeuristic *>(combine_.proto_) ;

    return (combine_.action_) ;
}


CbcGenCtlBlk::CGControl
CbcGenCtlBlk::getGreedyCover (CbcHeuristic *&gen, CbcModel *model,
                              bool alwaysCreate)

{
    if (greedyCover_.action_ != CbcGenCtlBlk::CGOff &&
            (greedyCover_.proto_ == 0 || alwaysCreate)) {
        if (greedyCover_.proto_) {
            delete greedyCover_.proto_ ;
        }
        greedyCover_.proto_ = new CbcHeuristicGreedyCover(*model) ;
    }
    gen = dynamic_cast<CbcHeuristic *>(greedyCover_.proto_) ;

    return (greedyCover_.action_) ;
}


CbcGenCtlBlk::CGControl
CbcGenCtlBlk::getGreedyEquality (CbcHeuristic *&gen, CbcModel *model,
                                 bool alwaysCreate)

{
    if (greedyEquality_.action_ != CbcGenCtlBlk::CGOff &&
            (greedyEquality_.proto_ == 0 || alwaysCreate)) {
        if (greedyEquality_.proto_) {
            delete greedyEquality_.proto_ ;
        }
        greedyEquality_.proto_ = new CbcHeuristicGreedyEquality(*model) ;
    }
    gen = dynamic_cast<CbcHeuristic *>(greedyEquality_.proto_) ;

    return (greedyEquality_.action_) ;
}


CbcGenCtlBlk::CGControl
CbcGenCtlBlk::getRounding (CbcHeuristic *&gen, CbcModel *model,
                           bool alwaysCreate)

{
    if (rounding_.action_ != CbcGenCtlBlk::CGOff &&
            (rounding_.proto_ == 0 || alwaysCreate)) {
        if (rounding_.proto_) {
            delete rounding_.proto_ ;
        }
        rounding_.proto_ = new CbcRounding(*model) ;
    }
    gen = dynamic_cast<CbcHeuristic *>(rounding_.proto_) ;

    return (rounding_.action_) ;
}


CbcGenCtlBlk::CGControl
CbcGenCtlBlk::getTreeLocal (CbcTreeLocal *&localTree, CbcModel *model,
                            bool alwaysCreate)

{
    if (localTree_.action_ != CbcGenCtlBlk::CGOff &&
            (localTree_.proto_ == 0 || alwaysCreate)) {
        if (localTree_.proto_) {
            delete localTree_.proto_ ;
        }
        localTree_.proto_ =
            new CbcTreeLocal(model, localTree_.soln_, localTree_.range_,
                             localTree_.typeCuts_, localTree_.maxDiverge_,
                             localTree_.timeLimit_, localTree_.nodeLimit_,
                             localTree_.refine_) ;
    }
    localTree = localTree_.proto_ ;

    return (localTree_.action_) ;
}

/*
  A bunch of little translation helper routines leading up to a version of
  setBaBStatus that figures it all out given a CbcModel and BACWhere code.
  This translation needs to be centralised to avoid sprinkling magic numbers
  all through the code.

  Be a bit careful with the translation routines --- they aren't sensitive to
  where the search stopped.
*/

CbcGenCtlBlk::BACMajor CbcGenCtlBlk::translateMajor (int status)

{
    switch (status) {
    case -1: {
        return (BACNotRun) ;
    }
    case  0: {
        return (BACFinish) ;
    }
    case  1: {
        return (BACStop) ;
    }
    case  2: {
        return (BACAbandon) ;
    }
    case  5: {
        return (BACUser) ;
    }
    default: {
        return (BACInvalid) ;
    }
    }
}


CbcGenCtlBlk::BACMinor CbcGenCtlBlk::translateMinor (int status)

{
    switch (status) {
    case -1: {
        return (BACmInvalid) ;
    }
    case  0: {
        return (BACmFinish) ;
    }
    case  1: {
        return (BACmInfeas) ;
    }
    case  2: {
        return (BACmGap) ;
    }
    case  3: {
        return (BACmNodeLimit) ;
    }
    case  4: {
        return (BACmTimeLimit) ;
    }
    case  5: {
        return (BACmUser) ;
    }
    case  6: {
        return (BACmSolnLimit) ;
    }
    case  7: {
        return (BACmUbnd) ;
    }
    default: {
        return (BACmOther) ;
    }
    }
}

/*
  A bit different --- given an OSI, use its interrogation functions to choose
  an appropriate BACMinor code. Not everything matches up, eh?
*/
CbcGenCtlBlk::BACMinor
CbcGenCtlBlk::translateMinor (const OsiSolverInterface *osi)

{
    if (osi->isProvenOptimal()) {
        return (BACmFinish) ;
    } else if (osi->isProvenPrimalInfeasible()) {
        return (BACmInfeas) ;
    } else if (osi->isProvenDualInfeasible()) {
        return (BACmUbnd) ;
    } else {
        return (BACmOther) ;
    }
}


/*
  A routine to set the bab_ status block given a CbcModel and an indication
  of where we're at in the search. Really, this is just a big mapping from
  CbcModel codes to CbcGeneric codes.
*/

void CbcGenCtlBlk::setBaBStatus (const CbcModel *model, BACWhere where,
                                 bool haveAnswer,
                                 OsiSolverInterface *answerSolver)

{
    CbcGenCtlBlk::BACMajor major ;
    CbcGenCtlBlk::BACMinor minor ;

    major = translateMajor(model->status()) ;

    if (where == CbcGenCtlBlk::BACwBareRoot ||
            where == CbcGenCtlBlk::BACwIPPRelax) {
        minor = translateMinor(model->solver()) ;
    } else {
        minor = translateMinor(model->secondaryStatus()) ;
    }

    setBaBStatus(major, minor, where, haveAnswer, answerSolver) ;

    return ;
}


/*
  Last, but not least, a routine to print the result.
*/

void CbcGenCtlBlk::printBaBStatus ()

{
    std::cout << "BAC result: stopped " ;

    switch (bab_.where_) {
    case BACwNotStarted: {
        std::cout << "before root relaxation" ;
        break ;
    }
    case BACwBareRoot: {
        std::cout << "after root relaxation" ;
        break ;
    }
    case BACwIPP: {
        std::cout << "after integer preprocessing" ;
        break ;
    }
    case BACwIPPRelax: {
        std::cout << "after solving preprocessed relaxation" ;
        break ;
    }
    case BACwBAC: {
        std::cout << "after branch-and-cut" ;
        break ;
    }
    default: {
        std::cout << "!!invalid phase code!!" ;
        break ;
    }
    }

    std::cout << std::endl << "    Branch-and-cut " ;

    switch (bab_.majorStatus_) {
    case BACNotRun: {
        std::cout << "never got started" ;
        break ;
    }
    case BACFinish: {
        std::cout << "finished" ;
        break ;
    }
    case BACStop: {
        std::cout << "stopped on a limit" ;
        break ;
    }
    case BACAbandon: {
        std::cout << "was abandoned" ;
        break ;
    }
    case BACUser: {
        std::cout << "stopped due to a user event" ;
        break ;
    }
    default: {
        std::cout << "!!invalid major status code!!" ;
        break ;
    }
    }

    std::cout << "; minor status is " ;

    switch (bab_.minorStatus_) {
    case BACmFinish: {
        std::cout << "optimal" ;
        break ;
    }
    case BACmInfeas: {
        std::cout << "infeasible" ;
        break ;
    }
    case BACmUbnd: {
        std::cout << "unbounded" ;
        break ;
    }
    case BACmGap: {
        std::cout << "reached specified integrality gap." ;
        break ;
    }
    case BACmNodeLimit: {
        std::cout << "reached node limit" ;
        break ;
    }
    case BACmTimeLimit: {
        std::cout << "reached time limit" ;
        break ;
    }
    case BACmSolnLimit: {
        std::cout << "reached limit on number of solutions" ;
        break ;
    }
    case BACmUser: {
        std::cout << "stopped due to a user event" ;
        break ;
    }
    case BACmOther: {
        std::cout << "other" ;
        break ;
    }
    default: {
        std::cout << "!!invalid minor status code!!" ;
        break ;
    }
    }

    std::cout << "." << std::endl ;
}

