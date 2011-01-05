/*
  Copyright (C) 2007, Lou Hafer, International Business Machines Corporation
  and others.  All Rights Reserved.

  This code is licensed under the terms of the Eclipse Public License (EPL).

  $Id$
*/
/*
  This file is part of cbc-generic.
*/

#include <iostream>

#include "CoinTime.hpp"

#include "OsiSolverInterface.hpp"
#include "OsiChooseVariable.hpp"

#include "CglPreProcess.hpp"

#include "CbcModel.hpp"
#include "CbcCutGenerator.hpp"
#include "CbcBranchActual.hpp"
#include "CbcStrategy.hpp"

#include "CbcGenCtlBlk.hpp"
#include "CbcGenParam.hpp"
#include "CbcGenCbcParam.hpp"

#define CBC_TRACK_SOLVERS 1
// #define COIN_CBC_VERBOSITY 5

/*
  The support functions for the main branch-and-cut action routine.
*/

namespace {

char svnid[] = "$Id: CbcGenBaB.cpp 1173 2009-06-04 09:44:10Z forrest $" ;

/*
  A hack to fix variables based on reduced cost prior to branch-and-cut. Note
  that we're *not* looking at the integrality gap here. Given the reduced costs
  of the root relaxation, we're simply placing a bet that variables with really
  unfavourable reduced costs that are at their most favourable bound in the
  root relaxation will never move from that bound.

  For the standard OsiSolverInterface, this requires a bit of effort as the
  solution and bounds arrays are const and the functions to change them have
  incompatible interfaces.
*/

void reducedCostHack (OsiSolverInterface *osi, double threshold)

{
    int numCols = osi->getNumCols() ;
    int i ;
    const double *lower = osi->getColLower() ;
    const double *upper = osi->getColUpper() ;
    const double *solution = osi->getColSolution() ;
    const double *dj = osi->getReducedCost() ;
    /*
      First task: scan the columns looking for variables that are at their
      favourable bound and have reduced cost that exceeds the threshold. Remember
      the column index and the value.
    */
    double *chgBnds = new double [numCols] ;
    int *chgCols = new int [numCols] ;

    int numFixed = 0 ;
    for (i = 0 ; i < numCols ; i++) {
        if (osi->isInteger(i)) {
            double value = solution[i] ;
            if (value < lower[i] + 1.0e-5 && dj[i] > threshold) {
                chgCols[numFixed] = i ;
                chgBnds[numFixed] = lower[i] ;
                numFixed++ ;
            } else if (value > upper[i] - 1.0e-5 && dj[i] < -threshold) {
                chgCols[numFixed] = i ;
                chgBnds[numFixed] = upper[i] ;
                numFixed++ ;
            }
        }
    }
    /*
      Second task: For variables that we want to fix, we need to:
        * Prepare an array with the new lower and upper bounds for variables that
          will be fixed. setColSetBounds requires an array with column indices and
          an array with new values for both bounds.
        * Set the correct value in a copy of the current solution. setColSolution
          requires a complete solution.
    */
    if (numFixed > 0) {
        double *newSoln = CoinCopyOfArray(solution, numCols) ;
        double *newBnds = new double [2*numFixed] ;
        double *bndPtr = &newBnds[0] ;
        for (i = 0 ; i < numFixed ; i++) {
            int j = chgCols[i] ;
            double val = chgBnds[i] ;
            *bndPtr++ = val ;
            *bndPtr++ = val ;
            newSoln[j] = val ;
        }
        osi->setColSetBounds(&chgCols[0], &chgCols[numFixed], &newBnds[0]) ;
        osi->setColSolution(&newSoln[0]) ;

        std::cout
            << "Reduced cost fixing prior to B&C: " << numFixed
            << " columns fixed." << std::endl ;

        delete[] newSoln ;
        delete[] newBnds ;
    }

    delete[] chgBnds ;
    delete[] chgCols ;

    return ;
}

/*
  Helper routine to solve a continuous relaxation and print something
  intelligent when the result is other than optimal. Returns true if the
  result is optimal, false otherwise.
*/

bool solveRelaxation (CbcModel *model)

{
    OsiSolverInterface *osi = model->solver() ;

    model->initialSolve() ;

    if (!(osi->isProvenOptimal())) {
        bool reason = false ;
        if (osi->isProvenPrimalInfeasible()) {
            std::cout
                << "Continuous relaxation is primal infeasible." << std::endl ;
            reason = true ;
        }
        if (osi->isProvenDualInfeasible()) {
            std::cout
                << "Continuous relaxation is dual infeasible." << std::endl ;
            reason = true ;
        }
        if (osi->isIterationLimitReached()) {
            std::cout
                << "Continuous solver reached iteration limit." << std::endl ;
            reason = true ;
        }
        if (osi->isAbandoned()) {
            std::cout
                << "Continuous solver abandoned the problem." << std::endl ;
            reason = true ;
        }
        if (reason == false) {
            std::cout
                << "Continuous solver failed for unknown reason." << std::endl ;
        }
        return (false) ;
    }

    return (true) ;
}


/*
  Helper routine to establish a priority vector.
*/

void setupPriorities (CbcModel *model, CbcGenCtlBlk::BPControl how)

{
    int numCols = model->getNumCols() ;
    int *sort = new int[numCols] ;
    double *dsort = new double[numCols] ;
    int *priority = new int[numCols] ;
    const double *objective = model->getObjCoefficients() ;
    int iColumn ;
    int n = 0 ;
    bool priorityOK = true ;

    for (iColumn = 0 ; iColumn < numCols ; iColumn++) {
        if (model->isInteger(iColumn)) {
            sort[n] = n ;
            if (how == CbcGenCtlBlk::BPCost) {
                dsort[n++] = -objective[iColumn] ;
            } else if (how == CbcGenCtlBlk::BPOrder) {
                dsort[n++] = iColumn ;
            } else {
                std::cerr
                    << "setupPriorities: Unrecognised priority specification."
                    << std::endl ;
                priorityOK = false ;
            }
        }
    }

    if (priorityOK) {
        CoinSort_2(dsort, dsort + n, sort) ;

        int level = 0 ;
        double last = -1.0e100 ;
        for (int i = 0 ; i < n ; i++) {
            int iPut = sort[i] ;
            if (dsort[i] != last) {
                level++ ;
                last = dsort[i] ;
            }
            priority[iPut] = level ;
        }

        model->passInPriorities(priority, false) ;
    }

    delete [] priority ;
    delete [] sort ;
    delete [] dsort ;

    return ;
}


/*
  Helper routine to install a batch of heuristics. Each call to getXXXHeuristic
  will return a pointer to the heuristic object in gen iff the heuristic is
  enabled.
*/

void installHeuristics (CbcGenCtlBlk *ctlBlk, CbcModel *model)

{
    CbcGenCtlBlk::CGControl action ;
    CbcHeuristic *gen ;
    CbcTreeLocal *localTree ;
    /*
      FPump goes first because it only works before there's a solution.
    */
    action = ctlBlk->getFPump(gen, model) ;
    if (action != CbcGenCtlBlk::CGOff) {
        model->addHeuristic(gen, "FPump") ;
    }
    action = ctlBlk->getRounding(gen, model) ;
    if (action != CbcGenCtlBlk::CGOff) {
        model->addHeuristic(gen, "Rounding") ;
    }
    action = ctlBlk->getCombine(gen, model) ;
    if (action != CbcGenCtlBlk::CGOff) {
        model->addHeuristic(gen, "Combine") ;
    }
    action = ctlBlk->getGreedyCover(gen, model) ;
    if (action != CbcGenCtlBlk::CGOff) {
        model->addHeuristic(gen, "GCov") ;
    }
    action = ctlBlk->getGreedyEquality(gen, model) ;
    if (action != CbcGenCtlBlk::CGOff) {
        model->addHeuristic(gen, "GEq") ;
    }
    /*
      This one's a bit different. We acquire the local tree and install it in the
      model.
    */
    action = ctlBlk->getTreeLocal(localTree, model) ;
    if (action != CbcGenCtlBlk::CGOff) {
        model->passInTreeHandler(*localTree) ;
    }

    return ;
}


/*
  Helper routine to install cut generators.

  I need to install the new lift-and-project generator (LandP). Also need to
  figure out stored cuts.
*/

void installCutGenerators (CbcGenCtlBlk *ctlBlk, CbcModel *model)

{
    int switches[20] ;
    int genCnt = 0 ;
    CbcGenCtlBlk::CGControl action ;
    CglCutGenerator *gen ;

    /*
      The magic numbers for the howOften parameter that determines how often the
      generator is invoked. -100 is disabled, -99 is root only, -98 will stay
      active only so long as it generates cuts that improve the objective. A value
      1 <= k <= 90 means the generator will be called every k nodes. If k is
      negative, then it can be switched off if unproductive. If k is positive,
      it'll carry on regardless.
    */
    int howOften[CbcGenCtlBlk::CGMarker] ;
    howOften[CbcGenCtlBlk::CGOff] = -100 ;
    howOften[CbcGenCtlBlk::CGOn] = -1 ;
    howOften[CbcGenCtlBlk::CGRoot] = -99 ;
    howOften[CbcGenCtlBlk::CGIfMove] = -98 ;
    howOften[CbcGenCtlBlk::CGForceOn] = 1 ;
    howOften[CbcGenCtlBlk::CGForceBut] = 1 ;

    /*
      A negative value for rowCuts means that the specified actions happen only at
      the root.
    */
    action = ctlBlk->getProbing(gen) ;
    if (action != CbcGenCtlBlk::CGOff) {
        if (action == CbcGenCtlBlk::CGForceBut) {
            CglProbing *probingGen = dynamic_cast<CglProbing *>(gen) ;
            probingGen->setRowCuts(-3) ;
        }
        model->addCutGenerator(gen, howOften[action], "Probing") ;
        switches[genCnt++] = 0 ;
    }
    action = ctlBlk->getGomory(gen) ;
    if (action != CbcGenCtlBlk::CGOff) {
        model->addCutGenerator(gen, howOften[action], "Gomory") ;
        switches[genCnt++] = -1 ;
    }
    action = ctlBlk->getKnapsack(gen) ;
    if (action != CbcGenCtlBlk::CGOff) {
        model->addCutGenerator(gen, howOften[action], "Knapsack") ;
        switches[genCnt++] = 0 ;
    }
    action = ctlBlk->getRedSplit(gen) ;
    if (action != CbcGenCtlBlk::CGOff) {
        model->addCutGenerator(gen, howOften[action], "RedSplit") ;
        switches[genCnt++] = 1 ;
    }
    action = ctlBlk->getClique(gen) ;
    if (action != CbcGenCtlBlk::CGOff) {
        model->addCutGenerator(gen, howOften[action], "Clique") ;
        switches[genCnt++] = 0 ;
    }
    action = ctlBlk->getMir(gen) ;
    if (action != CbcGenCtlBlk::CGOff) {
        model->addCutGenerator(gen, howOften[action], "MIR2") ;
        switches[genCnt++] = -1 ;
    }
    action = ctlBlk->getFlow(gen) ;
    if (action != CbcGenCtlBlk::CGOff) {
        model->addCutGenerator(gen, howOften[action], "Flow") ;
        switches[genCnt++] = 1 ;
    }
    action = ctlBlk->getTwomir(gen) ;
    if (action != CbcGenCtlBlk::CGOff) {
        model->addCutGenerator(gen, howOften[action], "2-MIR") ;
        switches[genCnt++] = 1 ;
    }
    /*
      Set control parameters on cut generators. cutDepth says `use this generator
      when (depth in tree) mod cutDepth == 0'. setSwitchOffIfLessThan says `switch
      this generator off if the number of cuts at the root is less than the given
      value'. Sort of. I need to document the magic numbers for howOften , etc.
    */
    genCnt = model->numberCutGenerators() ;
    int iGen ;
    for (iGen = 0 ; iGen < genCnt ; iGen++) {
        CbcCutGenerator *generator = model->cutGenerator(iGen) ;
        int howOften = generator->howOften() ;
        if (howOften == -98 || howOften == -99) {
            generator->setSwitchOffIfLessThan(switches[iGen]) ;
        }
        generator->setTiming(true) ;
        int cutDepth = ctlBlk->getCutDepth() ;
        if (cutDepth >= 0) {
            generator->setWhatDepth(cutDepth) ;
        }
    }
    /*
      Now some additional control parameters that affect cut generation activity.

      Minimum drop is the minimum objective degradation required to continue with
      cut passes.  We want at least .05 unless the objective is tiny, in which
      case we'll drop down to a floor of .0001.
    */
    {
        double objFrac = fabs(model->getMinimizationObjValue()) * .001 + .0001 ;
        double minDrop = CoinMin(.05, objFrac) ;
        model->setMinimumDrop(minDrop) ;
    }
    /*
      Set the maximum number of rounds of cut generation at the root and at nodes
      in the tree. If the value is positive, cut generation will terminate early
      if the objective degradation doesn't meet the minimum drop requirement. If
      the value is negatie, minimum drop is not considered.

      At the root, for small problems, push for 100 passes (really we're betting
      that we'll stop because no cuts were generated). For medium size problems,
      the same but say we can quit if we're not achieving the minimum drop.  For
      big problems, cut the number of rounds to 20.  The user may have expressed
      an opinion; if so, it's already set.

      Once we're in the tree, aim for one pass per activation.
    */
    if (ctlBlk->setByUser_[CbcCbcParam::CUTPASS] == false) {
        int numCols = model->getNumCols() ;
        if (numCols < 500)
            model->setMaximumCutPassesAtRoot(-100) ;
        else if (numCols < 5000)
            model->setMaximumCutPassesAtRoot(100) ;
        else
            model->setMaximumCutPassesAtRoot(20) ;
    }

    model->setMaximumCutPasses(1) ;

    return ;
}

/*
  Install `objects' (integers, SOS sets, etc.) in the OSI. Cribbed from
  CoinSolve 061216 and subjected to moderate rewriting. A substantial amount
  of code that's only relevant for AMPL has been deleted. We're only supporting
  OsiObjects in cbc-generic.
*/

void setupObjects (OsiSolverInterface *osi,
                   bool didIPP, CglPreProcess *ippObj)

{
    int numInts = osi->getNumIntegers() ;
    int numObjs = osi->numberObjects() ;
    /*
      Does this OSI have defined objects already? If not, we'd best define the
      basic integer objects.
    */
    if (numInts == 0 || numObjs == 0) {
        osi->findIntegers(false) ;
        numInts = osi->getNumIntegers() ;
        numObjs = osi->numberObjects() ;
    }
    /*
      If we did preprocessing and discovered SOS sets, create SOS objects and
      install them in the OSI. The priority of SOS objects is set so that larger
      sets have higher (lower numeric value) priority. The priority of the
      original objects is reset to be lower than the priority of any SOS object.
      Since the SOS objects are copied into the OSI, we need to delete our
      originals once they've been installed.

      It's not clear to me that this is the right thing to do, particularly if
      the OSI comes equipped with complex objects.  -- lh, 061216 --
    */
    if (didIPP && ippObj->numberSOS()) {
        OsiObject **oldObjs = osi->objects() ;
        int numCols = osi->getNumCols() ;

        for (int iObj = 0 ; iObj < numObjs ; iObj++) {
            oldObjs[iObj]->setPriority(numCols + 1) ;
        }

        int numSOS = ippObj->numberSOS() ;
        OsiObject **sosObjs = new OsiObject *[numSOS] ;
        const int *starts = ippObj->startSOS() ;
        const int *which = ippObj->whichSOS() ;
        const int *type = ippObj->typeSOS() ;
        const double *weight = ippObj->weightSOS() ;
        int iSOS ;
        for (iSOS = 0 ; iSOS < numSOS ; iSOS++) {
            int iStart = starts[iSOS] ;
            int sosLen = starts[iSOS+1] - iStart ;
            sosObjs[iSOS] =
                new OsiSOS(osi, sosLen, which + iStart, weight + iStart, type[iSOS]) ;
            sosObjs[iSOS]->setPriority(numCols - sosLen) ;
        }
        osi->addObjects(numSOS, sosObjs) ;

        for (iSOS = 0 ; iSOS < numSOS ; iSOS++)
            delete sosObjs[iSOS] ;
        delete [] sosObjs ;
    }

    return ;
}

} // end local namespace


namespace CbcGenParamUtils {

/*
  Run branch-and-cut.
*/

int doBaCParam (CoinParam *param)

{
    assert (param != 0) ;
    CbcGenParam *genParam = dynamic_cast<CbcGenParam *>(param) ;
    assert (genParam != 0) ;
    CbcGenCtlBlk *ctlBlk = genParam->obj() ;
    assert (ctlBlk != 0) ;
    CbcModel *model = ctlBlk->model_ ;
    assert (model != 0) ;
    /*
      Setup to return nonfatal/fatal error (1/-1) by default.
    */
    int retval ;
    if (CoinParamUtils::isInteractive()) {
        retval = 1 ;
    } else {
        retval = -1 ;
    }
    ctlBlk->setBaBStatus(CbcGenCtlBlk::BACAbandon, CbcGenCtlBlk::BACmInvalid,
                         CbcGenCtlBlk::BACwNotStarted, false, 0) ;
    /*
      We ain't gonna do squat without a good model.
    */
    if (!ctlBlk->goodModel_) {
        std::cout << "** Current model not valid!" << std::endl ;
        return (retval) ;
    }
    /*
      Start the clock ticking.
    */
    double time1 = CoinCpuTime() ;
    double time2 ;
    /*
      Create a clone of the model which we can modify with impunity. Extract
      the underlying solver for convenient access.
    */
    CbcModel babModel(*model) ;
    OsiSolverInterface *babSolver = babModel.solver() ;
    assert (babSolver != 0) ;
# if CBC_TRACK_SOLVERS > 0
    std::cout
        << "doBaCParam: initial babSolver is "
        << std::hex << babSolver << std::dec
        << ", log level " << babSolver->messageHandler()->logLevel()
        << "." << std::endl ;
# endif
    /*
      Solve the root relaxation. Bail unless it solves to optimality.
    */
    if (!solveRelaxation(&babModel)) {
        ctlBlk->setBaBStatus(&babModel, CbcGenCtlBlk::BACwBareRoot) ;
        return (0) ;
    }
# if COIN_CBC_VERBOSITY > 0
    std::cout
        << "doBaCParam: initial relaxation z = "
        << babSolver->getObjValue() << "." << std::endl ;
# endif
    /*
      Are we up for fixing variables based on reduced cost alone?
    */
    if (ctlBlk->djFix_.action_ == true) {
        reducedCostHack(babSolver, ctlBlk->djFix_.threshold_) ;
    }
    /*
      Time to consider preprocessing. We'll do a bit of setup before getting to
      the meat of the issue.

      preIppSolver will hold a clone of the unpreprocessed constraint system.
      We'll need it when we postprocess. ippSolver holds the preprocessed
      constraint system.  Again, we clone it and give the clone to babModel for
      B&C. Presumably we need an unmodified copy of the preprocessed system to
      do postprocessing, but the copy itself is hidden inside the preprocess
      object.
    */
    OsiSolverInterface *preIppSolver = 0 ;
    CglPreProcess ippObj ;
    bool didIPP = false ;

    int numberChanged = 0 ;
    int numberOriginalColumns = babSolver->getNumCols() ;
    CbcGenCtlBlk::IPPControl ippAction = ctlBlk->getIPPAction() ;

    if (!(ippAction == CbcGenCtlBlk::IPPOff ||
            ippAction == CbcGenCtlBlk::IPPStrategy)) {
        double timeLeft = babModel.getMaximumSeconds() ;
        preIppSolver = babSolver->clone() ;
        OsiSolverInterface *ippSolver ;
#   if CBC_TRACK_SOLVERS > 0
        std::cout
            << "doBaCParam: clone made prior to IPP is "
            << std::hex << preIppSolver << std::dec
            << ", log level " << preIppSolver->messageHandler()->logLevel()
            << "." << std::endl ;
#   endif

        preIppSolver->setHintParam(OsiDoInBranchAndCut, true, OsiHintDo) ;
        ippObj.messageHandler()->setLogLevel(babModel.logLevel()) ;

        CglProbing probingGen ;
        probingGen.setUsingObjective(true) ;
        probingGen.setMaxPass(3) ;
        probingGen.setMaxProbeRoot(preIppSolver->getNumCols()) ;
        probingGen.setMaxElements(100) ;
        probingGen.setMaxLookRoot(50) ;
        probingGen.setRowCuts(3) ;
        ippObj.addCutGenerator(&probingGen) ;
        /*
          For preProcessNonDefault, the 2nd parameter controls the conversion of
          clique and SOS constraints. 0 does nothing, -1 converts <= to ==, and
          2 and 3 form SOS sets under strict and not-so-strict conditions,
          respectively.
        */
        int convert = 0 ;
        if (ippAction == CbcGenCtlBlk::IPPEqual) {
            convert = -1 ;
        } else if (ippAction == CbcGenCtlBlk::IPPEqualAll) {
            convert = -2 ;
        } else if (ippAction == CbcGenCtlBlk::IPPSOS) {
            convert = 2 ;
        } else if (ippAction == CbcGenCtlBlk::IPPTrySOS) {
            convert = 3 ;
        }

        ippSolver = ippObj.preProcessNonDefault(*preIppSolver, convert, 10) ;
#   if CBC_TRACK_SOLVERS > 0
        std::cout
            << "doBaCParam: solver returned from IPP is "
            << std::hex << ippSolver << std::dec ;
        if (ippSolver) {
            std::cout
                << ", log level " << ippSolver->messageHandler()->logLevel() ;
        }
        std::cout << "." << std::endl ;
#   endif
        /*
          ippSolver == 0 is success of a sort --- integer preprocess has found the
          problem to be infeasible or unbounded. Need to think about how to indicate
          status.
        */
        if (!ippSolver) {
            std::cout
                << "Integer preprocess says infeasible or unbounded" << std::endl ;
            delete preIppSolver ;
            ctlBlk->setBaBStatus(&babModel, CbcGenCtlBlk::BACwIPP) ;
            return (0) ;
        }
#   if COIN_CBC_VERBOSITY > 0
        else {
            std::cout
                << "After integer preprocessing, model has "
                << ippSolver->getNumRows()
                << " rows, " << ippSolver->getNumCols() << " columns, and "
                << ippSolver->getNumElements() << " elements." << std::endl ;
        }
#   endif

        preIppSolver->setHintParam(OsiDoInBranchAndCut, false, OsiHintDo) ;
        ippSolver->setHintParam(OsiDoInBranchAndCut, false, OsiHintDo) ;

        if (ippAction == CbcGenCtlBlk::IPPSave) {
            ippSolver->writeMps("presolved", "mps", 1.0) ;
            std::cout
                << "Integer preprocessed model written to `presolved.mps' "
                << "as minimisation problem." << std::endl ;
        }

        OsiSolverInterface *osiTmp = ippSolver->clone() ;
        babModel.assignSolver(osiTmp) ;
        babSolver = babModel.solver() ;
#   if CBC_TRACK_SOLVERS > 0
        std::cout
            << "doBaCParam: clone of IPP solver passed to babModel is "
            << std::hex << babSolver << std::dec
            << ", log level " << babSolver->messageHandler()->logLevel()
            << "." << std::endl ;
#   endif
        if (!solveRelaxation(&babModel)) {
            delete preIppSolver ;
            ctlBlk->setBaBStatus(&babModel, CbcGenCtlBlk::BACwIPPRelax) ;
            return (0) ;
        }
#   if COIN_CBC_VERBOSITY > 0
        std::cout
            << "doBaCParam: presolved relaxation z = "
            << babSolver->getObjValue() << "." << std::endl ;
#   endif
        babModel.setMaximumSeconds(timeLeft - (CoinCpuTime() - time1)) ;
        didIPP = true ;
    }
    /*
      At this point, babModel and babSolver hold the constraint system we'll use
      for B&C (either the original system or the preprocessed system) and we have
      a solution to the lp relaxation.

      If we're using the COSTSTRATEGY option, set up priorities here and pass
      them to the babModel.
    */
    if (ctlBlk->priorityAction_ != CbcGenCtlBlk::BPOff) {
        setupPriorities(&babModel, ctlBlk->priorityAction_) ;
    }
    /*
      Install heuristics and cutting planes.
    */
    installHeuristics(ctlBlk, &babModel) ;
    installCutGenerators(ctlBlk, &babModel) ;
    /*
      Set up status print frequency for babModel.
    */
    if (babModel.getNumCols() > 2000 || babModel.getNumRows() > 1500 ||
            babModel.messageHandler()->logLevel() > 1)
        babModel.setPrintFrequency(100) ;
    /*
      If we've read in a known good solution for debugging, activate the row cut
      debugger.
    */
    if (ctlBlk->debugSol_.values_) {
        if (ctlBlk->debugSol_.numCols_ == babModel.getNumCols()) {
            babSolver->activateRowCutDebugger(ctlBlk->debugSol_.values_) ;
        } else {
            std::cout
                << "doBaCParam: debug file has incorrect number of columns."
                << std::endl ;
        }
    }
    /*
      Set ratio-based integrality gap, if specified by user.
    */
    if (ctlBlk->setByUser_[CbcCbcParam::GAPRATIO] == true) {
        double obj = babSolver->getObjValue() ;
        double gapRatio = babModel.getDblParam(CbcModel::CbcAllowableFractionGap) ;
        double gap = gapRatio * (1.0e-5 + fabs(obj)) ;
        babModel.setAllowableGap(gap) ;
        std::cout
            << "doBaCParam: Continuous objective = " << obj
            << ", so allowable gap set to " << gap << std::endl ;
    }
    /*
      A bit of mystery code. As best I can figure, setSpecialOptions(2) suppresses
      the removal of warm start information when checkSolution runs an lp to check
      a solution. John's comment, ``probably faster to use a basis to get integer
      solutions'' makes some sense in this context. Didn't try to track down
      moreMipOptions just yet.
    */
    babModel.setSpecialOptions(babModel.specialOptions() | 2) ;
    /*
      { int ndx = whichParam(MOREMIPOPTIONS,numberParameters,parameters) ;
        int moreMipOptions = parameters[ndx].intValue() ;
        if (moreMipOptions >= 0)
        { printf("more mip options %d\n",moreMipOptions);
          babModel.setSearchStrategy(moreMipOptions); } }
    */
    /*
      Begin the final run-up to branch-and-cut.

      Make sure that objects are set up in the solver. It's possible that whoever
      loaded the model into the solver also set up objects. But it's also
      entirely likely that none exist to this point (and interesting to note that
      IPP doesn't need to know anything about objects).
    */
    setupObjects(babSolver, didIPP, &ippObj) ;
    /*
      Set the branching method. We can't do this until we establish objects,
      because the constructor will set up arrays based on the number of objects,
      and there's no provision to set this information after creation. Arguably not
      good --- it'd be nice to set this in the prototype model that's cloned for
      this routine. In CoinSolve, shadowPriceMode is handled with the TESTOSI
      option.
    */
    OsiChooseStrong strong(babSolver) ;
    strong.setNumberBeforeTrusted(babModel.numberBeforeTrust()) ;
    strong.setNumberStrong(babModel.numberStrong()) ;
    strong.setShadowPriceMode(ctlBlk->chooseStrong_.shadowPriceMode_) ;
    CbcBranchDefaultDecision decision ;
    decision.setChooseMethod(strong) ;
    babModel.setBranchingMethod(decision) ;
    /*
      Here I've deleted a huge block of code that deals with external priorities,
      branch direction, pseudocosts, and solution. (PRIORITYIN) Also a block of
      code that generates C++ code.
    */
    /*
      Set up strategy for branch-and-cut. Note that the integer code supplied to
      setupPreProcessing is *not* compatible with the IPPAction enum. But at least
      it's documented. See desiredPreProcess_ in CbcStrategyDefault. `1' is
      accidentally equivalent to IPPOn.
    */

    if (ippAction == CbcGenCtlBlk::IPPStrategy) {
        CbcStrategyDefault strategy(true, 5, 5) ;
        strategy.setupPreProcessing(1) ;
        babModel.setStrategy(strategy) ;
    }
    /*
      Yes! At long last, we're ready for the big call. Do branch and cut. In
      general, the solver used to return the solution will not be the solver we
      passed in, so reset babSolver here.
    */
    int statistics = (ctlBlk->printOpt_ > 0) ? ctlBlk->printOpt_ : 0 ;
# if CBC_TRACK_SOLVERS > 0
    std::cout
        << "doBaCParam: solver at call to branchAndBound is "
        << std::hex << babModel.solver() << std::dec
        << ", log level " << babModel.solver()->messageHandler()->logLevel()
        << "." << std::endl ;
# endif
    babModel.branchAndBound(statistics) ;
    babSolver = babModel.solver() ;
# if CBC_TRACK_SOLVERS > 0
    std::cout
        << "doBaCParam: solver at return from branchAndBound is "
        << std::hex << babModel.solver() << std::dec
        << ", log level " << babModel.solver()->messageHandler()->logLevel()
        << "." << std::endl ;
# endif
    /*
      Write out solution to preprocessed model.
    */
    if (ctlBlk->debugCreate_ == "createAfterPre" &&
            babModel.bestSolution()) {
        CbcGenParamUtils::saveSolution(babSolver, "debug.file") ;
    }
    /*
      Print some information about branch-and-cut.
    */
# if COIN_CBC_VERBOSITY > 0
    std::cout
        << "Cuts at root node changed objective from "
        << babModel.getContinuousObjective()
        << " to " << babModel.rootObjectiveAfterCuts() << std::endl ;

    for (int iGen = 0 ; iGen < babModel.numberCutGenerators() ; iGen++) {
        CbcCutGenerator *generator = babModel.cutGenerator(iGen) ;
        std::cout
            << generator->cutGeneratorName() << " was tried "
            << generator->numberTimesEntered() << " times and created "
            << generator->numberCutsInTotal() << " cuts of which "
            << generator->numberCutsActive()
            << " were active after adding rounds of cuts" ;
        if (generator->timing()) {
            std::cout << " ( " << generator->timeInCutGenerator() << " seconds)" ;
        }
        std::cout << std::endl ;
    }
# endif

    time2 = CoinCpuTime();
    ctlBlk->totalTime_ += time2 - time1;
    /*
      If we performed integer preprocessing, time to back it out.
    */
    if (ippAction != CbcGenCtlBlk::IPPOff) {
#   if CBC_TRACK_SOLVERS > 0
        std::cout
            << "doBaCParam: solver passed to IPP postprocess is "
            << std::hex << babSolver << std::dec << "." << std::endl ;
#   endif
        ippObj.postProcess(*babSolver);
        babModel.assignSolver(preIppSolver) ;
        babSolver = babModel.solver() ;
#   if CBC_TRACK_SOLVERS > 0
        std::cout
            << "doBaCParam: solver in babModel after IPP postprocess is "
            << std::hex << babSolver << std::dec << "." << std::endl ;
#   endif
    }
    /*
      Write out postprocessed solution to debug file, if requested.
    */
    if (ctlBlk->debugCreate_ == "create" && babModel.bestSolution()) {
        CbcGenParamUtils::saveSolution(babSolver, "debug.file") ;
    }
    /*
      If we have a good solution, detach the solver with the answer. Fill in the
      rest of the status information for the benefit of the wider world.
    */
    bool keepAnswerSolver = false ;
    OsiSolverInterface *answerSolver = 0 ;
    if (babModel.bestSolution()) {
        babModel.setModelOwnsSolver(false) ;
        keepAnswerSolver = true ;
        answerSolver = babSolver ;
    }
    ctlBlk->setBaBStatus(&babModel, CbcGenCtlBlk::BACwBAC,
                         keepAnswerSolver, answerSolver) ;
    /*
      And one last bit of information & statistics.
    */
    ctlBlk->printBaBStatus() ;
    std::cout << "    " ;
    if (keepAnswerSolver) {
        std::cout
            << "objective " << babModel.getObjValue() << "; " ;
    }
    std::cout
        << babModel.getNodeCount() << " nodes and "
        << babModel.getIterationCount() << " iterations - took "
        << time2 - time1 << " seconds" << std::endl ;

    return (0) ;
}

} // end namespace CbcGenParamutils

