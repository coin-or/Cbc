/* $Id$ */
// Copyright (C) 2002, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

#if defined(_MSC_VER)
// Turn off compiler warning about long names
#pragma warning(disable : 4786)
#include <windows.h> // for Sleep()
#ifdef small
#undef small
#endif
#else
#include <unistd.h> // for usleep()
#endif

#include "CbcConfig.h"
#ifdef COIN_HAS_NTY
#include "CbcSymmetry.hpp"
#endif
//#define DEBUG_SOLUTION
#ifdef DEBUG_SOLUTION
#define COIN_DETAIL
#endif
#include <string>
//#define CBC_DEBUG 1
//#define CHECK_CUT_COUNTS
//#define CHECK_NODE
//#define CBC_CHECK_BASIS
#include <cassert>
#include <cfloat>
#define CUTS
#include "OsiSolverInterface.hpp"
#include "OsiChooseVariable.hpp"
#include "OsiAuxInfo.hpp"
#include "OsiSolverBranch.hpp"
#include "CoinWarmStartBasis.hpp"
#include "CoinTime.hpp"
#include "CbcModel.hpp"
#include "CbcNode.hpp"
#include "CbcStatistics.hpp"
#include "CbcStrategy.hpp"
#include "CbcBranchActual.hpp"
#include "CbcBranchDynamic.hpp"
#include "OsiRowCut.hpp"
#include "OsiRowCutDebugger.hpp"
#include "OsiCuts.hpp"
#include "CbcCountRowCut.hpp"
#include "CbcFeasibilityBase.hpp"
#include "CbcMessage.hpp"
#ifdef COIN_HAS_CLP
#include "OsiClpSolverInterface.hpp"
#include "ClpSimplexOther.hpp"
#include "ClpSolve.hpp"
#include "ClpDualRowSteepest.hpp"
#include "ClpPrimalColumnPivot.hpp"
#endif
using namespace std;
#include "CglCutGenerator.hpp"

CbcNode::CbcNode()
  : nodeInfo_(NULL)
  , objectiveValue_(1.0e100)
  , guessedObjectiveValue_(1.0e100)
  , sumInfeasibilities_(0.0)
  , branch_(NULL)
  , depth_(-1)
  , numberUnsatisfied_(0)
  , nodeNumber_(-1)
  , state_(0)
{
#ifdef CHECK_NODE
  printf("CbcNode %p Constructor\n", this);
#endif
}
// Print
void CbcNode::print() const
{
  printf("number %d obj %g depth %d sumun %g nunsat %d state %d\n",
    nodeNumber_, objectiveValue_, depth_, sumInfeasibilities_, numberUnsatisfied_, state_);
}
CbcNode::CbcNode(CbcModel *model,
  CbcNode *lastNode)
  : nodeInfo_(NULL)
  , objectiveValue_(1.0e100)
  , guessedObjectiveValue_(1.0e100)
  , sumInfeasibilities_(0.0)
  , branch_(NULL)
  , depth_(-1)
  , numberUnsatisfied_(0)
  , nodeNumber_(-1)
  , state_(0)
{
#ifdef CHECK_NODE
  printf("CbcNode %p Constructor from model\n", this);
#endif
  model->setObjectiveValue(this, lastNode);

  if (lastNode) {
    if (lastNode->nodeInfo_) {
      lastNode->nodeInfo_->increment();
    }
  }
  nodeNumber_ = model->getNodeCount();
}

#define CBC_NEW_CREATEINFO
#ifdef CBC_NEW_CREATEINFO

/*
  New createInfo, with basis manipulation hidden inside mergeBasis. Allows
  solvers to override and carry over all information from one basis to
  another.
*/

void CbcNode::createInfo(CbcModel *model,
  CbcNode *lastNode,
  const CoinWarmStartBasis *lastws,
  const double *lastLower, const double *lastUpper,
  int numberOldActiveCuts, int numberNewCuts)

{
  OsiSolverInterface *solver = model->solver();
  CbcStrategy *strategy = model->strategy();
  /*
      The root --- no parent. Create full basis and bounds information.
    */
  if (!lastNode) {
    if (!strategy)
      nodeInfo_ = new CbcFullNodeInfo(model, solver->getNumRows());
    else
      nodeInfo_ = strategy->fullNodeInfo(model, solver->getNumRows());
  } else {
    /*
          Not the root. Create an edit from the parent's basis & bound information.
          This is not quite as straightforward as it seems. We need to reintroduce
          cuts we may have dropped out of the basis, in the correct position, because
          this whole process is strictly positional. Start by grabbing the current
          basis.
        */
    bool mustDeleteBasis;
    const CoinWarmStartBasis *ws = dynamic_cast< const CoinWarmStartBasis * >(solver->getPointerToWarmStart(mustDeleteBasis));
    assert(ws != NULL); // make sure not volume
    //int numberArtificials = lastws->getNumArtificial();
    int numberColumns = solver->getNumCols();
    int numberRowsAtContinuous = model->numberRowsAtContinuous();
    int currentNumberCuts = model->currentNumberCuts();
#ifdef CBC_CHECK_BASIS
    std::cout
      << "Before expansion: orig " << numberRowsAtContinuous
      << ", old " << numberOldActiveCuts
      << ", new " << numberNewCuts
      << ", current " << currentNumberCuts << "." << std::endl;
    ws->print();
#endif
    /*
          Clone the basis and resize it to hold the structural constraints, plus
          all the cuts: old cuts, both active and inactive (currentNumberCuts),
          and new cuts (numberNewCuts). This will become the expanded basis.
        */
    CoinWarmStartBasis *expanded = dynamic_cast< CoinWarmStartBasis * >(ws->clone());
    int iCompact = numberRowsAtContinuous + numberOldActiveCuts + numberNewCuts;
    // int nPartial = numberRowsAtContinuous+currentNumberCuts;
    int iFull = numberRowsAtContinuous + currentNumberCuts + numberNewCuts;
    // int maxBasisLength = ((iFull+15)>>4)+((numberColumns+15)>>4);
    // printf("l %d full %d\n",maxBasisLength,iFull);
    expanded->resize(iFull, numberColumns);
#ifdef CBC_CHECK_BASIS
    std::cout
      << "\tFull basis " << iFull << " rows, "
      << numberColumns << " columns; compact "
      << iCompact << " rows." << std::endl;
#endif
    /*
          Now flesh out the expanded basis. The clone already has the
          correct status information for the variables and for the structural
          (numberRowsAtContinuous) constraints. Any indices beyond nPartial must be
          cuts created while processing this node --- they can be copied en bloc
          into the correct position in the expanded basis. The space reserved for
          xferRows is a gross overestimate.
        */
    CoinWarmStartBasis::XferVec xferRows;
    xferRows.reserve(iFull - numberRowsAtContinuous + 1);
    if (numberNewCuts) {
      xferRows.push_back(
        CoinWarmStartBasis::XferEntry(iCompact - numberNewCuts,
          iFull - numberNewCuts, numberNewCuts));
    }
    /*
          From nPartial down, record the entries we want to copy from the current
          basis (the entries for the active cuts; non-zero in the list returned
          by addedCuts). Fill the expanded basis with entries showing a status of
          basic for the deactivated (loose) cuts.
        */
    CbcCountRowCut **cut = model->addedCuts();
    iFull -= (numberNewCuts + 1);
    iCompact -= (numberNewCuts + 1);
    int runLen = 0;
    CoinWarmStartBasis::XferEntry entry(-1, -1, -1);
    while (iFull >= numberRowsAtContinuous) {
      for (; iFull >= numberRowsAtContinuous && cut[iFull - numberRowsAtContinuous]; iFull--)
        runLen++;
      if (runLen) {
        iCompact -= runLen;
        entry.first = iCompact + 1;
        entry.second = iFull + 1;
        entry.third = runLen;
        runLen = 0;
        xferRows.push_back(entry);
      }
      for (; iFull >= numberRowsAtContinuous && !cut[iFull - numberRowsAtContinuous]; iFull--)
        expanded->setArtifStatus(iFull, CoinWarmStartBasis::basic);
    }
    /*
          Finally, call mergeBasis to copy over entries from the current basis to
          the expanded basis. Since we cloned the expanded basis from the active basis
          and haven't changed the number of variables, only row status entries need
          to be copied.
        */
    expanded->mergeBasis(ws, &xferRows, 0);

#ifdef CBC_CHECK_BASIS
    std::cout << "Expanded basis:" << std::endl;
    expanded->print();
    std::cout << "Diffing against:" << std::endl;
    lastws->print();
#endif
    assert(expanded->getNumArtificial() >= lastws->getNumArtificial());
#ifdef CLP_INVESTIGATE
    if (!expanded->fullBasis()) {
      int iFull = numberRowsAtContinuous + currentNumberCuts + numberNewCuts;
      printf("cont %d old %d new %d current %d full inc %d full %d\n",
        numberRowsAtContinuous, numberOldActiveCuts, numberNewCuts,
        currentNumberCuts, iFull, iFull - numberNewCuts);
    }
#endif

    /*
          Now that we have two bases in proper positional correspondence, creating
          the actual diff is dead easy.

          Note that we're going to compare the expanded basis here to the stripped
          basis (lastws) produced by addCuts. It doesn't affect the correctness (the
          diff process has no knowledge of the meaning of an entry) but it does
          mean that we'll always generate a whack of diff entries because the expanded
          basis is considerably larger than the stripped basis.
        */
    CoinWarmStartDiff *basisDiff = expanded->generateDiff(lastws);

    /*
          Diff the bound vectors. It's assumed the number of structural variables
          is not changing. For branching objects that change bounds on integer
          variables, we should see at least one bound change as a consequence
          of applying the branch that generated this subproblem from its parent.
          This need not hold for other types of branching objects (hyperplane
          branches, for example).
        */
    const double *lower = solver->getColLower();
    const double *upper = solver->getColUpper();

    double *boundChanges = new double[2 * numberColumns];
    int *variables = new int[2 * numberColumns];
    int numberChangedBounds = 0;

    int i;
    for (i = 0; i < numberColumns; i++) {
      if (lower[i] != lastLower[i]) {
        variables[numberChangedBounds] = i;
        boundChanges[numberChangedBounds++] = lower[i];
      }
      if (upper[i] != lastUpper[i]) {
        variables[numberChangedBounds] = i | 0x80000000;
        boundChanges[numberChangedBounds++] = upper[i];
      }
#ifdef CBC_DEBUG
      if (lower[i] != lastLower[i]) {
        std::cout
          << "lower on " << i << " changed from "
          << lastLower[i] << " to " << lower[i] << std::endl;
      }
      if (upper[i] != lastUpper[i]) {
        std::cout
          << "upper on " << i << " changed from "
          << lastUpper[i] << " to " << upper[i] << std::endl;
      }
#endif
    }
#ifdef CBC_DEBUG
    std::cout << numberChangedBounds << " changed bounds." << std::endl;
#endif
    //if (lastNode->branchingObject()->boundBranch())
    //assert (numberChangedBounds);
    /*
          Hand the lot over to the CbcPartialNodeInfo constructor, then clean up and
          return.
        */
    if (!strategy) {
      delete nodeInfo_;
      nodeInfo_ = new CbcPartialNodeInfo(lastNode->nodeInfo_, this, numberChangedBounds,
        variables, boundChanges, basisDiff);
    } else {
      nodeInfo_ = strategy->partialNodeInfo(model, lastNode->nodeInfo_, this,
        numberChangedBounds, variables, boundChanges,
        basisDiff);
    }
    delete basisDiff;
    delete[] boundChanges;
    delete[] variables;
    delete expanded;
    if (mustDeleteBasis)
      delete ws;
  }
  // Set node number
  nodeInfo_->setNodeNumber(model->getNodeCount2());
  state_ |= 2; // say active
}

#else // CBC_NEW_CREATEINFO

/*
  Original createInfo, with bare manipulation of basis vectors. Fails if solver
  maintains additional information in basis.
*/

void CbcNode::createInfo(CbcModel *model,
  CbcNode *lastNode,
  const CoinWarmStartBasis *lastws,
  const double *lastLower, const double *lastUpper,
  int numberOldActiveCuts, int numberNewCuts)
{
  OsiSolverInterface *solver = model->solver();
  CbcStrategy *strategy = model->strategy();
  /*
      The root --- no parent. Create full basis and bounds information.
    */
  if (!lastNode) {
    if (!strategy)
      nodeInfo_ = new CbcFullNodeInfo(model, solver->getNumRows());
    else
      nodeInfo_ = strategy->fullNodeInfo(model, solver->getNumRows());
  }
  /*
      Not the root. Create an edit from the parent's basis & bound information.
      This is not quite as straightforward as it seems. We need to reintroduce
      cuts we may have dropped out of the basis, in the correct position, because
      this whole process is strictly positional. Start by grabbing the current
      basis.
    */
  else {
    bool mustDeleteBasis;
    const CoinWarmStartBasis *ws = dynamic_cast< const CoinWarmStartBasis * >(solver->getPointerToWarmStart(mustDeleteBasis));
    assert(ws != NULL); // make sure not volume
    //int numberArtificials = lastws->getNumArtificial();
    int numberColumns = solver->getNumCols();

    const double *lower = solver->getColLower();
    const double *upper = solver->getColUpper();

    int i;
    /*
        Create a clone and resize it to hold all the structural constraints, plus
        all the cuts: old cuts, both active and inactive (currentNumberCuts), and
        new cuts (numberNewCuts).

        TODO: You'd think that the set of constraints (logicals) in the expanded
        basis should match the set represented in lastws. At least, that's
        what I thought. But at the point I first looked hard at this bit of
        code, it turned out that lastws was the stripped basis produced at
        the end of addCuts(), rather than the raw basis handed back by
        addCuts1(). The expanded basis here is equivalent to the raw basis of
        addCuts1(). I said ``whoa, that's not good, I must have introduced a
        bug'' and went back to John's code to see where I'd gone wrong.
        And discovered the same `error' in his code.

        After a bit of thought, my conclusion is that correctness is not
        affected by whether lastws is the stripped or raw basis. The diffs
        have no semantics --- just a set of changes that need to be made
        to convert lastws into expanded. I think the only effect is that we
        store a lot more diffs (everything in expanded that's not covered by
        the stripped basis). But I need to give this more thought. There
        may well be some subtle error cases.

        In the mean time, I've twiddled addCuts() to set lastws to the raw
        basis. Makes me (Lou) less nervous to compare apples to apples.
        */
    CoinWarmStartBasis *expanded = dynamic_cast< CoinWarmStartBasis * >(ws->clone());
    int numberRowsAtContinuous = model->numberRowsAtContinuous();
    int iFull = numberRowsAtContinuous + model->currentNumberCuts() + numberNewCuts;
    //int numberArtificialsNow = iFull;
    //int maxBasisLength = ((iFull+15)>>4)+((numberColumns+15)>>4);
    //printf("l %d full %d\n",maxBasisLength,iFull);
    if (expanded)
      expanded->resize(iFull, numberColumns);
#ifdef CBC_CHECK_BASIS
    printf("Before expansion: orig %d, old %d, new %d, current %d\n",
      numberRowsAtContinuous, numberOldActiveCuts, numberNewCuts,
      model->currentNumberCuts());
    ws->print();
#endif
    /*
        Now fill in the expanded basis. Any indices beyond nPartial must
        be cuts created while processing this node --- they can be copied directly
        into the expanded basis. From nPartial down, pull the status of active cuts
        from ws, interleaving with a B entry for the deactivated (loose) cuts.
        */
    int numberDropped = model->currentNumberCuts() - numberOldActiveCuts;
    int iCompact = iFull - numberDropped;
    CbcCountRowCut **cut = model->addedCuts();
    int nPartial = model->currentNumberCuts() + numberRowsAtContinuous;
    iFull--;
    for (; iFull >= nPartial; iFull--) {
      CoinWarmStartBasis::Status status = ws->getArtifStatus(--iCompact);
      //assert (status != CoinWarmStartBasis::basic); // may be permanent cut
      expanded->setArtifStatus(iFull, status);
    }
    for (; iFull >= numberRowsAtContinuous; iFull--) {
      if (cut[iFull - numberRowsAtContinuous]) {
        CoinWarmStartBasis::Status status = ws->getArtifStatus(--iCompact);
        // If no cut generator being used then we may have basic variables
        //if (model->getMaximumCutPasses()&&
        //  status == CoinWarmStartBasis::basic)
        //printf("cut basic\n");
        expanded->setArtifStatus(iFull, status);
      } else {
        expanded->setArtifStatus(iFull, CoinWarmStartBasis::basic);
      }
    }
#ifdef CBC_CHECK_BASIS
    printf("Expanded basis\n");
    expanded->print();
    printf("Diffing against\n");
    lastws->print();
#endif
    /*
        Now that we have two bases in proper positional correspondence, creating
        the actual diff is dead easy.
        */

    CoinWarmStartDiff *basisDiff = expanded->generateDiff(lastws);
    /*
        Diff the bound vectors. It's assumed the number of structural variables is
        not changing. Assuming that branching objects all involve integer variables,
        we should see at least one bound change as a consequence of processing this
        subproblem. Different types of branching objects could break this assertion.
        Not true at all - we have not applied current branch - JJF.
        */
    double *boundChanges = new double[2 * numberColumns];
    int *variables = new int[2 * numberColumns];
    int numberChangedBounds = 0;
    for (i = 0; i < numberColumns; i++) {
      if (lower[i] != lastLower[i]) {
        variables[numberChangedBounds] = i;
        boundChanges[numberChangedBounds++] = lower[i];
      }
      if (upper[i] != lastUpper[i]) {
        variables[numberChangedBounds] = i | 0x80000000;
        boundChanges[numberChangedBounds++] = upper[i];
      }
#ifdef CBC_DEBUG
      if (lower[i] != lastLower[i])
        printf("lower on %d changed from %g to %g\n",
          i, lastLower[i], lower[i]);
      if (upper[i] != lastUpper[i])
        printf("upper on %d changed from %g to %g\n",
          i, lastUpper[i], upper[i]);
#endif
    }
#ifdef CBC_DEBUG
    printf("%d changed bounds\n", numberChangedBounds);
#endif
    //if (lastNode->branchingObject()->boundBranch())
    //assert (numberChangedBounds);
    /*
        Hand the lot over to the CbcPartialNodeInfo constructor, then clean up and
        return.
        */
    if (!strategy)
      nodeInfo_ = new CbcPartialNodeInfo(lastNode->nodeInfo_, this, numberChangedBounds,
        variables, boundChanges, basisDiff);
    else
      nodeInfo_ = strategy->partialNodeInfo(model, lastNode->nodeInfo_, this, numberChangedBounds,
        variables, boundChanges, basisDiff);
    delete basisDiff;
    delete[] boundChanges;
    delete[] variables;
    delete expanded;
    if (mustDeleteBasis)
      delete ws;
  }
  // Set node number
  nodeInfo_->setNodeNumber(model->getNodeCount2());
  state_ |= 2; // say active
}

#endif // CBC_NEW_CREATEINFO
/*
  The routine scans through the object list of the model looking for objects
  that indicate infeasibility. It tests each object using strong branching
  and selects the one with the least objective degradation.  A corresponding
  branching object is left attached to lastNode.

  If strong branching is disabled, a candidate object is chosen essentially
  at random (whatever object ends up in pos'n 0 of the candidate array).

  If a branching candidate is found to be monotone, bounds are set to fix the
  variable and the routine immediately returns (the caller is expected to
  reoptimize).

  If a branching candidate is found to result in infeasibility in both
  directions, the routine immediately returns an indication of infeasibility.

  Returns:  0	both branch directions are feasible
  -1	branching variable is monotone
  -2	infeasible

  Original comments:
  Here could go cuts etc etc
  For now just fix on objective from strong branching.
*/

int CbcNode::chooseBranch(CbcModel *model, CbcNode *lastNode, int numberPassesLeft)

{
  if (lastNode)
    depth_ = lastNode->depth_ + 1;
  else
    depth_ = 0;
  delete branch_;
  branch_ = NULL;
  OsiSolverInterface *solver = model->solver();
  // Mark variables which need to be clean
  char *cleanVariables = NULL;
#ifdef COIN_HAS_CLP
  OsiClpSolverInterface *osiclp = dynamic_cast< OsiClpSolverInterface * >(solver);
  int saveClpOptions = 0;
  if (osiclp) {
    // for faster hot start
    saveClpOptions = osiclp->specialOptions();
    osiclp->setSpecialOptions(saveClpOptions | 8192);
    if ((model->moreSpecialOptions2() & 32768) != 0) {
      cleanVariables = model->setupCleanVariables(); // for odd ints/sos etc
    }
  }
#else
  OsiSolverInterface *osiclp = NULL;
#endif
  double saveObjectiveValue = solver->getObjValue();
  double objectiveValue = CoinMax(solver->getObjSense() * saveObjectiveValue, objectiveValue_);
  const double *lower = solver->getColLower();
  const double *upper = solver->getColUpper();
  // See what user thinks
  int anyAction = model->problemFeasibility()->feasible(model, 0);
  if (anyAction) {
    // will return -2 if infeasible , 0 if treat as integer
    return anyAction - 1;
  }
  double integerTolerance = model->getDblParam(CbcModel::CbcIntegerTolerance);
  // point to useful information
  OsiBranchingInformation usefulInfo = model->usefulInformation();
  // and modify
  usefulInfo.depth_ = depth_;
  int i;
  bool beforeSolution = model->getSolutionCount() == 0;
  int numberStrong = model->numberStrong();
  // switch off strong if hotstart
  const double *hotstartSolution = model->hotstartSolution();
  const int *hotstartPriorities = model->hotstartPriorities();
  int numberObjects = model->numberObjects();
  int numberColumns = model->getNumCols();
  double *saveUpper = new double[numberColumns];
  double *saveLower = new double[numberColumns];
  for (i = 0; i < numberColumns; i++) {
    saveLower[i] = lower[i];
    saveUpper[i] = upper[i];
  }

  // Save solution in case heuristics need good solution later

  double *saveSolution = new double[numberColumns];
  memcpy(saveSolution, solver->getColSolution(), numberColumns * sizeof(double));
  model->reserveCurrentSolution(saveSolution);
  if (hotstartSolution) {
    numberStrong = 0;
    if ((model->moreSpecialOptions() & 1024) != 0 || true) {
      int nBad = 0;
      int nUnsat = 0;
      int nDiff = 0;
      for (int i = 0; i < numberObjects; i++) {
        OsiObject *object = model->modifiableObject(i);
        const CbcSimpleInteger *thisOne = dynamic_cast< const CbcSimpleInteger * >(object);
        if (thisOne) {
          int iColumn = thisOne->columnNumber();
          double targetValue = hotstartSolution[iColumn];
          double value = saveSolution[iColumn];
          if (fabs(value - floor(value + 0.5)) > 1.0e-6) {
            nUnsat++;
#ifdef CLP_INVESTIGATE
            printf("H %d is %g target %g\n", iColumn, value, targetValue);
#endif
          } else if (fabs(targetValue - value) > 1.0e-6) {
            nDiff++;
          }
          if (targetValue < saveLower[iColumn] || targetValue > saveUpper[iColumn]) {
#ifdef CLP_INVESTIGATE
            printf("%d has target %g and current bounds %g and %g\n",
              iColumn, targetValue, saveLower[iColumn], saveUpper[iColumn]);
#endif
            nBad++;
          }
        }
      }
#ifdef CLP_INVESTIGATE
      printf("Hot %d unsatisfied, %d outside limits, %d different\n",
        nUnsat, nBad, nDiff);
#endif
      if (nBad) {
        // switch off as not possible
        hotstartSolution = NULL;
        model->setHotstartSolution(NULL, NULL);
        usefulInfo.hotstartSolution_ = NULL;
      }
    }
  }
  int numberStrongDone = 0;
  int numberUnfinished = 0;
  int numberStrongInfeasible = 0;
  int numberStrongIterations = 0;
  int saveNumberStrong = numberStrong;
  bool checkFeasibility = numberObjects > model->numberIntegers();
  int maximumStrong = CoinMax(CoinMin(numberStrong, numberObjects), 1);
  /*
      Get a branching decision object. Use the default decision criteria unless
      the user has loaded a decision method into the model.
    */
  CbcBranchDecision *decision = model->branchingMethod();
  CbcDynamicPseudoCostBranchingObject *dynamicBranchingObject = dynamic_cast< CbcDynamicPseudoCostBranchingObject * >(decision);
  if (!decision || dynamicBranchingObject)
    decision = new CbcBranchDefaultDecision();
  decision->initialize(model);
  CbcStrongInfo *choice = new CbcStrongInfo[maximumStrong];
  // May go round twice if strong branching fixes all local candidates
  bool finished = false;
  double estimatedDegradation = 0.0;
  while (!finished) {
    finished = true;
    // Some objects may compute an estimate of best solution from here
    estimatedDegradation = 0.0;
    //int numberIntegerInfeasibilities=0; // without odd ones
    numberStrongDone = 0;
    numberUnfinished = 0;
    numberStrongInfeasible = 0;
    numberStrongIterations = 0;

    // We may go round this loop twice (only if we think we have solution)
    for (int iPass = 0; iPass < 2; iPass++) {

      // compute current state
      //int numberObjectInfeasibilities; // just odd ones
      //model->feasibleSolution(
      //                      numberIntegerInfeasibilities,
      //                      numberObjectInfeasibilities);
      // Some objects may compute an estimate of best solution from here
      estimatedDegradation = 0.0;
      numberUnsatisfied_ = 0;
      // initialize sum of "infeasibilities"
      sumInfeasibilities_ = 0.0;
      int bestPriority = COIN_INT_MAX;
      /*
              Scan for branching objects that indicate infeasibility. Choose the best
              maximumStrong candidates, using priority as the first criteria, then
              integer infeasibility.

              The algorithm is to fill the choice array with a set of good candidates (by
              infeasibility) with priority bestPriority.  Finding a candidate with
              priority better (less) than bestPriority flushes the choice array. (This
              serves as initialization when the first candidate is found.)

              A new candidate is added to choices only if its infeasibility exceeds the
              current max infeasibility (mostAway). When a candidate is added, it
              replaces the candidate with the smallest infeasibility (tracked by
              iSmallest).
            */
      int iSmallest = 0;
      double mostAway = 1.0e-100;
      for (i = 0; i < maximumStrong; i++)
        choice[i].possibleBranch = NULL;
      numberStrong = 0;
      bool canDoOneHot = false;
      for (i = 0; i < numberObjects; i++) {
        OsiObject *object = model->modifiableObject(i);
        int preferredWay;
        double infeasibility = object->infeasibility(&usefulInfo, preferredWay);
        int priorityLevel = object->priority();
        if (hotstartSolution) {
          // we are doing hot start
          const CbcSimpleInteger *thisOne = dynamic_cast< const CbcSimpleInteger * >(object);
          if (thisOne) {
            int iColumn = thisOne->columnNumber();
            bool canDoThisHot = true;
            double targetValue = hotstartSolution[iColumn];
            if (saveUpper[iColumn] > saveLower[iColumn]) {
              double value = saveSolution[iColumn];
              // clean
              value = CoinMin(value, saveUpper[iColumn]);
              value = CoinMax(value, saveLower[iColumn]);
              if (hotstartPriorities)
                priorityLevel = hotstartPriorities[iColumn];
              //double originalLower = thisOne->originalLower();
              //double originalUpper = thisOne->originalUpper();
              // switch off if not possible
              if (targetValue >= saveLower[iColumn] && targetValue <= saveUpper[iColumn]) {
                /* priority outranks rest always if negative
                                   otherwise can be downgraded if at correct level.
                                   Infeasibility may be increased to choose 1.0 values first.
                                   choose one near wanted value
                                */
                if (fabs(value - targetValue) > integerTolerance) {
                  //if (infeasibility>0.01)
                  //infeasibility = fabs(1.0e6-fabs(value-targetValue));
                  //else
                  infeasibility = fabs(value - targetValue);
                  //if (targetValue==1.0)
                  //infeasibility += 1.0;
                  if (value > targetValue) {
                    preferredWay = -1;
                  } else {
                    preferredWay = 1;
                  }
                  priorityLevel = CoinAbs(priorityLevel);
                } else if (priorityLevel < 0) {
                  priorityLevel = CoinAbs(priorityLevel);
                  if (targetValue == saveLower[iColumn]) {
                    infeasibility = integerTolerance + 1.0e-12;
                    preferredWay = -1;
                  } else if (targetValue == saveUpper[iColumn]) {
                    infeasibility = integerTolerance + 1.0e-12;
                    preferredWay = 1;
                  } else {
                    // can't
                    priorityLevel += 10000000;
                    canDoThisHot = false;
                  }
                } else {
                  priorityLevel += 10000000;
                  canDoThisHot = false;
                }
              } else {
                // switch off if not possible
                canDoThisHot = false;
              }
              if (canDoThisHot)
                canDoOneHot = true;
            } else if (targetValue < saveLower[iColumn] || targetValue > saveUpper[iColumn]) {
            }
          } else {
            priorityLevel += 10000000;
          }
        }
        if (infeasibility) {
          // Increase estimated degradation to solution
          estimatedDegradation += CoinMin(object->upEstimate(), object->downEstimate());
          numberUnsatisfied_++;
          sumInfeasibilities_ += infeasibility;
          // Better priority? Flush choices.
          if (priorityLevel < bestPriority) {
            int j;
            iSmallest = 0;
            for (j = 0; j < maximumStrong; j++) {
              choice[j].upMovement = 0.0;
              delete choice[j].possibleBranch;
              choice[j].possibleBranch = NULL;
            }
            bestPriority = priorityLevel;
            mostAway = 1.0e-100;
            numberStrong = 0;
          } else if (priorityLevel > bestPriority) {
            continue;
          }
          // Check for suitability based on infeasibility.
          if (infeasibility > mostAway) {
            //add to list
            choice[iSmallest].upMovement = infeasibility;
            delete choice[iSmallest].possibleBranch;
            CbcObject *obj = dynamic_cast< CbcObject * >(object);
            assert(obj);
            choice[iSmallest].possibleBranch = obj->createCbcBranch(solver, &usefulInfo, preferredWay);
            numberStrong = CoinMax(numberStrong, iSmallest + 1);
            // Save which object it was
            choice[iSmallest].objectNumber = i;
            int j;
            iSmallest = -1;
            mostAway = 1.0e50;
            for (j = 0; j < maximumStrong; j++) {
              if (choice[j].upMovement < mostAway) {
                mostAway = choice[j].upMovement;
                iSmallest = j;
              }
            }
          }
        }
      }
      if (!canDoOneHot && hotstartSolution) {
        // switch off as not possible
        hotstartSolution = NULL;
        model->setHotstartSolution(NULL, NULL);
        usefulInfo.hotstartSolution_ = NULL;
      }
      if (numberUnsatisfied_) {
        // some infeasibilities - go to next steps
#ifdef CLP_INVESTIGATE
        if (hotstartSolution) {
          int k = choice[0].objectNumber;
          OsiObject *object = model->modifiableObject(k);
          const CbcSimpleInteger *thisOne = dynamic_cast< const CbcSimpleInteger * >(object);
          assert(thisOne);
          int iColumn = thisOne->columnNumber();
          double targetValue = hotstartSolution[iColumn];
          double value = saveSolution[iColumn];
          printf("Branch on %d has target %g (value %g) and current bounds %g and %g\n",
            iColumn, targetValue, value, saveLower[iColumn], saveUpper[iColumn]);
        }
#endif
        break;
      } else if (!iPass) {
        // looks like a solution - get paranoid
        bool roundAgain = false;
        // get basis
        CoinWarmStartBasis *ws = dynamic_cast< CoinWarmStartBasis * >(solver->getWarmStart());
        if (!ws)
          break;
        for (i = 0; i < numberColumns; i++) {
          double value = saveSolution[i];
          if (value < lower[i]) {
            saveSolution[i] = lower[i];
            roundAgain = true;
            ws->setStructStatus(i, CoinWarmStartBasis::atLowerBound);
          } else if (value > upper[i]) {
            saveSolution[i] = upper[i];
            roundAgain = true;
            ws->setStructStatus(i, CoinWarmStartBasis::atUpperBound);
          }
        }
        if (roundAgain && saveNumberStrong) {
          // restore basis
          solver->setWarmStart(ws);
          delete ws;
          solver->resolve();
          memcpy(saveSolution, solver->getColSolution(), numberColumns * sizeof(double));
          model->reserveCurrentSolution(saveSolution);
          if (!solver->isProvenOptimal()) {
            // infeasible
            anyAction = -2;
            break;
          }
        } else {
          delete ws;
          break;
        }
      }
    }
    /* Some solvers can do the strong branching calculations faster if
           they do them all at once.  At present only Clp does for ordinary
           integers but I think this coding would be easy to modify
        */
    bool allNormal = true; // to say if we can do fast strong branching
    // Say which one will be best
    int bestChoice = 0;
    double worstInfeasibility = 0.0;
    for (i = 0; i < numberStrong; i++) {
      choice[i].numIntInfeasUp = numberUnsatisfied_;
      choice[i].numIntInfeasDown = numberUnsatisfied_;
      choice[i].fix = 0; // say not fixed
      if (!dynamic_cast< const CbcSimpleInteger * >(model->object(choice[i].objectNumber)))
        allNormal = false; // Something odd so lets skip clever fast branching
      if (!model->object(choice[i].objectNumber)->boundBranch())
        numberStrong = 0; // switch off
      if (choice[i].possibleBranch->numberBranches() > 2)
        numberStrong = 0; // switch off
      // Do best choice in case switched off
      if (choice[i].upMovement > worstInfeasibility) {
        worstInfeasibility = choice[i].upMovement;
        bestChoice = i;
      }
    }
    //if (!model->parentModel())
    //solver->writeMps("query");
    // If we have hit max time don't do strong branching
    bool hitMaxTime = (model->getCurrentSeconds() > model->getDblParam(CbcModel::CbcMaximumSeconds));
    // also give up if we are looping round too much
    if (hitMaxTime || numberPassesLeft <= 0)
      numberStrong = 0;
    /*
          Is strong branching enabled? If so, set up and do it. Otherwise, we'll
          fall through to simple branching.

          Setup for strong branching involves saving the current basis (for restoration
          afterwards) and setting up for hot starts.
        */
    if (numberStrong && saveNumberStrong) {

      bool solveAll = false; // set true to say look at all even if some fixed (experiment)
      solveAll = true;
      // worth trying if too many times
      // Save basis
      CoinWarmStart *ws = solver->getWarmStart();
      // save limit
      int saveLimit;
      solver->getIntParam(OsiMaxNumIterationHotStart, saveLimit);
      if (beforeSolution && saveLimit < 100)
        solver->setIntParam(OsiMaxNumIterationHotStart, 100); // go to end
#ifdef COIN_HAS_CLP
      /* If we are doing all strong branching in one go then we create new arrays
               to store information.  If clp NULL then doing old way.
               Going down -
               outputSolution[2*i] is final solution.
               outputStuff[2*i] is status (0 - finished, 1 infeas, other unknown
               outputStuff[2*i+numberStrong] is number iterations
               On entry newUpper[i] is new upper bound, on exit obj change
               Going up -
               outputSolution[2*i+1] is final solution.
               outputStuff[2*i+1] is status (0 - finished, 1 infeas, other unknown
               outputStuff[2*i+1+numberStrong] is number iterations
            On entry newLower[i] is new lower bound, on exit obj change
            */
      ClpSimplex *clp = NULL;
      double *newLower = NULL;
      double *newUpper = NULL;
      double **outputSolution = NULL;
      int *outputStuff = NULL;
      // Go back to normal way if user wants it
      if (osiclp && (osiclp->specialOptions() & 16) != 0 && osiclp->specialOptions() > 0)
        allNormal = false;
      if (osiclp && !allNormal) {
        // say do fast
        int easy = 1;
        osiclp->setHintParam(OsiDoInBranchAndCut, true, OsiHintDo, &easy);
      }
      if (osiclp && allNormal) {
        clp = osiclp->getModelPtr();
        // Clp - do a different way
        newLower = new double[numberStrong];
        newUpper = new double[numberStrong];
        outputSolution = new double *[2 * numberStrong];
        outputStuff = new int[4 * numberStrong];
        int *which = new int[numberStrong];
        int startFinishOptions;
        int specialOptions = osiclp->specialOptions();
        int clpOptions = clp->specialOptions();
        int returnCode = 0;
#define CRUNCH
#ifdef CRUNCH
        // Crunch down problem
        int numberRows = clp->numberRows();
        // Use dual region
        double *rhs = clp->dualRowSolution();
        int *whichRow = new int[3 * numberRows];
        int *whichColumn = new int[2 * numberColumns];
        int nBound;
        ClpSimplex *small = static_cast< ClpSimplexOther * >(clp)->crunch(rhs, whichRow, whichColumn, nBound, true);
        if (!small) {
          anyAction = -2;
          //printf("XXXX Inf by inspection\n");
          delete[] whichColumn;
          whichColumn = NULL;
          delete[] whichRow;
          whichRow = NULL;
          break;
        } else {
          clp = small;
        }
#else
        int saveLogLevel = clp->logLevel();
        int saveMaxIts = clp->maximumIterations();
#endif
        clp->setLogLevel(0);
        if ((specialOptions & 1) == 0) {
          startFinishOptions = 0;
          clp->setSpecialOptions(clpOptions | (64 | 1024));
        } else {
          startFinishOptions = 1 + 2 + 4;
          //startFinishOptions=1+4; // for moment re-factorize
          if ((specialOptions & 4) == 0)
            clp->setSpecialOptions(clpOptions | (64 | 128 | 512 | 1024 | 4096));
          else
            clp->setSpecialOptions(clpOptions | (64 | 128 | 512 | 1024 | 2048 | 4096));
        }
        // User may want to clean up before strong branching
        if ((clp->specialOptions() & 32) != 0) {
          clp->primal(1);
          if (clp->numberIterations())
            model->messageHandler()->message(CBC_ITERATE_STRONG, *model->messagesPointer())
              << clp->numberIterations()
              << CoinMessageEol;
        }
        clp->setMaximumIterations(saveLimit);
#ifdef CRUNCH
        int *backColumn = whichColumn + numberColumns;
#endif
        for (i = 0; i < numberStrong; i++) {
          int iObject = choice[i].objectNumber;
          const OsiObject *object = model->object(iObject);
          const CbcSimpleInteger *simple = static_cast< const CbcSimpleInteger * >(object);
          int iSequence = simple->columnNumber();
          newLower[i] = ceil(saveSolution[iSequence]);
          newUpper[i] = floor(saveSolution[iSequence]);
#ifdef CRUNCH
          iSequence = backColumn[iSequence];
          assert(iSequence >= 0);
#endif
          which[i] = iSequence;
          outputSolution[2 * i] = new double[numberColumns];
          outputSolution[2 * i + 1] = new double[numberColumns];
        }
        //clp->writeMps("bad");
        returnCode = clp->strongBranching(numberStrong, which,
          newLower, newUpper, outputSolution,
          outputStuff, outputStuff + 2 * numberStrong, !solveAll, false,
          startFinishOptions);
#ifndef CRUNCH
        clp->setSpecialOptions(clpOptions); // restore
        clp->setMaximumIterations(saveMaxIts);
        clp->setLogLevel(saveLogLevel);
#endif
        if (returnCode == -2) {
          // bad factorization!!!
          // Doing normal way
          // Mark hot start
          solver->markHotStart();
          clp = NULL;
        } else {
#ifdef CRUNCH
          // extract solution
          //bool checkSol=true;
          for (i = 0; i < numberStrong; i++) {
            int iObject = choice[i].objectNumber;
            const OsiObject *object = model->object(iObject);
            const CbcSimpleInteger *simple = static_cast< const CbcSimpleInteger * >(object);
            int iSequence = simple->columnNumber();
            which[i] = iSequence;
            double *sol = outputSolution[2 * i];
            double *sol2 = outputSolution[2 * i + 1];
            //bool x=true;
            //bool x2=true;
            for (int iColumn = numberColumns - 1; iColumn >= 0; iColumn--) {
              int jColumn = backColumn[iColumn];
              if (jColumn >= 0) {
                sol[iColumn] = sol[jColumn];
                sol2[iColumn] = sol2[jColumn];
              } else {
                sol[iColumn] = saveSolution[iColumn];
                sol2[iColumn] = saveSolution[iColumn];
              }
            }
          }
#endif
        }
#ifdef CRUNCH
        delete[] whichColumn;
        delete[] whichRow;
        delete small;
#endif
        delete[] which;
      } else {
        // Doing normal way
        // Mark hot start
        solver->markHotStart();
      }
#else /* COIN_HAS_CLP */

      OsiSolverInterface *clp = NULL;
      double **outputSolution = NULL;
      int *outputStuff = NULL;
      double *newLower = NULL;
      double *newUpper = NULL;

      solver->markHotStart();

#endif /* COIN_HAS_CLP */
      /*
              Open a loop to do the strong branching LPs. For each candidate variable,
              solve an LP with the variable forced down, then up. If a direction turns
              out to be infeasible or monotonic (i.e., over the dual objective cutoff),
              force the objective change to be big (1.0e100). If we determine the problem
              is infeasible, or find a monotone variable, escape the loop.

              TODO: The `restore bounds' part might be better encapsulated as an
            unbranch() method. Branching objects more exotic than simple integers
            or cliques might not restrict themselves to variable bounds.

              TODO: Virtuous solvers invalidate the current solution (or give bogus
            results :-) when the bounds are changed out from under them. So we
            need to do all the work associated with finding a new solution before
            restoring the bounds.
            */
      for (i = 0; i < numberStrong; i++) {
        double objectiveChange;
        double newObjectiveValue = 1.0e100;
        // status is 0 finished, 1 infeasible and other
        int iStatus;
        /*
                  Try the down direction first. (Specify the initial branching alternative as
                  down with a call to way(-1). Each subsequent call to branch() performs the
                  specified branch and advances the branch object state to the next branch
                  alternative.)
                */
        if (!clp) {
          choice[i].possibleBranch->way(-1);
          choice[i].possibleBranch->branch();
          bool feasible = true;
          if (checkFeasibility) {
            // check branching did not make infeasible
            int iColumn;
            int numberColumns = solver->getNumCols();
            const double *columnLower = solver->getColLower();
            const double *columnUpper = solver->getColUpper();
            for (iColumn = 0; iColumn < numberColumns; iColumn++) {
              if (columnLower[iColumn] > columnUpper[iColumn] + 1.0e-5)
                feasible = false;
            }
          }
          if (feasible) {
            solver->solveFromHotStart();
            if ((model->moreSpecialOptions2() & 32768) != 0 && solver->isProvenOptimal()) {
              // If any small values re-do
              model->cleanBounds(solver, cleanVariables);
            }
            numberStrongDone++;
            numberStrongIterations += solver->getIterationCount();
            /*
                        We now have an estimate of objective degradation that we can use for strong
                        branching. If we're over the cutoff, the variable is monotone up.
                        If we actually made it to optimality, check for a solution, and if we have
                        a good one, call setBestSolution to process it. Note that this may reduce the
                        cutoff, so we check again to see if we can declare this variable monotone.
                        */
            if (solver->isProvenOptimal())
              iStatus = 0; // optimal
            else if (solver->isIterationLimitReached()
              && !solver->isDualObjectiveLimitReached())
              iStatus = 2; // unknown
            else
              iStatus = 1; // infeasible
            newObjectiveValue = solver->getObjSense() * solver->getObjValue();
            choice[i].numItersDown = solver->getIterationCount();
          } else {
            iStatus = 1; // infeasible
            newObjectiveValue = 1.0e100;
            choice[i].numItersDown = 0;
          }
        } else {
          iStatus = outputStuff[2 * i];
          choice[i].numItersDown = outputStuff[2 * numberStrong + 2 * i];
          numberStrongDone++;
          numberStrongIterations += choice[i].numItersDown;
          newObjectiveValue = objectiveValue + newUpper[i];
          solver->setColSolution(outputSolution[2 * i]);
        }
        objectiveChange = CoinMax(newObjectiveValue - objectiveValue_, 0.0);
        if (!iStatus) {
          choice[i].finishedDown = true;
          if (newObjectiveValue >= model->getCutoff()) {
            objectiveChange = 1.0e100; // say infeasible
            numberStrongInfeasible++;
          } else {
            // See if integer solution
            if (model->feasibleSolution(choice[i].numIntInfeasDown,
                  choice[i].numObjInfeasDown)
              && model->problemFeasibility()->feasible(model, -1) >= 0) {
              model->setBestSolution(CBC_STRONGSOL,
                newObjectiveValue,
                solver->getColSolution());
              // only needed for odd solvers
              newObjectiveValue = solver->getObjSense() * solver->getObjValue();
              objectiveChange = CoinMax(newObjectiveValue - objectiveValue_, 0.0);
              model->setLastHeuristic(NULL);
              model->incrementUsed(solver->getColSolution());
              if (newObjectiveValue >= model->getCutoff()) { //  *new* cutoff
                objectiveChange = 1.0e100;
                numberStrongInfeasible++;
              }
            }
          }
        } else if (iStatus == 1) {
          objectiveChange = 1.0e100;
          numberStrongInfeasible++;
        } else {
          // Can't say much as we did not finish
          choice[i].finishedDown = false;
          numberUnfinished++;
        }
        choice[i].downMovement = objectiveChange;

        // restore bounds
        if (!clp) {
          for (int j = 0; j < numberColumns; j++) {
            if (saveLower[j] != lower[j])
              solver->setColLower(j, saveLower[j]);
            if (saveUpper[j] != upper[j])
              solver->setColUpper(j, saveUpper[j]);
          }
        }
        //printf("Down on %d, status is %d, obj %g its %d cost %g finished %d inf %d infobj %d\n",
        //     choice[i].objectNumber,iStatus,newObjectiveValue,choice[i].numItersDown,
        //     choice[i].downMovement,choice[i].finishedDown,choice[i].numIntInfeasDown,
        //     choice[i].numObjInfeasDown);

        // repeat the whole exercise, forcing the variable up
        if (!clp) {
          bool feasible = true;
          // If odd branching then maybe just one possibility
          if (choice[i].possibleBranch->numberBranchesLeft() > 0) {
            choice[i].possibleBranch->branch();
            if (checkFeasibility) {
              // check branching did not make infeasible
              int iColumn;
              int numberColumns = solver->getNumCols();
              const double *columnLower = solver->getColLower();
              const double *columnUpper = solver->getColUpper();
              for (iColumn = 0; iColumn < numberColumns; iColumn++) {
                if (columnLower[iColumn] > columnUpper[iColumn] + 1.0e-5)
                  feasible = false;
              }
            }
          } else {
            // second branch infeasible
            feasible = false;
          }
          if (feasible) {
            solver->solveFromHotStart();
            if ((model->moreSpecialOptions2() & 32768) != 0 && solver->isProvenOptimal()) {
              // If any small values re-do
              model->cleanBounds(solver, cleanVariables);
            }
            numberStrongDone++;
            numberStrongIterations += solver->getIterationCount();
            /*
                        We now have an estimate of objective degradation that we can use for strong
                        branching. If we're over the cutoff, the variable is monotone up.
                        If we actually made it to optimality, check for a solution, and if we have
                        a good one, call setBestSolution to process it. Note that this may reduce the
                        cutoff, so we check again to see if we can declare this variable monotone.
                        */
            if (solver->isProvenOptimal())
              iStatus = 0; // optimal
            else if (solver->isIterationLimitReached()
              && !solver->isDualObjectiveLimitReached())
              iStatus = 2; // unknown
            else
              iStatus = 1; // infeasible
            newObjectiveValue = solver->getObjSense() * solver->getObjValue();
            choice[i].numItersUp = solver->getIterationCount();
          } else {
            iStatus = 1; // infeasible
            newObjectiveValue = 1.0e100;
            choice[i].numItersDown = 0;
          }
        } else {
          iStatus = outputStuff[2 * i + 1];
          choice[i].numItersUp = outputStuff[2 * numberStrong + 2 * i + 1];
          numberStrongDone++;
          numberStrongIterations += choice[i].numItersUp;
          newObjectiveValue = objectiveValue + newLower[i];
          solver->setColSolution(outputSolution[2 * i + 1]);
        }
        objectiveChange = CoinMax(newObjectiveValue - objectiveValue_, 0.0);
        if (!iStatus) {
          choice[i].finishedUp = true;
          if (newObjectiveValue >= model->getCutoff()) {
            objectiveChange = 1.0e100; // say infeasible
            numberStrongInfeasible++;
          } else {
            // See if integer solution
            if (model->feasibleSolution(choice[i].numIntInfeasUp,
                  choice[i].numObjInfeasUp)
              && model->problemFeasibility()->feasible(model, -1) >= 0) {
              model->setBestSolution(CBC_STRONGSOL,
                newObjectiveValue,
                solver->getColSolution());
              // only needed for odd solvers
              newObjectiveValue = solver->getObjSense() * solver->getObjValue();
              objectiveChange = CoinMax(newObjectiveValue - objectiveValue_, 0.0);
              model->setLastHeuristic(NULL);
              model->incrementUsed(solver->getColSolution());
              if (newObjectiveValue >= model->getCutoff()) { //  *new* cutoff
                objectiveChange = 1.0e100;
                numberStrongInfeasible++;
              }
            }
          }
        } else if (iStatus == 1) {
          objectiveChange = 1.0e100;
          numberStrongInfeasible++;
        } else {
          // Can't say much as we did not finish
          choice[i].finishedUp = false;
          numberUnfinished++;
        }
        choice[i].upMovement = objectiveChange;

        // restore bounds
        if (!clp) {
          for (int j = 0; j < numberColumns; j++) {
            if (saveLower[j] != lower[j])
              solver->setColLower(j, saveLower[j]);
            if (saveUpper[j] != upper[j])
              solver->setColUpper(j, saveUpper[j]);
          }
        }

        //printf("Up on %d, status is %d, obj %g its %d cost %g finished %d inf %d infobj %d\n",
        //     choice[i].objectNumber,iStatus,newObjectiveValue,choice[i].numItersUp,
        //     choice[i].upMovement,choice[i].finishedUp,choice[i].numIntInfeasUp,
        //     choice[i].numObjInfeasUp);

        /*
                  End of evaluation for this candidate variable. Possibilities are:
                  * Both sides below cutoff; this variable is a candidate for branching.
                  * Both sides infeasible or above the objective cutoff: no further action
                  here. Break from the evaluation loop and assume the node will be purged
                  by the caller.
                  * One side below cutoff: Install the branch (i.e., fix the variable). Break
                  from the evaluation loop and assume the node will be reoptimised by the
                  caller.
                */
        // reset
        choice[i].possibleBranch->resetNumberBranchesLeft();
        if (choice[i].upMovement < 1.0e100) {
          if (choice[i].downMovement < 1.0e100) {
            // feasible - no action
          } else {
            // up feasible, down infeasible
            anyAction = -1;
            //printf("Down infeasible for choice %d sequence %d\n",i,
            // model->object(choice[i].objectNumber)->columnNumber());
            if (!solveAll) {
              choice[i].possibleBranch->way(1);
              choice[i].possibleBranch->branch();
              break;
            } else {
              choice[i].fix = 1;
            }
          }
        } else {
          if (choice[i].downMovement < 1.0e100) {
            // down feasible, up infeasible
            anyAction = -1;
            //printf("Up infeasible for choice %d sequence %d\n",i,
            // model->object(choice[i].objectNumber)->columnNumber());
            if (!solveAll) {
              choice[i].possibleBranch->way(-1);
              choice[i].possibleBranch->branch();
              break;
            } else {
              choice[i].fix = -1;
            }
          } else {
            // neither side feasible
            anyAction = -2;
            //printf("Both infeasible for choice %d sequence %d\n",i,
            // model->object(choice[i].objectNumber)->columnNumber());
            break;
          }
        }
        bool hitMaxTime = (model->getCurrentSeconds() > model->getDblParam(CbcModel::CbcMaximumSeconds));
        if (hitMaxTime) {
          numberStrong = i + 1;
          break;
        }
      }
      if (!clp) {
        // Delete the snapshot
        solver->unmarkHotStart();
      } else {
        delete[] newLower;
        delete[] newUpper;
        delete[] outputStuff;
        int i;
        for (i = 0; i < 2 * numberStrong; i++)
          delete[] outputSolution[i];
        delete[] outputSolution;
      }
      solver->setIntParam(OsiMaxNumIterationHotStart, saveLimit);
      // restore basis
      solver->setWarmStart(ws);
      // Unless infeasible we will carry on
      // But we could fix anyway
      if (anyAction == -1 && solveAll) {
        // apply and take off
        for (i = 0; i < numberStrong; i++) {
          if (choice[i].fix) {
            choice[i].possibleBranch->way(choice[i].fix);
            choice[i].possibleBranch->branch();
          }
        }
        bool feasible = true;
        if (checkFeasibility) {
          // check branching did not make infeasible
          int iColumn;
          int numberColumns = solver->getNumCols();
          const double *columnLower = solver->getColLower();
          const double *columnUpper = solver->getColUpper();
          for (iColumn = 0; iColumn < numberColumns; iColumn++) {
            if (columnLower[iColumn] > columnUpper[iColumn] + 1.0e-5)
              feasible = false;
          }
        }
        if (feasible) {
          // can do quick optimality check
          int easy = 2;
          solver->setHintParam(OsiDoInBranchAndCut, true, OsiHintDo, &easy);
          solver->resolve();
          solver->setHintParam(OsiDoInBranchAndCut, true, OsiHintDo, NULL);
          feasible = solver->isProvenOptimal();
        }
        if (feasible) {
          memcpy(saveSolution, solver->getColSolution(), numberColumns * sizeof(double));
          model->reserveCurrentSolution(saveSolution);
          memcpy(saveLower, solver->getColLower(), numberColumns * sizeof(double));
          memcpy(saveUpper, solver->getColUpper(), numberColumns * sizeof(double));
          // Clean up all candidates whih are fixed
          int numberLeft = 0;
          for (i = 0; i < numberStrong; i++) {
            CbcStrongInfo thisChoice = choice[i];
            choice[i].possibleBranch = NULL;
            const OsiObject *object = model->object(thisChoice.objectNumber);
            double infeasibility = object->checkInfeasibility(&usefulInfo);
            if (!infeasibility) {
              // take out
              delete thisChoice.possibleBranch;
            } else {
              choice[numberLeft++] = thisChoice;
            }
          }
          numberStrong = numberLeft;
          for (; i < maximumStrong; i++) {
            delete choice[i].possibleBranch;
            choice[i].possibleBranch = NULL;
          }
          // If all fixed then round again
          if (!numberLeft) {
            finished = false;
            numberStrong = 0;
            saveNumberStrong = 0;
            maximumStrong = 1;
          } else {
            anyAction = 0;
          }
          // If these two uncommented then different action
          anyAction = -1;
          finished = true;
          //printf("some fixed but continuing %d left\n",numberLeft);
        } else {
          anyAction = -2; // say infeasible
        }
      }
      delete ws;
      //int numberNodes = model->getNodeCount();
      // update number of strong iterations etc
      model->incrementStrongInfo(numberStrongDone, numberStrongIterations,
        anyAction == -2 ? 0 : numberStrongInfeasible, anyAction == -2);

      /*
              anyAction >= 0 indicates that strong branching didn't produce any monotone
              variables. Sift through the candidates for the best one.

              QUERY: Setting numberNodes looks to be a distributed noop. numberNodes is
              local to this code block. Perhaps should be numberNodes_ from model?
              Unclear what this calculation is doing.
            */
      if (anyAction >= 0) {

        // get average cost per iteration and assume stopped ones
        // would stop after 50% more iterations at average cost??? !!! ???
        double averageCostPerIteration = 0.0;
        double totalNumberIterations = 1.0;
        int smallestNumberInfeasibilities = COIN_INT_MAX;
        for (i = 0; i < numberStrong; i++) {
          totalNumberIterations += choice[i].numItersDown + choice[i].numItersUp;
          averageCostPerIteration += choice[i].downMovement + choice[i].upMovement;
          smallestNumberInfeasibilities = CoinMin(CoinMin(choice[i].numIntInfeasDown,
                                                    choice[i].numIntInfeasUp),
            smallestNumberInfeasibilities);
        }
        //if (smallestNumberInfeasibilities>=numberIntegerInfeasibilities)
        //numberNodes=1000000; // switch off search for better solution
        averageCostPerIteration /= totalNumberIterations;
        // all feasible - choose best bet

        // New method does all at once so it can be more sophisticated
        // in deciding how to balance actions.
        // But it does need arrays
        double *changeUp = new double[numberStrong];
        int *numberInfeasibilitiesUp = new int[numberStrong];
        double *changeDown = new double[numberStrong];
        int *numberInfeasibilitiesDown = new int[numberStrong];
        CbcBranchingObject **objects = new CbcBranchingObject *[numberStrong];
        for (i = 0; i < numberStrong; i++) {
          int iColumn = choice[i].possibleBranch->variable();
          model->messageHandler()->message(CBC_STRONG, *model->messagesPointer())
            << i << iColumn
            << choice[i].downMovement << choice[i].numIntInfeasDown
            << choice[i].upMovement << choice[i].numIntInfeasUp
            << choice[i].possibleBranch->value()
            << CoinMessageEol;
          changeUp[i] = choice[i].upMovement;
          numberInfeasibilitiesUp[i] = choice[i].numIntInfeasUp;
          changeDown[i] = choice[i].downMovement;
          numberInfeasibilitiesDown[i] = choice[i].numIntInfeasDown;
          objects[i] = choice[i].possibleBranch;
        }
        int whichObject = decision->bestBranch(objects, numberStrong, numberUnsatisfied_,
          changeUp, numberInfeasibilitiesUp,
          changeDown, numberInfeasibilitiesDown,
          objectiveValue_);
        // move branching object and make sure it will not be deleted
        if (whichObject >= 0) {
          branch_ = objects[whichObject];
          if (model->messageHandler()->logLevel() > 3)
            printf("Choosing column %d\n", choice[whichObject].possibleBranch->variable());
          choice[whichObject].possibleBranch = NULL;
        }
        delete[] changeUp;
        delete[] numberInfeasibilitiesUp;
        delete[] changeDown;
        delete[] numberInfeasibilitiesDown;
        delete[] objects;
      }
#ifdef COIN_HAS_CLP
      if (osiclp && !allNormal) {
        // back to normal
        osiclp->setHintParam(OsiDoInBranchAndCut, true, OsiHintDo, NULL);
      }
#endif
    }
    /*
          Simple branching. Probably just one, but we may have got here
          because of an odd branch e.g. a cut
        */
    else {
      // not strong
      // C) create branching object
      branch_ = choice[bestChoice].possibleBranch;
      choice[bestChoice].possibleBranch = NULL;
    }
  }
  // Set guessed solution value
  guessedObjectiveValue_ = objectiveValue_ + estimatedDegradation;
  //printf ("Node %d depth %d unsatisfied %d sum %g obj %g guess %g\n",
  //	    model->getNodeCount(),depth_,numberUnsatisfied_,
  //	    sumInfeasibilities_,objectiveValue_,guessedObjectiveValue_);
  /*
      Cleanup, then we're outta here.
    */
  if (!model->branchingMethod() || dynamicBranchingObject)
    delete decision;

  for (i = 0; i < maximumStrong; i++)
    delete choice[i].possibleBranch;
  delete[] choice;
  delete[] saveLower;
  delete[] saveUpper;

  // restore solution
  solver->setColSolution(saveSolution);
  delete[] saveSolution;
#ifdef COIN_HAS_CLP
  delete[] cleanVariables;
  if (osiclp)
    osiclp->setSpecialOptions(saveClpOptions);
#endif
  return anyAction;
}

/*
  Version for dynamic pseudo costs.

  **** For now just return if anything odd
  later allow even if odd

  The routine scans through the object list of the model looking for objects
  that indicate infeasibility. It tests each object using strong branching
  and selects the one with the least objective degradation.  A corresponding
  branching object is left attached to lastNode.
  This version gives preference in evaluation to variables which
  have not been evaluated many times.  It also uses numberStrong
  to say give up if last few tries have not changed incumbent.
  See Achterberg, Koch and Martin.

  If strong branching is disabled, a candidate object is chosen essentially
  at random (whatever object ends up in pos'n 0 of the candidate array).

  If a branching candidate is found to be monotone, bounds are set to fix the
  variable and the routine immediately returns (the caller is expected to
  reoptimize).

  If a branching candidate is found to result in infeasibility in both
  directions, the routine immediately returns an indication of infeasibility.

  Returns:  0	both branch directions are feasible
  -1	branching variable is monotone
  -2	infeasible
  -3   Use another method

  For now just fix on objective from strong branching.
*/

int CbcNode::chooseDynamicBranch(CbcModel *model, CbcNode *lastNode,
  OsiSolverBranch *& /*branches*/,
  int numberPassesLeft)

{
  if (lastNode)
    depth_ = lastNode->depth_ + 1;
  else
    depth_ = 0;
  // Go to other choose if hot start
  if (model->hotstartSolution() && (((model->moreSpecialOptions() & 1024) == 0) || true))
    return -3;
  delete branch_;
  branch_ = NULL;
  OsiSolverInterface *solver = model->solver();
  //#define CHECK_DEBUGGER_PATH
#ifdef CHECK_DEBUGGER_PATH
  bool onOptimalPath = false;
  if ((model->specialOptions() & 1) != 0) {
    const OsiRowCutDebugger *debugger = solver->getRowCutDebugger();
    if (debugger) {
      onOptimalPath = true;
    }
  }
#endif
  // get information on solver type
  const OsiAuxInfo *auxInfo = solver->getAuxiliaryInfo();
  const OsiBabSolver *auxiliaryInfo = dynamic_cast< const OsiBabSolver * >(auxInfo);
  if (!auxiliaryInfo) {
    // use one from CbcModel
    auxiliaryInfo = model->solverCharacteristics();
  }
  int numberObjects = model->numberObjects();
  // If very odd set of objects then use older chooseBranch
  bool useOldWay = false;
  // point to useful information
  OsiBranchingInformation usefulInfo = model->usefulInformation();
  if (numberObjects > model->numberIntegers()) {
    for (int i = model->numberIntegers(); i < numberObjects; i++) {
      OsiObject *object = model->modifiableObject(i);
      CbcObject *obj = dynamic_cast< CbcObject * >(object);
      if (!obj || !obj->optionalObject()) {
        double infeasibility = object->checkInfeasibility(&usefulInfo);
        if (infeasibility) {
          useOldWay = true;
          break;
        }
      } else {
        obj->initializeForBranching(model);
      }
    }
  }
  if ((model->specialOptions() & 128) != 0)
    useOldWay = false; // allow
  // For now return if not simple
  if (useOldWay)
    return -3;
  // Modify useful info
  usefulInfo.depth_ = depth_;
  if ((model->specialOptions() & 128) != 0) {
    // SOS - shadow prices
    int numberRows = solver->getNumRows();
    const double *pi = usefulInfo.pi_;
    double sumPi = 0.0;
    for (int i = 0; i < numberRows; i++)
      sumPi += fabs(pi[i]);
    sumPi /= static_cast< double >(numberRows);
    // and scale back
    sumPi *= 0.01;
    usefulInfo.defaultDual_ = sumPi; // switch on
    int numberColumns = solver->getNumCols();
    int size = CoinMax(numberColumns, 2 * numberRows);
    usefulInfo.usefulRegion_ = new double[size];
    CoinZeroN(usefulInfo.usefulRegion_, size);
    usefulInfo.indexRegion_ = new int[size];
    // pi may change
    usefulInfo.pi_ = CoinCopyOfArray(usefulInfo.pi_, numberRows);
  }
  assert(auxiliaryInfo);
  double cutoff = model->getCutoff();
  const double *lower = solver->getColLower();
  const double *upper = solver->getColUpper();
  // See if user thinks infeasible
  int anyAction = model->problemFeasibility()->feasible(model, 0);
  if (anyAction) {
    // will return -2 if infeasible , 0 if treat as integer
    return anyAction - 1;
  }
  int i;
  int saveStateOfSearch = model->stateOfSearch() % 10;
  int numberStrong = model->numberStrong();
  /* Ranging is switched off.
       The idea is that you can find out the effect of one iteration
       on each unsatisfied variable cheaply.  Then use this
       if you have not got much else to go on.
    */
  //#define RANGING
#ifdef RANGING
  // must have clp
#ifndef COIN_HAS_CLP
#warning("Ranging switched off as not Clp");
#undef RANGING
#endif
  // Pass number
  int kPass = 0;
  int numberRows = solver->getNumRows();
#endif
  int numberColumns = model->getNumCols();
  double *saveUpper = new double[numberColumns];
  double *saveLower = new double[numberColumns];
  for (i = 0; i < numberColumns; i++) {
    saveLower[i] = lower[i];
    saveUpper[i] = upper[i];
  }

  // Save solution in case heuristics need good solution later

  double *saveSolution = new double[numberColumns];
  memcpy(saveSolution, solver->getColSolution(), numberColumns * sizeof(double));
  model->reserveCurrentSolution(saveSolution);
  const double *hotstartSolution = model->hotstartSolution();
  const int *hotstartPriorities = model->hotstartPriorities();
  double integerTolerance = model->getDblParam(CbcModel::CbcIntegerTolerance);
  if (hotstartSolution) {
    if ((model->moreSpecialOptions() & 1024) != 0 || true) {
      int nBad = 0;
      int nUnsat = 0;
      int nDiff = 0;
      for (int i = 0; i < numberObjects; i++) {
        OsiObject *object = model->modifiableObject(i);
        const CbcSimpleInteger *thisOne = dynamic_cast< const CbcSimpleInteger * >(object);
        if (thisOne) {
          int iColumn = thisOne->columnNumber();
          double targetValue = hotstartSolution[iColumn];
          double value = saveSolution[iColumn];
          if (fabs(value - floor(value + 0.5)) > 1.0e-6) {
            nUnsat++;
#ifdef CLP_INVESTIGATE
            printf("H %d is %g target %g\n", iColumn, value, targetValue);
#endif
          } else if (fabs(targetValue - value) > 1.0e-6) {
            nDiff++;
          }
          if (targetValue < saveLower[iColumn] || targetValue > saveUpper[iColumn]) {
#ifdef CLP_INVESTIGATE
            printf("%d has target %g and current bounds %g and %g\n",
              iColumn, targetValue, saveLower[iColumn], saveUpper[iColumn]);
#endif
            nBad++;
          }
        }
      }
#ifdef CLP_INVESTIGATE
      printf("Hot %d unsatisfied, %d outside limits, %d different\n",
        nUnsat, nBad, nDiff);
#endif
      if (nBad) {
        // switch off as not possible
        hotstartSolution = NULL;
        model->setHotstartSolution(NULL, NULL);
        usefulInfo.hotstartSolution_ = NULL;
      }
    }
  }
  /*
      Get a branching decision object. Use the default dynamic decision criteria unless
      the user has loaded a decision method into the model.
    */
  CbcBranchDecision *decision = model->branchingMethod();
  if (!decision)
    decision = new CbcBranchDynamicDecision();
  int xMark = 0;
  // Get arrays to sort
  double *sort = new double[numberObjects];
  int *whichObject = new int[numberObjects];
#ifdef RANGING
  int xPen = 0;
  int *objectMark = new int[2 * numberObjects + 1];
#endif
  // Arrays with movements
  double *upEstimate = new double[numberObjects];
  double *downEstimate = new double[numberObjects];
  double estimatedDegradation = 0.0;
  int numberNodes = model->getNodeCount();
  int saveLogLevel = model->logLevel();
#ifdef JJF_ZERO
  if ((numberNodes % 500) == 0) {
    model->setLogLevel(6);
    // Get average up and down costs
    double averageUp = 0.0;
    double averageDown = 0.0;
    int numberUp = 0;
    int numberDown = 0;
    int i;
    for (i = 0; i < numberObjects; i++) {
      OsiObject *object = model->modifiableObject(i);
      CbcSimpleIntegerDynamicPseudoCost *dynamicObject = dynamic_cast< CbcSimpleIntegerDynamicPseudoCost * >(object);
      assert(dynamicObject);
      int numberUp2 = 0;
      int numberDown2 = 0;
      double up = 0.0;
      double down = 0.0;
      if (dynamicObject->numberTimesUp()) {
        numberUp++;
        averageUp += dynamicObject->upDynamicPseudoCost();
        numberUp2 += dynamicObject->numberTimesUp();
        up = dynamicObject->upDynamicPseudoCost();
      }
      if (dynamicObject->numberTimesDown()) {
        numberDown++;
        averageDown += dynamicObject->downDynamicPseudoCost();
        numberDown2 += dynamicObject->numberTimesDown();
        down = dynamicObject->downDynamicPseudoCost();
      }
      if (numberUp2 || numberDown2)
        printf("col %d - up %d times cost %g, - down %d times cost %g\n",
          dynamicObject->columnNumber(), numberUp2, up, numberDown2, down);
    }
    if (numberUp)
      averageUp /= static_cast< double >(numberUp);
    else
      averageUp = 1.0;
    if (numberDown)
      averageDown /= static_cast< double >(numberDown);
    else
      averageDown = 1.0;
    printf("total - up %d vars average %g, - down %d vars average %g\n",
      numberUp, averageUp, numberDown, averageDown);
  }
#endif
  int numberBeforeTrust = model->numberBeforeTrust();
  // May go round twice if strong branching fixes all local candidates
  bool finished = false;
  int numberToFix = 0;
  // Mark variables which need to be clean
  char *cleanVariables = NULL;
#ifdef COIN_HAS_CLP
  OsiClpSolverInterface *osiclp = dynamic_cast< OsiClpSolverInterface * >(solver);
  int saveClpOptions = 0;
  if (osiclp) {
    if ((model->moreSpecialOptions2() & 32768) != 0) {
      cleanVariables = model->setupCleanVariables(); // for odd ints/sos etc
    }
    // for faster hot start
    saveClpOptions = osiclp->specialOptions();
    osiclp->setSpecialOptions(saveClpOptions | 8192);
  }
#else
  OsiSolverInterface *osiclp = NULL;
#endif
  //const CglTreeProbingInfo * probingInfo = NULL; //model->probingInfo();
  // Old code left in with DEPRECATED_STRATEGY
  assert(model->searchStrategy() == -1 || model->searchStrategy() == 1 || model->searchStrategy() == 2);
#ifdef DEPRECATED_STRATEGY
  int saveSearchStrategy2 = model->searchStrategy();
#endif
  // Get average up and down costs
  {
    double averageUp = 0.0;
    double averageDown = 0.0;
    int numberUp = 0;
    int numberDown = 0;
    int i;
    for (i = 0; i < numberObjects; i++) {
      OsiObject *object = model->modifiableObject(i);
      CbcSimpleIntegerDynamicPseudoCost *dynamicObject = dynamic_cast< CbcSimpleIntegerDynamicPseudoCost * >(object);
      if (dynamicObject) {
        if (dynamicObject->numberTimesUp()) {
          numberUp++;
          averageUp += dynamicObject->upDynamicPseudoCost();
        }
        if (dynamicObject->numberTimesDown()) {
          numberDown++;
          averageDown += dynamicObject->downDynamicPseudoCost();
        }
      }
    }
    if (numberUp)
      averageUp /= static_cast< double >(numberUp);
    else
      averageUp = 1.0;
    if (numberDown)
      averageDown /= static_cast< double >(numberDown);
    else
      averageDown = 1.0;
    for (i = 0; i < numberObjects; i++) {
      OsiObject *object = model->modifiableObject(i);
      CbcSimpleIntegerDynamicPseudoCost *dynamicObject = dynamic_cast< CbcSimpleIntegerDynamicPseudoCost * >(object);
      if (dynamicObject) {
        if (!dynamicObject->numberTimesUp())
          dynamicObject->setUpDynamicPseudoCost(averageUp);
        if (!dynamicObject->numberTimesDown())
          dynamicObject->setDownDynamicPseudoCost(averageDown);
      }
    }
  }
  /*
      1 strong
      2 no strong
      3 strong just before solution
      4 no strong just before solution
      5 strong first time or before solution
      6 strong first time
    */
  int useShadow = model->moreSpecialOptions() & 7;
  if (useShadow > 2) {
    if (model->getSolutionCount()) {
      if (numberNodes || useShadow < 5) {
        useShadow = 0;
        // zap pseudo shadow prices
        model->pseudoShadow(-1);
        // and switch off
        model->setMoreSpecialOptions(model->moreSpecialOptions() & (~1023));
      } else {
        useShadow = 1;
      }
    } else if (useShadow < 5) {
      useShadow -= 2;
    } else {
      useShadow = 1;
    }
  }
  if (useShadow) {
    // pseudo shadow prices
    model->pseudoShadow((model->moreSpecialOptions() >> 3) & 63);
  }
#ifdef DEPRECATED_STRATEGY
  { // in for tabbing
  }
  else if (saveSearchStrategy2 < 1999)
  {
    // pseudo shadow prices
    model->pseudoShadow(NULL, NULL);
  }
  else if (saveSearchStrategy2 < 2999)
  {
    // leave old ones
  }
  else if (saveSearchStrategy2 < 3999)
  {
    // pseudo shadow prices at root
    if (!numberNodes)
      model->pseudoShadow(NULL, NULL);
  }
  else
  {
    abort();
  }
  if (saveSearchStrategy2 >= 0)
    saveSearchStrategy2 = saveSearchStrategy2 % 1000;
  if (saveSearchStrategy2 == 999)
    saveSearchStrategy2 = -1;
  int saveSearchStrategy = saveSearchStrategy2 < 99 ? saveSearchStrategy2 : saveSearchStrategy2 - 100;
#endif //DEPRECATED_STRATEGY
  int numberNotTrusted = 0;
  int numberStrongDone = 0;
  int numberUnfinished = 0;
  int numberStrongInfeasible = 0;
  int numberStrongIterations = 0;
  int strongType = 0;
#define DO_ALL_AT_ROOT
#ifdef DO_ALL_AT_ROOT
  int saveSatisfiedVariables = 0;
  int saveNumberToDo = 0;
#endif
  // so we can save lots of stuff
  CbcStrongInfo choice;
  memset(&choice, 0, sizeof(CbcStrongInfo));
  CbcDynamicPseudoCostBranchingObject *choiceObject = NULL;
  if (model->allDynamic()) {
    CbcSimpleIntegerDynamicPseudoCost *object = NULL;
    choiceObject = new CbcDynamicPseudoCostBranchingObject(model, 0, -1, 0.5, object);
  }
  choice.possibleBranch = choiceObject;
  numberPassesLeft = CoinMax(numberPassesLeft, 2);
  /* How dogged to be in strong branching
       0 - default
       1 - go to end on first time
       2 - always go to end
     */
  int goToEndInStrongBranching = (model->moreSpecialOptions2() & (3 * 8192)) >> 13;
#ifdef COIN_HAS_NTY
  // 1 after, 2 strong, 3 until depth 5
  int orbitOption = (model->moreSpecialOptions2() & (128 | 256)) >> 7;
#endif
  //#define DEBUG_SOLUTION
#ifdef DEBUG_SOLUTION
  bool onOptimalPath = false;
  if ((model->specialOptions() & 1) != 0) {
    const OsiRowCutDebugger *debugger = model->continuousSolver()->getRowCutDebugger();
    if (debugger) {
      const OsiRowCutDebugger *debugger2 = model->solver()->getRowCutDebugger();
      printf("On optimal in CbcNode %s\n", debugger2 ? "" : "but bad cuts");
      onOptimalPath = true;
    }
  }
#endif
  while (!finished) {
    numberPassesLeft--;
    finished = true;
    decision->initialize(model);
    // Some objects may compute an estimate of best solution from here
    estimatedDegradation = 0.0;
    numberToFix = 0;
    int numberToDo = 0;
    int iBestNot = -1;
    int iBestGot = -1;
    double best = 0.0;
    numberNotTrusted = 0;
    numberStrongDone = 0;
    numberUnfinished = 0;
    numberStrongInfeasible = 0;
    numberStrongIterations = 0;
#ifdef RANGING
    int *which = objectMark + numberObjects + 1;
    int neededPenalties;
    int optionalPenalties;
#endif
    // We may go round this loop three times (only if we think we have solution)
    for (int iPass = 0; iPass < 3; iPass++) {

      // Some objects may compute an estimate of best solution from here
      estimatedDegradation = 0.0;
      numberUnsatisfied_ = 0;
      // initialize sum of "infeasibilities"
      sumInfeasibilities_ = 0.0;
      int bestPriority = COIN_INT_MAX;
#ifdef JJF_ZERO
      int number01 = 0;
      const cliqueEntry *entry = NULL;
      const int *toZero = NULL;
      const int *toOne = NULL;
      const int *backward = NULL;
      int numberUnsatisProbed = 0;
      int numberUnsatisNotProbed = 0; // 0-1
      if (probingInfo) {
        number01 = probingInfo->numberIntegers();
        entry = probingInfo->fixEntries();
        toZero = probingInfo->toZero();
        toOne = probingInfo->toOne();
        backward = probingInfo->backward();
        if (!toZero[number01] || number01 < numberObjects || true) {
          // no info
          probingInfo = NULL;
        }
      }
#endif
      /*
              Scan for branching objects that indicate infeasibility. Choose candidates
              using priority as the first criteria, then integer infeasibility.

              The algorithm is to fill the array with a set of good candidates (by
              infeasibility) with priority bestPriority.  Finding a candidate with
              priority better (less) than bestPriority flushes the choice array. (This
              serves as initialization when the first candidate is found.)

            */
      numberToDo = 0;
#ifdef RANGING
      neededPenalties = 0;
      optionalPenalties = numberObjects;
#endif
      iBestNot = -1;
      double bestNot = 0.0;
      iBestGot = -1;
      best = 0.0;
      /* Problem type as set by user or found by analysis.  This will be extended
            0 - not known
            1 - Set partitioning <=
            2 - Set partitioning ==
            3 - Set covering
            4 - all +- 1 or all +1 and odd
            */
      int problemType = model->problemType();
      bool canDoOneHot = false;
      for (i = 0; i < numberObjects; i++) {
        OsiObject *object = model->modifiableObject(i);
        CbcSimpleIntegerDynamicPseudoCost *dynamicObject = dynamic_cast< CbcSimpleIntegerDynamicPseudoCost * >(object);
        double infeasibility = object->checkInfeasibility(&usefulInfo);
        int priorityLevel = object->priority();
        if (hotstartSolution) {
          // we are doing hot start
          const CbcSimpleInteger *thisOne = dynamic_cast< const CbcSimpleInteger * >(object);
          if (thisOne) {
            int iColumn = thisOne->columnNumber();
            bool canDoThisHot = true;
            double targetValue = hotstartSolution[iColumn];
            if (saveUpper[iColumn] > saveLower[iColumn]) {
              double value = saveSolution[iColumn];
              if (hotstartPriorities)
                priorityLevel = hotstartPriorities[iColumn];
              //double originalLower = thisOne->originalLower();
              //double originalUpper = thisOne->originalUpper();
              // switch off if not possible
              if (targetValue >= saveLower[iColumn] && targetValue <= saveUpper[iColumn]) {
                /* priority outranks rest always if negative
                                   otherwise can be downgraded if at correct level.
                                   Infeasibility may be increased to choose 1.0 values first.
                                   choose one near wanted value
                                */
                if (fabs(value - targetValue) > integerTolerance) {
                  //if (infeasibility>0.01)
                  //infeasibility = fabs(1.0e6-fabs(value-targetValue));
                  //else
                  infeasibility = fabs(value - targetValue);
                  priorityLevel = CoinAbs(priorityLevel);
                } else if (priorityLevel < 0) {
                  priorityLevel = CoinAbs(priorityLevel);
                  if (targetValue == saveLower[iColumn] || targetValue == saveUpper[iColumn]) {
                    infeasibility = integerTolerance + 1.0e-12;
                  } else {
                    // can't
                    priorityLevel += 10000000;
                    canDoThisHot = false;
                  }
                } else {
                  priorityLevel += 10000000;
                  canDoThisHot = false;
                }
              } else {
                // switch off if not possible
                canDoThisHot = false;
              }
              if (canDoThisHot)
                canDoOneHot = true;
            } else if (targetValue < saveLower[iColumn] || targetValue > saveUpper[iColumn]) {
            }
          } else {
            priorityLevel += 10000000;
          }
        }
#define ZERO_ONE 0
#define ZERO_FAKE 1.0e20;
#if ZERO_ONE == 1
        // branch on 0-1 first (temp)
        if (fabs(saveSolution[dynamicObject->columnNumber()]) < 1.0)
          priorityLevel--;
#endif
#if ZERO_ONE == 2
        if (fabs(saveSolution[dynamicObject->columnNumber()]) < 1.0)
          infeasibility *= ZERO_FAKE;
#endif
        if (infeasibility) {
          int iColumn = numberColumns + i;
          bool gotDown = false;
          int numberThisDown = 0;
          bool gotUp = false;
          int numberThisUp = 0;
          double downGuess = object->downEstimate();
          double upGuess = object->upEstimate();
          if (dynamicObject) {
            // Use this object's numberBeforeTrust
            int numberBeforeTrustThis = dynamicObject->numberBeforeTrust();
            iColumn = dynamicObject->columnNumber();
            gotDown = false;
            numberThisDown = dynamicObject->numberTimesDown();
            if (numberThisDown >= numberBeforeTrustThis)
              gotDown = true;
            gotUp = false;
            numberThisUp = dynamicObject->numberTimesUp();
            if (numberThisUp >= numberBeforeTrustThis)
              gotUp = true;
            if (!depth_ && false) {
              // try closest to 0.5
              double part = saveSolution[iColumn] - floor(saveSolution[iColumn]);
              infeasibility = fabs(0.5 - part);
            }
            if (problemType > 0 && problemType < 4 && false) {
              // try closest to 0.5
              double part = saveSolution[iColumn] - floor(saveSolution[iColumn]);
              infeasibility = 0.5 - fabs(0.5 - part);
            }
#ifdef JJF_ZERO
            if (probingInfo) {
              int iSeq = backward[iColumn];
              assert(iSeq >= 0);
              infeasibility = 1.0 + (toZero[iSeq + 1] - toZero[iSeq]) + 5.0 * CoinMin(toOne[iSeq] - toZero[iSeq], toZero[iSeq + 1] - toOne[iSeq]);
              if (toZero[iSeq + 1] > toZero[iSeq]) {
                numberUnsatisProbed++;
              } else {
                numberUnsatisNotProbed++;
              }
            }
#endif
          } else {
            // see if SOS
            CbcSOS *sosObject = dynamic_cast< CbcSOS * >(object);
            if (sosObject) {
              gotDown = false;
              numberThisDown = sosObject->numberTimesDown();
              if (numberThisDown >= numberBeforeTrust)
                gotDown = true;
              gotUp = false;
              numberThisUp = sosObject->numberTimesUp();
              if (numberThisUp >= numberBeforeTrust)
                gotUp = true;
            } else {
              gotDown = true;
              numberThisDown = 999999;
              downGuess = 1.0e20;
              gotUp = true;
              numberThisUp = 999999;
              upGuess = 1.0e20;
              numberPassesLeft = 0;
            }
          }
          // Increase estimated degradation to solution
          estimatedDegradation += CoinMin(downGuess, upGuess);
          downEstimate[i] = downGuess;
          upEstimate[i] = upGuess;
          numberUnsatisfied_++;
          sumInfeasibilities_ += infeasibility;
          // Better priority? Flush choices.
          if (priorityLevel < bestPriority) {
            numberToDo = 0;
            bestPriority = priorityLevel;
            iBestGot = -1;
            best = 0.0;
            numberNotTrusted = 0;
#ifdef RANGING
            neededPenalties = 0;
            optionalPenalties = numberObjects;
#endif
          } else if (priorityLevel > bestPriority) {
            continue;
          }
          if (!gotUp || !gotDown)
            numberNotTrusted++;
          // Check for suitability based on infeasibility.
          if ((gotDown && gotUp) && numberStrong > 0) {
            sort[numberToDo] = -infeasibility;
            if (infeasibility > best) {
              best = infeasibility;
              iBestGot = numberToDo;
            }
#ifdef RANGING
            if (dynamicObject) {
              objectMark[--optionalPenalties] = numberToDo;
              which[optionalPenalties] = iColumn;
            }
#endif
          } else {
#ifdef RANGING
            if (dynamicObject) {
              objectMark[neededPenalties] = numberToDo;
              which[neededPenalties++] = iColumn;
            }
#endif
            sort[numberToDo] = -10.0 * infeasibility;
            if (!(numberThisUp + numberThisDown))
              sort[numberToDo] *= 100.0; // make even more likely
            if (iColumn < numberColumns) {
              double part = saveSolution[iColumn] - floor(saveSolution[iColumn]);
              if (1.0 - fabs(part - 0.5) > bestNot) {
                iBestNot = numberToDo;
                bestNot = 1.0 - fabs(part - 0.5);
              }
            } else {
              // SOS
              if (-sort[numberToDo] > bestNot) {
                iBestNot = numberToDo;
                bestNot = -sort[numberToDo];
              }
            }
          }
          if (model->messageHandler()->logLevel() > 3) {
            printf("%d (%d) down %d %g up %d %g - infeas %g - sort %g solution %g\n",
              i, iColumn, numberThisDown, object->downEstimate(), numberThisUp, object->upEstimate(),
              infeasibility, sort[numberToDo], saveSolution[iColumn]);
          }
          whichObject[numberToDo++] = i;
        } else {
          // for debug
          downEstimate[i] = -1.0;
          upEstimate[i] = -1.0;
        }
      }
      if (numberUnsatisfied_) {
        //if (probingInfo&&false)
        //printf("nunsat %d, %d probed, %d other 0-1\n",numberUnsatisfied_,
        // numberUnsatisProbed,numberUnsatisNotProbed);
        // some infeasibilities - go to next steps
        if (!canDoOneHot && hotstartSolution) {
          // switch off as not possible
          hotstartSolution = NULL;
          model->setHotstartSolution(NULL, NULL);
          usefulInfo.hotstartSolution_ = NULL;
        }
        break;
      } else if (!iPass) {
        // may just need resolve
        model->resolve(NULL, 11, saveSolution, saveLower, saveUpper);
        double newObjValue = solver->getObjSense() * solver->getObjValue();
        objectiveValue_ = CoinMax(objectiveValue_, newObjValue);
        if (!solver->isProvenOptimal()) {
          // infeasible
          anyAction = -2;
          break;
        }
        // Double check looks OK - just look at rows with all integers
        if (model->allDynamic()) {
          double *solution = CoinCopyOfArray(saveSolution, numberColumns);
          for (int i = 0; i < numberColumns; i++) {
            if (model->isInteger(i))
              solution[i] = floor(solution[i] + 0.5);
          }
          int numberRows = solver->getNumRows();
          double *rowActivity = new double[numberRows];
          CoinZeroN(rowActivity, numberRows);
          solver->getMatrixByCol()->times(solution, rowActivity);
          //const double * element = model->solver()->getMatrixByCol()->getElements();
          const int *row = model->solver()->getMatrixByCol()->getIndices();
          const CoinBigIndex *columnStart = model->solver()->getMatrixByCol()->getVectorStarts();
          const int *columnLength = model->solver()->getMatrixByCol()->getVectorLengths();
          int nFree = 0;
          int nFreeNon = 0;
          int nFixedNon = 0;
          double mostAway = 0.0;
          int whichAway = -1;
          const double *columnLower = solver->getColLower();
          const double *columnUpper = solver->getColUpper();
          for (int i = 0; i < numberColumns; i++) {
            if (!model->isInteger(i)) {
              // mark rows as flexible
              CoinBigIndex start = columnStart[i];
              CoinBigIndex end = start + columnLength[i];
              for (CoinBigIndex j = start; j < end; j++) {
                int iRow = row[j];
                rowActivity[iRow] = COIN_DBL_MAX;
              }
            } else if (columnLower[i] < columnUpper[i]) {
              double solutionValue = saveSolution[i];
              if (fabs(solution[i] - solutionValue) > integerTolerance && (solutionValue - columnLower[i]) > integerTolerance && (columnUpper[i] - solutionValue) > integerTolerance) {
                nFreeNon++;
                if (fabs(solution[i] - saveSolution[i]) > mostAway) {
                  mostAway = fabs(solution[i] - saveSolution[i]);
                  whichAway = i;
                }
              } else {
                nFree++;
              }
            } else if (solution[i] != saveSolution[i]) {
              nFixedNon++;
            }
          }
          const double *lower = solver->getRowLower();
          const double *upper = solver->getRowUpper();
          bool satisfied = true;
          for (int i = 0; i < numberRows; i++) {
            double value = rowActivity[i];
            if (value != COIN_DBL_MAX) {
              if (value > upper[i] + 1.0e-5 || value < lower[i] - 1.0e-5) {
                satisfied = false;
              }
            }
          }
          delete[] rowActivity;
          if (!satisfied) {
#ifdef CLP_INVESTIGATE
            printf("%d free ok %d free off target %d fixed off target\n",
              nFree, nFreeNon, nFixedNon);
#endif
            if (nFreeNon) {
              // try branching on these
              delete branch_;
              numberUnsatisfied_ = nFreeNon;
              for (int i = 0; i < numberObjects; i++) {
                OsiObject *object = model->modifiableObject(i);
                CbcSimpleIntegerDynamicPseudoCost *obj = dynamic_cast< CbcSimpleIntegerDynamicPseudoCost * >(object);
                assert(obj);
                int iColumn = obj->columnNumber();
                if (iColumn == whichAway) {
                  int preferredWay = (saveSolution[iColumn] > solution[iColumn])
                    ? -1
                    : +1;
                  usefulInfo.integerTolerance_ = 0.0;
                  branch_ = obj->createCbcBranch(solver, &usefulInfo, preferredWay);
                  iBestGot = i;
                  break;
                }
              }
              anyAction = 0;
              delete[] solution;
              break;
            }
          }
          delete[] solution;
        }
      } else if (iPass == 1) {
        // looks like a solution - get paranoid
        bool roundAgain = false;
        // get basis
        CoinWarmStartBasis *ws = dynamic_cast< CoinWarmStartBasis * >(solver->getWarmStart());
        if (!ws)
          break;
        double tolerance;
        solver->getDblParam(OsiPrimalTolerance, tolerance);
        for (i = 0; i < numberColumns; i++) {
          double value = saveSolution[i];
          if (value < lower[i] - tolerance) {
            saveSolution[i] = lower[i];
            roundAgain = true;
            ws->setStructStatus(i, CoinWarmStartBasis::atLowerBound);
          } else if (value > upper[i] + tolerance) {
            saveSolution[i] = upper[i];
            roundAgain = true;
            ws->setStructStatus(i, CoinWarmStartBasis::atUpperBound);
          }
        }
        if (roundAgain) {
          // restore basis
          solver->setWarmStart(ws);
          solver->setColSolution(saveSolution);
          delete ws;
          bool takeHint;
          OsiHintStrength strength;
          solver->getHintParam(OsiDoDualInResolve, takeHint, strength);
          solver->setHintParam(OsiDoDualInResolve, false, OsiHintDo);
          model->resolve(NULL, 11, saveSolution, saveLower, saveUpper);
          double newObjValue = solver->getObjSense() * solver->getObjValue();
          objectiveValue_ = CoinMax(objectiveValue_, newObjValue);
          solver->setHintParam(OsiDoDualInResolve, takeHint, strength);
          if (!solver->isProvenOptimal()) {
            // infeasible
            anyAction = -2;
            break;
          }
        } else {
          delete ws;
          break;
        }
      }
    }
    if (anyAction == -2) {
      break;
    }
    // skip if solution
    if (!numberUnsatisfied_)
      break;
    int skipAll = (numberNotTrusted == 0 || numberToDo == 1) ? 1 : 0;
    bool doneHotStart = false;
    //DEPRECATED_STRATEGYint searchStrategy = saveSearchStrategy>=0 ? (saveSearchStrategy%10) : -1;
    int searchStrategy = model->searchStrategy();
    // But adjust depending on ratio of iterations
    if (searchStrategy > 0) {
      if (numberBeforeTrust >= /*5*/ 10 && numberBeforeTrust <= 10) {
        if (searchStrategy != 2) {
          assert(searchStrategy == 1);
          if (depth_ > 5) {
            int numberIterations = model->getIterationCount();
            int numberStrongIterations = model->numberStrongIterations();
            if (numberStrongIterations > numberIterations + 10000) {
              searchStrategy = 2;
              skipAll = 1;
            } else if (numberStrongIterations * 4 + 1000 < numberIterations) {
              searchStrategy = 3;
              skipAll = 0;
            }
          } else {
            searchStrategy = 3;
            skipAll = 0;
          }
        }
      }
    }
    // worth trying if too many times
    // Save basis
    CoinWarmStart *ws = NULL;
    // save limit
    int saveLimit = 0;
    solver->getIntParam(OsiMaxNumIterationHotStart, saveLimit);
    if (!numberPassesLeft)
      skipAll = 1;
    if (!skipAll) {
      ws = solver->getWarmStart();
      int limit = 100;
      if (!saveStateOfSearch && saveLimit < limit && saveLimit == 100)
        solver->setIntParam(OsiMaxNumIterationHotStart, limit);
    }
    // Say which one will be best
    int whichChoice = 0;
    int bestChoice;
    if (iBestGot >= 0)
      bestChoice = iBestGot;
    else
      bestChoice = iBestNot;
    assert(bestChoice >= 0);
    // If we have hit max time don't do strong branching
    bool hitMaxTime = (model->getCurrentSeconds() > model->getDblParam(CbcModel::CbcMaximumSeconds));
    // also give up if we are looping round too much
    if (hitMaxTime || numberPassesLeft <= 0 || useShadow == 2) {
      int iObject = whichObject[bestChoice];
      OsiObject *object = model->modifiableObject(iObject);
      int preferredWay;
      object->infeasibility(&usefulInfo, preferredWay);
      CbcObject *obj = dynamic_cast< CbcObject * >(object);
      assert(obj);
      branch_ = obj->createCbcBranch(solver, &usefulInfo, preferredWay);
      {
        CbcBranchingObject *branchObj = dynamic_cast< CbcBranchingObject * >(branch_);
        assert(branchObj);
        branchObj->way(preferredWay);
      }
      delete ws;
      ws = NULL;
      break;
    } else {
      // say do fast
      int easy = 1;
      if (!skipAll)
        solver->setHintParam(OsiDoInBranchAndCut, true, OsiHintDo, &easy);
      int iDo;
#define RESET_BOUNDS
#ifdef RANGING
      bool useRanging = model->allDynamic() && !skipAll;
      if (useRanging) {
        double currentObjective = solver->getObjValue() * solver->getObjSense();
        double gap = cutoff - currentObjective;
        // relax a bit
        gap *= 1.0000001;
        gap = CoinMax(1.0e-5, gap);
        // off penalties if too much
        double needed = neededPenalties;
        needed *= numberRows;
        if (numberNodes) {
          if (needed > 1.0e6) {
            neededPenalties = 0;
          } else if (gap < 1.0e5) {
            // maybe allow some not needed
            int extra = static_cast< int >((1.0e6 - needed) / numberRows);
            int nStored = numberObjects - optionalPenalties;
            extra = CoinMin(extra, nStored);
            for (int i = 0; i < extra; i++) {
              objectMark[neededPenalties] = objectMark[optionalPenalties + i];
              which[neededPenalties++] = which[optionalPenalties + i];
              ;
            }
          }
        }
        if (osiclp && neededPenalties) {
          assert(!doneHotStart);
          xPen += neededPenalties;
          which--;
          which[0] = neededPenalties;
          osiclp->passInRanges(which);
          // Mark hot start and get ranges
          if (kPass) {
            // until can work out why solution can go funny
            int save = osiclp->specialOptions();
            osiclp->setSpecialOptions(save | 256);
            solver->markHotStart();
#ifdef RESET_BOUNDS
            memcpy(saveLower, solver->getColLower(), solver->getNumCols() * sizeof(double));
            memcpy(saveUpper, solver->getColUpper(), solver->getNumCols() * sizeof(double));
#endif
            osiclp->setSpecialOptions(save);
          } else {
            solver->markHotStart();
#ifdef RESET_BOUNDS
            memcpy(saveLower, solver->getColLower(), solver->getNumCols() * sizeof(double));
            memcpy(saveUpper, solver->getColUpper(), solver->getNumCols() * sizeof(double));
#endif
          }
          doneHotStart = true;
          xMark++;
          kPass++;
          osiclp->passInRanges(NULL);
          const double *downCost = osiclp->upRange();
          const double *upCost = osiclp->downRange();
          bool problemFeasible = true;
          int numberFixed = 0;
          for (int i = 0; i < neededPenalties; i++) {
            int j = objectMark[i];
            int iObject = whichObject[j];
            OsiObject *object = model->modifiableObject(iObject);
            CbcSimpleIntegerDynamicPseudoCost *dynamicObject = dynamic_cast< CbcSimpleIntegerDynamicPseudoCost * >(object);
            // Use this object's numberBeforeTrust
            int numberBeforeTrustThis = dynamicObject->numberBeforeTrust();
            int iSequence = dynamicObject->columnNumber();
            double value = saveSolution[iSequence];
            value -= floor(value);
            double upPenalty = CoinMin(upCost[i], 1.0e110) * (1.0 - value);
            double downPenalty = CoinMin(downCost[i], 1.0e110) * value;
            int numberThisDown = dynamicObject->numberTimesDown();
            int numberThisUp = dynamicObject->numberTimesUp();
            if (!numberBeforeTrustThis) {
              // override
              downEstimate[iObject] = downPenalty;
              upEstimate[iObject] = upPenalty;
              double min1 = CoinMin(downEstimate[iObject],
                upEstimate[iObject]);
              double max1 = CoinMax(downEstimate[iObject],
                upEstimate[iObject]);
              min1 = 0.8 * min1 + 0.2 * max1;
              sort[j] = -min1;
            } else if (numberThisDown < numberBeforeTrustThis || numberThisUp < numberBeforeTrustThis) {
              double invTrust = 1.0 / static_cast< double >(numberBeforeTrustThis);
              if (numberThisDown < numberBeforeTrustThis) {
                double fraction = numberThisDown * invTrust;
                downEstimate[iObject] = fraction * downEstimate[iObject] + (1.0 - fraction) * downPenalty;
              }
              if (numberThisUp < numberBeforeTrustThis) {
                double fraction = numberThisUp * invTrust;
                upEstimate[iObject] = fraction * upEstimate[iObject] + (1.0 - fraction) * upPenalty;
              }
              double min1 = CoinMin(downEstimate[iObject],
                upEstimate[iObject]);
              double max1 = CoinMax(downEstimate[iObject],
                upEstimate[iObject]);
              min1 = 0.8 * min1 + 0.2 * max1;
              min1 *= 10.0;
              if (!(numberThisDown + numberThisUp))
                min1 *= 100.0;
              sort[j] = -min1;
            }
            // seems unreliable
            if (false && CoinMax(downPenalty, upPenalty) > gap) {
              COIN_DETAIL_PRINT(printf("gap %g object %d has down range %g, up %g\n",
                gap, i, downPenalty, upPenalty));
              //sort[j] -= 1.0e50; // make more likely to be chosen
              int number;
              if (downPenalty > gap) {
                number = dynamicObject->numberTimesDown();
                if (upPenalty > gap)
                  problemFeasible = false;
                CbcBranchingObject *branch = dynamicObject->createCbcBranch(solver, &usefulInfo, 1);
                //branch->fix(solver,saveLower,saveUpper,1);
                delete branch;
              } else {
                number = dynamicObject->numberTimesUp();
                CbcBranchingObject *branch = dynamicObject->createCbcBranch(solver, &usefulInfo, 1);
                //branch->fix(solver,saveLower,saveUpper,-1);
                delete branch;
              }
              if (number >= numberBeforeTrustThis)
                dynamicObject->setNumberBeforeTrust(CoinMin(number + 1, 5 * numberBeforeTrust));
              numberFixed++;
            }
            if (!numberNodes)
              COIN_DETAIL_PRINT(printf("%d pen down ps %g -> %g up ps %g -> %g\n",
                iObject, downPenalty, downPenalty, upPenalty, upPenalty));
          }
          if (numberFixed && problemFeasible) {
            assert(doneHotStart);
            solver->unmarkHotStart();
            model->resolve(NULL, 11, saveSolution, saveLower, saveUpper);
#ifdef CHECK_DEBUGGER_PATH
            if ((model->specialOptions() & 1) != 0 && onOptimalPath) {
              const OsiRowCutDebugger *debugger = solver->getRowCutDebugger();
              if (!debugger) {
                printf("Strong branching down on %d went off optimal path\n", iObject);
                abort();
              }
            }
#endif
            double newObjValue = solver->getObjSense() * solver->getObjValue();
            objectiveValue_ = CoinMax(objectiveValue_, newObjValue);
            solver->markHotStart();
#ifdef RESET_BOUNDS
            memcpy(saveLower, solver->getColLower(), solver->getNumCols() * sizeof(double));
            memcpy(saveUpper, solver->getColUpper(), solver->getNumCols() * sizeof(double));
#endif
            problemFeasible = solver->isProvenOptimal();
          }
          if (!problemFeasible) {
            COIN_DETAIL_PRINT(fprintf(stdout, "both ways infeas on ranging - code needed\n"));
            anyAction = -2;
            if (!choiceObject) {
              delete choice.possibleBranch;
              choice.possibleBranch = NULL;
            }
            //printf("Both infeasible for choice %d sequence %d\n",i,
            // model->object(choice.objectNumber)->columnNumber());
            // Delete the snapshot
            solver->unmarkHotStart();
            // back to normal
            solver->setHintParam(OsiDoInBranchAndCut, true, OsiHintDo, NULL);
            // restore basis
            solver->setWarmStart(ws);
            doneHotStart = false;
            delete ws;
            ws = NULL;
            break;
          }
        }
      }
#endif /* RANGING */
      {
        int numberIterations = model->getIterationCount();
        //numberIterations += (model->numberExtraIterations()>>2);
        const int *strongInfo = model->strongInfo();
        //int numberDone = strongInfo[0]-strongInfo[3];
        int numberFixed = strongInfo[1] - strongInfo[4];
        int numberInfeasible = strongInfo[2] - strongInfo[5];
        assert(!strongInfo[3]);
        assert(!strongInfo[4]);
        assert(!strongInfo[5]);
        int numberStrongIterations = model->numberStrongIterations();
        int numberRows = solver->getNumRows();
        if (numberStrongIterations > numberIterations + CoinMin(100, 10 * numberRows) && depth_ >= 4 && numberNodes > 100) {
          if (20 * numberInfeasible + 4 * numberFixed < numberNodes) {
            // Say never do
            if (numberBeforeTrust == 10)
              skipAll = -1;
          }
        }
      }
      // make sure best will be first
      if (iBestGot >= 0)
        sort[iBestGot] = -COIN_DBL_MAX;
      // Actions 0 - exit for repeat, 1 resolve and try old choice,2 exit for continue
      if (anyAction)
        numberToDo = 0; // skip as we will be trying again
      // Sort
      CoinSort_2(sort, sort + numberToDo, whichObject);
      // Change in objective opposite infeasible
      double worstFeasible = 0.0;
      // Just first if strong off
      if (!numberStrong)
        numberToDo = CoinMin(numberToDo, 1);
      if (searchStrategy == 2)
        numberToDo = CoinMin(numberToDo, 20);
      iDo = 0;
      int saveLimit2;
      solver->getIntParam(OsiMaxNumIterationHotStart, saveLimit2);
      int numberTest = numberNotTrusted > 0 ? numberStrong : (numberStrong + 1) / 2;
      if (searchStrategy == 3) {
        // Previously decided we need strong
        numberTest = numberStrong;
      }
      // Try nearly always off
      if (skipAll >= 0) {
        if (searchStrategy < 2) {
          //if ((numberNodes%20)!=0) {
          if ((model->specialOptions() & 8) == 0) {
            numberTest = 0;
          }
          //} else {
          //numberTest=2*numberStrong;
          //skipAll=0;
          //}
        }
      } else {
        // Just take first
        numberTest = 1;
      }
      int testDepth = (skipAll >= 0) ? 8 : 4;
      if (depth_ < testDepth && numberStrong) {
        if (searchStrategy != 2) {
          int numberRows = solver->getNumRows();
          // whether to do this or not is important - think
          if (numberRows < 300 || numberRows + numberColumns < 2500) {
            if (depth_ < 7)
              numberStrong = CoinMin(3 * numberStrong, numberToDo);
            if (!depth_)
              numberStrong = CoinMin(6 * numberStrong, numberToDo);
          }
          numberTest = numberStrong;
          skipAll = 0;
        }
      }
      // Do at least 5 strong
      if (numberColumns < 1000 && (depth_ < 15 || numberNodes < 1000000))
        numberTest = CoinMax(numberTest, 5);
      if ((model->specialOptions() & 8) == 0) {
        if (skipAll) {
          numberTest = 0;
        }
      } else {
        // do 5 as strong is fixing
        numberTest = CoinMax(numberTest, 5);
      }
      // see if switched off
      if (skipAll < 0) {
        numberTest = 0;
      }
      int realMaxHotIterations = 999999;
      if (skipAll < 0)
        numberToDo = 1;
      strongType = 0;
#ifdef DO_ALL_AT_ROOT
      if (model->strongStrategy()) {
        int iStrategy = model->strongStrategy();
        int kDepth = iStrategy / 100;
        if (kDepth)
          iStrategy -= 100 * kDepth;
        else
          kDepth = 5;
        double objValue = solver->getObjSense() * solver->getObjValue();
        double bestPossible = model->getBestPossibleObjValue();
        bestPossible += 1.0e-7 * (1.0 + fabs(bestPossible));
        int jStrategy = iStrategy / 10;
        if (jStrategy) {
          if ((jStrategy & 1) != 0 && !depth_)
            strongType = 2;
          else if ((jStrategy & 2) != 0 && depth_ <= kDepth)
            strongType = 2;
          else if ((jStrategy & 4) != 0 && objValue < bestPossible)
            strongType = 2;
          iStrategy -= 10 * jStrategy;
        }
        if (!strongType) {
          if ((iStrategy & 1) != 0 && !depth_)
            strongType = 1;
          else if ((iStrategy & 2) != 0 && depth_ <= kDepth)
            strongType = 1;
          else if ((iStrategy & 4) != 0 && objValue < bestPossible)
            strongType = 1;
        }
        saveNumberToDo = numberToDo;
        if (strongType == 2) {
          // add in satisfied
          const int *integerVariable = model->integerVariable();
          int numberIntegers = model->numberIntegers();
          if (numberIntegers == numberObjects) {
            numberToDo = 0;
            for (int i = 0; i < numberIntegers; i++) {
              int iColumn = integerVariable[i];
              if (saveUpper[iColumn] > saveLower[iColumn]) {
                whichObject[numberToDo++] = i;
              }
            }
            saveSatisfiedVariables = numberToDo - saveNumberToDo;
          } else {
            strongType = 1;
          }
        }
        if (strongType) {
          numberTest = numberToDo;
          numberStrong = numberToDo;
          skipAll = 0;
          searchStrategy = 0;
          solver->setIntParam(OsiMaxNumIterationHotStart, 100000);
          //printf("Strong branching type %d\n",strongType);
        }
      }
#endif
#ifdef COIN_HAS_NTY
      const int *orbits = NULL;
#endif
#ifdef COIN_HAS_NTY
      if (orbitOption == 2 /* was >1*/) {
        CbcSymmetry *symmetryInfo = model->symmetryInfo();
        CbcNodeInfo *infoX = lastNode ? lastNode->nodeInfo() : NULL;
        bool worthTrying = false;
        if (infoX) {
          CbcNodeInfo *info = infoX;
          for (int i = 0; i < NTY_BAD_DEPTH; i++) {
            if (!info->parent()) {
              worthTrying = true;
              break;
            }
            info = info->parent();
            if (info->symmetryWorked()) {
              worthTrying = true;
              break;
            }
          }
        } else {
          worthTrying = true;
        }
        if (symmetryInfo && worthTrying) {
          symmetryInfo->ChangeBounds(solver->getColLower(),
            solver->getColUpper(),
            solver->getNumCols(), false);
          symmetryInfo->Compute_Symmetry();
          symmetryInfo->fillOrbits();
          orbits = symmetryInfo->whichOrbit();
          int iColumn = -1;
          if (orbits && symmetryInfo->numberUsefulOrbits()) {
            bool doBranch = true;
            int numberUsefulOrbits = symmetryInfo->numberUsefulOrbits();
            if (numberUsefulOrbits < 2) {
              assert(numberUsefulOrbits);
              double largest = -1.0;
              for (int i = 0; i < numberColumns; i++) {
                if (orbits[i] >= 0) {
                  if (saveSolution[i] > largest) {
                    largest = saveSolution[i];
                    iColumn = i;
                  }
                }
              }
            } else {
#if COIN_HAS_NTY2 == 1
              // take largest
              int iOrbit = symmetryInfo->largestOrbit(solver->getColLower(),
                solver->getColUpper());
              double largest = -1.0;
              for (int i = 0; i < numberColumns; i++) {
                if (orbits[i] == iOrbit) {
                  if (saveSolution[i] > largest) {
                    largest = saveSolution[i];
                    iColumn = i;
                  }
                }
              }
#endif
              if (orbitOption == 2) {
                // strong
                int nDo = 0;
                const double *lower = solver->getColLower();
                const double *upper = solver->getColUpper();
                const int *integerVariable = model->integerVariable();
                for (int iOrbit = 0; iOrbit < numberUsefulOrbits; iOrbit++) {
                  double distance = 1.0;
                  int iColumn = -1;
                  int numberIntegers = model->numberIntegers();
                  for (int j = 0; j < numberIntegers; j++) {
                    int i = integerVariable[j];
                    if (orbits[i] == iOrbit && lower[i] == 0.0 && upper[i] == 1.0) {
                      double away = fabs(saveSolution[i] - 0.5);
                      if (away < distance && away < 0.4999) {
                        distance = away;
                        iColumn = j;
                      }
                    }
                  }
                  if (iColumn >= 0)
                    whichObject[nDo++] = iColumn;
                }
                if (nDo)
                  numberToDo = nDo;
                doBranch = false;
              } else if (orbitOption == 3) {
                // subset
                int nDo = 0;
                for (int iDo = 0; iDo < numberToDo; iDo++) {
                  int iObject = whichObject[iDo];
                  OsiObject *object = model->modifiableObject(iObject);
                  CbcSimpleIntegerDynamicPseudoCost *dynamicObject = dynamic_cast< CbcSimpleIntegerDynamicPseudoCost * >(object);
                  int iColumn = dynamicObject ? dynamicObject->columnNumber() : -1;
                  if (iColumn < 0 || orbits[iColumn] >= 0)
                    whichObject[nDo++] = whichObject[iDo];
                }
                assert(nDo);
                //printf("nDo %d\n",nDo);
                numberToDo = nDo;
                doBranch = false;
                /* need NULL as if two in same orbit and strong branching fixes
			 then we may be in trouble.
			 Strong option should be OK as only one in set done.
		       */
                orbits = NULL;
              }
            }
            if (doBranch) {
              orbitOption = 0;
              branch_ = new CbcOrbitalBranchingObject(model, iColumn, 1, 0, NULL);
              if (infoX)
                infoX->setSymmetryWorked();
              numberToDo = 0;
            }
          }
        }
      }
#endif
      for (iDo = 0; iDo < numberToDo; iDo++) {
        int iObject = whichObject[iDo];
        OsiObject *object = model->modifiableObject(iObject);
        CbcSimpleIntegerDynamicPseudoCost *dynamicObject = dynamic_cast< CbcSimpleIntegerDynamicPseudoCost * >(object);
        int iColumn = dynamicObject ? dynamicObject->columnNumber() : numberColumns + iObject;
        int preferredWay;
        double infeasibility = object->infeasibility(&usefulInfo, preferredWay);
        bool feasibleSolution = false;
        double predictedChange = 0.0;
        // may have become feasible
        if (!infeasibility) {
          if (strongType != 2 || solver->getColLower()[iColumn] == solver->getColUpper()[iColumn])
            continue;
        }
#ifndef NDEBUG
        if (iColumn < numberColumns) {
          const double *solution = model->testSolution();
          assert(saveSolution[iColumn] == solution[iColumn]);
        }
#endif
        CbcSimpleInteger *obj = dynamic_cast< CbcSimpleInteger * >(object);
        if (obj) {
          if (choiceObject) {
            obj->fillCreateBranch(choiceObject, &usefulInfo, preferredWay);
            choiceObject->setObject(dynamicObject);
          } else {
            choice.possibleBranch = obj->createCbcBranch(solver, &usefulInfo, preferredWay);
          }
        } else {
          CbcObject *obj = dynamic_cast< CbcObject * >(object);
          assert(obj);
          choice.possibleBranch = obj->createCbcBranch(solver, &usefulInfo, preferredWay);
        }
        // Save which object it was
        choice.objectNumber = iObject;
        choice.numIntInfeasUp = numberUnsatisfied_;
        choice.numIntInfeasDown = numberUnsatisfied_;
        if (strongType != 2) {
          choice.upMovement = upEstimate[iObject];
          choice.downMovement = downEstimate[iObject];
        } else {
          choice.upMovement = 0.1;
          choice.downMovement = 0.1;
        }
        assert(choice.upMovement >= 0.0);
        assert(choice.downMovement >= 0.0);
        choice.fix = 0; // say not fixed
        // see if can skip strong branching
        int canSkip = choice.possibleBranch->fillStrongInfo(choice);
        if ((numberTest <= 0 || skipAll)) {
          if (iDo > 20) {
            if (!choiceObject) {
              delete choice.possibleBranch;
              choice.possibleBranch = NULL;
            }
            break; // give up anyway
          }
        }
        if (model->messageHandler()->logLevel() > 3 && numberBeforeTrust && dynamicObject)
          dynamicObject->print(1, choice.possibleBranch->value());
        if (strongType)
          canSkip = 0;
        if (skipAll < 0)
          canSkip = 1;
        if (!canSkip) {
          if (!doneHotStart) {
            // Mark hot start
            doneHotStart = true;
            solver->markHotStart();
#ifdef RESET_BOUNDS
            memcpy(saveLower, solver->getColLower(), solver->getNumCols() * sizeof(double));
            memcpy(saveUpper, solver->getColUpper(), solver->getNumCols() * sizeof(double));
#endif
            if (!solver->isProvenOptimal()) {
              skipAll = -2;
              canSkip = 1;
            }
            xMark++;
          }
        }
        if (!canSkip) {
          numberTest--;
          // just do a few
          if (searchStrategy == 2)
            solver->setIntParam(OsiMaxNumIterationHotStart, 10);
          double objectiveChange;
          double newObjectiveValue = 1.0e100;
          int j;
#ifdef COIN_HAS_CLP
          int saveMaxHotIts = 0;
          int saveOsiClpOptions = 0;
          if (osiclp && goToEndInStrongBranching) {
            /* How dogged to be in strong branching
			 0 - default
			 1 - go to end on first time
			 2 - always go to end
		      */
            osiclp->getIntParam(OsiMaxNumIterationHotStart, saveMaxHotIts);
            saveOsiClpOptions = osiclp->specialOptions();
            if (goToEndInStrongBranching == 2 || (dynamicObject && dynamicObject->numberTimesBranched() == 0)) {
              if (osiclp->getNumRows() < 200 || goToEndInStrongBranching == 2) {
                osiclp->setIntParam(OsiMaxNumIterationHotStart,
                  10 * (osiclp->getNumRows() + numberColumns));
                osiclp->setSpecialOptions(saveOsiClpOptions & (~32));
              }
            }
          }
#endif
          // status is 0 finished, 1 infeasible and other
          int iStatus;
          /*
                      Try the down direction first. (Specify the initial branching alternative as
                      down with a call to way(-1). Each subsequent call to branch() performs the
                      specified branch and advances the branch object state to the next branch
                      alternative.)
                    */
          choice.possibleBranch->way(-1);
          predictedChange = choice.possibleBranch->branch();
#ifdef COIN_HAS_NTY
          if (orbits) {
            // can fix all in orbit
            int fixOrbit = orbits[iObject];
            if (fixOrbit >= 0) {
              //printf("fixing all in orbit %d for column %d\n",fixOrbit,iObject);
              for (int i = 0; i < numberColumns; i++) {
                if (orbits[i] == fixOrbit)
                  solver->setColUpper(i, 0.0);
              }
            }
          }
#endif
          solver->solveFromHotStart();
          if ((model->moreSpecialOptions2() & 32768) != 0 && solver->isProvenOptimal()) {
            // If any small values re-do
            model->cleanBounds(solver, cleanVariables);
          }
          bool needHotStartUpdate = false;
          numberStrongDone++;
          numberStrongIterations += solver->getIterationCount();
          /*
                      We now have an estimate of objective degradation that we can use for strong
                      branching. If we're over the cutoff, the variable is monotone up.
                      If we actually made it to optimality, check for a solution, and if we have
                      a good one, call setBestSolution to process it. Note that this may reduce the
                      cutoff, so we check again to see if we can declare this variable monotone.
                    */
          if (solver->isProvenOptimal())
            iStatus = 0; // optimal
          else if (solver->isIterationLimitReached()
            && !solver->isDualObjectiveLimitReached()) {
            iStatus = 2; // unknown
          } else {
            iStatus = 1; // infeasible
#ifdef CONFLICT_CUTS
#undef CONFLICT_CUTS
            //#define CONFLICT_CUTS 2
#endif
#ifdef CONFLICT_CUTS
#ifdef COIN_HAS_CLP
            if (osiclp && (model->moreSpecialOptions() & 4194304) != 0) {
              const CbcFullNodeInfo *topOfTree = model->topOfTree();
              if (topOfTree) {
#if CONFLICT_CUTS == 2
                OsiRowCut *cut = osiclp->smallModelCut(topOfTree->lower(),
                  topOfTree->upper(),
                  model->numberRowsAtContinuous(),
                  model->whichGenerator());
#else
                OsiRowCut *cut = osiclp->modelCut(topOfTree->lower(),
                  topOfTree->upper(),
                  model->numberRowsAtContinuous(),
                  model->whichGenerator(), 0);
#endif
                if (cut) {
                  if (model->messageHandler()->logLevel() > 1)
                    printf("Conflict cut found in strong branching (%d elements)\n",
                      cut->row().getNumElements());
                  //cut->print();
                  if ((model->specialOptions() & 1) != 0) {
                    const OsiRowCutDebugger *debugger = model->continuousSolver()->getRowCutDebugger();
                    if (debugger) {
                      if (debugger->invalidCut(*cut)) {
                        model->continuousSolver()->applyRowCuts(1, cut);
                        model->continuousSolver()->writeMps("bad");
                      }
                      CoinAssert(!debugger->invalidCut(*cut));
                    }
                  }
                  model->makeGlobalCut(cut);
                }
              }
            }
#endif
#endif
          }
          // say infeasible if branch says so
          if (predictedChange == COIN_DBL_MAX)
            iStatus = 1;
          if (iStatus != 2 && solver->getIterationCount() > realMaxHotIterations)
            numberUnfinished++;
          newObjectiveValue = solver->getObjSense() * solver->getObjValue();
          choice.numItersDown = solver->getIterationCount();
          objectiveChange = CoinMax(newObjectiveValue - objectiveValue_, 0.0);
          // Update branching information if wanted
          CbcBranchingObject *cbcobj = dynamic_cast< CbcBranchingObject * >(choice.possibleBranch);
          if (cbcobj) {
            CbcObject *object = cbcobj->object();
            assert(object);
            CbcObjectUpdateData update = object->createUpdateInformation(solver, this, cbcobj);
            update.objectNumber_ = choice.objectNumber;
            model->addUpdateInformation(update);
          } else {
            decision->updateInformation(solver, this);
          }
          if (!iStatus) {
            choice.finishedDown = true;
            if (newObjectiveValue >= cutoff) {
              objectiveChange = 1.0e100; // say infeasible
              numberStrongInfeasible++;
            } else {
#define CBCNODE_TIGHTEN_BOUNDS
#ifdef CBCNODE_TIGHTEN_BOUNDS
              // Can we tighten bounds?
              if (iColumn < numberColumns && cutoff < 1.0e20
                && objectiveChange > 1.0e-5) {
                double value = saveSolution[iColumn];
                double down = value - floor(value - integerTolerance);
                double changePer = objectiveChange / (down + 1.0e-7);
                double distance = (cutoff - objectiveValue_) / changePer;
                distance += 1.0e-3;
                if (distance < 5.0) {
                  double newLower = ceil(value - distance);
                  if (newLower > saveLower[iColumn]) {
                    //printf("Could increase lower bound on %d from %g to %g\n",
                    //   iColumn,saveLower[iColumn],newLower);
                    saveLower[iColumn] = newLower;
                    solver->setColLower(iColumn, newLower);
                  }
                }
              }
#endif
              // See if integer solution
              feasibleSolution = model->feasibleSolution(choice.numIntInfeasDown,
                choice.numObjInfeasDown);
              if (feasibleSolution
                && model->problemFeasibility()->feasible(model, -1) >= 0) {
                if (auxiliaryInfo->solutionAddsCuts()) {
                  needHotStartUpdate = true;
                  solver->unmarkHotStart();
                }
                model->setLogLevel(saveLogLevel);
                model->setBestSolution(CBC_STRONGSOL,
                  newObjectiveValue,
                  solver->getColSolution());
                if (needHotStartUpdate) {
                  model->resolve(NULL, 11, saveSolution, saveLower, saveUpper);
                  newObjectiveValue = solver->getObjSense() * solver->getObjValue();
                  objectiveValue_ = CoinMax(objectiveValue_, newObjectiveValue);
                  objectiveChange = CoinMax(newObjectiveValue - objectiveValue_, 0.0);
                  model->feasibleSolution(choice.numIntInfeasDown,
                    choice.numObjInfeasDown);
                }
                model->setLastHeuristic(NULL);
                model->incrementUsed(solver->getColSolution());
                cutoff = model->getCutoff();
                if (newObjectiveValue >= cutoff) { //  *new* cutoff
                  objectiveChange = 1.0e100;
                  numberStrongInfeasible++;
                }
              }
            }
          } else if (iStatus == 1) {
            choice.finishedDown = true;
            objectiveChange = COIN_DBL_MAX;
            numberStrongInfeasible++;
          } else {
            // Can't say much as we did not finish
            choice.finishedDown = false;
            numberUnfinished++;
          }
          choice.downMovement = objectiveChange;

          // restore bounds
          for (j = 0; j < numberColumns; j++) {
            if (saveLower[j] != lower[j])
              solver->setColLower(j, saveLower[j]);
            if (saveUpper[j] != upper[j])
              solver->setColUpper(j, saveUpper[j]);
          }
          if (needHotStartUpdate) {
            needHotStartUpdate = false;
            model->resolve(NULL, 11, saveSolution, saveLower, saveUpper);
#ifdef CHECK_DEBUGGER_PATH
            if ((model->specialOptions() & 1) != 0 && onOptimalPath) {
              const OsiRowCutDebugger *debugger = solver->getRowCutDebugger();
              if (!debugger) {
                printf("Strong branching down on %d went off optimal path\n", iObject);
                model->solver()->writeMps("query");
                abort();
              }
            }
#endif
            double newObjValue = solver->getObjSense() * solver->getObjValue();
            objectiveValue_ = CoinMax(objectiveValue_, newObjValue);
            //we may again have an integer feasible solution
            int numberIntegerInfeasibilities;
            int numberObjectInfeasibilities;
            if (model->feasibleSolution(
                  numberIntegerInfeasibilities,
                  numberObjectInfeasibilities)) {
#ifdef BONMIN
              //In this case node has become integer feasible, let us exit the loop
              std::cout << "Node has become integer feasible" << std::endl;
              numberUnsatisfied_ = 0;
              break;
#endif
              double objValue = solver->getObjValue();
              model->setLogLevel(saveLogLevel);
              model->setBestSolution(CBC_STRONGSOL,
                objValue,
                solver->getColSolution());
              model->resolve(NULL, 11, saveSolution, saveLower, saveUpper);
              double newObjValue = solver->getObjSense() * solver->getObjValue();
              objectiveValue_ = CoinMax(objectiveValue_, newObjValue);
              cutoff = model->getCutoff();
            }
            solver->markHotStart();
#ifdef RESET_BOUNDS
            memcpy(saveLower, solver->getColLower(), solver->getNumCols() * sizeof(double));
            memcpy(saveUpper, solver->getColUpper(), solver->getNumCols() * sizeof(double));
#endif
            if (!solver->isProvenOptimal()) {
              skipAll = -2;
              canSkip = 1;
            }
            xMark++;
          }
#if 0 //def DO_ALL_AT_ROOT
                    if (strongType)
                        printf("Down on %d, status is %d, obj %g its %d cost %g finished %d inf %d infobj %d\n",
                               choice.objectNumber, iStatus, newObjectiveValue, choice.numItersDown,
                               choice.downMovement, choice.finishedDown, choice.numIntInfeasDown,
                               choice.numObjInfeasDown);
#endif

          // repeat the whole exercise, forcing the variable up
          predictedChange = choice.possibleBranch->branch();
          solver->solveFromHotStart();
#ifdef COIN_HAS_CLP
          if (osiclp && goToEndInStrongBranching) {
            osiclp->setIntParam(OsiMaxNumIterationHotStart, saveMaxHotIts);
            osiclp->setSpecialOptions(saveOsiClpOptions);
          }
#endif
          if ((model->moreSpecialOptions2() & 32768) != 0 && solver->isProvenOptimal()) {
            // If any small values re-do
            model->cleanBounds(solver, cleanVariables);
          }
          numberStrongDone++;
          numberStrongIterations += solver->getIterationCount();
          /*
                      We now have an estimate of objective degradation that we can use for strong
                      branching. If we're over the cutoff, the variable is monotone up.
                      If we actually made it to optimality, check for a solution, and if we have
                      a good one, call setBestSolution to process it. Note that this may reduce the
                      cutoff, so we check again to see if we can declare this variable monotone.
                    */
          if (solver->isProvenOptimal())
            iStatus = 0; // optimal
          else if (solver->isIterationLimitReached()
            && !solver->isDualObjectiveLimitReached()) {
            iStatus = 2; // unknown
          } else {
            iStatus = 1; // infeasible
#ifdef CONFLICT_CUTS
#ifdef COIN_HAS_CLP
            if (osiclp && (model->moreSpecialOptions() & 4194304) != 0) {
              const CbcFullNodeInfo *topOfTree = model->topOfTree();
              if (topOfTree) {
#if CONFLICT_CUTS == 2
                OsiRowCut *cut = osiclp->smallModelCut(topOfTree->lower(),
                  topOfTree->upper(),
                  model->numberRowsAtContinuous(),
                  model->whichGenerator());
#else
                OsiRowCut *cut = osiclp->modelCut(topOfTree->lower(),
                  topOfTree->upper(),
                  model->numberRowsAtContinuous(),
                  model->whichGenerator(), 0);
#endif
                if (cut) {
                  //printf("XXXXXX found conflict cut in strong branching\n");
                  //cut->print();
                  if ((model->specialOptions() & 1) != 0) {
                    const OsiRowCutDebugger *debugger = model->continuousSolver()->getRowCutDebugger();
                    if (debugger) {
                      if (debugger->invalidCut(*cut)) {
                        model->continuousSolver()->applyRowCuts(1, cut);
                        model->continuousSolver()->writeMps("bad");
                      }
                      CoinAssert(!debugger->invalidCut(*cut));
                    }
                  }
                  model->makeGlobalCut(cut);
                }
              }
            }
#endif
#endif
          }
          // say infeasible if branch says so
          if (predictedChange == COIN_DBL_MAX)
            iStatus = 1;
          if (iStatus != 2 && solver->getIterationCount() > realMaxHotIterations)
            numberUnfinished++;
          newObjectiveValue = solver->getObjSense() * solver->getObjValue();
          choice.numItersUp = solver->getIterationCount();
          objectiveChange = CoinMax(newObjectiveValue - objectiveValue_, 0.0);
          // Update branching information if wanted
          cbcobj = dynamic_cast< CbcBranchingObject * >(choice.possibleBranch);
          if (cbcobj) {
            CbcObject *object = cbcobj->object();
            assert(object);
            CbcObjectUpdateData update = object->createUpdateInformation(solver, this, cbcobj);
            update.objectNumber_ = choice.objectNumber;
            model->addUpdateInformation(update);
          } else {
            decision->updateInformation(solver, this);
          }
          if (!iStatus) {
            choice.finishedUp = true;
            if (newObjectiveValue >= cutoff) {
              objectiveChange = 1.0e100; // say infeasible
              numberStrongInfeasible++;
            } else {
#ifdef CBCNODE_TIGHTEN_BOUNDS
              // Can we tighten bounds?
              if (iColumn < numberColumns && cutoff < 1.0e20
                && objectiveChange > 1.0e-5) {
                double value = saveSolution[iColumn];
                double up = ceil(value + integerTolerance) - value;
                double changePer = objectiveChange / (up + 1.0e-7);
                double distance = (cutoff - objectiveValue_) / changePer;
                distance += 1.0e-3;
                if (distance < 5.0) {
                  double newUpper = floor(value + distance);
                  if (newUpper < saveUpper[iColumn]) {
                    //printf("Could decrease upper bound on %d from %g to %g\n",
                    //   iColumn,saveUpper[iColumn],newUpper);
                    saveUpper[iColumn] = newUpper;
                    solver->setColUpper(iColumn, newUpper);
                  }
                }
              }
#endif
              // See if integer solution
              feasibleSolution = model->feasibleSolution(choice.numIntInfeasUp,
                choice.numObjInfeasUp);
              if (feasibleSolution
                && model->problemFeasibility()->feasible(model, -1) >= 0) {
#ifdef BONMIN
                std::cout << "Node has become integer feasible" << std::endl;
                numberUnsatisfied_ = 0;
                break;
#endif
                if (auxiliaryInfo->solutionAddsCuts()) {
                  needHotStartUpdate = true;
                  solver->unmarkHotStart();
                }
                model->setLogLevel(saveLogLevel);
                model->setBestSolution(CBC_STRONGSOL,
                  newObjectiveValue,
                  solver->getColSolution());
                if (choice.finishedDown) {
                  double cutoff = model->getCutoff();
                  double downObj = objectiveValue_
                    + choice.downMovement;
                  if (downObj >= cutoff) {
                    choice.downMovement = 1.0e100;
                    numberStrongInfeasible++;
                  }
                }
                if (needHotStartUpdate) {
                  model->resolve(NULL, 11, saveSolution, saveLower, saveUpper);
#ifdef CHECK_DEBUGGER_PATH
                  if ((model->specialOptions() & 1) != 0 && onOptimalPath) {
                    const OsiRowCutDebugger *debugger = solver->getRowCutDebugger();
                    if (!debugger) {
                      printf("Strong branching up on %d went off optimal path\n", iObject);
                      abort();
                    }
                  }
#endif
                  newObjectiveValue = solver->getObjSense() * solver->getObjValue();
                  objectiveValue_ = CoinMax(objectiveValue_, newObjectiveValue);
                  objectiveChange = CoinMax(newObjectiveValue - objectiveValue_, 0.0);
                  model->feasibleSolution(choice.numIntInfeasDown,
                    choice.numObjInfeasDown);
                }
                model->setLastHeuristic(NULL);
                model->incrementUsed(solver->getColSolution());
                cutoff = model->getCutoff();
                if (newObjectiveValue >= cutoff) { //  *new* cutoff
                  objectiveChange = 1.0e100;
                  numberStrongInfeasible++;
                }
              }
            }
          } else if (iStatus == 1) {
            choice.finishedUp = true;
            objectiveChange = COIN_DBL_MAX;
            numberStrongInfeasible++;
          } else {
            // Can't say much as we did not finish
            choice.finishedUp = false;
            numberUnfinished++;
          }
          choice.upMovement = objectiveChange;

          // restore bounds
          for (j = 0; j < numberColumns; j++) {
            if (saveLower[j] != lower[j])
              solver->setColLower(j, saveLower[j]);
            if (saveUpper[j] != upper[j])
              solver->setColUpper(j, saveUpper[j]);
          }
          if (needHotStartUpdate) {
            needHotStartUpdate = false;
            model->resolve(NULL, 11, saveSolution, saveLower, saveUpper);
#ifdef CHECK_DEBUGGER_PATH
            if ((model->specialOptions() & 1) != 0 && onOptimalPath) {
              const OsiRowCutDebugger *debugger = solver->getRowCutDebugger();
              if (!debugger) {
                printf("Strong branching up on %d went off optimal path\n", iObject);
                abort();
              }
            }
#endif
            double newObjValue = solver->getObjSense() * solver->getObjValue();
            objectiveValue_ = CoinMax(objectiveValue_, newObjValue);
            //we may again have an integer feasible solution
            int numberIntegerInfeasibilities;
            int numberObjectInfeasibilities;
            if (model->feasibleSolution(
                  numberIntegerInfeasibilities,
                  numberObjectInfeasibilities)) {
              double objValue = solver->getObjValue();
              model->setLogLevel(saveLogLevel);
              model->setBestSolution(CBC_STRONGSOL,
                objValue,
                solver->getColSolution());
              model->resolve(NULL, 11, saveSolution, saveLower, saveUpper);
#ifdef CHECK_DEBUGGER_PATH
              if ((model->specialOptions() & 1) != 0 && onOptimalPath) {
                const OsiRowCutDebugger *debugger = solver->getRowCutDebugger();
                if (!debugger) {
                  printf("Strong branching up on %d went off optimal path\n", iObject);
                  abort();
                }
              }
#endif
              double newObjValue = solver->getObjSense() * solver->getObjValue();
              objectiveValue_ = CoinMax(objectiveValue_, newObjValue);
              cutoff = model->getCutoff();
            }
            solver->markHotStart();
#ifdef RESET_BOUNDS
            memcpy(saveLower, solver->getColLower(), solver->getNumCols() * sizeof(double));
            memcpy(saveUpper, solver->getColUpper(), solver->getNumCols() * sizeof(double));
#endif
            if (!solver->isProvenOptimal()) {
              skipAll = -2;
              canSkip = 1;
            }
            xMark++;
          }

#if 0 //def DO_ALL_AT_ROOT
                    if (strongType)
                        printf("Up on %d, status is %d, obj %g its %d cost %g finished %d inf %d infobj %d\n",
                               choice.objectNumber, iStatus, newObjectiveValue, choice.numItersUp,
                               choice.upMovement, choice.finishedUp, choice.numIntInfeasUp,
                               choice.numObjInfeasUp);
#endif
        }

        solver->setIntParam(OsiMaxNumIterationHotStart, saveLimit2);
        /*
                  End of evaluation for this candidate variable. Possibilities are:
                  * Both sides below cutoff; this variable is a candidate for branching.
                  * Both sides infeasible or above the objective cutoff: no further action
                  here. Break from the evaluation loop and assume the node will be purged
                  by the caller.
                  * One side below cutoff: Install the branch (i.e., fix the variable). Break
                  from the evaluation loop and assume the node will be reoptimised by the
                  caller.
                */
        // reset
        choice.possibleBranch->resetNumberBranchesLeft();
        if (choice.upMovement < 1.0e100) {
          if (choice.downMovement < 1.0e100) {
            // In case solution coming in was odd
            choice.upMovement = CoinMax(0.0, choice.upMovement);
            choice.downMovement = CoinMax(0.0, choice.downMovement);
#if ZERO_ONE == 2
            // branch on 0-1 first (temp)
            if (fabs(choice.possibleBranch->value()) < 1.0) {
              choice.upMovement *= ZERO_FAKE;
              choice.downMovement *= ZERO_FAKE;
            }
#endif
            // feasible - see which best
            if (!canSkip) {
              if (model->messageHandler()->logLevel() > 3)
                printf("sort %g downest %g upest %g ", sort[iDo], downEstimate[iObject],
                  upEstimate[iObject]);
              model->messageHandler()->message(CBC_STRONG, *model->messagesPointer())
                << iObject << iColumn
                << choice.downMovement << choice.numIntInfeasDown
                << choice.upMovement << choice.numIntInfeasUp
                << choice.possibleBranch->value()
                << CoinMessageEol;
            }
            int betterWay = 0;
            // If was feasible (extra strong branching) skip
            if (infeasibility) {
              CbcBranchingObject *branchObj = dynamic_cast< CbcBranchingObject * >(branch_);
              if (branch_)
                assert(branchObj);
              betterWay = decision->betterBranch(choice.possibleBranch,
                branchObj,
                choice.upMovement,
                choice.numIntInfeasUp,
                choice.downMovement,
                choice.numIntInfeasDown);
            }
            if (betterWay) {
              // C) create branching object
              if (choiceObject) {
                delete branch_;
                branch_ = choice.possibleBranch->clone();
              } else {
                delete branch_;
                branch_ = choice.possibleBranch;
                choice.possibleBranch = NULL;
              }
              {
                CbcBranchingObject *branchObj = dynamic_cast< CbcBranchingObject * >(branch_);
                assert(branchObj);
                //branchObj->way(preferredWay);
                branchObj->way(betterWay);
              }
              bestChoice = choice.objectNumber;
              whichChoice = iDo;
              if (numberStrong <= 1) {
                delete ws;
                ws = NULL;
                break;
              }
            } else {
              if (!choiceObject) {
                delete choice.possibleBranch;
                choice.possibleBranch = NULL;
              }
              if (iDo >= 2 * numberStrong) {
                delete ws;
                ws = NULL;
                break;
              }
              if (!dynamicObject || dynamicObject->numberTimesUp() > 1) {
                if (iDo - whichChoice >= numberStrong) {
                  if (!choiceObject) {
                    delete choice.possibleBranch;
                    choice.possibleBranch = NULL;
                  }
                  break; // give up
                }
              } else {
                if (iDo - whichChoice >= 2 * numberStrong) {
                  delete ws;
                  ws = NULL;
                  if (!choiceObject) {
                    delete choice.possibleBranch;
                    choice.possibleBranch = NULL;
                  }
                  break; // give up
                }
              }
            }
          } else {
            // up feasible, down infeasible
            anyAction = -1;
            worstFeasible = CoinMax(worstFeasible, choice.upMovement);
            model->messageHandler()->message(CBC_STRONG, *model->messagesPointer())
              << iObject << iColumn
              << choice.downMovement << choice.numIntInfeasDown
              << choice.upMovement << choice.numIntInfeasUp
              << choice.possibleBranch->value()
              << CoinMessageEol;
            //printf("Down infeasible for choice %d sequence %d\n",i,
            // model->object(choice.objectNumber)->columnNumber());
            choice.fix = 1;
            numberToFix++;
            choice.possibleBranch->fix(solver, saveLower, saveUpper, 1);
            if (!choiceObject) {
              delete choice.possibleBranch;
              choice.possibleBranch = NULL;
            } else {
              //choiceObject = new CbcDynamicPseudoCostBranchingObject(*choiceObject);
              choice.possibleBranch = choiceObject;
            }
            assert(doneHotStart);
            solver->unmarkHotStart();
            model->resolve(NULL, 11, saveSolution, saveLower, saveUpper);
#ifdef CHECK_DEBUGGER_PATH
            if ((model->specialOptions() & 1) != 0 && onOptimalPath) {
              const OsiRowCutDebugger *debugger = solver->getRowCutDebugger();
              if (!debugger) {
                printf("Strong branching down on %d went off optimal path\n", iObject);
                abort();
              }
            }
#endif
            double newObjValue = solver->getObjSense() * solver->getObjValue();
            objectiveValue_ = CoinMax(objectiveValue_, newObjValue);
            bool goneInfeasible = (!solver->isProvenOptimal() || solver->isDualObjectiveLimitReached());
            solver->markHotStart();
#ifdef RESET_BOUNDS
            memcpy(saveLower, solver->getColLower(), solver->getNumCols() * sizeof(double));
            memcpy(saveUpper, solver->getColUpper(), solver->getNumCols() * sizeof(double));
#endif
            if (!solver->isProvenOptimal()) {
              skipAll = -2;
              canSkip = 1;
            }
            xMark++;
            // may be infeasible (if other way stopped on iterations)
            if (goneInfeasible) {
              // neither side feasible
              anyAction = -2;
              if (!choiceObject) {
                delete choice.possibleBranch;
                choice.possibleBranch = NULL;
              }
              //printf("Both infeasible for choice %d sequence %d\n",i,
              // model->object(choice.objectNumber)->columnNumber());
              delete ws;
              ws = NULL;
              break;
            }
          }
        } else {
          if (choice.downMovement < 1.0e100) {
            // down feasible, up infeasible
            anyAction = -1;
            worstFeasible = CoinMax(worstFeasible, choice.downMovement);
            model->messageHandler()->message(CBC_STRONG, *model->messagesPointer())
              << iObject << iColumn
              << choice.downMovement << choice.numIntInfeasDown
              << choice.upMovement << choice.numIntInfeasUp
              << choice.possibleBranch->value()
              << CoinMessageEol;
            choice.fix = -1;
            numberToFix++;
            choice.possibleBranch->fix(solver, saveLower, saveUpper, -1);
            if (!choiceObject) {
              delete choice.possibleBranch;
              choice.possibleBranch = NULL;
            } else {
              //choiceObject = new CbcDynamicPseudoCostBranchingObject(*choiceObject);
              choice.possibleBranch = choiceObject;
            }
            assert(doneHotStart);
            solver->unmarkHotStart();
            model->resolve(NULL, 11, saveSolution, saveLower, saveUpper);
#ifdef CHECK_DEBUGGER_PATH
            if ((model->specialOptions() & 1) != 0 && onOptimalPath) {
              const OsiRowCutDebugger *debugger = solver->getRowCutDebugger();
              if (!debugger) {
                printf("Strong branching down on %d went off optimal path\n", iObject);
                solver->writeMps("query");
                abort();
              }
            }
#endif
            double newObjValue = solver->getObjSense() * solver->getObjValue();
            objectiveValue_ = CoinMax(objectiveValue_, newObjValue);
            bool goneInfeasible = (!solver->isProvenOptimal() || solver->isDualObjectiveLimitReached());
            solver->markHotStart();
#ifdef RESET_BOUNDS
            memcpy(saveLower, solver->getColLower(), solver->getNumCols() * sizeof(double));
            memcpy(saveUpper, solver->getColUpper(), solver->getNumCols() * sizeof(double));
#endif
            if (!solver->isProvenOptimal()) {
              skipAll = -2;
              canSkip = 1;
            }
            xMark++;
            // may be infeasible (if other way stopped on iterations)
            if (goneInfeasible) {
              // neither side feasible
              anyAction = -2;
              if (!choiceObject) {
                delete choice.possibleBranch;
                choice.possibleBranch = NULL;
              }
              delete ws;
              ws = NULL;
              break;
            }
          } else {
            // neither side feasible
            anyAction = -2;
            if (!choiceObject) {
              delete choice.possibleBranch;
              choice.possibleBranch = NULL;
            }
            delete ws;
            ws = NULL;
            break;
          }
        }
        // Check max time
        hitMaxTime = (model->getCurrentSeconds() > model->getDblParam(CbcModel::CbcMaximumSeconds));
        if (hitMaxTime) {
          // make sure rest are fast
          for (int jDo = iDo + 1; jDo < numberToDo; jDo++) {
            int iObject = whichObject[iDo];
            OsiObject *object = model->modifiableObject(iObject);
            CbcSimpleIntegerDynamicPseudoCost *dynamicObject = dynamic_cast< CbcSimpleIntegerDynamicPseudoCost * >(object);
            if (dynamicObject)
              dynamicObject->setNumberBeforeTrust(0);
          }
          numberTest = 0;
        }
        if (!choiceObject) {
          delete choice.possibleBranch;
        }
      }
      if (model->messageHandler()->logLevel() > 3) {
        if (anyAction == -2) {
          printf("infeasible\n");
        } else if (anyAction == -1) {
          printf("%d fixed AND choosing %d iDo %d iChosenWhen %d numberToDo %d\n", numberToFix, bestChoice,
            iDo, whichChoice, numberToDo);
        } else {
          int iObject = whichObject[whichChoice];
          OsiObject *object = model->modifiableObject(iObject);
          CbcSimpleIntegerDynamicPseudoCost *dynamicObject = dynamic_cast< CbcSimpleIntegerDynamicPseudoCost * >(object);
          if (dynamicObject) {
            int iColumn = dynamicObject->columnNumber();
            printf("choosing %d (column %d) iChosenWhen %d numberToDo %d\n", bestChoice,
              iColumn, whichChoice, numberToDo);
          }
        }
      }
      if (doneHotStart) {
        // Delete the snapshot
        solver->unmarkHotStart();
        // back to normal
        solver->setHintParam(OsiDoInBranchAndCut, true, OsiHintDo, NULL);
        // restore basis
        solver->setWarmStart(ws);
      }
      solver->setIntParam(OsiMaxNumIterationHotStart, saveLimit);
      // Unless infeasible we will carry on
      // But we could fix anyway
      if (numberToFix && !hitMaxTime) {
        if (anyAction != -2) {
          // apply and take off
          bool feasible = true;
          // can do quick optimality check
          int easy = 2;
          solver->setHintParam(OsiDoInBranchAndCut, true, OsiHintDo, &easy);
          model->resolve(NULL, 11, saveSolution, saveLower, saveUpper);
#ifdef CHECK_DEBUGGER_PATH
          if ((model->specialOptions() & 1) != 0 && onOptimalPath) {
            const OsiRowCutDebugger *debugger = solver->getRowCutDebugger();
            if (!debugger) {
              printf("Strong branching went off optimal path\n");
              abort();
            }
          }
#endif
          double newObjValue = solver->getObjSense() * solver->getObjValue();
          objectiveValue_ = CoinMax(objectiveValue_, newObjValue);
          solver->setHintParam(OsiDoInBranchAndCut, true, OsiHintDo, NULL);
          feasible = solver->isProvenOptimal();
          if (feasible) {
            anyAction = 0;
          } else {
            anyAction = -2;
            finished = true;
          }
        }
      }
      // If  fixed then round again
      // See if candidate still possible
      if (branch_) {
        const OsiObject *object = model->object(bestChoice);
        double infeasibility = object->checkInfeasibility(&usefulInfo);
        if (!infeasibility) {
          // take out
          delete branch_;
          branch_ = NULL;
        } else {
          // get preferred way
          int preferredWay;
          object->infeasibility(&usefulInfo, preferredWay);
          CbcBranchingObject *branchObj = dynamic_cast< CbcBranchingObject * >(branch_);
          assert(branchObj);
          branchObj->way(preferredWay);
#ifdef CBCNODE_TIGHTEN_BOUNDS
          bool fixed = branchObj->tighten(solver);
          if (fixed) {
            //printf("Variable now fixed!\n");
            // take out
            delete branch_;
            branch_ = NULL;
          }
#endif
        }
      }
      if (!branch_ && anyAction != -2 && !hitMaxTime) {
        finished = false;
      }
      delete ws;
    }
  }
  // update number of strong iterations etc
  model->incrementStrongInfo(numberStrongDone, numberStrongIterations,
    anyAction == -2 ? 0 : numberToFix, anyAction == -2);
  if (model->searchStrategy() == -1) {
#ifndef COIN_DEVELOP
    if (solver->messageHandler()->logLevel() > 1)
#endif
      printf("%d strong, %d iters, %d inf, %d not finished, %d not trusted\n",
        numberStrongDone, numberStrongIterations, numberStrongInfeasible, numberUnfinished,
        numberNotTrusted);
    // decide what to do
    int strategy = 1;
    if (((numberUnfinished * 4 > numberStrongDone && numberStrongInfeasible * 40 < numberStrongDone) || numberStrongInfeasible < 0) && model->numberStrong() < 10 && model->numberBeforeTrust() <= 20 && model->numberObjects() > CoinMax(1000, solver->getNumRows())) {
      strategy = 2;
#ifdef COIN_DEVELOP
      //if (model->logLevel()>1)
      printf("going to strategy 2\n");
#endif
      // Weaken
      model->setNumberStrong(2);
      model->setNumberBeforeTrust(1);
      model->synchronizeNumberBeforeTrust();
    }
    if (numberNodes)
      strategy = 1; // should only happen after hot start
    model->setSearchStrategy(strategy);
  } else if (numberStrongDone) {
    //printf("%d strongB, %d iters, %d inf, %d not finished, %d not trusted\n",
    //   numberStrongDone,numberStrongIterations,numberStrongInfeasible,numberUnfinished,
    //   numberNotTrusted);
  }
  if (model->searchStrategy() == 1 && numberNodes > 500 && numberNodes < -510) {
#ifndef COIN_DEVELOP
    if (solver->messageHandler()->logLevel() > 1)
#endif
      printf("after %d nodes - %d strong, %d iters, %d inf, %d not finished, %d not trusted\n",
        numberNodes, numberStrongDone, numberStrongIterations, numberStrongInfeasible, numberUnfinished,
        numberNotTrusted);
    // decide what to do
    if (numberUnfinished * 10 > numberStrongDone + 1 || !numberStrongInfeasible) {
      COIN_DETAIL_PRINT(printf("going to strategy 2\n"));
      // Weaken
      model->setNumberStrong(2);
      model->setNumberBeforeTrust(1);
      model->synchronizeNumberBeforeTrust();
      model->setSearchStrategy(2);
    }
  }
  if (numberUnfinished * 10 < numberStrongDone && model->numberStrongIterations() * 20 < model->getIterationCount() && !auxiliaryInfo->solutionAddsCuts()) {
    //printf("increasing trust\n");
    model->synchronizeNumberBeforeTrust(2);
  }

  // Set guessed solution value
  guessedObjectiveValue_ = objectiveValue_ + estimatedDegradation;
  int kColumn = -1;
  if (branch_) {
    CbcObject *obj = (dynamic_cast< CbcBranchingObject * >(branch_))->object();
    CbcSimpleInteger *branchObj = dynamic_cast< CbcSimpleInteger * >(obj);
    if (branchObj) {
      kColumn = branchObj->columnNumber();
    }
  }
#ifdef COIN_HAS_NTY
  if (orbitOption && kColumn >= 0) {
    CbcSymmetry *symmetryInfo = model->symmetryInfo();
    CbcNodeInfo *infoX = lastNode ? lastNode->nodeInfo() : NULL;
    bool worthTrying = false;
    if (infoX) {
      CbcNodeInfo *info = infoX;
      for (int i = 0; i < NTY_BAD_DEPTH; i++) {
        if (!info->parent()) {
          worthTrying = true;
          break;
        }
        info = info->parent();
        if (info->symmetryWorked()) {
          worthTrying = true;
          break;
        }
      }
    } else {
      worthTrying = true;
    }
    if (orbitOption == 3 && depth_ > 5)
      worthTrying = false;
    if (symmetryInfo && worthTrying) {
      if ((orbitOption & 1) == 1) {
        symmetryInfo->ChangeBounds(solver->getColLower(),
          solver->getColUpper(),
          solver->getNumCols(), false);
        symmetryInfo->Compute_Symmetry();
        symmetryInfo->fillOrbits();
      }
      const int *orbits = symmetryInfo->whichOrbit();
      if (orbits && orbits[kColumn] >= 0) {
        int numberUsefulOrbits = symmetryInfo->numberUsefulOrbits();
        if (solver->messageHandler()->logLevel() > 1)
          printf("Orbital Branching on %d - way %d n %d\n", kColumn, way(), numberUsefulOrbits);
        if (numberUsefulOrbits < 1000 || orbitOption == 3) {
          delete branch_;
          branch_ = new CbcOrbitalBranchingObject(model, kColumn, 1, 0, NULL);
          if (infoX)
            infoX->setSymmetryWorked();
        }
      }
    }
  }
#endif
  if (model->logLevel() > 1)
    printf("Node %d depth %d unsatisfied %d sum %g obj %g guess %g branching on %d\n",
      model->getNodeCount(), depth_, numberUnsatisfied_,
      sumInfeasibilities_, objectiveValue_, guessedObjectiveValue_,
      kColumn);
#ifdef DO_ALL_AT_ROOT
  if (strongType) {
    char general[200];
    if (strongType == 1)
      sprintf(general, "Strong branching on all %d unsatisfied, %d iterations (depth %d)\n",
        saveNumberToDo, numberStrongIterations, depth_);
    else
      sprintf(general, "Strong branching on all %d unfixed variables (%d unsatisfied), %d iterations (depth %d)\n",
        saveNumberToDo + saveSatisfiedVariables, saveNumberToDo, numberStrongIterations, depth_);
    model->messageHandler()->message(CBC_FPUMP2, model->messages())
      << general << CoinMessageEol;
  }
#endif
#ifdef DEBUG_SOLUTION
  if (onOptimalPath && anyAction == -2) {
    printf("Gone off optimal path in CbcNode\n");
    assert(!onOptimalPath || anyAction != -2);
  }
#endif
  delete[] cleanVariables;
  /*
      Cleanup, then we're finished
    */
  if (!model->branchingMethod())
    delete decision;

  delete choiceObject;
  delete[] sort;
  delete[] whichObject;
#ifdef RANGING
  delete[] objectMark;
#endif
  delete[] saveLower;
  delete[] saveUpper;
  delete[] upEstimate;
  delete[] downEstimate;
#ifdef COIN_HAS_CLP
  if (osiclp) {
    osiclp->setSpecialOptions(saveClpOptions);
  }
#endif
  // restore solution
  solver->setColSolution(saveSolution);
  model->reserveCurrentSolution(saveSolution);
  delete[] saveSolution;
  model->setStateOfSearch(saveStateOfSearch);
  model->setLogLevel(saveLogLevel);
  // delete extra regions
  if (usefulInfo.usefulRegion_) {
    delete[] usefulInfo.usefulRegion_;
    delete[] usefulInfo.indexRegion_;
    delete[] usefulInfo.pi_;
    usefulInfo.usefulRegion_ = NULL;
    usefulInfo.indexRegion_ = NULL;
    usefulInfo.pi_ = NULL;
  }
  useShadow = model->moreSpecialOptions() & 7;
  if ((useShadow == 5 && model->getSolutionCount()) || useShadow == 6) {
    // zap pseudo shadow prices
    model->pseudoShadow(-1);
    // and switch off
    model->setMoreSpecialOptions(model->moreSpecialOptions() & (~1023));
  }
  return anyAction;
}
// 0 is down, 1 is up
typedef struct {
  double initialValue; // initial value
  double upLowerBound; // Lower bound when going up
  double downUpperBound; // Upper bound when going down
  double movement[2]; // cost  (and initial away from feasible)
  double sumModified[2]; // Sum of integer changes
  int modified[2]; // Number integers changed
  int numIntInfeas[2]; // without odd ones
  int numObjInfeas[2]; // just odd ones
  bool finished[2]; // true if solver finished
  int numIters[2]; // number of iterations in solver (-1 if never solved)
  double *integerSolution; // output if thinks integer solution
#ifdef COIN_HAS_CLP
  ClpDualRowSteepest *steepest;
#endif
  int columnNumber; // Which column it is
} StrongInfo;
typedef struct {
  double integerTolerance;
  double *originalSolution;
  CoinWarmStart *ws;
  double *newObjective;
#ifdef COIN_HAS_CLP
  ClpDualRowSteepest *dualRowPivot;
  ClpPrimalColumnPivot *primalColumnPivot;
#endif
  int *back;
  int solveType;
} StrongStaticInfo;
typedef struct {
  StrongStaticInfo *staticInfo;
  StrongInfo *choice;
  OsiSolverInterface *solver;
  double *tempSolution;
  CoinWarmStart *tempBasis;
  int whichChoice;
} StrongBundle;
/* return 1 if possible solution (for solveType 100 if infeasible)
   2 set if down was infeasible
   4 set if up was infeasible
 */
int solveAnalyze(void *info)
{
  StrongBundle *bundle = reinterpret_cast< StrongBundle * >(info);
  StrongInfo *choice = bundle->choice;
  StrongStaticInfo *staticInfo = bundle->staticInfo;
  OsiSolverInterface *solver = bundle->solver;
  int solveType = staticInfo->solveType;
  if (solveType == 77) {
    return 0;
  }
  const double *saveSolution = staticInfo->originalSolution;
  int iColumn = choice->columnNumber;
  const int *back = staticInfo->back;
  double newObjectiveValue = 1.0e100;
  double integerTolerance = staticInfo->integerTolerance;
  double bestSolutionValue = COIN_DBL_MAX;
  int returnStatus = 0;
  // status is 0 finished, 1 infeasible and other
  int iStatus;
  /*
    Try the down direction first. (Specify the initial branching alternative as
    down with a call to way(-1). Each subsequent call to branch() performs the
    specified branch and advances the branch object state to the next branch
    alternative.)
  */
  for (int iWay = 0; iWay < 2; iWay++) {
    if (choice->numIters[iWay] == 0) {
      int numberColumns = solver->getNumCols();
      if (solveType != 100) {
        double saveBound;
        if (iWay == 0) {
          saveBound = solver->getColUpper()[iColumn];
          solver->setColUpper(iColumn, choice->downUpperBound);
        } else {
          saveBound = solver->getColLower()[iColumn];
          solver->setColLower(iColumn, choice->upLowerBound);
        }
        if ((solveType & 2) == 0) {
          solver->solveFromHotStart();
        } else {
          // restore basis
          solver->setWarmStart(staticInfo->ws);
#ifdef COIN_HAS_CLP
          if (staticInfo->dualRowPivot) {
            OsiClpSolverInterface *osiclp = dynamic_cast< OsiClpSolverInterface * >(solver);
            ClpSimplex *simplex = osiclp->getModelPtr();
            simplex->setDualRowPivotAlgorithm(*staticInfo->dualRowPivot);
            //simplex->dualRowPivot()->saveWeights(simplex,4);
            simplex->setWhatsChanged(ALL_SAME_EXCEPT_COLUMN_BOUNDS);
            simplex->dual(0, 5);
          } else {
#endif
            solver->resolve();
#ifdef COIN_HAS_CLP
          }
#endif
        }
        if (iWay == 0)
          solver->setColUpper(iColumn, saveBound);
        else
          solver->setColLower(iColumn, saveBound);
        /*
	  We now have an estimate of objective degradation that we can use for strong
	  branching. If we're over the cutoff, the variable is monotone up.
	  If we actually made it to optimality, check for a solution, and if we have
	  a good one, call setBestSolution to process it. Note that this may reduce the
	  cutoff, so we check again to see if we can declare this variable monotone.
	*/
        if (solver->isProvenOptimal()) {
          iStatus = 0; // optimal
        } else if (solver->isIterationLimitReached()
          && !solver->isDualObjectiveLimitReached()) {
          iStatus = 2; // unknown
        } else {
          iStatus = 1; // infeasible
        }
        newObjectiveValue = solver->getObjSense() * solver->getObjValue();
        choice->numIters[iWay] = solver->getIterationCount();
        // Look at interaction
        const double *thisSolution = solver->getColSolution();
        int numberModified = 0;
        double sumModified = 0.0;
        int numberInfeas = 0;
        for (int i = 0; i < numberColumns; i++) {
          if (back[i] >= 0) {
            double value = thisSolution[i];
            if (iColumn != i) {
              double difference = fabs(saveSolution[i] - value);
              if (difference > integerTolerance) {
                numberModified++;
                sumModified += difference;
              }
            }
            if (fabs(value - floor(value + 0.5)) > integerTolerance)
              numberInfeas++;
            ;
          }
        }
        choice->numIntInfeas[iWay] = numberInfeas;
        choice->sumModified[iWay] = sumModified;
        choice->modified[iWay] = numberModified;
        if (!iStatus) {
          choice->finished[iWay] = true;
          if (!numberInfeas) {
            returnStatus = 1;
            if (!choice->integerSolution) {
              bestSolutionValue = newObjectiveValue;
              choice->integerSolution = CoinCopyOfArray(thisSolution, numberColumns);
              ;
            } else if (bestSolutionValue > newObjectiveValue) {
              memcpy(choice->integerSolution, thisSolution, numberColumns * sizeof(double));
            }
          }
        } else if (iStatus == 1) {
          newObjectiveValue = 1.0e100;
        } else {
          // Can't say much as we did not finish
          choice->finished[iWay] = false;
        }
        choice->movement[iWay] = newObjectiveValue;
      } else {
#ifdef COIN_HAS_CLP
        OsiClpSolverInterface *osiclp = dynamic_cast< OsiClpSolverInterface * >(solver);
        ClpSimplex *simplex = osiclp ? osiclp->getModelPtr() : NULL;
#endif
        // doing continuous and general integer
        solver->setColSolution(staticInfo->originalSolution);
        solver->setWarmStart(staticInfo->ws);
        double saveBound;
        double newBound;
        if (iWay == 0) {
          saveBound = solver->getColUpper()[iColumn];
          solver->setColUpper(iColumn, choice->downUpperBound);
          newBound = choice->downUpperBound;
        } else {
          saveBound = solver->getColLower()[iColumn];
          solver->setColLower(iColumn, choice->upLowerBound);
          newBound = choice->upLowerBound;
        }
#if 0 //def COIN_HAS_CLP
	if (simplex) {
	  // set solution to new bound (if basic will be recomputed)
	  simplex->primalColumnSolution()[iColumn]=newBound;
	}
#endif
        solver->setHintParam(OsiDoDualInResolve, true, OsiHintDo);
#define PRINT_ANALYZE 0
#if PRINT_ANALYZE > 0
        osiclp->getModelPtr()->setLogLevel(1);
        solver->setHintParam(OsiDoReducePrint, false, OsiHintTry);
#endif
        solver->resolve();
        if (iWay == 0) {
#if PRINT_ANALYZE > 0
          printf("column %d down original %g <= %g <= %g upper now %g - result %s\n",
            iColumn, solver->getColLower()[iColumn],
            staticInfo->originalSolution[iColumn], saveBound,
            newBound, solver->isProvenOptimal() ? "ok" : "infeas");
#endif
          solver->setColUpper(iColumn, saveBound);
        } else {
#if PRINT_ANALYZE > 0
          printf("column %d up original %g <= %g <= %g lower now %g - result %s\n",
            iColumn, saveBound, staticInfo->originalSolution[iColumn],
            solver->getColUpper()[iColumn],
            newBound, solver->isProvenOptimal() ? "ok" : "infeas");
#endif
          solver->setColLower(iColumn, saveBound);
        }
        choice->numIters[iWay] = solver->getIterationCount();
        if (solver->isProvenOptimal()) {
          //printf("Way %d - all way %d iterations - column %d\n",
          //	 iWay,solver->getIterationCount(),iColumn);
          // can go all way
          choice->movement[iWay] = newBound;
        } else {
          // zero objective
          double offset;
          solver->getDblParam(OsiObjOffset, offset);
          solver->setDblParam(OsiObjOffset, 0.0);
          solver->setObjective(staticInfo->newObjective + numberColumns);
          if (iWay == 0) {
            solver->setObjCoeff(iColumn, 1.0);
          } else {
            solver->setObjCoeff(iColumn, -1.0);
          }
          solver->setColSolution(staticInfo->originalSolution);
          solver->setWarmStart(staticInfo->ws);
          solver->setHintParam(OsiDoDualInResolve, false, OsiHintDo);
          solver->resolve();
          //printf("Way %d - first solve %d iterations, second %d - column %d\n",
          //	 iWay,choice->numIters[iWay],solver->getIterationCount(),iColumn);
          choice->movement[iWay] = solver->getColSolution()[iColumn];
          choice->numIters[iWay] += solver->getIterationCount();
#if PRINT_ANALYZE > 0
          if (iWay == 0) {
            printf("column %d down can get to %g - result %s\n",
              iColumn, solver->getColSolution()[iColumn], solver->isProvenOptimal() ? "ok" : "infeas");
          } else {
            printf("column %d up can get to %g - result %s\n",
              iColumn, solver->getColSolution()[iColumn], solver->isProvenOptimal() ? "ok" : "infeas");
          }
#endif
          // reset objective
          solver->setDblParam(OsiObjOffset, offset);
          solver->setObjective(staticInfo->newObjective);
          if (!solver->isProvenOptimal()) {
#ifdef COIN_HAS_CLP
            OsiClpSolverInterface *osiclp = dynamic_cast< OsiClpSolverInterface * >(solver);
            ClpSimplex *simplex = osiclp->getModelPtr();
            double sum = simplex->sumPrimalInfeasibilities();
            sum /= static_cast< double >(simplex->numberPrimalInfeasibilities());
            if (sum > 1.0e-3) {
#endif
              choice->modified[0] = 1;
              returnStatus = 1;
              solver->writeMps("bad", "mps");
              abort();
#ifdef COIN_HAS_CLP
            }
#endif
          }
        }
        //solver->setObjCoeff(iColumn,0.0);
      }
    }
  }
  return returnStatus;
}
#ifdef THREADS_IN_ANALYZE
void *cbc_parallelManager(void *stuff)
{
  CoinPthreadStuff *driver = reinterpret_cast< CoinPthreadStuff * >(stuff);
  int whichThread = driver->whichThread();
  CoinThreadInfo *threadInfo = driver->threadInfoPointer(whichThread);
  threadInfo->status = -1;
  int *which = threadInfo->stuff;
  pthread_barrier_wait(driver->barrierPointer());
#if 0
  int status=-1;
  while (status!=100)
    status=timedWait(driver,1000,2);
  pthread_cond_signal(driver->conditionPointer(1));
  pthread_mutex_unlock(driver->mutexPointer(1,whichThread));
#endif
  // so now mutex_ is locked
  int whichLocked = 0;
  while (true) {
    pthread_mutex_t *mutexPointer = driver->mutexPointer(whichLocked, whichThread);
    // wait
    //printf("Child waiting for %d - status %d %d %d\n",
    //	   whichLocked,lockedX[0],lockedX[1],lockedX[2]);
#ifdef DETAIL_THREAD
    printf("thread %d about to lock mutex %d\n", whichThread, whichLocked);
#endif
    pthread_mutex_lock(mutexPointer);
    whichLocked++;
    if (whichLocked == 3)
      whichLocked = 0;
    int unLock = whichLocked + 1;
    if (unLock == 3)
      unLock = 0;
    //printf("child pointer %p status %d\n",threadInfo,threadInfo->status);
    assert(threadInfo->status >= 0);
    if (threadInfo->status == 1000)
      pthread_exit(NULL);
    int type = threadInfo->status;
    int &returnCode = which[0];
    int iPass = which[1];
    //CoinIndexedVector * array;
    //double dummy;
    switch (type) {
      // dummy
    case 0:
      break;
    case 1:
      returnCode = solveAnalyze(threadInfo->extraInfo);
      threadInfo->stuff[3] = 0;
      break;
    case 100:
      // initialization
      break;
    }
    threadInfo->status = (type != 1) ? -1 : -2;
#ifdef DETAIL_THREAD
    printf("thread %d about to unlock mutex %d\n", whichThread, unLock);
#endif
    pthread_mutex_unlock(driver->mutexPointer(unLock, whichThread));
  }
}
#endif
int CbcNode::analyze(CbcModel *model, double *results)
{
#define COIN_DETAIL
  int i;
  int numberIterationsAllowed = model->numberAnalyzeIterations();
  int numberColumns = model->getNumCols();
  int numberRows = model->getNumRows();
  int numberObjects = model->numberObjects();
  int numberIntegers = model->numberIntegers();
  int numberLookIntegers = 0;
  int highestPriority = COIN_INT_MAX;
  int *back = new int[numberColumns];
  const int *integerVariable = model->integerVariable();
  for (i = 0; i < numberIntegers; i++) {
    highestPriority = CoinMin(highestPriority, model->modifiableObject(i)->priority());
  }
  for (i = 0; i < numberColumns; i++)
    back[i] = -1;
  for (i = 0; i < numberIntegers; i++) {
    int iColumn = integerVariable[i];
    back[iColumn] = i;
    if (model->modifiableObject(i)->priority() == highestPriority) {
      numberLookIntegers++;
    } else {
      back[iColumn] = i + numberColumns;
    }
  }
  /*
    0 - just look 
    0 (1) bit - use to set priorities 
    1 (2) bit - look at bounds on all variables and more iterations
    2 (4) bit - do threaded (if parallelMode()==1 then not repeatable if any fixed)
    3 (8) bit - 
    4 (16) bit - do even if m*n>1,000,000
    5 (32) bit - printing time
    6 (64) bit - save mps file
  */
  int solveType;
  char general[200];
  if (numberIterationsAllowed > 0) {
    solveType = 0;
  } else {
    solveType = -numberIterationsAllowed;
    if ((solveType & 16) == 0) {
      double size = numberRows;
      size *= numberLookIntegers;
      if (size > 1000000) {
        if ((solveType & 32) != 0)
          model->messageHandler()->message(CBC_GENERAL, *model->messagesPointer())
            << "Skipping analyze as problem too large"
            << CoinMessageEol;
        return 0;
      }
    }
    sprintf(general, "Analyze options %d", solveType);
    model->messageHandler()->message(CBC_GENERAL, *model->messagesPointer())
      << general
      << CoinMessageEol;
    if ((solveType & 1) != 0)
      model->messageHandler()->message(CBC_GENERAL, *model->messagesPointer())
        << "Using to set priorities (probably bad idea)"
        << CoinMessageEol;
    if ((solveType & 2) != 0)
      model->messageHandler()->message(CBC_GENERAL, *model->messagesPointer())
        << "Use more iterations and look at continuous/general integer variables"
        << CoinMessageEol;
    if ((solveType & 4) != 0)
      model->messageHandler()->message(CBC_GENERAL, *model->messagesPointer())
        << "Use threads"
        << CoinMessageEol;
    if ((solveType & 32) != 0)
      model->messageHandler()->message(CBC_GENERAL, *model->messagesPointer())
        << "32 switches on more printing, (16 bit allows large problems)"
        << CoinMessageEol;
  }
  OsiSolverInterface *solver = model->solver();
  objectiveValue_ = solver->getObjSense() * solver->getObjValue();
  const double *lower = solver->getColLower();
  const double *upper = solver->getColUpper();
  const double *dj = solver->getReducedCost();
  // What results is
  double *newLower = results;
  double *objLower = newLower + numberIntegers;
  double *newUpper = objLower + numberIntegers;
  double *objUpper = newUpper + numberIntegers;
  double *interAction = objUpper + numberIntegers;
  for (i = 0; i < numberIntegers; i++) {
    int iColumn = integerVariable[i];
    newLower[i] = lower[iColumn];
    objLower[i] = -COIN_DBL_MAX;
    newUpper[i] = upper[iColumn];
    objUpper[i] = -COIN_DBL_MAX;
    interAction[i] = 0.0;
  }
  double *objMovement = new double[2 * numberIntegers];
  memset(objMovement, 0, 2 * numberIntegers * sizeof(double));
  double *saveUpper = new double[numberColumns];
  double *saveLower = new double[numberColumns];
  // Save solution in case heuristics need good solution later

  double *saveSolution = new double[numberColumns];
  memcpy(saveSolution, solver->getColSolution(), numberColumns * sizeof(double));
  model->reserveCurrentSolution(saveSolution);
  for (i = 0; i < numberColumns; i++) {
    saveLower[i] = lower[i];
    saveUpper[i] = upper[i];
  }
  // Get arrays to sort
  double *sort = new double[numberObjects];
  int *whichObject = new int[numberObjects];
  int numberToFix = 0;
  int numberToDo = 0;
  double integerTolerance = model->getDblParam(CbcModel::CbcIntegerTolerance);
  // point to useful information
  OsiBranchingInformation usefulInfo = model->usefulInformation();
  // and modify
  usefulInfo.depth_ = depth_;

  // compute current state
  int numberObjectInfeasibilities; // just odd ones
  int numberIntegerInfeasibilities;
  model->feasibleSolution(
    numberIntegerInfeasibilities,
    numberObjectInfeasibilities);
  if (solveType) {
    if ((solveType & 2) == 0)
      numberIterationsAllowed = 200 * numberIntegerInfeasibilities;
    else
      numberIterationsAllowed = COIN_INT_MAX;
  }
  int saveAllowed = numberIterationsAllowed;
#ifdef COIN_HAS_CLP
  OsiClpSolverInterface *osiclp = dynamic_cast< OsiClpSolverInterface * >(solver);
  int saveClpOptions = 0;
  bool fastIterations = (model->specialOptions() & 8) != 0;
  if (osiclp) {
    saveClpOptions = osiclp->specialOptions();
    // for faster hot start
    if (fastIterations)
      osiclp->setSpecialOptions(saveClpOptions | 8192);
    else
      osiclp->setSpecialOptions(saveClpOptions | 2048); // switch off crunch
  }
#else
  bool fastIterations = false;
#endif
  /*
    Scan for branching objects that indicate infeasibility.
    
    The algorithm is to fill the array with a set of good candidates (by
    infeasibility).
    
  */
  numberToDo = 0;
  for (i = 0; i < numberObjects; i++) {
    OsiObject *object = model->modifiableObject(i);
    CbcSimpleIntegerDynamicPseudoCost *dynamicObject = dynamic_cast< CbcSimpleIntegerDynamicPseudoCost * >(object);
    if (!dynamicObject)
      continue;
    if (dynamicObject->priority() != highestPriority)
      continue;
    double infeasibility = object->checkInfeasibility(&usefulInfo);
    int iColumn = dynamicObject->columnNumber();
    if (saveUpper[iColumn] == saveLower[iColumn])
      continue;
    if (infeasibility)
      sort[numberToDo] = -1.0e10 - infeasibility;
    else
      sort[numberToDo] = -fabs(dj[iColumn]);
    whichObject[numberToDo++] = i;
  }
  // Save basis
  CoinWarmStart *ws = solver->getWarmStart();
  int saveLimit;
  solver->getIntParam(OsiMaxNumIterationHotStart, saveLimit);
  int targetIterations = CoinMax(500, numberIterationsAllowed / numberObjects);
  if (saveLimit < targetIterations)
    solver->setIntParam(OsiMaxNumIterationHotStart, targetIterations);
  if ((solveType & 2) == 0) {
    // Mark hot start
    solver->markHotStart();
  }
  solver->setHintParam(OsiDoDualInResolve, true, OsiHintDo);
  // Sort
  CoinSort_2(sort, sort + numberToDo, whichObject);
  double *currentSolution = model->currentSolution();
  double objMin = 1.0e50;
  double objMax = -1.0e50;
  bool needResolve = false;
  int maxChoices = 1;
  int currentChoice = 0;
  int numberThreads = 0;
  bool doAtEnd = false;
  if (model->parallelMode() && (solveType & 4) != 0) {
    numberThreads = model->getNumberThreads();
    sprintf(general, "Using %d threads in analysis\n", numberThreads);
    model->messageHandler()->message(CBC_GENERAL, *model->messagesPointer())
      << general
      << CoinMessageEol;
    if (model->parallelMode() == 1) {
      maxChoices = numberThreads;
    } else {
      maxChoices = numberToDo;
      if ((solveType & 2) != 0)
        maxChoices = numberColumns;
      doAtEnd = true;
    }
  }
  StrongInfo *choices = new StrongInfo[maxChoices];
  StrongStaticInfo staticInfo;
  int numberBundles = CoinMax(1, numberThreads);
  StrongBundle *bundles = new StrongBundle[numberBundles];
  /*
    0 - available - no need to look at results
    1 - not available
    2 - available - need to look at results
  */
#ifndef NUMBER_THREADS
#define NUMBER_THREADS 4
#endif
  int status[NUMBER_THREADS];
  memset(status, 0, sizeof(status));
  memset(&staticInfo, 0, sizeof(staticInfo));
  staticInfo.solveType = solveType;
  staticInfo.originalSolution = saveSolution;
  staticInfo.back = back;
  staticInfo.ws = ws;
  staticInfo.integerTolerance = integerTolerance;
  double time1 = model->getCurrentSeconds();
#define DO_STEEPEST_SERIAL 1
#ifdef COIN_HAS_CLP
  if (osiclp && (solveType & 2) != 0 && (!numberThreads || DO_STEEPEST_SERIAL)) {
    ClpSimplex *simplex = osiclp->getModelPtr();
    simplex->setLogLevel(0);
    simplex->dual(0, 1);
    ClpDualRowPivot *dualRowPivot = simplex->dualRowPivot();
    ClpDualRowSteepest *steep = dynamic_cast< ClpDualRowSteepest * >(dualRowPivot);
    if (steep) {
      staticInfo.dualRowPivot = new ClpDualRowSteepest(*steep);
      staticInfo.dualRowPivot->setMode(1); // full steepest edge
      simplex->spareIntArray_[0] = 0;
      simplex->spareIntArray_[1] = numberRows;
      staticInfo.dualRowPivot->saveWeights(simplex, 7);
    }
  }
#endif
  for (int i = 0; i < numberBundles; i++) {
    memset(bundles + i, 0, sizeof(StrongBundle));
    bundles[i].staticInfo = &staticInfo;
  }
#if defined(THREADS_IN_ANALYZE) && defined(COIN_HAS_CLP)
#define USE_STRONG_THREADS
  CoinPthreadStuff threadInfo(numberThreads, cbc_parallelManager);
  int threadNeedsRefreshing[NUMBER_THREADS];
  for (int i = 0; i < numberThreads; i++) {
    threadInfo.threadInfo_[i].extraInfo2 = solver->clone();
    threadNeedsRefreshing[i] = 0;
  }
#ifdef COIN_HAS_CLP
  int numberSteepThreads = 0;
  int step = numberThreads ? (numberRows + numberThreads - 1) / numberThreads : 0;
  int first = 0;
  for (int i = 0; i < numberThreads; i++) {
    if (osiclp && (solveType & 2) != 0 && !DO_STEEPEST_SERIAL) {
      OsiSolverInterface *solver = reinterpret_cast< OsiSolverInterface * >(threadInfo.threadInfo_[i].extraInfo2);
      OsiClpSolverInterface *osiclp = dynamic_cast< OsiClpSolverInterface * >(solver);
      ClpSimplex *simplex = osiclp->getModelPtr();
      simplex->setLogLevel(0);
      simplex->dual(0, 1);
      ClpDualRowPivot *dualRowPivot = simplex->dualRowPivot();
      ClpDualRowSteepest *steep = dynamic_cast< ClpDualRowSteepest * >(dualRowPivot);
      if (steep) {
        numberSteepThreads = numberThreads;
        ClpDualRowSteepest *dualRowPivot = new ClpDualRowSteepest(*steep);
        dualRowPivot->setMode(1); // full steepest edge
        simplex->spareIntArray_[0] = 0;
        simplex->spareIntArray_[1] = numberRows;
        simplex->spareIntArray_[0] = first;
        simplex->spareIntArray_[1] = CoinMin(first + step, numberRows);
        first += step;
        if (i == 0)
          staticInfo.dualRowPivot = dualRowPivot;
        choices[i].steepest = dualRowPivot;
        dualRowPivot->saveWeights(simplex, 7);
      }
    }
  }
  if (numberSteepThreads && false) {
    int numberDone = 0;
    int iDo = 0;
    staticInfo.solveType = 200;
    while (numberDone < numberSteepThreads) {
      int iThread;
      threadInfo.waitParallelTask(1, iThread, iDo < numberToDo);
      int threadStatus = 1 + threadInfo.threadInfo_[iThread].status;
      iThread = iThread;
      if (threadStatus == 0 && iDo < numberSteepThreads) {
        StrongInfo &choice = choices[iThread];
        StrongBundle &bundle = bundles[iThread];
        bundle.whichChoice = iThread;
        memset(&choice, 0, sizeof(StrongInfo));
        iDo++; //started this one
        bundle.choice = &choice;
        bundle.solver = solver;
        bundle.solver = reinterpret_cast< OsiSolverInterface * >(threadInfo.threadInfo_[iThread].extraInfo2);
        threadStatus = 0;
#ifdef DETAIL_THREAD
        printf("Starting steep task on thread %d\n",
          choice.iThread);
#endif
        threadInfo.startParallelTask(1, iThread, &bundle);
      }
      if (!threadStatus) {
#ifdef _MSC_VER
        Sleep(1);
#else
        usleep(1000);
#endif
        continue;
      }
      if (threadStatus) {
        numberDone++;
        // say available
        threadInfo.sayIdle(iThread);
      }
      staticInfo.solveType = solveType;
    }
    OsiSolverInterface *solver0 = reinterpret_cast< OsiSolverInterface * >(threadInfo.threadInfo_[0].extraInfo2);
    CoinIndexedVector *savedWeights0 = staticInfo.dualRowPivot->savedWeights();
    int *index0 = savedWeights0->getIndices();
    double *weight0 = savedWeights0->denseVector();
    int step = (numberRows + numberSteepThreads - 1) / numberSteepThreads;
    int first = step;
    //memset(weight0+first,0,(numberRows-first)*sizeof(double));
    for (int i = 1; i < numberSteepThreads; i++) {
      int n = CoinMin(step, numberRows - first);
      CoinIndexedVector *savedWeights = choices[i].steepest->savedWeights();
      int *index = savedWeights->getIndices();
      double *weight = savedWeights->denseVector();
      memcpy(index0 + first, index + first, n * sizeof(int));
      memcpy(weight0 + first, weight + first, n * sizeof(double));
      first += step;
      delete choices[i].steepest;
      choices[i].steepest = NULL;
    }
    //for (int j=0;j<numberRows;j++)
    //weight0[j]=1.0;
  }
#endif
#endif
  double bestSolutionValue = model->getMinimizationObjValue();
  double *bestSolution = NULL;
  double cutoff;
  solver->getDblParam(OsiDualObjectiveLimit, cutoff);
  double maxMovement = 2.0 * (cutoff - objectiveValue_) + 1.0e-6;
  /*
    Now calculate the cost forcing the variable up and down.
  */
  int iDo = 0;
  int iDone = -1;
  int numberDone = 0;
  int iThread = 0;
  int threadStatus = 0;
  int whenPrint = (numberToDo + 9) / 10;
  while (numberDone < numberToDo) {
    if ((solveType & 32) != 0 && numberDone == whenPrint) {
      whenPrint += (numberToDo + 9) / 10;
      sprintf(general, "%d variables looked at - %d changed by initial strong branching (%.2f seconds - %d iterations)", numberDone, numberToFix, model->getCurrentSeconds() - time1, saveAllowed - numberIterationsAllowed);
      model->messageHandler()->message(CBC_GENERAL, *model->messagesPointer())
        << general
        << CoinMessageEol;
    }
#ifdef USE_STRONG_THREADS
    if (numberThreads) {
      threadInfo.waitParallelTask(1, iThread, iDo < numberToDo);
      threadStatus = 1 + threadInfo.threadInfo_[iThread].status;
      if (!doAtEnd)
        currentChoice = iThread;
      if (threadNeedsRefreshing[iThread]) {
        OsiSolverInterface *solver = reinterpret_cast< OsiSolverInterface * >(threadInfo.threadInfo_[iThread].extraInfo2);
        if ((threadNeedsRefreshing[iThread] & 1) != 0)
          solver->setColLower(saveLower);
        if ((threadNeedsRefreshing[iThread] & 2) != 0)
          solver->setColUpper(saveUpper);
        threadNeedsRefreshing[iThread] = 0;
      }
    }
#endif
    if (threadStatus == 0 && iDo < numberToDo) {
      StrongInfo &choice = choices[currentChoice];
      StrongBundle &bundle = bundles[iThread];
      bundle.whichChoice = currentChoice;
      memset(&choice, 0, sizeof(StrongInfo));
      int iObject = whichObject[iDo];
      iDo++; //started this one
      OsiObject *object = model->modifiableObject(iObject);
      CbcSimpleIntegerDynamicPseudoCost *dynamicObject = dynamic_cast< CbcSimpleIntegerDynamicPseudoCost * >(object);
      int iColumn = dynamicObject->columnNumber();
      double value = currentSolution[iColumn];
      double nearest = floor(value + 0.5);
      double lowerValue = floor(value);
      bool satisfied = false;
      if (fabs(value - nearest) <= integerTolerance || value < saveLower[iColumn] || value > saveUpper[iColumn]) {
        satisfied = true;
        if (nearest < saveUpper[iColumn]) {
          lowerValue = nearest;
        } else {
          lowerValue = nearest - 1;
        }
      }
      double upperValue = lowerValue + 1.0;
      // Save which object it was
      choice.columnNumber = iColumn;
      choice.initialValue = value;
      choice.upLowerBound = upperValue;
      choice.downUpperBound = lowerValue;
      choice.numIntInfeas[1] = numberUnsatisfied_;
      choice.numIntInfeas[0] = numberUnsatisfied_;
      choice.movement[0] = 0.0;
      choice.movement[1] = 0.0;
      choice.numIters[0] = 0;
      choice.numIters[1] = 0;
      if (fabs(value - lowerValue) <= integerTolerance)
        choice.numIters[0] = -1; // mark as not done
      if (fabs(value - upperValue) <= integerTolerance)
        choice.numIters[1] = -1; // mark as not done
      bundle.choice = &choice;
      bundle.solver = solver;
#ifdef USE_STRONG_THREADS
      if (numberThreads) {
        bundle.solver = reinterpret_cast< OsiSolverInterface * >(threadInfo.threadInfo_[iThread].extraInfo2);
        threadStatus = 0;
#ifdef DETAIL_THREAD
        printf("Starting task for column %d on thread %d\n",
          choice.columnNumber, iThread);
#endif
        threadInfo.startParallelTask(1, iThread, &bundle);
      } else {
#endif
        threadStatus = 2;
        solveAnalyze(&bundle);
#ifdef USE_STRONG_THREADS
      }
#endif
    }
    if (!threadStatus) {
#ifdef _MSC_VER
      Sleep(1);
#else
      usleep(1000);
#endif
      continue;
    }
    if (threadStatus) {
      int whichChoice = bundles[iThread].whichChoice;
      StrongInfo &choice = choices[whichChoice];
      int iColumn = choice.columnNumber;
      if (choice.integerSolution) {
        double *foundSolution = choice.integerSolution;
        solver->setColSolution(foundSolution);
        // See if integer solution
        int numberInfeas = 0;
        int numberOddInfeas = 0;
        if (model->feasibleSolution(numberInfeas, numberOddInfeas)
          && model->problemFeasibility()->feasible(model, -1) >= 0) {
          double newObjectiveValue;
          solver->getDblParam(OsiObjOffset, newObjectiveValue);
          newObjectiveValue = -newObjectiveValue;
          const double *cost = solver->getObjCoefficients();
          for (int i = 0; i < numberColumns; i++)
            newObjectiveValue += cost[i] * foundSolution[i];
          if (newObjectiveValue < bestSolutionValue) {
            if (doAtEnd) {
              if (!bestSolution)
                bestSolution = CoinCopyOfArray(foundSolution, numberColumns);
              else
                memcpy(bestSolution, foundSolution, numberColumns * sizeof(double));
              bestSolutionValue = newObjectiveValue;
            } else {
              model->setBestSolution(CBC_STRONGSOL,
                newObjectiveValue,
                foundSolution);
              model->setLastHeuristic(NULL);
              model->incrementUsed(solver->getColSolution());
              bestSolutionValue = model->getMinimizationObjValue();
            }
          }
        }
        delete[] foundSolution;
      }
      for (int iWay = 0; iWay < 2; iWay++) {
        numberIterationsAllowed -= choice.numIters[iWay];
        choice.movement[iWay] -= objectiveValue_;
      }
      // If objective goes above certain amount we can set bound
      int jInt = back[iColumn];
      OsiObject *object = model->modifiableObject(jInt);
      CbcSimpleIntegerDynamicPseudoCost *dynamicObject = dynamic_cast< CbcSimpleIntegerDynamicPseudoCost * >(object);
      if (dynamicObject) {
        if (choice.numIters[0] >= 0) {
          dynamicObject->addToSumDownCost(CoinMin(choice.movement[0], maxMovement));
          dynamicObject->addToSumDownChange(choice.initialValue - choice.downUpperBound);
        }
        if (choice.numIters[1] >= 0) {
          dynamicObject->addToSumUpCost(CoinMin(choice.movement[1], maxMovement));
          dynamicObject->addToSumUpChange(choice.upLowerBound - choice.initialValue);
        }
      }
      newLower[jInt] = choice.upLowerBound;
      if (choice.finished[0])
        objLower[jInt] = choice.movement[0] + objectiveValue_;
      else
        objLower[jInt] = objectiveValue_;
      newUpper[jInt] = choice.downUpperBound;
      if (choice.finished[1])
        objUpper[jInt] = choice.movement[1] + objectiveValue_;
      else
        objUpper[jInt] = objectiveValue_;
      objMin = CoinMin(CoinMin(objLower[jInt], objUpper[jInt]), objMin);
      objMovement[2 * jInt] = choice.movement[0];
      objMovement[2 * jInt + 1] = choice.movement[1];
      double sumModified = choice.modified[0] + choice.modified[1] + 1.0e-15 * (choice.sumModified[0] + choice.sumModified[1]);
      if (choice.numIters[0] >= 0 && choice.numIters[1] >= 0)
        sumModified *= 0.6;
      interAction[jInt] = sumModified;
      /*
	End of evaluation for this candidate variable. Possibilities are:
	* Both sides below cutoff; this variable is a candidate for branching.
	* Both sides infeasible or above the objective cutoff: no further action
	here. Break from the evaluation loop and assume the node will be purged
	by the caller.
	* One side below cutoff: Install the branch (i.e., fix the variable). Break
	from the evaluation loop and assume the node will be reoptimised by the
	caller.
      */
      threadStatus = 0;
      currentChoice++;
      numberDone++;
#ifdef USE_STRONG_THREADS
      // say available
      if (numberThreads) {
        threadInfo.sayIdle(iThread);
      }
#endif
      if (doAtEnd)
        continue;
      if (choice.movement[1] < 1.0e100) {
        if (choice.movement[0] < 1.0e100) {
          objMax = CoinMax(CoinMax(objLower[jInt], objUpper[jInt]), objMax);
          // In case solution coming in was odd
          choice.movement[1] = CoinMax(0.0, choice.movement[1]);
          choice.movement[0] = CoinMax(0.0, choice.movement[0]);
          // feasible -
          model->messageHandler()->message(CBC_STRONG, *model->messagesPointer())
            << iColumn << iColumn
            << choice.movement[0] << choice.numIntInfeas[0]
            << choice.movement[1] << choice.numIntInfeas[1]
            << choice.initialValue
            << CoinMessageEol;
        } else {
          // up feasible, down infeasible
          needResolve = true;
          numberToFix++;
          saveLower[iColumn] = choice.upLowerBound;
          solver->setColLower(iColumn, choice.upLowerBound);
#ifdef USE_STRONG_THREADS
          for (int i = 0; i < numberThreads; i++) {
            threadNeedsRefreshing[i] |= 1;
          }
#endif
        }
      } else {
        if (choice.movement[0] < 1.0e100) {
          // down feasible, up infeasible
          needResolve = true;
          numberToFix++;
          saveUpper[iColumn] = choice.downUpperBound;
          solver->setColUpper(iColumn, choice.downUpperBound);
#ifdef USE_STRONG_THREADS
          for (int i = 0; i < numberThreads; i++) {
            threadNeedsRefreshing[i] |= 2;
          }
#endif
        } else {
          // neither side feasible
          COIN_DETAIL_PRINT(printf("Both infeasible for choice %d sequence %d\n", i,
            model->object(choice.objectNumber)->columnNumber()));
          //solver->writeMps("bad");
          numberToFix = -1;
          break;
        }
      }
      if (numberIterationsAllowed <= 0)
        break;
      if (currentChoice == maxChoices)
        currentChoice = 0;
    }
    //printf("obj %d, col %d, down %g up %g value %g\n",iObject,iColumn,
    //     choice.downMovement,choice.upMovement,value);
  }
  // Do at end if deterministic
  if (doAtEnd) {
    if (bestSolution) {
      model->setBestSolution(CBC_STRONGSOL,
        bestSolutionValue,
        bestSolution);
      model->setLastHeuristic(NULL);
      model->incrementUsed(solver->getColSolution());
      delete[] bestSolution;
    }
    for (int iDo = 0; iDo < numberLookIntegers; iDo++) {
      StrongInfo &choice = choices[iDo];
      int iColumn = choice.columnNumber;
      int iObject = iColumn;
      int jInt = back[iColumn];
      double value = choice.initialValue;
      double lowerValue = choice.downUpperBound;
      double upperValue = choice.upLowerBound;
      if (choice.movement[1] < 1.0e100) {
        if (choice.movement[0] < 1.0e100) {
          objMax = CoinMax(CoinMax(objLower[jInt], objUpper[jInt]), objMax);
          // In case solution coming in was odd
          choice.movement[1] = CoinMax(0.0, choice.movement[1]);
          choice.movement[0] = CoinMax(0.0, choice.movement[0]);
          // feasible -
          model->messageHandler()->message(CBC_STRONG, *model->messagesPointer())
            << iObject << iColumn
            << choice.movement[0] << choice.numIntInfeas[0]
            << choice.movement[1] << choice.numIntInfeas[1]
            << value
            << CoinMessageEol;
        } else {
          // up feasible, down infeasible
          numberToFix++;
          saveLower[iColumn] = upperValue;
          solver->setColLower(iColumn, upperValue);
        }
      } else {
        if (choice.movement[0] < 1.0e100) {
          // down feasible, up infeasible
          needResolve = true;
          numberToFix++;
          saveUpper[iColumn] = lowerValue;
          solver->setColUpper(iColumn, lowerValue);
        } else {
          // neither side feasible
          COIN_DETAIL_PRINT(printf("Both infeasible for choice %d sequence %d\n", i,
            model->object(choice.objectNumber)->columnNumber()));
          //solver->writeMps("bad");
          numberToFix = -1;
          break;
        }
      }
    }
  }
  if (false) {
    const double *lower = solver->getColLower();
    const double *upper = solver->getColUpper();
    for (int i = 0; i < numberColumns; i++) {
      if (lower[i] != saveLower[i] || upper[i] != saveUpper[i]) {
        printf("%d changed- saved %g,%g now %g,%g\n", i, saveLower[i], saveUpper[i],
          lower[i], upper[i]);
      }
    }
  }
  COIN_DETAIL_PRINT(printf("Best possible solution %g, can fix more if solution of %g found - looked at %d variables in %d iterations\n",
    objMin, objMax, iDo, model->numberAnalyzeIterations() - numberIterationsAllowed));
  model->setNumberAnalyzeIterations(numberIterationsAllowed);
  if (numberToFix > 0) {
    sprintf(general, "%d variable bounds modified by initial strong branching (%.2f seconds - %d iterations)", numberToFix, model->getCurrentSeconds() - time1, saveAllowed - numberIterationsAllowed);
  } else if (numberToFix < 0) {
    sprintf(general, "initial strong branching found to be infeasible (%.2f seconds - %d iterations)", model->getCurrentSeconds() - time1, saveAllowed - numberIterationsAllowed);
  } else if ((solveType & 32) != 0) {
    sprintf(general, "No variables fixed by initial strong branching (%.2f seconds - %d iterations)", model->getCurrentSeconds() - time1, saveAllowed - numberIterationsAllowed);
  } else {
    general[0] = '\0';
  }
  if (general[0] != '\0')
    model->messageHandler()->message(CBC_GENERAL, *model->messagesPointer())
      << general
      << CoinMessageEol;
  double smallestEffect = COIN_DBL_MAX;
  double largestEffect = 0.0;
  for (i = 0; i < numberIntegers; i++) {
    int iColumn = integerVariable[i];
    if (back[iColumn] >= numberColumns)
      continue;
    smallestEffect = CoinMin(smallestEffect, interAction[i]);
    largestEffect = CoinMax(largestEffect, interAction[i]);
  }
  double groupValue[11];
  int groupCounts[11] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
  groupValue[10] = largestEffect;
  for (int i = 0; i < 10; i++)
    groupValue[i] = smallestEffect + i * 0.1 * (largestEffect - smallestEffect);
  sprintf(general, "Looked at %d integer variables - smallest interaction %g",
    numberLookIntegers, smallestEffect);
  model->messageHandler()->message((solveType & 32) == 0 ? CBC_FPUMP2 : CBC_FPUMP1,
    *model->messagesPointer())
    << general
    << CoinMessageEol;
  for (int i = 0; i < numberIntegers; i++) {
    int iColumn = integerVariable[i];
    if (back[iColumn] >= numberColumns)
      continue;
    double value = interAction[i];
    int j;
    for (j = 0; j < 11; j++) {
      if (value <= groupValue[j] || j == 10)
        break;
    }
    groupCounts[j]++;
  }
  general[0] = '\0';
  for (int i = 0; i < 11; i++)
    sprintf(general + strlen(general), "%d <= %g ", groupCounts[i], groupValue[i]);
  model->messageHandler()->message((solveType & 32) == 0 ? CBC_FPUMP2 : CBC_FPUMP1,
    *model->messagesPointer())
    << general
    << CoinMessageEol;
  smallestEffect = COIN_DBL_MAX;
  largestEffect = 0.0;
  int numberChanged = 0;
  int numberZeroMoved = 0;
  for (i = 0; i < numberIntegers; i++) {
    int iColumn = integerVariable[i];
    if (back[iColumn] >= numberColumns)
      continue;
    for (int iWay = 0; iWay < 2; iWay++) {
      double value = objMovement[2 * i + iWay];
      if (value < 1.0e-7) {
        numberZeroMoved++;
      } else if (value < 1.0e50) {
        smallestEffect = CoinMin(smallestEffect, value);
        largestEffect = CoinMax(largestEffect, value);
      } else {
        numberChanged++;
      }
    }
  }
  memset(groupCounts, 0, sizeof(groupCounts));
  groupValue[10] = largestEffect;
  for (int i = 0; i < 10; i++)
    groupValue[i] = smallestEffect + i * 0.1 * (largestEffect - smallestEffect);
  sprintf(general, "Strong branching - %d bounds changed, %d zero objective changes and %d nonzero (smallest %g)",
    numberChanged, numberZeroMoved,
    2 * numberLookIntegers - numberChanged - numberZeroMoved, smallestEffect);
  model->messageHandler()->message((solveType & 32) == 0 ? CBC_FPUMP2 : CBC_FPUMP1,
    *model->messagesPointer())
    << general
    << CoinMessageEol;
  sprintf(general, "Breakdown ");
  for (i = 0; i < numberIntegers; i++) {
    int iColumn = integerVariable[i];
    if (back[iColumn] >= numberColumns)
      continue;
    for (int iWay = 0; iWay < 2; iWay++) {
      double value = objMovement[2 * i + iWay];
      int j;
      for (j = 0; j < 11; j++) {
        if (value <= groupValue[j] || j == 10)
          break;
      }
      groupCounts[j]++;
    }
  }
  for (int i = 0; i < 11; i++)
    sprintf(general + strlen(general), "%d <= %g ", groupCounts[i], groupValue[i]);
  model->messageHandler()->message((solveType & 32) == 0 ? CBC_FPUMP2 : CBC_FPUMP1,
    *model->messagesPointer())
    << general
    << CoinMessageEol;
  delete[] objMovement;
  if ((solveType & 2) == 0) {
    // Delete the snapshot
    solver->unmarkHotStart();
  }
  // back to normal
  solver->setHintParam(OsiDoInBranchAndCut, true, OsiHintDo, NULL);
  solver->setIntParam(OsiMaxNumIterationHotStart, saveLimit);
  // restore basis
  solver->setWarmStart(ws);
  // skip if infeasible
  if (numberToFix < 0)
    solveType = 0;
  int numberBoundsChanged = 0;
  if ((solveType & 16) == 0) {
    double size = numberRows;
    size *= numberColumns;
    if (size > 1000000) {
      if ((solveType & 32) != 0)
        model->messageHandler()->message(CBC_GENERAL, *model->messagesPointer())
          << "Skipping analyze on other columns as problem too large"
          << CoinMessageEol;
      solveType &= ~2;
    }
  }
  if ((solveType & 2) != 0) {
#ifdef COIN_HAS_CLP
    int saveOptions = osiclp ? osiclp->specialOptions() : 0;
    if (osiclp) {
      //ClpPrimalColumnPivot * primalColumnPivot=NULL;
      osiclp->setSpecialOptions(saveOptions | 2048); // off crunch
    }
#endif
    double *newLower = new double[2 * numberColumns];
    double *newUpper = newLower + numberColumns;
    // look at ints/all - should be parametrics - for now primal
    OsiSolverInterface *temp = solver->clone();
    // add constraint
    int *indices = reinterpret_cast< int * >(newUpper);
    double *obj = newLower;
    memcpy(obj, solver->getObjCoefficients(), numberColumns * sizeof(double));
    int n = 0;
    for (int i = 0; i < numberColumns; i++) {
      if (obj[i]) {
        indices[n] = i;
        obj[n++] = obj[i];
      }
    }
    if (n) {
      double cutoff = model->getCutoff();
      // relax a little bit
      cutoff += 1.0e-4;
      double offset;
      temp->getDblParam(OsiObjOffset, offset);
      temp->addRow(n, indices, obj, -COIN_DBL_MAX, CoinMin(cutoff, 1.0e25) + offset);
      temp->setDblParam(OsiObjOffset, 0.0);
#if defined(THREADS_IN_ANALYZE) && defined(COIN_HAS_CLP)
      for (int iThread = 0; iThread < numberThreads; iThread++) {
        OsiSolverInterface *solver = reinterpret_cast< OsiSolverInterface * >(threadInfo.threadInfo_[iThread].extraInfo2);
        solver->addRow(n, indices, obj, -COIN_DBL_MAX, CoinMin(cutoff, 1.0e25) + offset);
      }
#endif
    }
    //temp->setHintParam(OsiDoDualInResolve, false, OsiHintDo) ;
    temp->setHintParam(OsiDoReducePrint, true, OsiHintTry);
    temp->setDblParam(OsiDualObjectiveLimit, COIN_DBL_MAX);
    temp->resolve();
    {
      const double *lower = temp->getColLower();
      const double *upper = temp->getColUpper();
      for (int i = 0; i < numberColumns; i++) {
        assert(lower[i] == saveLower[i]);
        assert(upper[i] == saveUpper[i]);
      }
    }
    delete ws;
    ws = temp->getWarmStart();
    staticInfo.ws = ws;
    staticInfo.newObjective = new double[2 * numberColumns];
    memcpy(staticInfo.newObjective, solver->getObjCoefficients(), numberColumns * sizeof(double));
    memset(staticInfo.newObjective + numberColumns, 0, numberColumns * sizeof(double));
#if defined(THREADS_IN_ANALYZE) && defined(COIN_HAS_CLP)
    for (int iThread = 0; iThread < numberThreads; iThread++) {
      OsiSolverInterface *solver = reinterpret_cast< OsiSolverInterface * >(threadInfo.threadInfo_[iThread].extraInfo2);
      solver->setObjective(newLower);
      solver->setDblParam(OsiDualObjectiveLimit, COIN_DBL_MAX);
      threadNeedsRefreshing[iThread] = 3;
    }
#endif
    for (int i = 0; i < numberColumns; i++) {
      newLower[i] = lower[i];
      newUpper[i] = upper[i];
    }
    double *thisSolution = CoinCopyOfArray(temp->getColSolution(), numberColumns);
    double primalTolerance;
    solver->getDblParam(OsiPrimalTolerance, primalTolerance);
    iDo = 0;
    iDone = -1;
    numberDone = 0;
    int iThread = 0;
    threadStatus = 0;
    currentChoice = 0;
    staticInfo.solveType = 100; //mark for analyze
    staticInfo.originalSolution = thisSolution;
    whenPrint = (numberColumns + 9) / 10;
    while (numberDone < numberColumns) {
      if ((solveType & 32) != 0 && numberDone == whenPrint) {
        whenPrint += (numberColumns + 9) / 10;
        sprintf(general, "%d variables looked at - %d bounds changed by secondary solves (%.2f seconds - %d iterations)", numberDone, numberBoundsChanged, model->getCurrentSeconds() - time1, saveAllowed - numberIterationsAllowed);
        model->messageHandler()->message(CBC_GENERAL, *model->messagesPointer())
          << general
          << CoinMessageEol;
      }
#ifdef USE_STRONG_THREADS
      if (numberThreads) {
        threadInfo.waitParallelTask(1, iThread, iDo < numberColumns);
        threadStatus = 1 + threadInfo.threadInfo_[iThread].status;
        currentChoice = iThread;
        if (threadNeedsRefreshing[iThread]) {
          OsiSolverInterface *solver = reinterpret_cast< OsiSolverInterface * >(threadInfo.threadInfo_[iThread].extraInfo2);
          if ((threadNeedsRefreshing[iThread] & 1) != 0)
            solver->setColLower(saveLower);
          if ((threadNeedsRefreshing[iThread] & 2) != 0)
            solver->setColUpper(saveUpper);
          threadNeedsRefreshing[iThread] = 0;
        }
      }
#endif
      if (threadStatus == 0 && iDo < numberColumns) {
        StrongInfo &choice = choices[currentChoice];
        StrongBundle &bundle = bundles[currentChoice];
        bundle.whichChoice = currentChoice;
        memset(&choice, 0, sizeof(StrongInfo));
        int iColumn = iDo;
        iDo++; //started this one
        int typeSolve = 0;
        if (thisSolution[iColumn] > newLower[iColumn] + integerTolerance)
          typeSolve = 1;
        if (thisSolution[iColumn] < newUpper[iColumn] - integerTolerance)
          typeSolve += 2;
        if (typeSolve && back[iColumn] >= 0) {
          if (thisSolution[iColumn] < newLower[iColumn] + 0.9999)
            typeSolve &= ~1;
          if (thisSolution[iColumn] > newUpper[iColumn] - 0.9999)
            typeSolve &= ~2;
          if (temp->isBinary(iColumn))
            typeSolve = 0; // already done
        }
        if (typeSolve == 0 || newUpper[iColumn] == newLower[iColumn]) {
#ifdef USE_STRONG_THREADS
          // say available
          if (numberThreads) {
            threadInfo.sayIdle(iThread);
          }
#endif
          numberDone++;
          continue;
        }
        // Save which object it was
        choice.columnNumber = iColumn;
        choice.initialValue = thisSolution[iColumn];
        choice.movement[0] = COIN_DBL_MAX;
        choice.movement[1] = -COIN_DBL_MAX;
        choice.upLowerBound = newUpper[iColumn];
        choice.downUpperBound = newLower[iColumn];
        if ((typeSolve & 1) == 0)
          choice.numIters[0] = -1; // mark as not done
        if ((typeSolve & 2) == 0)
          choice.numIters[1] = -1; // mark as not done
        bundle.choice = &choice;
        bundle.solver = temp;
#ifdef USE_STRONG_THREADS
        if (numberThreads) {
          bundle.solver = reinterpret_cast< OsiSolverInterface * >(threadInfo.threadInfo_[iThread].extraInfo2);
          threadStatus = 0;
#ifdef DETAIL_THREAD
          printf("Starting task for column %d on thread %d\n",
            choice.columnNumber, iThread);
#endif
          threadInfo.startParallelTask(1, iThread, &bundle);
        } else {
#endif
          threadStatus = 2;
          solveAnalyze(&bundle);
#ifdef USE_STRONG_THREADS
        }
#endif
      }
      if (threadStatus) {
        int whichChoice = bundles[iThread].whichChoice;
        StrongInfo &choice = choices[whichChoice];
        int iColumn = choice.columnNumber;
        if (choice.modified[0]) {
          numberToFix = -numberColumns - 1;
        }
        double gotLower = COIN_DBL_MAX;
        double gotUpper = -COIN_DBL_MAX;
        if (choice.numIters[0] >= 0) {
          // go down
          double value = choice.movement[0];
          if (value > newLower[iColumn] + 100.0 * integerTolerance) {
            if (back[iColumn] >= 0)
              value = ceil(value - integerTolerance);
            else
              value = CoinMax(newLower[iColumn], value - 1.0e-5 - 1.0e-8 * fabs(value));
            if (value > newLower[iColumn] + 1.0e-8 * (1.0 + fabs(value))) {
              sprintf(general, "Secondary analysis solve increases lower bound on %d from %g to %g%s",
                iColumn, newUpper[iColumn], value, (back[iColumn] >= 0) ? "(integer)" : "");
              model->messageHandler()->message(CBC_FPUMP2, *model->messagesPointer())
                << general
                << CoinMessageEol;
              numberBoundsChanged++;
              if (value > newUpper[iColumn] - primalTolerance) {
                value = newUpper[iColumn];
                if (value > newUpper[iColumn] + 10.0 * primalTolerance) {
                  // infeasible
                  numberToFix = -numberColumns - 1;
                }
              }
              gotLower = value;
            }
          }
        }
        if (choice.numIters[1] >= 0) {
          // go up
          double value = choice.movement[1];
          if (value < newUpper[iColumn] - 100.0 * integerTolerance) {
            if (back[iColumn] >= 0)
              value = floor(value + integerTolerance);
            else
              value = CoinMin(newUpper[iColumn], value + 1.0e-5 + 1.0e-8 * fabs(value));
            if (value < newUpper[iColumn] - 1.0e-8 * (1.0 + fabs(value))) {
              sprintf(general, "Secondary analysis solve decreases upper bound on %d from %g to %g%s",
                iColumn, newUpper[iColumn], value, (back[iColumn] >= 0) ? "(integer)" : "");
              model->messageHandler()->message(CBC_FPUMP2, *model->messagesPointer())
                << general
                << CoinMessageEol;
              numberBoundsChanged++;
              if (value < newLower[iColumn] + primalTolerance) {
                value = newLower[iColumn];
                if (value < newLower[iColumn] - 10.0 * primalTolerance) {
                  // infeasible
                  numberToFix = -numberColumns - 1;
                }
              }
              gotUpper = value;
            }
          }
        }
        if (gotLower != COIN_DBL_MAX) {
          newLower[iColumn] = gotLower;
          temp->setColLower(iColumn, gotLower);
          if (!doAtEnd)
            solver->setColLower(iColumn, gotLower);
        }
        if (gotUpper != -COIN_DBL_MAX) {
          gotUpper = CoinMax(gotUpper, newLower[iColumn]);
          newUpper[iColumn] = gotUpper;
          temp->setColUpper(iColumn, gotUpper);
          if (!doAtEnd)
            solver->setColUpper(iColumn, gotUpper);
        }
#if 0
	if ((model->specialOptions()&1) != 0) {
	  const OsiRowCutDebugger *debugger = solver->getRowCutDebugger() ;
	  if (!debugger) {
	    abort();
	  } else {
	    printf("still ok\n");
	  }
	}
#endif
        threadStatus = 0;
        currentChoice++;
        numberDone++;
        for (int iWay = 0; iWay < 2; iWay++) {
          if (choice.numIters[iWay] > 0)
            numberIterationsAllowed -= choice.numIters[iWay];
        }
        if (currentChoice == maxChoices)
          currentChoice = 0;
#ifdef USE_STRONG_THREADS
        // say available
        if (numberThreads) {
          threadInfo.sayIdle(iThread);
        }
#endif
      }
    }
    delete[] thisSolution;
    delete temp;
    delete[] newLower;
#ifdef COIN_HAS_CLP
    if (osiclp) {
      //ClpPrimalColumnPivot * primalColumnPivot=NULL;
      osiclp->setSpecialOptions(saveOptions);
    }
#endif
  }
  delete[] staticInfo.newObjective;
#ifdef COIN_HAS_CLP
  if (osiclp) {
    delete staticInfo.dualRowPivot;
    delete staticInfo.primalColumnPivot;
    ClpSimplex *simplex = osiclp->getModelPtr();
    ClpDualRowPivot *dualRowPivot = simplex->dualRowPivot();
    ClpDualRowSteepest *steep = dynamic_cast< ClpDualRowSteepest * >(dualRowPivot);
    if (steep)
      steep->setMode(3);
  }
#endif
  if ((solveType & 64) != 0) {
    OsiSolverInterface *temp = solver->clone();
    int numberRows = solver->getNumRows();
    int numberContinuousRows = model->numberRowsAtContinuous();
    int *del = new int[numberRows - numberContinuousRows];
    for (int i = numberContinuousRows; i < numberRows; i++)
      del[i - numberContinuousRows] = i;
    temp->deleteRows(numberRows - numberContinuousRows, del);
    delete[] del;
#ifdef COIN_HAS_CLP
    if (!osiclp) {
#endif
      solver->writeMps("analyzed");
      temp->writeMps("analyzed2");
#ifdef COIN_HAS_CLP
    } else {
      OsiClpSolverInterface *osiclp2 = dynamic_cast< OsiClpSolverInterface * >(temp);
      osiclp->getModelPtr()->writeMps("analyzed.mps", 2, 1);
      osiclp2->getModelPtr()->writeMps("analyzed2.mps", 2, 1);
    }
#endif
    delete temp;
    model->messageHandler()->message(CBC_GENERAL, *model->messagesPointer())
      << "Models saved on 'analyzed' and 'analyzed2'"
      << CoinMessageEol;
  }
  delete[] choices;
  for (int i = 0; i < numberBundles; i++) {
    delete[] bundles[i].tempSolution;
    delete bundles[i].tempBasis;
  }
  delete[] bundles;
#ifdef USE_STRONG_THREADS
  if (numberThreads) {
    threadInfo.waitAllTasks();
    for (int i = 0; i < numberThreads; i++) {
      delete reinterpret_cast< OsiSolverInterface * >(threadInfo.threadInfo_[i].extraInfo2);
    }
  }
#endif
  delete ws;

  delete[] sort;
  delete[] whichObject;
  delete[] saveLower;
  delete[] saveUpper;
  delete[] back;
  // restore solution
  solver->setColSolution(saveSolution);
#ifdef COIN_HAS_CLP
  if (osiclp)
    osiclp->setSpecialOptions(saveClpOptions);
#endif
  delete[] saveSolution;
  solver->resolve();
  if (numberToFix < 0 && !solver->isProvenOptimal()) {
    // infeasible
    model->messageHandler()->message(CBC_GENERAL, *model->messagesPointer())
      << "Analysis shows problem to be infeasible"
      << CoinMessageEol;
    return numberToFix;
  }
  if (numberBoundsChanged) {
    sprintf(general, "%d bounds changed by secondary solves (%.2f seconds - %d iterations)",
      numberBoundsChanged, model->getCurrentSeconds() - time1, saveAllowed - numberIterationsAllowed);
    model->messageHandler()->message(CBC_GENERAL, *model->messagesPointer())
      << general
      << CoinMessageEol;
  } else if ((solveType & 32) != 0) {
    sprintf(general, "No bounds changed by secondary solves (%.2f seconds - %d iterations)",
      model->getCurrentSeconds() - time1, saveAllowed - numberIterationsAllowed);
    model->messageHandler()->message(CBC_GENERAL, *model->messagesPointer())
      << general
      << CoinMessageEol;
  }
  model->reserveCurrentSolution(solver->getColSolution());
  if ((solveType & 1) != 0) {
    if (groupCounts[0] * 4 > numberIntegers) {
      // change priority on this group
      int generalPriority = -10000000;
      for (int i = 0; i < numberIntegers; i++) {
        OsiObject *object = model->modifiableObject(i);
        CbcSimpleInteger *integerObject = dynamic_cast< CbcSimpleInteger * >(object);
        if (!integerObject)
          continue;
        generalPriority = CoinMax(generalPriority, integerObject->priority());
      }
      for (int i = 0; i < numberIntegers; i++) {
        OsiObject *object = model->modifiableObject(i);
        CbcSimpleInteger *integerObject = dynamic_cast< CbcSimpleInteger * >(object);
        if (!integerObject)
          continue;
        if (!interAction[i] && integerObject->priority() == generalPriority)
          integerObject->setPriority(generalPriority + 1);
      }
    }
  }
  return numberToFix;
}

CbcNode::CbcNode(const CbcNode &rhs)
  : CoinTreeNode(rhs)
{
#ifdef CHECK_NODE
  printf("CbcNode %p Constructor from rhs %p\n", this, &rhs);
#endif
  if (rhs.nodeInfo_)
    nodeInfo_ = rhs.nodeInfo_->clone();
  else
    nodeInfo_ = NULL;
  objectiveValue_ = rhs.objectiveValue_;
  guessedObjectiveValue_ = rhs.guessedObjectiveValue_;
  sumInfeasibilities_ = rhs.sumInfeasibilities_;
  if (rhs.branch_)
    branch_ = rhs.branch_->clone();
  else
    branch_ = NULL;
  depth_ = rhs.depth_;
  numberUnsatisfied_ = rhs.numberUnsatisfied_;
  nodeNumber_ = rhs.nodeNumber_;
  state_ = rhs.state_;
  if (nodeInfo_)
    assert((state_ & 2) != 0);
  else
    assert((state_ & 2) == 0);
}

CbcNode &
CbcNode::operator=(const CbcNode &rhs)
{
  if (this != &rhs) {
    delete nodeInfo_;
    if (rhs.nodeInfo_)
      nodeInfo_ = rhs.nodeInfo_->clone();
    else
      nodeInfo_ = NULL;
    objectiveValue_ = rhs.objectiveValue_;
    guessedObjectiveValue_ = rhs.guessedObjectiveValue_;
    sumInfeasibilities_ = rhs.sumInfeasibilities_;
    if (rhs.branch_)
      branch_ = rhs.branch_->clone();
    else
      branch_ = NULL,
      depth_ = rhs.depth_;
    numberUnsatisfied_ = rhs.numberUnsatisfied_;
    nodeNumber_ = rhs.nodeNumber_;
    state_ = rhs.state_;
    if (nodeInfo_)
      assert((state_ & 2) != 0);
    else
      assert((state_ & 2) == 0);
  }
  return *this;
}
CbcNode::~CbcNode()
{
#ifdef CHECK_NODE
  if (nodeInfo_) {
    printf("CbcNode %p Destructor nodeInfo %p (%d)\n",
      this, nodeInfo_, nodeInfo_->numberPointingToThis());
    //assert(nodeInfo_->numberPointingToThis()>=0);
  } else {
    printf("CbcNode %p Destructor nodeInfo %p (?)\n",
      this, nodeInfo_);
  }
#endif
  if (nodeInfo_) {
    // was if (nodeInfo_&&(state_&2)!=0) {
    nodeInfo_->nullOwner();
    int numberToDelete = nodeInfo_->numberBranchesLeft();
    //    CbcNodeInfo * parent = nodeInfo_->parent();
    //assert (nodeInfo_->numberPointingToThis()>0);
    if (nodeInfo_->decrement(numberToDelete) == 0 || (state_ & 2) == 0) {
      if ((state_ & 2) == 0)
        nodeInfo_->nullParent();
      delete nodeInfo_;
    } else {
      //printf("node %p nodeinfo %p parent %p\n",this,nodeInfo_,nodeInfo_->parent());
      // anyway decrement parent
      //if (parent)
      ///parent->decrement(1);
    }
  }
  delete branch_;
}
// Decrement  active cut counts
void CbcNode::decrementCuts(int change)
{
  if (nodeInfo_)
    assert((state_ & 2) != 0);
  else
    assert((state_ & 2) == 0);
  if (nodeInfo_) {
    nodeInfo_->decrementCuts(change);
  }
}
void CbcNode::decrementParentCuts(CbcModel *model, int change)
{
  if (nodeInfo_)
    assert((state_ & 2) != 0);
  else
    assert((state_ & 2) == 0);
  if (nodeInfo_) {
    nodeInfo_->decrementParentCuts(model, change);
  }
}

/*
  Initialize reference counts (numberPointingToThis, numberBranchesLeft_)
  in the attached nodeInfo_.
*/
void CbcNode::initializeInfo()
{
  assert(nodeInfo_ && branch_);
  nodeInfo_->initializeInfo(branch_->numberBranches());
  assert((state_ & 2) != 0);
  assert(nodeInfo_->numberBranchesLeft() == branch_->numberBranchesLeft());
}
// Nulls out node info
void CbcNode::nullNodeInfo()
{
  nodeInfo_ = NULL;
  // say not active
  state_ &= ~2;
}

int CbcNode::branch(OsiSolverInterface *solver)
{
  double changeInGuessed;
  assert(nodeInfo_->numberBranchesLeft() == branch_->numberBranchesLeft());
  if (!solver)
    changeInGuessed = branch_->branch();
  else
    changeInGuessed = branch_->branch(solver);
  guessedObjectiveValue_ += changeInGuessed;
  //#define PRINTIT
#ifdef PRINTIT
  int numberLeft = nodeInfo_->numberBranchesLeft();
  CbcNodeInfo *parent = nodeInfo_->parent();
  int parentNodeNumber = -1;
  CbcBranchingObject *object1 = dynamic_cast< CbcBranchingObject * >(branch_);
  //OsiObject * object = object1->
  //int sequence = object->columnNumber);
  int id = -1;
  double value = 0.0;
  if (object1) {
    id = object1->variable();
    value = object1->value();
  }
  printf("id %d value %g objvalue %g\n", id, value, objectiveValue_);
  if (parent)
    parentNodeNumber = parent->nodeNumber();
  printf("Node number %d, %s, way %d, depth %d, parent node number %d\n",
    nodeInfo_->nodeNumber(), (numberLeft == 2) ? "leftBranch" : "rightBranch",
    way(), depth_, parentNodeNumber);
  assert(parentNodeNumber != nodeInfo_->nodeNumber());
#endif
  return nodeInfo_->branchedOn();
}
/* Active arm of the attached OsiBranchingObject.

   In the simplest instance, coded -1 for the down arm of the branch, +1 for
   the up arm. But see OsiBranchingObject::way()
   Use nodeInfo--.numberBranchesLeft_ to see how active

   Except that there is no OsiBranchingObject::way(), and this'll fail in any
   event because we have various OsiXXXBranchingObjects which aren't descended
   from CbcBranchingObjects. I think branchIndex() is the appropriate
   equivalent, but could be wrong. (lh, 061220)

   071212: I'm finally getting back to cbc-generic and rescuing a lot of my
   annotation from branches/devel (which was killed in summer). I'm going to
   put back an assert(obj) just to see what happens. It's still present as of
   the most recent change to CbcNode (r833).

   080104: Yep, we can arrive here with an OsiBranchingObject. Removed the
   assert, it's served its purpose.

   080226: John finally noticed this problem and added a way() method to the
   OsiBranchingObject hierarchy. Removing my workaround.

*/
int CbcNode::way() const
{
  if (branch_) {
    CbcBranchingObject *obj = dynamic_cast< CbcBranchingObject * >(branch_);
    if (obj) {
      return obj->way();
    } else {
      OsiTwoWayBranchingObject *obj2 = dynamic_cast< OsiTwoWayBranchingObject * >(branch_);
      assert(obj2);
      return obj2->way();
    }
  } else {
    return 0;
  }
}
/* Create a branching object for the node

   The routine scans the object list of the model and selects a set of
   unsatisfied objects as candidates for branching. The candidates are
   evaluated, and an appropriate branch object is installed.

   The numberPassesLeft is decremented to stop fixing one variable each time
   and going on and on (e.g. for stock cutting, air crew scheduling)

   If evaluation determines that an object is monotone or infeasible,
   the routine returns immediately. In the case of a monotone object,
   the branch object has already been called to modify the model.

   Return value:
   <ul>
   <li>  0: A branching object has been installed
   <li> -1: A monotone object was discovered
   <li> -2: An infeasible object was discovered
   </ul>
   Branch state:
   <ul>
   <li> -1: start
   <li> -1: A monotone object was discovered
   <li> -2: An infeasible object was discovered
   </ul>
*/
int CbcNode::chooseOsiBranch(CbcModel *model,
  CbcNode *lastNode,
  OsiBranchingInformation *usefulInfo,
  int branchState)
{
  int returnStatus = 0;
  if (lastNode)
    depth_ = lastNode->depth_ + 1;
  else
    depth_ = 0;
  OsiSolverInterface *solver = model->solver();
  objectiveValue_ = solver->getObjValue() * solver->getObjSense();
  usefulInfo->objectiveValue_ = objectiveValue_;
  usefulInfo->depth_ = depth_;
  const double *saveInfoSol = usefulInfo->solution_;
  double *saveSolution = new double[solver->getNumCols()];
  memcpy(saveSolution, solver->getColSolution(), solver->getNumCols() * sizeof(double));
  usefulInfo->solution_ = saveSolution;
  OsiChooseVariable *choose = model->branchingMethod()->chooseMethod();
  int numberUnsatisfied = -1;
  if (branchState < 0) {
    // initialize
    // initialize sum of "infeasibilities"
    sumInfeasibilities_ = 0.0;
    numberUnsatisfied = choose->setupList(usefulInfo, true);
    numberUnsatisfied_ = numberUnsatisfied;
    branchState = 0;
    if (numberUnsatisfied_ < 0) {
      // infeasible
      delete[] saveSolution;
      return -2;
    }
  }
  // unset best
  int best = -1;
  choose->setBestObjectIndex(-1);
  if (numberUnsatisfied) {
    if (branchState > 0 || !choose->numberOnList()) {
      // we need to return at once - don't do strong branching or anything
      if (choose->numberOnList() || !choose->numberStrong()) {
        best = choose->candidates()[0];
        choose->setBestObjectIndex(best);
      } else {
        // nothing on list - need to try again - keep any solution
        numberUnsatisfied = choose->setupList(usefulInfo, false);
        numberUnsatisfied_ = numberUnsatisfied;
        if (numberUnsatisfied) {
          best = choose->candidates()[0];
          choose->setBestObjectIndex(best);
        }
      }
    } else {
      // carry on with strong branching or whatever
      int returnCode = choose->chooseVariable(solver, usefulInfo, true);
      // update number of strong iterations etc
      model->incrementStrongInfo(choose->numberStrongDone(), choose->numberStrongIterations(),
        returnCode == -1 ? 0 : choose->numberStrongFixed(), returnCode == -1);
      if (returnCode > 1) {
        // has fixed some
        returnStatus = -1;
      } else if (returnCode == -1) {
        // infeasible
        returnStatus = -2;
      } else if (returnCode == 0) {
        // normal
        returnStatus = 0;
        numberUnsatisfied = 1;
      } else {
        // ones on list satisfied - double check
        numberUnsatisfied = choose->setupList(usefulInfo, false);
        numberUnsatisfied_ = numberUnsatisfied;
        if (numberUnsatisfied) {
          best = choose->candidates()[0];
          choose->setBestObjectIndex(best);
        }
      }
    }
  }
  delete branch_;
  branch_ = NULL;
  guessedObjectiveValue_ = COIN_DBL_MAX; //objectiveValue_; // for now
  if (!returnStatus) {
    if (numberUnsatisfied) {
      // create branching object
      const OsiObject *obj = model->solver()->object(choose->bestObjectIndex());
      //const OsiSolverInterface * solver = usefulInfo->solver_;
      branch_ = obj->createBranch(model->solver(), usefulInfo, obj->whichWay());
    }
  }
  usefulInfo->solution_ = saveInfoSol;
  delete[] saveSolution;
  // may have got solution
  if (choose->goodSolution()
    && model->problemFeasibility()->feasible(model, -1) >= 0) {
    // yes
    double objValue = choose->goodObjectiveValue();
    model->setBestSolution(CBC_STRONGSOL,
      objValue,
      choose->goodSolution());
    model->setLastHeuristic(NULL);
    model->incrementUsed(choose->goodSolution());
    choose->clearGoodSolution();
  }
  return returnStatus;
}
int CbcNode::chooseClpBranch(CbcModel *model,
  CbcNode *lastNode)
{
  assert(lastNode);
  depth_ = lastNode->depth_ + 1;
  delete branch_;
  branch_ = NULL;
  OsiSolverInterface *solver = model->solver();
  //double saveObjectiveValue = solver->getObjValue();
  //double objectiveValue = CoinMax(solver->getObjSense()*saveObjectiveValue,objectiveValue_);
  const double *lower = solver->getColLower();
  const double *upper = solver->getColUpper();
  // point to useful information
  OsiBranchingInformation usefulInfo = model->usefulInformation();
  // and modify
  usefulInfo.depth_ = depth_;
  int i;
  //bool beforeSolution = model->getSolutionCount()==0;
  int numberObjects = model->numberObjects();
  int numberColumns = model->getNumCols();
  double *saveUpper = new double[numberColumns];
  double *saveLower = new double[numberColumns];

  // Save solution in case heuristics need good solution later

  double *saveSolution = new double[numberColumns];
  memcpy(saveSolution, solver->getColSolution(), numberColumns * sizeof(double));
  model->reserveCurrentSolution(saveSolution);
  for (i = 0; i < numberColumns; i++) {
    saveLower[i] = lower[i];
    saveUpper[i] = upper[i];
  }
  // Save basis
  CoinWarmStart *ws = solver->getWarmStart();
  numberUnsatisfied_ = 0;
  // initialize sum of "infeasibilities"
  sumInfeasibilities_ = 0.0;
  // Note looks as if off end (hidden one)
  OsiObject *object = model->modifiableObject(numberObjects);
  CbcGeneralDepth *thisOne = dynamic_cast< CbcGeneralDepth * >(object);
  assert(thisOne);
  OsiClpSolverInterface *clpSolver
    = dynamic_cast< OsiClpSolverInterface * >(solver);
  assert(clpSolver);
  ClpSimplex *simplex = clpSolver->getModelPtr();
  int preferredWay;
  double infeasibility = object->infeasibility(&usefulInfo, preferredWay);
  if (thisOne->whichSolution() >= 0) {
    ClpNode *nodeInfo = NULL;
    if ((model->moreSpecialOptions() & 33554432) == 0) {
      nodeInfo = thisOne->nodeInfo(thisOne->whichSolution());
      nodeInfo->applyNode(simplex, 2);
    } else {
      // from diving
      CbcSubProblem **nodes = reinterpret_cast< CbcSubProblem ** >(model->temporaryPointer());
      assert(nodes);
      int numberDo = thisOne->numberNodes() - 1;
      for (int iNode = 0; iNode < numberDo; iNode++)
        nodes[iNode]->apply(solver, 1);
      nodes[numberDo]->apply(solver, 9 + 16);
    }
    int saveLogLevel = simplex->logLevel();
    simplex->setLogLevel(0);
    simplex->dual();
    simplex->setLogLevel(saveLogLevel);
    double cutoff = model->getCutoff();
    bool goodSolution = true;
    if (simplex->status()) {
      //simplex->writeMps("bad7.mps",2);
      if (nodeInfo) {
        if (nodeInfo->objectiveValue() > cutoff - 1.0e-2)
          goodSolution = false;
        else
          assert(!simplex->status());
      } else {
        // debug diving
        assert(!simplex->status());
      }
    }
    if (goodSolution) {
      double newObjectiveValue = solver->getObjSense() * solver->getObjValue();
      // See if integer solution
      int numInf;
      int numInf2;
      bool gotSol = model->feasibleSolution(numInf, numInf2);
      if (!gotSol) {
        COIN_DETAIL_PRINT(printf("numinf %d\n", numInf));
        double *sol = simplex->primalColumnSolution();
        for (int i = 0; i < numberColumns; i++) {
          if (simplex->isInteger(i)) {
            double value = floor(sol[i] + 0.5);
            if (fabs(value - sol[i]) > 1.0e-7) {
              COIN_DETAIL_PRINT(printf("%d value %g\n", i, sol[i]));
              if (fabs(value - sol[i]) < 1.0e-3) {
                sol[i] = value;
              }
            }
          }
        }
        simplex->writeMps("bad8.mps", 2);
        bool gotSol = model->feasibleSolution(numInf, numInf2);
        if (!gotSol)
          assert(gotSol);
      }
      model->setBestSolution(CBC_STRONGSOL,
        newObjectiveValue,
        solver->getColSolution());
      model->setLastHeuristic(NULL);
      model->incrementUsed(solver->getColSolution());
    }
  }
  // restore bounds
  {
    for (int j = 0; j < numberColumns; j++) {
      if (saveLower[j] != lower[j])
        solver->setColLower(j, saveLower[j]);
      if (saveUpper[j] != upper[j])
        solver->setColUpper(j, saveUpper[j]);
    }
  }
  // restore basis
  solver->setWarmStart(ws);
  delete ws;
  int anyAction;
  //#define CHECK_PATH
#ifdef CHECK_PATH
  extern int gotGoodNode_Z;
  if (gotGoodNode_Z >= 0)
    printf("good node %d %g\n", gotGoodNode_Z, infeasibility);
#endif
  if (infeasibility > 0.0) {
    if (infeasibility == COIN_DBL_MAX) {
      anyAction = -2; // infeasible
    } else {
      branch_ = thisOne->createCbcBranch(solver, &usefulInfo, preferredWay);
      if (branch_) {
        // Set to first one (and change when re-pushing)
        CbcGeneralBranchingObject *branch = dynamic_cast< CbcGeneralBranchingObject * >(branch_);
        branch->state(objectiveValue_, sumInfeasibilities_,
          numberUnsatisfied_, 0);
        branch->setNode(this);
        anyAction = 0;
      } else {
        anyAction = -2; // mark as infeasible
      }
    }
  } else {
    anyAction = -1;
  }
#ifdef CHECK_PATH
  gotGoodNode_Z = -1;
#endif
  // Set guessed solution value
  guessedObjectiveValue_ = objectiveValue_ + 1.0e-5;
  delete[] saveLower;
  delete[] saveUpper;

  // restore solution
  solver->setColSolution(saveSolution);
  delete[] saveSolution;
  return anyAction;
}
/* Double checks in case node can change its mind!
   Returns objective value
   Can change objective etc */
double
CbcNode::checkIsCutoff(double cutoff)
{
  branch_->checkIsCutoff(cutoff);
  return objectiveValue_;
}

/* vi: softtabstop=2 shiftwidth=2 expandtab tabstop=2
*/
