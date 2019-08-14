/* $Id$ */
// Copyright (C) 2003, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

#if defined(_MSC_VER)
// Turn off compiler warning about long names
#pragma warning(disable : 4786)
#endif
#include "CbcConfig.h"
#include <cassert>
#include <cstdlib>
#include <cmath>
#include <cfloat>

#ifdef COIN_HAS_CLP
#include "OsiClpSolverInterface.hpp"
#else
#include "OsiSolverInterface.hpp"
#endif
//#define CGL_DEBUG 1
#ifdef CGL_DEBUG
#include "OsiRowCutDebugger.hpp"
#endif
#include "CbcModel.hpp"
#include "CbcMessage.hpp"
#include "CbcCutGenerator.hpp"
#include "CbcBranchDynamic.hpp"
#include "CglProbing.hpp"
#include "CoinTime.hpp"
#ifdef CBC_THREAD
// need time on a thread by thread basis
#include <pthread.h>
#include <time.h>
#endif

// Default Constructor
CbcCutGenerator::CbcCutGenerator()
  : timeInCutGenerator_(0.0)
  , model_(NULL)
  , generator_(NULL)
  , generatorName_(NULL)
  , whenCutGenerator_(-1)
  , whenCutGeneratorInSub_(-100)
  , switchOffIfLessThan_(0)
  , depthCutGenerator_(-1)
  , depthCutGeneratorInSub_(-1)
  , inaccuracy_(0)
  , numberTimes_(0)
  , numberCuts_(0)
  , numberElements_(0)
  , numberColumnCuts_(0)
  , numberCutsActive_(0)
  , numberCutsAtRoot_(0)
  , numberActiveCutsAtRoot_(0)
  , numberShortCutsAtRoot_(0)
  , switches_(1)
  , maximumTries_(-1)
{
}
// Normal constructor
CbcCutGenerator::CbcCutGenerator(CbcModel *model, CglCutGenerator *generator,
  int howOften, const char *name,
  bool normal, bool atSolution,
  bool infeasible, int howOftenInSub,
  int whatDepth, int whatDepthInSub,
  int switchOffIfLessThan)
  : timeInCutGenerator_(0.0)
  , depthCutGenerator_(whatDepth)
  , depthCutGeneratorInSub_(whatDepthInSub)
  , inaccuracy_(0)
  , numberTimes_(0)
  , numberCuts_(0)
  , numberElements_(0)
  , numberColumnCuts_(0)
  , numberCutsActive_(0)
  , numberCutsAtRoot_(0)
  , numberActiveCutsAtRoot_(0)
  , numberShortCutsAtRoot_(0)
  , switches_(1)
  , maximumTries_(-1)
{
  if (howOften < -1900) {
    setGlobalCuts(true);
    howOften += 2000;
  } else if (howOften < -900) {
    setGlobalCutsAtRoot(true);
    howOften += 1000;
  }
  model_ = model;
  generator_ = generator->clone();
  generator_->refreshSolver(model_->solver());
  setNeedsOptimalBasis(generator_->needsOptimalBasis());
  whenCutGenerator_ = howOften;
  whenCutGeneratorInSub_ = howOftenInSub;
  switchOffIfLessThan_ = switchOffIfLessThan;
  if (name)
    generatorName_ = CoinStrdup(name);
  else
    generatorName_ = CoinStrdup("Unknown");
  setNormal(normal);
  setAtSolution(atSolution);
  setWhenInfeasible(infeasible);
}

// Copy constructor
CbcCutGenerator::CbcCutGenerator(const CbcCutGenerator &rhs)
{
  model_ = rhs.model_;
  generator_ = rhs.generator_->clone();
  //generator_->refreshSolver(model_->solver());
  whenCutGenerator_ = rhs.whenCutGenerator_;
  whenCutGeneratorInSub_ = rhs.whenCutGeneratorInSub_;
  switchOffIfLessThan_ = rhs.switchOffIfLessThan_;
  depthCutGenerator_ = rhs.depthCutGenerator_;
  depthCutGeneratorInSub_ = rhs.depthCutGeneratorInSub_;
  generatorName_ = CoinStrdup(rhs.generatorName_);
  switches_ = rhs.switches_;
  maximumTries_ = rhs.maximumTries_;
  timeInCutGenerator_ = rhs.timeInCutGenerator_;
  savedCuts_ = rhs.savedCuts_;
  inaccuracy_ = rhs.inaccuracy_;
  numberTimes_ = rhs.numberTimes_;
  numberCuts_ = rhs.numberCuts_;
  numberElements_ = rhs.numberElements_;
  numberColumnCuts_ = rhs.numberColumnCuts_;
  numberCutsActive_ = rhs.numberCutsActive_;
  numberCutsAtRoot_ = rhs.numberCutsAtRoot_;
  numberActiveCutsAtRoot_ = rhs.numberActiveCutsAtRoot_;
  numberShortCutsAtRoot_ = rhs.numberShortCutsAtRoot_;
}

// Assignment operator
CbcCutGenerator &
CbcCutGenerator::operator=(const CbcCutGenerator &rhs)
{
  if (this != &rhs) {
    delete generator_;
    free(generatorName_);
    model_ = rhs.model_;
    generator_ = rhs.generator_->clone();
    generator_->refreshSolver(model_->solver());
    whenCutGenerator_ = rhs.whenCutGenerator_;
    whenCutGeneratorInSub_ = rhs.whenCutGeneratorInSub_;
    switchOffIfLessThan_ = rhs.switchOffIfLessThan_;
    depthCutGenerator_ = rhs.depthCutGenerator_;
    depthCutGeneratorInSub_ = rhs.depthCutGeneratorInSub_;
    generatorName_ = CoinStrdup(rhs.generatorName_);
    switches_ = rhs.switches_;
    maximumTries_ = rhs.maximumTries_;
    timeInCutGenerator_ = rhs.timeInCutGenerator_;
    savedCuts_ = rhs.savedCuts_;
    inaccuracy_ = rhs.inaccuracy_;
    numberTimes_ = rhs.numberTimes_;
    numberCuts_ = rhs.numberCuts_;
    numberElements_ = rhs.numberElements_;
    numberColumnCuts_ = rhs.numberColumnCuts_;
    numberCutsActive_ = rhs.numberCutsActive_;
    numberCutsAtRoot_ = rhs.numberCutsAtRoot_;
    numberActiveCutsAtRoot_ = rhs.numberActiveCutsAtRoot_;
    numberShortCutsAtRoot_ = rhs.numberShortCutsAtRoot_;
  }
  return *this;
}

// Destructor
CbcCutGenerator::~CbcCutGenerator()
{
  free(generatorName_);
  delete generator_;
}

/* This is used to refresh any inforamtion.
   It also refreshes the solver in the cut generator
   in case generator wants to do some work
*/
void CbcCutGenerator::refreshModel(CbcModel *model)
{
  model_ = model;
  // added test - helps if generator not thread safe
  if (whenCutGenerator_ != -100)
    generator_->refreshSolver(model_->solver());
}
/* Generate cuts for the model data contained in si.
   The generated cuts are inserted into and returned in the
   collection of cuts cs.
*/
bool CbcCutGenerator::generateCuts(OsiCuts &cs, int fullScan, OsiSolverInterface *solver, CbcNode *node)
{
  /*
	  Make some decisions about whether we'll generate cuts. First convert
	  whenCutGenerator_ to a set of canonical values for comparison to the node
	  count.

		 0 <	mod 1000000, with a result of 0 forced to 1
	   -99 <= <= 0	convert to 1
	  -100 =	Off, period
	*/
  int depth;
  if (node)
    depth = node->depth();
  else
    depth = 0;
  int howOften = whenCutGenerator_;
  if (dynamic_cast< CglProbing * >(generator_)) {
    if (howOften == -100 && model_->doCutsNow(3)) {
      howOften = 1; // do anyway
    }
  }
  if (howOften == -100)
    return false;
  int pass = model_->getCurrentPassNumber() - 1;
  if (maximumTries_ > 0) {
    // howOften means what it says
    if ((pass % howOften) != 0 || depth)
      return false;
    else
      howOften = 1;
  }
  if (howOften > 0)
    howOften = howOften % 1000000;
  else
    howOften = 1;
  if (!howOften)
    howOften = 1;
  bool returnCode = false;
  //OsiSolverInterface * solver = model_->solver();
  // Reset cuts on first pass
  if (!pass)
    savedCuts_ = OsiCuts();
  /*
	  Determine if we should generate cuts based on node count.
	*/
  bool doThis = (model_->getNodeCount() % howOften) == 0;
  /*
	  If the user has provided a depth specification, it will override the node
	  count specification.
	*/
  if (depthCutGenerator_ > 0) {
    doThis = (depth % depthCutGenerator_) == 0;
    if (depth < depthCutGenerator_)
      doThis = true; // and also at top of tree
  }
  /*
	  A few magic numbers ...

	  The distinction between -100 and 100 for howOften is that we can override 100
	  with fullScan. -100 means no cuts, period. As does the magic number -200 for
	  whenCutGeneratorInSub_.
	*/

  // But turn off if 100
  if (howOften == 100)
    doThis = false;
  // Switch off if special setting
  if (whenCutGeneratorInSub_ == -200 && model_->parentModel()) {
    fullScan = 0;
    doThis = false;
  }
  if (fullScan || doThis) {
    CoinThreadRandom *randomNumberGenerator = NULL;
#ifdef COIN_HAS_CLP
    {
      OsiClpSolverInterface *clpSolver
        = dynamic_cast< OsiClpSolverInterface * >(solver);
      if (clpSolver)
        randomNumberGenerator = clpSolver->getModelPtr()->randomNumberGenerator();
    }
#endif
    double time1 = 0.0;
    //#undef CBC_THREAD
#if defined(_MSC_VER) || defined(__MSVCRT__) || !defined(CBC_THREAD)
    if (timing()) 
      time1 = CoinCpuTime();
#else
    struct timespec currTime;
    clockid_t threadClockId;
    if (timing()) {
      if (!model_->getNumberThreads()) {
	time1 = CoinCpuTime();
      } else {
	// Get thread clock Id
	pthread_getcpuclockid(pthread_self(), &threadClockId);
	// Using thread clock Id get the clock time
	clock_gettime(threadClockId, &currTime);
	time1 = static_cast<double>(currTime.tv_sec)
	  +1.0e-9*static_cast<double>(currTime.tv_nsec);
      }
    }
#endif
    //#define CBC_DEBUG
    int numberRowCutsBefore = cs.sizeRowCuts();
    int numberColumnCutsBefore = cs.sizeColCuts();
#ifdef JJF_ZERO
    int cutsBefore = cs.sizeCuts();
#endif
    CglTreeInfo info;
    info.level = depth;
    info.pass = pass;
    info.formulation_rows = model_->numberRowsAtContinuous();
    info.inTree = node != NULL;
    if (model_->parentModel()) {
      info.parentSolver = model_->parentModel()->continuousSolver();
      // indicate if doing full search
      info.hasParent = ((model_->specialOptions() & 67108864) == 0) ? 1 : 2;
    } else {
      info.hasParent = 0;
      info.parentSolver = NULL;
    }
    info.originalColumns = model_->originalColumns();
    info.randomNumberGenerator = randomNumberGenerator;
    info.options = (globalCutsAtRoot()) ? 8 : 0;
    if (ineffectualCuts())
      info.options |= 32;
    if (globalCuts())
      info.options |= 16;
    if (fullScan < 0)
      info.options |= 128;
    if (whetherInMustCallAgainMode())
      info.options |= 1024;
    // See if we want alternate set of cuts
    if ((model_->moreSpecialOptions() & 16384) != 0)
      info.options |= 256;
    if (model_->parentModel())
      info.options |= 512;
    // above had &&!model_->parentModel()&&depth<2)
    incrementNumberTimesEntered();
    CglProbing *generator = dynamic_cast< CglProbing * >(generator_);
    //if (!depth&&!pass)
    //printf("Cut generator %s when %d\n",generatorName_,whenCutGenerator_);
    if (!generator) {
      // Pass across model information in case it could be useful
      //void * saveData = solver->getApplicationData();
      //solver->setApplicationData(model_);
      generator_->generateCuts(*solver, cs, info);
      //solver->setApplicationData(saveData);
    } else {
      // Probing - return tight column bounds
      CglTreeProbingInfo *info2 = model_->probingInfo();
      bool doCuts = false;
      if (info2 && !depth) {
        info2->options = (globalCutsAtRoot()) ? 8 : 0;
        info2->level = depth;
        info2->pass = pass;
        info2->formulation_rows = model_->numberRowsAtContinuous();
        info2->inTree = node != NULL;
        if (model_->parentModel()) {
          info2->parentSolver = model_->parentModel()->continuousSolver();
          // indicate if doing full search
          info2->hasParent = ((model_->specialOptions() & 67108864) == 0) ? 1 : 2;
        } else {
          info2->hasParent = 0;
          info2->parentSolver = NULL;
        }
        info2->originalColumns = model_->originalColumns();
        info2->randomNumberGenerator = randomNumberGenerator;
        generator->generateCutsAndModify(*solver, cs, info2);
        doCuts = true;
      } else if (depth) {
        /* The idea behind this is that probing may work in a different
                   way deep in tree.  So every now and then try various
                   combinations to see what works.
                */
#define TRY_NOW_AND_THEN
#ifdef TRY_NOW_AND_THEN
        if ((numberTimes_ == 200 || (numberTimes_ > 200 && (numberTimes_ % 2000) == 0))
          && !model_->parentModel() && info.formulation_rows > 200) {
          /* In tree, every now and then try various combinations
                       maxStack, maxProbe (last 5 digits)
                       123 is special and means CglProbing will try and
                       be intelligent.
                    */
          int test[] = {
            100123,
            199999,
            200123,
            299999,
            500123,
            599999,
            1000123,
            1099999,
            2000123,
            2099999
          };
          int n = static_cast< int >(sizeof(test) / sizeof(int));
          int saveStack = generator->getMaxLook();
          int saveNumber = generator->getMaxProbe();
          int kr1 = 0;
          int kc1 = 0;
          int bestStackTree = -1;
          int bestNumberTree = -1;
          for (int i = 0; i < n; i++) {
            //OsiCuts cs2 = cs;
            int stack = test[i] / 100000;
            int number = test[i] - 100000 * stack;
            generator->setMaxLook(stack);
            generator->setMaxProbe(number);
            int numberRowCutsBefore = cs.sizeRowCuts();
            int numberColumnCutsBefore = cs.sizeColCuts();
            generator_->generateCuts(*solver, cs, info);
            int numberRowCuts = cs.sizeRowCuts() - numberRowCutsBefore;
            int numberColumnCuts = cs.sizeColCuts() - numberColumnCutsBefore;
#ifdef CLP_INVESTIGATE
            if (numberRowCuts < kr1 || numberColumnCuts < kc1)
              printf("Odd ");
#endif
            if (numberRowCuts > kr1 || numberColumnCuts > kc1) {
#ifdef CLP_INVESTIGATE
              printf("*** ");
#endif
              kr1 = numberRowCuts;
              kc1 = numberColumnCuts;
              bestStackTree = stack;
              bestNumberTree = number;
              doCuts = true;
            }
#ifdef CLP_INVESTIGATE
            printf("maxStack %d number %d gives %d row cuts and %d column cuts\n",
              stack, number, numberRowCuts, numberColumnCuts);
#endif
          }
          generator->setMaxLook(saveStack);
          generator->setMaxProbe(saveNumber);
          if (bestStackTree > 0) {
            generator->setMaxLook(bestStackTree);
            generator->setMaxProbe(bestNumberTree);
#ifdef CLP_INVESTIGATE
            printf("RRNumber %d -> %d, stack %d -> %d\n",
              saveNumber, bestNumberTree, saveStack, bestStackTree);
#endif
          } else {
            // no good
            generator->setMaxLook(0);
#ifdef CLP_INVESTIGATE
            printf("RRSwitching off number %d -> %d, stack %d -> %d\n",
              saveNumber, saveNumber, saveStack, 1);
#endif
          }
        }
#endif
        if (generator->getMaxLook() > 0 && !doCuts) {
          generator->generateCutsAndModify(*solver, cs, &info);
          doCuts = true;
        }
      } else {
        // at root - don't always do
        if (pass < 15 || (pass & 1) == 0) {
          generator->generateCutsAndModify(*solver, cs, &info);
          doCuts = true;
        }
      }
      if (doCuts && generator->tightLower()) {
        // probing may have tightened bounds - check
        const double *tightLower = generator->tightLower();
        const double *lower = solver->getColLower();
        const double *tightUpper = generator->tightUpper();
        const double *upper = solver->getColUpper();
        const double *solution = solver->getColSolution();
        int j;
        int numberColumns = solver->getNumCols();
        double primalTolerance = 1.0e-8;
        const char *tightenBounds = generator->tightenBounds();
#ifdef CGL_DEBUG
        const OsiRowCutDebugger *debugger = solver->getRowCutDebugger();
        if (debugger && debugger->onOptimalPath(*solver)) {
          printf("On optimal path CbcCut\n");
          int nCols = solver->getNumCols();
          int i;
          const double *optimal = debugger->optimalSolution();
          const double *objective = solver->getObjCoefficients();
          double objval1 = 0.0, objval2 = 0.0;
          for (i = 0; i < nCols; i++) {
#if CGL_DEBUG > 1
            printf("%d %g %g %g %g\n", i, lower[i], solution[i], upper[i], optimal[i]);
#endif
            objval1 += solution[i] * objective[i];
            objval2 += optimal[i] * objective[i];
            assert(optimal[i] >= lower[i] - 1.0e-5 && optimal[i] <= upper[i] + 1.0e-5);
            assert(optimal[i] >= tightLower[i] - 1.0e-5 && optimal[i] <= tightUpper[i] + 1.0e-5);
          }
          printf("current obj %g, integer %g\n", objval1, objval2);
        }
#endif
        bool feasible = true;
        if ((model_->getThreadMode() & 2) == 0) {
          for (j = 0; j < numberColumns; j++) {
            if (solver->isInteger(j)) {
              if (tightUpper[j] < upper[j]) {
                double nearest = floor(tightUpper[j] + 0.5);
                //assert (fabs(tightUpper[j]-nearest)<1.0e-5); may be infeasible
                solver->setColUpper(j, nearest);
                if (nearest < solution[j] - primalTolerance)
                  returnCode = true;
              }
              if (tightLower[j] > lower[j]) {
                double nearest = floor(tightLower[j] + 0.5);
                //assert (fabs(tightLower[j]-nearest)<1.0e-5); may be infeasible
                solver->setColLower(j, nearest);
                if (nearest > solution[j] + primalTolerance)
                  returnCode = true;
              }
            } else {
              if (upper[j] > lower[j]) {
                if (tightUpper[j] == tightLower[j]) {
                  // fix
                  //if (tightLower[j]!=lower[j])
                  solver->setColLower(j, tightLower[j]);
                  //if (tightUpper[j]!=upper[j])
                  solver->setColUpper(j, tightUpper[j]);
                  if (tightLower[j] > solution[j] + primalTolerance || tightUpper[j] < solution[j] - primalTolerance)
                    returnCode = true;
                } else if (tightenBounds && tightenBounds[j]) {
                  solver->setColLower(j, CoinMax(tightLower[j], lower[j]));
                  solver->setColUpper(j, CoinMin(tightUpper[j], upper[j]));
                  if (tightLower[j] > solution[j] + primalTolerance || tightUpper[j] < solution[j] - primalTolerance)
                    returnCode = true;
                }
              }
            }
            if (upper[j] < lower[j] - 1.0e-3) {
              feasible = false;
              break;
            }
          }
        } else {
          CoinPackedVector lbs;
          CoinPackedVector ubs;
          int numberChanged = 0;
          bool ifCut = false;
          for (j = 0; j < numberColumns; j++) {
            if (solver->isInteger(j)) {
              if (tightUpper[j] < upper[j]) {
                double nearest = floor(tightUpper[j] + 0.5);
                //assert (fabs(tightUpper[j]-nearest)<1.0e-5); may be infeasible
                ubs.insert(j, nearest);
                numberChanged++;
                if (nearest < solution[j] - primalTolerance)
                  ifCut = true;
              }
              if (tightLower[j] > lower[j]) {
                double nearest = floor(tightLower[j] + 0.5);
                //assert (fabs(tightLower[j]-nearest)<1.0e-5); may be infeasible
                lbs.insert(j, nearest);
                numberChanged++;
                if (nearest > solution[j] + primalTolerance)
                  ifCut = true;
              }
            } else {
              if (upper[j] > lower[j]) {
                if (tightUpper[j] == tightLower[j]) {
                  // fix
                  lbs.insert(j, tightLower[j]);
                  ubs.insert(j, tightUpper[j]);
                  if (tightLower[j] > solution[j] + primalTolerance || tightUpper[j] < solution[j] - primalTolerance)
                    ifCut = true;
                } else if (tightenBounds && tightenBounds[j]) {
                  lbs.insert(j, CoinMax(tightLower[j], lower[j]));
                  ubs.insert(j, CoinMin(tightUpper[j], upper[j]));
                  if (tightLower[j] > solution[j] + primalTolerance || tightUpper[j] < solution[j] - primalTolerance)
                    ifCut = true;
                }
              }
            }
            if (upper[j] < lower[j] - 1.0e-3) {
              feasible = false;
              break;
            }
          }
          if (numberChanged) {
            OsiColCut cc;
            cc.setUbs(ubs);
            cc.setLbs(lbs);
            if (ifCut) {
              cc.setEffectiveness(100.0);
            } else {
              cc.setEffectiveness(1.0e-5);
            }
            cs.insert(cc);
          }
        }
        if (!feasible) {
          // not feasible -add infeasible cut
          OsiRowCut rc;
          rc.setLb(COIN_DBL_MAX);
          rc.setUb(0.0);
          cs.insert(rc);
        }
      }
      //if (!solver->basisIsAvailable())
      //returnCode=true;
      if (!returnCode) {
        // bounds changed but still optimal
#ifdef COIN_HAS_CLP
        OsiClpSolverInterface *clpSolver
          = dynamic_cast< OsiClpSolverInterface * >(solver);
        if (clpSolver) {
          clpSolver->setLastAlgorithm(2);
        }
#endif
      }
#ifdef JJF_ZERO
      // Pass across info to pseudocosts
      char *mark = new char[numberColumns];
      memset(mark, 0, numberColumns);
      int nLook = generator->numberThisTime();
      const int *lookedAt = generator->lookedAt();
      const int *fixedDown = generator->fixedDown();
      const int *fixedUp = generator->fixedUp();
      for (j = 0; j < nLook; j++)
        mark[lookedAt[j]] = 1;
      int numberObjects = model_->numberObjects();
      for (int i = 0; i < numberObjects; i++) {
        CbcSimpleIntegerDynamicPseudoCost *obj1 = dynamic_cast< CbcSimpleIntegerDynamicPseudoCost * >(model_->modifiableObject(i));
        if (obj1) {
          int iColumn = obj1->columnNumber();
          if (mark[iColumn])
            obj1->setProbingInformation(fixedDown[iColumn], fixedUp[iColumn]);
        }
      }
      delete[] mark;
#endif
    }
    CbcCutModifier *modifier = model_->cutModifier();
    if (modifier) {
      int numberRowCutsAfter = cs.sizeRowCuts();
      int k;
      int nOdd = 0;
      //const OsiSolverInterface * solver = model_->solver();
      for (k = numberRowCutsAfter - 1; k >= numberRowCutsBefore; k--) {
        OsiRowCut &thisCut = cs.rowCut(k);
        int returnCode = modifier->modify(solver, thisCut);
        if (returnCode) {
          nOdd++;
          if (returnCode == 3)
            cs.eraseRowCut(k);
        }
      }
      if (nOdd)
        COIN_DETAIL_PRINT(printf("Cut generator %s produced %d cuts of which %d were modified\n",
          generatorName_, numberRowCutsAfter - numberRowCutsBefore, nOdd));
    }
    {
      // make all row cuts without test for duplicate
      int numberRowCutsAfter = cs.sizeRowCuts();
      int k;
#ifdef CGL_DEBUG
      const OsiRowCutDebugger *debugger = solver->getRowCutDebugger();
#endif
      //#define WEAKEN_CUTS 1
#ifdef WEAKEN_CUTS
      const double *lower = solver->getColLower();
      const double *upper = solver->getColUpper();
      const double *solution = solver->getColSolution();
#endif
      for (k = numberRowCutsBefore; k < numberRowCutsAfter; k++) {
        OsiRowCut *thisCut = cs.rowCutPtr(k);
#ifdef WEAKEN_CUTS
        // weaken cut if coefficients not integer

        double lb = thisCut->lb();
        double ub = thisCut->ub();
        if (lb < -1.0e100 || ub > 1.0e100) {
          // normal cut
          CoinPackedVector rpv = thisCut->row();
          const int n = rpv.getNumElements();
          const int *indices = rpv.getIndices();
          const double *elements = rpv.getElements();
          double bound = 0.0;
          double sum = 0.0;
          bool integral = true;
          int nInteger = 0;
          for (int k = 0; k < n; k++) {
            double value = fabs(elements[k]);
            int column = indices[k];
            sum += value;
            if (value != floor(value + 0.5))
              integral = false;
            if (solver->isInteger(column)) {
              nInteger++;
              double largerBound = CoinMax(fabs(lower[column]),
                fabs(upper[column]));
              double solutionBound = fabs(solution[column]) + 10.0;
              bound += CoinMin(largerBound, solutionBound);
            }
          }
#if WEAKEN_CUTS == 1
          // leave if all 0-1
          if (nInteger == bound)
            integral = true;
#endif
          if (!integral) {
            double weakenBy = 1.0e-7 * (bound + sum);
#if WEAKEN_CUTS > 2
            weakenBy *= 10.0;
#endif
            if (lb < -1.0e100)
              thisCut->setUb(ub + weakenBy);
            else
              thisCut->setLb(lb - weakenBy);
          }
        }
#endif
#ifdef CGL_DEBUG
        if (debugger && debugger->onOptimalPath(*solver)) {
#if CGL_DEBUG > 1
          const double *optimal = debugger->optimalSolution();
          CoinPackedVector rpv = thisCut->row();
          const int n = rpv.getNumElements();
          const int *indices = rpv.getIndices();
          const double *elements = rpv.getElements();

          double lb = thisCut->lb();
          double ub = thisCut->ub();
          double sum = 0.0;

          for (int k = 0; k < n; k++) {
            int column = indices[k];
            sum += optimal[column] * elements[k];
          }
          // is it nearly violated
          if (sum > ub - 1.0e-8 || sum < lb + 1.0e-8) {
            double violation = CoinMax(sum - ub, lb - sum);
            std::cout << generatorName_ << " cut with " << n
                      << " coefficients, nearly cuts off known solutions by " << violation
                      << ", lo=" << lb << ", ub=" << ub << std::endl;
            for (int k = 0; k < n; k++) {
              int column = indices[k];
              std::cout << "( " << column << " , " << elements[k] << " ) ";
              if ((k % 4) == 3)
                std::cout << std::endl;
            }
            std::cout << std::endl;
            std::cout << "Non zero solution values are" << std::endl;
            int j = 0;
            for (int k = 0; k < n; k++) {
              int column = indices[k];
              if (fabs(optimal[column]) > 1.0e-9) {
                std::cout << "( " << column << " , " << optimal[column] << " ) ";
                if ((j % 4) == 3)
                  std::cout << std::endl;
                j++;
              }
            }
            std::cout << std::endl;
          }
#endif
          assert(!debugger->invalidCut(*thisCut));
          if (debugger->invalidCut(*thisCut))
            abort();
        }
#endif
        thisCut->mutableRow().setTestForDuplicateIndex(false);
      }
    }
    // Add in saved cuts if violated
    if (false && !depth) {
      const double *solution = solver->getColSolution();
      double primalTolerance = 1.0e-7;
      int numberCuts = savedCuts_.sizeRowCuts();
      for (int k = numberCuts - 1; k >= 0; k--) {
        const OsiRowCut *thisCut = savedCuts_.rowCutPtr(k);
        double sum = 0.0;
        int n = thisCut->row().getNumElements();
        const int *column = thisCut->row().getIndices();
        const double *element = thisCut->row().getElements();
        assert(n);
        for (int i = 0; i < n; i++) {
          double value = element[i];
          sum += value * solution[column[i]];
        }
        if (sum > thisCut->ub() + primalTolerance) {
          sum = sum - thisCut->ub();
        } else if (sum < thisCut->lb() - primalTolerance) {
          sum = thisCut->lb() - sum;
        } else {
          sum = 0.0;
        }
        if (sum) {
          // add to candidates and take out here
          cs.insert(*thisCut);
          savedCuts_.eraseRowCut(k);
        }
      }
    }
    if (!atSolution()) {
      int numberRowCutsAfter = cs.sizeRowCuts();
      int k;
      int nEls = 0;
      int nCuts = numberRowCutsAfter - numberRowCutsBefore;
      // Remove NULL cuts!
      int nNull = 0;
      const double *solution = solver->getColSolution();
      bool feasible = true;
      double primalTolerance = 1.0e-7;
      int shortCut = (depth) ? -1 : generator_->maximumLengthOfCutInTree();
      for (k = numberRowCutsAfter - 1; k >= numberRowCutsBefore; k--) {
        const OsiRowCut *thisCut = cs.rowCutPtr(k);
        double sum = 0.0;
        if (thisCut->lb() <= thisCut->ub()) {
          int n = thisCut->row().getNumElements();
          if (n <= shortCut)
            numberShortCutsAtRoot_++;
          const int *column = thisCut->row().getIndices();
          const double *element = thisCut->row().getElements();
          if (n <= 0) {
            // infeasible cut - give up
            feasible = false;
            break;
          }
          nEls += n;
          for (int i = 0; i < n; i++) {
            double value = element[i];
            sum += value * solution[column[i]];
          }
          if (sum > thisCut->ub() + primalTolerance) {
            sum = sum - thisCut->ub();
          } else if (sum < thisCut->lb() - primalTolerance) {
            sum = thisCut->lb() - sum;
          } else {
            sum = 0.0;
            cs.eraseRowCut(k);
            nNull++;
          }
        }
      }
      //if (nNull)
      //printf("%s has %d cuts and %d elements - %d null!\n",generatorName_,
      //       nCuts,nEls,nNull);
      numberRowCutsAfter = cs.sizeRowCuts();
      nCuts = numberRowCutsAfter - numberRowCutsBefore;
      nEls = 0;
      for (k = numberRowCutsBefore; k < numberRowCutsAfter; k++) {
        const OsiRowCut *thisCut = cs.rowCutPtr(k);
        int n = thisCut->row().getNumElements();
        nEls += n;
      }
      //printf("%s has %d cuts and %d elements\n",generatorName_,
      //     nCuts,nEls);
      CoinBigIndex nElsNow = solver->getMatrixByCol()->getNumElements();
      int numberColumns = solver->getNumCols();
      int numberRows = solver->getNumRows();
      //double averagePerRow = static_cast<double>(nElsNow)/
      //static_cast<double>(numberRows);
      CoinBigIndex nAdd;
      CoinBigIndex nAdd2;
      CoinBigIndex nReasonable;
      if (!model_->parentModel() && depth < 2) {
        if (inaccuracy_ < 3) {
          nAdd = 10000;
          if (pass > 0 && numberColumns > -500)
            nAdd = CoinMin(nAdd, nElsNow + 2 * numberRows);
        } else {
          nAdd = 10000;
          if (pass > 0)
            nAdd = CoinMin(nAdd, nElsNow + 2 * numberRows);
        }
        nAdd2 = 5 * numberColumns;
        nReasonable = CoinMax(nAdd2, nElsNow / 8 + nAdd);
        if (!depth && !pass) {
          // allow more
          nAdd += nElsNow / 2;
          nAdd2 += nElsNow / 2;
          nReasonable += nElsNow / 2;
        }
        //if (!depth&&ineffectualCuts())
        //nReasonable *= 2;
      } else {
        nAdd = 200;
        nAdd2 = 2 * numberColumns;
        nReasonable = CoinMax(nAdd2, nElsNow / 8 + nAdd);
      }
      //#define UNS_WEIGHT 0.1
#ifdef UNS_WEIGHT
      const double *colLower = solver->getColLower();
      const double *colUpper = solver->getColUpper();
#endif
      if (/*nEls>CoinMax(nAdd2,nElsNow/8+nAdd)*/ nCuts && feasible) {
        //printf("need to remove cuts\n");
        // just add most effective
#ifndef JJF_ONE
        CoinBigIndex nDelete = nEls - nReasonable;

        nElsNow = nEls;
        double *sort = new double[nCuts];
        int *which = new int[nCuts];
        // For parallel cuts
        double *element2 = new double[numberColumns];
        //#define USE_OBJECTIVE 2
#ifdef USE_OBJECTIVE
        const double *objective = solver->getObjCoefficients();
#if USE_OBJECTIVE > 1
        double objNorm = 0.0;
        for (int i = 0; i < numberColumns; i++)
          objNorm += objective[i] * objective[i];
        if (objNorm)
          objNorm = 1.0 / sqrt(objNorm);
        else
          objNorm = 1.0;
        objNorm *= 0.01; // downgrade
#endif
#endif
        CoinZeroN(element2, numberColumns);
        for (k = numberRowCutsBefore; k < numberRowCutsAfter; k++) {
          const OsiRowCut *thisCut = cs.rowCutPtr(k);
          double sum = 0.0;
          if (thisCut->lb() <= thisCut->ub()) {
            int n = thisCut->row().getNumElements();
            const int *column = thisCut->row().getIndices();
            const double *element = thisCut->row().getElements();
            assert(n);
#ifdef UNS_WEIGHT
            double normU = 0.0;
            double norm = 1.0e-3;
            int nU = 0;
            for (int i = 0; i < n; i++) {
              double value = element[i];
              int iColumn = column[i];
              double solValue = solution[iColumn];
              sum += value * solValue;
              value *= value;
              norm += value;
              if (solValue > colLower[iColumn] + 1.0e-6 && solValue < colUpper[iColumn] - 1.0e-6) {
                normU += value;
                nU++;
              }
            }
#ifdef JJF_ZERO
            int nS = n - nU;
            if (numberColumns > 20000) {
              if (nS > 50) {
                double ratio = 50.0 / nS;
                normU /= ratio;
              }
            }
#endif
            norm += UNS_WEIGHT * (normU - norm);
#else
            double norm = 1.0e-3;
#ifdef USE_OBJECTIVE
            double obj = 0.0;
#endif
            for (int i = 0; i < n; i++) {
              int iColumn = column[i];
              double value = element[i];
              sum += value * solution[iColumn];
              norm += value * value;
#ifdef USE_OBJECTIVE
              obj += value * objective[iColumn];
#endif
            }
#endif
            if (sum > thisCut->ub()) {
              sum = sum - thisCut->ub();
            } else if (sum < thisCut->lb()) {
              sum = thisCut->lb() - sum;
            } else {
              sum = 0.0;
            }
#ifdef USE_OBJECTIVE
            if (sum) {
#if USE_OBJECTIVE == 1
              obj = CoinMax(1.0e-6, fabs(obj));
              norm = sqrt(obj * norm);
              //sum += fabs(obj)*invObjNorm;
              //printf("sum %g norm %g normobj %g invNorm %g mod %g\n",
              //     sum,norm,obj,invObjNorm,obj*invObjNorm);
              // normalize
              sum /= sqrt(norm);
#else
              // normalize
              norm = 1.0 / sqrt(norm);
              sum = (sum + objNorm * obj) * norm;
#endif
            }
#else
            // normalize
            sum /= sqrt(norm);
#endif
            //sum /= pow(norm,0.3);
            // adjust for length
            //sum /= pow(reinterpret_cast<double>(n),0.2);
            //sum /= sqrt((double) n);
            // randomize
            //double randomNumber =
            //model_->randomNumberGenerator()->randomDouble();
            //sum *= (0.5+randomNumber);
          } else {
            // keep
            sum = COIN_DBL_MAX;
          }
          sort[k - numberRowCutsBefore] = sum;
          which[k - numberRowCutsBefore] = k;
        }
        CoinSort_2(sort, sort + nCuts, which);
        // Now see which ones are too similar
        int nParallel = 0;
        double testValue = (depth > 1) ? 0.99 : 0.999999;
        for (k = 0; k < nCuts; k++) {
          int j = which[k];
          const OsiRowCut *thisCut = cs.rowCutPtr(j);
          if (thisCut->lb() > thisCut->ub())
            break; // cut is infeasible
          int n = thisCut->row().getNumElements();
          const int *column = thisCut->row().getIndices();
          const double *element = thisCut->row().getElements();
          assert(n);
          double norm = 0.0;
          double lb = thisCut->lb();
          double ub = thisCut->ub();
          for (int i = 0; i < n; i++) {
            double value = element[i];
            element2[column[i]] = value;
            norm += value * value;
          }
          int kkk = CoinMin(nCuts, k + 5);
          for (int kk = k + 1; kk < kkk; kk++) {
            int jj = which[kk];
            const OsiRowCut *thisCut2 = cs.rowCutPtr(jj);
            if (thisCut2->lb() > thisCut2->ub())
              break; // cut is infeasible
            int nB = thisCut2->row().getNumElements();
            const int *columnB = thisCut2->row().getIndices();
            const double *elementB = thisCut2->row().getElements();
            assert(nB);
            double normB = 0.0;
            double product = 0.0;
            for (int i = 0; i < nB; i++) {
              double value = elementB[i];
              normB += value * value;
              product += value * element2[columnB[i]];
            }
            if (product > 0.0 && product * product > testValue * norm * normB) {
              bool parallel = true;
              double lbB = thisCut2->lb();
              double ubB = thisCut2->ub();
              if ((lb < -1.0e20 && lbB > -1.0e20) || (lbB < -1.0e20 && lb > -1.0e20))
                parallel = false;
              double tolerance;
              tolerance = CoinMax(fabs(lb), fabs(lbB)) + 1.0e-6;
              if (fabs(lb - lbB) > tolerance)
                parallel = false;
              if ((ub > 1.0e20 && ubB < 1.0e20) || (ubB > 1.0e20 && ub < 1.0e20))
                parallel = false;
              tolerance = CoinMax(fabs(ub), fabs(ubB)) + 1.0e-6;
              if (fabs(ub - ubB) > tolerance)
                parallel = false;
              if (parallel) {
                nParallel++;
                sort[k] = 0.0;
                break;
              }
            }
          }
          for (int i = 0; i < n; i++) {
            element2[column[i]] = 0.0;
          }
        }
        delete[] element2;
        CoinSort_2(sort, sort + nCuts, which);
        k = 0;
        while (nDelete > 0 || !sort[k]) {
          int iCut = which[k];
          const OsiRowCut *thisCut = cs.rowCutPtr(iCut);
          int n = thisCut->row().getNumElements();
          // may be best, just to save if short
          if (false && n && sort[k]) {
            // add to saved cuts
            savedCuts_.insert(*thisCut);
          }
          nDelete -= n;
          k++;
          if (k >= nCuts)
            break;
        }
        std::sort(which, which + k);
        k--;
        for (; k >= 0; k--) {
          cs.eraseRowCut(which[k]);
        }
        delete[] sort;
        delete[] which;
        numberRowCutsAfter = cs.sizeRowCuts();
#else
        double *norm = new double[nCuts];
        int *which = new int[2 * nCuts];
        double *score = new double[nCuts];
        double *ortho = new double[nCuts];
        int nIn = 0;
        int nOut = nCuts;
        // For parallel cuts
        double *element2 = new double[numberColumns];
        const double *objective = solver->getObjCoefficients();
        double objNorm = 0.0;
        for (int i = 0; i < numberColumns; i++)
          objNorm += objective[i] * objective[i];
        if (objNorm)
          objNorm = 1.0 / sqrt(objNorm);
        else
          objNorm = 1.0;
        objNorm *= 0.1; // weight of 0.1
        CoinZeroN(element2, numberColumns);
        int numberRowCuts = numberRowCutsAfter - numberRowCutsBefore;
        int iBest = -1;
        double best = 0.0;
        int nPossible = 0;
        double testValue = (depth > 1) ? 0.7 : 0.5;
        for (k = 0; k < numberRowCuts; k++) {
          const OsiRowCut *thisCut = cs.rowCutPtr(k + numberRowCutsBefore);
          double sum = 0.0;
          if (thisCut->lb() <= thisCut->ub()) {
            int n = thisCut->row().getNumElements();
            const int *column = thisCut->row().getIndices();
            const double *element = thisCut->row().getElements();
            assert(n);
            double normThis = 1.0e-6;
            double obj = 0.0;
            for (int i = 0; i < n; i++) {
              int iColumn = column[i];
              double value = element[i];
              sum += value * solution[iColumn];
              normThis += value * value;
              obj += value * objective[iColumn];
            }
            if (sum > thisCut->ub()) {
              sum = sum - thisCut->ub();
            } else if (sum < thisCut->lb()) {
              sum = thisCut->lb() - sum;
            } else {
              sum = 0.0;
            }
            if (sum) {
              normThis = 1.0 / sqrt(normThis);
              norm[k] = normThis;
              sum *= normThis;
              obj *= normThis;
              score[k] = sum + obj * objNorm;
              ortho[k] = 1.0;
            }
          } else {
            // keep and discard others
            nIn = 1;
            which[0] = k;
            for (int j = 0; j < numberRowCuts; j++) {
              if (j != k)
                which[nOut++] = j;
            }
            iBest = -1;
            break;
          }
          if (sum) {
            if (score[k] > best) {
              best = score[k];
              iBest = nPossible;
            }
            which[nPossible++] = k;
          } else {
            which[nOut++] = k;
          }
        }
        while (iBest >= 0) {
          int kBest = which[iBest];
          int j = which[nIn];
          which[iBest] = j;
          which[nIn++] = kBest;
          const OsiRowCut *thisCut = cs.rowCutPtr(kBest + numberRowCutsBefore);
          int n = thisCut->row().getNumElements();
          nReasonable -= n;
          if (nReasonable <= 0) {
            for (k = nIn; k < nPossible; k++)
              which[nOut++] = which[k];
            break;
          }
          // Now see which ones are too similar and choose next
          iBest = -1;
          best = 0.0;
          int nOld = nPossible;
          nPossible = nIn;
          const int *column = thisCut->row().getIndices();
          const double *element = thisCut->row().getElements();
          assert(n);
          double normNew = norm[kBest];
          for (int i = 0; i < n; i++) {
            double value = element[i];
            element2[column[i]] = value;
          }
          for (int j = nIn; j < nOld; j++) {
            k = which[j];
            const OsiRowCut *thisCut2 = cs.rowCutPtr(k + numberRowCutsBefore);
            int nB = thisCut2->row().getNumElements();
            const int *columnB = thisCut2->row().getIndices();
            const double *elementB = thisCut2->row().getElements();
            assert(nB);
            double normB = norm[k];
            double product = 0.0;
            for (int i = 0; i < nB; i++) {
              double value = elementB[i];
              product += value * element2[columnB[i]];
            }
            double orthoScore = 1.0 - product * normNew * normB;
            if (orthoScore >= testValue) {
              ortho[k] = CoinMin(orthoScore, ortho[k]);
              double test = score[k] + ortho[k];
              if (test > best) {
                best = score[k];
                iBest = nPossible;
              }
              which[nPossible++] = k;
            } else {
              which[nOut++] = k;
            }
          }
          for (int i = 0; i < n; i++) {
            element2[column[i]] = 0.0;
          }
        }
        delete[] score;
        delete[] ortho;
        std::sort(which + nCuts, which + nOut);
        k = nOut - 1;
        for (; k >= nCuts; k--) {
          cs.eraseRowCut(which[k] + numberRowCutsBefore);
        }
        delete[] norm;
        delete[] which;
        numberRowCutsAfter = cs.sizeRowCuts();
#endif
      }
    }
#ifdef CBC_DEBUG
    {
      int numberRowCutsAfter = cs.sizeRowCuts();
      int k;
      int nBad = 0;
      for (k = numberRowCutsBefore; k < numberRowCutsAfter; k++) {
        OsiRowCut thisCut = cs.rowCut(k);
        if (thisCut.lb() > thisCut.ub() || thisCut.lb() > 1.0e8 || thisCut.ub() < -1.0e8)
          printf("cut from %s has bounds %g and %g!\n",
            generatorName_, thisCut.lb(), thisCut.ub());
        if (thisCut.lb() <= thisCut.ub()) {
          /* check size of elements.
                       We can allow smaller but this helps debug generators as it
                       is unsafe to have small elements */
          int n = thisCut.row().getNumElements();
          const int *column = thisCut.row().getIndices();
          const double *element = thisCut.row().getElements();
          assert(n);
          for (int i = 0; i < n; i++) {
            double value = element[i];
            if (fabs(value) <= 1.0e-12 || fabs(value) >= 1.0e20)
              nBad++;
          }
        }
        if (nBad)
          printf("Cut generator %s produced %d cuts of which %d had tiny or large elements\n",
            generatorName_, numberRowCutsAfter - numberRowCutsBefore, nBad);
      }
    }
#endif
    int numberRowCutsAfter = cs.sizeRowCuts();
    int numberColumnCutsAfter = cs.sizeColCuts();
    if (numberRowCutsBefore < numberRowCutsAfter) {
      for (int k = numberRowCutsBefore; k < numberRowCutsAfter; k++) {
        OsiRowCut thisCut = cs.rowCut(k);
        int n = thisCut.row().getNumElements();
        numberElements_ += n;
      }
#ifdef JJF_ZERO
      printf("generator %s generated %d row cuts\n",
        generatorName_, numberRowCutsAfter - numberRowCutsBefore);
#endif
      numberCuts_ += numberRowCutsAfter - numberRowCutsBefore;
    }
    if (numberColumnCutsBefore < numberColumnCutsAfter) {
#ifdef JJF_ZERO
      printf("generator %s generated %d column cuts\n",
        generatorName_, numberColumnCutsAfter - numberColumnCutsBefore);
#endif
      numberColumnCuts_ += numberColumnCutsAfter - numberColumnCutsBefore;
    }
    if (timing()) {
      // Using thread clock Id get the clock time
#if defined(_MSC_VER) || defined(__MSVCRT__) || !defined(CBC_THREAD)
      timeInCutGenerator_ += CoinCpuTime() - time1;
#else
      if (!model_->getNumberThreads()) {
	timeInCutGenerator_ += CoinCpuTime() - time1;
      } else {
	clock_gettime(threadClockId, &currTime);
	timeInCutGenerator_ += static_cast<double>(currTime.tv_sec) + 1.0e-9*
	  static_cast<double>(currTime.tv_nsec) - time1;
      }
#endif
    }
    // switch off if first time and no good
    if (node == NULL && !pass) {
      if (numberRowCutsAfter - numberRowCutsBefore
        < switchOffIfLessThan_ /*&& numberCuts_ < switchOffIfLessThan_*/) {
        // switch off
        maximumTries_ = 0;
        whenCutGenerator_ = -100;
        //whenCutGenerator_ = -100;
        //whenCutGeneratorInSub_ = -200;
      }
    }
    if (maximumTries_ > 0) {
      maximumTries_--;
      if (!maximumTries_)
        whenCutGenerator_ = -100;
    }
  }
  return returnCode;
}
void CbcCutGenerator::setHowOften(int howOften)
{

  if (howOften >= 1000000) {
    // leave Probing every SCANCUTS_PROBING
    howOften = howOften % 1000000;
    CglProbing *generator = dynamic_cast< CglProbing * >(generator_);

    if (generator && howOften > SCANCUTS_PROBING)
      howOften = SCANCUTS_PROBING + 1000000;
    else
      howOften += 1000000;
  }
  whenCutGenerator_ = howOften;
}
void CbcCutGenerator::setWhatDepth(int value)
{
  depthCutGenerator_ = value;
}
void CbcCutGenerator::setWhatDepthInSub(int value)
{
  depthCutGeneratorInSub_ = value;
}
// Add in statistics from other
void CbcCutGenerator::addStatistics(const CbcCutGenerator *other)
{
  // Time in cut generator
  timeInCutGenerator_ += other->timeInCutGenerator_;
  // Number times cut generator entered
  numberTimes_ += other->numberTimes_;
  // Total number of cuts added
  numberCuts_ += other->numberCuts_;
  // Total number of elements added
  numberElements_ += other->numberElements_;
  // Total number of column cuts added
  numberColumnCuts_ += other->numberColumnCuts_;
  // Total number of cuts active after (at end of n cut passes at each node)
  numberCutsActive_ += other->numberCutsActive_;
  // Number of cuts generated at root
  numberCutsAtRoot_ += other->numberCutsAtRoot_;
  // Number of cuts active at root
  numberActiveCutsAtRoot_ += other->numberActiveCutsAtRoot_;
  // Number of short cuts at root
  numberShortCutsAtRoot_ += other->numberShortCutsAtRoot_;
}
// Scale back statistics by factor
void CbcCutGenerator::scaleBackStatistics(int factor)
{
  // leave time
  // Number times cut generator entered
  numberTimes_ = (numberTimes_ + factor - 1) / factor;
  // Total number of cuts added
  numberCuts_ = (numberCuts_ + factor - 1) / factor;
  // Total number of elements added
  numberElements_ = (numberElements_ + factor - 1) / factor;
  // Total number of column cuts added
  numberColumnCuts_ = (numberColumnCuts_ + factor - 1) / factor;
  // Total number of cuts active after (at end of n cut passes at each node)
  numberCutsActive_ = (numberCutsActive_ + factor - 1) / factor;
  // Number of cuts generated at root
  numberCutsAtRoot_ = (numberCutsAtRoot_ + factor - 1) / factor;
  // Number of cuts active at root
  numberActiveCutsAtRoot_ = (numberActiveCutsAtRoot_ + factor - 1) / factor;
  // Number of short cuts at root
  numberShortCutsAtRoot_ = (numberShortCutsAtRoot_ + factor - 1) / factor;
}
// Create C++ lines to get to current state
void CbcCutGenerator::generateTuning(FILE *fp)
{
  fprintf(fp, "// Cbc tuning for generator %s\n", generatorName_);
  fprintf(fp, "   generator->setHowOften(%d);\n", whenCutGenerator_);
  fprintf(fp, "   generator->setSwitchOffIfLessThan(%d);\n", switchOffIfLessThan_);
  fprintf(fp, "   generator->setWhatDepth(%d);\n", depthCutGenerator_);
  fprintf(fp, "   generator->setInaccuracy(%d);\n", inaccuracy_);
  if (timing())
    fprintf(fp, "   generator->setTiming(true);\n");
  if (normal())
    fprintf(fp, "   generator->setNormal(true);\n");
  if (atSolution())
    fprintf(fp, "   generator->setAtSolution(true);\n");
  if (whenInfeasible())
    fprintf(fp, "   generator->setWhenInfeasible(true);\n");
  if (needsOptimalBasis())
    fprintf(fp, "   generator->setNeedsOptimalBasis(true);\n");
  if (mustCallAgain())
    fprintf(fp, "   generator->setMustCallAgain(true);\n");
  if (whetherToUse())
    fprintf(fp, "   generator->setWhetherToUse(true);\n");
}

/* vi: softtabstop=2 shiftwidth=2 expandtab tabstop=2
*/
