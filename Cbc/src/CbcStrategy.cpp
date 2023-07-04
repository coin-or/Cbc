/* $Id$ */
// Copyright (C) 2005, International Business Machines
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

#include "OsiSolverInterface.hpp"
#ifdef COIN_HAS_CLP
#include "OsiClpSolverInterface.hpp"
#endif
#include "CbcModel.hpp"
#include "CbcMessage.hpp"
#include "CbcStrategy.hpp"
#include "CbcCutGenerator.hpp"
#include "CbcBranchActual.hpp"
#include "CbcNode.hpp"
#include "CoinWarmStart.hpp"
#include "CglPreProcess.hpp"
// Cuts

#include "CglGomory.hpp"
#include "CglProbing.hpp"
#include "CglKnapsackCover.hpp"
#include "CglOddHole.hpp"
#include "CglClique.hpp"
#include "CglFlowCover.hpp"
#include "CglMixedIntegerRounding2.hpp"

// Heuristics

#include "CbcHeuristic.hpp"
#include "CbcHeuristicLocal.hpp"
#include "CbcHeuristicRINS.hpp"

// Default Constructor
CbcStrategy::CbcStrategy()
  : depth_(0)
  , preProcessState_(0)
  , process_(NULL)
{
}

// Destructor
CbcStrategy::~CbcStrategy()
{
  delete process_;
}
// Delete pre-processing object to save memory
void CbcStrategy::deletePreProcess()
{
  delete process_;
  process_ = NULL;
}
// Return a new Full node information pointer (descendant of CbcFullNodeInfo)
CbcNodeInfo *
CbcStrategy::fullNodeInfo(CbcModel *model, int numberRowsAtContinuous) const
{
  return new CbcFullNodeInfo(model, numberRowsAtContinuous);
}
// Return a new Partial node information pointer (descendant of CbcPartialNodeInfo)
CbcNodeInfo *
CbcStrategy::partialNodeInfo(CbcModel * /*model*/,
  CbcNodeInfo *parent, CbcNode *owner,
  int numberChangedBounds, const int *variables,
  const double *boundChanges,
  const CoinWarmStartDiff *basisDiff) const
{
  return new CbcPartialNodeInfo(parent, owner, numberChangedBounds, variables,
    boundChanges, basisDiff);
}
/* After a CbcModel::resolve this can return a status
   -1 no effect
   0 treat as optimal
   1 as 0 but do not do any more resolves (i.e. no more cuts)
   2 treat as infeasible
*/
int CbcStrategy::status(CbcModel * /*model*/, CbcNodeInfo * /*parent*/,
  int /*whereFrom*/)
{
  return -1;
}

// Default Constructor
CbcStrategyDefault::CbcStrategyDefault(int cutsOnlyAtRoot,
  int numberStrong,
  int numberBeforeTrust,
  int printLevel)
  : CbcStrategy()
  , cutsOnlyAtRoot_(cutsOnlyAtRoot)
  , numberStrong_(numberStrong)
  , numberBeforeTrust_(numberBeforeTrust)
  , printLevel_(printLevel)
  , desiredPreProcess_(0)
  , preProcessPasses_(0)
{
}

// Destructor
CbcStrategyDefault::~CbcStrategyDefault()
{
}

// Clone
CbcStrategy *
CbcStrategyDefault::clone() const
{
  return new CbcStrategyDefault(*this);
}

// Copy constructor
CbcStrategyDefault::CbcStrategyDefault(const CbcStrategyDefault &rhs)
  : CbcStrategy(rhs)
  , cutsOnlyAtRoot_(rhs.cutsOnlyAtRoot_)
  , numberStrong_(rhs.numberStrong_)
  , numberBeforeTrust_(rhs.numberBeforeTrust_)
  , printLevel_(rhs.printLevel_)
  , desiredPreProcess_(rhs.desiredPreProcess_)
  , preProcessPasses_(rhs.preProcessPasses_)
{
  setNested(rhs.getNested());
}

/*
  Set up cut generators. Will instantiate Probing, Gomory, Knapsack, Clique,
  FlowCover, and MIR2 generators. Probing should be the first in the vector
  of generators as it tightens bounds on continuous variables.

  Cut generators already installed will dominate cut generators instantiated
  here.

  There's a classic magic number overloaded parameter example here. The
  variable genFlags below is interpreted as single-bit flags to control
  whether a cut generator will be instantiated: Probing:1, Gomory:2,
  Knapsack:4, Clique:8, FlowCover:16, MIR2:32. Normally it's hardcoded to 63.
  If CBC_GENERATE_TEST is defined, and the model's node limit is set between
  190000 and 190064, genFlags is loaded with the low-order bits.
*/
void CbcStrategyDefault::setupCutGenerators(CbcModel &model)
{
  if (cutsOnlyAtRoot_ < 0)
    return; // no cuts wanted

  // Magic number overloaded parameter -- see comment at head.
  int genFlags = 63;
#ifdef CBC_GENERATE_TEST
  int nNodes = model.getMaximumNodes();
  if (nNodes >= 190000 && nNodes < 190064)
    genFlags = nNodes - 190000;
#endif

  CglProbing generator1;
  generator1.setUsingObjective(true);
  generator1.setMaxPass(1);
  generator1.setMaxPassRoot(1);
  // Number of unsatisfied variables to look at
  generator1.setMaxProbe(10);
  // How far to follow the consequences
  generator1.setMaxLook(10);
  // Only look at rows with fewer than this number of elements
  generator1.setMaxElements(200);
  generator1.setMaxElementsRoot(300);
  //generator1.setRowCuts(3);

  CglGomory generator2;
  // try larger limit
  generator2.setLimit(300);

  CglKnapsackCover generator3;

  //CglOddHole generator4;
  //generator4.setMinimumViolation(0.005);
  //generator4.setMinimumViolationPer(0.00002);
  // try larger limit
  //generator4.setMaximumEntries(200);

  CglClique generator5;
  generator5.setStarCliqueReport(false);
  generator5.setRowCliqueReport(false);

  CglMixedIntegerRounding2 mixedGen;
  CglFlowCover flowGen;

  /*
      Add in generators. Do not override generators already installed.
    */
  int setting = cutsOnlyAtRoot_ ? -99 : -1;
  int numberGenerators = model.numberCutGenerators();
  int iGenerator;
  bool found;
  found = false;
  for (iGenerator = 0; iGenerator < numberGenerators; iGenerator++) {
    CglCutGenerator *generator = model.cutGenerator(iGenerator)->generator();
    CglProbing *cgl = dynamic_cast< CglProbing * >(generator);
    if (cgl) {
      found = true;
      break;
    }
  }
  if (!found && (genFlags & 1) != 0)
    model.addCutGenerator(&generator1, setting, "Probing");
  found = false;

  for (iGenerator = 0; iGenerator < numberGenerators; iGenerator++) {
    CglCutGenerator *generator = model.cutGenerator(iGenerator)->generator();
    CglGomory *cgl = dynamic_cast< CglGomory * >(generator);
    if (cgl) {
      found = true;
      break;
    }
  }
  if (!found && (genFlags & 2) != 0)
    model.addCutGenerator(&generator2, setting, "Gomory");

  found = false;
  for (iGenerator = 0; iGenerator < numberGenerators; iGenerator++) {
    CglCutGenerator *generator = model.cutGenerator(iGenerator)->generator();
    CglKnapsackCover *cgl = dynamic_cast< CglKnapsackCover * >(generator);
    if (cgl) {
      found = true;
      break;
    }
  }
  if (!found && (genFlags & 4) != 0)
    model.addCutGenerator(&generator3, setting, "Knapsack");
  //model.addCutGenerator(&generator4,setting,"OddHole");

  found = false;
  for (iGenerator = 0; iGenerator < numberGenerators; iGenerator++) {
    CglCutGenerator *generator = model.cutGenerator(iGenerator)->generator();
    CglClique *cgl = dynamic_cast< CglClique * >(generator);
    if (cgl) {
      found = true;
      break;
    }
  }
  if (!found && (genFlags & 8) != 0)
    model.addCutGenerator(&generator5, setting, "Clique");

  found = false;
  for (iGenerator = 0; iGenerator < numberGenerators; iGenerator++) {
    CglCutGenerator *generator = model.cutGenerator(iGenerator)->generator();
    CglFlowCover *cgl = dynamic_cast< CglFlowCover * >(generator);
    if (cgl) {
      found = true;
      break;
    }
  }
  if (!found && (genFlags & 16) != 0)
    model.addCutGenerator(&flowGen, setting, "FlowCover");

  found = false;
  for (iGenerator = 0; iGenerator < numberGenerators; iGenerator++) {
    CglCutGenerator *generator = model.cutGenerator(iGenerator)->generator();
    CglMixedIntegerRounding2 *cgl = dynamic_cast< CglMixedIntegerRounding2 * >(generator);
    if (cgl) {
      found = true;
      break;
    }
  }
  if (!found && (genFlags & 32) != 0)
    model.addCutGenerator(&mixedGen, setting, "MixedIntegerRounding2");

  // Say we want timings
  int newNumberGenerators = model.numberCutGenerators();
  for (iGenerator = numberGenerators; iGenerator < newNumberGenerators; iGenerator++) {
    CbcCutGenerator *generator = model.cutGenerator(iGenerator);
    generator->setTiming(true);
  }

  // Caution! Undocumented magic numbers.
  int currentPasses = model.getMaximumCutPassesAtRoot();
  if (currentPasses >= 0) {
    if (model.getNumCols() < 5000)
      model.setMaximumCutPassesAtRoot(CoinMax(50, currentPasses)); // use minimum drop
    else
      model.setMaximumCutPassesAtRoot(CoinMax(20, currentPasses));
  } else {
    currentPasses = -currentPasses;
    if (model.getNumCols() < 500)
      model.setMaximumCutPassesAtRoot(-CoinMax(100, currentPasses)); // always do 100 if possible
    else
      model.setMaximumCutPassesAtRoot(-CoinMax(20, currentPasses));
  }
}
// Setup heuristics
void CbcStrategyDefault::setupHeuristics(CbcModel &model)
{
  // Allow rounding heuristic

  CbcRounding heuristic1(model);
  heuristic1.setHeuristicName("rounding");
  int numberHeuristics = model.numberHeuristics();
  int iHeuristic;
  bool found;
  found = false;
  for (iHeuristic = 0; iHeuristic < numberHeuristics; iHeuristic++) {
    CbcHeuristic *heuristic = model.heuristic(iHeuristic);
    CbcRounding *cgl = dynamic_cast< CbcRounding * >(heuristic);
    if (cgl) {
      found = true;
      break;
    }
  }
  if (!found)
    model.addHeuristic(&heuristic1);
#ifdef JJF_ZERO
  // Allow join solutions
  CbcHeuristicLocal heuristic2(model);
  heuristic2.setHeuristicName("join solutions");
  heuristic2.setSearchType(1);
  found = false;
  for (iHeuristic = 0; iHeuristic < numberHeuristics; iHeuristic++) {
    CbcHeuristic *heuristic = model.heuristic(iHeuristic);
    CbcHeuristicLocal *cgl = dynamic_cast< CbcHeuristicLocal * >(heuristic);
    if (cgl) {
      found = true;
      break;
    }
  }
  if (!found)
    model.addHeuristic(&heuristic2);
#endif
}
// Do printing stuff
void CbcStrategyDefault::setupPrinting(CbcModel &model, int modelLogLevel)
{
  if (!modelLogLevel) {
    model.solver()->setHintParam(OsiDoReducePrint, true, OsiHintTry);
    model.messageHandler()->setLogLevel(0);
    model.solver()->messageHandler()->setLogLevel(0);
  } else if (modelLogLevel == 1) {
    model.solver()->setHintParam(OsiDoReducePrint, true, OsiHintTry);
    model.messageHandler()->setLogLevel(1);
    model.solver()->messageHandler()->setLogLevel(0);
  } else {
    model.messageHandler()->setLogLevel(CoinMax(2, model.messageHandler()->logLevel()));
    model.solver()->messageHandler()->setLogLevel(CoinMax(1, model.solver()->messageHandler()->logLevel()));
    model.setPrintFrequency(CoinMin(50, model.printFrequency()));
  }
}

/*
 Aside from setting CbcModel::numberStrong_ and numberBeforeTrust, the big
 activity is integer preprocessing. Surely this code to do preprocessing
 duplicates code to do preprocessing up in the solver main routine. Most of the
 effort goes into manipulating SOS sets.
*/
// Other stuff e.g. strong branching
void CbcStrategyDefault::setupOther(CbcModel &model)
{
  // See if preprocessing wanted
  if (desiredPreProcess_) {
    delete process_;
    /*
          Inaccurate as of 080122 --- assignSolver (below) can now be instructed not to
          delete the existing solver when the preprocessed solver is assigned to the
          model. 'Course, we do need to hold on to a pointer somewhere, and that must
          be captured before this call.
        */
    // solver_ should have been cloned outside
    CglPreProcess *process = new CglPreProcess();
    // Pass in models message handler
    process->passInMessageHandler(model.messageHandler());
    OsiSolverInterface *solver = model.solver();
#ifdef COIN_HAS_CLP
    OsiClpSolverInterface *clpSolver = dynamic_cast< OsiClpSolverInterface * >(solver);
    if (clpSolver && false) {
      // see if all coefficients multiple of 0.01 (close enough)
      CoinPackedMatrix *matrix = clpSolver->getModelPtr()->matrix();
      double *element = matrix->getMutableElements();
      //const int * row = matrix->getIndices();
      const CoinBigIndex *columnStart = matrix->getVectorStarts();
      const int *columnLength = matrix->getVectorLengths();
#ifdef COIN_DETAIL
      int numberInt = 0;
#endif
      int numberNon = 0;
      int numberClose = 0;
      int numberColumns = clpSolver->getNumCols();
      int iColumn;
      for (iColumn = 0; iColumn < numberColumns; iColumn++) {
        for (CoinBigIndex j = columnStart[iColumn];
             j < columnStart[iColumn] + columnLength[iColumn]; j++) {
          //int iRow = row[j];
          double value1 = element[j];
          double value = fabs(value1);
          if (value > 1.0e7) {
            if (value != floor(value))
              numberNon++;
#ifdef COIN_DETAIL
            else
              numberInt++;
#endif
          } else {
            int iValue = static_cast< int >(100 * (value + 0.005));
            double value2 = iValue;
            if (value2 == 100.0 * value) {
#ifdef COIN_DETAIL
              numberInt++;
#endif
            } else if (fabs(value2 - 100.0 * value) < 1.0e-5) {
              numberClose++;
            } else {
              numberNon++;
            }
          }
        }
      }
      if (!numberNon && numberClose) {
        COIN_DETAIL_PRINT(printf("Tidying %d multiples of 0.01, %d close\n",
          numberInt, numberClose));
        for (iColumn = 0; iColumn < numberColumns; iColumn++) {
          for (CoinBigIndex j = columnStart[iColumn];
               j < columnStart[iColumn] + columnLength[iColumn]; j++) {
            //int iRow = row[j];
            double value1 = element[j];
            double value = fabs(value1);
            if (value < 1.0e7) {
              int iValue = static_cast< int >(100 * (value + 0.005));
              double value2 = iValue;
              if (value2 != 100.0 * value) {
                value2 *= 0.01;
                if (fabs(value - floor(value + 0.5)) <= 1.0e-7)
                  value2 = floor(value + 0.5);
                if (value1 < 0.0)
                  value2 = -value2;
                element[j] = value2;
              }
            }
          }
        }
      }
    }
#endif
    {
      // mark some columns as ineligible for presolve
      int numberColumns = solver->getNumCols();
      char *prohibited = new char[numberColumns];
      memset(prohibited, 0, numberColumns);
      int numberProhibited = 0;
      /*
              Create CbcSimpleInteger objects would be more accurate in the general
              case.  The `false' parameter says we won't delete existing objects.

              Only Clp will produce SOS objects in findIntegers (080122), and that's
              where a possible conversion can occur. If clp is holding OsiSOS objects,
              they'll be converted to CbcSOS objects.
            */
      // convert to Cbc integers
      model.findIntegers(false);
      int numberObjects = model.numberObjects();
      if (numberObjects) {
        OsiObject **objects = model.objects();
        for (int iObject = 0; iObject < numberObjects; iObject++) {
          CbcSOS *obj = dynamic_cast< CbcSOS * >(objects[iObject]);
          if (obj) {
            // SOS
            int n = obj->numberMembers();
            const int *which = obj->members();
            for (int i = 0; i < n; i++) {
              int iColumn = which[i];
              prohibited[iColumn] = 1;
              numberProhibited++;
            }
          }
        }
      }
      if (numberProhibited)
        process->passInProhibited(prohibited, numberColumns);
      delete[] prohibited;
    }
    int logLevel = model.messageHandler()->logLevel();
#ifdef COIN_HAS_CLP
    //OsiClpSolverInterface * clpSolver = dynamic_cast< OsiClpSolverInterface*> (solver);
    ClpSimplex *lpSolver = NULL;
    if (clpSolver) {
      if (clpSolver->messageHandler()->logLevel())
        clpSolver->messageHandler()->setLogLevel(1);
      if (logLevel > -1)
        clpSolver->messageHandler()->setLogLevel(CoinMin(logLevel, clpSolver->messageHandler()->logLevel()));
      lpSolver = clpSolver->getModelPtr();
      /// If user left factorization frequency then compute
      lpSolver->defaultFactorizationFrequency();
    }
#endif
    // Tell solver we are in Branch and Cut
    solver->setHintParam(OsiDoInBranchAndCut, true, OsiHintDo);
    // Default set of cut generators
    // Limited set that could reduce problem size (drop rows / fix values)
    CglProbing generator1;
    generator1.setUsingObjective(true);
    generator1.setMaxPass(1);
    generator1.setMaxPassRoot(1);
    generator1.setMaxProbeRoot(CoinMin(3000, solver->getNumCols()));
    generator1.setMaxProbeRoot(123);
    generator1.setMaxElements(100);
    generator1.setMaxElementsRoot(200);
    generator1.setMaxLookRoot(50);
    generator1.setRowCuts(3);
    //generator1.messageHandler()->setLogLevel(logLevel);
    // Not needed with pass in process->messageHandler()->setLogLevel(logLevel);
    // Add in generators
    process->addCutGenerator(&generator1);
    int translate[] = { 9999, 0, 2, -2, 3, 4, 4, 4 };
    OsiSolverInterface *solver2 = process->preProcessNonDefault(*solver,
      translate[desiredPreProcess_], preProcessPasses_, 6);
    // Tell solver we are not in Branch and Cut
    solver->setHintParam(OsiDoInBranchAndCut, false, OsiHintDo);
    if (solver2)
      solver2->setHintParam(OsiDoInBranchAndCut, false, OsiHintDo);
    bool feasible = true;
    if (!solver2) {
      feasible = false;
      //printf("Pre-processing says infeasible\n");
      delete process;
      preProcessState_ = -1;
      process_ = NULL;
    } else {
      // now tighten bounds
#ifdef COIN_HAS_CLP
      if (clpSolver) {
        // model has changed
        solver = model.solver();
        OsiClpSolverInterface *clpSolver = dynamic_cast< OsiClpSolverInterface * >(solver);
        ClpSimplex *lpSolver = clpSolver->getModelPtr();
        lpSolver->passInMessageHandler(solver->messageHandler());
        if (lpSolver->tightenPrimalBounds() == 0) {
          lpSolver->dual();
        } else {
          feasible = false;
        }
      }
#endif
      if (feasible) {
        preProcessState_ = 1;
        process_ = process;
        /* Note that original solver will be kept (with false)
                   and that final solver will also be kept.
                   This is for post-processing

		   Keep in mind when examining this that linear presolve does not
		   understand SOS.
                */
        OsiSolverInterface *solver3 = solver2->clone();
        model.assignSolver(solver3, false);
        if (process_->numberSOS()) {
          int numberSOS = process_->numberSOS();
          int numberIntegers = model.numberIntegers();
          /* model may not have created objects
                       If none then create
                       NOTE - put back to original column numbers as
                       CbcModel will pack down ALL as it doesn't know where from
                    */
          bool someObjects = model.numberObjects() > 0;
          if (!numberIntegers || !model.numberObjects()) {
            model.findIntegers(true);
            numberIntegers = model.numberIntegers();
          }
          OsiObject **oldObjects = model.objects();
          // Do sets and priorities
          OsiObject **objects = new OsiObject *[numberSOS];
          // set old objects to have low priority
          int numberOldObjects = model.numberObjects();
          int numberColumns = model.getNumCols();
          for (int iObj = 0; iObj < numberOldObjects; iObj++) {
            int oldPriority = oldObjects[iObj]->priority();
            oldObjects[iObj]->setPriority(numberColumns + oldPriority);
          }
          const int *starts = process_->startSOS();
          const int *which = process_->whichSOS();
          const int *type = process_->typeSOS();
          const double *weight = process_->weightSOS();
          int iSOS;
          for (iSOS = 0; iSOS < numberSOS; iSOS++) {
            int iStart = starts[iSOS];
            int n = starts[iSOS + 1] - iStart;
            objects[iSOS] = new CbcSOS(&model, n, which + iStart, weight + iStart,
              iSOS, type[iSOS]);
            // branch on long sets first
            objects[iSOS]->setPriority(numberColumns - n);
          }
          model.addObjects(numberSOS, objects);
          for (iSOS = 0; iSOS < numberSOS; iSOS++)
            delete objects[iSOS];
          delete[] objects;
          if (!someObjects) {
            // put back old column numbers
            const int *originalColumns = process_->originalColumns();
            // use reverse lookup to fake it
            int n = originalColumns[numberColumns - 1] + 1;
            int *fake = new int[n];
            int i;
            // This was wrong (now is correct) - so could never have been called
            abort();
            for (i = 0; i < n; i++)
              fake[i] = -1;
            for (i = 0; i < numberColumns; i++)
              fake[originalColumns[i]] = i;
            for (int iObject = 0; iObject < model.numberObjects(); iObject++) {
              // redo ids etc
              CbcSimpleInteger *obj = dynamic_cast< CbcSimpleInteger * >(model.modifiableObject(iObject));
              if (obj) {
                obj->resetSequenceEtc(n, fake);
              } else {
                // redo ids etc
                CbcObject *obj = dynamic_cast< CbcObject * >(model.modifiableObject(iObject));
                assert(obj);
                obj->redoSequenceEtc(&model, n, fake);
              }
            }
            delete[] fake;
          }
        }
      } else {
        //printf("Pre-processing says infeasible\n");
        delete process;
        preProcessState_ = -1;
        process_ = NULL;
      }
    }
  }
  model.setNumberStrong(numberStrong_);
  model.setNumberBeforeTrust(numberBeforeTrust_);
}

// Create C++ lines to get to current state
void CbcStrategyDefault::generateCpp(FILE *fp)
{
  fprintf(fp, "0#include \"CbcStrategy.hpp\"\n");
  fprintf(fp, "3  CbcStrategyDefault strategy(%s,%d,%d,%d);\n",
    cutsOnlyAtRoot_ ? "1" : "0",
    numberStrong_,
    numberBeforeTrust_,
    printLevel_);
  fprintf(fp, "3  strategy.setupPreProcessing(%d,%d);\n",
    desiredPreProcess_, preProcessPasses_);
}
// Default Constructor
CbcStrategyDefaultSubTree::CbcStrategyDefaultSubTree(CbcModel *parent,
  int cutsOnlyAtRoot,
  int numberStrong,
  int numberBeforeTrust,
  int printLevel)
  : CbcStrategy()
  , parentModel_(parent)
  , cutsOnlyAtRoot_(cutsOnlyAtRoot)
  , numberStrong_(numberStrong)
  , numberBeforeTrust_(numberBeforeTrust)
  , printLevel_(printLevel)
{
}

// Destructor
CbcStrategyDefaultSubTree::~CbcStrategyDefaultSubTree()
{
}

// Clone
CbcStrategy *
CbcStrategyDefaultSubTree::clone() const
{
  return new CbcStrategyDefaultSubTree(*this);
}

// Copy constructor
CbcStrategyDefaultSubTree::CbcStrategyDefaultSubTree(const CbcStrategyDefaultSubTree &rhs)
  : CbcStrategy(rhs)
  , parentModel_(rhs.parentModel_)
  , cutsOnlyAtRoot_(rhs.cutsOnlyAtRoot_)
  , numberStrong_(rhs.numberStrong_)
  , numberBeforeTrust_(rhs.numberBeforeTrust_)
  , printLevel_(rhs.printLevel_)
{
  setNested(rhs.getNested());
}

// Setup cut generators
void CbcStrategyDefaultSubTree::setupCutGenerators(CbcModel &model)
{
  // Set up some cut generators and defaults
  if (cutsOnlyAtRoot_ < 0)
    return; // no cuts wanted
  // Probing first as gets tight bounds on continuous

  CglProbing generator1;
  generator1.setUsingObjective(true);
  generator1.setMaxPass(1);
  // Number of unsatisfied variables to look at
  generator1.setMaxProbe(10);
  // How far to follow the consequences
  generator1.setMaxLook(10);
  // Only look at rows with fewer than this number of elements
  generator1.setMaxElements(200);
  //generator1.setRowCuts(3);

  CglGomory generator2;
  // try larger limit
  generator2.setLimit(300);

  CglKnapsackCover generator3;

  //CglOddHole generator4;
  //generator4.setMinimumViolation(0.005);
  //generator4.setMinimumViolationPer(0.00002);
  // try larger limit
  //generator4.setMaximumEntries(200);

  CglClique generator5;
  generator5.setStarCliqueReport(false);
  generator5.setRowCliqueReport(false);

  CglMixedIntegerRounding2 mixedGen;
  CglFlowCover flowGen;

  // Add in generators
  int setting = cutsOnlyAtRoot_ ? -99 : -1;
  int numberGenerators = model.numberCutGenerators();
  int numberParentGenerators = parentModel_->numberCutGenerators();
  int iGenerator;
  bool found;
  found = false;
  int howOften = 0;
  for (iGenerator = 0; iGenerator < numberParentGenerators; iGenerator++) {
    CglCutGenerator *generator = parentModel_->cutGenerator(iGenerator)->generator();
    CglProbing *cgl = dynamic_cast< CglProbing * >(generator);
    if (cgl) {
      found = true;
      howOften = parentModel_->cutGenerator(iGenerator)->howOften();
      break;
    }
  }

  if (found && (howOften >= -1 || howOften == -98)) {
    found = false;
    for (iGenerator = 0; iGenerator < numberGenerators; iGenerator++) {
      CglCutGenerator *generator = model.cutGenerator(iGenerator)->generator();
      CglProbing *cgl = dynamic_cast< CglProbing * >(generator);
      if (cgl) {
        found = true;
        break;
      }
    }
    if (!found) {
      if (howOften == -1)
        howOften = -98;
      else if (howOften == -98)
        howOften = -99;
      model.addCutGenerator(&generator1, setting, "Probing");
      CbcCutGenerator *generator = model.cutGenerator(numberGenerators);
      generator->setHowOften(howOften);
      numberGenerators++;
    }
  }
  found = false;
  for (iGenerator = 0; iGenerator < numberParentGenerators; iGenerator++) {
    CglCutGenerator *generator = parentModel_->cutGenerator(iGenerator)->generator();
    CglGomory *cgl = dynamic_cast< CglGomory * >(generator);
    if (cgl) {
      found = true;
      howOften = parentModel_->cutGenerator(iGenerator)->howOften();
      break;
    }
  }
  if (found && howOften >= 0) {
    found = false;
    for (iGenerator = 0; iGenerator < numberGenerators; iGenerator++) {
      CglCutGenerator *generator = model.cutGenerator(iGenerator)->generator();
      CglGomory *cgl = dynamic_cast< CglGomory * >(generator);
      if (cgl) {
        found = true;
        break;
      }
    }
    if (!found)
      model.addCutGenerator(&generator2, setting, "Gomory");
  }
  found = false;
  for (iGenerator = 0; iGenerator < numberParentGenerators; iGenerator++) {
    CglCutGenerator *generator = parentModel_->cutGenerator(iGenerator)->generator();
    CglKnapsackCover *cgl = dynamic_cast< CglKnapsackCover * >(generator);
    if (cgl) {
      found = true;
      howOften = parentModel_->cutGenerator(iGenerator)->howOften();
      break;
    }
  }
  if (found && howOften >= 0) {
    found = false;
    for (iGenerator = 0; iGenerator < numberGenerators; iGenerator++) {
      CglCutGenerator *generator = model.cutGenerator(iGenerator)->generator();
      CglKnapsackCover *cgl = dynamic_cast< CglKnapsackCover * >(generator);
      if (cgl) {
        found = true;
        break;
      }
    }
    if (!found)
      model.addCutGenerator(&generator3, setting, "Knapsack");
  }
  found = false;
  for (iGenerator = 0; iGenerator < numberParentGenerators; iGenerator++) {
    CglCutGenerator *generator = parentModel_->cutGenerator(iGenerator)->generator();
    CglClique *cgl = dynamic_cast< CglClique * >(generator);
    if (cgl) {
      found = true;
      howOften = parentModel_->cutGenerator(iGenerator)->howOften();
      break;
    }
  }
  if (found && howOften >= 0) {
    found = false;
    for (iGenerator = 0; iGenerator < numberGenerators; iGenerator++) {
      CglCutGenerator *generator = model.cutGenerator(iGenerator)->generator();
      CglClique *cgl = dynamic_cast< CglClique * >(generator);
      if (cgl) {
        found = true;
        break;
      }
    }
    if (!found)
      model.addCutGenerator(&generator5, setting, "Clique");
  }
  found = false;
  for (iGenerator = 0; iGenerator < numberParentGenerators; iGenerator++) {
    CglCutGenerator *generator = parentModel_->cutGenerator(iGenerator)->generator();
    CglFlowCover *cgl = dynamic_cast< CglFlowCover * >(generator);
    if (cgl) {
      found = true;
      howOften = parentModel_->cutGenerator(iGenerator)->howOften();
      break;
    }
  }
  if (found && howOften >= 0) {
    found = false;
    for (iGenerator = 0; iGenerator < numberGenerators; iGenerator++) {
      CglCutGenerator *generator = model.cutGenerator(iGenerator)->generator();
      CglFlowCover *cgl = dynamic_cast< CglFlowCover * >(generator);
      if (cgl) {
        found = true;
        break;
      }
    }
    if (!found)
      model.addCutGenerator(&flowGen, setting, "FlowCover");
    found = false;
  }
  for (iGenerator = 0; iGenerator < numberParentGenerators; iGenerator++) {
    CglCutGenerator *generator = parentModel_->cutGenerator(iGenerator)->generator();
    CglMixedIntegerRounding2 *cgl = dynamic_cast< CglMixedIntegerRounding2 * >(generator);
    if (cgl) {
      found = true;
      howOften = parentModel_->cutGenerator(iGenerator)->howOften();
      break;
    }
  }
  if (found && howOften >= 0) {
    found = false;
    for (iGenerator = 0; iGenerator < numberGenerators; iGenerator++) {
      CglCutGenerator *generator = model.cutGenerator(iGenerator)->generator();
      CglMixedIntegerRounding2 *cgl = dynamic_cast< CglMixedIntegerRounding2 * >(generator);
      if (cgl) {
        found = true;
        break;
      }
    }
    if (!found)
      model.addCutGenerator(&mixedGen, setting, "MixedIntegerRounding2");
  }
#ifdef JJF_ZERO
  // Say we want timings
  int newNumberGenerators = model.numberCutGenerators();
  for (iGenerator = numberGenerators; iGenerator < newNumberGenerators; iGenerator++) {
    CbcCutGenerator *generator = model.cutGenerator(iGenerator);
    generator->setTiming(true);
  }
#endif
  if (model.getNumCols() < -500)
    model.setMaximumCutPassesAtRoot(-100); // always do 100 if possible
  else if (model.getNumCols() < 5000)
    model.setMaximumCutPassesAtRoot(100); // use minimum drop
  else
    model.setMaximumCutPassesAtRoot(20);
}
// Setup heuristics
void CbcStrategyDefaultSubTree::setupHeuristics(CbcModel &model)
{
  // Allow rounding heuristic

  CbcRounding heuristic1(model);
  heuristic1.setHeuristicName("rounding");
  int numberHeuristics = model.numberHeuristics();
  int iHeuristic;
  bool found;
  found = false;
  for (iHeuristic = 0; iHeuristic < numberHeuristics; iHeuristic++) {
    CbcHeuristic *heuristic = model.heuristic(iHeuristic);
    CbcRounding *cgl = dynamic_cast< CbcRounding * >(heuristic);
    if (cgl) {
      found = true;
      break;
    }
  }
  if (!found)
    model.addHeuristic(&heuristic1);
  if ((model.moreSpecialOptions() & 32768) != 0) {
    // Allow join solutions
    CbcHeuristicLocal heuristic2(model);
    heuristic2.setHeuristicName("join solutions");
    //sheuristic2.setSearchType(1);
    found = false;
    for (iHeuristic = 0; iHeuristic < numberHeuristics; iHeuristic++) {
      CbcHeuristic *heuristic = model.heuristic(iHeuristic);
      CbcHeuristicLocal *cgl = dynamic_cast< CbcHeuristicLocal * >(heuristic);
      if (cgl) {
        found = true;
        break;
      }
    }
    if (!found)
      model.addHeuristic(&heuristic2);
    // Allow RINS
    CbcHeuristicRINS heuristic5(model);
    heuristic5.setHeuristicName("RINS");
    heuristic5.setFractionSmall(0.5);
    heuristic5.setDecayFactor(5.0);
    //heuristic5.setSearchType(1);
    found = false;
    for (iHeuristic = 0; iHeuristic < numberHeuristics; iHeuristic++) {
      CbcHeuristic *heuristic = model.heuristic(iHeuristic);
      CbcHeuristicLocal *cgl = dynamic_cast< CbcHeuristicLocal * >(heuristic);
      if (cgl) {
        found = true;
        break;
      }
    }
    if (!found)
      model.addHeuristic(&heuristic5);
  }
}
// Do printing stuff
void CbcStrategyDefaultSubTree::setupPrinting(CbcModel &model, int modelLogLevel)
{
  if (!modelLogLevel) {
    model.solver()->setHintParam(OsiDoReducePrint, true, OsiHintTry);
    model.messageHandler()->setLogLevel(0);
    model.solver()->messageHandler()->setLogLevel(0);
  } else if (modelLogLevel == 1) {
    model.solver()->setHintParam(OsiDoReducePrint, true, OsiHintTry);
    model.messageHandler()->setLogLevel(1);
    model.solver()->messageHandler()->setLogLevel(0);
  } else {
    model.messageHandler()->setLogLevel(2);
    model.solver()->messageHandler()->setLogLevel(1);
    model.setPrintFrequency(50);
  }
}
// Other stuff e.g. strong branching
void CbcStrategyDefaultSubTree::setupOther(CbcModel &model)
{
  model.setNumberStrong(numberStrong_);
  model.setNumberBeforeTrust(numberBeforeTrust_);
}
// For uniform setting of cut and heuristic options
void setCutAndHeuristicOptions(CbcModel &model)
{
  int numberGenerators = model.numberCutGenerators();
  int iGenerator;
  for (iGenerator = 0; iGenerator < numberGenerators; iGenerator++) {
    CglCutGenerator *generator = model.cutGenerator(iGenerator)->generator();
    CglProbing *cglProbing = dynamic_cast< CglProbing * >(generator);
    if (cglProbing) {
      cglProbing->setUsingObjective(1);
      cglProbing->setMaxPass(1);
      cglProbing->setMaxPassRoot(1);
      // Number of unsatisfied variables to look at
      cglProbing->setMaxProbe(10);
      cglProbing->setMaxProbeRoot(50);
      //cglProbing->setMaxProbeRoot(123);
      // How far to follow the consequences
      cglProbing->setMaxLook(5);
      cglProbing->setMaxLookRoot(50);
      cglProbing->setMaxLookRoot(10);
      // Only look at rows with fewer than this number of elements
      cglProbing->setMaxElements(200);
      cglProbing->setMaxElementsRoot(300);
      cglProbing->setRowCuts(3);
    }
#ifdef JJF_ZERO
    CglGomory *cglGomory = dynamic_cast< CglGomory * >(generator);
    if (cglGomory) {
      // try larger limit
      cglGomory->setLimitAtRoot(1000);
      cglGomory->setLimit(50);
    }
    CglKnapsackCover *cglKnapsackCover = dynamic_cast< CglKnapsackCover * >(generator);
    if (cglKnapsackCover) {
    }
#endif
  }
}

/* vi: softtabstop=2 shiftwidth=2 expandtab tabstop=2
*/
