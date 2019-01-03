/* $Id$ */
// Copyright (C) 2004, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

#include "CbcModel.hpp"
#include "CbcNode.hpp"
#include "CbcTreeLocal.hpp"
#include "CoinPackedMatrix.hpp"
#include "CoinTime.hpp"
#include "OsiRowCutDebugger.hpp"
#include <cassert>
#ifdef JJF_ZERO
// gdb doesn't always put breakpoints in this virtual function
// just stick xxxxxx() where you want to start
static void xxxxxx()
{
  printf("break\n");
}
#endif
CbcTreeLocal::CbcTreeLocal()
  : localNode_(NULL)
  , bestSolution_(NULL)
  , savedSolution_(NULL)
  , saveNumberSolutions_(0)
  , model_(NULL)
  , originalLower_(NULL)
  , originalUpper_(NULL)
  , range_(0)
  , typeCuts_(-1)
  , maxDiversification_(0)
  , diversification_(0)
  , nextStrong_(false)
  , rhs_(0.0)
  , savedGap_(0.0)
  , bestCutoff_(0.0)
  , timeLimit_(0)
  , startTime_(0)
  , nodeLimit_(0)
  , startNode_(-1)
  , searchType_(-1)
  , refine_(false)
{
}
/* Constructor with solution.
   range is upper bound on difference from given solution.
   maxDiversification is maximum number of diversifications to try
   timeLimit is seconds in subTree
   nodeLimit is nodes in subTree
*/
CbcTreeLocal::CbcTreeLocal(CbcModel *model, const double *solution,
  int range, int typeCuts, int maxDiversification,
  int timeLimit, int nodeLimit, bool refine)
  : localNode_(NULL)
  , bestSolution_(NULL)
  , savedSolution_(NULL)
  , saveNumberSolutions_(0)
  , model_(model)
  , originalLower_(NULL)
  , originalUpper_(NULL)
  , range_(range)
  , typeCuts_(typeCuts)
  , maxDiversification_(maxDiversification)
  , diversification_(0)
  , nextStrong_(false)
  , rhs_(0.0)
  , savedGap_(0.0)
  , bestCutoff_(0.0)
  , timeLimit_(timeLimit)
  , startTime_(0)
  , nodeLimit_(nodeLimit)
  , startNode_(-1)
  , searchType_(-1)
  , refine_(refine)
{

  OsiSolverInterface *solver = model_->solver();
  const double *lower = solver->getColLower();
  const double *upper = solver->getColUpper();
  //const double * solution = solver->getColSolution();
  //const double * objective = solver->getObjCoefficients();
  double primalTolerance;
  solver->getDblParam(OsiPrimalTolerance, primalTolerance);

  // Get increment
  model_->analyzeObjective();

  {
    // needed to sync cutoffs
    double value;
    solver->getDblParam(OsiDualObjectiveLimit, value);
    model_->setCutoff(value * solver->getObjSense());
  }
  bestCutoff_ = model_->getCutoff();
  // save current gap
  savedGap_ = model_->getDblParam(CbcModel::CbcAllowableGap);

  // make sure integers found
  model_->findIntegers(false);
  int numberIntegers = model_->numberIntegers();
  const int *integerVariable = model_->integerVariable();
  int i;
  double direction = solver->getObjSense();
  double newSolutionValue = 1.0e50;
  if (solution) {
    // copy solution
    solver->setColSolution(solution);
    newSolutionValue = direction * solver->getObjValue();
  }
  originalLower_ = new double[numberIntegers];
  originalUpper_ = new double[numberIntegers];
  bool all01 = true;
  int number01 = 0;
  for (i = 0; i < numberIntegers; i++) {
    int iColumn = integerVariable[i];
    originalLower_[i] = lower[iColumn];
    originalUpper_[i] = upper[iColumn];
    if (upper[iColumn] - lower[iColumn] > 1.5)
      all01 = false;
    else if (upper[iColumn] - lower[iColumn] == 1.0)
      number01++;
  }
  if (all01 && !typeCuts_)
    typeCuts_ = 1; // may as well so we don't have to deal with refine
  if (!number01 && !typeCuts_) {
    if (model_->messageHandler()->logLevel() > 1)
      printf("** No 0-1 variables and local search only on 0-1 - switching off\n");
    typeCuts_ = -1;
  } else {
    if (model_->messageHandler()->logLevel() > 1) {
      std::string type;
      if (all01) {
        printf("%d 0-1 variables normal local  cuts\n",
          number01);
      } else if (typeCuts_) {
        printf("%d 0-1 variables, %d other - general integer local cuts\n",
          number01, numberIntegers - number01);
      } else {
        printf("%d 0-1 variables, %d other - local cuts but just on 0-1 variables\n",
          number01, numberIntegers - number01);
      }
      printf("maximum diversifications %d, initial cutspace %d, max time %d seconds, max nodes %d\n",
        maxDiversification_, range_, timeLimit_, nodeLimit_);
    }
  }
  int numberColumns = model_->getNumCols();
  savedSolution_ = new double[numberColumns];
  memset(savedSolution_, 0, numberColumns * sizeof(double));
  if (solution) {
    rhs_ = range_;
    // Check feasible
    int goodSolution = createCut(solution, cut_);
    if (goodSolution >= 0) {
      for (i = 0; i < numberIntegers; i++) {
        int iColumn = integerVariable[i];
        double value = floor(solution[iColumn] + 0.5);
        // fix so setBestSolution will work
        solver->setColLower(iColumn, value);
        solver->setColUpper(iColumn, value);
      }
      model_->reserveCurrentSolution();
      // Create cut and get total gap
      if (newSolutionValue < bestCutoff_) {
        model_->setBestSolution(CBC_ROUNDING, newSolutionValue, solution);
        bestCutoff_ = model_->getCutoff();
        // save as best solution
        memcpy(savedSolution_, model_->bestSolution(), numberColumns * sizeof(double));
      }
      for (i = 0; i < numberIntegers; i++) {
        int iColumn = integerVariable[i];
        // restore bounds
        solver->setColLower(iColumn, originalLower_[i]);
        solver->setColUpper(iColumn, originalUpper_[i]);
      }
      // make sure can't stop on gap
      model_->setDblParam(CbcModel::CbcAllowableGap, -1.0e50);
    } else {
      model_ = NULL;
    }
  } else {
    // no solution
    rhs_ = 1.0e50;
    // make sure can't stop on gap
    model_->setDblParam(CbcModel::CbcAllowableGap, -1.0e50);
  }
}
CbcTreeLocal::~CbcTreeLocal()
{
  delete[] originalLower_;
  delete[] originalUpper_;
  delete[] bestSolution_;
  delete[] savedSolution_;
  delete localNode_;
}
// Copy constructor
CbcTreeLocal::CbcTreeLocal(const CbcTreeLocal &rhs)
  : CbcTree(rhs)
  , saveNumberSolutions_(rhs.saveNumberSolutions_)
  , model_(rhs.model_)
  , range_(rhs.range_)
  , typeCuts_(rhs.typeCuts_)
  , maxDiversification_(rhs.maxDiversification_)
  , diversification_(rhs.diversification_)
  , nextStrong_(rhs.nextStrong_)
  , rhs_(rhs.rhs_)
  , savedGap_(rhs.savedGap_)
  , bestCutoff_(rhs.bestCutoff_)
  , timeLimit_(rhs.timeLimit_)
  , startTime_(rhs.startTime_)
  , nodeLimit_(rhs.nodeLimit_)
  , startNode_(rhs.startNode_)
  , searchType_(rhs.searchType_)
  , refine_(rhs.refine_)
{
  cut_ = rhs.cut_;
  fixedCut_ = rhs.fixedCut_;
  if (rhs.localNode_)
    localNode_ = new CbcNode(*rhs.localNode_);
  else
    localNode_ = NULL;
  if (rhs.originalLower_) {
    int numberIntegers = model_->numberIntegers();
    originalLower_ = new double[numberIntegers];
    memcpy(originalLower_, rhs.originalLower_, numberIntegers * sizeof(double));
    originalUpper_ = new double[numberIntegers];
    memcpy(originalUpper_, rhs.originalUpper_, numberIntegers * sizeof(double));
  } else {
    originalLower_ = NULL;
    originalUpper_ = NULL;
  }
  if (rhs.bestSolution_) {
    int numberColumns = model_->getNumCols();
    bestSolution_ = new double[numberColumns];
    memcpy(bestSolution_, rhs.bestSolution_, numberColumns * sizeof(double));
  } else {
    bestSolution_ = NULL;
  }
  if (rhs.savedSolution_) {
    int numberColumns = model_->getNumCols();
    savedSolution_ = new double[numberColumns];
    memcpy(savedSolution_, rhs.savedSolution_, numberColumns * sizeof(double));
  } else {
    savedSolution_ = NULL;
  }
}
//----------------------------------------------------------------
// Assignment operator
//-------------------------------------------------------------------
CbcTreeLocal &
CbcTreeLocal::operator=(const CbcTreeLocal &rhs)
{
  if (this != &rhs) {
    CbcTree::operator=(rhs);
    saveNumberSolutions_ = rhs.saveNumberSolutions_;
    cut_ = rhs.cut_;
    fixedCut_ = rhs.fixedCut_;
    delete localNode_;
    if (rhs.localNode_)
      localNode_ = new CbcNode(*rhs.localNode_);
    else
      localNode_ = NULL;
    model_ = rhs.model_;
    range_ = rhs.range_;
    typeCuts_ = rhs.typeCuts_;
    maxDiversification_ = rhs.maxDiversification_;
    diversification_ = rhs.diversification_;
    nextStrong_ = rhs.nextStrong_;
    rhs_ = rhs.rhs_;
    savedGap_ = rhs.savedGap_;
    bestCutoff_ = rhs.bestCutoff_;
    timeLimit_ = rhs.timeLimit_;
    startTime_ = rhs.startTime_;
    nodeLimit_ = rhs.nodeLimit_;
    startNode_ = rhs.startNode_;
    searchType_ = rhs.searchType_;
    refine_ = rhs.refine_;
    delete[] originalLower_;
    delete[] originalUpper_;
    if (rhs.originalLower_) {
      int numberIntegers = model_->numberIntegers();
      originalLower_ = new double[numberIntegers];
      memcpy(originalLower_, rhs.originalLower_, numberIntegers * sizeof(double));
      originalUpper_ = new double[numberIntegers];
      memcpy(originalUpper_, rhs.originalUpper_, numberIntegers * sizeof(double));
    } else {
      originalLower_ = NULL;
      originalUpper_ = NULL;
    }
    delete[] bestSolution_;
    if (rhs.bestSolution_) {
      int numberColumns = model_->getNumCols();
      bestSolution_ = new double[numberColumns];
      memcpy(bestSolution_, rhs.bestSolution_, numberColumns * sizeof(double));
    } else {
      bestSolution_ = NULL;
    }
    delete[] savedSolution_;
    if (rhs.savedSolution_) {
      int numberColumns = model_->getNumCols();
      savedSolution_ = new double[numberColumns];
      memcpy(savedSolution_, rhs.savedSolution_, numberColumns * sizeof(double));
    } else {
      savedSolution_ = NULL;
    }
  }
  return *this;
}
// Clone
CbcTree *
CbcTreeLocal::clone() const
{
  return new CbcTreeLocal(*this);
}
// Pass in solution (so can be used after heuristic)
void CbcTreeLocal::passInSolution(const double *solution, double solutionValue)
{
  int numberColumns = model_->getNumCols();
  delete[] savedSolution_;
  savedSolution_ = new double[numberColumns];
  memcpy(savedSolution_, solution, numberColumns * sizeof(double));
  rhs_ = range_;
  // Check feasible
  int goodSolution = createCut(solution, cut_);
  if (goodSolution >= 0) {
    bestCutoff_ = CoinMin(solutionValue, model_->getCutoff());
  } else {
    model_ = NULL;
  }
}
// Return the top node of the heap
CbcNode *
CbcTreeLocal::top() const
{
#ifdef CBC_DEBUG
  int smallest = 9999999;
  int largest = -1;
  double smallestD = 1.0e30;
  double largestD = -1.0e30;
  int n = nodes_.size();
  for (int i = 0; i < n; i++) {
    int nn = nodes_[i]->nodeInfo()->nodeNumber();
    double dd = nodes_[i]->objectiveValue();
    largest = CoinMax(largest, nn);
    smallest = CoinMin(smallest, nn);
    largestD = CoinMax(largestD, dd);
    smallestD = CoinMin(smallestD, dd);
  }
  if (model_->messageHandler()->logLevel() > 1) {
    printf("smallest %d, largest %d, top %d\n", smallest, largest,
      nodes_.front()->nodeInfo()->nodeNumber());
    printf("smallestD %g, largestD %g, top %g\n", smallestD, largestD, nodes_.front()->objectiveValue());
  }
#endif
  return nodes_.front();
}

// Add a node to the heap
void CbcTreeLocal::push(CbcNode *x)
{
  if (typeCuts_ >= 0 && !nodes_.size() && searchType_ < 0) {
    startNode_ = model_->getNodeCount();
    // save copy of node
    localNode_ = new CbcNode(*x);

    if (cut_.row().getNumElements()) {
      // Add to global cuts
      // we came in with solution
      model_->makeGlobalCut(cut_);
      if (model_->messageHandler()->logLevel() > 1)
        printf("initial cut - rhs %g %g\n",
          cut_.lb(), cut_.ub());
      searchType_ = 1;
    } else {
      // stop on first solution
      searchType_ = 0;
    }
    startTime_ = static_cast< int >(CoinCpuTime());
    saveNumberSolutions_ = model_->getSolutionCount();
  }
  nodes_.push_back(x);
#ifdef CBC_DEBUG
  if (model_->messageHandler()->logLevel() > 0)
    printf("pushing node onto heap %d %x %x\n",
      x->nodeInfo()->nodeNumber(), x, x->nodeInfo());
#endif
  std::push_heap(nodes_.begin(), nodes_.end(), comparison_);
}

// Remove the top node from the heap
void CbcTreeLocal::pop()
{
  std::pop_heap(nodes_.begin(), nodes_.end(), comparison_);
  nodes_.pop_back();
}
// Test if empty - does work if so
bool CbcTreeLocal::empty()
{
  if (typeCuts_ < 0)
    return !nodes_.size();
  /* state -
       0 iterating
       1 subtree finished optimal solution for subtree found
       2 subtree finished and no solution found
       3 subtree exiting and solution found
       4 subtree exiting and no solution found
    */
  int state = 0;
  assert(searchType_ != 2);
  if (searchType_) {
    if (CoinCpuTime() - startTime_ > timeLimit_ || model_->getNodeCount() - startNode_ >= nodeLimit_) {
      state = 4;
    }
  } else {
    if (model_->getSolutionCount() > saveNumberSolutions_) {
      state = 4;
    }
  }
  if (!nodes_.size())
    state = 2;
  if (!state) {
    return false;
  }
  // Finished this phase
  int numberColumns = model_->getNumCols();
  if (model_->getSolutionCount() > saveNumberSolutions_) {
    if (model_->getCutoff() < bestCutoff_) {
      // Save solution
      if (!bestSolution_)
        bestSolution_ = new double[numberColumns];
      memcpy(bestSolution_, model_->bestSolution(), numberColumns * sizeof(double));
      bestCutoff_ = model_->getCutoff();
    }
    state--;
  }
  // get rid of all nodes (safe even if already done)
  double bestPossibleObjective;
  cleanTree(model_, -COIN_DBL_MAX, bestPossibleObjective);

  double increment = model_->getDblParam(CbcModel::CbcCutoffIncrement);
  if (model_->messageHandler()->logLevel() > 1)
    printf("local state %d after %d nodes and %d seconds, new solution %g, best solution %g, k was %g\n",
      state,
      model_->getNodeCount() - startNode_,
      static_cast< int >(CoinCpuTime()) - startTime_,
      model_->getCutoff() + increment, bestCutoff_ + increment, rhs_);
  saveNumberSolutions_ = model_->getSolutionCount();
  bool finished = false;
  bool lastTry = false;
  switch (state) {
  case 1:
    // solution found and subtree exhausted
    if (rhs_ > 1.0e30) {
      finished = true;
    } else {
      // find global cut and reverse
      reverseCut(1);
      searchType_ = 1; // first false
      rhs_ = range_; // reset range
      nextStrong_ = false;

      // save best solution in this subtree
      memcpy(savedSolution_, model_->bestSolution(), numberColumns * sizeof(double));
    }
    break;
  case 2:
    // solution not found and subtree exhausted
    if (rhs_ > 1.0e30) {
      finished = true;
    } else {
      // find global cut and reverse
      reverseCut(2);
      searchType_ = 1; // first false
      if (diversification_ < maxDiversification_) {
        if (nextStrong_) {
          diversification_++;
          // cut is valid so don't model_->setCutoff(1.0e50);
          searchType_ = 0;
        }
        nextStrong_ = true;
        rhs_ += range_ / 2;
      } else {
        // This will be last try (may hit max time)
        lastTry = true;
        if (!maxDiversification_)
          typeCuts_ = -1; // make sure can't start again
        model_->setCutoff(bestCutoff_);
        if (model_->messageHandler()->logLevel() > 1)
          printf("Exiting local search with current set of cuts\n");
        rhs_ = 1.0e100;
        // Can now stop on gap
        model_->setDblParam(CbcModel::CbcAllowableGap, savedGap_);
      }
    }
    break;
  case 3:
    // solution found and subtree not exhausted
    if (rhs_ < 1.0e30) {
      if (searchType_) {
        if (!typeCuts_ && refine_ && searchType_ == 1) {
          // We need to check we have best solution given these 0-1 values
          OsiSolverInterface *subSolver = model_->continuousSolver()->clone();
          CbcModel *subModel = model_->subTreeModel(subSolver);
          CbcTree normalTree;
          subModel->passInTreeHandler(normalTree);
          int numberIntegers = model_->numberIntegers();
          const int *integerVariable = model_->integerVariable();
          const double *solution = model_->bestSolution();
          int i;
          int numberColumns = model_->getNumCols();
          for (i = 0; i < numberIntegers; i++) {
            int iColumn = integerVariable[i];
            double value = floor(solution[iColumn] + 0.5);
            if (!typeCuts_ && originalUpper_[i] - originalLower_[i] > 1.0)
              continue; // skip as not 0-1
            if (originalLower_[i] == originalUpper_[i])
              continue;
            subSolver->setColLower(iColumn, value);
            subSolver->setColUpper(iColumn, value);
          }
          subSolver->initialSolve();
          // We can copy cutoff
          // But adjust
          subModel->setCutoff(model_->getCutoff() + model_->getDblParam(CbcModel::CbcCutoffIncrement) + 1.0e-6);
          subModel->setSolutionCount(0);
          assert(subModel->isProvenOptimal());
          if (!subModel->typePresolve()) {
            subModel->branchAndBound();
            if (subModel->status()) {
              model_->incrementSubTreeStopped();
            }
            //printf("%g %g %g %g\n",subModel->getCutoff(),model_->getCutoff(),
            //   subModel->getMinimizationObjValue(),model_->getMinimizationObjValue());
            double newCutoff = subModel->getMinimizationObjValue() - subModel->getDblParam(CbcModel::CbcCutoffIncrement);
            if (subModel->getSolutionCount()) {
              if (!subModel->status())
                assert(subModel->isProvenOptimal());
              memcpy(model_->bestSolution(), subModel->bestSolution(),
                numberColumns * sizeof(double));
              model_->setCutoff(newCutoff);
            }
          } else if (subModel->typePresolve() == 1) {
            CbcModel *model2 = subModel->integerPresolve(true);
            if (model2) {
              // Do complete search
              model2->branchAndBound();
              // get back solution
              subModel->originalModel(model2, false);
              if (model2->status()) {
                model_->incrementSubTreeStopped();
              }
              double newCutoff = model2->getMinimizationObjValue() - model2->getDblParam(CbcModel::CbcCutoffIncrement);
              if (model2->getSolutionCount()) {
                if (!model2->status())
                  assert(model2->isProvenOptimal());
                memcpy(model_->bestSolution(), subModel->bestSolution(),
                  numberColumns * sizeof(double));
                model_->setCutoff(newCutoff);
              }
              delete model2;
            } else {
              // infeasible - could just be - due to cutoff
            }
          } else {
            // too dangerous at present
            assert(subModel->typePresolve() != 2);
          }
          if (model_->getCutoff() < bestCutoff_) {
            // Save solution
            if (!bestSolution_)
              bestSolution_ = new double[numberColumns];
            memcpy(bestSolution_, model_->bestSolution(), numberColumns * sizeof(double));
            bestCutoff_ = model_->getCutoff();
          }
          delete subModel;
        }
        // we have done search to make sure best general solution
        searchType_ = 1;
        // Reverse cut weakly
        reverseCut(3, rhs_);
      } else {
        searchType_ = 1;
        // delete last cut
        deleteCut(cut_);
      }
    } else {
      searchType_ = 1;
    }
    // save best solution in this subtree
    memcpy(savedSolution_, model_->bestSolution(), numberColumns * sizeof(double));
    nextStrong_ = false;
    rhs_ = range_;
    break;
  case 4:
    // solution not found and subtree not exhausted
    if (maxDiversification_) {
      if (nextStrong_) {
        // Reverse cut weakly
        reverseCut(4, rhs_);
        model_->setCutoff(1.0e50);
        diversification_++;
        searchType_ = 0;
      } else {
        // delete last cut
        deleteCut(cut_);
        searchType_ = 1;
      }
      nextStrong_ = true;
      rhs_ += range_ / 2;
    } else {
      // special case when using as heuristic
      // Reverse cut weakly if lb -infinity
      reverseCut(4, rhs_);
      // This will be last try (may hit max time0
      lastTry = true;
      model_->setCutoff(bestCutoff_);
      if (model_->messageHandler()->logLevel() > 1)
        printf("Exiting local search with current set of cuts\n");
      rhs_ = 1.0e100;
      // Can now stop on gap
      model_->setDblParam(CbcModel::CbcAllowableGap, savedGap_);
      typeCuts_ = -1;
    }
    break;
  }
  if (rhs_ < 1.0e30 || lastTry) {
    int goodSolution = createCut(savedSolution_, cut_);
    if (goodSolution >= 0) {
      // Add to global cuts
      model_->makeGlobalCut(cut_);
      CbcRowCuts *global = model_->globalCuts();
      int n = global->sizeRowCuts();
      OsiRowCut *rowCut = global->rowCutPtr(n - 1);
      if (model_->messageHandler()->logLevel() > 1)
        printf("inserting cut - now %d cuts, rhs %g %g, cutspace %g, diversification %d\n",
          n, rowCut->lb(), rowCut->ub(), rhs_, diversification_);
      const OsiRowCutDebugger *debugger = model_->solver()->getRowCutDebuggerAlways();
      if (debugger) {
        if (debugger->invalidCut(*rowCut))
          printf("ZZZZTree Global cut - cuts off optimal solution!\n");
      }
      for (int i = 0; i < n; i++) {
        rowCut = global->rowCutPtr(i);
        if (model_->messageHandler()->logLevel() > 1)
          printf("%d - rhs %g %g\n",
            i, rowCut->lb(), rowCut->ub());
      }
    }
    // put back node
    startTime_ = static_cast< int >(CoinCpuTime());
    startNode_ = model_->getNodeCount();
    if (localNode_) {
      // save copy of node
      CbcNode *localNode2 = new CbcNode(*localNode_);
      // But localNode2 now owns cuts so swap
      //printf("pushing local node2 onto heap %d %x %x\n",localNode_->nodeNumber(),
      //   localNode_,localNode_->nodeInfo());
      nodes_.push_back(localNode_);
      localNode_ = localNode2;
      std::make_heap(nodes_.begin(), nodes_.end(), comparison_);
    }
  }
  return finished;
}
// We may have got an intelligent tree so give it one more chance
void CbcTreeLocal::endSearch()
{
  if (typeCuts_ >= 0) {
    // copy best solution to model
    int numberColumns = model_->getNumCols();
    if (bestSolution_ && bestCutoff_ < model_->getCutoff()) {
      memcpy(model_->bestSolution(), bestSolution_, numberColumns * sizeof(double));
      model_->setCutoff(bestCutoff_);
      // recompute objective value
      const double *objCoef = model_->getObjCoefficients();
      double objOffset = 0.0;
      model_->continuousSolver()->getDblParam(OsiObjOffset, objOffset);

      // Compute dot product of objCoef and colSol and then adjust by offset
      double objValue = -objOffset;
      for (int i = 0; i < numberColumns; i++)
        objValue += objCoef[i] * bestSolution_[i];
      model_->setMinimizationObjValue(objValue);
    }
    // Can now stop on gap
    model_->setDblParam(CbcModel::CbcAllowableGap, savedGap_);
  }
}
// Create cut
int CbcTreeLocal::createCut(const double *solution, OsiRowCut &rowCut)
{
  if (rhs_ > 1.0e20)
    return -1;
  OsiSolverInterface *solver = model_->solver();
  const double *rowLower = solver->getRowLower();
  const double *rowUpper = solver->getRowUpper();
  //const double * solution = solver->getColSolution();
  //const double * objective = solver->getObjCoefficients();
  double integerTolerance = model_->getDblParam(CbcModel::CbcIntegerTolerance);
  double primalTolerance;
  solver->getDblParam(OsiPrimalTolerance, primalTolerance);
  // relax
  primalTolerance *= 1000.0;

  int numberRows = model_->getNumRows();

  int numberIntegers = model_->numberIntegers();
  const int *integerVariable = model_->integerVariable();
  int i;

  // Check feasible

  double *rowActivity = new double[numberRows];
  memset(rowActivity, 0, numberRows * sizeof(double));
  solver->getMatrixByCol()->times(solution, rowActivity);
  int goodSolution = 0;
  // check was feasible
  for (i = 0; i < numberRows; i++) {
    if (rowActivity[i] < rowLower[i] - primalTolerance) {
      goodSolution = -1;
    } else if (rowActivity[i] > rowUpper[i] + primalTolerance) {
      goodSolution = -1;
    }
  }
  delete[] rowActivity;
  for (i = 0; i < numberIntegers; i++) {
    int iColumn = integerVariable[i];
    double value = solution[iColumn];
    if (fabs(floor(value + 0.5) - value) > integerTolerance) {
      goodSolution = -1;
    }
  }
  // zap cut
  if (goodSolution == 0) {
    // Create cut and get total gap
    CoinPackedVector cut;
    double rhs = rhs_;
    double maxValue = 0.0;
    for (i = 0; i < numberIntegers; i++) {
      int iColumn = integerVariable[i];
      double value = floor(solution[iColumn] + 0.5);
      /*
              typeCuts_ == 0 restricts to binary, 1 allows general integer. But we're
              still restricted to being up against a bound. Consider: the notion is that
              the cut restricts us to a k-neighbourhood. For binary variables, this
              amounts to k variables which change value. For general integer, we could
              end up with a single variable sucking up all of k (hence mu --- the
              variable must swing to its other bound to look like a movement of 1).  For
              variables in the middle of a range, we're talking about fabs(sol<j> - x<j>).
            */
      if (!typeCuts_ && originalUpper_[i] - originalLower_[i] > 1.0)
        continue; // skip as not 0-1
      if (originalLower_[i] == originalUpper_[i])
        continue;
      double mu = 1.0 / (originalUpper_[i] - originalLower_[i]);
      if (value == originalLower_[i]) {
        rhs += mu * originalLower_[i];
        cut.insert(iColumn, 1.0);
        maxValue += originalUpper_[i];
      } else if (value == originalUpper_[i]) {
        rhs -= mu * originalUpper_[i];
        cut.insert(iColumn, -1.0);
        maxValue += originalLower_[i];
      }
    }
    if (maxValue < rhs - primalTolerance) {
      if (model_->messageHandler()->logLevel() > 1)
        printf("slack cut\n");
      goodSolution = 1;
    }
    rowCut.setRow(cut);
    rowCut.setLb(-COIN_DBL_MAX);
    rowCut.setUb(rhs);
    rowCut.setGloballyValid();
    if (model_->messageHandler()->logLevel() > 1)
      printf("Cut size: %i Cut rhs: %g\n", cut.getNumElements(), rhs);
#ifdef CBC_DEBUG
    if (model_->messageHandler()->logLevel() > 0) {
      int k;
      for (k = 0; k < cut.getNumElements(); k++) {
        printf("%i    %g ", cut.getIndices()[k], cut.getElements()[k]);
        if ((k + 1) % 5 == 0)
          printf("\n");
      }
      if (k % 5 != 0)
        printf("\n");
    }
#endif
    return goodSolution;
  } else {
    if (model_->messageHandler()->logLevel() > 1)
      printf("Not a good solution\n");
    return -1;
  }
}
// Other side of last cut branch
void CbcTreeLocal::reverseCut(int state, double bias)
{
  // find global cut
  CbcRowCuts *global = model_->globalCuts();
  int n = global->sizeRowCuts();
  int i;
  OsiRowCut *rowCut = NULL;
  for (i = 0; i < n; i++) {
    rowCut = global->rowCutPtr(i);
    if (cut_ == *rowCut) {
      break;
    }
  }
  if (!rowCut) {
    // must have got here in odd way e.g. strong branching
    return;
  }
  if (rowCut->lb() > -1.0e10)
    return;
  // get smallest element
  double smallest = COIN_DBL_MAX;
  CoinPackedVector row = cut_.row();
  for (int k = 0; k < row.getNumElements(); k++)
    smallest = CoinMin(smallest, fabs(row.getElements()[k]));
  if (!typeCuts_ && !refine_) {
    // Reverse cut very very weakly
    if (state > 2)
      smallest = 0.0;
  }
  // replace by other way
  if (model_->messageHandler()->logLevel() > 1)
    printf("reverseCut - changing cut %d out of %d, old rhs %g %g ",
      i, n, rowCut->lb(), rowCut->ub());
  rowCut->setLb(rowCut->ub() + smallest - bias);
  rowCut->setUb(COIN_DBL_MAX);
  if (model_->messageHandler()->logLevel() > 1)
    printf("new rhs %g %g, bias %g smallest %g ",
      rowCut->lb(), rowCut->ub(), bias, smallest);
  const OsiRowCutDebugger *debugger = model_->solver()->getRowCutDebuggerAlways();
  if (debugger) {
    if (debugger->invalidCut(*rowCut))
      printf("ZZZZTree Global cut - cuts off optimal solution!\n");
  }
}
// Delete last cut branch
void CbcTreeLocal::deleteCut(OsiRowCut &cut)
{
  // find global cut
  CbcRowCuts *global = model_->globalCuts();
  int n = global->sizeRowCuts();
  int i;
  OsiRowCut *rowCut = NULL;
  for (i = 0; i < n; i++) {
    rowCut = global->rowCutPtr(i);
    if (cut == *rowCut) {
      break;
    }
  }
  assert(i < n);
  // delete last cut
  if (model_->messageHandler()->logLevel() > 1)
    printf("deleteCut - deleting cut %d out of %d, rhs %g %g\n",
      i, n, rowCut->lb(), rowCut->ub());
  global->eraseRowCut(i);
}
// Create C++ lines to get to current state
void CbcTreeLocal::generateCpp(FILE *fp)
{
  CbcTreeLocal other;
  fprintf(fp, "0#include \"CbcTreeLocal.hpp\"\n");
  fprintf(fp, "5  CbcTreeLocal localTree(cbcModel,NULL);\n");
  if (range_ != other.range_)
    fprintf(fp, "5  localTree.setRange(%d);\n", range_);
  if (typeCuts_ != other.typeCuts_)
    fprintf(fp, "5  localTree.setTypeCuts(%d);\n", typeCuts_);
  if (maxDiversification_ != other.maxDiversification_)
    fprintf(fp, "5  localTree.setMaxDiversification(%d);\n", maxDiversification_);
  if (timeLimit_ != other.timeLimit_)
    fprintf(fp, "5  localTree.setTimeLimit(%d);\n", timeLimit_);
  if (nodeLimit_ != other.nodeLimit_)
    fprintf(fp, "5  localTree.setNodeLimit(%d);\n", nodeLimit_);
  if (refine_ != other.refine_)
    fprintf(fp, "5  localTree.setRefine(%s);\n", refine_ ? "true" : "false");
  fprintf(fp, "5  cbcModel->passInTreeHandler(localTree);\n");
}

CbcTreeVariable::CbcTreeVariable()
  : localNode_(NULL)
  , bestSolution_(NULL)
  , savedSolution_(NULL)
  , saveNumberSolutions_(0)
  , model_(NULL)
  , originalLower_(NULL)
  , originalUpper_(NULL)
  , range_(0)
  , typeCuts_(-1)
  , maxDiversification_(0)
  , diversification_(0)
  , nextStrong_(false)
  , rhs_(0.0)
  , savedGap_(0.0)
  , bestCutoff_(0.0)
  , timeLimit_(0)
  , startTime_(0)
  , nodeLimit_(0)
  , startNode_(-1)
  , searchType_(-1)
  , refine_(false)
{
}
/* Constructor with solution.
   range is upper bound on difference from given solution.
   maxDiversification is maximum number of diversifications to try
   timeLimit is seconds in subTree
   nodeLimit is nodes in subTree
*/
CbcTreeVariable::CbcTreeVariable(CbcModel *model, const double *solution,
  int range, int typeCuts, int maxDiversification,
  int timeLimit, int nodeLimit, bool refine)
  : localNode_(NULL)
  , bestSolution_(NULL)
  , savedSolution_(NULL)
  , saveNumberSolutions_(0)
  , model_(model)
  , originalLower_(NULL)
  , originalUpper_(NULL)
  , range_(range)
  , typeCuts_(typeCuts)
  , maxDiversification_(maxDiversification)
  , diversification_(0)
  , nextStrong_(false)
  , rhs_(0.0)
  , savedGap_(0.0)
  , bestCutoff_(0.0)
  , timeLimit_(timeLimit)
  , startTime_(0)
  , nodeLimit_(nodeLimit)
  , startNode_(-1)
  , searchType_(-1)
  , refine_(refine)
{

  OsiSolverInterface *solver = model_->solver();
  const double *lower = solver->getColLower();
  const double *upper = solver->getColUpper();
  //const double * solution = solver->getColSolution();
  //const double * objective = solver->getObjCoefficients();
  double primalTolerance;
  solver->getDblParam(OsiPrimalTolerance, primalTolerance);

  // Get increment
  model_->analyzeObjective();

  {
    // needed to sync cutoffs
    double value;
    solver->getDblParam(OsiDualObjectiveLimit, value);
    model_->setCutoff(value * solver->getObjSense());
  }
  bestCutoff_ = model_->getCutoff();
  // save current gap
  savedGap_ = model_->getDblParam(CbcModel::CbcAllowableGap);

  // make sure integers found
  model_->findIntegers(false);
  int numberIntegers = model_->numberIntegers();
  const int *integerVariable = model_->integerVariable();
  int i;
  double direction = solver->getObjSense();
  double newSolutionValue = 1.0e50;
  if (solution) {
    // copy solution
    solver->setColSolution(solution);
    newSolutionValue = direction * solver->getObjValue();
  }
  originalLower_ = new double[numberIntegers];
  originalUpper_ = new double[numberIntegers];
  bool all01 = true;
  int number01 = 0;
  for (i = 0; i < numberIntegers; i++) {
    int iColumn = integerVariable[i];
    originalLower_[i] = lower[iColumn];
    originalUpper_[i] = upper[iColumn];
    if (upper[iColumn] - lower[iColumn] > 1.5)
      all01 = false;
    else if (upper[iColumn] - lower[iColumn] == 1.0)
      number01++;
  }
  if (all01 && !typeCuts_)
    typeCuts_ = 1; // may as well so we don't have to deal with refine
  if (!number01 && !typeCuts_) {
    if (model_->messageHandler()->logLevel() > 1)
      printf("** No 0-1 variables and local search only on 0-1 - switching off\n");
    typeCuts_ = -1;
  } else {
    if (model_->messageHandler()->logLevel() > 1) {
      std::string type;
      if (all01) {
        printf("%d 0-1 variables normal local  cuts\n",
          number01);
      } else if (typeCuts_) {
        printf("%d 0-1 variables, %d other - general integer local cuts\n",
          number01, numberIntegers - number01);
      } else {
        printf("%d 0-1 variables, %d other - local cuts but just on 0-1 variables\n",
          number01, numberIntegers - number01);
      }
      printf("maximum diversifications %d, initial cutspace %d, max time %d seconds, max nodes %d\n",
        maxDiversification_, range_, timeLimit_, nodeLimit_);
    }
  }
  int numberColumns = model_->getNumCols();
  savedSolution_ = new double[numberColumns];
  memset(savedSolution_, 0, numberColumns * sizeof(double));
  if (solution) {
    rhs_ = range_;
    // Check feasible
    int goodSolution = createCut(solution, cut_);
    if (goodSolution >= 0) {
      for (i = 0; i < numberIntegers; i++) {
        int iColumn = integerVariable[i];
        double value = floor(solution[iColumn] + 0.5);
        // fix so setBestSolution will work
        solver->setColLower(iColumn, value);
        solver->setColUpper(iColumn, value);
      }
      model_->reserveCurrentSolution();
      // Create cut and get total gap
      if (newSolutionValue < bestCutoff_) {
        model_->setBestSolution(CBC_ROUNDING, newSolutionValue, solution);
        bestCutoff_ = model_->getCutoff();
        // save as best solution
        memcpy(savedSolution_, model_->bestSolution(), numberColumns * sizeof(double));
      }
      for (i = 0; i < numberIntegers; i++) {
        int iColumn = integerVariable[i];
        // restore bounds
        solver->setColLower(iColumn, originalLower_[i]);
        solver->setColUpper(iColumn, originalUpper_[i]);
      }
      // make sure can't stop on gap
      model_->setDblParam(CbcModel::CbcAllowableGap, -1.0e50);
    } else {
      model_ = NULL;
    }
  } else {
    // no solution
    rhs_ = 1.0e50;
    // make sure can't stop on gap
    model_->setDblParam(CbcModel::CbcAllowableGap, -1.0e50);
  }
}
CbcTreeVariable::~CbcTreeVariable()
{
  delete[] originalLower_;
  delete[] originalUpper_;
  delete[] bestSolution_;
  delete[] savedSolution_;
  delete localNode_;
}
// Copy constructor
CbcTreeVariable::CbcTreeVariable(const CbcTreeVariable &rhs)
  : CbcTree(rhs)
  , saveNumberSolutions_(rhs.saveNumberSolutions_)
  , model_(rhs.model_)
  , range_(rhs.range_)
  , typeCuts_(rhs.typeCuts_)
  , maxDiversification_(rhs.maxDiversification_)
  , diversification_(rhs.diversification_)
  , nextStrong_(rhs.nextStrong_)
  , rhs_(rhs.rhs_)
  , savedGap_(rhs.savedGap_)
  , bestCutoff_(rhs.bestCutoff_)
  , timeLimit_(rhs.timeLimit_)
  , startTime_(rhs.startTime_)
  , nodeLimit_(rhs.nodeLimit_)
  , startNode_(rhs.startNode_)
  , searchType_(rhs.searchType_)
  , refine_(rhs.refine_)
{
  cut_ = rhs.cut_;
  fixedCut_ = rhs.fixedCut_;
  if (rhs.localNode_)
    localNode_ = new CbcNode(*rhs.localNode_);
  else
    localNode_ = NULL;
  if (rhs.originalLower_) {
    int numberIntegers = model_->numberIntegers();
    originalLower_ = new double[numberIntegers];
    memcpy(originalLower_, rhs.originalLower_, numberIntegers * sizeof(double));
    originalUpper_ = new double[numberIntegers];
    memcpy(originalUpper_, rhs.originalUpper_, numberIntegers * sizeof(double));
  } else {
    originalLower_ = NULL;
    originalUpper_ = NULL;
  }
  if (rhs.bestSolution_) {
    int numberColumns = model_->getNumCols();
    bestSolution_ = new double[numberColumns];
    memcpy(bestSolution_, rhs.bestSolution_, numberColumns * sizeof(double));
  } else {
    bestSolution_ = NULL;
  }
  if (rhs.savedSolution_) {
    int numberColumns = model_->getNumCols();
    savedSolution_ = new double[numberColumns];
    memcpy(savedSolution_, rhs.savedSolution_, numberColumns * sizeof(double));
  } else {
    savedSolution_ = NULL;
  }
}
//----------------------------------------------------------------
// Assignment operator
//-------------------------------------------------------------------
CbcTreeVariable &
CbcTreeVariable::operator=(const CbcTreeVariable &rhs)
{
  if (this != &rhs) {
    CbcTree::operator=(rhs);
    saveNumberSolutions_ = rhs.saveNumberSolutions_;
    cut_ = rhs.cut_;
    fixedCut_ = rhs.fixedCut_;
    delete localNode_;
    if (rhs.localNode_)
      localNode_ = new CbcNode(*rhs.localNode_);
    else
      localNode_ = NULL;
    model_ = rhs.model_;
    range_ = rhs.range_;
    typeCuts_ = rhs.typeCuts_;
    maxDiversification_ = rhs.maxDiversification_;
    diversification_ = rhs.diversification_;
    nextStrong_ = rhs.nextStrong_;
    rhs_ = rhs.rhs_;
    savedGap_ = rhs.savedGap_;
    bestCutoff_ = rhs.bestCutoff_;
    timeLimit_ = rhs.timeLimit_;
    startTime_ = rhs.startTime_;
    nodeLimit_ = rhs.nodeLimit_;
    startNode_ = rhs.startNode_;
    searchType_ = rhs.searchType_;
    refine_ = rhs.refine_;
    delete[] originalLower_;
    delete[] originalUpper_;
    if (rhs.originalLower_) {
      int numberIntegers = model_->numberIntegers();
      originalLower_ = new double[numberIntegers];
      memcpy(originalLower_, rhs.originalLower_, numberIntegers * sizeof(double));
      originalUpper_ = new double[numberIntegers];
      memcpy(originalUpper_, rhs.originalUpper_, numberIntegers * sizeof(double));
    } else {
      originalLower_ = NULL;
      originalUpper_ = NULL;
    }
    delete[] bestSolution_;
    if (rhs.bestSolution_) {
      int numberColumns = model_->getNumCols();
      bestSolution_ = new double[numberColumns];
      memcpy(bestSolution_, rhs.bestSolution_, numberColumns * sizeof(double));
    } else {
      bestSolution_ = NULL;
    }
    delete[] savedSolution_;
    if (rhs.savedSolution_) {
      int numberColumns = model_->getNumCols();
      savedSolution_ = new double[numberColumns];
      memcpy(savedSolution_, rhs.savedSolution_, numberColumns * sizeof(double));
    } else {
      savedSolution_ = NULL;
    }
  }
  return *this;
}
// Clone
CbcTree *
CbcTreeVariable::clone() const
{
  return new CbcTreeVariable(*this);
}
// Pass in solution (so can be used after heuristic)
void CbcTreeVariable::passInSolution(const double *solution, double solutionValue)
{
  int numberColumns = model_->getNumCols();
  delete[] savedSolution_;
  savedSolution_ = new double[numberColumns];
  memcpy(savedSolution_, solution, numberColumns * sizeof(double));
  rhs_ = range_;
  // Check feasible
  int goodSolution = createCut(solution, cut_);
  if (goodSolution >= 0) {
    bestCutoff_ = CoinMin(solutionValue, model_->getCutoff());
  } else {
    model_ = NULL;
  }
}
// Return the top node of the heap
CbcNode *
CbcTreeVariable::top() const
{
#ifdef CBC_DEBUG
  int smallest = 9999999;
  int largest = -1;
  double smallestD = 1.0e30;
  double largestD = -1.0e30;
  int n = nodes_.size();
  for (int i = 0; i < n; i++) {
    int nn = nodes_[i]->nodeInfo()->nodeNumber();
    double dd = nodes_[i]->objectiveValue();
    largest = CoinMax(largest, nn);
    smallest = CoinMin(smallest, nn);
    largestD = CoinMax(largestD, dd);
    smallestD = CoinMin(smallestD, dd);
  }
  if (model_->messageHandler()->logLevel() > 1) {
    printf("smallest %d, largest %d, top %d\n", smallest, largest,
      nodes_.front()->nodeInfo()->nodeNumber());
    printf("smallestD %g, largestD %g, top %g\n", smallestD, largestD, nodes_.front()->objectiveValue());
  }
#endif
  return nodes_.front();
}

// Add a node to the heap
void CbcTreeVariable::push(CbcNode *x)
{
  if (typeCuts_ >= 0 && !nodes_.size() && searchType_ < 0) {
    startNode_ = model_->getNodeCount();
    // save copy of node
    localNode_ = new CbcNode(*x);

    if (cut_.row().getNumElements()) {
      // Add to global cuts
      // we came in with solution
      model_->makeGlobalCut(cut_);
      if (model_->messageHandler()->logLevel() > 1)
        printf("initial cut - rhs %g %g\n",
          cut_.lb(), cut_.ub());
      searchType_ = 1;
    } else {
      // stop on first solution
      searchType_ = 0;
    }
    startTime_ = static_cast< int >(CoinCpuTime());
    saveNumberSolutions_ = model_->getSolutionCount();
  }
  nodes_.push_back(x);
#ifdef CBC_DEBUG
  if (model_->messageHandler()->logLevel() > 0)
    printf("pushing node onto heap %d %x %x\n",
      x->nodeInfo()->nodeNumber(), x, x->nodeInfo());
#endif
  std::push_heap(nodes_.begin(), nodes_.end(), comparison_);
}

// Remove the top node from the heap
void CbcTreeVariable::pop()
{
  std::pop_heap(nodes_.begin(), nodes_.end(), comparison_);
  nodes_.pop_back();
}
// Test if empty - does work if so
bool CbcTreeVariable::empty()
{
  if (typeCuts_ < 0)
    return !nodes_.size();
  /* state -
       0 iterating
       1 subtree finished optimal solution for subtree found
       2 subtree finished and no solution found
       3 subtree exiting and solution found
       4 subtree exiting and no solution found
    */
  int state = 0;
  assert(searchType_ != 2);
  if (searchType_) {
    if (CoinCpuTime() - startTime_ > timeLimit_ || model_->getNodeCount() - startNode_ >= nodeLimit_) {
      state = 4;
    }
  } else {
    if (model_->getSolutionCount() > saveNumberSolutions_) {
      state = 4;
    }
  }
  if (!nodes_.size())
    state = 2;
  if (!state) {
    return false;
  }
  // Finished this phase
  int numberColumns = model_->getNumCols();
  if (model_->getSolutionCount() > saveNumberSolutions_) {
    if (model_->getCutoff() < bestCutoff_) {
      // Save solution
      if (!bestSolution_)
        bestSolution_ = new double[numberColumns];
      memcpy(bestSolution_, model_->bestSolution(), numberColumns * sizeof(double));
      bestCutoff_ = model_->getCutoff();
    }
    state--;
  }
  // get rid of all nodes (safe even if already done)
  double bestPossibleObjective;
  cleanTree(model_, -COIN_DBL_MAX, bestPossibleObjective);

  double increment = model_->getDblParam(CbcModel::CbcCutoffIncrement);
  if (model_->messageHandler()->logLevel() > 1)
    printf("local state %d after %d nodes and %d seconds, new solution %g, best solution %g, k was %g\n",
      state,
      model_->getNodeCount() - startNode_,
      static_cast< int >(CoinCpuTime()) - startTime_,
      model_->getCutoff() + increment, bestCutoff_ + increment, rhs_);
  saveNumberSolutions_ = model_->getSolutionCount();
  bool finished = false;
  bool lastTry = false;
  switch (state) {
  case 1:
    // solution found and subtree exhausted
    if (rhs_ > 1.0e30) {
      finished = true;
    } else {
      // find global cut and reverse
      reverseCut(1);
      searchType_ = 1; // first false
      rhs_ = range_; // reset range
      nextStrong_ = false;

      // save best solution in this subtree
      memcpy(savedSolution_, model_->bestSolution(), numberColumns * sizeof(double));
    }
    break;
  case 2:
    // solution not found and subtree exhausted
    if (rhs_ > 1.0e30) {
      finished = true;
    } else {
      // find global cut and reverse
      reverseCut(2);
      searchType_ = 1; // first false
      if (diversification_ < maxDiversification_) {
        if (nextStrong_) {
          diversification_++;
          // cut is valid so don't model_->setCutoff(1.0e50);
          searchType_ = 0;
        }
        nextStrong_ = true;
        rhs_ += range_ / 2;
      } else {
        // This will be last try (may hit max time)
        lastTry = true;
        if (!maxDiversification_)
          typeCuts_ = -1; // make sure can't start again
        model_->setCutoff(bestCutoff_);
        if (model_->messageHandler()->logLevel() > 1)
          printf("Exiting local search with current set of cuts\n");
        rhs_ = 1.0e100;
        // Can now stop on gap
        model_->setDblParam(CbcModel::CbcAllowableGap, savedGap_);
      }
    }
    break;
  case 3:
    // solution found and subtree not exhausted
    if (rhs_ < 1.0e30) {
      if (searchType_) {
        if (!typeCuts_ && refine_ && searchType_ == 1) {
          // We need to check we have best solution given these 0-1 values
          OsiSolverInterface *subSolver = model_->continuousSolver()->clone();
          CbcModel *subModel = model_->subTreeModel(subSolver);
          CbcTree normalTree;
          subModel->passInTreeHandler(normalTree);
          int numberIntegers = model_->numberIntegers();
          const int *integerVariable = model_->integerVariable();
          const double *solution = model_->bestSolution();
          int i;
          int numberColumns = model_->getNumCols();
          for (i = 0; i < numberIntegers; i++) {
            int iColumn = integerVariable[i];
            double value = floor(solution[iColumn] + 0.5);
            if (!typeCuts_ && originalUpper_[i] - originalLower_[i] > 1.0)
              continue; // skip as not 0-1
            if (originalLower_[i] == originalUpper_[i])
              continue;
            subSolver->setColLower(iColumn, value);
            subSolver->setColUpper(iColumn, value);
          }
          subSolver->initialSolve();
          // We can copy cutoff
          // But adjust
          subModel->setCutoff(model_->getCutoff() + model_->getDblParam(CbcModel::CbcCutoffIncrement) + 1.0e-6);
          subModel->setSolutionCount(0);
          assert(subModel->isProvenOptimal());
          if (!subModel->typePresolve()) {
            subModel->branchAndBound();
            if (subModel->status()) {
              model_->incrementSubTreeStopped();
            }
            //printf("%g %g %g %g\n",subModel->getCutoff(),model_->getCutoff(),
            //   subModel->getMinimizationObjValue(),model_->getMinimizationObjValue());
            double newCutoff = subModel->getMinimizationObjValue() - subModel->getDblParam(CbcModel::CbcCutoffIncrement);
            if (subModel->getSolutionCount()) {
              if (!subModel->status())
                assert(subModel->isProvenOptimal());
              memcpy(model_->bestSolution(), subModel->bestSolution(),
                numberColumns * sizeof(double));
              model_->setCutoff(newCutoff);
            }
          } else if (subModel->typePresolve() == 1) {
            CbcModel *model2 = subModel->integerPresolve(true);
            if (model2) {
              // Do complete search
              model2->branchAndBound();
              // get back solution
              subModel->originalModel(model2, false);
              if (model2->status()) {
                model_->incrementSubTreeStopped();
              }
              double newCutoff = model2->getMinimizationObjValue() - model2->getDblParam(CbcModel::CbcCutoffIncrement);
              if (model2->getSolutionCount()) {
                if (!model2->status())
                  assert(model2->isProvenOptimal());
                memcpy(model_->bestSolution(), subModel->bestSolution(),
                  numberColumns * sizeof(double));
                model_->setCutoff(newCutoff);
              }
              delete model2;
            } else {
              // infeasible - could just be - due to cutoff
            }
          } else {
            // too dangerous at present
            assert(subModel->typePresolve() != 2);
          }
          if (model_->getCutoff() < bestCutoff_) {
            // Save solution
            if (!bestSolution_)
              bestSolution_ = new double[numberColumns];
            memcpy(bestSolution_, model_->bestSolution(), numberColumns * sizeof(double));
            bestCutoff_ = model_->getCutoff();
          }
          delete subModel;
        }
        // we have done search to make sure best general solution
        searchType_ = 1;
        // Reverse cut weakly
        reverseCut(3, rhs_);
      } else {
        searchType_ = 1;
        // delete last cut
        deleteCut(cut_);
      }
    } else {
      searchType_ = 1;
    }
    // save best solution in this subtree
    memcpy(savedSolution_, model_->bestSolution(), numberColumns * sizeof(double));
    nextStrong_ = false;
    rhs_ = range_;
    break;
  case 4:
    // solution not found and subtree not exhausted
    if (maxDiversification_) {
      if (nextStrong_) {
        // Reverse cut weakly
        reverseCut(4, rhs_);
        model_->setCutoff(1.0e50);
        diversification_++;
        searchType_ = 0;
      } else {
        // delete last cut
        deleteCut(cut_);
        searchType_ = 1;
      }
      nextStrong_ = true;
      rhs_ += range_ / 2;
    } else {
      // special case when using as heuristic
      // Reverse cut weakly if lb -infinity
      reverseCut(4, rhs_);
      // This will be last try (may hit max time0
      lastTry = true;
      model_->setCutoff(bestCutoff_);
      if (model_->messageHandler()->logLevel() > 1)
        printf("Exiting local search with current set of cuts\n");
      rhs_ = 1.0e100;
      // Can now stop on gap
      model_->setDblParam(CbcModel::CbcAllowableGap, savedGap_);
      typeCuts_ = -1;
    }
    break;
  }
  if (rhs_ < 1.0e30 || lastTry) {
    int goodSolution = createCut(savedSolution_, cut_);
    if (goodSolution >= 0) {
      // Add to global cuts
      model_->makeGlobalCut(cut_);
      CbcRowCuts *global = model_->globalCuts();
      int n = global->sizeRowCuts();
      OsiRowCut *rowCut = global->rowCutPtr(n - 1);
      if (model_->messageHandler()->logLevel() > 1)
        printf("inserting cut - now %d cuts, rhs %g %g, cutspace %g, diversification %d\n",
          n, rowCut->lb(), rowCut->ub(), rhs_, diversification_);
      const OsiRowCutDebugger *debugger = model_->solver()->getRowCutDebuggerAlways();
      if (debugger) {
        if (debugger->invalidCut(*rowCut))
          printf("ZZZZTree Global cut - cuts off optimal solution!\n");
      }
      for (int i = 0; i < n; i++) {
        rowCut = global->rowCutPtr(i);
        if (model_->messageHandler()->logLevel() > 1)
          printf("%d - rhs %g %g\n",
            i, rowCut->lb(), rowCut->ub());
      }
    }
    // put back node
    startTime_ = static_cast< int >(CoinCpuTime());
    startNode_ = model_->getNodeCount();
    if (localNode_) {
      // save copy of node
      CbcNode *localNode2 = new CbcNode(*localNode_);
      // But localNode2 now owns cuts so swap
      //printf("pushing local node2 onto heap %d %x %x\n",localNode_->nodeNumber(),
      //   localNode_,localNode_->nodeInfo());
      nodes_.push_back(localNode_);
      localNode_ = localNode2;
      std::make_heap(nodes_.begin(), nodes_.end(), comparison_);
    }
  }
  return finished;
}
// We may have got an intelligent tree so give it one more chance
void CbcTreeVariable::endSearch()
{
  if (typeCuts_ >= 0) {
    // copy best solution to model
    int numberColumns = model_->getNumCols();
    if (bestSolution_ && bestCutoff_ < model_->getCutoff()) {
      memcpy(model_->bestSolution(), bestSolution_, numberColumns * sizeof(double));
      model_->setCutoff(bestCutoff_);
      // recompute objective value
      const double *objCoef = model_->getObjCoefficients();
      double objOffset = 0.0;
      model_->continuousSolver()->getDblParam(OsiObjOffset, objOffset);

      // Compute dot product of objCoef and colSol and then adjust by offset
      double objValue = -objOffset;
      for (int i = 0; i < numberColumns; i++)
        objValue += objCoef[i] * bestSolution_[i];
      model_->setMinimizationObjValue(objValue);
    }
    // Can now stop on gap
    model_->setDblParam(CbcModel::CbcAllowableGap, savedGap_);
  }
}
// Create cut
int CbcTreeVariable::createCut(const double *solution, OsiRowCut &rowCut)
{
  if (rhs_ > 1.0e20)
    return -1;
  OsiSolverInterface *solver = model_->solver();
  const double *rowLower = solver->getRowLower();
  const double *rowUpper = solver->getRowUpper();
  //const double * solution = solver->getColSolution();
  //const double * objective = solver->getObjCoefficients();
  double integerTolerance = model_->getDblParam(CbcModel::CbcIntegerTolerance);
  double primalTolerance;
  solver->getDblParam(OsiPrimalTolerance, primalTolerance);
  // relax
  primalTolerance *= 1000.0;

  int numberRows = model_->getNumRows();

  int numberIntegers = model_->numberIntegers();
  const int *integerVariable = model_->integerVariable();
  int i;

  // Check feasible

  double *rowActivity = new double[numberRows];
  memset(rowActivity, 0, numberRows * sizeof(double));
  solver->getMatrixByCol()->times(solution, rowActivity);
  int goodSolution = 0;
  // check was feasible
  for (i = 0; i < numberRows; i++) {
    if (rowActivity[i] < rowLower[i] - primalTolerance) {
      goodSolution = -1;
    } else if (rowActivity[i] > rowUpper[i] + primalTolerance) {
      goodSolution = -1;
    }
  }
  delete[] rowActivity;
  for (i = 0; i < numberIntegers; i++) {
    int iColumn = integerVariable[i];
    double value = solution[iColumn];
    if (fabs(floor(value + 0.5) - value) > integerTolerance) {
      goodSolution = -1;
    }
  }
  // zap cut
  if (goodSolution == 0) {
    // Create cut and get total gap
    CoinPackedVector cut;
    double rhs = rhs_;
    double maxValue = 0.0;
    for (i = 0; i < numberIntegers; i++) {
      int iColumn = integerVariable[i];
      double value = floor(solution[iColumn] + 0.5);
      if (!typeCuts_ && originalUpper_[i] - originalLower_[i] > 1.0)
        continue; // skip as not 0-1
      if (originalLower_[i] == originalUpper_[i])
        continue;
      double mu = 1.0 / (originalUpper_[i] - originalLower_[i]);
      if (value == originalLower_[i]) {
        rhs += mu * originalLower_[i];
        cut.insert(iColumn, 1.0);
        maxValue += originalUpper_[i];
      } else if (value == originalUpper_[i]) {
        rhs -= mu * originalUpper_[i];
        cut.insert(iColumn, -1.0);
        maxValue += originalLower_[i];
      }
    }
    if (maxValue < rhs - primalTolerance) {
      if (model_->messageHandler()->logLevel() > 1)
        printf("slack cut\n");
      goodSolution = 1;
    }
    rowCut.setRow(cut);
    rowCut.setLb(-COIN_DBL_MAX);
    rowCut.setUb(rhs);
    rowCut.setGloballyValid();
    if (model_->messageHandler()->logLevel() > 1)
      printf("Cut size: %i Cut rhs: %g\n", cut.getNumElements(), rhs);
#ifdef CBC_DEBUG
    if (model_->messageHandler()->logLevel() > 0) {
      int k;
      for (k = 0; k < cut.getNumElements(); k++) {
        printf("%i    %g ", cut.getIndices()[k], cut.getElements()[k]);
        if ((k + 1) % 5 == 0)
          printf("\n");
      }
      if (k % 5 != 0)
        printf("\n");
    }
#endif
    return goodSolution;
  } else {
    if (model_->messageHandler()->logLevel() > 1)
      printf("Not a good solution\n");
    return -1;
  }
}
// Other side of last cut branch
void CbcTreeVariable::reverseCut(int state, double bias)
{
  // find global cut
  CbcRowCuts *global = model_->globalCuts();
  int n = global->sizeRowCuts();
  int i;
  OsiRowCut *rowCut = NULL;
  for (i = 0; i < n; i++) {
    rowCut = global->rowCutPtr(i);
    if (cut_ == *rowCut) {
      break;
    }
  }
  if (!rowCut) {
    // must have got here in odd way e.g. strong branching
    return;
  }
  if (rowCut->lb() > -1.0e10)
    return;
  // get smallest element
  double smallest = COIN_DBL_MAX;
  CoinPackedVector row = cut_.row();
  for (int k = 0; k < row.getNumElements(); k++)
    smallest = CoinMin(smallest, fabs(row.getElements()[k]));
  if (!typeCuts_ && !refine_) {
    // Reverse cut very very weakly
    if (state > 2)
      smallest = 0.0;
  }
  // replace by other way
  if (model_->messageHandler()->logLevel() > 1)
    printf("reverseCut - changing cut %d out of %d, old rhs %g %g ",
      i, n, rowCut->lb(), rowCut->ub());
  rowCut->setLb(rowCut->ub() + smallest - bias);
  rowCut->setUb(COIN_DBL_MAX);
  if (model_->messageHandler()->logLevel() > 1)
    printf("new rhs %g %g, bias %g smallest %g ",
      rowCut->lb(), rowCut->ub(), bias, smallest);
  const OsiRowCutDebugger *debugger = model_->solver()->getRowCutDebuggerAlways();
  if (debugger) {
    if (debugger->invalidCut(*rowCut))
      printf("ZZZZTree Global cut - cuts off optimal solution!\n");
  }
}
// Delete last cut branch
void CbcTreeVariable::deleteCut(OsiRowCut &cut)
{
  // find global cut
  CbcRowCuts *global = model_->globalCuts();
  int n = global->sizeRowCuts();
  int i;
  OsiRowCut *rowCut = NULL;
  for (i = 0; i < n; i++) {
    rowCut = global->rowCutPtr(i);
    if (cut == *rowCut) {
      break;
    }
  }
  assert(i < n);
  // delete last cut
  if (model_->messageHandler()->logLevel() > 1)
    printf("deleteCut - deleting cut %d out of %d, rhs %g %g\n",
      i, n, rowCut->lb(), rowCut->ub());
  global->eraseRowCut(i);
}
// Create C++ lines to get to current state
void CbcTreeVariable::generateCpp(FILE *fp)
{
  CbcTreeVariable other;
  fprintf(fp, "0#include \"CbcTreeVariable.hpp\"\n");
  fprintf(fp, "5  CbcTreeVariable variableTree(cbcModel,NULL);\n");
  if (range_ != other.range_)
    fprintf(fp, "5  variableTree.setRange(%d);\n", range_);
  if (typeCuts_ != other.typeCuts_)
    fprintf(fp, "5  variableTree.setTypeCuts(%d);\n", typeCuts_);
  if (maxDiversification_ != other.maxDiversification_)
    fprintf(fp, "5  variableTree.setMaxDiversification(%d);\n", maxDiversification_);
  if (timeLimit_ != other.timeLimit_)
    fprintf(fp, "5  variableTree.setTimeLimit(%d);\n", timeLimit_);
  if (nodeLimit_ != other.nodeLimit_)
    fprintf(fp, "5  variableTree.setNodeLimit(%d);\n", nodeLimit_);
  if (refine_ != other.refine_)
    fprintf(fp, "5  variableTree.setRefine(%s);\n", refine_ ? "true" : "false");
  fprintf(fp, "5  cbcModel->passInTreeHandler(variableTree);\n");
}

/* vi: softtabstop=2 shiftwidth=2 expandtab tabstop=2
*/
