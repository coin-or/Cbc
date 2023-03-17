// $Id: CbcHeuristicDW.cpp 1899 2013-04-09 18:12:08Z stefan $
// Copyright (C) 2006, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

#if defined(_MSC_VER)
// Turn off compiler warning about long names
#pragma warning(disable : 4786)
#endif
#include <cassert>
#include <cstdlib>
#include <cmath>
#include <cfloat>

#include "OsiSolverInterface.hpp"
#include "OsiClpSolverInterface.hpp"
#include "CbcModel.hpp"
#include "CbcMessage.hpp"
#include "CbcHeuristicDW.hpp"
#include "CbcStrategy.hpp"
#include "ClpPresolve.hpp"
#include "CglProbing.hpp"

static int dummyCallBack(CbcHeuristicDW * /*heuristic*/,
  CbcModel * /*thisModel*/, int /*whereFrom*/)
{
  return 0;
}

// Default Constructor
CbcHeuristicDW::CbcHeuristicDW()
  : CbcHeuristic()
{
  setDefaults();
}

// Constructor with model - assumed before cuts

CbcHeuristicDW::CbcHeuristicDW(CbcModel &model, int keepContinuous)
  : CbcHeuristic(model)
{
  setDefaults();
  functionPointer_ = dummyCallBack;
  assert(model.solver());
  solver_ = model.solver()->clone();
  findStructure();
}
/* Constructor with model - assumed before cuts
 */
CbcHeuristicDW::CbcHeuristicDW(CbcModel &model,
  int callBack(CbcHeuristicDW *currentHeuristic,
    CbcModel *thisModel,
    int whereFrom),
  int keepContinuous)
  : CbcHeuristic(model)
{
  setDefaults();
  functionPointer_ = callBack;
  assert(model.solver());
  solver_ = model.solver()->clone();
  findStructure();
}
// Set default values
void CbcHeuristicDW::setDefaults()
{
  targetObjective_ = -COIN_DBL_MAX;
  bestObjective_ = COIN_DBL_MAX;
  lastObjective_ = COIN_DBL_MAX;
  fullDWEverySoOften_ = 0;
  numberPasses_ = 0;
  howOften_ = 100;
  decayFactor_ = 0.5;
  functionPointer_ = NULL;
  solver_ = NULL;
  dwSolver_ = NULL;
  bestSolution_ = NULL;
  continuousSolution_ = NULL;
  fixedDj_ = NULL;
  saveLower_ = NULL;
  saveUpper_ = NULL;
  random_ = NULL;
  affinity_ = NULL;
  weights_ = NULL;
  objectiveDW_ = NULL;
  numberColumnsDW_ = NULL;
  whichRowBlock_ = NULL;
  whichColumnBlock_ = NULL;
  dwBlock_ = NULL;
  backwardRow_ = NULL;
  rowsInBlock_ = NULL;
  columnsInBlock_ = NULL;
  startRowBlock_ = NULL;
  startColumnBlock_ = NULL;
  intsInBlock_ = NULL;
  fingerPrint_ = NULL;
  fullDWEverySoOften_ = 0;
  numberPasses_ = 0;
  numberBadPasses_ = COIN_INT_MAX;
  maximumDW_ = 0;
  numberDW_ = 0;
  numberDWTimes_ = 0;
  sizeFingerPrint_ = 0;
  numberMasterColumns_ = 0;
  numberMasterRows_ = 0;
  numberBlocks_ = 0;
  keepContinuous_ = 0;
  phase_ = 0;
  pass_ = 0;
  nNeededBase_ = 200;
  nNodesBase_ = 500;
  nNeeded_ = nNeededBase_;
  nNodes_ = nNodesBase_;
  solveState_ = 0;
}
// Guts of copy
void CbcHeuristicDW::gutsOfCopy(const CbcHeuristicDW &rhs)
{
  targetObjective_ = rhs.targetObjective_;
  bestObjective_ = rhs.bestObjective_;
  lastObjective_ = rhs.lastObjective_;
  fullDWEverySoOften_ = rhs.fullDWEverySoOften_;
  numberPasses_ = rhs.numberPasses_;
  numberBadPasses_ = rhs.numberBadPasses_;
  howOften_ = rhs.howOften_;
  decayFactor_ = rhs.decayFactor_;
  fullDWEverySoOften_ = rhs.fullDWEverySoOften_;
  numberPasses_ = rhs.numberPasses_;
  maximumDW_ = rhs.maximumDW_;
  numberDW_ = rhs.numberDW_;
  numberDWTimes_ = rhs.numberDWTimes_;
  sizeFingerPrint_ = rhs.sizeFingerPrint_;
  numberMasterColumns_ = rhs.numberMasterColumns_;
  numberMasterRows_ = rhs.numberMasterRows_;
  numberBlocks_ = rhs.numberBlocks_;
  keepContinuous_ = rhs.keepContinuous_;
  phase_ = rhs.phase_;
  pass_ = rhs.pass_;
  nNeededBase_ = rhs.nNeededBase_;
  nNodesBase_ = rhs.nNodesBase_;
  nNeeded_ = rhs.nNeeded_;
  nNodes_ = rhs.nNodes_;
  solveState_ = rhs.solveState_;
  functionPointer_ = rhs.functionPointer_;
  if (rhs.solver_)
    solver_ = rhs.solver_->clone();
  else
    solver_ = NULL;
  if (rhs.dwSolver_)
    dwSolver_ = rhs.dwSolver_->clone();
  else
    dwSolver_ = NULL;
  if (rhs.saveLower_) {
    int numberColumns = solver_->getNumCols();
    int numberRows = solver_->getNumRows();
    saveLower_ = CoinCopyOfArray(rhs.saveLower_, numberColumns);
    saveUpper_ = CoinCopyOfArray(rhs.saveUpper_, numberColumns);
    whichColumnBlock_ = CoinCopyOfArray(rhs.whichColumnBlock_, numberColumns);
    columnsInBlock_ = CoinCopyOfArray(rhs.columnsInBlock_, numberColumns);
    whichRowBlock_ = CoinCopyOfArray(rhs.whichRowBlock_, numberRows);
    rowsInBlock_ = CoinCopyOfArray(rhs.rowsInBlock_, numberRows);
    if (rhs.affinity_)
      affinity_ = CoinCopyOfArray(rhs.affinity_, numberBlocks_ * numberBlocks_);
    else
      affinity_ = NULL;
    backwardRow_ = CoinCopyOfArray(rhs.backwardRow_, numberRows);
    startRowBlock_ = CoinCopyOfArray(rhs.startRowBlock_, numberBlocks_ + 1);
    startColumnBlock_ = CoinCopyOfArray(rhs.startColumnBlock_, numberBlocks_ + 1);
    intsInBlock_ = CoinCopyOfArray(rhs.intsInBlock_, numberBlocks_);
  } else {
    saveLower_ = NULL;
    saveUpper_ = NULL;
    affinity_ = NULL;
    whichRowBlock_ = NULL;
    whichColumnBlock_ = NULL;
    backwardRow_ = NULL;
    rowsInBlock_ = NULL;
    columnsInBlock_ = NULL;
    startRowBlock_ = NULL;
    startColumnBlock_ = NULL;
    intsInBlock_ = NULL;
  }
  if (rhs.weights_) {
    assert(maximumDW_);
    weights_ = CoinCopyOfArray(rhs.weights_, maximumDW_);
    random_ = CoinCopyOfArray(rhs.random_, numberMasterRows_);
    dwBlock_ = CoinCopyOfArray(rhs.dwBlock_, maximumDW_);
    fingerPrint_ = CoinCopyOfArray(rhs.fingerPrint_,
      sizeFingerPrint_ * maximumDW_);
    objectiveDW_ = CoinCopyOfArray(rhs.objectiveDW_, numberDWTimes_);
    numberColumnsDW_ = CoinCopyOfArray(rhs.numberColumnsDW_, numberDWTimes_);
  } else {
    random_ = NULL;
    weights_ = NULL;
    objectiveDW_ = NULL;
    numberColumnsDW_ = NULL;
    dwBlock_ = NULL;
    fingerPrint_ = NULL;
  }
  if (rhs.bestSolution_) {
    int numberColumns = solver_->getNumCols();
    bestSolution_ = CoinCopyOfArray(rhs.bestSolution_, numberColumns);
  } else {
    bestSolution_ = NULL;
  }
  if (rhs.continuousSolution_) {
    int numberColumns = solver_->getNumCols();
    continuousSolution_ = CoinCopyOfArray(rhs.continuousSolution_, numberColumns);
  } else {
    continuousSolution_ = NULL;
  }
  if (rhs.fixedDj_) {
    int numberColumns = solver_->getNumCols();
    fixedDj_ = CoinCopyOfArray(rhs.fixedDj_, numberColumns);
  } else {
    fixedDj_ = NULL;
  }
}
// Guts of delete
void CbcHeuristicDW::gutsOfDelete()
{
  delete solver_;
  delete dwSolver_;
  delete[] bestSolution_;
  delete[] continuousSolution_;
  delete[] fixedDj_;
  delete[] saveLower_;
  delete[] saveUpper_;
  delete[] random_;
  delete[] affinity_;
  delete[] weights_;
  delete[] objectiveDW_;
  delete[] numberColumnsDW_;
  delete[] whichRowBlock_;
  delete[] whichColumnBlock_;
  delete[] dwBlock_;
  delete[] backwardRow_;
  delete[] rowsInBlock_;
  delete[] columnsInBlock_;
  delete[] startRowBlock_;
  delete[] startColumnBlock_;
  delete[] intsInBlock_;
  delete[] fingerPrint_;
  //functionPointer_ = NULL;
  solver_ = NULL;
  dwSolver_ = NULL;
  bestSolution_ = NULL;
  continuousSolution_ = NULL;
  fixedDj_ = NULL;
  saveLower_ = NULL;
  saveUpper_ = NULL;
  random_ = NULL;
  affinity_ = NULL;
  weights_ = NULL;
  objectiveDW_ = NULL;
  numberColumnsDW_ = NULL;
  whichRowBlock_ = NULL;
  whichColumnBlock_ = NULL;
  dwBlock_ = NULL;
  backwardRow_ = NULL;
  rowsInBlock_ = NULL;
  columnsInBlock_ = NULL;
  startRowBlock_ = NULL;
  startColumnBlock_ = NULL;
  intsInBlock_ = NULL;
  fingerPrint_ = NULL;
  numberBlocks_ = 0;
}

// Destructor
CbcHeuristicDW::~CbcHeuristicDW()
{
  gutsOfDelete();
}

// Clone
CbcHeuristic *
CbcHeuristicDW::clone() const
{
  return new CbcHeuristicDW(*this);
}

// Assignment operator
CbcHeuristicDW &
CbcHeuristicDW::operator=(const CbcHeuristicDW &rhs)
{
  if (this != &rhs) {
    CbcHeuristic::operator=(rhs);
    gutsOfDelete();
    gutsOfCopy(rhs);
  }
  return *this;
}

// Create C++ lines to get to current state
void CbcHeuristicDW::generateCpp(FILE *fp)
{
  abort();
}

// Copy constructor
CbcHeuristicDW::CbcHeuristicDW(const CbcHeuristicDW &rhs)
  : CbcHeuristic(rhs)
{
  gutsOfCopy(rhs);
}
// Resets stuff if model changes
void CbcHeuristicDW::resetModel(CbcModel *model)
{
  if (model_ && numberBlocks_ && model->getNumCols() != model->getNumCols())
    abort();
  model_ = model;
}
/*
  First tries setting a variable to better value.  If feasible then
  tries setting others.  If not feasible then tries swaps
  Returns 1 if solution, 0 if not */
int CbcHeuristicDW::solution(double &solutionValue,
  double *betterSolution)
{
  numCouldRun_++;
  int returnCode = 0;
  const double *bestSolutionIn = model_->bestSolution();
  if (!bestSolutionIn && !bestSolution_)
    return 0; // No solution found yet
  if (numberBlocks_ < 3)
    return 0; // no point
#ifdef HEURISTIC_INFORM
  printf("Entering heuristic %s - nRuns %d numCould %d when %d\n",
    heuristicName(), numRuns_, numCouldRun_, when_);
#endif
  if (bestSolutionIn && objectiveValue(bestSolutionIn) < bestObjective_ - 1.0e-5)
    passInSolution(bestSolutionIn);
  char dwPrint[200];
  sprintf(dwPrint, "Before PASSes objective is %g",
    bestObjective_);
  lastObjective_ = bestObjective_;
  model_->messageHandler()->message(CBC_FPUMP1, model_->messages())
    << dwPrint
    << CoinMessageEol;
  double startTime = CoinCpuTime();
  double startTimeElapsed = CoinGetTimeOfDay();
  CoinWarmStart *basis = NULL;
  lastObjective_ = COIN_DBL_MAX;
  int passesToDW = dwSolver_ ? 0 : -1;
  bool goodSolution = true;
  int numberColumns = solver_->getNumCols();
  int logLevel = model_->messageHandler()->logLevel();
  // For moment just OsiClp
  OsiClpSolverInterface *solver = dynamic_cast< OsiClpSolverInterface * >(solver_);
  ClpSimplex *simplex = solver->getModelPtr();
  double *columnLower = simplex->columnLower();
  double *columnUpper = simplex->columnUpper();
  const double *cost = solver->getObjCoefficients();
  const double *dj = solver->getReducedCost();
  assert(solver);
  if (!continuousSolution_) {
    bool takeHint;
    OsiHintStrength strength;
    solver->getHintParam(OsiDoDualInResolve, takeHint, strength);
    solver->setHintParam(OsiDoDualInResolve, false, OsiHintDo);
    solver->resolve();
    solver->setHintParam(OsiDoDualInResolve, takeHint, OsiHintDo);
    continuousSolution_ = CoinCopyOfArray(solver->getColSolution(),
      numberColumns);
  }
  // data arrays
  // Block order for choosing (after sort)
  int *whichBlock = new int[8 * numberBlocks_];
  memset(whichBlock, 0, 8 * numberBlocks_ * sizeof(int));
  // Count of number of times block chosen
  int *doneBlock = whichBlock + numberBlocks_;
  // Pass at which block last used
  int *whenBlock = doneBlock + numberBlocks_;
  // Number of Big Djs' (? artificial costs) in block
  int *bigDjBlock = whenBlock + numberBlocks_;
  // Number of times block has helped improve solution
  int *goodBlock = bigDjBlock + numberBlocks_;
  int *priorityBlock = goodBlock + numberBlocks_;
  int *orderBlock = priorityBlock + numberBlocks_;
  // block can be fixed if nothing in master rows, maybe always same as continuous
  int *fixedBlock = orderBlock + numberBlocks_;
  // Mixture of stuff to sort blocks on
  double *blockSort = new double[4 * numberBlocks_];
  // Reduced cost (d sub j was old notation) contribution
  double *blockDj = blockSort + numberBlocks_;
  // Difference between current best and continuous solutions
  double *blockDiff = blockDj + numberBlocks_;
  // Difference between current best and continuous solutions (just integers)
  double *blockDiffInt = blockDiff + numberBlocks_;
  delete[] fixedDj_;
  fixedDj_ = CoinCopyOfArray(dj, numberColumns);
  int numberImproving = 0;
  int *whenBetter = new int[numberPasses_];
  double *improvement = new double[numberPasses_];
  // First has number
  int **improvingBlocks = new int *[numberPasses_];
  int numberBlocksUsed = numberBlocks_;
  // Get basic priority order
  for (int i = 0; i < numberBlocks_; i++) {
    whenBlock[i] = -intsInBlock_[i];
    orderBlock[i] = i;
  }
  CoinSort_2(whenBlock, whenBlock + numberBlocks_, orderBlock);
  for (int i = 0; i < numberBlocks_; i++) {
    whenBlock[i] = -100;
    doneBlock[i] = 0;
    whichBlock[i] = i;
  }
  bestObjective_ = objectiveValue(bestSolution_);
  double startTime2 = CoinCpuTime();
  double startTime2Elapsed = CoinGetTimeOfDay();
  // make sure all master columns freed up
  for (int iColumn = 0; iColumn < numberColumns; ++iColumn) {
    if (whichColumnBlock_[iColumn] < 0) {
      columnLower[iColumn] = saveLower_[iColumn];
      columnUpper[iColumn] = saveUpper_[iColumn];
    }
  }
  int numberNoMaster = 0;
  int numberSameAsContinuous = 0;
  int numberSameAsContinuousJustInts = 0;
  // Column copy
  //const double * element = solver->getMatrixByCol()->getElements();
  const int *row = solver->getMatrixByCol()->getIndices();
  const CoinBigIndex *columnStart = solver->getMatrixByCol()->getVectorStarts();
  const int *columnLength = solver->getMatrixByCol()->getVectorLengths();
  for (int iBlock = 0; iBlock < numberBlocks_; iBlock++) {
    int nElInMaster = 0;
    int numberDifferentContinuous = 0;
    int numberDifferentContinuousJustInts = 0;
    int start = startColumnBlock_[iBlock];
    int end = startColumnBlock_[iBlock + 1];
    for (int i = start; i < end; i++) {
      int iColumn = columnsInBlock_[i];
      for (CoinBigIndex j = columnStart[iColumn];
           j < columnStart[iColumn] + columnLength[iColumn]; j++) {
        int iRow = row[j];
        iRow = backwardRow_[iRow];
        if (iRow >= 0)
          nElInMaster++;
      }
      if (fabs(bestSolution_[iColumn] - continuousSolution_[iColumn]) < 1.0e-5) {
        numberDifferentContinuous++;
        if (solver->isInteger(iColumn))
          numberDifferentContinuousJustInts++;
      }
    }
    if (!nElInMaster) {
      fixedBlock[iBlock] = 10;
      numberNoMaster++;
    } else if (!numberDifferentContinuous) {
      fixedBlock[iBlock] = 2;
      numberSameAsContinuous++;
    } else if (!numberDifferentContinuousJustInts) {
      fixedBlock[iBlock] = 1;
      numberSameAsContinuousJustInts++;
    }
  }
  if (numberNoMaster) {
    sprintf(dwPrint, "*** %d blocks have no elements in master - can be solved separately",
      numberNoMaster);
    model_->messageHandler()->message(CBC_FPUMP1, model_->messages())
      << dwPrint
      << CoinMessageEol;
  }
  sprintf(dwPrint, "With initial best solution %d blocks were same as continuous, %d when just looking at integers",
    numberSameAsContinuous, numberSameAsContinuousJustInts);
  model_->messageHandler()->message(CBC_FPUMP1, model_->messages())
    << dwPrint
    << CoinMessageEol;
  int lastGoodPass = 0;
  for (pass_ = 0; pass_ < numberPasses_; pass_++) {
    double endTime2 = CoinCpuTime();
    double endTime2Elapsed = CoinGetTimeOfDay();
#ifndef SCALE_FACTOR
#define SCALE_FACTOR 1.0
#endif
    if (pass_) {
      sprintf(dwPrint, "PASS %d changed objective from %g to %g in %g seconds (%g elapsed) - total %g (%g elapsed) - current needed %d nodes %d - %s",
        pass_, lastObjective_ * SCALE_FACTOR,
        bestObjective_ * SCALE_FACTOR, endTime2 - startTime2,
        endTime2Elapsed - startTime2Elapsed,
        endTime2 - startTime, endTime2Elapsed - startTimeElapsed,
        nNeeded_, nNodes_,
        lastObjective_ > bestObjective_ + 1.0e-3 ? "improving" : "");
      model_->messageHandler()->message(CBC_FPUMP1, model_->messages())
        << dwPrint
        << CoinMessageEol;
    }
    if ((pass_ % 10) == 9) {
      for (int iImp = CoinMax(1, numberImproving - 10); iImp < numberImproving; iImp++) {
        int *blocks = improvingBlocks[iImp];
        int nBlocks = blocks[0];
        blocks++;
        sprintf(dwPrint, "Pass %d improved objective by %g using %d blocks - ",
          whenBetter[iImp], improvement[iImp], nBlocks);
        model_->messageHandler()->message(CBC_FPUMP1, model_->messages())
          << dwPrint
          << CoinMessageEol;
        if (logLevel > 1) {
          for (int i = 0; i < nBlocks; i++)
            printf("%d ", blocks[i]);
          printf("\n");
        }
      }
      if (logLevel > 1) {
        int *count = new int[numberImproving + 1];
        memset(count, 0, (numberImproving + 1) * sizeof(int));
        for (int i = 0; i < numberBlocks_; i++)
          count[goodBlock[i]]++;
        for (int i = 0; i < numberImproving; i++) {
          if (count[i])
            printf("%d blocks were involved in improvement %d times\n",
              count[i], i);
        }
        delete[] count;
      }
    }
    startTime2 = CoinCpuTime();
    startTime2Elapsed = CoinGetTimeOfDay();
    if (pass_ - lastGoodPass > numberBadPasses_) {
      sprintf(dwPrint, "Exiting on lack of progress");
      model_->messageHandler()->message(CBC_FPUMP1, model_->messages())
        << dwPrint
        << CoinMessageEol;
      break;
    }
    if (model_->getNodeCount() >= model_->getMaximumNodes() || model_->maximumSecondsReached()) {
      sprintf(dwPrint, "Exiting on time or interrupt");
      model_->messageHandler()->message(CBC_FPUMP1, model_->messages())
        << dwPrint
        << CoinMessageEol;
      break;
    }
    if (bestObjective_ >= lastObjective_ - 1.0e-3) {
      // what now
      // 0 - fine, 1 can't be better, 2 max node
      //assert(solveState);
      if (solveState_ < 2) {
        // more in
        sprintf(dwPrint, "No improvement - think we need more variables ");
        model_->messageHandler()->message(CBC_FPUMP1, model_->messages())
          << dwPrint
          << CoinMessageEol;
        nNeeded_ += nNeeded_ / 10;
        nNeeded_ = CoinMin(nNeeded_, 800);
        nNodes_ = nNodesBase_;
        (*(functionPointer_))(this, NULL, 6);
      } else {
        // more nodes fewer in
        sprintf(dwPrint, "No improvement - stopped on nodes - think we need more nodes ");
        model_->messageHandler()->message(CBC_FPUMP1, model_->messages())
          << dwPrint
          << CoinMessageEol;
        if (phase_) {
          nNodes_ += nNodes_ / 5;
          nNodes_ = CoinMin(nNodes_, 1000);
        }
        nNeeded_ -= nNeeded_ / 20;
        nNeeded_ = CoinMax(nNeeded_, 50);
        (*(functionPointer_))(this, NULL, 7);
      }
    } else {
      // improving
      (*(functionPointer_))(this, NULL, 4);
      solveState_ = 0;
      //lastObjective_=bestObjective_;
      if (phase_) {
        //nNeededBase_ += nNeededBase_/50;
        //nNodesBase_ += nNodesBase_/50;
      }
      nNeeded_ -= nNeeded_ / 10;
      nNeeded_ = CoinMax(nNeededBase_, nNeeded_);
      nNodes_ = nNodesBase_;
      (*(functionPointer_))(this, NULL, 8);
    }
    sprintf(dwPrint, "new needed %d, nodes %d", nNeeded_, nNodes_);
    model_->messageHandler()->message(CBC_FPUMP1, model_->messages())
      << dwPrint
      << CoinMessageEol;
    for (int i = 0; i < numberColumns; ++i) {
      if (solver_->isInteger(i)) {
        double value = floor(bestSolution_[i] + 0.5);
        columnLower[i] = value;
        columnUpper[i] = value;
      } else {
        columnLower[i] = saveLower_[i];
        columnUpper[i] = saveUpper_[i];
      }
    }
    if (goodSolution) {
      lastGoodPass = pass_;
      int lastNumberDW = numberDW_;
      solver->setColSolution(bestSolution_);
      solver->setHintParam(OsiDoReducePrint, false, OsiHintDo, 0);
      if (basis) {
        solver->setWarmStart(basis);
        delete basis;
        basis = NULL;
      }
      solver->resolve();
      solver->setHintParam(OsiDoReducePrint, true, OsiHintDo, 0);
      if (solver->getObjValue() < lastObjective_ - 1.0e-5) {
        memcpy(bestSolution_, solver->getColSolution(),
          numberColumns * sizeof(double));
        bestObjective_ = solver->getObjValue();
        int *blocks = new int[numberBlocksUsed + 1];
        blocks[0] = numberBlocksUsed;
        memcpy(blocks + 1, whichBlock, numberBlocksUsed * sizeof(int));
        improvingBlocks[numberImproving] = blocks;
        whenBetter[numberImproving] = pass_;
        improvement[numberImproving] = lastObjective_ - bestObjective_;
        numberImproving++;
        lastObjective_ = bestObjective_;
        if (pass_) {
          // update good
          for (int i = 0; i < numberBlocksUsed; i++)
            goodBlock[whichBlock[i]]++;
        }
      }
      goodSolution = false;
      if (fullDWEverySoOften_ > 0) {
        addDW(bestSolution_, numberBlocksUsed,
          whichBlock);
      }
      if (passesToDW == 0) {
        passesToDW = fullDWEverySoOften_;
        const double *duals = solver->getRowPrice();
        double *bestSolution2 = CoinCopyOfArray(bestSolution_,
          numberColumns);
        // Column copy
        const double *element = solver->getMatrixByCol()->getElements();
        const int *row = solver->getMatrixByCol()->getIndices();
        const CoinBigIndex *columnStart = solver->getMatrixByCol()->getVectorStarts();
        const int *columnLength = solver->getMatrixByCol()->getVectorLengths();
        int numberUsed = 0;
        for (int iBlock = 0; iBlock < numberBlocks_; iBlock++) {
          int start = startColumnBlock_[iBlock];
          int end = startColumnBlock_[iBlock + 1];
          ClpSimplex *tempModel = new ClpSimplex(solver->getModelPtr(),
            startRowBlock_[iBlock + 1] - startRowBlock_[iBlock],
            rowsInBlock_ + startRowBlock_[iBlock],
            end - start,
            columnsInBlock_ + startColumnBlock_[iBlock]);
          tempModel->setLogLevel(0);
          tempModel->setDualObjectiveLimit(COIN_DBL_MAX);
          double *objectiveX = tempModel->objective();
          double *columnLowerX = tempModel->columnLower();
          double *columnUpperX = tempModel->columnUpper();
          for (int i = start; i < end; i++) {
            int jColumn = i - start;
            int iColumn = columnsInBlock_[i];
            columnLowerX[jColumn] = CoinMax(saveLower_[iColumn], -1.0e12);
            columnUpperX[jColumn] = CoinMin(saveUpper_[iColumn], 1.0e12);
            if (solver->isInteger(iColumn))
              tempModel->setInteger(jColumn);
            double cost = objectiveX[jColumn];
            for (CoinBigIndex j = columnStart[iColumn];
                 j < columnStart[iColumn] + columnLength[iColumn]; j++) {
              int iRow = row[j];
              double elementValue = element[j];
              if (backwardRow_[iRow] >= 0) {
                cost -= elementValue * duals[iRow];
              }
            }
            objectiveX[jColumn] = cost;
          }
          OsiClpSolverInterface solverX(tempModel, true);
          CbcModel modelX(solverX);
          modelX.setLogLevel(1);
          modelX.setMoreSpecialOptions2(57);
          // need to stop after solutions and nodes
          //modelX.setMaximumNodes(nNodes_);
          modelX.setMaximumSolutions(1);
          modelX.branchAndBound();
          const double *bestSolutionX = modelX.bestSolution();
          if (bestSolutionX) {
            whichBlock[numberUsed++] = iBlock;
            for (int i = start; i < end; i++) {
              int jColumn = i - start;
              int iColumn = columnsInBlock_[i];
              bestSolution2[iColumn] = bestSolutionX[jColumn];
            }
          }
        }
        addDW(bestSolution2, numberUsed, whichBlock);
        if (!pass_ && false) {
          // see if gives a solution
          for (int i = 0; i < numberColumns; ++i) {
            if (solver_->isInteger(i)) {
              double value = floor(bestSolution2[i] + 0.5);
              columnLower[i] = value;
              columnUpper[i] = value;
            } else {
              columnLower[i] = saveLower_[i];
              columnUpper[i] = saveUpper_[i];
            }
          }
          solver_->resolve();
          if (solver_->isProvenOptimal()) {
            printf("DW1 sol %g\n", solver->getObjValue());
          }
        }
        // now try purer DW
        bool takeHint;
        OsiHintStrength strength;
        dwSolver_->getHintParam(OsiDoDualInResolve, takeHint, strength);
        dwSolver_->setHintParam(OsiDoDualInResolve, false, OsiHintDo);
        dwSolver_->resolve();
        dwSolver_->setHintParam(OsiDoDualInResolve, takeHint, OsiHintDo);
        duals = dwSolver_->getRowPrice();
        numberUsed = 0;
        for (int iBlock = 0; iBlock < numberBlocks_; iBlock++) {
          int start = startColumnBlock_[iBlock];
          int end = startColumnBlock_[iBlock + 1];
          ClpSimplex *tempModel = new ClpSimplex(solver->getModelPtr(),
            startRowBlock_[iBlock + 1] - startRowBlock_[iBlock],
            rowsInBlock_ + startRowBlock_[iBlock],
            end - start,
            columnsInBlock_ + startColumnBlock_[iBlock]);
          tempModel->setLogLevel(0);
          tempModel->setDualObjectiveLimit(COIN_DBL_MAX);
          double *objectiveX = tempModel->objective();
          double *columnLowerX = tempModel->columnLower();
          double *columnUpperX = tempModel->columnUpper();
          double convexityDual = duals[numberMasterRows_ + iBlock];
          for (int i = start; i < end; i++) {
            int jColumn = i - start;
            int iColumn = columnsInBlock_[i];
            columnLowerX[jColumn] = CoinMax(saveLower_[iColumn], -1.0e12);
            columnUpperX[jColumn] = CoinMin(saveUpper_[iColumn], 1.0e12);
            if (solver->isInteger(iColumn))
              tempModel->setInteger(jColumn);
            double cost = objectiveX[jColumn];
            for (CoinBigIndex j = columnStart[iColumn];
                 j < columnStart[iColumn] + columnLength[iColumn]; j++) {
              int iRow = row[j];
              double elementValue = element[j];
              if (backwardRow_[iRow] >= 0) {
                // duals are from dw
                cost -= elementValue * duals[backwardRow_[iRow]];
              }
            }
            objectiveX[jColumn] = cost;
          }
          OsiClpSolverInterface solverX(tempModel, true);
          solverX.initialSolve();
          double cObj = solverX.getObjValue();
          CbcModel modelX(solverX);
          modelX.setLogLevel(1);
          modelX.setMoreSpecialOptions2(57);
          modelX.setMaximumSolutions(1);
          modelX.branchAndBound();
          sprintf(dwPrint, "Block %d contobj %g intobj %g convdual %g",
            iBlock, cObj, modelX.getObjValue(), convexityDual);
          model_->messageHandler()->message(CBC_FPUMP2, model_->messages())
            << dwPrint
            << CoinMessageEol;
          const double *bestSolutionX = modelX.bestSolution();
          if (bestSolutionX) {
            whichBlock[numberUsed++] = iBlock;
            for (int i = start; i < end; i++) {
              int iColumn = columnsInBlock_[i];
              bestSolution2[iColumn] = bestSolutionX[i - start];
            }
          }
        }
        addDW(bestSolution2, numberUsed, whichBlock);
        if (!pass_ && false) {
          // see if gives a solution
          for (int i = 0; i < numberColumns; ++i) {
            if (solver_->isInteger(i)) {
              double value = floor(bestSolution2[i] + 0.5);
              columnLower[i] = value;
              columnUpper[i] = value;
            } else {
              columnLower[i] = saveLower_[i];
              columnUpper[i] = saveUpper_[i];
            }
          }
          solver_->resolve();
          if (solver_->isProvenOptimal()) {
            printf("DW sol %g\n", solver->getObjValue());
          }
        }
        delete[] bestSolution2;
      }
      passesToDW--;
      if (numberDW_ > lastNumberDW) {
        intArray_ = NULL;
        doubleArray_ = NULL;
        (*(functionPointer_))(this, NULL, 3);
      }
    }
    for (int i = 0; i < numberBlocks_; i++) {
      blockDj[i] = 0.0;
      blockDiff[i] = 0.0;
      blockDiffInt[i] = 0.0;
      bigDjBlock[i] = 0;
    }
    for (int i = 0; i < numberColumns; ++i) {
      int kBlock = whichColumnBlock_[i];
      if (kBlock >= 0) {
        columnLower[i] = bestSolution_[i];
        columnUpper[i] = bestSolution_[i];
      }
    }
    for (int iBlock = 0; iBlock < numberBlocks_; iBlock++) {
      for (int i = startColumnBlock_[iBlock];
           i < startColumnBlock_[iBlock + 1]; i++) {
        int iColumn = columnsInBlock_[i];
        if (continuousSolution_) {
          blockDiff[iBlock] += fabs((bestSolution_[iColumn] - continuousSolution_[iColumn]))
            * (fabs(cost[iColumn]) + 1.0e-5);
          if (solver->isInteger(iColumn))
            blockDiffInt[iBlock] += fabs((bestSolution_[iColumn] - continuousSolution_[iColumn]))
              * (fabs(cost[iColumn]) + 1.0e-5);
        }
        if (solver->isInteger(iColumn)) {
          if (bestSolution_[iColumn] < saveUpper_[iColumn] - 1.0e-1) {
            if (fixedDj_[iColumn] < -1.0e-5)
              blockDj[iBlock] += fixedDj_[iColumn];
            if (fixedDj_[iColumn] < -1.0e4)
              bigDjBlock[iBlock]++;
          }
          if (bestSolution_[iColumn] > saveLower_[iColumn] + 1.0e-1) {
            if (fixedDj_[iColumn] > 1.0e-5)
              blockDj[iBlock] -= fixedDj_[iColumn];
            if (fixedDj_[iColumn] > 1.0e4)
              bigDjBlock[iBlock]++;
          }
        }
      }
    }
    // Get average dj and difference
    int numberInDj = 0;
    double averageDj = 1.0e-12;
    int numberInDiff = 0;
    double averageDiff = 1.0e-12;
    for (int i = 0; i < numberBlocks_; i++) {
      assert(blockDj[i] <= 0.0);
      if (blockDj[i] < 0.0) {
        numberInDj++;
        averageDj -= blockDj[i];
      }
      assert(blockDiff[i] >= 0.0);
      if (blockDiff[i] > 0.0) {
        numberInDiff++;
        averageDiff += blockDiff[i];
      }
    }
    if (numberInDj)
      averageDj /= static_cast< double >(numberInDj);
    if (numberInDiff)
      averageDiff /= static_cast< double >(numberInDiff);
    double ratioDiff = averageDj / averageDiff;
    // downplay
    ratioDiff *= 1.0e-3;
    for (int i = 0; i < numberBlocks_; i++) {
      whichBlock[i] = i;
      assert(intsInBlock_[i]);
      blockSort[i] = blockDj[i];
#define TRY_ADJUST 3
#if TRY_ADJUST == 1
      blockSort[i] *= CoinDrand48() + 5.0e-1;
#elif TRY_ADJUST == 2
      blockSort[i] -= ratioDiff * blockDiff[i];
#elif TRY_ADJUST == 3
      blockSort[i] -= ratioDiff * blockDiff[i];
      blockSort[i] *= 0.1 * CoinDrand48() + 0.9;
#endif
      if (phase_ == 99)
        blockSort[i] -= 2.0 * averageDj * goodBlock[i];
      blockSort[i] /= static_cast< double >(intsInBlock_[i]);
      //blockSort[i] /= sqrt(static_cast<double>(intsInBlock_[i]));
      if (doneBlock[i]) {
        blockSort[i] /= static_cast< double >(doneBlock[i] + 1);
        if (whenBlock[i] > pass_ - 10)
          blockSort[i] += 1.0e2 * averageDj;
      }
    }
    CoinSort_2(blockSort, blockSort + numberBlocks_, whichBlock);
    // allow user to modify
    intArray_ = whichBlock;
    doubleArray_ = blockSort;
    (*(functionPointer_))(this, NULL, 1);
    int numberBlocksIn = 0;
    for (int iTry = 0; iTry < 2; iTry++) {
      int nFreed = 0;
      int nBigDjBlock = 0;
      numberBlocksUsed = 0;
      for (int i = 0; i < numberBlocks_; i++) {
        int iBlock = whichBlock[i];
        bool skipBlock = false;
        assert(iBlock >= 0);
        if ((doneBlock[iBlock] && !phase_) || !blockSort[i]) {
          //printf("already done block %d - dj %g\n",iBlock,blockDj[i]);
          skipBlock = true;
        } else if (bigDjBlock[iBlock]) {
          nBigDjBlock++;
          if (nBigDjBlock > 20 && !phase_) {
            skipBlock = true;
          }
        }
        int returnCode = (*(functionPointer_))(this, NULL, 5);
        if (returnCode < 0)
          skipBlock = true;
        else if (returnCode > 0)
          skipBlock = false;
        if (skipBlock) {
          whichBlock[i] -= 1000000;
        } else {
          // free up
          sprintf(dwPrint, "freeing block %d (already freed %d times) - dj %g, diff %g, diffint %g - %d columns (%d integer)",
            iBlock, doneBlock[iBlock], blockDj[iBlock],
            blockDiff[iBlock], blockDiffInt[iBlock],
            startColumnBlock_[iBlock + 1] - startColumnBlock_[iBlock],
            intsInBlock_[iBlock]);
          model_->messageHandler()->message(CBC_FPUMP2, model_->messages())
            << dwPrint
            << CoinMessageEol;
          numberBlocksIn++;
          doneBlock[iBlock]++;
          whenBlock[iBlock] = pass_;
          nFreed += intsInBlock_[iBlock];
          for (int j = startColumnBlock_[iBlock];
               j < startColumnBlock_[iBlock + 1]; j++) {
            int iColumn = columnsInBlock_[j];
            columnLower[iColumn] = saveLower_[iColumn];
            columnUpper[iColumn] = saveUpper_[iColumn];
          }
          if (!numberBlocksUsed && affinity_) {
            // re-sort rest using affinity
            const unsigned short int *aff = affinity_ + iBlock * numberBlocks_;
            for (int j = i + 1; j < numberBlocks_; j++) {
              int jBlock = whichBlock[j];
              blockSort[j] -= 1.0e2 * aff[jBlock];
            }
            CoinSort_2(blockSort + i + 1, blockSort + numberBlocks_, whichBlock + i + 1);
          }
          numberBlocksUsed = i + 1;
          if (nFreed >= nNeeded_ && numberBlocksIn > 3)
            break;
        }
      }
      sprintf(dwPrint, "%d big dj blocks found", nBigDjBlock);
      model_->messageHandler()->message(CBC_FPUMP2, model_->messages())
        << dwPrint
        << CoinMessageEol;
      if (nFreed)
        break;
      phase_ = 1; // round again
      sprintf(dwPrint, "Changing phase");
      model_->messageHandler()->message(CBC_FPUMP1, model_->messages())
        << dwPrint
        << CoinMessageEol;
      for (int i = 0; i < numberBlocks_; i++)
        whichBlock[i] += 1000000;
      //nNeeded=500; // allow more
    }
    if (numberBlocksIn == numberBlocks_)
      numberPasses_ = 1; // no point
    // pack down blocks
    int n = 0;
    for (int i = 0; i < numberBlocksUsed; i++) {
      if (whichBlock[i] >= 0)
        whichBlock[n++] = whichBlock[i];
    }
    numberBlocksUsed = n;
    int nFixedInts = 0;
    int nFixedContinuous = 0;
    int nFreeInts = 0;
    int nFreeContinuous = 0;
    int nFreeMaster = 0;
    int nFixedMaster = 0;
    for (int i = 0; i < numberColumns; ++i) {
      int kBlock = whichColumnBlock_[i];
      if (columnUpper[i] > columnLower[i]) {
        if (kBlock >= 0) {
          if (solver->isInteger(i))
            nFreeInts++;
          else
            nFreeContinuous++;
        } else {
          nFreeMaster++;
        }
      } else {
        if (kBlock >= 0) {
          if (solver->isInteger(i))
            nFixedInts++;
          else
            nFixedContinuous++;
        } else {
          nFixedMaster++;
        }
      }
    }
    sprintf(dwPrint, "Fixed %d ints, %d c, %d m - free %d, %d, %d",
      nFixedInts, nFixedContinuous, nFixedMaster,
      nFreeInts, nFreeContinuous, nFreeMaster);
    model_->messageHandler()->message(CBC_FPUMP1, model_->messages())
      << dwPrint
      << CoinMessageEol;
    // But free up and then fix again
    for (int i = 0; i < numberColumns; ++i) {
      columnLower[i] = saveLower_[i];
      columnUpper[i] = saveUpper_[i];
      int kBlock = whichColumnBlock_[i];
      if (kBlock >= 0 && solver->isInteger(i) && whenBlock[kBlock] != pass_) {
        columnLower[i] = bestSolution_[i];
        columnUpper[i] = bestSolution_[i];
      }
    }
    bool takeHint;
    OsiHintStrength strength;
    solver->getHintParam(OsiDoDualInResolve, takeHint, strength);
    solver->setHintParam(OsiDoDualInResolve, false, OsiHintDo);
    solver->messageHandler()->setLogLevel(1);
    solver->setHintParam(OsiDoReducePrint, false, OsiHintDo, 0);
    solver->resolve();
    solver->setHintParam(OsiDoReducePrint, true, OsiHintDo, 0);
    //solver->messageHandler()->setLogLevel(0) ;
    solver->setHintParam(OsiDoDualInResolve, takeHint, strength);
    if (solver->getObjValue() > bestObjective_ + 1.0e-5 * (1.0 + fabs(bestObjective_))) {
      // trouble
      if (logLevel > 1) {
        for (int i = 0; i < numberBlocks_; i++) {
          if (whenBlock[i] == pass_) {
            printf("Block %d free\n", i);
          }
        }
      }
      solver->writeMps("bad", "mps");
      const double *lower = solver->getColLower();
      const double *upper = solver->getColUpper();
      if (logLevel > 1) {
        printf("best obj %g\n", objectiveValue(bestSolution_));
        for (int i = 0; i < numberColumns; ++i) {
          double value = bestSolution_[i];
          if (value < lower[i] - 1.0e-5 || value > upper[i] + 1.0e-5)
            printf("column %d (block %d) %g %g <= %g <= %g %g\n",
              i, whichColumnBlock_[i], saveLower_[i],
              lower[i], value, upper[i], saveUpper_[i]);
        }
      }
      //abort();
      sprintf(dwPrint, "**** Continuous below best - bad solution passed in?!");
      model_->messageHandler()->message(CBC_FPUMP1, model_->messages())
        << dwPrint
        << CoinMessageEol;
      // give up
      break;
    }
    const double *tempSol = solver->getColSolution();
    // But free up and then fix again
    for (int i = 0; i < numberColumns; ++i) {
      int kBlock = whichColumnBlock_[i];
      if (kBlock >= 0 && !solver->isInteger(i) && whenBlock[kBlock] != pass_) {
        columnLower[i] = tempSol[i];
        columnUpper[i] = tempSol[i];
      }
    }
    solver->getHintParam(OsiDoDualInResolve, takeHint, strength);
    solver->setHintParam(OsiDoDualInResolve, false, OsiHintDo);
    solver->messageHandler()->setLogLevel(1);
    solver->setHintParam(OsiDoReducePrint, false, OsiHintDo, 0);
    solver->resolve();
    solver->setHintParam(OsiDoReducePrint, true, OsiHintDo, 0);
    //solver->messageHandler()->setLogLevel(0) ;
    solver->setHintParam(OsiDoDualInResolve, takeHint, strength);
    //solver->messageHandler()->setLogLevel(0) ;
    //lp->setLogLevel(0);
    ClpPresolve pinfo;
    // fix small infeasibilities
    pinfo.setPresolveActions(pinfo.presolveActions() | 0x4000);
    int numberPasses = 2; // can change this
    ClpSimplex *model2 = pinfo.presolvedModel(*simplex, 1.0e-8,
      true,
      numberPasses, true);
    if (!model2) {
      abort();
    } else {
      sprintf(dwPrint, "Reduced model has %d rows and %d columns",
        model2->numberRows(), model2->numberColumns());
      model_->messageHandler()->message(CBC_FPUMP1, model_->messages())
        << dwPrint
        << CoinMessageEol;
      //model2->setLogLevel(0);
      OsiClpSolverInterface solver2(model2);
      solver2.setWarmStart(NULL);
      //solver2.messageHandler()->setLogLevel(0) ;
      CbcModel model(solver2);
      model.setMaximumNodes(nNodes_);
      //model.setMaximumSolutions(2);
      model.solver()->setHintParam(OsiDoReducePrint, false, OsiHintDo, 0);
      model.solver()->setHintParam(OsiDoPresolveInInitial, false, OsiHintDo, 0);
      //model.solver()->setHintParam(OsiDoPresolveInResolve, false, OsiHintDo, 0) ;
      model.initialSolve();
      if (bestObjective_ > model.solver()->getObjValue() + 1.0e-1) {
        int nFix = 0;
#ifdef HOT_START
        // Set up hot start
        const double *dj = model.solver()->getReducedCost();
        const double *lower = model.solver()->getColLower();
        const double *upper = model.solver()->getColUpper();
        const double *solution = model.solver()->getColSolution();
        double gap = CoinMax(bestObjective_ - model.solver()->getObjValue(),
          1.0e-3);
        int numberColumns2 = model.solver()->getNumCols();
        const int *originalColumns = pinfo.originalColumns();
        double *hot = new double[2 * numberColumns2];
        int *hotPriorities = new int[numberColumns2 * 2];
        double *hotWeight = hot + numberColumns2;
        memset(hot, 0, numberColumns2 * sizeof(double));
        int *sort = hotPriorities + numberColumns2;
        for (int i = 0; i < numberColumns2; i++) {
          hotWeight[i] = COIN_DBL_MAX;
          sort[i] = i;
          if (solver2.isInteger(i)) {
            int iColumn = originalColumns[i];
            hot[i] = bestSolution_[iColumn];
            hotWeight[i] = -fabs(fixedDj_[iColumn]);
            if (bestSolution_[i] < saveLower_[iColumn] + 1.0e-6) {
              if (fixedDj_[i] > 0.0)
                hotWeight[i] = fixedDj_[iColumn];
            } else if (bestSolution_[i] > saveUpper_[iColumn] - 1.0e-6) {
              if (fixedDj_[i] < 0.0)
                hotWeight[i] = -fixedDj_[iColumn];
            }
            if (solution[i] < saveLower_[iColumn] + 1.0e-6) {
              if (dj[i] > gap) {
                solver2.setColUpper(i, saveLower_[iColumn]);
                nFix++;
              }
            } else if (solution[i] > saveUpper_[iColumn] - 1.0e-6) {
              if (-dj[i] > gap) {
                solver2.setColLower(i, saveUpper_[iColumn]);
                nFix++;
              }
            }
          }
        }
        CoinSort_2(hotWeight, hotWeight + numberColumns2, sort);
        for (int i = 0; i < numberColumns2; i++) {
          hotPriorities[sort[i]] = i + 1;
        }
        model.setHotstartSolution(hot, hotPriorities);
        delete[] hot;
        delete[] hotPriorities;
#endif
        if (logLevel > 1) {
          if (nFix)
            printf("Fixed another %d integers\n", nFix);
        }
        {
          // priorities
          memset(priorityBlock, 0, numberBlocks_ * sizeof(int));
          for (int i = 0; i < numberBlocks_; i++) {
            int iBlock = whichBlock[i];
            if (iBlock >= 0)
              priorityBlock[iBlock] = i + 1;
          }
#if 1
          assert(numberBlocks_ < 4000);
          // But make sure we do one block before next
          for (int i = 0; i < numberBlocks_; i++) {
            priorityBlock[i] += 4000 * orderBlock[i];
          }
#else
          // just do ones with many first
          for (int i = 0; i < numberBlocks; i++) {
            priorityBlock[i] = 10000 - intsInBlock_[i];
          }
#endif
          int numberColumns2 = model.solver()->getNumCols();
          const int *original = pinfo.originalColumns();
          int *priorities = new int[numberColumns2];
          int n = 0;
          for (int i = 0; i < numberColumns2; i++) {
            if (solver2.isInteger(i)) {
              int iBlock = whichColumnBlock_[original[i]];
              priorities[n++] = priorityBlock[iBlock];
            }
          }
          model.passInPriorities(priorities, false);
          delete[] priorities;
        }
        CglProbing probingGen;
        probingGen.setUsingObjective(1);
        probingGen.setMaxPass(1);
        probingGen.setMaxPassRoot(1);
        // Number of unsatisfied variables to look at
        probingGen.setMaxProbe(10);
        probingGen.setMaxProbeRoot(50);
        // How far to follow the consequences
        probingGen.setMaxLook(10);
        probingGen.setMaxLookRoot(50);
        probingGen.setMaxLookRoot(10);
        // Only look at rows with fewer than this number of elements
        probingGen.setMaxElements(200);
        probingGen.setMaxElementsRoot(300);
        probingGen.setRowCuts(0);
        model.addCutGenerator(&probingGen, 5, "Probing");
        model.setNumberThreads(model_->getNumberThreads());
        model.solver()->setHintParam(OsiDoReducePrint, true, OsiHintDo, 0);
        model.setPrintingMode(1);
        model.setCutoff(lastObjective_ - model_->getCutoffIncrement());
        model.setLogLevel(1);
        intArray_ = NULL;
        doubleArray_ = NULL;
        (*(functionPointer_))(this, &model, 2);
        model.branchAndBound();
        if (logLevel > 1) {
          printf("After B&B status %d objective %g\n", model.status(),
            model.getMinimizationObjValue());
        }
        int modelStatus = model.status();
        model.solver()->setHintParam(OsiDoReducePrint, false, OsiHintDo, 0);
        if (model.bestSolution() && model.getMinimizationObjValue() < lastObjective_) {
          const int *original = pinfo.originalColumns();
          int n2 = solver2.getNumCols();
          const double *bestSolution2 = model.bestSolution();
          for (int i = 0; i < n2; i++) {
            int iColumn = original[i];
            if (solver2.isInteger(i)) {
              double value = floor(bestSolution2[i] + 0.5);
              if (fabs(bestSolution2[i] - value) > 1.0e-5) {
                if (logLevel > 1) {
                  printf("bad %d %g\n", i, bestSolution2[i]);
                }
              } else {
                solver->setColLower(iColumn, value);
                solver->setColUpper(iColumn, value);
              }
            }
          }
          pinfo.postsolve(true);
          delete model2;
          bool takeHint;
          OsiHintStrength strength;
          solver->getHintParam(OsiDoDualInResolve, takeHint, strength);
          solver->setHintParam(OsiDoDualInResolve, false, OsiHintDo);
          solver->messageHandler()->setLogLevel(1);
          solver->setHintParam(OsiDoReducePrint, false, OsiHintDo, 0);
          solver->resolve();
          if (basis)
            delete basis;
          basis = solver->getWarmStart();
          memcpy(fixedDj_, solver->getReducedCost(),
            numberColumns * sizeof(double));
          solver->setHintParam(OsiDoReducePrint, true, OsiHintDo, 0);
          //solver->messageHandler()->setLogLevel(0) ;
          solver->setHintParam(OsiDoDualInResolve, takeHint, strength);
          //solver->
          CbcModel model(*solver);
          model.setNumberBeforeTrust(-1);
          model.setNumberStrong(0);
          model.setMoreSpecialOptions2(57);
          model.branchAndBound();
          if (model.getMinimizationObjValue() < lastObjective_) {
            memcpy(bestSolution_, model.bestSolution(),
              numberColumns * sizeof(double));
            bestObjective_ = model.getObjValue();
            if (logLevel > 1) {
              for (int i = 0; i < numberColumns; i++) {
                if (simplex->isInteger(i)) {
                  if (fabs(bestSolution_[i] - floor(bestSolution_[i] + 0.5)) > 1.0e-5) {
                    printf("bad after %d %g\n", i, bestSolution_[i]);
                  }
                }
              }
            }
            goodSolution = true;
          }
        } else if (modelStatus != 0) {
          // stopped on nodes
          // 0 - fine, 1 can't be better, 2 max node
          solveState_ = 2;
          pinfo.destroyPresolve();
          delete model2;
        } else {
          // complete search
          solveState_ = 1;
          pinfo.destroyPresolve();
          delete model2;
        }
      } else {
        // can't be better
        // 0 - fine, 1 can't be better, 2 max node
        solveState_ = 1;
        pinfo.destroyPresolve();
        delete model2;
      }
    }
  }
  delete[] whichBlock;
  delete[] blockSort;
  delete[] whenBetter;
  delete[] improvement;
  for (int i = 0; i < numberImproving; i++)
    delete[] improvingBlocks[i];
  delete[] improvingBlocks;
  delete basis;
  if (bestObjective_ < solutionValue) {
    solutionValue = bestObjective_;
    memcpy(betterSolution, bestSolution_,
      solver_->getNumCols() * sizeof(double));
  }
  return returnCode;
}
// update model
void CbcHeuristicDW::setModel(CbcModel *model)
{
  if (model != model_) {
    gutsOfDelete();
    model_ = model;
    assert(model->solver());
    solver_ = model->solver()->clone();
    findStructure();
  }
}
// Find structure
void CbcHeuristicDW::findStructure()
{
  int numberRows = solver_->getNumRows();
  int numberColumns = solver_->getNumCols();
  // look for DW
  int *blockStart = new int[3 * (numberRows + numberColumns) + 1];
  int *columnBlock = blockStart + numberRows;
  int *nextColumn = columnBlock + numberColumns;
  int *blockCount = nextColumn + numberColumns;
  int *blockEls = blockCount + numberRows + 1;
  int *countIntegers = blockEls + numberRows;
  memset(countIntegers, 0, numberColumns * sizeof(int));
  int direction[2] = { -1, 1 };
  int bestBreak = -1;
  double bestValue = 0.0;
  int iPass = 0;
  int halfway = (numberRows + 1) / 2;
  int firstMaster = -1;
  int lastMaster = -2;
  char dwPrint[200];
  // Column copy
  const CoinPackedMatrix *columnCopy = solver_->getMatrixByCol();
  //const double * element = columnCopy->getElements();
  const int *row = columnCopy->getIndices();
  const CoinBigIndex *columnStart = columnCopy->getVectorStarts();
  const int *columnLength = columnCopy->getVectorLengths();
  // Row copy
  const CoinPackedMatrix *rowCopy = solver_->getMatrixByRow();
  const int *column = rowCopy->getIndices();
  const int *rowLength = rowCopy->getVectorLengths();
  const CoinBigIndex *rowStart = rowCopy->getVectorStarts();
  //const double * elementByRow = rowCopy->getElements();
  while (iPass < 2) {
    int increment = direction[iPass];
    int start = increment > 0 ? 0 : numberRows - 1;
    int stop = increment > 0 ? numberRows : -1;
    int numberBlocks = 0;
    int thisBestBreak = -1;
    double thisBestValue = COIN_DBL_MAX;
    int numberRowsDone = 0;
    int numberMarkedColumns = 0;
    int maximumBlockSize = 0;
    for (int i = 0; i < numberRows + 2 * numberColumns; i++)
      blockStart[i] = -1;
    for (int i = 0; i < numberRows + 1; i++)
      blockCount[i] = 0;
    for (int iRow = start; iRow != stop; iRow += increment) {
      int iBlock = -1;
      for (CoinBigIndex j = rowStart[iRow]; j < rowStart[iRow] + rowLength[iRow]; j++) {
        int iColumn = column[j];
        int whichColumnBlock = columnBlock[iColumn];
        if (whichColumnBlock >= 0) {
          // column marked
          if (iBlock < 0) {
            // put row in that block
            iBlock = whichColumnBlock;
          } else if (iBlock != whichColumnBlock) {
            // merge
            blockCount[iBlock] += blockCount[whichColumnBlock];
            blockCount[whichColumnBlock] = 0;
            int jColumn = blockStart[whichColumnBlock];
            while (jColumn >= 0) {
              columnBlock[jColumn] = iBlock;
              iColumn = jColumn;
              jColumn = nextColumn[jColumn];
            }
            nextColumn[iColumn] = blockStart[iBlock];
            blockStart[iBlock] = blockStart[whichColumnBlock];
            blockStart[whichColumnBlock] = -1;
          }
        }
      }
      int n = numberMarkedColumns;
      if (iBlock < 0) {
        //new block
        if (rowLength[iRow]) {
          numberBlocks++;
          iBlock = numberBlocks;
          int jColumn = column[rowStart[iRow]];
          columnBlock[jColumn] = iBlock;
          blockStart[iBlock] = jColumn;
          numberMarkedColumns++;
          for (CoinBigIndex j = rowStart[iRow] + 1; j < rowStart[iRow] + rowLength[iRow]; j++) {
            int iColumn = column[j];
            columnBlock[iColumn] = iBlock;
            numberMarkedColumns++;
            nextColumn[jColumn] = iColumn;
            jColumn = iColumn;
          }
          blockCount[iBlock] = numberMarkedColumns - n;
        } else {
          // empty
          iBlock = numberRows;
        }
      } else {
        // put in existing block
        int jColumn = blockStart[iBlock];
        for (CoinBigIndex j = rowStart[iRow]; j < rowStart[iRow] + rowLength[iRow]; j++) {
          int iColumn = column[j];
          assert(columnBlock[iColumn] < 0 || columnBlock[iColumn] == iBlock);
          if (columnBlock[iColumn] < 0) {
            columnBlock[iColumn] = iBlock;
            numberMarkedColumns++;
            nextColumn[iColumn] = jColumn;
            jColumn = iColumn;
          }
        }
        blockStart[iBlock] = jColumn;
        blockCount[iBlock] += numberMarkedColumns - n;
      }
      maximumBlockSize = CoinMax(maximumBlockSize, blockCount[iBlock]);
      numberRowsDone++;
      if (thisBestValue * numberRowsDone > maximumBlockSize && numberRowsDone > halfway) {
        thisBestBreak = iRow;
        thisBestValue = static_cast< double >(maximumBlockSize) / static_cast< double >(numberRowsDone);
      }
    }
    if (thisBestBreak == stop)
      thisBestValue = COIN_DBL_MAX;
    iPass++;
    if (iPass == 1) {
      bestBreak = thisBestBreak;
      bestValue = thisBestValue;
    } else {
      if (bestValue < thisBestValue) {
        firstMaster = 0;
        lastMaster = bestBreak;
      } else {
        firstMaster = thisBestBreak; // ? +1
        lastMaster = numberRows;
      }
    }
  }
  numberBlocks_ = 0;
  if (firstMaster < lastMaster) {
    sprintf(dwPrint, "%d master rows %d <= < %d", lastMaster - firstMaster,
      firstMaster, lastMaster);
    model_->messageHandler()->message(CBC_FPUMP1, model_->messages())
      << dwPrint
      << CoinMessageEol;
    for (int i = 0; i < numberRows + 2 * numberColumns; i++)
      blockStart[i] = -1;
    for (int i = firstMaster; i < lastMaster; i++)
      blockStart[i] = -2;
    int firstRow = 0;
    int numberBlocks = -1;
    while (true) {
      for (; firstRow < numberRows; firstRow++) {
        if (blockStart[firstRow] == -1)
          break;
      }
      if (firstRow == numberRows)
        break;
      int nRows = 0;
      numberBlocks++;
      int numberStack = 1;
      blockCount[0] = firstRow;
      while (numberStack) {
        int iRow = blockCount[--numberStack];
        for (CoinBigIndex j = rowStart[iRow]; j < rowStart[iRow] + rowLength[iRow]; j++) {
          int iColumn = column[j];
          int iBlock = columnBlock[iColumn];
          if (iBlock < 0) {
            columnBlock[iColumn] = numberBlocks;
            for (CoinBigIndex k = columnStart[iColumn];
                 k < columnStart[iColumn] + columnLength[iColumn]; k++) {
              int jRow = row[k];
              int rowBlock = blockStart[jRow];
              if (rowBlock == -1) {
                nRows++;
                blockStart[jRow] = numberBlocks;
                blockCount[numberStack++] = jRow;
              }
            }
          }
        }
      }
      if (!nRows) {
        // empty!!
        numberBlocks--;
      }
      firstRow++;
    }
    // adjust
    numberBlocks++;
    for (int i = 0; i < numberBlocks; i++) {
      blockCount[i] = 0;
      nextColumn[i] = 0;
    }
    int numberEmpty = 0;
    int numberMaster = 0;
    memset(blockEls, 0, numberBlocks * sizeof(int));
    for (int iRow = 0; iRow < numberRows; iRow++) {
      int iBlock = blockStart[iRow];
      if (iBlock >= 0) {
        blockCount[iBlock]++;
        blockEls[iBlock] += rowLength[iRow];
      } else {
        if (iBlock == -2)
          numberMaster++;
        else
          numberEmpty++;
      }
    }
    int numberEmptyColumns = 0;
    int numberMasterColumns = 0;
    int numberMasterIntegers = 0;
    for (int iColumn = 0; iColumn < numberColumns; iColumn++) {
      int iBlock = columnBlock[iColumn];
      bool isInteger = (solver_->isInteger(iColumn));
      if (iBlock >= 0) {
        nextColumn[iBlock]++;
        if (isInteger)
          countIntegers[iBlock]++;
      } else {
        if (isInteger)
          numberMasterIntegers++;
        if (columnLength[iColumn])
          numberMasterColumns++;
        else
          numberEmptyColumns++;
      }
    }
    int largestRows = 0;
    int largestColumns = 0;
    for (int i = 0; i < numberBlocks; i++) {
      if (blockCount[i] + nextColumn[i] > largestRows + largestColumns) {
        largestRows = blockCount[i];
        largestColumns = nextColumn[i];
      }
    }
    bool useful = true;
    if (numberMaster > halfway || largestRows * 3 > numberRows)
      useful = false;
    sprintf(dwPrint, "%s %d blocks (largest %d,%d), %d master rows (%d empty) out of %d, %d master columns (%d empty, %d integer) out of %d",
      useful ? "**Useful" : "NoGood",
      numberBlocks, largestRows, largestColumns, numberMaster, numberEmpty, numberRows,
      numberMasterColumns, numberEmptyColumns, numberMasterIntegers,
      numberColumns);
    model_->messageHandler()->message(CBC_FPUMP1, model_->messages())
      << dwPrint
      << CoinMessageEol;
    // columnBlock is columnBlock and blockStart is rowBlock
    // See if we want to compress
    if (!keepContinuous_) {
      // use blockEls
      int newNumber = 0;
      for (int i = 0; i < numberBlocks; i++) {
        sprintf(dwPrint, "Block %d has %d rows and %d columns (%d elements, %d integers)",
          i, blockCount[i], nextColumn[i], blockEls[i], countIntegers[i]);
        model_->messageHandler()->message(CBC_FPUMP2, model_->messages())
          << dwPrint
          << CoinMessageEol;
        if (countIntegers[i]) {
          blockEls[i] = newNumber;
          newNumber++;
        } else {
          blockEls[i] = -1;
        }
      }
      for (int i = 0; i < numberRows; i++) {
        int iBlock = blockStart[i];
        if (iBlock >= 0)
          blockStart[i] = blockEls[iBlock];
      }
      for (int i = 0; i < numberColumns; i++) {
        int iBlock = columnBlock[i];
        if (iBlock >= 0)
          columnBlock[i] = blockEls[iBlock];
      }
      if (newNumber < numberBlocks) {
        sprintf(dwPrint, "Number of blocks reduced from %d to %d",
          numberBlocks, newNumber);
        model_->messageHandler()->message(CBC_FPUMP1, model_->messages())
          << dwPrint
          << CoinMessageEol;
      }
      numberBlocks = newNumber;
    }
    // now set up structures
    numberBlocks_ = numberBlocks;
    // so callBack can modify
    whichRowBlock_ = blockStart;
    whichColumnBlock_ = columnBlock;
    intArray_ = NULL;
    doubleArray_ = NULL;
    (*(functionPointer_))(this, NULL, 0);

    saveLower_ = CoinCopyOfArray(solver_->getColLower(), numberColumns);
    saveUpper_ = CoinCopyOfArray(solver_->getColUpper(), numberColumns);
    startRowBlock_ = new int[numberBlocks_ + 2];
    backwardRow_ = new int[numberRows];
    rowsInBlock_ = new int[numberRows];
    whichRowBlock_ = new int[numberRows];
    startColumnBlock_ = new int[numberBlocks_ + 2];
    columnsInBlock_ = new int[numberColumns];
    whichColumnBlock_ = new int[numberColumns];
    intsInBlock_ = new int[numberBlocks_];
    // use for counts
    memset(rowsInBlock_, 0, numberBlocks_ * sizeof(int));
    numberMasterRows_ = 0;
    for (int i = 0; i < numberRows; i++) {
      int iBlock = blockStart[i];
      if (iBlock >= 0) {
        rowsInBlock_[iBlock]++;
        whichRowBlock_[i] = iBlock;
        backwardRow_[i] = -1;
      } else {
        whichRowBlock_[i] = -1;
        backwardRow_[i] = numberMasterRows_;
        numberMasterRows_++;
      }
    }
    memset(columnsInBlock_, 0, numberBlocks_ * sizeof(int));
    memset(intsInBlock_, 0, numberBlocks_ * sizeof(int));
    numberMasterColumns_ = 0;
    for (int i = 0; i < numberColumns; i++) {
      int iBlock = columnBlock[i];
      if (iBlock >= 0) {
        columnsInBlock_[iBlock]++;
        whichColumnBlock_[i] = iBlock;
        if (solver_->isInteger(i))
          intsInBlock_[iBlock]++;
      } else {
        whichColumnBlock_[i] = -1;
        numberMasterColumns_++;
      }
    }
    // starts
    int nRow = 0;
    int nColumn = 0;
    int maxIntsInBlock = 0;
    for (int i = 0; i < numberBlocks_; i++) {
      maxIntsInBlock = CoinMax(maxIntsInBlock, intsInBlock_[i]);
      startRowBlock_[i] = nRow;
      startColumnBlock_[i] = nColumn;
      nRow += rowsInBlock_[i];
      nColumn += columnsInBlock_[i];
    }
    // may not be used - but set anyway
    sizeFingerPrint_ = (maxIntsInBlock + 31) / 32;
    startRowBlock_[numberBlocks_] = nRow;
    startColumnBlock_[numberBlocks_] = nColumn;
    startRowBlock_[numberBlocks_ + 1] = numberRows;
    startColumnBlock_[numberBlocks_ + 1] = numberColumns;
    // do lists
    for (int i = 0; i < numberRows; ++i) {
      int iBlock = whichRowBlock_[i];
      if (iBlock < 0)
        iBlock = numberBlocks_;
      int k = startRowBlock_[iBlock];
      startRowBlock_[iBlock] = k + 1;
      rowsInBlock_[k] = i;
    }
    for (int i = numberBlocks + 1; i > 0; i--)
      startRowBlock_[i] = startRowBlock_[i - 1];
    startRowBlock_[0] = 0;
    for (int i = 0; i < numberColumns; ++i) {
      int iBlock = whichColumnBlock_[i];
      if (iBlock < 0)
        iBlock = numberBlocks_;
      int k = startColumnBlock_[iBlock];
      startColumnBlock_[iBlock] = k + 1;
      columnsInBlock_[k] = i;
    }
    for (int i = numberBlocks + 1; i > 0; i--)
      startColumnBlock_[i] = startColumnBlock_[i - 1];
    startColumnBlock_[0] = 0;
    if (numberBlocks_ < 10000) {
      affinity_ = new unsigned short[numberBlocks_ * numberBlocks_];
      // compute space needed
      int *build = new int[numberMasterRows_];
      memset(build, 0, numberMasterRows_ * sizeof(int));
      int nSpace = 0;
      for (int iBlock = 0; iBlock < numberBlocks_; iBlock++) {
        int start = startColumnBlock_[iBlock];
        int end = startColumnBlock_[iBlock + 1];
        for (int i = start; i < end; i++) {
          int iColumn = columnsInBlock_[i];
          for (CoinBigIndex j = columnStart[iColumn];
               j < columnStart[iColumn] + columnLength[iColumn]; j++) {
            int iRow = row[j];
            iRow = backwardRow_[iRow];
            if (iRow >= 0)
              build[iRow]++;
          }
        }
        for (int i = 0; i < numberMasterRows_; i++) {
          int value = build[i];
          if (value) {
            build[i] = 0;
            nSpace++;
          }
        }
      }
      // get arrays
      int *starts = new int[numberBlocks_ + 1 + 2 * nSpace];
      memset(affinity_, 0, numberBlocks_ * numberBlocks_ * sizeof(unsigned short));
      // fill arrays
      int *rowM = starts + numberBlocks_ + 1;
      int *sumM = rowM + nSpace;
      nSpace = 0;
      starts[0] = 0;
      for (int iBlock = 0; iBlock < numberBlocks_; iBlock++) {
        int start = startColumnBlock_[iBlock];
        int end = startColumnBlock_[iBlock + 1];
        for (int i = start; i < end; i++) {
          int iColumn = columnsInBlock_[i];
          for (CoinBigIndex j = columnStart[iColumn];
               j < columnStart[iColumn] + columnLength[iColumn]; j++) {
            int iRow = row[j];
            iRow = backwardRow_[iRow];
            if (iRow >= 0)
              build[iRow]++;
          }
        }
        for (int i = 0; i < numberMasterRows_; i++) {
          int value = build[i];
          if (value) {
            build[i] = 0;
            sumM[nSpace] = value;
            rowM[nSpace++] = i;
          }
        }
        starts[iBlock + 1] = nSpace;
      }
      for (int iBlock = 0; iBlock < numberBlocks_; iBlock++) {
        int startI = starts[iBlock];
        int endI = starts[iBlock + 1];
        if (endI == startI)
          continue;
        for (int jBlock = iBlock + 1; jBlock < numberBlocks_; jBlock++) {
          int startJ = starts[jBlock];
          int endJ = starts[jBlock + 1];
          if (endJ == startJ)
            continue;
          double sum = 0.0;
          int i = startI;
          int j = startJ;
          int rowI = rowM[i];
          int rowJ = rowM[j];
          while (rowI != COIN_INT_MAX && rowJ != COIN_INT_MAX) {
            if (rowI < rowJ) {
              i++;
              if (i < endI)
                rowI = rowM[i];
              else
                rowI = COIN_INT_MAX;
            } else if (rowI > rowJ) {
              j++;
              if (j < endJ)
                rowJ = rowM[j];
              else
                rowJ = COIN_INT_MAX;
            } else {
              // bias ????????
              sum += sumM[i] * sumM[j];
              i++;
              if (i < endI)
                rowI = rowM[i];
              else
                rowI = COIN_INT_MAX;
              j++;
              if (j < endJ)
                rowJ = rowM[j];
              else
                rowJ = COIN_INT_MAX;
            }
          }
          if (sum > 65535)
            sum = 65535;
          unsigned short value = static_cast< unsigned short >(sum);
          affinity_[iBlock * numberBlocks + jBlock] = value;
          affinity_[jBlock * numberBlocks + iBlock] = value;
        }
      }
      // statistics
      int nTotalZero = 0;
      int base = 0;
      for (int iBlock = 0; iBlock < numberBlocks_; iBlock++) {
        //int aff = 0;
        int nZero = 0;
        for (int jBlock = 0; jBlock < numberBlocks_; jBlock++) {
          if (iBlock != jBlock) {
            if (affinity_[base + jBlock])
              ;//aff += affinity_[base + jBlock];
            else
              nZero++;
          }
        }
        //printf("Block %d has affinity %d but zero with %d blocks",
        //     iBlock,aff,nZero);
        nTotalZero += nZero;
        base += numberBlocks;
      }
      sprintf(dwPrint, "Total not affinity %d - average %g%%",
        nTotalZero, 100.0 * (static_cast< double >(nTotalZero) / (numberBlocks * numberBlocks)));
      model_->messageHandler()->message(CBC_FPUMP1, model_->messages())
        << dwPrint
        << CoinMessageEol;

      delete[] starts;
      delete[] build;
    } else {
      sprintf(dwPrint, "Too many blocks - no affinity");
      model_->messageHandler()->message(CBC_FPUMP1, model_->messages())
        << dwPrint
        << CoinMessageEol;
    }
    if (fullDWEverySoOften_ > 0) {
      setupDWStructures();
    }
  }
  delete[] blockStart;
}
// Add DW proposals
int CbcHeuristicDW::addDW(const double *solution, int numberBlocksUsed,
  const int *whichBlocks)
{
  char dwPrint[200];
  if (numberDW_ + numberBlocksUsed > maximumDW_) {
    // extend
    int n = maximumDW_ + 5 * numberBlocks_;
    double *weightsT = new double[n];
    int *dwBlockT = new int[n];
    unsigned int *fingerT = new unsigned int[n * sizeFingerPrint_];
    memcpy(weightsT, weights_, numberDW_ * sizeof(double));
    memcpy(dwBlockT, dwBlock_, numberDW_ * sizeof(int));
    memcpy(fingerT, fingerPrint_,
      numberDW_ * sizeFingerPrint_ * sizeof(unsigned int));
    delete[] weights_;
    weights_ = weightsT;
    delete[] dwBlock_;
    dwBlock_ = dwBlockT;
    delete[] fingerPrint_;
    fingerPrint_ = fingerT;
    maximumDW_ = n;
  }
  //int numberColumns = solver_->getNumCols();
  //int numberRows = solver_->getNumRows();
  // get space to add elements
#define MAX_ADD 100000
  CoinBigIndex *startsDW = new CoinBigIndex[numberBlocks_ + 1 + MAX_ADD];
  int *rowDW = reinterpret_cast< int * >(startsDW + numberBlocks_ + 1);
  double *elementDW = new double[MAX_ADD + 3 * numberBlocks_ + numberMasterRows_];
  double *newCost = elementDW + MAX_ADD;
  double *newLower = newCost + numberBlocks_;
  double *newUpper = newLower + numberBlocks_;
  double *build = newUpper + numberBlocks_;
  memset(build, 0, numberMasterRows_ * sizeof(double));
  int nAdd = 0;
  int nTotalAdded = 0;
  int nEls = 0;
  startsDW[0] = 0;
  // Column copy
  const double *element = solver_->getMatrixByCol()->getElements();
  const int *row = solver_->getMatrixByCol()->getIndices();
  const CoinBigIndex *columnStart = solver_->getMatrixByCol()->getVectorStarts();
  const int *columnLength = solver_->getMatrixByCol()->getVectorLengths();
  const double *objective = solver_->getObjCoefficients();
  for (int jBlock = 0; jBlock < numberBlocksUsed; jBlock++) {
    double thisWeight = 0.0;
    double thisWeightC = 0.0;
    double thisCost = 0.0;
    int nElInMaster = 0;
    int nElIntInMaster = 0;
    int nElIntInMaster1 = 0;
    int iBlock = whichBlocks[jBlock];
    int start = startColumnBlock_[iBlock];
    int end = startColumnBlock_[iBlock + 1];
    unsigned int *finger = fingerPrint_ + sizeFingerPrint_ * (nAdd + numberDW_);
    memset(finger, 0, sizeFingerPrint_ * sizeof(unsigned int));
    int iBit = 0;
    for (int i = start; i < end; i++) {
      int iColumn = columnsInBlock_[i];
      bool isInteger = solver_->isInteger(iColumn);
      double value = solution[iColumn];
      if (isInteger) {
        if (value > 1.0e-6)
          finger[0] |= 1 << iBit;
        iBit++;
        if (iBit == 32) {
          finger++;
          iBit++;
        }
      }
      thisCost += value * objective[iColumn];
      for (CoinBigIndex j = columnStart[iColumn];
           j < columnStart[iColumn] + columnLength[iColumn]; j++) {
        int iRow = row[j];
        iRow = backwardRow_[iRow];
        if (iRow >= 0) {
          nElInMaster++;
          double elementValue = element[j];
          build[iRow] += value * elementValue;
          if (isInteger) {
            nElIntInMaster++;
            if (value)
              nElIntInMaster1++;
            thisWeight += random_[iRow] * value * elementValue;
          } else {
            value = 0.0001 * floor(value * 10000.0 + 0.5);
          }
          thisWeightC += random_[iRow] * value * elementValue;
        }
      }
    }
    // see if already in
    sprintf(dwPrint, "block %d nel %d nelInt %d nelInt1 %d - weight %g (%g)",
      iBlock, nElInMaster, nElIntInMaster, nElIntInMaster1,
      thisWeight, thisWeightC);
    model_->messageHandler()->message(CBC_FPUMP2, model_->messages())
      << dwPrint
      << CoinMessageEol;
    int iProposal;
    for (iProposal = 0; iProposal < numberDW_; iProposal++) {
      if (iBlock == dwBlock_[iProposal] && weights_[iProposal] == thisWeightC)
        break;
    }
    if (iProposal < numberDW_) {
      sprintf(dwPrint, "above looks like duplicate");
      model_->messageHandler()->message(CBC_FPUMP2, model_->messages())
        << dwPrint
        << CoinMessageEol;
      memset(build, 0, numberMasterRows_ * sizeof(double));
      //iProposal=numberDW;
    }
    if (iProposal == numberDW_) {
      // new
      // could build in stages
      assert(nEls + numberMasterRows_ < MAX_ADD);
      for (int i = 0; i < numberMasterRows_; i++) {
        double value = build[i];
        if (value) {
          build[i] = 0.0;
          if (fabs(value) > 1.0e-10) {
            elementDW[nEls] = value;
            rowDW[nEls++] = i;
          }
        }
      }
      // convexity
      elementDW[nEls] = 1.0;
      rowDW[nEls++] = numberMasterRows_ + iBlock;
      weights_[nAdd + numberDW_] = thisWeightC;
      dwBlock_[nAdd + numberDW_] = iBlock;
      newLower[nAdd] = 0.0;
      newUpper[nAdd] = 1.0;
      newCost[nAdd++] = thisCost;
      startsDW[nAdd] = nEls;
    }
    if (nEls + numberMasterRows_ > MAX_ADD) {
      sprintf(dwPrint, "Adding %d proposals with %d elements - out of room",
        nAdd, nEls);
      model_->messageHandler()->message(CBC_FPUMP1, model_->messages())
        << dwPrint
        << CoinMessageEol;
      dwSolver_->addCols(nAdd, startsDW, rowDW, elementDW, newLower,
        newUpper, newCost);
      numberDW_ += nAdd;
      nTotalAdded += nAdd;
      nAdd = 0;
      nEls = 0;
    }
  }
  if (nAdd) {
    sprintf(dwPrint, "Adding %d proposals with %d elements",
      nAdd, nEls);
    model_->messageHandler()->message(CBC_FPUMP1, model_->messages())
      << dwPrint
      << CoinMessageEol;
    dwSolver_->addCols(nAdd, startsDW, rowDW, elementDW, newLower,
      newUpper, newCost);
    nTotalAdded += nAdd;
    numberDW_ += nAdd;
  }
  delete[] startsDW;
  delete[] elementDW;
  if (nTotalAdded) {
    double *objs = new double[numberDWTimes_ + 1];
    memcpy(objs, objectiveDW_, numberDWTimes_ * sizeof(double));
    delete[] objectiveDW_;
    objectiveDW_ = objs;
    int *temp = new int[numberDWTimes_ + 1];
    memcpy(temp, numberColumnsDW_, numberDWTimes_ * sizeof(int));
    delete[] numberColumnsDW_;
    numberColumnsDW_ = temp;
    numberColumnsDW_[numberDWTimes_] = dwSolver_->getNumCols();
    objectiveDW_[numberDWTimes_++] = objectiveValue(solution);
  }
  return nTotalAdded;
}
// Pass in a solution
void CbcHeuristicDW::passInSolution(const double *solution)
{
  if (fullDWEverySoOften_ > 0) {
    int *which = new int[numberBlocks_];
    for (int i = 0; i < numberBlocks_; i++)
      which[i] = i;
    addDW(solution, numberBlocks_, which);
    delete[] which;
  }
  if (objectiveValue(solution) < bestObjective_ - 1.0e-5) {
    bestObjective_ = objectiveValue(solution);
    int numberColumns = solver_->getNumCols();
    if (!bestSolution_)
      bestSolution_ = new double[numberColumns];
    memcpy(bestSolution_, solution, numberColumns * sizeof(double));
  }
}
// Pass in continuous solution
void CbcHeuristicDW::passInContinuousSolution(const double *solution)
{
  int numberColumns = solver_->getNumCols();
  if (!continuousSolution_)
    continuousSolution_ = new double[numberColumns];
  memcpy(continuousSolution_, solution, numberColumns * sizeof(double));
}
// Objective value (could also check validity)
double
CbcHeuristicDW::objectiveValue(const double *solution)
{
  // compute objective value
  double objOffset = 0.0;
  solver_->getDblParam(OsiObjOffset, objOffset);
  double objectiveValue = -objOffset;
  int numberColumns = solver_->getNumCols();
  const double *objective = solver_->getObjCoefficients();
  int logLevel = model_->messageHandler()->logLevel();
  for (int i = 0; i < numberColumns; i++) {
    double value = solution[i];
    if (logLevel > 1) {
      if (solver_->isInteger(i)) {
        if (fabs(value - floor(value + 0.5)) > 1.0e-7)
          printf("Bad integer value for %d of %g\n", i, value);
      }
    }
    objectiveValue += objective[i] * value;
  }
  return objectiveValue;
}
// Objective value when whichDw created
double
CbcHeuristicDW::objectiveValueWhen(int whichDW) const
{
  if (whichDW >= numberDWTimes_)
    return COIN_DBL_MAX;
  else
    return objectiveDW_[whichDW];
}
// Number of columns in DW
int CbcHeuristicDW::numberColumnsDW(int whichDW) const
{
  if (whichDW >= numberDWTimes_)
    return COIN_INT_MAX;
  else
    return numberColumnsDW_[whichDW];
}
// DW model (user must delete)
OsiSolverInterface *
CbcHeuristicDW::DWModel(int whichDW) const
{
  if (whichDW >= numberDWTimes_)
    return NULL;
  OsiSolverInterface *newSolver = dwSolver_->clone();
  int numberColumns2 = newSolver->getNumCols();
  int numberColumns = numberColumnsDW_[whichDW];
  if (numberColumns < numberColumns2) {
    int *del = new int[numberColumns2 - numberColumns];
    for (int i = numberColumns; i < numberColumns2; i++)
      del[i - numberColumns] = i;
    newSolver->deleteCols(numberColumns2 - numberColumns, del);
    delete[] del;
  }
  // Set all to integer that need setting
  for (int i = numberMasterColumns_; i < numberColumns; i++) {
    newSolver->setContinuous(i);
  }
  int numberDW = numberColumns - numberMasterColumns_;
  for (int iBlock = 0; iBlock < numberBlocks_; iBlock++) {
    bool allSame = true;
    unsigned int *finger = fingerPrint_;
    unsigned int *fingerTest = NULL;
    for (int i = 0; i < numberDW; i++) {
      if (dwBlock_[i] == iBlock) {
        if (fingerTest) {
          for (int j = 0; j < sizeFingerPrint_; j++) {
            if (finger[j] != fingerTest[j]) {
              allSame = false;
              break;
            }
          }
          if (!allSame)
            break;
        } else {
          fingerTest = finger;
        }
      }
      finger += sizeFingerPrint_;
    }
    if (!allSame) {
      // Set all to integer that need setting
      for (int i = 0; i < numberDW; i++) {
        if (iBlock == dwBlock_[i]) {
          int iColumn = numberMasterColumns_ + i;
          newSolver->setInteger(iColumn);
        }
      }
    }
  }
  //newSolver->writeMps("dw","mps");
  return newSolver;
}
/* DW Proposal actions
   fullDWEverySoOften -
   0 - off
   k - every k times solution gets better
*/
void CbcHeuristicDW::setProposalActions(int fullDWEverySoOften)
{
  fullDWEverySoOften_ = fullDWEverySoOften;
  if (fullDWEverySoOften_ > 0 && !random_)
    setupDWStructures();
}
// Set up DW structure
void CbcHeuristicDW::setupDWStructures()
{
  char dwPrint[200];
  random_ = new double[numberMasterRows_];
  for (int i = 0; i < numberMasterRows_; i++)
    random_[i] = CoinDrand48();
  weights_ = new double[numberBlocks_];
  dwBlock_ = new int[numberBlocks_];
  fingerPrint_ = new unsigned int[numberBlocks_ * sizeFingerPrint_];
  // create dwSolver
  int numberColumns = solver_->getNumCols();
  int numberRows = solver_->getNumRows();
  int *tempRow = new int[numberRows + numberColumns];
  int *tempColumn = tempRow + numberRows;
  int numberMasterRows = 0;
  for (int i = 0; i < numberRows; i++) {
    int iBlock = whichRowBlock_[i];
    if (iBlock < 0)
      tempRow[numberMasterRows++] = i;
  }
  int numberMasterColumns = 0;
  for (int i = 0; i < numberColumns; i++) {
    int iBlock = whichColumnBlock_[i];
    if (iBlock < 0)
      tempColumn[numberMasterColumns++] = i;
  }
  OsiClpSolverInterface *solver = dynamic_cast< OsiClpSolverInterface * >(solver_);
  ClpSimplex *tempModel = new ClpSimplex(solver->getModelPtr(),
    numberMasterRows, tempRow,
    numberMasterColumns, tempColumn);
  // add convexity constraints
  double *rhs = new double[numberBlocks_];
  for (int i = 0; i < numberBlocks_; i++)
    rhs[i] = 1.0;
  tempModel->addRows(numberBlocks_, rhs, rhs, NULL, NULL, NULL);
  delete[] rhs;
  OsiClpSolverInterface *clpSolver = new OsiClpSolverInterface(tempModel, true);
  clpSolver->getModelPtr()->setDualObjectiveLimit(COIN_DBL_MAX);
  dwSolver_ = clpSolver;
  sprintf(dwPrint, "DW model has %d master rows, %d master columns and %d convexity rows",
    numberMasterRows, numberMasterColumns, numberBlocks_);
  model_->messageHandler()->message(CBC_FPUMP1, model_->messages())
    << dwPrint
    << CoinMessageEol;
  // do master integers
  for (int i = 0; i < numberMasterColumns; i++) {
    int iColumn = tempColumn[i];
    if (solver->isInteger(iColumn))
      dwSolver_->setInteger(i);
  }
  delete[] tempRow;
}

/* vi: softtabstop=2 shiftwidth=2 expandtab tabstop=2
*/
