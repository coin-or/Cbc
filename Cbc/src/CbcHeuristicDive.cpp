/* $Id$ */
// Copyright (C) 2008, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

#if defined(_MSC_VER)
// Turn off compiler warning about long names
#pragma warning(disable : 4786)
#endif

#include "CbcStrategy.hpp"
#include "CbcModel.hpp"
#include "CbcSubProblem.hpp"
#include "CbcSimpleInteger.hpp"
#include "OsiAuxInfo.hpp"
#include "CoinTime.hpp"

#ifdef COIN_HAS_CLP
#include "OsiClpSolverInterface.hpp"
#endif
#include "CbcHeuristicDive.hpp"

//#define DIVE_FIX_BINARY_VARIABLES
//#define DIVE_DEBUG
#ifdef DIVE_DEBUG
#define DIVE_PRINT = 2
#endif

// Default Constructor
CbcHeuristicDive::CbcHeuristicDive()
  : CbcHeuristic()
{
  // matrix and row copy will automatically be empty
  downLocks_ = NULL;
  upLocks_ = NULL;
  downArray_ = NULL;
  upArray_ = NULL;
  priority_ = NULL;
  percentageToFix_ = 0.2;
  maxIterations_ = 100;
  maxSimplexIterations_ = 10000;
  maxSimplexIterationsAtRoot_ = 1000000;
  maxTime_ = 600;
  whereFrom_ = 255 - 2 - 16 + 256;
  decayFactor_ = 1.0;
  smallObjective_ = 1.0e-10;
}

// Constructor from model
CbcHeuristicDive::CbcHeuristicDive(CbcModel &model)
  : CbcHeuristic(model)
{
  downLocks_ = NULL;
  upLocks_ = NULL;
  downArray_ = NULL;
  upArray_ = NULL;
  priority_ = NULL;
  // Get a copy of original matrix
  assert(model.solver());
  // model may have empty matrix - wait until setModel
  const CoinPackedMatrix *matrix = model.solver()->getMatrixByCol();
  if (matrix) {
    matrix_ = *matrix;
    matrixByRow_ = *model.solver()->getMatrixByRow();
    validate();
  }
  percentageToFix_ = 0.2;
  maxTime_ = 600;
  smallObjective_ = 1.0e-10;
  maxIterations_ = 100;
  maxSimplexIterations_ = 10000;
  maxSimplexIterationsAtRoot_ = 1000000;
  whereFrom_ = 255 - 2 - 16 + 256;
  decayFactor_ = 1.0;
  smallObjective_ = 1.0e-10;
}

// Destructor
CbcHeuristicDive::~CbcHeuristicDive()
{
  delete[] downLocks_;
  delete[] upLocks_;
  delete[] priority_;
  assert(!downArray_);
}

// Create C++ lines to get to current state
void CbcHeuristicDive::generateCpp(FILE *fp, const char *heuristic)
{
  // hard coded as CbcHeuristic virtual
  CbcHeuristic::generateCpp(fp, heuristic);
  if (percentageToFix_ != 0.2)
    fprintf(fp, "3  %s.setPercentageToFix(%.f);\n", heuristic, percentageToFix_);
  else
    fprintf(fp, "4  %s.setPercentageToFix(%.f);\n", heuristic, percentageToFix_);
  if (maxIterations_ != 100)
    fprintf(fp, "3  %s.setMaxIterations(%d);\n", heuristic, maxIterations_);
  else
    fprintf(fp, "4  %s.setMaxIterations(%d);\n", heuristic, maxIterations_);
  if (maxSimplexIterations_ != 10000)
    fprintf(fp, "3  %s.setMaxSimplexIterations(%d);\n", heuristic, maxSimplexIterations_);
  else
    fprintf(fp, "4  %s.setMaxSimplexIterations(%d);\n", heuristic, maxSimplexIterations_);
  if (maxTime_ != 600)
    fprintf(fp, "3  %s.setMaxTime(%.2f);\n", heuristic, maxTime_);
  else
    fprintf(fp, "4  %s.setMaxTime(%.2f);\n", heuristic, maxTime_);
}

// Copy constructor
CbcHeuristicDive::CbcHeuristicDive(const CbcHeuristicDive &rhs)
  : CbcHeuristic(rhs)
  , matrix_(rhs.matrix_)
  , matrixByRow_(rhs.matrixByRow_)
  , percentageToFix_(rhs.percentageToFix_)
  , maxTime_(rhs.maxTime_)
  , smallObjective_(rhs.smallObjective_)
  , maxIterations_(rhs.maxIterations_)
  , maxSimplexIterations_(rhs.maxSimplexIterations_)
  , maxSimplexIterationsAtRoot_(rhs.maxSimplexIterationsAtRoot_)
{
  downArray_ = NULL;
  upArray_ = NULL;
  if (rhs.downLocks_) {
    int numberIntegers = model_->numberIntegers();
    downLocks_ = CoinCopyOfArray(rhs.downLocks_, numberIntegers);
    upLocks_ = CoinCopyOfArray(rhs.upLocks_, numberIntegers);
    priority_ = CoinCopyOfArray(rhs.priority_, numberIntegers);
  } else {
    downLocks_ = NULL;
    upLocks_ = NULL;
    priority_ = NULL;
  }
}

// Assignment operator
CbcHeuristicDive &
CbcHeuristicDive::operator=(const CbcHeuristicDive &rhs)
{
  if (this != &rhs) {
    CbcHeuristic::operator=(rhs);
    matrix_ = rhs.matrix_;
    matrixByRow_ = rhs.matrixByRow_;
    percentageToFix_ = rhs.percentageToFix_;
    maxIterations_ = rhs.maxIterations_;
    maxSimplexIterations_ = rhs.maxSimplexIterations_;
    maxSimplexIterationsAtRoot_ = rhs.maxSimplexIterationsAtRoot_;
    maxTime_ = rhs.maxTime_;
    smallObjective_ = rhs.smallObjective_;
    delete[] downLocks_;
    delete[] upLocks_;
    delete[] priority_;
    if (rhs.downLocks_) {
      int numberIntegers = model_->numberIntegers();
      downLocks_ = CoinCopyOfArray(rhs.downLocks_, numberIntegers);
      upLocks_ = CoinCopyOfArray(rhs.upLocks_, numberIntegers);
      priority_ = CoinCopyOfArray(rhs.priority_, numberIntegers);
    } else {
      downLocks_ = NULL;
      upLocks_ = NULL;
      priority_ = NULL;
    }
  }
  return *this;
}

// Resets stuff if model changes
void CbcHeuristicDive::resetModel(CbcModel *model)
{
  model_ = model;
  assert(model_->solver());
  // Get a copy of original matrix
  const CoinPackedMatrix *matrix = model_->solver()->getMatrixByCol();
  // model may have empty matrix - wait until setModel
  if (matrix) {
    matrix_ = *matrix;
    matrixByRow_ = *model->solver()->getMatrixByRow();
    validate();
  }
  setPriorities();
}

// update model
void CbcHeuristicDive::setModel(CbcModel *model)
{
  model_ = model;
  assert(model_->solver());
  // Get a copy of original matrix
  const CoinPackedMatrix *matrix = model_->solver()->getMatrixByCol();
  if (matrix) {
    matrix_ = *matrix;
    matrixByRow_ = *model->solver()->getMatrixByRow();
    // make sure model okay for heuristic
    validate();
  }
  setPriorities();
}

// Sets priorities if any
void CbcHeuristicDive::setPriorities()
{
  delete[] priority_;
  assert(model_);
  priority_ = NULL;
  if (!model_->objects())
    return;
  bool gotPriorities = false;
  int numberIntegers = model_->numberIntegers();
  int priority1 = -COIN_INT_MAX;
  int priority2 = COIN_INT_MAX;
  smallObjective_ = 0.0;
  const double *objective = model_->solver()->getObjCoefficients();
  int numberObjects = model_->numberObjects();
  for (int i = 0; i < numberObjects; i++) {
    OsiObject *object = model_->modifiableObject(i);
    const CbcSimpleInteger *thisOne = dynamic_cast< const CbcSimpleInteger * >(object);
    if (!thisOne)
      continue; // Not integer
    int iColumn = thisOne->columnNumber();
    smallObjective_ += objective[iColumn];
    int level = thisOne->priority();
    priority1 = CoinMax(priority1, level);
    priority2 = CoinMin(priority2, level);
    if (thisOne->preferredWay() != 0)
      gotPriorities = true;
  }
  smallObjective_ = CoinMax(1.0e-10, 1.0e-5 * (smallObjective_ / numberIntegers));
  if (gotPriorities || priority1 > priority2) {
    priority_ = new PriorityType[numberIntegers];
    int nInteger = 0;
    for (int i = 0; i < numberObjects; i++) {
      OsiObject *object = model_->modifiableObject(i);
      const CbcSimpleInteger *thisOne = dynamic_cast< const CbcSimpleInteger * >(object);
      if (!thisOne)
        continue; // Not integer
      int level = thisOne->priority() - priority2;
      assert(level < (1 << 29));
      assert(nInteger < numberIntegers);
      priority_[nInteger].priority = static_cast< unsigned int >(level);
      int direction = 0;
      if (thisOne->preferredWay() < 0)
        direction = 1;
      else if (thisOne->preferredWay() > 0)
        direction = 1 | 1;
      // at present don't try other way is not used
      priority_[nInteger++].direction = static_cast< unsigned char >(direction);
    }
    assert(nInteger == numberIntegers);
  }
}

bool CbcHeuristicDive::canHeuristicRun()
{
  if (model_->bestSolution() || model_->getNodeCount()) {
    if (when_ == 3 || (when_ == 4 && numberSolutionsFound_))
      return false;
  }
  return shouldHeurRun_randomChoice();
}

inline bool compareBinaryVars(const PseudoReducedCost obj1,
  const PseudoReducedCost obj2)
{
  return obj1.pseudoRedCost > obj2.pseudoRedCost;
}

// inner part of dive
int CbcHeuristicDive::solution(double &solutionValue, int &numberNodes,
  int &numberCuts, OsiRowCut **cuts,
  CbcSubProblem **&nodes,
  double *newSolution)
{
#if DIVE_PRINT
  int nRoundInfeasible = 0;
  int nRoundFeasible = 0;
  printf("Entering %s - fix %.1f%% maxTime %.2f maxPasses %d - max iterations %d (at root %d) - when to do %d\n",
    heuristicName_.c_str(), percentageToFix_ * 100.0, maxTime_, maxIterations_,
    maxSimplexIterations_, maxSimplexIterationsAtRoot_, when());
#endif
  int reasonToStop = 0;
  double time1 = CoinCpuTime();
  int numberSimplexIterations = 0;
  int maxSimplexIterations = (model_->getNodeCount()) ? maxSimplexIterations_
                                                      : maxSimplexIterationsAtRoot_;
  int maxIterationsInOneSolve = (maxSimplexIterations < 1000000) ? 1000 : 10000;
  // but can't be exactly coin_int_max
  maxSimplexIterations = CoinMin(maxSimplexIterations, COIN_INT_MAX >> 3);
  bool fixGeneralIntegers = false;
  //int maxIterations = maxIterations_;
  int saveSwitches = switches_;
  if ((maxIterations_ % 10) != 0) {
    int digit = maxIterations_ % 10;
    //maxIterations -= digit;
    switches_ |= 65536;
    if ((digit & 3) != 0)
      fixGeneralIntegers = true;
  }

  OsiSolverInterface *solver = cloneBut(6); // was model_->solver()->clone();
#ifdef COIN_HAS_CLP
  OsiClpSolverInterface *clpSolver
    = dynamic_cast< OsiClpSolverInterface * >(solver);
  if (clpSolver) {
    ClpSimplex *clpSimplex = clpSolver->getModelPtr();
    int oneSolveIts = clpSimplex->maximumIterations();
    oneSolveIts = CoinMin(1000 + 2 * (clpSimplex->numberRows() + clpSimplex->numberColumns()), oneSolveIts);
    if (maxSimplexIterations > 1000000)
      maxIterationsInOneSolve = oneSolveIts;
    clpSimplex->setMaximumIterations(oneSolveIts);
    if (!nodes) {
      // say give up easily
      clpSimplex->setMoreSpecialOptions(clpSimplex->moreSpecialOptions() | 64);
    } else {
      // get ray
      int specialOptions = clpSimplex->specialOptions();
      specialOptions &= ~0x3100000;
      specialOptions |= 32;
      clpSimplex->setSpecialOptions(specialOptions);
      clpSolver->setSpecialOptions(clpSolver->specialOptions() | 1048576);
      if ((model_->moreSpecialOptions() & 16777216) != 0) {
        // cutoff is constraint
        clpSolver->setDblParam(OsiDualObjectiveLimit, COIN_DBL_MAX);
      }
    }
  }
#endif
  const double *lower = solver->getColLower();
  const double *upper = solver->getColUpper();
  const double *rowLower = solver->getRowLower();
  const double *rowUpper = solver->getRowUpper();
  const double *solution = solver->getColSolution();
  const double *objective = solver->getObjCoefficients();
  double integerTolerance = model_->getDblParam(CbcModel::CbcIntegerTolerance);
  double primalTolerance;
  solver->getDblParam(OsiPrimalTolerance, primalTolerance);

  int numberRows = matrix_.getNumRows();
  assert(numberRows <= solver->getNumRows());
  int numberIntegers = model_->numberIntegers();
  const int *integerVariable = model_->integerVariable();
  double direction = solver->getObjSense(); // 1 for min, -1 for max
  double newSolutionValue = direction * solver->getObjValue();
  int returnCode = 0;
  // Column copy
  const double *element = matrix_.getElements();
  const int *row = matrix_.getIndices();
  const CoinBigIndex *columnStart = matrix_.getVectorStarts();
  const int *columnLength = matrix_.getVectorLengths();
#ifdef DIVE_FIX_BINARY_VARIABLES
  // Row copy
  const double *elementByRow = matrixByRow_.getElements();
  const int *column = matrixByRow_.getIndices();
  const CoinBigIndex *rowStart = matrixByRow_.getVectorStarts();
  const int *rowLength = matrixByRow_.getVectorLengths();
#endif

  // Get solution array for heuristic solution
  int numberColumns = solver->getNumCols();
  memcpy(newSolution, solution, numberColumns * sizeof(double));

  // vectors to store the latest variables fixed at their bounds
  int *columnFixed = new int[numberIntegers + numberColumns];
  int *back = columnFixed + numberIntegers;
  double *originalBound = new double[numberIntegers + 2 * numberColumns];
  double *lowerBefore = originalBound + numberIntegers;
  double *upperBefore = lowerBefore + numberColumns;
  memcpy(lowerBefore, lower, numberColumns * sizeof(double));
  memcpy(upperBefore, upper, numberColumns * sizeof(double));
  double *lastDjs = newSolution + numberColumns;
  bool *fixedAtLowerBound = new bool[numberIntegers];
  PseudoReducedCost *candidate = new PseudoReducedCost[numberIntegers];
  double *random = new double[numberIntegers];

  int maxNumberAtBoundToFix = static_cast< int >(floor(percentageToFix_ * numberIntegers));
  assert(!maxNumberAtBoundToFix || !nodes);

  // count how many fractional variables
  int numberFractionalVariables = 0;
  for (int i = 0; i < numberColumns; i++)
    back[i] = -1;
  for (int i = 0; i < numberIntegers; i++) {
    random[i] = randomNumberGenerator_.randomDouble() + 0.3;
    int iColumn = integerVariable[i];
    if (!isHeuristicInteger(solver, iColumn))
      continue;
    back[iColumn] = i;
    double value = newSolution[iColumn];
    // clean
    value = CoinMin(value, upperBefore[iColumn]);
    value = CoinMax(value, lowerBefore[iColumn]);
    newSolution[iColumn] = value;
    if (fabs(floor(value + 0.5) - value) > integerTolerance) {
      numberFractionalVariables++;
    }
  }

  const double *reducedCost = NULL;
  // See if not NLP
  if (!model_->solverCharacteristics() || model_->solverCharacteristics()->reducedCostsAccurate())
    reducedCost = solver->getReducedCost();

  int iteration = 0;
  int numberAtBoundFixed = 0;
#if DIVE_PRINT > 1
  int numberGeneralFixed = 0; // fixed as satisfied but not at bound
  int numberReducedCostFixed = 0;
#endif
  while (numberFractionalVariables) {
    iteration++;

    // initialize any data
    initializeData();

    // select a fractional variable to bound
    int bestColumn = -1;
    int bestRound; // -1 rounds down, +1 rounds up
    bool canRound = selectVariableToBranch(solver, newSolution,
      bestColumn, bestRound);
    // if the solution is not trivially roundable, we don't try to round;
    // if the solution is trivially roundable, we try to round. However,
    // if the rounded solution is worse than the current incumbent,
    // then we don't round and proceed normally. In this case, the
    // bestColumn will be a trivially roundable variable
    if (canRound) {
      // check if by rounding all fractional variables
      // we get a solution with an objective value
      // better than the current best integer solution
      double delta = 0.0;
      for (int i = 0; i < numberIntegers; i++) {
        int iColumn = integerVariable[i];
        if (!isHeuristicInteger(solver, iColumn))
          continue;
        double value = newSolution[iColumn];
        if (fabs(floor(value + 0.5) - value) > integerTolerance) {
          assert(downLocks_[i] == 0 || upLocks_[i] == 0);
          double obj = objective[iColumn];
          if (downLocks_[i] == 0 && upLocks_[i] == 0) {
            if (direction * obj >= 0.0)
              delta += (floor(value) - value) * obj;
            else
              delta += (ceil(value) - value) * obj;
          } else if (downLocks_[i] == 0)
            delta += (floor(value) - value) * obj;
          else
            delta += (ceil(value) - value) * obj;
        }
      }
      if (direction * (solver->getObjValue() + delta) < solutionValue) {
#if DIVE_PRINT
        nRoundFeasible++;
#endif
        if (!nodes || bestColumn < 0) {
          // Round all the fractional variables
          for (int i = 0; i < numberIntegers; i++) {
            int iColumn = integerVariable[i];
            if (!isHeuristicInteger(solver, iColumn))
              continue;
            double value = newSolution[iColumn];
            if (fabs(floor(value + 0.5) - value) > integerTolerance) {
              assert(downLocks_[i] == 0 || upLocks_[i] == 0);
              if (downLocks_[i] == 0 && upLocks_[i] == 0) {
                if (direction * objective[iColumn] >= 0.0)
                  newSolution[iColumn] = floor(value);
                else
                  newSolution[iColumn] = ceil(value);
              } else if (downLocks_[i] == 0)
                newSolution[iColumn] = floor(value);
              else
                newSolution[iColumn] = ceil(value);
            }
          }
          break;
        } else {
          // can't round if going to use in branching
          int i;
          for (i = 0; i < numberIntegers; i++) {
            int iColumn = integerVariable[i];
            if (!isHeuristicInteger(solver, iColumn))
              continue;
            double value = newSolution[bestColumn];
            if (fabs(floor(value + 0.5) - value) > integerTolerance) {
              if (iColumn == bestColumn) {
                assert(downLocks_[i] == 0 || upLocks_[i] == 0);
                double obj = objective[bestColumn];
                if (downLocks_[i] == 0 && upLocks_[i] == 0) {
                  if (direction * obj >= 0.0)
                    bestRound = -1;
                  else
                    bestRound = 1;
                } else if (downLocks_[i] == 0)
                  bestRound = -1;
                else
                  bestRound = 1;
                break;
              }
            }
          }
        }
      }
#if DIVE_PRINT
      else
        nRoundInfeasible++;
#endif
    }

    // do reduced cost fixing
#if DIVE_PRINT > 1
    numberReducedCostFixed = reducedCostFix(solver);
#else
    reducedCostFix(solver);
#endif

    numberAtBoundFixed = 0;
#if DIVE_PRINT > 1
    numberGeneralFixed = 0; // fixed as satisfied but not at bound
#endif
#ifdef DIVE_FIX_BINARY_VARIABLES
    // fix binary variables based on pseudo reduced cost
    if (binVarIndex_.size()) {
      int cnt = 0;
      int n = static_cast< int >(binVarIndex_.size());
      for (int j = 0; j < n; j++) {
        int iColumn1 = binVarIndex_[j];
        double value = newSolution[iColumn1];
        if (fabs(value) <= integerTolerance && lower[iColumn1] != upper[iColumn1]) {
          double maxPseudoReducedCost = 0.0;
#ifdef DIVE_DEBUG
          std::cout << "iColumn1 = " << iColumn1 << ", value = " << value << std::endl;
#endif
          int iRow = vbRowIndex_[j];
          double chosenValue = 0.0;
          for (int k = rowStart[iRow]; k < rowStart[iRow] + rowLength[iRow]; k++) {
            int iColumn2 = column[k];
#ifdef DIVE_DEBUG
            std::cout << "iColumn2 = " << iColumn2 << std::endl;
#endif
            if (iColumn1 != iColumn2) {
              double pseudoReducedCost = fabs(reducedCost[iColumn2] * elementByRow[k]);
#ifdef DIVE_DEBUG
              int k2;
              for (k2 = rowStart[iRow]; k2 < rowStart[iRow] + rowLength[iRow]; k2++) {
                if (column[k2] == iColumn1)
                  break;
              }
              std::cout << "reducedCost[" << iColumn2 << "] = "
                        << reducedCost[iColumn2]
                        << ", elementByRow[" << iColumn2 << "] = " << elementByRow[k]
                        << ", elementByRow[" << iColumn1 << "] = " << elementByRow[k2]
                        << ", pseudoRedCost = " << pseudoReducedCost
                        << std::endl;
#endif
              if (pseudoReducedCost > maxPseudoReducedCost)
                maxPseudoReducedCost = pseudoReducedCost;
            } else {
              // save value
              chosenValue = fabs(elementByRow[k]);
            }
          }
          assert(chosenValue);
          maxPseudoReducedCost /= chosenValue;
#ifdef DIVE_DEBUG
          std::cout << ", maxPseudoRedCost = " << maxPseudoReducedCost << std::endl;
#endif
          candidate[cnt].var = iColumn1;
          candidate[cnt++].pseudoRedCost = maxPseudoReducedCost;
        }
      }
#ifdef DIVE_DEBUG
      std::cout << "candidates for rounding = " << cnt << std::endl;
#endif
      std::sort(candidate, candidate + cnt, compareBinaryVars);
      for (int i = 0; i < cnt; i++) {
        int iColumn = candidate[i].var;
        if (numberAtBoundFixed < maxNumberAtBoundToFix) {
          columnFixed[numberAtBoundFixed] = iColumn;
          originalBound[numberAtBoundFixed] = upper[iColumn];
          fixedAtLowerBound[numberAtBoundFixed] = true;
          solver->setColUpper(iColumn, lower[iColumn]);
          numberAtBoundFixed++;
          if (numberAtBoundFixed == maxNumberAtBoundToFix)
            break;
        }
      }
    }
#endif

    // fix other integer variables that are at their bounds
    int cnt = 0;
#ifdef GAP
    double gap = 1.0e30;
#endif
    int fixPriority = COIN_INT_MAX;
    if (reducedCost && true) {
#ifndef JJF_ONE
      cnt = fixOtherVariables(solver, solution, candidate, random);
      if (priority_) {
        for (int i = 0; i < cnt; i++) {
          int iColumn = candidate[i].var;
          if (upper[iColumn] > lower[iColumn]) {
            int j = back[iColumn];
            fixPriority = CoinMin(fixPriority, static_cast< int >(priority_[j].priority));
          }
        }
      }
#else
#ifdef GAP
      double cutoff = model_->getCutoff();
      if (cutoff < 1.0e20 && false) {
        double direction = solver->getObjSense();
        gap = cutoff - solver->getObjValue() * direction;
        gap *= 0.1; // Fix more if plausible
        double tolerance;
        solver->getDblParam(OsiDualTolerance, tolerance);
        if (gap <= 0.0)
          gap = tolerance;
        gap += 100.0 * tolerance;
      }
      int nOverGap = 0;
#endif
      int numberFree = 0;
      int numberFixed = 0;
      for (int i = 0; i < numberIntegers; i++) {
        int iColumn = integerVariable[i];
        if (!isHeuristicInteger(solver, iColumn))
          continue;
        if (upper[iColumn] > lower[iColumn]) {
          numberFree++;
          if (priority_) {
            fixPriority = CoinMin(fixPriority, static_cast< int >(priority_[i].priority));
          }
          double value = newSolution[iColumn];
          if (fabs(floor(value + 0.5) - value) <= integerTolerance) {
            candidate[cnt].var = iColumn;
            candidate[cnt++].pseudoRedCost = fabs(reducedCost[iColumn] * random[i]);
#ifdef GAP
            if (fabs(reducedCost[iColumn]) > gap)
              nOverGap++;
#endif
          }
        } else {
          numberFixed++;
        }
      }
#ifdef GAP
      int nLeft = maxNumberAtBoundToFix - numberAtBoundFixed;
#ifdef CLP_INVESTIGATE4
      printf("cutoff %g obj %g nover %d - %d free, %d fixed\n",
        cutoff, solver->getObjValue(), nOverGap, numberFree, numberFixed);
#endif
      if (nOverGap > nLeft && true) {
        nOverGap = CoinMin(nOverGap, nLeft + maxNumberAtBoundToFix / 2);
        maxNumberAtBoundToFix += nOverGap - nLeft;
      }
#else
#ifdef CLP_INVESTIGATE4
      printf("cutoff %g obj %g - %d free, %d fixed\n",
        model_->getCutoff(), solver->getObjValue(), numberFree, numberFixed);
#endif
#endif
#endif
    } else {
      for (int i = 0; i < numberIntegers; i++) {
        int iColumn = integerVariable[i];
        if (!isHeuristicInteger(solver, iColumn))
          continue;
        if (upper[iColumn] > lower[iColumn]) {
          if (priority_) {
            fixPriority = CoinMin(fixPriority, static_cast< int >(priority_[i].priority));
          }
          double value = newSolution[iColumn];
          if (fabs(floor(value + 0.5) - value) <= integerTolerance) {
            candidate[cnt].var = iColumn;
            candidate[cnt++].pseudoRedCost = numberIntegers - i;
          }
        }
      }
    }
    std::sort(candidate, candidate + cnt, compareBinaryVars);
    // If getting on fix all
    if (iteration * 3 > maxIterations_ * 2)
      fixPriority = COIN_INT_MAX;
    for (int i = 0; i < cnt; i++) {
      int iColumn = candidate[i].var;
      if (upper[iColumn] > lower[iColumn]) {
        double value = newSolution[iColumn];
        if (fabs(floor(value + 0.5) - value) <= integerTolerance && numberAtBoundFixed < maxNumberAtBoundToFix) {
          // fix the variable at one of its bounds
          if (fabs(lower[iColumn] - value) <= integerTolerance || fixGeneralIntegers) {
            if (priority_) {
              int j = back[iColumn];
              if (priority_[j].priority > fixPriority)
                continue; // skip - only fix ones at high priority
              int thisRound = static_cast< int >(priority_[j].direction);
              if ((thisRound & 1) != 0) {
                // for now force way
                if ((thisRound & 2) != 0)
                  continue;
              }
            }
            if (fabs(lower[iColumn] - value) <= integerTolerance) {
              columnFixed[numberAtBoundFixed] = iColumn;
              originalBound[numberAtBoundFixed] = upper[iColumn];
              fixedAtLowerBound[numberAtBoundFixed] = true;
              solver->setColUpper(iColumn, lower[iColumn]);
            } else {
              // fix to interior value
#if DIVE_PRINT > 1
              numberGeneralFixed++;
#endif
              double fixValue = floor(value + 0.5);
              columnFixed[numberAtBoundFixed] = iColumn;
              originalBound[numberAtBoundFixed] = upper[iColumn];
              fixedAtLowerBound[numberAtBoundFixed] = true;
              solver->setColUpper(iColumn, fixValue);
              numberAtBoundFixed++;
              columnFixed[numberAtBoundFixed] = iColumn;
              originalBound[numberAtBoundFixed] = lower[iColumn];
              fixedAtLowerBound[numberAtBoundFixed] = false;
              solver->setColLower(iColumn, fixValue);
            }
            //if (priority_)
            //printf("fixing %d (priority %d) to lower bound of %g\n",
            //	 iColumn,priority_[back[iColumn]].priority,lower[iColumn]);
            numberAtBoundFixed++;
          } else if (fabs(upper[iColumn] - value) <= integerTolerance) {
            if (priority_) {
              int j = back[iColumn];
              if (priority_[j].priority > fixPriority)
                continue; // skip - only fix ones at high priority
              int thisRound = static_cast< int >(priority_[j].direction);
              if ((thisRound & 1) != 0) {
                // for now force way
                if ((thisRound & 2) == 0)
                  continue;
              }
            }
            columnFixed[numberAtBoundFixed] = iColumn;
            originalBound[numberAtBoundFixed] = lower[iColumn];
            fixedAtLowerBound[numberAtBoundFixed] = false;
            solver->setColLower(iColumn, upper[iColumn]);
            //if (priority_)
            //printf("fixing %d (priority %d) to upper bound of %g\n",
            //	 iColumn,priority_[back[iColumn]].priority,upper[iColumn]);
            numberAtBoundFixed++;
          }
          if (numberAtBoundFixed == maxNumberAtBoundToFix)
            break;
        }
      }
    }

    double originalBoundBestColumn;
    double bestColumnValue;
    int whichWay;
    if (bestColumn >= 0) {
      bestColumnValue = newSolution[bestColumn];
      if (bestRound < 0) {
        originalBoundBestColumn = upper[bestColumn];
        solver->setColUpper(bestColumn, floor(bestColumnValue));
#ifdef DIVE_DEBUG
        if (priority_) {
          printf("setting %d (priority %d) upper bound to %g (%g)\n",
            bestColumn, priority_[back[bestColumn]].priority, floor(bestColumnValue), bestColumnValue);
        }
#endif
        whichWay = 0;
      } else {
        originalBoundBestColumn = lower[bestColumn];
        solver->setColLower(bestColumn, ceil(bestColumnValue));
#ifdef DIVE_DEBUG
        if (priority_) {
          printf("setting %d (priority %d) lower bound to %g (%g)\n",
            bestColumn, priority_[back[bestColumn]].priority, ceil(bestColumnValue), bestColumnValue);
        }
#endif
        whichWay = 1;
      }
    } else {
      break;
    }
    int originalBestRound = bestRound;
    int saveModelOptions = model_->specialOptions();
    while (1) {

      model_->setSpecialOptions(saveModelOptions | 2048);
      solver->resolve();
      numberSimplexIterations += solver->getIterationCount();
#if DIVE_PRINT > 1
      int numberFractionalVariables = 0;
      double sumFractionalVariables = 0.0;
      int numberFixed = 0;
      for (int i = 0; i < numberIntegers; i++) {
        int iColumn = integerVariable[i];
        if (!isHeuristicInteger(solver, iColumn))
          continue;
        double value = newSolution[iColumn];
        double away = fabs(floor(value + 0.5) - value);
        if (away > integerTolerance) {
          numberFractionalVariables++;
          sumFractionalVariables += away;
        }
        if (upper[iColumn] == lower[iColumn])
          numberFixed++;
      }
      printf("pass %d obj %g %s its %d total %d fixed %d +(%d,%d)", iteration,
        solver->getObjValue(), solver->isProvenOptimal() ? "opt" : "infeasible",
        solver->getIterationCount(),
        solver->getIterationCount() + numberSimplexIterations,
        numberFixed,
        numberReducedCostFixed,
        numberAtBoundFixed - numberGeneralFixed);
      if (solver->isProvenOptimal()) {
        printf(" - %d at bound, %d away (sum %g)\n",
          numberIntegers - numberFixed - numberFractionalVariables,
          numberFractionalVariables, sumFractionalVariables);
      } else {
        printf("\n");
        if (fixGeneralIntegers) {
          int digit = maxIterations_ % 10;
          if (digit == 1) {
            // switch off for now
            switches_ = saveSwitches;
            fixGeneralIntegers = false;
          } else if (digit == 2) {
            // switch off always
            switches_ = saveSwitches;
            fixGeneralIntegers = false;
            maxIterations_ -= digit;
          }
        }
      }
#else
      if (!solver->isProvenOptimal()) {
        if (fixGeneralIntegers) {
          int digit = maxIterations_ % 10;
          if (digit == 1) {
            // switch off for now
            switches_ = saveSwitches;
            fixGeneralIntegers = false;
          } else if (digit == 2) {
            // switch off always
            switches_ = saveSwitches;
            fixGeneralIntegers = false;
            maxIterations_ -= digit;
          }
        }
      }
#endif
      model_->setSpecialOptions(saveModelOptions);
      if (!solver->isAbandoned() && !solver->isIterationLimitReached()) {
        //numberSimplexIterations += solver->getIterationCount();
      } else {
        numberSimplexIterations = maxSimplexIterations + 1;
        reasonToStop += 100;
        break;
      }

      if (!solver->isProvenOptimal()) {
        if (nodes) {
          if (solver->isProvenPrimalInfeasible()) {
            if (maxSimplexIterationsAtRoot_ != COIN_INT_MAX) {
              // stop now
              printf("stopping on first infeasibility\n");
              break;
            } else if (cuts) {
              // can do conflict cut
              printf("could do intermediate conflict cut\n");
              bool localCut;
              OsiRowCut *cut = model_->conflictCut(solver, localCut);
              if (cut) {
                if (!localCut) {
                  model_->makePartialCut(cut, solver);
                  cuts[numberCuts++] = cut;
                } else {
                  delete cut;
                }
              }
            }
          } else {
            reasonToStop += 10;
            break;
          }
        }
        if (numberAtBoundFixed > 0) {
          // Remove the bound fix for variables that were at bounds
          for (int i = 0; i < numberAtBoundFixed; i++) {
            int iColFixed = columnFixed[i];
            if (fixedAtLowerBound[i])
              solver->setColUpper(iColFixed, originalBound[i]);
            else
              solver->setColLower(iColFixed, originalBound[i]);
          }
          numberAtBoundFixed = 0;
        } else if (bestRound == originalBestRound) {
          bestRound *= (-1);
          whichWay |= 2;
          if (bestRound < 0) {
            solver->setColLower(bestColumn, originalBoundBestColumn);
            solver->setColUpper(bestColumn, floor(bestColumnValue));
          } else {
            solver->setColLower(bestColumn, ceil(bestColumnValue));
            solver->setColUpper(bestColumn, originalBoundBestColumn);
          }
        } else
          break;
      } else
        break;
    }

    if (!solver->isProvenOptimal() || direction * solver->getObjValue() >= solutionValue) {
      reasonToStop += 1;
    } else if (iteration > maxIterations_) {
      reasonToStop += 2;
    } else if (CoinCpuTime() - time1 > maxTime_) {
      reasonToStop += 3;
    } else if (numberSimplexIterations > maxSimplexIterations) {
      reasonToStop += 4;
      // also switch off
#if DIVE_PRINT
      printf("switching off diving as too many iterations %d, %d allowed\n",
        numberSimplexIterations, maxSimplexIterations);
#endif
      when_ = 0;
    } else if (solver->getIterationCount() > maxIterationsInOneSolve && iteration > 3 && !nodes) {
      reasonToStop += 5;
      // also switch off
#if DIVE_PRINT
      printf("switching off diving one iteration took %d iterations (total %d)\n",
        solver->getIterationCount(), numberSimplexIterations);
#endif
      when_ = 0;
    }

    memcpy(newSolution, solution, numberColumns * sizeof(double));
    numberFractionalVariables = 0;
    double sumFractionalVariables = 0.0;
    for (int i = 0; i < numberIntegers; i++) {
      int iColumn = integerVariable[i];
      if (!isHeuristicInteger(solver, iColumn))
        continue;
      double value = newSolution[iColumn];
      double away = fabs(floor(value + 0.5) - value);
      if (away > integerTolerance) {
        numberFractionalVariables++;
        sumFractionalVariables += away;
      }
    }
    if (nodes) {
      // save information
      //branchValues[numberNodes]=bestColumnValue;
      //statuses[numberNodes]=whichWay+(bestColumn<<2);
      //bases[numberNodes]=solver->getWarmStart();
      ClpSimplex *simplex = clpSolver->getModelPtr();
      CbcSubProblem *sub = new CbcSubProblem(clpSolver, lowerBefore, upperBefore,
        simplex->statusArray(), numberNodes);
      nodes[numberNodes] = sub;
      // other stuff
      sub->branchValue_ = bestColumnValue;
      sub->problemStatus_ = whichWay;
      sub->branchVariable_ = bestColumn;
      sub->objectiveValue_ = simplex->objectiveValue();
      sub->sumInfeasibilities_ = sumFractionalVariables;
      sub->numberInfeasibilities_ = numberFractionalVariables;
      printf("DiveNode %d column %d way %d bvalue %g obj %g\n",
        numberNodes, sub->branchVariable_, sub->problemStatus_,
        sub->branchValue_, sub->objectiveValue_);
      numberNodes++;
      if (solver->isProvenOptimal()) {
        memcpy(lastDjs, solver->getReducedCost(), numberColumns * sizeof(double));
        memcpy(lowerBefore, lower, numberColumns * sizeof(double));
        memcpy(upperBefore, upper, numberColumns * sizeof(double));
      }
    }
    if (!numberFractionalVariables || reasonToStop)
      break;
  }
  if (nodes) {
    printf("Exiting dive for reason %d\n", reasonToStop);
    if (reasonToStop > 1) {
      printf("problems in diving\n");
      int whichWay = nodes[numberNodes - 1]->problemStatus_;
      CbcSubProblem *sub;
      if ((whichWay & 2) == 0) {
        // leave both ways
        sub = new CbcSubProblem(*nodes[numberNodes - 1]);
        nodes[numberNodes++] = sub;
      } else {
        sub = nodes[numberNodes - 1];
      }
      if ((whichWay & 1) == 0)
        sub->problemStatus_ = whichWay | 1;
      else
        sub->problemStatus_ = whichWay & ~1;
    }
    if (!numberNodes) {
      // was good at start! - create fake
      clpSolver->resolve();
      numberSimplexIterations += clpSolver->getIterationCount();
      ClpSimplex *simplex = clpSolver->getModelPtr();
      CbcSubProblem *sub = new CbcSubProblem(clpSolver, lowerBefore, upperBefore,
        simplex->statusArray(), numberNodes);
      nodes[numberNodes] = sub;
      // other stuff
      sub->branchValue_ = 0.0;
      sub->problemStatus_ = 0;
      sub->branchVariable_ = -1;
      sub->objectiveValue_ = simplex->objectiveValue();
      sub->sumInfeasibilities_ = 0.0;
      sub->numberInfeasibilities_ = 0;
      printf("DiveNode %d column %d way %d bvalue %g obj %g\n",
        numberNodes, sub->branchVariable_, sub->problemStatus_,
        sub->branchValue_, sub->objectiveValue_);
      numberNodes++;
      assert(solver->isProvenOptimal());
    }
    nodes[numberNodes - 1]->problemStatus_ |= 256 * reasonToStop;
    // use djs as well
    if (solver->isProvenPrimalInfeasible() && cuts) {
      // can do conflict cut and re-order
      printf("could do final conflict cut\n");
      bool localCut;
      OsiRowCut *cut = model_->conflictCut(solver, localCut);
      if (cut) {
        printf("cut - need to use conflict and previous djs\n");
        if (!localCut) {
          model_->makePartialCut(cut, solver);
          cuts[numberCuts++] = cut;
        } else {
          delete cut;
        }
      } else {
        printf("bad conflict - just use previous djs\n");
      }
    }
  }

  // re-compute new solution value
  double objOffset = 0.0;
  solver->getDblParam(OsiObjOffset, objOffset);
  newSolutionValue = -objOffset;
  for (int i = 0; i < numberColumns; i++)
    newSolutionValue += objective[i] * newSolution[i];
  newSolutionValue *= direction;
  //printf("new solution value %g %g\n",newSolutionValue,solutionValue);
  if (newSolutionValue < solutionValue && !reasonToStop) {
    double *rowActivity = new double[numberRows];
    memset(rowActivity, 0, numberRows * sizeof(double));
    // paranoid check
    memset(rowActivity, 0, numberRows * sizeof(double));
    for (int i = 0; i < numberColumns; i++) {
      CoinBigIndex j;
      double value = newSolution[i];
      if (value) {
        for (j = columnStart[i];
             j < columnStart[i] + columnLength[i]; j++) {
          int iRow = row[j];
          rowActivity[iRow] += value * element[j];
        }
      }
    }
    // check was approximately feasible
    bool feasible = true;
    for (int i = 0; i < numberRows; i++) {
      if (rowActivity[i] < rowLower[i]) {
        if (rowActivity[i] < rowLower[i] - 1000.0 * primalTolerance)
          feasible = false;
      } else if (rowActivity[i] > rowUpper[i]) {
        if (rowActivity[i] > rowUpper[i] + 1000.0 * primalTolerance)
          feasible = false;
      }
    }
    for (int i = 0; i < numberIntegers; i++) {
      int iColumn = integerVariable[i];
      if (!isHeuristicInteger(solver, iColumn))
        continue;
      double value = newSolution[iColumn];
      if (fabs(floor(value + 0.5) - value) > integerTolerance) {
        feasible = false;
        break;
      }
    }
    if (feasible) {
      // new solution
      solutionValue = newSolutionValue;
      //printf("** Solution of %g found by CbcHeuristicDive\n",newSolutionValue);
      //if (cuts)
      //clpSolver->getModelPtr()->writeMps("good8.mps", 2);
      returnCode = 1;
    } else {
      // Can easily happen
      //printf("Debug CbcHeuristicDive giving bad solution\n");
    }
    delete[] rowActivity;
  }

#if DIVE_PRINT
  std::cout << heuristicName_
            << " nRoundInfeasible = " << nRoundInfeasible
            << ", nRoundFeasible = " << nRoundFeasible
            << ", returnCode = " << returnCode
            << ", reasonToStop = " << reasonToStop
            << ", simplexIts = " << numberSimplexIterations
            << ", iterations = " << iteration << std::endl;
#endif

  delete[] columnFixed;
  delete[] originalBound;
  delete[] fixedAtLowerBound;
  delete[] candidate;
  delete[] random;
  delete[] downArray_;
  downArray_ = NULL;
  delete[] upArray_;
  upArray_ = NULL;
  delete solver;
  switches_ = saveSwitches;
  return returnCode;
}
// See if diving will give better solution
// Sets value of solution
// Returns 1 if solution, 0 if not
int CbcHeuristicDive::solution(double &solutionValue,
  double *betterSolution)
{
  int nodeCount = model_->getNodeCount();
  if (feasibilityPumpOptions_ > 0 && (nodeCount % feasibilityPumpOptions_) != 0)
    return 0;
  ++numCouldRun_;

  // test if the heuristic can run
  if (!canHeuristicRun())
    return 0;

#ifdef JJF_ZERO
  // See if to do
  if (!when() || (when() % 10 == 1 && model_->phase() != 1) || (when() % 10 == 2 && (model_->phase() != 2 && model_->phase() != 3)))
    return 0; // switched off
#endif
#ifdef HEURISTIC_INFORM
  printf("Entering heuristic %s - nRuns %d numCould %d when %d\n",
    heuristicName(), numRuns_, numCouldRun_, when_);
#endif
#ifdef DIVE_DEBUG
  std::cout << "solutionValue = " << solutionValue << std::endl;
#endif
  // Get solution array for heuristic solution
  int numberColumns = model_->solver()->getNumCols();
  double *newSolution = CoinCopyOfArray(model_->solver()->getColSolution(),
    numberColumns);
  int numberCuts = 0;
  int numberNodes = -1;
  CbcSubProblem **nodes = NULL;
  int returnCode = solution(solutionValue, numberNodes, numberCuts,
    NULL, nodes,
    newSolution);
  if (returnCode == 1)
    memcpy(betterSolution, newSolution, numberColumns * sizeof(double));

  delete[] newSolution;
  return returnCode;
}
/* returns 0 if no solution, 1 if valid solution
   with better objective value than one passed in
   also returns list of nodes
   This does Fractional Diving
*/
int CbcHeuristicDive::fathom(CbcModel *model, int &numberNodes,
  CbcSubProblem **&nodes)
{
  double solutionValue = model->getCutoff();
  numberNodes = 0;
  // Get solution array for heuristic solution
  int numberColumns = model_->solver()->getNumCols();
  double *newSolution = new double[4 * numberColumns];
  double *lastDjs = newSolution + numberColumns;
  double *originalLower = lastDjs + numberColumns;
  double *originalUpper = originalLower + numberColumns;
  memcpy(originalLower, model_->solver()->getColLower(),
    numberColumns * sizeof(double));
  memcpy(originalUpper, model_->solver()->getColUpper(),
    numberColumns * sizeof(double));
  int numberCuts = 0;
  OsiRowCut **cuts = NULL; //new OsiRowCut * [maxIterations_];
  nodes = new CbcSubProblem *[maxIterations_ + 2];
  int returnCode = solution(solutionValue, numberNodes, numberCuts,
    cuts, nodes,
    newSolution);

  if (returnCode == 1) {
    // copy to best solution ? or put in solver
    printf("Solution from heuristic fathom\n");
  }
  int numberFeasibleNodes = numberNodes;
  if (returnCode != 1)
    numberFeasibleNodes--;
  if (numberFeasibleNodes > 0) {
    CoinWarmStartBasis *basis = nodes[numberFeasibleNodes - 1]->status_;
    //double * sort = new double [numberFeasibleNodes];
    //int * whichNode = new int [numberFeasibleNodes];
    //int numberNodesNew=0;
    // use djs on previous unless feasible
    for (int iNode = 0; iNode < numberFeasibleNodes; iNode++) {
      CbcSubProblem *sub = nodes[iNode];
      double branchValue = sub->branchValue_;
      int iStatus = sub->problemStatus_;
      int iColumn = sub->branchVariable_;
      bool secondBranch = (iStatus & 2) != 0;
      bool branchUp;
      if (!secondBranch)
        branchUp = (iStatus & 1) != 0;
      else
        branchUp = (iStatus & 1) == 0;
      double djValue = lastDjs[iColumn];
      sub->djValue_ = fabs(djValue);
      if (!branchUp && floor(branchValue) == originalLower[iColumn]
        && basis->getStructStatus(iColumn) == CoinWarmStartBasis::atLowerBound) {
        if (djValue > 0.0) {
          // naturally goes to LB
          printf("ignoring branch down on %d (node %d) from value of %g - branch was %s - dj %g\n",
            iColumn, iNode, branchValue, secondBranch ? "second" : "first",
            djValue);
          sub->problemStatus_ |= 4;
          //} else {
          // put on list
          //sort[numberNodesNew]=djValue;
          //whichNode[numberNodesNew++]=iNode;
        }
      } else if (branchUp && ceil(branchValue) == originalUpper[iColumn]
        && basis->getStructStatus(iColumn) == CoinWarmStartBasis::atUpperBound) {
        if (djValue < 0.0) {
          // naturally goes to UB
          printf("ignoring branch up on %d (node %d) from value of %g - branch was %s - dj %g\n",
            iColumn, iNode, branchValue, secondBranch ? "second" : "first",
            djValue);
          sub->problemStatus_ |= 4;
          //} else {
          // put on list
          //sort[numberNodesNew]=-djValue;
          //whichNode[numberNodesNew++]=iNode;
        }
      }
    }
    // use conflict to order nodes
    for (int iCut = 0; iCut < numberCuts; iCut++) {
    }
    //CoinSort_2(sort,sort+numberNodesNew,whichNode);
    // create nodes
    // last node will have one way already done
  }
  for (int iCut = 0; iCut < numberCuts; iCut++) {
    delete cuts[iCut];
  }
  delete[] cuts;
  delete[] newSolution;
  return returnCode;
}

// Validate model i.e. sets when_ to 0 if necessary (may be NULL)
void CbcHeuristicDive::validate()
{
  if (model_ && (when() % 100) < 10) {
    if (model_->numberIntegers() != model_->numberObjects() && (model_->numberObjects() || (model_->specialOptions() & 1024) == 0)) {
      int numberOdd = 0;
      for (int i = 0; i < model_->numberObjects(); i++) {
        if (!model_->object(i)->canDoHeuristics())
          numberOdd++;
      }
      if (numberOdd)
        setWhen(0);
    }
  }

  int numberIntegers = model_->numberIntegers();
  const int *integerVariable = model_->integerVariable();
  delete[] downLocks_;
  delete[] upLocks_;
  downLocks_ = new unsigned short[numberIntegers];
  upLocks_ = new unsigned short[numberIntegers];
  // Column copy
  const double *element = matrix_.getElements();
  const int *row = matrix_.getIndices();
  const CoinBigIndex *columnStart = matrix_.getVectorStarts();
  OsiSolverInterface *solver = model_->solver();
  const int *columnLength = matrix_.getVectorLengths();
  const double *rowLower = solver->getRowLower();
  const double *rowUpper = solver->getRowUpper();
  for (int i = 0; i < numberIntegers; i++) {
    int iColumn = integerVariable[i];
    if (!isHeuristicInteger(solver, iColumn))
      continue;
    int down = 0;
    int up = 0;
    if (columnLength[iColumn] > 65535) {
      setWhen(0);
      break; // unlikely to work
    }
    for (CoinBigIndex j = columnStart[iColumn];
         j < columnStart[iColumn] + columnLength[iColumn]; j++) {
      int iRow = row[j];
      if (rowLower[iRow] > -1.0e20 && rowUpper[iRow] < 1.0e20) {
        up++;
        down++;
      } else if (element[j] > 0.0) {
        if (rowUpper[iRow] < 1.0e20)
          up++;
        else
          down++;
      } else {
        if (rowLower[iRow] > -1.0e20)
          up++;
        else
          down++;
      }
    }
    downLocks_[i] = static_cast< unsigned short >(down);
    upLocks_[i] = static_cast< unsigned short >(up);
  }

#ifdef DIVE_FIX_BINARY_VARIABLES
  selectBinaryVariables();
#endif
}

// Select candidate binary variables for fixing
void CbcHeuristicDive::selectBinaryVariables()
{
  // Row copy
  const double *elementByRow = matrixByRow_.getElements();
  const int *column = matrixByRow_.getIndices();
  const CoinBigIndex *rowStart = matrixByRow_.getVectorStarts();
  const int *rowLength = matrixByRow_.getVectorLengths();

  const int numberRows = matrixByRow_.getNumRows();
  const int numberCols = matrixByRow_.getNumCols();

  OsiSolverInterface *solver = model_->solver();
  const double *lower = solver->getColLower();
  const double *upper = solver->getColUpper();
  const double *rowLower = solver->getRowLower();
  const double *rowUpper = solver->getRowUpper();

  //  const char * integerType = model_->integerType();

  //  const int numberIntegers = model_->numberIntegers();
  //  const int * integerVariable = model_->integerVariable();
  const double *objective = solver->getObjCoefficients();

  // vector to store the row number of variable bound rows
  int *rowIndexes = new int[numberCols];
  memset(rowIndexes, -1, numberCols * sizeof(int));

  for (int i = 0; i < numberRows; i++) {
    int positiveBinary = -1;
    int negativeBinary = -1;
    int nPositiveOther = 0;
    int nNegativeOther = 0;
    for (CoinBigIndex k = rowStart[i]; k < rowStart[i] + rowLength[i]; k++) {
      int iColumn = column[k];
      if (isHeuristicInteger(solver, iColumn) && lower[iColumn] == 0.0 && upper[iColumn] == 1.0 && objective[iColumn] == 0.0 && elementByRow[k] > 0.0 && positiveBinary < 0)
        positiveBinary = iColumn;
      else if (isHeuristicInteger(solver, iColumn) && lower[iColumn] == 0.0 && upper[iColumn] == 1.0 && objective[iColumn] == 0.0 && elementByRow[k] < 0.0 && negativeBinary < 0)
        negativeBinary = iColumn;
      else if ((elementByRow[k] > 0.0 && lower[iColumn] >= 0.0) || (elementByRow[k] < 0.0 && upper[iColumn] <= 0.0))
        nPositiveOther++;
      else if ((elementByRow[k] > 0.0 && lower[iColumn] <= 0.0) || (elementByRow[k] < 0.0 && upper[iColumn] >= 0.0))
        nNegativeOther++;
      if (nPositiveOther > 0 && nNegativeOther > 0)
        break;
    }
    int binVar = -1;
    if (positiveBinary >= 0 && (negativeBinary >= 0 || nNegativeOther > 0) && nPositiveOther == 0 && rowLower[i] == 0.0 && rowUpper[i] > 0.0)
      binVar = positiveBinary;
    else if (negativeBinary >= 0 && (positiveBinary >= 0 || nPositiveOther > 0) && nNegativeOther == 0 && rowLower[i] < 0.0 && rowUpper[i] == 0.0)
      binVar = negativeBinary;
    if (binVar >= 0) {
      if (rowIndexes[binVar] == -1)
        rowIndexes[binVar] = i;
      else if (rowIndexes[binVar] >= 0)
        rowIndexes[binVar] = -2;
    }
  }

  for (int j = 0; j < numberCols; j++) {
    if (rowIndexes[j] >= 0) {
      binVarIndex_.push_back(j);
      vbRowIndex_.push_back(rowIndexes[j]);
    }
  }

#ifdef DIVE_DEBUG
  std::cout << "number vub Binary = " << binVarIndex_.size() << std::endl;
#endif

  delete[] rowIndexes;
}

/*
  Perform reduced cost fixing on integer variables.

  The variables in question are already nonbasic at bound. We're just nailing
  down the current situation.
*/

int CbcHeuristicDive::reducedCostFix(OsiSolverInterface *solver)

{
  //return 0; // temp
#ifndef JJF_ONE
  if (!model_->solverCharacteristics()->reducedCostsAccurate())
    return 0; //NLP
#endif
  double cutoff = model_->getCutoff();
  if (cutoff > 1.0e20)
    return 0;
#ifdef DIVE_DEBUG
  std::cout << "cutoff = " << cutoff << std::endl;
#endif
  double direction = solver->getObjSense();
  double gap = cutoff - solver->getObjValue() * direction;
  gap *= 0.5; // Fix more
  double tolerance;
  solver->getDblParam(OsiDualTolerance, tolerance);
  if (gap <= 0.0)
    gap = tolerance; //return 0;
  gap += 100.0 * tolerance;
  double integerTolerance = model_->getDblParam(CbcModel::CbcIntegerTolerance);

  const double *lower = solver->getColLower();
  const double *upper = solver->getColUpper();
  const double *solution = solver->getColSolution();
  const double *reducedCost = solver->getReducedCost();

  int numberIntegers = model_->numberIntegers();
  const int *integerVariable = model_->integerVariable();

  int numberFixed = 0;

#ifdef COIN_HAS_CLP
  OsiClpSolverInterface *clpSolver
    = dynamic_cast< OsiClpSolverInterface * >(solver);
  ClpSimplex *clpSimplex = NULL;
  if (clpSolver)
    clpSimplex = clpSolver->getModelPtr();
#endif
  for (int i = 0; i < numberIntegers; i++) {
    int iColumn = integerVariable[i];
    if (!isHeuristicInteger(solver, iColumn))
      continue;
    double djValue = direction * reducedCost[iColumn];
    if (upper[iColumn] - lower[iColumn] > integerTolerance) {
      if (solution[iColumn] < lower[iColumn] + integerTolerance && djValue > gap) {
#ifdef COIN_HAS_CLP
        // may just have been fixed before
        if (clpSimplex) {
          if (clpSimplex->getColumnStatus(iColumn) == ClpSimplex::basic) {
#ifdef COIN_DEVELOP
            printf("DJfix %d has status of %d, dj of %g gap %g, bounds %g %g\n",
              iColumn, clpSimplex->getColumnStatus(iColumn),
              djValue, gap, lower[iColumn], upper[iColumn]);
#endif
          } else {
            //assert(clpSimplex->getColumnStatus(iColumn) == ClpSimplex::atLowerBound || clpSimplex->getColumnStatus(iColumn) == ClpSimplex::isFixed);
          }
        }
#endif
        solver->setColUpper(iColumn, lower[iColumn]);
        numberFixed++;
      } else if (solution[iColumn] > upper[iColumn] - integerTolerance && -djValue > gap) {
#ifdef COIN_HAS_CLP
        // may just have been fixed before
        if (clpSimplex) {
          if (clpSimplex->getColumnStatus(iColumn) == ClpSimplex::basic) {
#ifdef COIN_DEVELOP
            printf("DJfix %d has status of %d, dj of %g gap %g, bounds %g %g\n",
              iColumn, clpSimplex->getColumnStatus(iColumn),
              djValue, gap, lower[iColumn], upper[iColumn]);
#endif
          } else {
            //assert(clpSimplex->getColumnStatus(iColumn) == ClpSimplex::atUpperBound || clpSimplex->getColumnStatus(iColumn) == ClpSimplex::isFixed);
          }
        }
#endif
        solver->setColLower(iColumn, upper[iColumn]);
        numberFixed++;
      }
    }
  }
  return numberFixed;
}
// Fix other variables at bounds
int CbcHeuristicDive::fixOtherVariables(OsiSolverInterface *solver,
  const double *solution,
  PseudoReducedCost *candidate,
  const double *random)
{
  const double *lower = solver->getColLower();
  const double *upper = solver->getColUpper();
  double integerTolerance = model_->getDblParam(CbcModel::CbcIntegerTolerance);
  double primalTolerance;
  solver->getDblParam(OsiPrimalTolerance, primalTolerance);

  int numberIntegers = model_->numberIntegers();
  const int *integerVariable = model_->integerVariable();
  const double *reducedCost = solver->getReducedCost();
  // fix other integer variables that are at their bounds
  int cnt = 0;
#ifdef GAP
  double direction = solver->getObjSense(); // 1 for min, -1 for max
  double gap = 1.0e30;
#endif
#ifdef GAP
  double cutoff = model_->getCutoff();
  if (cutoff < 1.0e20 && false) {
    double direction = solver->getObjSense();
    gap = cutoff - solver->getObjValue() * direction;
    gap *= 0.1; // Fix more if plausible
    double tolerance;
    solver->getDblParam(OsiDualTolerance, tolerance);
    if (gap <= 0.0)
      gap = tolerance;
    gap += 100.0 * tolerance;
  }
  int nOverGap = 0;
#endif
#ifdef CLP_INVESTIGATE4
  int numberFree = 0;
  int numberFixedAlready = 0;
#endif
  for (int i = 0; i < numberIntegers; i++) {
    int iColumn = integerVariable[i];
    if (!isHeuristicInteger(solver, iColumn))
      continue;
    if (upper[iColumn] > lower[iColumn]) {
#ifdef CLP_INVESTIGATE4
      numberFree++;
#endif
      double value = solution[iColumn];
      if (fabs(floor(value + 0.5) - value) <= integerTolerance) {
        candidate[cnt].var = iColumn;
        candidate[cnt++].pseudoRedCost = fabs(reducedCost[iColumn] * random[i]);
#ifdef GAP
        if (fabs(reducedCost[iColumn]) > gap)
          nOverGap++;
#endif
      }
    } else {
#ifdef CLP_INVESTIGATE4
      numberFixedAlready++;
#endif
    }
  }
#ifdef GAP
  int nLeft = maxNumberToFix - numberFixedAlready;
#ifdef CLP_INVESTIGATE4
  printf("cutoff %g obj %g nover %d - %d free, %d fixed\n",
    cutoff, solver->getObjValue(), nOverGap, numberFree,
    numberFixedAlready);
#endif
  if (nOverGap > nLeft && true) {
    nOverGap = CoinMin(nOverGap, nLeft + maxNumberToFix / 2);
    maxNumberToFix += nOverGap - nLeft;
  }
#else
#ifdef CLP_INVESTIGATE4
  printf("cutoff %g obj %g - %d free, %d fixed\n",
    model_->getCutoff(), solver->getObjValue(), numberFree,
    numberFixedAlready);
#endif
#endif
  return cnt;
}

/* vi: softtabstop=2 shiftwidth=2 expandtab tabstop=2
*/
