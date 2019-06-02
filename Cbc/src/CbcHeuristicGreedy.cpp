/* $Id$ */
// Copyright (C) 2005, International Business Machines
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
#include "CbcModel.hpp"
#include "CbcStrategy.hpp"
#include "CbcHeuristicGreedy.hpp"
#include "CoinSort.hpp"
#include "CglPreProcess.hpp"
// Default Constructor
CbcHeuristicGreedyCover::CbcHeuristicGreedyCover()
  : CbcHeuristic()
{
  // matrix  will automatically be empty
  originalNumberRows_ = 0;
  algorithm_ = 0;
  numberTimes_ = 100;
}

// Constructor from model
CbcHeuristicGreedyCover::CbcHeuristicGreedyCover(CbcModel &model)
  : CbcHeuristic(model)
{
  gutsOfConstructor(&model);
  algorithm_ = 0;
  numberTimes_ = 100;
  whereFrom_ = 1;
}

// Destructor
CbcHeuristicGreedyCover::~CbcHeuristicGreedyCover()
{
}

// Clone
CbcHeuristic *
CbcHeuristicGreedyCover::clone() const
{
  return new CbcHeuristicGreedyCover(*this);
}
// Guts of constructor from a CbcModel
void CbcHeuristicGreedyCover::gutsOfConstructor(CbcModel *model)
{
  model_ = model;
  // Get a copy of original matrix
  assert(model->solver());
  if (model->solver()->getNumRows()) {
    matrix_ = *model->solver()->getMatrixByCol();
  }
  originalNumberRows_ = model->solver()->getNumRows();
}
// Create C++ lines to get to current state
void CbcHeuristicGreedyCover::generateCpp(FILE *fp)
{
  CbcHeuristicGreedyCover other;
  fprintf(fp, "0#include \"CbcHeuristicGreedy.hpp\"\n");
  fprintf(fp, "3  CbcHeuristicGreedyCover heuristicGreedyCover(*cbcModel);\n");
  CbcHeuristic::generateCpp(fp, "heuristicGreedyCover");
  if (algorithm_ != other.algorithm_)
    fprintf(fp, "3  heuristicGreedyCover.setAlgorithm(%d);\n", algorithm_);
  else
    fprintf(fp, "4  heuristicGreedyCover.setAlgorithm(%d);\n", algorithm_);
  if (numberTimes_ != other.numberTimes_)
    fprintf(fp, "3  heuristicGreedyCover.setNumberTimes(%d);\n", numberTimes_);
  else
    fprintf(fp, "4  heuristicGreedyCover.setNumberTimes(%d);\n", numberTimes_);
  fprintf(fp, "3  cbcModel->addHeuristic(&heuristicGreedyCover);\n");
}

// Copy constructor
CbcHeuristicGreedyCover::CbcHeuristicGreedyCover(const CbcHeuristicGreedyCover &rhs)
  : CbcHeuristic(rhs)
  , matrix_(rhs.matrix_)
  , originalNumberRows_(rhs.originalNumberRows_)
  , algorithm_(rhs.algorithm_)
  , numberTimes_(rhs.numberTimes_)
{
}

// Assignment operator
CbcHeuristicGreedyCover &
CbcHeuristicGreedyCover::operator=(const CbcHeuristicGreedyCover &rhs)
{
  if (this != &rhs) {
    CbcHeuristic::operator=(rhs);
    matrix_ = rhs.matrix_;
    originalNumberRows_ = rhs.originalNumberRows_;
    algorithm_ = rhs.algorithm_;
    numberTimes_ = rhs.numberTimes_;
  }
  return *this;
}
// Returns 1 if solution, 0 if not
int CbcHeuristicGreedyCover::solution(double &solutionValue,
  double *betterSolution)
{
  numCouldRun_++;
  if (!model_)
    return 0;
  // See if to do
  if (!when() || (when() == 1 && model_->phase() != 1))
    return 0; // switched off
  if (model_->getNodeCount() > numberTimes_)
    return 0;
#ifdef HEURISTIC_INFORM
  printf("Entering heuristic %s - nRuns %d numCould %d when %d\n",
    heuristicName(), numRuns_, numCouldRun_, when_);
#endif
  // See if at root node
  bool atRoot = model_->getNodeCount() == 0;
  int passNumber = model_->getCurrentPassNumber();
  if (atRoot && passNumber > 1)
    return 0;
  OsiSolverInterface *solver = model_->solver();
  const double *columnLower = solver->getColLower();
  const double *columnUpper = solver->getColUpper();
  // And original upper bounds in case we want to use them
  const double *originalUpper = model_->continuousSolver()->getColUpper();
  // But not if algorithm says so
  if ((algorithm_ % 10) == 0)
    originalUpper = columnUpper;
  const double *rowLower = solver->getRowLower();
  const double *solution = solver->getColSolution();
  const double *objective = solver->getObjCoefficients();
  double integerTolerance = model_->getDblParam(CbcModel::CbcIntegerTolerance);
  double primalTolerance;
  solver->getDblParam(OsiPrimalTolerance, primalTolerance);

  // This is number of rows when matrix was passed in
  int numberRows = originalNumberRows_;
  if (!numberRows)
    return 0; // switched off

  numRuns_++;
  assert(numberRows == matrix_.getNumRows());
  int iRow, iColumn;
  double direction = solver->getObjSense();
  double offset;
  solver->getDblParam(OsiObjOffset, offset);
  double newSolutionValue = -offset;
  int returnCode = 0;

  // Column copy
  const double *element = matrix_.getElements();
  const int *row = matrix_.getIndices();
  const CoinBigIndex *columnStart = matrix_.getVectorStarts();
  const int *columnLength = matrix_.getVectorLengths();

  // Get solution array for heuristic solution
  int numberColumns = solver->getNumCols();
  double *newSolution = new double[numberColumns];
  double *rowActivity = new double[numberRows];
  memset(rowActivity, 0, numberRows * sizeof(double));
  bool allOnes = true;
  // Get rounded down solution
  for (iColumn = 0; iColumn < numberColumns; iColumn++) {
    CoinBigIndex j;
    double value = solution[iColumn];
    if (isHeuristicInteger(solver, iColumn)) {
      // Round down integer
      if (fabs(floor(value + 0.5) - value) < integerTolerance) {
        value = floor(CoinMax(value + 1.0e-3, columnLower[iColumn]));
      } else {
        value = CoinMax(floor(value), columnLower[iColumn]);
      }
    }
    // make sure clean
    value = CoinMin(value, columnUpper[iColumn]);
    value = CoinMax(value, columnLower[iColumn]);
    newSolution[iColumn] = value;
    double cost = direction * objective[iColumn];
    newSolutionValue += value * cost;
    for (j = columnStart[iColumn];
         j < columnStart[iColumn] + columnLength[iColumn]; j++) {
      int iRow = row[j];
      rowActivity[iRow] += value * element[j];
      if (element[j] != 1.0)
        allOnes = false;
    }
  }
  // See if we round up
  bool roundup = ((algorithm_ % 100) != 0);
  if (roundup && allOnes) {
    // Get rounded up solution
    for (iColumn = 0; iColumn < numberColumns; iColumn++) {
      CoinBigIndex j;
      double value = solution[iColumn];
      if (isHeuristicInteger(solver, iColumn)) {
        // but round up if no activity
        if (roundup && value >= 0.499999 && !newSolution[iColumn]) {
          bool choose = true;
          for (j = columnStart[iColumn];
               j < columnStart[iColumn] + columnLength[iColumn]; j++) {
            int iRow = row[j];
            if (rowActivity[iRow]) {
              choose = false;
              break;
            }
          }
          if (choose) {
            newSolution[iColumn] = 1.0;
            double cost = direction * objective[iColumn];
            newSolutionValue += cost;
            for (j = columnStart[iColumn];
                 j < columnStart[iColumn] + columnLength[iColumn]; j++) {
              int iRow = row[j];
              rowActivity[iRow] += 1.0;
            }
          }
        }
      }
    }
  }
  // Get initial list
  int *which = new int[numberColumns];
  for (iColumn = 0; iColumn < numberColumns; iColumn++)
    which[iColumn] = iColumn;
  int numberLook = numberColumns;
  // See if we want to perturb more
  double perturb = ((algorithm_ % 10) == 0) ? 0.1 : 0.25;
  // Keep going round until a solution
  while (true) {
    // Get column with best ratio
    int bestColumn = -1;
    double bestRatio = COIN_DBL_MAX;
    double bestStepSize = 0.0;
    int newNumber = 0;
    for (int jColumn = 0; jColumn < numberLook; jColumn++) {
      int iColumn = which[jColumn];
      CoinBigIndex j;
      double value = newSolution[iColumn];
      double cost = direction * objective[iColumn];
      if (isHeuristicInteger(solver, iColumn)) {
        // use current upper or original upper
        if (value + 0.99 < originalUpper[iColumn]) {
          double sum = 0.0;
          int numberExact = 0;
          for (j = columnStart[iColumn];
               j < columnStart[iColumn] + columnLength[iColumn]; j++) {
            int iRow = row[j];
            double gap = rowLower[iRow] - rowActivity[iRow];
            double elementValue = allOnes ? 1.0 : element[j];
            if (gap > 1.0e-7) {
              sum += CoinMin(elementValue, gap);
              if (fabs(elementValue - gap) < 1.0e-7)
                numberExact++;
            }
          }
          // could bias if exact
          if (sum > 0.0) {
            // add to next time
            which[newNumber++] = iColumn;
            double ratio = (cost / sum) * (1.0 + perturb * randomNumberGenerator_.randomDouble());
            // If at root choose first
            if (atRoot)
              ratio = iColumn;
            if (ratio < bestRatio) {
              bestRatio = ratio;
              bestColumn = iColumn;
              bestStepSize = 1.0;
            }
          }
        }
      } else {
        // continuous
        if (value < columnUpper[iColumn]) {
          // Go through twice - first to get step length
          double step = 1.0e50;
          for (j = columnStart[iColumn];
               j < columnStart[iColumn] + columnLength[iColumn]; j++) {
            int iRow = row[j];
            if (rowActivity[iRow] < rowLower[iRow] - 1.0e-10 && element[j] * step + rowActivity[iRow] >= rowLower[iRow]) {
              step = (rowLower[iRow] - rowActivity[iRow]) / element[j];
              ;
            }
          }
          // now ratio
          if (step < 1.0e50) {
            // add to next time
            which[newNumber++] = iColumn;
            assert(step > 0.0);
            double sum = 0.0;
            for (j = columnStart[iColumn];
                 j < columnStart[iColumn] + columnLength[iColumn]; j++) {
              int iRow = row[j];
              double newActivity = element[j] * step + rowActivity[iRow];
              if (rowActivity[iRow] < rowLower[iRow] - 1.0e-10 && newActivity >= rowLower[iRow] - 1.0e-12) {
                sum += element[j];
              }
            }
            assert(sum > 0.0);
            double ratio = (cost / sum) * (1.0 + perturb * randomNumberGenerator_.randomDouble());
            if (ratio < bestRatio) {
              bestRatio = ratio;
              bestColumn = iColumn;
              bestStepSize = step;
            }
          }
        }
      }
    }
    if (bestColumn < 0)
      break; // we have finished
    // Increase chosen column
    newSolution[bestColumn] += bestStepSize;
    double cost = direction * objective[bestColumn];
    newSolutionValue += bestStepSize * cost;
    for (CoinBigIndex j = columnStart[bestColumn];
         j < columnStart[bestColumn] + columnLength[bestColumn]; j++) {
      int iRow = row[j];
      rowActivity[iRow] += bestStepSize * element[j];
    }
  }
  delete[] which;
  if (newSolutionValue < solutionValue) {
    // check feasible
    memset(rowActivity, 0, numberRows * sizeof(double));
    for (iColumn = 0; iColumn < numberColumns; iColumn++) {
      CoinBigIndex j;
      double value = newSolution[iColumn];
      if (value) {
        for (j = columnStart[iColumn];
             j < columnStart[iColumn] + columnLength[iColumn]; j++) {
          int iRow = row[j];
          rowActivity[iRow] += value * element[j];
        }
      }
    }
    // check was approximately feasible
    bool feasible = true;
    for (iRow = 0; iRow < numberRows; iRow++) {
      if (rowActivity[iRow] < rowLower[iRow]) {
        if (rowActivity[iRow] < rowLower[iRow] - 10.0 * primalTolerance)
          feasible = false;
      }
    }
    if (feasible) {
      // new solution
      memcpy(betterSolution, newSolution, numberColumns * sizeof(double));
      solutionValue = newSolutionValue;
      //printf("** Solution of %g found by rounding\n",newSolutionValue);
      returnCode = 1;
    } else {
      // Can easily happen
      //printf("Debug CbcHeuristicGreedyCover giving bad solution\n");
    }
  }
  delete[] newSolution;
  delete[] rowActivity;
  return returnCode;
}
// update model
void CbcHeuristicGreedyCover::setModel(CbcModel *model)
{
  gutsOfConstructor(model);
  validate();
}
// Resets stuff if model changes
void CbcHeuristicGreedyCover::resetModel(CbcModel *model)
{
  gutsOfConstructor(model);
}
// Validate model i.e. sets when_ to 0 if necessary (may be NULL)
void CbcHeuristicGreedyCover::validate()
{
  if (model_ && when() < 10) {
    if (model_->numberIntegers() != model_->numberObjects() && (model_->numberObjects() || (model_->specialOptions() & 1024) == 0)) {
      int numberOdd = 0;
      for (int i = 0; i < model_->numberObjects(); i++) {
        if (!model_->object(i)->canDoHeuristics())
          numberOdd++;
      }
      if (numberOdd)
        setWhen(0);
    }
    // Only works if costs positive, coefficients positive and all rows G
    OsiSolverInterface *solver = model_->solver();
    const double *columnLower = solver->getColLower();
    const double *rowUpper = solver->getRowUpper();
    const double *objective = solver->getObjCoefficients();
    double direction = solver->getObjSense();

    int numberRows = solver->getNumRows();
    int numberColumns = solver->getNumCols();
    // Column copy
    matrix_.setDimensions(numberRows, numberColumns);
    const double *element = matrix_.getElements();
    const CoinBigIndex *columnStart = matrix_.getVectorStarts();
    const int *columnLength = matrix_.getVectorLengths();
    bool good = true;
    for (int iRow = 0; iRow < numberRows; iRow++) {
      if (rowUpper[iRow] < 1.0e30)
        good = false;
    }
    for (int iColumn = 0; iColumn < numberColumns; iColumn++) {
      if (objective[iColumn] * direction < 0.0)
        good = false;
      if (columnLower[iColumn] < 0.0)
        good = false;
      CoinBigIndex j;
      for (j = columnStart[iColumn];
           j < columnStart[iColumn] + columnLength[iColumn]; j++) {
        if (element[j] < 0.0)
          good = false;
      }
    }
    if (!good)
      setWhen(0); // switch off
  }
}
// Default Constructor
CbcHeuristicGreedyEquality::CbcHeuristicGreedyEquality()
  : CbcHeuristic()
{
  // matrix  will automatically be empty
  fraction_ = 1.0; // no branch and bound
  originalNumberRows_ = 0;
  algorithm_ = 0;
  numberTimes_ = 100;
  whereFrom_ = 1;
}

// Constructor from model
CbcHeuristicGreedyEquality::CbcHeuristicGreedyEquality(CbcModel &model)
  : CbcHeuristic(model)
{
  // Get a copy of original matrix
  gutsOfConstructor(&model);
  fraction_ = 1.0; // no branch and bound
  algorithm_ = 0;
  numberTimes_ = 100;
  whereFrom_ = 1;
}

// Destructor
CbcHeuristicGreedyEquality::~CbcHeuristicGreedyEquality()
{
}

// Clone
CbcHeuristic *
CbcHeuristicGreedyEquality::clone() const
{
  return new CbcHeuristicGreedyEquality(*this);
}
// Guts of constructor from a CbcModel
void CbcHeuristicGreedyEquality::gutsOfConstructor(CbcModel *model)
{
  model_ = model;
  // Get a copy of original matrix
  assert(model->solver());
  if (model->solver()->getNumRows()) {
    matrix_ = *model->solver()->getMatrixByCol();
  }
  originalNumberRows_ = model->solver()->getNumRows();
}
// Create C++ lines to get to current state
void CbcHeuristicGreedyEquality::generateCpp(FILE *fp)
{
  CbcHeuristicGreedyEquality other;
  fprintf(fp, "0#include \"CbcHeuristicGreedy.hpp\"\n");
  fprintf(fp, "3  CbcHeuristicGreedyEquality heuristicGreedyEquality(*cbcModel);\n");
  CbcHeuristic::generateCpp(fp, "heuristicGreedyEquality");
  if (algorithm_ != other.algorithm_)
    fprintf(fp, "3  heuristicGreedyEquality.setAlgorithm(%d);\n", algorithm_);
  else
    fprintf(fp, "4  heuristicGreedyEquality.setAlgorithm(%d);\n", algorithm_);
  if (fraction_ != other.fraction_)
    fprintf(fp, "3  heuristicGreedyEquality.setFraction(%g);\n", fraction_);
  else
    fprintf(fp, "4  heuristicGreedyEquality.setFraction(%g);\n", fraction_);
  if (numberTimes_ != other.numberTimes_)
    fprintf(fp, "3  heuristicGreedyEquality.setNumberTimes(%d);\n", numberTimes_);
  else
    fprintf(fp, "4  heuristicGreedyEquality.setNumberTimes(%d);\n", numberTimes_);
  fprintf(fp, "3  cbcModel->addHeuristic(&heuristicGreedyEquality);\n");
}

// Copy constructor
CbcHeuristicGreedyEquality::CbcHeuristicGreedyEquality(const CbcHeuristicGreedyEquality &rhs)
  : CbcHeuristic(rhs)
  , matrix_(rhs.matrix_)
  , fraction_(rhs.fraction_)
  , originalNumberRows_(rhs.originalNumberRows_)
  , algorithm_(rhs.algorithm_)
  , numberTimes_(rhs.numberTimes_)
{
}

// Assignment operator
CbcHeuristicGreedyEquality &
CbcHeuristicGreedyEquality::operator=(const CbcHeuristicGreedyEquality &rhs)
{
  if (this != &rhs) {
    CbcHeuristic::operator=(rhs);
    matrix_ = rhs.matrix_;
    fraction_ = rhs.fraction_;
    originalNumberRows_ = rhs.originalNumberRows_;
    algorithm_ = rhs.algorithm_;
    numberTimes_ = rhs.numberTimes_;
  }
  return *this;
}
// Returns 1 if solution, 0 if not
int CbcHeuristicGreedyEquality::solution(double &solutionValue,
  double *betterSolution)
{
  numCouldRun_++;
  if (!model_)
    return 0;
  // See if to do
  if (!when() || (when() == 1 && model_->phase() != 1))
    return 0; // switched off
  if (model_->getNodeCount() > numberTimes_)
    return 0;
  // See if at root node
  bool atRoot = model_->getNodeCount() == 0;
  int passNumber = model_->getCurrentPassNumber();
  if (atRoot && passNumber > 1)
    return 0;
  OsiSolverInterface *solver = model_->solver();
  const double *columnLower = solver->getColLower();
  const double *columnUpper = solver->getColUpper();
  // And original upper bounds in case we want to use them
  const double *originalUpper = model_->continuousSolver()->getColUpper();
  // But not if algorithm says so
  if ((algorithm_ % 10) == 0)
    originalUpper = columnUpper;
  const double *rowLower = solver->getRowLower();
  const double *rowUpper = solver->getRowUpper();
  const double *solution = solver->getColSolution();
  const double *objective = solver->getObjCoefficients();
  double integerTolerance = model_->getDblParam(CbcModel::CbcIntegerTolerance);
  double primalTolerance;
  solver->getDblParam(OsiPrimalTolerance, primalTolerance);

  // This is number of rows when matrix was passed in
  int numberRows = originalNumberRows_;
  if (!numberRows)
    return 0; // switched off
  numRuns_++;

  assert(numberRows == matrix_.getNumRows());
  int iRow, iColumn;
  double direction = solver->getObjSense();
  double offset;
  solver->getDblParam(OsiObjOffset, offset);
  double newSolutionValue = -offset;
  int returnCode = 0;

  // Column copy
  const double *element = matrix_.getElements();
  const int *row = matrix_.getIndices();
  const CoinBigIndex *columnStart = matrix_.getVectorStarts();
  const int *columnLength = matrix_.getVectorLengths();

  // Get solution array for heuristic solution
  int numberColumns = solver->getNumCols();
  double *newSolution = new double[numberColumns];
  double *rowActivity = new double[numberRows];
  memset(rowActivity, 0, numberRows * sizeof(double));
  double rhsNeeded = 0;
  for (iRow = 0; iRow < numberRows; iRow++)
    rhsNeeded += rowUpper[iRow];
  rhsNeeded *= fraction_;
  bool allOnes = true;
  // Get rounded down solution
  for (iColumn = 0; iColumn < numberColumns; iColumn++) {
    CoinBigIndex j;
    double value = solution[iColumn];
    if (isHeuristicInteger(solver, iColumn)) {
      // Round down integer
      if (fabs(floor(value + 0.5) - value) < integerTolerance) {
        value = floor(CoinMax(value + 1.0e-3, columnLower[iColumn]));
      } else {
        value = CoinMax(floor(value), columnLower[iColumn]);
      }
    }
    // make sure clean
    value = CoinMin(value, columnUpper[iColumn]);
    value = CoinMax(value, columnLower[iColumn]);
    newSolution[iColumn] = value;
    double cost = direction * objective[iColumn];
    newSolutionValue += value * cost;
    for (j = columnStart[iColumn];
         j < columnStart[iColumn] + columnLength[iColumn]; j++) {
      int iRow = row[j];
      rowActivity[iRow] += value * element[j];
      rhsNeeded -= value * element[j];
      if (element[j] != 1.0)
        allOnes = false;
    }
  }
  // See if we round up
  bool roundup = ((algorithm_ % 100) != 0);
  if (roundup && allOnes) {
    // Get rounded up solution
    for (iColumn = 0; iColumn < numberColumns; iColumn++) {
      CoinBigIndex j;
      double value = solution[iColumn];
      if (isHeuristicInteger(solver, iColumn)) {
        // but round up if no activity
        if (roundup && value >= 0.6 && !newSolution[iColumn]) {
          bool choose = true;
          for (j = columnStart[iColumn];
               j < columnStart[iColumn] + columnLength[iColumn]; j++) {
            int iRow = row[j];
            if (rowActivity[iRow]) {
              choose = false;
              break;
            }
          }
          if (choose) {
            newSolution[iColumn] = 1.0;
            double cost = direction * objective[iColumn];
            newSolutionValue += cost;
            for (j = columnStart[iColumn];
                 j < columnStart[iColumn] + columnLength[iColumn]; j++) {
              int iRow = row[j];
              rowActivity[iRow] += 1.0;
              rhsNeeded -= 1.0;
            }
          }
        }
      }
    }
  }
  // Get initial list
  int *which = new int[numberColumns];
  for (iColumn = 0; iColumn < numberColumns; iColumn++)
    which[iColumn] = iColumn;
  int numberLook = numberColumns;
  // See if we want to perturb more
  double perturb = ((algorithm_ % 10) == 0) ? 0.1 : 0.25;
  // Keep going round until a solution
  while (true) {
    // Get column with best ratio
    int bestColumn = -1;
    double bestRatio = COIN_DBL_MAX;
    double bestStepSize = 0.0;
    int newNumber = 0;
    for (int jColumn = 0; jColumn < numberLook; jColumn++) {
      int iColumn = which[jColumn];
      CoinBigIndex j;
      double value = newSolution[iColumn];
      double cost = direction * objective[iColumn];
      if (isHeuristicInteger(solver, iColumn)) {
        // use current upper or original upper
        if (value + 0.9999 < originalUpper[iColumn]) {
          double movement = 1.0;
          double sum = 0.0;
          for (j = columnStart[iColumn];
               j < columnStart[iColumn] + columnLength[iColumn]; j++) {
            int iRow = row[j];
            double gap = rowUpper[iRow] - rowActivity[iRow];
            double elementValue = allOnes ? 1.0 : element[j];
            sum += elementValue;
            if (movement * elementValue > gap) {
              movement = gap / elementValue;
            }
          }
          if (movement > 0.999999) {
            // add to next time
            which[newNumber++] = iColumn;
            double ratio = (cost / sum) * (1.0 + perturb * randomNumberGenerator_.randomDouble());
            // If at root
            if (atRoot) {
              if (fraction_ == 1.0)
                ratio = iColumn; // choose first
              else
                ratio = -solution[iColumn]; // choose largest
            }
            if (ratio < bestRatio) {
              bestRatio = ratio;
              bestColumn = iColumn;
              bestStepSize = 1.0;
            }
          }
        }
      } else {
        // continuous
        if (value < columnUpper[iColumn]) {
          double movement = 1.0e50;
          double sum = 0.0;
          for (j = columnStart[iColumn];
               j < columnStart[iColumn] + columnLength[iColumn]; j++) {
            int iRow = row[j];
            if (element[j] * movement + rowActivity[iRow] > rowUpper[iRow]) {
              movement = (rowUpper[iRow] - rowActivity[iRow]) / element[j];
              ;
            }
            sum += element[j];
          }
          // now ratio
          if (movement > 1.0e-7) {
            // add to next time
            which[newNumber++] = iColumn;
            double ratio = (cost / sum) * (1.0 + perturb * randomNumberGenerator_.randomDouble());
            if (ratio < bestRatio) {
              bestRatio = ratio;
              bestColumn = iColumn;
              bestStepSize = movement;
            }
          }
        }
      }
    }
    if (bestColumn < 0)
      break; // we have finished
    // Increase chosen column
    newSolution[bestColumn] += bestStepSize;
    double cost = direction * objective[bestColumn];
    newSolutionValue += bestStepSize * cost;
    for (CoinBigIndex j = columnStart[bestColumn];
         j < columnStart[bestColumn] + columnLength[bestColumn]; j++) {
      int iRow = row[j];
      rowActivity[iRow] += bestStepSize * element[j];
      rhsNeeded -= bestStepSize * element[j];
    }
    if (rhsNeeded < 1.0e-8)
      break;
  }
  delete[] which;
  if (fraction_ < 1.0 && rhsNeeded < 1.0e-8 && newSolutionValue < solutionValue) {
    // do branch and cut
    // fix all nonzero
    OsiSolverInterface *newSolver = model_->continuousSolver()->clone();
    for (iColumn = 0; iColumn < numberColumns; iColumn++) {
      if (isHeuristicInteger(newSolver, iColumn))
        newSolver->setColLower(iColumn, newSolution[iColumn]);
    }
    int returnCode = smallBranchAndBound(newSolver, 200, newSolution, newSolutionValue,
      solutionValue, "CbcHeuristicGreedy");
    if (returnCode < 0)
      returnCode = 0; // returned on size
    if ((returnCode & 2) != 0) {
      // could add cut
      returnCode &= ~2;
    }
    rhsNeeded = 1.0 - returnCode;
    delete newSolver;
  }
  if (newSolutionValue < solutionValue && rhsNeeded < 1.0e-8) {
    // check feasible
    memset(rowActivity, 0, numberRows * sizeof(double));
    for (iColumn = 0; iColumn < numberColumns; iColumn++) {
      CoinBigIndex j;
      double value = newSolution[iColumn];
      if (value) {
        for (j = columnStart[iColumn];
             j < columnStart[iColumn] + columnLength[iColumn]; j++) {
          int iRow = row[j];
          rowActivity[iRow] += value * element[j];
        }
      }
    }
    // check was approximately feasible
    bool feasible = true;
    for (iRow = 0; iRow < numberRows; iRow++) {
      if (rowActivity[iRow] < rowLower[iRow]) {
        if (rowActivity[iRow] < rowLower[iRow] - 10.0 * primalTolerance)
          feasible = false;
      }
    }
    if (feasible) {
      // new solution
      memcpy(betterSolution, newSolution, numberColumns * sizeof(double));
      solutionValue = newSolutionValue;
      returnCode = 1;
    }
  }
  delete[] newSolution;
  delete[] rowActivity;
  if (atRoot && fraction_ == 1.0) {
    // try quick search
    fraction_ = 0.4;
    int newCode = this->solution(solutionValue, betterSolution);
    if (newCode)
      returnCode = 1;
    fraction_ = 1.0;
  }
  return returnCode;
}
// update model
void CbcHeuristicGreedyEquality::setModel(CbcModel *model)
{
  gutsOfConstructor(model);
  validate();
}
// Resets stuff if model changes
void CbcHeuristicGreedyEquality::resetModel(CbcModel *model)
{
  gutsOfConstructor(model);
}
// Validate model i.e. sets when_ to 0 if necessary (may be NULL)
void CbcHeuristicGreedyEquality::validate()
{
  if (model_ && when() < 10) {
    if (model_->numberIntegers() != model_->numberObjects())
      setWhen(0);
    // Only works if costs positive, coefficients positive and all rows E or L
    // And if values are integer
    OsiSolverInterface *solver = model_->solver();
    const double *columnLower = solver->getColLower();
    const double *rowUpper = solver->getRowUpper();
    const double *rowLower = solver->getRowLower();
    const double *objective = solver->getObjCoefficients();
    double direction = solver->getObjSense();

    int numberRows = solver->getNumRows();
    int numberColumns = solver->getNumCols();
    matrix_.setDimensions(numberRows, numberColumns);
    // Column copy
    const double *element = matrix_.getElements();
    const CoinBigIndex *columnStart = matrix_.getVectorStarts();
    const int *columnLength = matrix_.getVectorLengths();
    bool good = true;
    for (int iRow = 0; iRow < numberRows; iRow++) {
      if (rowUpper[iRow] > 1.0e30)
        good = false;
      if (rowLower[iRow] > 0.0 && rowLower[iRow] != rowUpper[iRow])
        good = false;
      if (floor(rowUpper[iRow] + 0.5) != rowUpper[iRow])
        good = false;
    }
    for (int iColumn = 0; iColumn < numberColumns; iColumn++) {
      if (objective[iColumn] * direction < 0.0)
        good = false;
      if (columnLower[iColumn] < 0.0)
        good = false;
      CoinBigIndex j;
      for (j = columnStart[iColumn];
           j < columnStart[iColumn] + columnLength[iColumn]; j++) {
        if (element[j] < 0.0)
          good = false;
        if (floor(element[j] + 0.5) != element[j])
          good = false;
      }
    }
    if (!good)
      setWhen(0); // switch off
  }
}

// Default Constructor
CbcHeuristicGreedySOS::CbcHeuristicGreedySOS()
  : CbcHeuristic()
{
  originalRhs_ = NULL;
  // matrix  will automatically be empty
  originalNumberRows_ = 0;
  algorithm_ = 0;
  numberTimes_ = 100;
}

// Constructor from model
CbcHeuristicGreedySOS::CbcHeuristicGreedySOS(CbcModel &model)
  : CbcHeuristic(model)
{
  gutsOfConstructor(&model);
  algorithm_ = 2;
  numberTimes_ = 100;
  whereFrom_ = 1;
}

// Destructor
CbcHeuristicGreedySOS::~CbcHeuristicGreedySOS()
{
  delete[] originalRhs_;
}

// Clone
CbcHeuristic *
CbcHeuristicGreedySOS::clone() const
{
  return new CbcHeuristicGreedySOS(*this);
}
// Guts of constructor from a CbcModel
void CbcHeuristicGreedySOS::gutsOfConstructor(CbcModel *model)
{
  model_ = model;
  // Get a copy of original matrix
  assert(model->solver());
  if (model->solver()->getNumRows()) {
    matrix_ = *model->solver()->getMatrixByCol();
  }
  originalNumberRows_ = model->solver()->getNumRows();
  originalRhs_ = new double[originalNumberRows_];
}
// Create C++ lines to get to current state
void CbcHeuristicGreedySOS::generateCpp(FILE *fp)
{
  CbcHeuristicGreedySOS other;
  fprintf(fp, "0#include \"CbcHeuristicGreedy.hpp\"\n");
  fprintf(fp, "3  CbcHeuristicGreedySOS heuristicGreedySOS(*cbcModel);\n");
  CbcHeuristic::generateCpp(fp, "heuristicGreedySOS");
  if (algorithm_ != other.algorithm_)
    fprintf(fp, "3  heuristicGreedySOS.setAlgorithm(%d);\n", algorithm_);
  else
    fprintf(fp, "4  heuristicGreedySOS.setAlgorithm(%d);\n", algorithm_);
  if (numberTimes_ != other.numberTimes_)
    fprintf(fp, "3  heuristicGreedySOS.setNumberTimes(%d);\n", numberTimes_);
  else
    fprintf(fp, "4  heuristicGreedySOS.setNumberTimes(%d);\n", numberTimes_);
  fprintf(fp, "3  cbcModel->addHeuristic(&heuristicGreedySOS);\n");
}

// Copy constructor
CbcHeuristicGreedySOS::CbcHeuristicGreedySOS(const CbcHeuristicGreedySOS &rhs)
  : CbcHeuristic(rhs)
  , matrix_(rhs.matrix_)
  , originalNumberRows_(rhs.originalNumberRows_)
  , algorithm_(rhs.algorithm_)
  , numberTimes_(rhs.numberTimes_)
{
  originalRhs_ = CoinCopyOfArray(rhs.originalRhs_, originalNumberRows_);
}

// Assignment operator
CbcHeuristicGreedySOS &
CbcHeuristicGreedySOS::operator=(const CbcHeuristicGreedySOS &rhs)
{
  if (this != &rhs) {
    CbcHeuristic::operator=(rhs);
    matrix_ = rhs.matrix_;
    originalNumberRows_ = rhs.originalNumberRows_;
    algorithm_ = rhs.algorithm_;
    numberTimes_ = rhs.numberTimes_;
    delete[] originalRhs_;
    originalRhs_ = CoinCopyOfArray(rhs.originalRhs_, originalNumberRows_);
  }
  return *this;
}
// Returns 1 if solution, 0 if not
int CbcHeuristicGreedySOS::solution(double &solutionValue,
  double *betterSolution)
{
  numCouldRun_++;
  if (!model_)
    return 0;
  // See if to do
  if (!when() || (when() == 1 && model_->phase() != 1))
    return 0; // switched off
  if (model_->getNodeCount() > numberTimes_)
    return 0;
  // See if at root node
  bool atRoot = model_->getNodeCount() == 0;
  int passNumber = model_->getCurrentPassNumber();
  if (atRoot && passNumber > 1)
    return 0;
  OsiSolverInterface *solver = model_->solver();
  int numberColumns = solver->getNumCols();
  // This is number of rows when matrix was passed in
  int numberRows = originalNumberRows_;
  if (!numberRows)
    return 0; // switched off

  const double *columnLower = solver->getColLower();
  const double *columnUpper = solver->getColUpper();
  // modified rhs
  double *rhs = CoinCopyOfArray(originalRhs_, numberRows);
  // Column copy
  const double *element = matrix_.getElements();
  const int *row = matrix_.getIndices();
  const CoinBigIndex *columnStart = matrix_.getVectorStarts();
  const int *columnLength = matrix_.getVectorLengths();
  int *sosRow = new int[numberColumns];
  char *sos = new char[numberRows];
  memset(sos, 'a', numberRows);
  int nonSOS = 0;
  // If bit set then use current
  if ((algorithm_ & 1) != 0) {
    const CoinPackedMatrix *matrix = solver->getMatrixByCol();
    element = matrix->getElements();
    row = matrix->getIndices();
    columnStart = matrix->getVectorStarts();
    columnLength = matrix->getVectorLengths();
    //rhs = new double [numberRows];
    const double *rowLower = solver->getRowLower();
    const double *rowUpper = solver->getRowUpper();
    bool good = true;
    for (int iRow = 0; iRow < numberRows; iRow++) {
      if (rowLower[iRow] == 1.0 && rowUpper[iRow] == 1.0) {
        // SOS
        rhs[iRow] = -1.0;
      } else if (rowLower[iRow] > 0.0 && rowUpper[iRow] < 1.0e10) {
        good = false;
      } else if (rowUpper[iRow] < 0.0) {
        good = false;
      } else if (rowUpper[iRow] < 1.0e10) {
        rhs[iRow] = rowUpper[iRow];
        if (rhs[iRow] < 0)
          sos[iRow] = 0; // can't be SOS
      } else {
        rhs[iRow] = rowLower[iRow];
        if (rhs[iRow] < 0)
          sos[iRow] = 0; // can't be SOS
      }
    }
    for (int iColumn = 0; iColumn < numberColumns; iColumn++) {
      if (!columnLength[iColumn])
        continue;
      if (columnLower[iColumn] < 0.0 || columnUpper[iColumn] > 1.0)
        good = false;
      CoinBigIndex j;
      int nSOS = 0;
      int iSOS = -1;
      if (!isHeuristicInteger(solver, iColumn))
        good = false;
      for (j = columnStart[iColumn];
           j < columnStart[iColumn] + columnLength[iColumn]; j++) {
        if (element[j] < 0.0)
          good = false;
        int iRow = row[j];
        if (rhs[iRow] == -1.0 && sos[iRow] == 'a') {
          if (element[j] != 1.0)
            good = false;
          iSOS = iRow;
          nSOS++;
        }
      }
      if (nSOS > 1)
        good = false;
      else if (!nSOS)
        nonSOS++;
      sosRow[iColumn] = iSOS;
    }
    if (!good) {
      delete[] sosRow;
      delete[] rhs;
      delete[] sos;
      setWhen(0); // switch off
      return 0;
    }
  } else {
    abort(); // not allowed yet
  }
  const double *solution = solver->getColSolution();
  const double *objective = solver->getObjCoefficients();
  const double *rowLower = solver->getRowLower();
  const double *rowUpper = solver->getRowUpper();
  double integerTolerance = model_->getDblParam(CbcModel::CbcIntegerTolerance);
  double primalTolerance;
  solver->getDblParam(OsiPrimalTolerance, primalTolerance);

  numRuns_++;
  assert(numberRows == matrix_.getNumRows());
  // set up linked list for sets
  int *firstGub = new int[numberRows];
  int *nextGub = new int[numberColumns];
  int iRow, iColumn;
  double direction = solver->getObjSense();
  double *slackCost = new double[numberRows];
  double *modifiedCost = CoinCopyOfArray(objective, numberColumns);
  for (int iRow = 0; iRow < numberRows; iRow++) {
    slackCost[iRow] = 1.0e30;
    firstGub[iRow] = -1;
  }
  // Take off cost of gub slack
  for (int iColumn = 0; iColumn < numberColumns; iColumn++) {
    nextGub[iColumn] = -1;
    int iRow = sosRow[iColumn];
    if (columnLength[iColumn] == 1 && iRow >= 0) {
      // SOS slack
      double cost = direction * objective[iColumn];
      assert(rhs[iRow] < 0.0);
      slackCost[iRow] = CoinMin(slackCost[iRow], cost);
    }
  }
  double offset2 = 0.0;
  for (int iRow = 0; iRow < numberRows; iRow++) {
    if (sos[iRow] == 'a') {
      // row is possible
      sos[iRow] = 0;
      if (rhs[iRow] < 0.0) {
        sos[iRow] = 1;
        rhs[iRow] = 1.0;
      } else if (rhs[iRow] != rowUpper[iRow]) {
        // G row
        sos[iRow] = -1;
      }
      if (slackCost[iRow] == 1.0e30) {
        slackCost[iRow] = 0.0;
      } else {
        offset2 += slackCost[iRow];
        sos[iRow] = 2;
      }
    }
  }
  for (int iColumn = 0; iColumn < numberColumns; iColumn++) {
    double cost = direction * modifiedCost[iColumn];
    CoinBigIndex j;
    for (j = columnStart[iColumn];
         j < columnStart[iColumn] + columnLength[iColumn]; j++) {
      int iRow = row[j];
      if (sos[iRow] > 0) {
        cost -= slackCost[iRow];
        if (firstGub[iRow] < 0) {
          firstGub[iRow] = iColumn;
        } else {
          int jColumn = firstGub[iRow];
          while (nextGub[jColumn] >= 0)
            jColumn = nextGub[jColumn];
          nextGub[jColumn] = iColumn;
        }
        // Only in one sos
        break;
      }
    }
    modifiedCost[iColumn] = cost;
  }
  delete[] slackCost;
  double offset;
  solver->getDblParam(OsiObjOffset, offset);
  double newSolutionValue = -offset + offset2;
  int returnCode = 0;

  // Get solution array for heuristic solution
  double *newSolution = new double[numberColumns];
  double *rowActivity = new double[numberRows];
  double *contribution = new double[numberColumns];
  int *which = new int[numberColumns];
  double *newSolution0 = new double[numberColumns];
  if ((algorithm_ & (2 | 4)) == 0) {
    // get solution as small as possible
    for (iColumn = 0; iColumn < numberColumns; iColumn++)
      newSolution0[iColumn] = columnLower[iColumn];
  } else {
    // Get rounded down solution
    for (iColumn = 0; iColumn < numberColumns; iColumn++) {
      double value = solution[iColumn];
      // Round down integer
      if (fabs(floor(value + 0.5) - value) < integerTolerance) {
        value = floor(CoinMax(value + 1.0e-3, columnLower[iColumn]));
      } else {
        value = CoinMax(floor(value), columnLower[iColumn]);
      }
      // make sure clean
      value = CoinMin(value, columnUpper[iColumn]);
      value = CoinMax(value, columnLower[iColumn]);
      newSolution0[iColumn] = value;
    }
  }
  double *rowWeight = new double[numberRows];
  for (int i = 0; i < numberRows; i++)
    rowWeight[i] = 1.0;
  double costBias = 0.0;
  int nPass = ((algorithm_ & 4) != 0) ? 1 : 10;
  for (int iPass = 0; iPass < nPass; iPass++) {
    newSolutionValue = -offset + offset2;
    memcpy(newSolution, newSolution0, numberColumns * sizeof(double));
    // get row activity
    memset(rowActivity, 0, numberRows * sizeof(double));
    for (iColumn = 0; iColumn < numberColumns; iColumn++) {
      CoinBigIndex j;
      double value = newSolution[iColumn];
      for (j = columnStart[iColumn];
           j < columnStart[iColumn] + columnLength[iColumn]; j++) {
        int iRow = row[j];
        rowActivity[iRow] += value * element[j];
      }
    }
    if (!rowWeight) {
      rowWeight = CoinCopyOfArray(rowActivity, numberRows);
    }
    for (iColumn = 0; iColumn < numberColumns; iColumn++) {
      CoinBigIndex j;
      double value = newSolution[iColumn];
      double cost = modifiedCost[iColumn];
      double forSort = 1.0e-24;
      bool hasSlack = false;
      bool willFit = true;
      bool gRow = false;
      newSolutionValue += value * cost;
      cost += 1.0e-12;
      for (j = columnStart[iColumn];
           j < columnStart[iColumn] + columnLength[iColumn]; j++) {
        int iRow = row[j];
        int type = sos[iRow];
        double gap = rhs[iRow] - rowActivity[iRow] + 1.0e-8;
        switch (type) {
        case -1:
          // G row
          gRow = true;
#if 0
	    if (rhs[iRow]>rowWeight[iRow]||(algorithm_&(2|4))!=0)
	      forSort += element[j];
	    else
	      forSort += 0.1*element[j];
#else
          forSort += rowWeight[iRow] * element[j];
#endif
          break;
        case 0:
          // L row
          if (gap < element[j]) {
            willFit = false;
          } else {
            forSort += element[j];
          }
          break;
        case 1:
          // SOS without slack
          if (gap < element[j]) {
            willFit = false;
          }
          break;
        case 2:
          // SOS with slack
          hasSlack = true;
          if (gap < element[j]) {
            willFit = false;
          }
          break;
        }
      }
      bool isSlack = hasSlack && (columnLength[iColumn] == 1);
      if (forSort < 1.0e-24)
        forSort = 1.0e-12;
      if ((algorithm_ & 4) != 0 && forSort > 1.0e-24)
        forSort = 1.0;
      // Use smallest cost if will fit
      if (willFit && (hasSlack || gRow) && value == 0.0 && columnUpper[iColumn]) {
        if (hasSlack && !gRow) {
          if (cost > 1.0e-12) {
            forSort = 2.0e30;
          } else if (cost == 1.0e-12) {
            if (!isSlack)
              forSort = 1.0e29;
            else
              forSort = 1.0e28;
          } else {
            forSort = cost / forSort;
          }
        } else {
          if (!gRow || true)
            forSort = (cost + costBias) / forSort;
          else
            forSort = 1.0e-12 / forSort;
        }
      } else {
        // put at end
        forSort = 1.0e30;
      }
      which[iColumn] = iColumn;
      contribution[iColumn] = forSort;
    }
    CoinSort_2(contribution, contribution + numberColumns, which);
    // Go through columns
    int nAdded = 0;
    int nSlacks = 0;
    for (int jColumn = 0; jColumn < numberColumns; jColumn++) {
      if (contribution[jColumn] >= 1.0e30)
        break;
      int iColumn = which[jColumn];
      double value = newSolution[iColumn];
      if (value)
        continue;
      bool possible = true;
      CoinBigIndex j;
      for (j = columnStart[iColumn];
           j < columnStart[iColumn] + columnLength[iColumn]; j++) {
        int iRow = row[j];
        if (sos[iRow] > 0 && rowActivity[iRow]) {
          possible = false;
        } else {
          double gap = rhs[iRow] - rowActivity[iRow] + 1.0e-8;
          if (gap < element[j] && sos[iRow] >= 0)
            possible = false;
        }
      }
      if (possible) {
        //#define REPORT 1
#ifdef REPORT
        if ((nAdded % 1000) == 0) {
          double gap = 0.0;
          for (int i = 0; i < numberRows; i++) {
            if (rowUpper[i] > 1.0e20)
              gap += CoinMax(rowLower[i] - rowActivity[i], 0.0);
          }
          if (gap)
            printf("after %d added gap %g - %d slacks\n",
              nAdded, gap, nSlacks);
        }
#endif
        nAdded++;
        if (columnLength[iColumn] == 1)
          nSlacks++;
        // Increase chosen column
        newSolution[iColumn] = 1.0;
        double cost = modifiedCost[iColumn];
        newSolutionValue += cost;
        for (CoinBigIndex j = columnStart[iColumn];
             j < columnStart[iColumn] + columnLength[iColumn]; j++) {
          int iRow = row[j];
          rowActivity[iRow] += element[j];
        }
      }
    }
#ifdef REPORT
    {
      double under = 0.0;
      double over = 0.0;
      double gap = 0.0;
      int nUnder = 0;
      int nOver = 0;
      int nGap = 0;
      for (iRow = 0; iRow < numberRows; iRow++) {
        if (rowActivity[iRow] < rowLower[iRow] - 10.0 * primalTolerance) {
          double value = rowLower[iRow] - rowActivity[iRow];
#if REPORT > 1
          printf("below on %d is %g - activity %g lower %g\n",
            iRow, value, rowActivity[iRow], rowLower[iRow]);
#endif
          under += value;
          nUnder++;
        } else if (rowActivity[iRow] > rowUpper[iRow] + 10.0 * primalTolerance) {
          double value = rowActivity[iRow] - rowUpper[iRow];
#if REPORT > 1
          printf("above on %d is %g - activity %g upper %g\n",
            iRow, value, rowActivity[iRow], rowUpper[iRow]);
#endif
          over += value;
          nOver++;
        } else {
          double value = rowActivity[iRow] - rowLower[iRow];
          if (value && value < 1.0e20) {
#if REPORT > 1
            printf("gap on %d is %g - activity %g lower %g\n",
              iRow, value, rowActivity[iRow], rowLower[iRow]);
#endif
            gap += value;
            nGap++;
          }
        }
      }
      printf("final under %g (%d) - over %g (%d) - free %g (%d) - %d added - solvalue %g\n",
        under, nUnder, over, nOver, gap, nGap, nAdded, newSolutionValue);
    }
#endif
    double gap = 0.0;
    double over = 0.0;
    int nL = 0;
    int nG = 0;
    int nUnder = 0;
    for (iRow = 0; iRow < numberRows; iRow++) {
      if (rowLower[iRow] < -1.0e20)
        nL++;
      if (rowUpper[iRow] > 1.0e20)
        nG++;
      if (rowActivity[iRow] < rowLower[iRow] - 10.0 * primalTolerance) {
        gap += rowLower[iRow] - rowActivity[iRow];
        nUnder++;
        rowWeight[iRow] *= 1.1;
      } else if (rowActivity[iRow] > rowUpper[iRow] + 10.0 * primalTolerance) {
        gap += rowActivity[iRow] - rowUpper[iRow];
      } else {
        over += rowActivity[iRow] - rowLower[iRow];
        //rowWeight[iRow] *= 0.9;
      }
    }
    if (nG && !nL) {
      // can we fix
      // get list of columns which can go down without making
      // things much worse
      int nPossible = 0;
      int nEasyDown = 0;
      int nSlackDown = 0;
      for (int iColumn = 0; iColumn < numberColumns; iColumn++) {
        if (newSolution[iColumn] && columnUpper[iColumn] > columnLower[iColumn]) {
          bool canGoDown = true;
          bool under = false;
          int iSos = -1;
          for (CoinBigIndex j = columnStart[iColumn];
               j < columnStart[iColumn] + columnLength[iColumn]; j++) {
            int iRow = row[j];
            if (sos[iRow] < 0) {
              double over = rowActivity[iRow] - rowLower[iRow];
              if (over >= 0.0 && element[j] > over + 1.0e-12) {
                canGoDown = false;
                break;
              } else if (over < 0.0) {
                under = true;
              }
            } else {
              iSos = iRow;
            }
          }
          if (canGoDown) {
            if (!under) {
              if (iSos >= 0) {
                // find cheapest
                double cheapest = modifiedCost[iColumn];
                int iCheapest = -1;
                int jColumn = firstGub[iSos];
                assert(jColumn >= 0);
                while (jColumn >= 0) {
                  if (modifiedCost[jColumn] < cheapest) {
                    cheapest = modifiedCost[jColumn];
                    iCheapest = jColumn;
                  }
                  jColumn = nextGub[jColumn];
                }
                if (iCheapest >= 0) {
                  // Decrease column
                  newSolution[iColumn] = 0.0;
                  newSolutionValue -= modifiedCost[iColumn];
                  for (CoinBigIndex j = columnStart[iColumn];
                       j < columnStart[iColumn] + columnLength[iColumn]; j++) {
                    int iRow = row[j];
                    rowActivity[iRow] -= element[j];
                  }
                  // Increase chosen column
                  newSolution[iCheapest] = 1.0;
                  newSolutionValue += modifiedCost[iCheapest];
                  for (CoinBigIndex j = columnStart[iCheapest];
                       j < columnStart[iCheapest] + columnLength[iCheapest]; j++) {
                    int iRow = row[j];
                    rowActivity[iRow] += element[j];
                  }
                  nEasyDown++;
                  if (columnLength[iColumn] > 1) {
                    //printf("%d is easy down\n",iColumn);
                  } else {
                    nSlackDown++;
                  }
                }
              } else if (modifiedCost[iColumn] > 0.0) {
                // easy down
                // Decrease column
                newSolution[iColumn] = 0.0;
                newSolutionValue -= modifiedCost[iColumn];
                for (CoinBigIndex j = columnStart[iColumn];
                     j < columnStart[iColumn] + columnLength[iColumn]; j++) {
                  int iRow = row[j];
                  rowActivity[iRow] -= element[j];
                }
                nEasyDown++;
              }
            } else {
              which[nPossible++] = iColumn;
            }
          }
        }
      }
#ifdef REPORT
      printf("%d possible down, %d easy down of which %d are slacks\n",
        nPossible, nEasyDown, nSlackDown);
#endif
      double *needed = new double[numberRows];
      for (int i = 0; i < numberRows; i++) {
        double value = rowLower[i] - rowActivity[i];
        if (value < 1.0e-8)
          value = 0.0;
        needed[i] = value;
      }
      if (gap && /*nUnder==1 &&*/ nonSOS) {
        double *weight = new double[numberColumns];
        int *sort = new int[numberColumns];
        // look at ones not in set
        int nPossible = 0;
        for (int iColumn = 0; iColumn < numberColumns; iColumn++) {
          if (!newSolution[iColumn] && columnUpper[iColumn] > columnLower[iColumn]) {
            int iSos = -1;
            double value = 0.0;
            for (CoinBigIndex j = columnStart[iColumn];
                 j < columnStart[iColumn] + columnLength[iColumn]; j++) {
              int iRow = row[j];
              if (sos[iRow] < 0) {
                if (needed[iRow])
                  value += CoinMin(element[j] / needed[iRow], 1.0);
              } else {
                iSos = iRow;
              }
            }
            if (value && iSos < 0) {
              weight[nPossible] = -value;
              sort[nPossible++] = iColumn;
            }
          }
        }
        CoinSort_2(weight, weight + nPossible, sort);
        for (int i = 0; i < nPossible; i++) {
          int iColumn = sort[i];
          double helps = 0.0;
          for (CoinBigIndex j = columnStart[iColumn];
               j < columnStart[iColumn] + columnLength[iColumn]; j++) {
            int iRow = row[j];
            if (needed[iRow])
              helps += CoinMin(needed[iRow], element[j]);
          }
          if (helps) {
            newSolution[iColumn] = 1.0;
            newSolutionValue += modifiedCost[iColumn];
            for (CoinBigIndex j = columnStart[iColumn];
                 j < columnStart[iColumn] + columnLength[iColumn]; j++) {
              int iRow = row[j];
              rowActivity[iRow] += element[j];
              if (needed[iRow]) {
                needed[iRow] -= element[j];
                if (needed[iRow] < 1.0e-8)
                  needed[iRow] = 0.0;
              }
            }
            gap -= helps;
#ifdef REPORT
            {
              double gap2 = 0.0;
              for (iRow = 0; iRow < numberRows; iRow++) {
                if (rowActivity[iRow] < rowLower[iRow] - 10.0 * primalTolerance) {
                  gap2 += rowLower[iRow] - rowActivity[iRow];
                }
              }
              printf("estimated gap (nonsos) %g - computed %g\n",
                gap, gap2);
            }
#endif
            if (gap < 1.0e-12)
              break;
          }
        }
        delete[] weight;
        delete[] sort;
      }
      if (gap && nPossible /*&&nUnder==1*/ && true && model_->bestSolution()) {
        double *weight = new double[numberColumns];
        int *sort = new int[numberColumns];
        // look at ones in sets
        const double *goodSolution = model_->bestSolution();
        int nPossible = 0;
        double largestWeight = 0.0;
        for (int iColumn = 0; iColumn < numberColumns; iColumn++) {
          if (!newSolution[iColumn] && goodSolution[iColumn] && columnUpper[iColumn] > columnLower[iColumn]) {
            int iSos = -1;
            double value = 0.0;
            for (CoinBigIndex j = columnStart[iColumn];
                 j < columnStart[iColumn] + columnLength[iColumn]; j++) {
              int iRow = row[j];
              if (sos[iRow] < 0) {
                if (needed[iRow])
                  value += CoinMin(element[j] / needed[iRow], 1.0);
              } else {
                iSos = iRow;
              }
            }
            if (value && iSos >= 0) {
              // see if value bigger than current
              int jColumn = firstGub[iSos];
              assert(jColumn >= 0);
              while (jColumn >= 0) {
                if (newSolution[jColumn])
                  break;
                jColumn = nextGub[jColumn];
              }
              assert(jColumn >= 0);
              double value2 = 0.0;
              for (CoinBigIndex j = columnStart[jColumn];
                   j < columnStart[jColumn] + columnLength[jColumn]; j++) {
                int iRow = row[j];
                if (needed[iRow])
                  value2 += CoinMin(element[j] / needed[iRow], 1.0);
              }
              if (value > value2) {
                weight[nPossible] = -(value - value2);
                largestWeight = CoinMax(largestWeight, (value - value2));
                sort[nPossible++] = iColumn;
              }
            }
          }
        }
        if (nPossible) {
          double *temp = new double[numberRows];
          int *which2 = new int[numberRows];
          memset(temp, 0, numberRows * sizeof(double));
          // modify so ones just more than gap best
          if (largestWeight > gap && nUnder == 1) {
            double offset = 4 * largestWeight;
            for (int i = 0; i < nPossible; i++) {
              double value = -weight[i];
              if (value > gap - 1.0e-12)
                weight[i] = -(offset - (value - gap));
            }
          }
          CoinSort_2(weight, weight + nPossible, sort);
          for (int i = 0; i < nPossible; i++) {
            int iColumn = sort[i];
            int n = 0;
            // find jColumn
            int iSos = -1;
            for (CoinBigIndex j = columnStart[iColumn];
                 j < columnStart[iColumn] + columnLength[iColumn]; j++) {
              int iRow = row[j];
              temp[iRow] = element[j];
              which2[n++] = iRow;
              if (sos[iRow] >= 0) {
                iSos = iRow;
              }
            }
            int jColumn = firstGub[iSos];
            assert(jColumn >= 0);
            while (jColumn >= 0) {
              if (newSolution[jColumn])
                break;
              jColumn = nextGub[jColumn];
            }
            assert(jColumn >= 0);
            for (CoinBigIndex j = columnStart[jColumn];
                 j < columnStart[jColumn] + columnLength[jColumn]; j++) {
              int iRow = row[j];
              if (!temp[iRow])
                which2[n++] = iRow;
              temp[iRow] -= element[j];
            }
            double helps = 0.0;
            for (int i = 0; i < n; i++) {
              int iRow = which2[i];
              double newValue = rowActivity[iRow] + temp[iRow];
              if (temp[iRow] > 1.0e-8) {
                if (rowActivity[iRow] < rowLower[iRow] - 1.0e-8) {
                  helps += CoinMin(temp[iRow],
                    rowLower[iRow] - rowActivity[iRow]);
                }
              } else if (temp[iRow] < -1.0e-8) {
                if (newValue < rowLower[iRow] - 1.0e-12) {
                  helps -= CoinMin(-temp[iRow],
                    1.0 * (rowLower[iRow] - newValue));
                }
              }
            }
            if (helps > 0.0) {
              newSolution[iColumn] = 1.0;
              newSolution[jColumn] = 0.0;
              newSolutionValue += modifiedCost[iColumn] - modifiedCost[jColumn];
              for (int i = 0; i < n; i++) {
                int iRow = which2[i];
                double newValue = rowActivity[iRow] + temp[iRow];
                rowActivity[iRow] = newValue;
                temp[iRow] = 0.0;
              }
              gap -= helps;
#ifdef REPORT
              {
                double gap2 = 0.0;
                for (iRow = 0; iRow < numberRows; iRow++) {
                  if (rowActivity[iRow] < rowLower[iRow] - 10.0 * primalTolerance) {
                    gap2 += rowLower[iRow] - rowActivity[iRow];
                  }
                }
                printf("estimated gap %g - computed %g\n",
                  gap, gap2);
              }
#endif
              if (gap < 1.0e-8)
                break;
            } else {
              for (int i = 0; i < n; i++)
                temp[which2[i]] = 0.0;
            }
          }
          delete[] which2;
          delete[] temp;
        }
        delete[] weight;
        delete[] sort;
      }
      delete[] needed;
    }
#ifdef REPORT
    {
      double gap = 0.0;
      double over = 0.0;
      for (iRow = 0; iRow < numberRows; iRow++) {
        if (rowActivity[iRow] < rowLower[iRow] - 10.0 * primalTolerance) {
          double value = rowLower[iRow] - rowActivity[iRow];
#if REPORT > 1
          printf("below on %d is %g - activity %g lower %g\n",
            iRow, value, rowActivity[iRow], rowLower[iRow]);
#endif
          gap += value;
        } else if (rowActivity[iRow] > rowUpper[iRow] + 10.0 * primalTolerance) {
          double value = rowActivity[iRow] - rowUpper[iRow];
#if REPORT > 1
          printf("above on %d is %g - activity %g upper %g\n",
            iRow, value, rowActivity[iRow], rowUpper[iRow]);
#endif
          gap += value;
        } else {
          double value = rowActivity[iRow] - rowLower[iRow];
          if (value) {
#if REPORT > 1
            printf("over on %d is %g - activity %g lower %g\n",
              iRow, value, rowActivity[iRow], rowLower[iRow]);
#endif
            over += value;
          }
        }
      }
      printf("modified final gap %g - over %g - %d added - solvalue %g\n",
        gap, over, nAdded, newSolutionValue);
    }
#endif
    if (!gap) {
      break;
    } else {
      if (iPass == 0) {
        costBias = 10.0 * newSolutionValue / static_cast< double >(nAdded);
      } else {
        costBias *= 10.0;
      }
    }
  }
  delete[] newSolution0;
  delete[] rowWeight;
  delete[] sos;
  delete[] firstGub;
  delete[] nextGub;
  if (newSolutionValue < solutionValue) {
    // check feasible
    memset(rowActivity, 0, numberRows * sizeof(double));
    for (iColumn = 0; iColumn < numberColumns; iColumn++) {
      CoinBigIndex j;
      double value = newSolution[iColumn];
      if (value) {
        for (j = columnStart[iColumn];
             j < columnStart[iColumn] + columnLength[iColumn]; j++) {
          int iRow = row[j];
          rowActivity[iRow] += value * element[j];
        }
      }
    }
    // check was approximately feasible
    bool feasible = true;
    for (iRow = 0; iRow < numberRows; iRow++) {
      if (rowActivity[iRow] < rowLower[iRow]) {
        if (rowActivity[iRow] < rowLower[iRow] - 10.0 * primalTolerance)
          feasible = false;
      } else if (rowActivity[iRow] > rowUpper[iRow]) {
        if (rowActivity[iRow] > rowUpper[iRow] + 10.0 * primalTolerance)
          feasible = false;
      }
    }
    if (feasible) {
      // new solution
      memcpy(betterSolution, newSolution, numberColumns * sizeof(double));
      solutionValue = newSolutionValue;
      //printf("** Solution of %g found by rounding\n",newSolutionValue);
      returnCode = 1;
    } else {
      // Can easily happen
      //printf("Debug CbcHeuristicGreedySOS giving bad solution\n");
    }
  }
  delete[] sosRow;
  delete[] newSolution;
  delete[] rowActivity;
  delete[] modifiedCost;
  delete[] contribution;
  delete[] which;
  delete[] rhs;
  return returnCode;
}
// update model
void CbcHeuristicGreedySOS::setModel(CbcModel *model)
{
  delete[] originalRhs_;
  gutsOfConstructor(model);
  validate();
}
// Resets stuff if model changes
void CbcHeuristicGreedySOS::resetModel(CbcModel *model)
{
  delete[] originalRhs_;
  gutsOfConstructor(model);
}
// Validate model i.e. sets when_ to 0 if necessary (may be NULL)
void CbcHeuristicGreedySOS::validate()
{
  if (model_ && when() < 10) {
    if (model_->numberIntegers() != model_->numberObjects() && (model_->numberObjects() || (model_->specialOptions() & 1024) == 0)) {
      int numberOdd = 0;
      for (int i = 0; i < model_->numberObjects(); i++) {
        if (!model_->object(i)->canDoHeuristics())
          numberOdd++;
      }
      if (numberOdd)
        setWhen(0);
    }
    // Only works if coefficients positive and all rows L/G or SOS
    OsiSolverInterface *solver = model_->solver();
    const double *columnUpper = solver->getColUpper();
    const double *columnLower = solver->getColLower();
    const double *rowLower = solver->getRowLower();
    const double *rowUpper = solver->getRowUpper();

    int numberRows = solver->getNumRows();
    // Column copy
    const double *element = matrix_.getElements();
    const int *row = matrix_.getIndices();
    const CoinBigIndex *columnStart = matrix_.getVectorStarts();
    const int *columnLength = matrix_.getVectorLengths();
    bool good = true;
    assert(originalRhs_);
    for (int iRow = 0; iRow < numberRows; iRow++) {
      if (rowLower[iRow] == 1.0 && rowUpper[iRow] == 1.0) {
        // SOS
        originalRhs_[iRow] = -1.0;
      } else if (rowLower[iRow] > 0.0 && rowUpper[iRow] < 1.0e10) {
        good = false;
      } else if (rowUpper[iRow] < 0.0) {
        good = false;
      } else if (rowUpper[iRow] < 1.0e10) {
        originalRhs_[iRow] = rowUpper[iRow];
      } else {
        originalRhs_[iRow] = rowLower[iRow];
      }
    }
    int numberColumns = solver->getNumCols();
    for (int iColumn = 0; iColumn < numberColumns; iColumn++) {
      if (!columnLength[iColumn])
        continue;
      if (columnLower[iColumn] < 0.0 || columnUpper[iColumn] > 1.0)
        good = false;
      CoinBigIndex j;
      int nSOS = 0;
      if (!isHeuristicInteger(solver, iColumn))
        good = false;
      for (j = columnStart[iColumn];
           j < columnStart[iColumn] + columnLength[iColumn]; j++) {
        if (element[j] < 0.0)
          good = false;
        int iRow = row[j];
        if (originalRhs_[iRow] == -1.0) {
          if (element[j] != 1.0)
            good = false;
          nSOS++;
        }
      }
      if (nSOS > 1)
        good = false;
    }
    if (!good)
      setWhen(0); // switch off
  }
}

/* vi: softtabstop=2 shiftwidth=2 expandtab tabstop=2
*/
