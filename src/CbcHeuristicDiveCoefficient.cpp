// Copyright (C) 2008, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

#if defined(_MSC_VER)
// Turn off compiler warning about long names
#pragma warning(disable : 4786)
#endif

//#define PRINT_DEBUG

#include "CbcHeuristicDiveCoefficient.hpp"
#include "CbcStrategy.hpp"

// Default Constructor
CbcHeuristicDiveCoefficient::CbcHeuristicDiveCoefficient()
  : CbcHeuristicDive()
{
  whereFrom_ |= 16 * (1 + 256);
}

// Constructor from model
CbcHeuristicDiveCoefficient::CbcHeuristicDiveCoefficient(CbcModel &model)
  : CbcHeuristicDive(model)
{
  whereFrom_ |= 16 * (1 + 256);
}

// Destructor
CbcHeuristicDiveCoefficient::~CbcHeuristicDiveCoefficient()
{
}

// Clone
CbcHeuristicDiveCoefficient *
CbcHeuristicDiveCoefficient::clone() const
{
  return new CbcHeuristicDiveCoefficient(*this);
}

// Create C++ lines to get to current state
void CbcHeuristicDiveCoefficient::generateCpp(FILE *fp)
{
  CbcHeuristicDiveCoefficient other;
  fprintf(fp, "0#include \"CbcHeuristicDiveCoefficient.hpp\"\n");
  fprintf(fp, "3  CbcHeuristicDiveCoefficient heuristicDiveCoefficient(*cbcModel);\n");
  CbcHeuristic::generateCpp(fp, "heuristicDiveCoefficient");
  fprintf(fp, "3  cbcModel->addHeuristic(&heuristicDiveCoefficient);\n");
}

// Copy constructor
CbcHeuristicDiveCoefficient::CbcHeuristicDiveCoefficient(const CbcHeuristicDiveCoefficient &rhs)
  : CbcHeuristicDive(rhs)
{
}

// Assignment operator
CbcHeuristicDiveCoefficient &
CbcHeuristicDiveCoefficient::operator=(const CbcHeuristicDiveCoefficient &rhs)
{
  if (this != &rhs) {
    CbcHeuristicDive::operator=(rhs);
  }
  return *this;
}

bool CbcHeuristicDiveCoefficient::selectVariableToBranch(OsiSolverInterface *solver,
  const double *newSolution,
  int &bestColumn,
  int &bestRound)
{
  int numberIntegers = model_->numberIntegers();
  const int *integerVariable = model_->integerVariable();
  double integerTolerance = model_->getDblParam(CbcModel::CbcIntegerTolerance);
  const double * lower = solver->getColLower();
  const double * upper = solver->getColUpper();

  bestColumn = -1;
  bestRound = -1; // -1 rounds down, +1 rounds up
  double bestFraction = COIN_DBL_MAX;
  int bestLocks = COIN_INT_MAX;
  bool allTriviallyRoundableSoFar = true;
  int bestPriority = COIN_INT_MAX;
  for (int i = 0; i < numberIntegers; i++) {
    int iColumn = integerVariable[i];
    if (!isHeuristicInteger(solver, iColumn))
      continue;
    double value = newSolution[iColumn];
    // deal with tolerance problems
    if (value<lower[iColumn] || value>upper[iColumn])
      continue;
    double fraction = value - floor(value);
    int round = 0;
    if (fabs(floor(value + 0.5) - value) > integerTolerance) {
      int nDownLocks = downLocks_[i];
      int nUpLocks = upLocks_[i];
      if (allTriviallyRoundableSoFar || (nDownLocks > 0 && nUpLocks > 0)) {

        if (allTriviallyRoundableSoFar && nDownLocks > 0 && nUpLocks > 0) {
          allTriviallyRoundableSoFar = false;
          bestFraction = COIN_DBL_MAX;
          bestLocks = COIN_INT_MAX;
        }

        // the variable cannot be rounded
        int nLocks = nDownLocks;
        if (nDownLocks < nUpLocks)
          round = -1;
        else if (nDownLocks > nUpLocks) {
          round = 1;
          fraction = 1.0 - fraction;
          nLocks = nUpLocks;
        } else if (fraction < 0.5)
          round = -1;
        else {
          round = 1;
          fraction = 1.0 - fraction;
          nLocks = nUpLocks;
        }

        // if variable is not binary, penalize it
        if (!solver->isBinary(iColumn))
          fraction *= 1000.0;

        // if priorities then use
        if (priority_) {
          int thisRound = static_cast< int >(priority_[i].direction);
          if ((thisRound & 1) != 0)
            round = ((thisRound & 2) == 0) ? -1 : +1;
          if (priority_[i].priority > bestPriority) {
            nLocks = COIN_INT_MAX;
          } else if (priority_[i].priority < bestPriority) {
            bestPriority = static_cast< int >(priority_[i].priority);
            bestLocks = COIN_INT_MAX;
          }
        }
        if (nLocks < bestLocks || (nLocks == bestLocks && fraction < bestFraction)) {
          bestColumn = iColumn;
          bestLocks = nLocks;
          bestFraction = fraction;
          bestRound = round;
        }
      }
    }
  }
  return allTriviallyRoundableSoFar;
}

/* vi: softtabstop=2 shiftwidth=2 expandtab tabstop=2
*/
