/* $Id$ */
// Copyright (C) 2008, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

#if defined(_MSC_VER)
// Turn off compiler warning about long names
#pragma warning(disable : 4786)
#endif

#include "CbcHeuristicDiveLineSearch.hpp"
#include "CbcStrategy.hpp"

// Default Constructor
CbcHeuristicDiveLineSearch::CbcHeuristicDiveLineSearch()
  : CbcHeuristicDive()
{
}

// Constructor from model
CbcHeuristicDiveLineSearch::CbcHeuristicDiveLineSearch(CbcModel &model)
  : CbcHeuristicDive(model)
{
}

// Destructor
CbcHeuristicDiveLineSearch::~CbcHeuristicDiveLineSearch()
{
}

// Clone
CbcHeuristicDiveLineSearch *
CbcHeuristicDiveLineSearch::clone() const
{
  return new CbcHeuristicDiveLineSearch(*this);
}

// Create C++ lines to get to current state
void CbcHeuristicDiveLineSearch::generateCpp(FILE *fp)
{
  CbcHeuristicDiveLineSearch other;
  fprintf(fp, "0#include \"CbcHeuristicDiveLineSearch.hpp\"\n");
  fprintf(fp, "3  CbcHeuristicDiveLineSearch heuristicDiveLineSearch(*cbcModel);\n");
  CbcHeuristic::generateCpp(fp, "heuristicDiveLineSearch");
  fprintf(fp, "3  cbcModel->addHeuristic(&heuristicDiveLineSearch);\n");
}

// Copy constructor
CbcHeuristicDiveLineSearch::CbcHeuristicDiveLineSearch(const CbcHeuristicDiveLineSearch &rhs)
  : CbcHeuristicDive(rhs)
{
}

// Assignment operator
CbcHeuristicDiveLineSearch &
CbcHeuristicDiveLineSearch::operator=(const CbcHeuristicDiveLineSearch &rhs)
{
  if (this != &rhs) {
    CbcHeuristicDive::operator=(rhs);
  }
  return *this;
}

bool CbcHeuristicDiveLineSearch::selectVariableToBranch(OsiSolverInterface *solver,
  const double *newSolution,
  int &bestColumn,
  int &bestRound)
{
  int numberIntegers = model_->numberIntegers();
  const int *integerVariable = model_->integerVariable();
  double integerTolerance = model_->getDblParam(CbcModel::CbcIntegerTolerance);

  // get the LP relaxation solution at the root node
  double *rootNodeLPSol = model_->continuousSolution();

  bestColumn = -1;
  bestRound = -1; // -1 rounds down, +1 rounds up
  double bestRelDistance = COIN_DBL_MAX;
  bool allTriviallyRoundableSoFar = true;
  int bestPriority = COIN_INT_MAX;
  for (int i = 0; i < numberIntegers; i++) {
    int iColumn = integerVariable[i];
    if (!isHeuristicInteger(solver, iColumn))
      continue;
    double rootValue = rootNodeLPSol[iColumn];
    double value = newSolution[iColumn];
    double fraction = value - floor(value);
    int round = 0;
    if (fabs(floor(value + 0.5) - value) > integerTolerance) {
      if (allTriviallyRoundableSoFar || (downLocks_[i] > 0 && upLocks_[i] > 0)) {

        if (allTriviallyRoundableSoFar && downLocks_[i] > 0 && upLocks_[i] > 0) {
          allTriviallyRoundableSoFar = false;
          bestRelDistance = COIN_DBL_MAX;
        }

        double relDistance;
        if (value < rootValue) {
          round = -1;
          relDistance = fraction / (rootValue - value);
        } else if (value > rootValue) {
          round = 1;
          relDistance = (1.0 - fraction) / (value - rootValue);
        } else {
          round = -1;
          relDistance = COIN_DBL_MAX;
        }

        // if variable is not binary, penalize it
        if (!solver->isBinary(iColumn))
          relDistance *= 1000.0;

        // if priorities then use
        if (priority_) {
          int thisRound = static_cast< int >(priority_[i].direction);
          if ((thisRound & 1) != 0)
            round = ((thisRound & 2) == 0) ? -1 : +1;
          if (priority_[i].priority > bestPriority) {
            relDistance = COIN_DBL_MAX;
          } else if (priority_[i].priority < bestPriority) {
            bestPriority = static_cast< int >(priority_[i].priority);
            bestRelDistance = COIN_DBL_MAX;
          }
        }
        if (relDistance < bestRelDistance) {
          bestColumn = iColumn;
          bestRelDistance = relDistance;
          bestRound = round;
        }
      }
    }
  }
  return allTriviallyRoundableSoFar;
}

/* vi: softtabstop=2 shiftwidth=2 expandtab tabstop=2
*/
