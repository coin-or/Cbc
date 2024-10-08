// Copyright (C) 2008, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

#if defined(_MSC_VER)
// Turn off compiler warning about long names
#pragma warning(disable : 4786)
#endif

#include "CbcHeuristicDiveVectorLength.hpp"
#include "CbcStrategy.hpp"

// Default Constructor
CbcHeuristicDiveVectorLength::CbcHeuristicDiveVectorLength()
  : CbcHeuristicDive()
{
}

// Constructor from model
CbcHeuristicDiveVectorLength::CbcHeuristicDiveVectorLength(CbcModel &model)
  : CbcHeuristicDive(model)
{
}

// Destructor
CbcHeuristicDiveVectorLength::~CbcHeuristicDiveVectorLength()
{
}

// Clone
CbcHeuristicDiveVectorLength *
CbcHeuristicDiveVectorLength::clone() const
{
  return new CbcHeuristicDiveVectorLength(*this);
}

// Create C++ lines to get to current state
void CbcHeuristicDiveVectorLength::generateCpp(FILE *fp)
{
  CbcHeuristicDiveVectorLength other;
  fprintf(fp, "0#include \"CbcHeuristicDiveVectorLength.hpp\"\n");
  fprintf(fp, "3  CbcHeuristicDiveVectorLength heuristicDiveVectorLength(*cbcModel);\n");
  CbcHeuristic::generateCpp(fp, "heuristicDiveVectorLength");
  fprintf(fp, "3  cbcModel->addHeuristic(&heuristicDiveVectorLength);\n");
}

// Copy constructor
CbcHeuristicDiveVectorLength::CbcHeuristicDiveVectorLength(const CbcHeuristicDiveVectorLength &rhs)
  : CbcHeuristicDive(rhs)
{
}

// Assignment operator
CbcHeuristicDiveVectorLength &
CbcHeuristicDiveVectorLength::operator=(const CbcHeuristicDiveVectorLength &rhs)
{
  if (this != &rhs) {
    CbcHeuristicDive::operator=(rhs);
  }
  return *this;
}

bool CbcHeuristicDiveVectorLength::selectVariableToBranch(OsiSolverInterface *solver,
  const double *newSolution,
  int &bestColumn,
  int &bestRound)
{
  const double *objective = solver->getObjCoefficients();
  double direction = solver->getObjSense(); // 1 for min, -1 for max

  const int *columnLength = matrix_.getVectorLengths();
  int numberIntegers = model_->numberIntegers();
  const int *integerVariable = model_->integerVariable();
  double integerTolerance = model_->getDblParam(CbcModel::CbcIntegerTolerance);

  bestColumn = -1;
  bestRound = -1; // -1 rounds down, +1 rounds up
  double bestScore = COIN_DBL_MAX;
  bool allTriviallyRoundableSoFar = true;
  int bestPriority = COIN_INT_MAX;
  for (int i = 0; i < numberIntegers; i++) {
    int iColumn = integerVariable[i];
    if (!isHeuristicInteger(solver, iColumn))
      continue;
    double value = newSolution[iColumn];
    double fraction = value - floor(value);
    int round = 0;
    if (fabs(floor(value + 0.5) - value) > integerTolerance) {
      if (allTriviallyRoundableSoFar || (downLocks_[i] > 0 && upLocks_[i] > 0)) {

        if (allTriviallyRoundableSoFar && downLocks_[i] > 0 && upLocks_[i] > 0) {
          allTriviallyRoundableSoFar = false;
          bestScore = COIN_DBL_MAX;
        }

        // the variable cannot be rounded
        double obj = direction * objective[iColumn];
        if (obj > smallObjective_) {
          round = 1; // round up
        } else if (obj < -smallObjective_) {
          round = -1; // round down
        } else {
          if (fraction < 0.4)
            round = -1;
          else
            round = 1;
        }
        double objDelta;
        if (round == 1)
          objDelta = (1.0 - fraction) * std::max(obj, smallObjective_);
        else
          objDelta = -fraction * std::min(obj, -smallObjective_);

        // we want the smaller score
        double score = objDelta / (static_cast< double >(columnLength[iColumn]) + 1.0);

        // if variable is not binary, penalize it
        if (!solver->isBinary(iColumn))
          score *= 1000.0;

        // if priorities then use
        if (priority_) {
          int thisRound = static_cast< int >(priority_[i].direction);
          if ((thisRound & 1) != 0)
            round = ((thisRound & 2) == 0) ? -1 : +1;
          if (priority_[i].priority > bestPriority) {
            score = COIN_DBL_MAX;
          } else if (priority_[i].priority < bestPriority) {
            bestPriority = static_cast< int >(priority_[i].priority);
            bestScore = COIN_DBL_MAX;
          }
        }
        if (score < bestScore) {
          bestColumn = iColumn;
          bestScore = score;
          bestRound = round;
        }
      }
    }
  }
  return allTriviallyRoundableSoFar;
}

/* vi: softtabstop=2 shiftwidth=2 expandtab tabstop=2
*/
