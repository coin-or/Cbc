// Copyright (C) 2026, Haroldo Gambini Santos and others.
// All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

#if defined(_MSC_VER)
#pragma warning(disable : 4786)
#endif

#include "CbcHeuristicDiveConfigurable.hpp"
#include "CbcModel.hpp"
#include "CbcStrategy.hpp"
#include "CoinStaticConflictGraph.hpp"
#include "OsiSolverInterface.hpp"

// Default Constructor
CbcHeuristicDiveConfigurable::CbcHeuristicDiveConfigurable()
  : CbcHeuristicDive()
  , minFractionality_(0.01)
  , weightFractionality_(1.0)
  , powerFractionality_(1.0)
  , weightLocks_(0.1)
  , powerLocks_(0.5)
  , weightConflict_(0.0)
  , powerConflict_(0.5)
  , weightObjCoeff_(0.0)
  , powerObjCoeff_(0.2)
  , weightNonzeros_(0.0)
  , powerNonzeros_(0.25)
  , useConflictForDirection_(true)
  , preferObjectiveDirection_(false)
  , randomFactor_(0.0)
  , targetFractionality_(0.5)
  , fixCount_(0)
{
  whereFrom_ |= 16 * (1 + 256);
}

// Constructor from model
CbcHeuristicDiveConfigurable::CbcHeuristicDiveConfigurable(CbcModel &model)
  : CbcHeuristicDive(model)
  , minFractionality_(0.01)
  , weightFractionality_(1.0)
  , powerFractionality_(1.0)
  , weightLocks_(0.1)
  , powerLocks_(0.5)
  , weightConflict_(0.0)
  , powerConflict_(0.5)
  , weightObjCoeff_(0.0)
  , powerObjCoeff_(0.2)
  , weightNonzeros_(0.0)
  , powerNonzeros_(0.25)
  , useConflictForDirection_(true)
  , preferObjectiveDirection_(false)
  , randomFactor_(0.0)
  , targetFractionality_(0.5)
  , fixCount_(0)
{
  whereFrom_ |= 16 * (1 + 256);
}

// Copy constructor
CbcHeuristicDiveConfigurable::CbcHeuristicDiveConfigurable(
  const CbcHeuristicDiveConfigurable &rhs)
  : CbcHeuristicDive(rhs)
  , minFractionality_(rhs.minFractionality_)
  , weightFractionality_(rhs.weightFractionality_)
  , powerFractionality_(rhs.powerFractionality_)
  , weightLocks_(rhs.weightLocks_)
  , powerLocks_(rhs.powerLocks_)
  , weightConflict_(rhs.weightConflict_)
  , powerConflict_(rhs.powerConflict_)
  , weightObjCoeff_(rhs.weightObjCoeff_)
  , powerObjCoeff_(rhs.powerObjCoeff_)
  , weightNonzeros_(rhs.weightNonzeros_)
  , powerNonzeros_(rhs.powerNonzeros_)
  , useConflictForDirection_(rhs.useConflictForDirection_)
  , preferObjectiveDirection_(rhs.preferObjectiveDirection_)
  , randomFactor_(rhs.randomFactor_)
  , targetFractionality_(rhs.targetFractionality_)
  , fixCount_(rhs.fixCount_)
{
}

// Destructor
CbcHeuristicDiveConfigurable::~CbcHeuristicDiveConfigurable()
{
}

// Clone
CbcHeuristicDiveConfigurable *
CbcHeuristicDiveConfigurable::clone() const
{
  return new CbcHeuristicDiveConfigurable(*this);
}

// Assignment operator
CbcHeuristicDiveConfigurable &
CbcHeuristicDiveConfigurable::operator=(const CbcHeuristicDiveConfigurable &rhs)
{
  if (this != &rhs) {
    CbcHeuristicDive::operator=(rhs);
    minFractionality_ = rhs.minFractionality_;
    weightFractionality_ = rhs.weightFractionality_;
    powerFractionality_ = rhs.powerFractionality_;
    weightLocks_ = rhs.weightLocks_;
    powerLocks_ = rhs.powerLocks_;
    weightConflict_ = rhs.weightConflict_;
    powerConflict_ = rhs.powerConflict_;
    weightObjCoeff_ = rhs.weightObjCoeff_;
    powerObjCoeff_ = rhs.powerObjCoeff_;
    weightNonzeros_ = rhs.weightNonzeros_;
    powerNonzeros_ = rhs.powerNonzeros_;
    useConflictForDirection_ = rhs.useConflictForDirection_;
    preferObjectiveDirection_ = rhs.preferObjectiveDirection_;
    randomFactor_ = rhs.randomFactor_;
    targetFractionality_ = rhs.targetFractionality_;
    fixCount_ = rhs.fixCount_;
  }
  return *this;
}

void CbcHeuristicDiveConfigurable::initializeData()
{
  // If fixCount_ is set, override percentageToFix_ to produce that count
  if (fixCount_ > 0 && model_->numberIntegers() > 0) {
    percentageToFix_ = static_cast<double>(fixCount_) / model_->numberIntegers();
  }
}

// Create C++ lines to get to current state
void CbcHeuristicDiveConfigurable::generateCpp(FILE *fp)
{
  CbcHeuristicDiveConfigurable other;
  fprintf(fp, "0#include \"CbcHeuristicDiveConfigurable.hpp\"\n");
  fprintf(fp, "3  CbcHeuristicDiveConfigurable heuristicDiveConfigurable(*cbcModel);\n");
  CbcHeuristic::generateCpp(fp, "heuristicDiveConfigurable");
  fprintf(fp, "3  cbcModel->addHeuristic(&heuristicDiveConfigurable);\n");
}

bool CbcHeuristicDiveConfigurable::selectVariableToBranch(
  OsiSolverInterface *solver,
  const double *newSolution,
  int &bestColumn,
  int &bestRound)
{
  int numberIntegers = model_->numberIntegers();
  const int *integerVariable = model_->integerVariable();
  double integerTolerance = model_->getDblParam(CbcModel::CbcIntegerTolerance);
  const double *lower = solver->getColLower();
  const double *upper = solver->getColUpper();
  const double *objective = solver->getObjCoefficients();
  double direction = solver->getObjSenseInCbc();

  // Get conflict graph (may be NULL)
  const CoinStaticConflictGraph *cgraph = solver->getCGraph();
  int numCols = solver->getNumCols();

  // Column lengths for nonzeros criterion
  const CoinPackedMatrix *matrix = solver->getMatrixByCol();
  const int *columnLength = matrix ? matrix->getVectorLengths() : nullptr;

  bestColumn = -1;
  bestRound = -1;
  double bestScore = -1.0;
  bool allTriviallyRoundable = true;
  // Fallback: track best candidate with relaxed threshold
  int fallbackColumn = -1;
  int fallbackRound = -1;
  double fallbackScore = -1.0;

  for (int i = 0; i < numberIntegers; i++) {
    int iColumn = integerVariable[i];
    if (!isHeuristicInteger(solver, iColumn))
      continue;
    double value = newSolution[iColumn];
    if (value < lower[iColumn] || value > upper[iColumn])
      continue;

    double fraction = value - floor(value);
    double distFromInt = std::min(fraction, 1.0 - fraction);

    // Must be fractional at all
    if (distFromInt < integerTolerance)
      continue;

    int nDownLocks = downLocks_[i];
    int nUpLocks = upLocks_[i];

    // Trivially roundable check
    if (nDownLocks == 0 || nUpLocks == 0) {
      if (allTriviallyRoundable)
        continue; // skip trivially roundable when scoring non-trivial
      else
        continue; // already found non-trivial candidates
    }

    // First non-trivially-roundable variable resets the search
    if (allTriviallyRoundable) {
      allTriviallyRoundable = false;
      bestScore = -1.0;
      fallbackScore = -1.0;
    }

    // Skip near-integer variables (but track as fallback)
    bool belowThreshold = (distFromInt < minFractionality_);

    // --- Compute score ---
    double score = 0.0;

    if (weightFractionality_ > 0.0) {
      // Score based on proximity to target fractionality
      // targetFractionality_=0.5: prefer most fractional (closest to 0.5)
      // targetFractionality_=0.6: prefer vars near 0.6 (biased toward rounding up)
      double frac = value - floor(value);
      double proximity = 1.0 - fabs(frac - targetFractionality_);
      double fracScore = (powerFractionality_ == 1.0)
        ? proximity
        : std::pow(std::max(proximity, 0.01), powerFractionality_);
      score += weightFractionality_ * fracScore;
    }

    if (weightLocks_ > 0.0) {
      int minLocks = std::min(nDownLocks, nUpLocks);
      double lockScore = 1.0 / (1.0 + minLocks);
      if (powerLocks_ != 1.0)
        lockScore = std::pow(lockScore, powerLocks_);
      score += weightLocks_ * lockScore;
    }

    if (weightConflict_ > 0.0 && cgraph && solver->isBinary(iColumn)) {
      size_t d1 = cgraph->degree(static_cast<size_t>(iColumn));
      size_t d0 = cgraph->degree(static_cast<size_t>(iColumn + numCols));
      double cs = static_cast<double>(d0 + d1);
      if (cs > 0.0) {
        if (powerConflict_ != 1.0)
          cs = std::pow(cs, powerConflict_);
        score += weightConflict_ * cs;
      }
    }

    if (weightObjCoeff_ > 0.0) {
      double absObj = fabs(objective[iColumn]);
      if (absObj > 0.0) {
        double objScore = (powerObjCoeff_ == 1.0)
          ? absObj
          : std::pow(absObj, powerObjCoeff_);
        score += weightObjCoeff_ * objScore;
      }
    }

    if (weightNonzeros_ > 0.0 && columnLength) {
      int nz = columnLength[iColumn];
      if (nz > 0) {
        double nzScore = (powerNonzeros_ == 1.0)
          ? static_cast<double>(nz)
          : std::pow(static_cast<double>(nz), powerNonzeros_);
        score += weightNonzeros_ * nzScore;
      }
    }

    // Determine rounding direction
    int round = 0;
    if (nDownLocks < nUpLocks)
      round = -1;
    else if (nDownLocks > nUpLocks)
      round = 1;
    else {
      if (preferObjectiveDirection_)
        round = (direction * objective[iColumn] >= 0.0) ? -1 : 1;
      else
        round = (fraction < 0.5) ? -1 : 1;
    }

    // Conflict-graph direction override (only for ties)
    if (useConflictForDirection_ && cgraph && solver->isBinary(iColumn)
      && nDownLocks == nUpLocks) {
      size_t d1 = cgraph->degree(static_cast<size_t>(iColumn));
      size_t d0 = cgraph->degree(static_cast<size_t>(iColumn + numCols));
      if (d1 > d0) round = -1;
      else if (d0 > d1) round = 1;
    }

    if (priority_) {
      int thisRound = static_cast<int>(priority_[i].direction);
      if ((thisRound & 1) != 0)
        round = ((thisRound & 2) == 0) ? -1 : +1;
    }

    if (belowThreshold) {
      // Track as fallback only
      if (score > fallbackScore) {
        fallbackScore = score;
        fallbackColumn = iColumn;
        fallbackRound = round;
      }
    } else {
      // Add random noise for diversification
      if (randomFactor_ > 0.0)
        score *= (1.0 + randomFactor_ * (randomNumberGenerator_.randomDouble() - 0.5));
      // Primary candidate
      if (score > bestScore) {
        bestScore = score;
        bestColumn = iColumn;
        bestRound = round;
      }
    }
  }

  // If no variable above threshold, use fallback (relaxed threshold)
  if (bestColumn < 0 && fallbackColumn >= 0) {
    bestColumn = fallbackColumn;
    bestRound = fallbackRound;
  }

  return allTriviallyRoundable;
}

/* vi: softtabstop=2 shiftwidth=2 expandtab tabstop=2
*/
