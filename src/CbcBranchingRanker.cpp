// Copyright (C) 2025, Haroldo Gambini Santos and others.
// All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

#include "CbcBranchingRanker.hpp"

#include <algorithm>
#include <cmath>

// --- Construction ----------------------------------------------------------

CbcBranchingRanker::CbcBranchingRanker()
  : weightConflict_(0.0)
  , formula_(CONFLICT_MIN)
  , scalingPowerTrusted_(0.5)
  , scalingPowerUntrusted_(1.0)
  , nBoostsApplied_(0)
  , nZeroScore_(0)
  , headerPrinted_(false)
{
}

CbcBranchingRanker::CbcBranchingRanker(const CbcBranchingRanker &rhs)
  : weightConflict_(rhs.weightConflict_)
  , formula_(rhs.formula_)
  , scalingPowerTrusted_(rhs.scalingPowerTrusted_)
  , scalingPowerUntrusted_(rhs.scalingPowerUntrusted_)
  , nBoostsApplied_(0)
  , nZeroScore_(0)
  , headerPrinted_(false)
{
}

CbcBranchingRanker &CbcBranchingRanker::operator=(const CbcBranchingRanker &rhs)
{
  if (this != &rhs) {
    weightConflict_ = rhs.weightConflict_;
    formula_ = rhs.formula_;
    scalingPowerTrusted_ = rhs.scalingPowerTrusted_;
    scalingPowerUntrusted_ = rhs.scalingPowerUntrusted_;
  }
  return *this;
}

// --- Core computation ------------------------------------------------------

double CbcBranchingRanker::conflictScore(std::size_t d0, std::size_t d1) const
{
  switch (formula_) {
  case CONFLICT_SUM:
    return static_cast< double >(d0 + d1);
  case CONFLICT_PRODUCT:
    return std::sqrt(static_cast< double >(d0) * static_cast< double >(d1));
  case CONFLICT_MIN:
  default:
    return static_cast< double >(std::min(d0, d1));
  }
}

double CbcBranchingRanker::applyConflictBoost(double sortKey, std::size_t d0,
  std::size_t d1, bool trusted) const
{
  if (weightConflict_ == 0.0)
    return sortKey;

  double cs = conflictScore(d0, d1);
  if (cs <= 0.0) {
    ++nZeroScore_;
    return sortKey;
  }

  const double power = trusted ? scalingPowerTrusted_ : scalingPowerUntrusted_;
  const double scaledCs = (power == 1.0) ? cs : std::pow(cs, power);

  // sortKey is negative; multiplying by (1 + positive) makes it more negative
  // (higher priority in the sort).
  ++nBoostsApplied_;
  return sortKey * (1.0 + weightConflict_ * scaledCs);
}

// --- Diagnostics -----------------------------------------------------------

const char *CbcBranchingRanker::formulaName() const
{
  switch (formula_) {
  case CONFLICT_SUM:     return "sum";
  case CONFLICT_PRODUCT: return "product";
  case CONFLICT_MIN:
  default:               return "min";
  }
}

void CbcBranchingRanker::resetCounters() const
{
  nBoostsApplied_ = 0;
  nZeroScore_     = 0;
  headerPrinted_  = false;
}

/* vi: softtabstop=2 shiftwidth=2 expandtab tabstop=2
*/
