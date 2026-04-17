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
  , weightRange_(0.0)
  , scalingPowerRangeTrusted_(0.5)
  , scalingPowerRangeUntrusted_(1.0)
  , weightNonzeros_(0.0)
  , scalingPowerNzTrusted_(0.25)
  , scalingPowerNzUntrusted_(0.5)
  , nBoostsApplied_(0)
  , nZeroScore_(0)
  , nRangeBoostsApplied_(0)
  , nNzBoostsApplied_(0)
  , headerPrinted_(false)
{
}

CbcBranchingRanker::CbcBranchingRanker(const CbcBranchingRanker &rhs)
  : weightConflict_(rhs.weightConflict_)
  , formula_(rhs.formula_)
  , scalingPowerTrusted_(rhs.scalingPowerTrusted_)
  , scalingPowerUntrusted_(rhs.scalingPowerUntrusted_)
  , weightRange_(rhs.weightRange_)
  , scalingPowerRangeTrusted_(rhs.scalingPowerRangeTrusted_)
  , scalingPowerRangeUntrusted_(rhs.scalingPowerRangeUntrusted_)
  , weightNonzeros_(rhs.weightNonzeros_)
  , scalingPowerNzTrusted_(rhs.scalingPowerNzTrusted_)
  , scalingPowerNzUntrusted_(rhs.scalingPowerNzUntrusted_)
  , nBoostsApplied_(0)
  , nZeroScore_(0)
  , nRangeBoostsApplied_(0)
  , nNzBoostsApplied_(0)
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
    weightRange_ = rhs.weightRange_;
    scalingPowerRangeTrusted_ = rhs.scalingPowerRangeTrusted_;
    scalingPowerRangeUntrusted_ = rhs.scalingPowerRangeUntrusted_;
    weightNonzeros_ = rhs.weightNonzeros_;
    scalingPowerNzTrusted_ = rhs.scalingPowerNzTrusted_;
    scalingPowerNzUntrusted_ = rhs.scalingPowerNzUntrusted_;
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
  nBoostsApplied_      = 0;
  nZeroScore_          = 0;
  nRangeBoostsApplied_ = 0;
  nNzBoostsApplied_    = 0;
  headerPrinted_       = false;
}

double CbcBranchingRanker::applyNonzerosBoost(double sortKey, int nz,
  bool trusted) const
{
  if (weightNonzeros_ == 0.0 || nz <= 0)
    return sortKey;

  // Raw score = nz raised to a slow-growing power (default 0.25 = 4th root).
  // No normalization needed: 4th-root grows so slowly it stays in a small range
  // (nz=1→1.0, nz=16→2.0, nz=100→3.16) and the weight keeps the boost mild.
  const double power = trusted ? scalingPowerNzTrusted_ : scalingPowerNzUntrusted_;
  const double score = std::pow(static_cast< double >(nz), power);

  ++nNzBoostsApplied_;
  return sortKey * (1.0 + weightNonzeros_ * score);
}

double CbcBranchingRanker::applyRangeBoost(double sortKey, double lb, double ub,
  bool trusted) const
{
  if (weightRange_ == 0.0)
    return sortKey;

  const double range = ub - lb;
  if (range <= 0.0)
    return sortKey;

  // Score = 1/range: binary [0,1] → 1.0, [0,9] → 0.111, [0,99] → 0.010.
  // Naturally in (0, 1] for integers (minimum range = 1).
  const double score = 1.0 / range;
  const double power = trusted ? scalingPowerRangeTrusted_ : scalingPowerRangeUntrusted_;
  const double scaledScore = (power == 1.0) ? score : std::pow(score, power);

  ++nRangeBoostsApplied_;
  return sortKey * (1.0 + weightRange_ * scaledScore);
}

/* vi: softtabstop=2 shiftwidth=2 expandtab tabstop=2
*/
