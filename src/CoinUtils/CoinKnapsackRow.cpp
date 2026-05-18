/* -*- mode: C++; tab-width: 2; indent-tabs-mode: nil; -*-
 *
 * This file is part of the COIN-OR CoinUtils package
 *
 * @file   CoinKnapsackRow.cpp
 * @brief  Helper class to analyze a single constraint row and store it in knapsack form, i.e.
 *  only binary variables and positive coefficients.
 *
 * Copyright (C) 2025
 * All rights reserved.
 *
 * This code is licensed under the terms of the Eclipse Public License (EPL).
 */

#include "CoinKnapsackRow.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <limits>

#include "CoinPragma.hpp"

namespace {
inline void update_two_largest(double val, double (&v)[2])
{
  if (val > v[0]) {
    v[1] = v[0];
    v[0] = val;
  } else if (val > v[1]) {
    v[1] = val;
  }
}

inline void update_two_smallest(double val, double (&v)[2])
{
  if (val < v[0]) {
    v[1] = v[0];
    v[0] = val;
  } else if (val < v[1]) {
    v[1] = val;
  }
}
} // namespace

CoinKnapsackRow::CoinKnapsackRow(
  size_t numCols,
  const char *colType,
  const double *colLB,
  const double *colUB,
  double primalTolerance,
  double infinity)
  : colType_(colType)
  , colLB_(colLB)
  , colUB_(colUB)
  , columns_(new CoinTerm[numCols])
  , nzs_(0)
  , numCols_(numCols)
  , fixedVariables_(new int[numCols])
  , nFixed_(0)
  , rhs_(std::numeric_limits< double >::quiet_NaN())
  , primalTolerance_(primalTolerance)
  , infinity_(infinity)
{
}

CoinKnapsackRow::~CoinKnapsackRow()
{
  delete[] columns_;
  delete[] fixedVariables_;
}

//! Reset cached values so the object can process another row.
void CoinKnapsackRow::resetComputedValues()
{
  rhs_ = std::numeric_limits< double >::quiet_NaN();
  nzs_ = nFixed_ = 0;
  twoLargest_[0] = twoLargest_[1] = -infinity_;
  twoSmallest_[0] = twoSmallest_[1] = +infinity_;
}

//! Sort stored columns from smallest to largest coefficient.
void CoinKnapsackRow::sortColumns()
{
  std::sort(columns_, columns_ + nzs_, CoinTerm::AscendingValueThenIndex);
}

bool CoinKnapsackRow::isExplicitClique() const
{
  if (nzs_ < 2) {
    return false;
  }
  if (!std::isfinite(rhs_) || rhs_ < 0.0 || rhs_ >= infinity_) {
    return false;
  }
  if (!std::isfinite(twoSmallest_[0]) || !std::isfinite(twoSmallest_[1])) {
    return false;
  }
  return (twoSmallest_[0] + twoSmallest_[1]) > (rhs_ + primalTolerance_);
}

/*! Analyse a single constraint row, to reduce it to a
  knapsack structure. Uses the tolerance/infinity
  provided at construction time. */
void CoinKnapsackRow::processRow(
  const int *idx,
  const double *coef,
  size_t nz,
  char sense,
  double multiplier,
  double rhs)
{
  resetComputedValues();

  const double mult = multiplier;
  rhs_ = mult * rhs;

  // discount fixed variables from RHS
  for (size_t j = 0; j < nz; ++j) {
    const size_t idxCol = static_cast< size_t >(idx[j]);
    if (colLB_[idxCol] == colUB_[idxCol]) {
      const double coefCol = coef[j] * mult;
      rhs_ -= coefCol * colLB_[idxCol];
    }
  }

  for (size_t j = 0; j < nz; ++j) {
    const int idxCol = idx[j];

    // fixed variables already considered
    if (colLB_[idxCol] == colUB_[idxCol]) {
      continue;
    }

    const char cType = colType_[idxCol];

    if (cType == CoinColumnType::SemiContinuous || cType == CoinColumnType::SemiInteger) {
      // not handled yet
      rhs_ = std::numeric_limits< double >::quiet_NaN();
      return;
    }

    const double coefCol = coef[j] * mult;

    if (cType == CoinColumnType::Continuous || cType == CoinColumnType::GeneralInteger) {
      if (coefCol < 0.0) {
        if (colUB_[idxCol] >= infinity_) {
          rhs_ = std::numeric_limits< double >::quiet_NaN();
          return;
        }
        rhs_ -= coefCol * colUB_[idxCol];
      } else if (coefCol > 0.0) {
        if (colLB_[idxCol] <= -infinity_) {
          rhs_ = std::numeric_limits< double >::quiet_NaN();
          return;
        }
        rhs_ -= coefCol * colLB_[idxCol];
      }
      continue;
    }

#ifdef DEBUGCG
    assert(cType == CoinColumnType::Binary);
    assert(nzs_ < numCols_);
#endif

    CoinTerm &newTerm = columns_[nzs_++];
    if (coefCol >= 0.0) {
      newTerm.index = static_cast< int >(idxCol);
      newTerm.value = coefCol;
    } else {
      newTerm.index = static_cast< int >(idxCol + numCols_);
      newTerm.value = -coefCol;
      rhs_ += newTerm.value;
    }

    update_two_largest(newTerm.value, twoLargest_);
    update_two_smallest(newTerm.value, twoSmallest_);
  }

  if (nzs_ == 0) {
    return;
  }

  for (size_t j = 0; j < nzs_; ++j) {
    if (columns_[j].value > rhs_ + primalTolerance_) {
#ifdef DEBUGCG
      assert(nFixed_ < numCols_);
#endif
      fixedVariables_[nFixed_++] = columns_[j].index;
    }
  }
}

void CoinKnapsackRow::copyColumnIndices(size_t *dest) const
{
  for (size_t i = 0; i < nzs_; ++i) {
    dest[i] = static_cast< size_t >(columns_[i].index);
  }
}

int CoinKnapsackRow::rowIterations(
  const char sense,
  const double rhs,
  const double range,
  double multiplier[],
  double adjustedRHS[])
{
  switch (sense) {
  case 'L':
    multiplier[0] = 1.0;
    adjustedRHS[0] = rhs;
    return 1;
  case 'G':
    multiplier[0] = -1.0;
    adjustedRHS[0] = rhs;
    return 1;
  case 'E':
    multiplier[0] = 1.0;
    adjustedRHS[0] = rhs;
    multiplier[1] = -1.0;
    adjustedRHS[1] = rhs;
    return 2;
  case 'R':
    multiplier[0] = 1.0;
    adjustedRHS[0] = rhs;
    multiplier[1] = -1.0;
    adjustedRHS[1] = (rhs - range);
    return 2;
  default:
    // unknown sense
    return 0;
  }
}
