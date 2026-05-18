/* -*- mode: C++; tab-width: 2; indent-tabs-mode: nil; -*-
 *
 * This file is part of the COIN-OR CoinUtils package
 *
 * @file   CoinBoundPropagation.cpp
 * @brief  Bound propagation for binary variables in a MILP.
 *
 * Copyright (C) 2025
 * All rights reserved.
 *
 * This code is licensed under the terms of the Eclipse Public License (EPL).
 */

#include "CoinBoundPropagation.hpp"

#include <algorithm>
#include <cmath>
#include <vector>

#include "CoinPragma.hpp"
#include "CoinKnapsackRow.hpp"
#include "CoinPackedMatrix.hpp"
#include "CoinTypes.h"

#ifdef COIN_BT_STATS
#include <chrono>
#endif

namespace {

/**
 * Quick single-pass pre-check: returns true if the row iteration described
 * by (multiplier, adjustedRhs) might produce a variable fixing or reveal
 * infeasibility.  Returns false when it is certain no such outcome is
 * possible, allowing processRow to be skipped entirely.
 *
 * This check mirrors the transformation in CoinKnapsackRow::processRow and
 * is always sound: it never returns false for a row that would actually yield
 * a fixing or detect infeasibility.
 */
static bool rowNeedsProcessing(
  const int *idx, const double *coef, size_t nz,
  double multiplier, double adjustedRhs,
  const double *colLB, const double *colUB, const char *colType,
  double primalTolerance, double infinity)
{
  double rhs = multiplier * adjustedRhs;
  double maxCoef = 0.0;
  // Track whether any binary with c > 0 has been encountered.  Once rhs >= maxCoef
  // >= 0 is established *before* any c>0 binary, remaining c<0 entries can only
  // widen the gap rhs-maxCoef (proven: for any c<0 update, rhs grows by |c| and
  // maxCoef grows by at most |c|, so rhs-maxCoef is non-decreasing when rhs>=0).
  // A c>0 binary increases maxCoef without touching rhs, so we must keep scanning.
  bool seenPositiveC = false;

  for (size_t j = 0; j < nz; ++j) {
    const int col = idx[j];
    const double c = coef[j] * multiplier;

    if (colLB[col] == colUB[col]) {
      rhs -= c * colLB[col];
      continue;
    }

    const char ct = colType[col];
    if (ct == CoinColumnType::SemiContinuous || ct == CoinColumnType::SemiInteger)
      return true; // processRow will abort on semi types; let it do so
    if (ct == CoinColumnType::Continuous || ct == CoinColumnType::GeneralInteger) {
      if (c < 0.0) {
        if (colUB[col] >= infinity)
          return false; // processRow would return unbounded (NaN), no fixings
        rhs -= c * colUB[col];
      } else if (c > 0.0) {
        if (colLB[col] <= -infinity)
          return false;
        rhs -= c * colLB[col];
      }
      continue;
    }

    // Binary variable: mirror the complementation logic in processRow.
    if (c >= 0.0) {
      seenPositiveC = true;
      maxCoef = std::max(maxCoef, c);
    } else {
      rhs += (-c); // complementation adds |c| to rhs (same as processRow)
      maxCoef = std::max(maxCoef, -c);
      // Early exit: once rhs >= maxCoef >= 0 and no c>0 binary has been seen,
      // all remaining entries are either non-binary (already handled) or c<0
      // binaries, which can only preserve rhs >= maxCoef.  No fixing possible.
      if (!seenPositiveC && rhs >= 0.0 && maxCoef <= rhs + primalTolerance)
        return false;
    }
  }

  // Potential infeasibility: processRow must be called so it can abort the
  // whole pass via the infeasible_ flag.
  if (rhs < -primalTolerance)
    return true;

  return maxCoef > rhs + primalTolerance;
}

} // namespace

CoinBoundPropagation::CoinBoundPropagation(
  int numCols,
  const char *colType,
  const double *colLB,
  const double *colUB,
  const CoinPackedMatrix *matrixByRow,
  const char *sense,
  const double *rowRHS,
  const double *rowRange,
  double primalTolerance,
  double infinity,
  int maxRowNz)
  : newBounds_()
  , infeasible_(false)
  , infeasibleRow_(-1)
  , infeasibleCol_(-1)
  , complete_(true)
{
  // Mutable copies of column bounds updated as fixings are propagated.
  std::vector< double > mutableLB(colLB, colLB + numCols);
  std::vector< double > mutableUB(colUB, colUB + numCols);
  double *mColLB = mutableLB.data();
  double *mColUB = mutableUB.data();

  // fixedTo[j] == -1: not yet fixed, 0: fixed to 0, 1: fixed to 1.
  std::vector< int > fixedTo(static_cast< size_t >(numCols), -1);

  const int *idxs = matrixByRow->getIndices();
  const double *coefs = matrixByRow->getElements();
  const CoinBigIndex *start = matrixByRow->getVectorStarts();
  const int *length = matrixByRow->getVectorLengths();
  const size_t nRows = static_cast< size_t >(matrixByRow->getNumRows());

  CoinKnapsackRow knapsackRow(
    static_cast< size_t >(numCols),
    colType, mColLB, mColUB,
    primalTolerance, infinity);

  double multipliers[2];
  double rhsAdjustments[2];

#ifdef COIN_BT_STATS
  rowStats_.reserve(nRows);
#endif

  for (size_t idxRow = 0; idxRow < nRows; ++idxRow) {
    const char rowSense = sense[idxRow];
    const CoinBigIndex rowStart = start[idxRow];
    const size_t rowLength = static_cast< size_t >(length[idxRow]);
    const int *rowIdxs = idxs + rowStart;
    const double *rowCoefs = coefs + rowStart;
    const double range = rowRange ? rowRange[idxRow] : 0.0;

    // Hard skip: row exceeds the user-supplied nonzero limit.
    // This is a heuristic — fixings and infeasibility from this row are lost.
    if (maxRowNz >= 0 && static_cast< int >(rowLength) > maxRowNz) {
      complete_ = false;
#ifdef COIN_BT_STATS
      rowStats_.push_back({idxRow, 0.0, 0, true});
#endif
      continue;
    }

#ifdef COIN_BT_STATS
    const auto statsT0 = std::chrono::high_resolution_clock::now();
    const size_t fixingsBefore = newBounds_.size();
#endif

    const int numIter = CoinKnapsackRow::rowIterations(
      rowSense, rowRHS[idxRow], range, multipliers, rhsAdjustments);

    bool rowSkipped = (numIter == 0);

    for (int it = 0; it < numIter; ++it) {
      // Mathematical pre-check: skip this iteration if no fixing is possible.
      // This is always sound — it never suppresses a real fixing or infeasibility.
      if (!rowNeedsProcessing(rowIdxs, rowCoefs, rowLength,
                               multipliers[it], rhsAdjustments[it],
                               mColLB, mColUB, colType,
                               primalTolerance, infinity)) {
        rowSkipped = true;
        continue;
      }

      knapsackRow.processRow(
        rowIdxs, rowCoefs, rowLength, rowSense,
        multipliers[it], rhsAdjustments[it]);

      // Skip rows that are unbounded or have no binary variables.
      if (knapsackRow.isUnbounded()) {
        rowSkipped = true;
        continue;
      }

      rowSkipped = false;

      const double rhs = knapsackRow.rhs();

      // Direct row infeasibility: effective RHS is finite but strictly negative.
      if (std::isfinite(rhs) && rhs < -primalTolerance) {
#ifdef COIN_BT_STATS
        {
          const auto statsT1 = std::chrono::high_resolution_clock::now();
          rowStats_.push_back({idxRow,
            std::chrono::duration< double >(statsT1 - statsT0).count(),
            newBounds_.size() - fixingsBefore, false});
        }
#endif
        infeasibleRow_ = static_cast< int >(idxRow);
        infeasible_ = true;
        return;
      }

      const size_t nFixed = knapsackRow.nFixedVariables();
      if (nFixed == 0)
        continue;

      const int *fixedVars = knapsackRow.fixedVariables();
      for (size_t fi = 0; fi < nFixed; ++fi) {
        const int rawIdx = fixedVars[fi];
        // rawIdx < numCols  → original variable must be 0
        // rawIdx >= numCols → complemented, so original must be 1
        const int origCol = (rawIdx < numCols) ? rawIdx : rawIdx - numCols;
        const int newVal = (rawIdx < numCols) ? 0 : 1;

        if (fixedTo[origCol] == -1) {
          // First time this variable is fixed.
          fixedTo[origCol] = newVal;
          const double lb = static_cast< double >(newVal);
          const double ub = static_cast< double >(newVal);
          newBounds_.push_back(
            std::make_pair(static_cast< size_t >(origCol),
              std::make_pair(lb, ub)));
          mColLB[origCol] = lb;
          mColUB[origCol] = ub;
        } else if (fixedTo[origCol] != newVal) {
#ifdef COIN_BT_STATS
          {
            const auto statsT1 = std::chrono::high_resolution_clock::now();
            rowStats_.push_back({idxRow,
              std::chrono::duration< double >(statsT1 - statsT0).count(),
              newBounds_.size() - fixingsBefore, false});
          }
#endif
          // Contradictory fixing: same variable implied to be both 0 and 1.
          infeasibleCol_ = origCol;
          infeasibleRow_ = static_cast< int >(idxRow); // row that found the conflict
          infeasible_ = true;
          return;
        }
        // If fixedTo[origCol] == newVal the fixing is a duplicate; skip.
      }
    } // row iterations

#ifdef COIN_BT_STATS
    {
      const auto statsT1 = std::chrono::high_resolution_clock::now();
      rowStats_.push_back({idxRow,
        std::chrono::duration< double >(statsT1 - statsT0).count(),
        newBounds_.size() - fixingsBefore,
        rowSkipped});
    }
#endif
  } // all rows
}
