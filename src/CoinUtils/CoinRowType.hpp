/* -*- mode: C++; tab-width: 2; indent-tabs-mode: nil; -*-
 *
 * @file   CoinRowType.hpp
 * @brief  Lightweight constraint-type classifier for profiling.
 *
 * Classifies a single constraint row into one of the standard MIP
 * row types (set packing, covering, knapsack, etc.) using the same
 * logic as OsiFeatures but without requiring an OsiSolverInterface.
 */

#ifndef COIN_ROW_TYPE_HPP
#define COIN_ROW_TYPE_HPP

#include <cmath>
#include <cstddef>
#include "CoinUtilsConfig.h"

enum CoinRowType {
  CoinRowPacking = 0,
  CoinRowPartitioning,
  CoinRowCovering,
  CoinRowCardinality,
  CoinRowInvKnapsack,
  CoinRowKnapsack,
  CoinRowIntKnapsack,
  CoinRowBinPacking,
  CoinRowFlowBin,
  CoinRowFlowMixed,
  CoinRowMixedBin,
  CoinRowGenInt,
  CoinRowSingleton,
  CoinRowOther,
  CoinRowTypeCount
};

inline const char *coinRowTypeName(CoinRowType t)
{
  static const char *names[] = {
    "Packing", "Partitioning", "Covering",
    "Cardinality", "InvKnapsack", "Knapsack", "IntKnapsack",
    "BinPacking", "FlowBin", "FlowMixed", "MixedBin",
    "GenInt", "Singleton", "Other"
  };
  return (t >= 0 && t < CoinRowTypeCount) ? names[t] : "Unknown";
}

/**
 * Classify a single constraint row.
 *
 * @param nz      number of non-zeros in the row
 * @param idx     column indices
 * @param coef    coefficients
 * @param sense   row sense ('L','G','E','R','N')
 * @param rhs     right-hand side
 * @param colType column types (CoinColumnType codes)
 * @param colLB   column lower bounds
 * @param colUB   column upper bounds
 */
inline CoinRowType classifyRow(
  int nz,
  const int *idx,
  const double *coef,
  char sense,
  double rhs,
  const char *colType,
  const double *colLB,
  const double *colUB)
{
  if (nz <= 0)
    return CoinRowOther;
  if (nz == 1)
    return CoinRowSingleton;

  int nBin = 0, nCont = 0, nInt = 0;
  double minV = coef[0], maxV = coef[0];
  int nPos = 0, nNeg = 0;
  bool allInt = true;

  for (int j = 0; j < nz; ++j) {
    double v = coef[j];
    if (v < minV) minV = v;
    if (v > maxV) maxV = v;
    if (v > 1e-16) ++nPos;
    else if (v < -1e-16) ++nNeg;
    if (allInt && fabs(v - round(v)) > 1e-16) allInt = false;

    int c = idx[j];
    char ct = colType[c];
    // binary: integer type with bounds in {0,1}
    if (ct == 1 && colLB[c] >= -1e-16 && colUB[c] <= 1.0 + 1e-16
      && (colUB[c] - colLB[c]) > 0.5)
      ++nBin;
    else if (ct == 2)
      ++nInt;
    else
      ++nCont;
  }

  bool rowBin = (nBin == nz);

  if (rowBin) {
    bool allPlusOnes = (fabs(minV - 1.0) <= 1e-16) && (fabs(maxV - 1.0) <= 1e-16);
    bool allMinusOnes = (fabs(minV + 1.0) <= 1e-16) && (fabs(maxV + 1.0) <= 1e-16);

    if (allPlusOnes || allMinusOnes) {
      char eSense = allMinusOnes
        ? (sense == 'L' ? 'G' : (sense == 'G' ? 'L' : sense))
        : sense;
      double eRhs = allMinusOnes ? -rhs : rhs;

      if (fabs(eRhs - 1.0) <= 1e-16) {
        if (eSense == 'E') return CoinRowPartitioning;
        if (eSense == 'L') return CoinRowPacking;
        if (eSense == 'G') return CoinRowCovering;
      }
      if (eRhs >= 1.99) {
        if (eSense == 'E') return CoinRowCardinality;
        if (eSense == 'L') return CoinRowInvKnapsack;
      }
    } else {
      if (rhs >= 1.1 && (maxV - minV >= 0.1) && nNeg == 0)
        return allInt ? CoinRowIntKnapsack : CoinRowKnapsack;
      if (rhs >= 1.1 && nNeg == 1 && nz >= 2)
        return CoinRowBinPacking;
    }

    if (nNeg >= 2 && nPos >= 2 && sense == 'E')
      return CoinRowFlowBin;
  } else {
    if (nNeg >= 2 && nPos >= 2 && sense == 'E')
      return CoinRowFlowMixed;
    if (nCont > 0)
      return CoinRowMixedBin;
    if (nInt > 0)
      return CoinRowGenInt;
  }

  return CoinRowOther;
}

#endif // COIN_ROW_TYPE_HPP
