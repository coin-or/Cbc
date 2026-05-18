/* -*- mode: C++; tab-width: 2; indent-tabs-mode: nil; -*-
 *
 * This file is part of the COIN-OR CoinUtils package
 *
 * @file   CoinKnapsackRow.hpp
 * @brief  Helper class to analyze a single constraint and store it in a knapack form.
 *
 * Copyright (C) 2025
 * All rights reserved.
 *
 * This code is licensed under the terms of the Eclipse Public License (EPL).
 */

#ifndef COIN_KNAPSACK_ROW_HPP
#define COIN_KNAPSACK_ROW_HPP

#include <cmath>
#include <cstddef>

#include "CoinUtilsConfig.h"
#include "CoinColumnType.hpp"
#include "CoinTerm.hpp"

/**
 * @brief Helper that transforms a constraint row into the knapsack
 *        representation.
 *
 * The class collects all binary variables (including complemented versions),
 * updates the effective right-hand side, and tracks simple statistics such as
 * the two largest/smallest coefficients and variables that should be fixed due
 * to their coefficients exceeding the row bound.
 */
class COINUTILSLIB_EXPORT CoinKnapsackRow {
public:

  /**
   * @brief Construct a knapsack helper bound to the global column metadata.
   *
   * The constructor does not copy column data; it only keeps the provided
   * pointers so that repeated calls to #processRow can access the column types
   * and bounds. The caller must ensure that the pointed-to arrays remain valid
   * for the lifetime of the helper.
   *
   * @param numCols Total number of structural columns in the model; used to
   *        validate indexes passed to #processRow.
   * @param colType Pointer to an array of length @p numCols describing the type
   *        (binary, integer, continuous, etc.) of each column.
   * @param colLB Pointer to an array with the column lower bounds (length
   *        @p numCols).
   * @param colUB Pointer to an array with the column upper bounds (length
   *        @p numCols).
   * @param primalTolerance Numerical tolerance used to decide whether constraints/bounds
   * are violated. Defaults to $1\times
   *        10^{-7}$.
   * @param infinity Numerical value representing practical infinity when
   *        checking row feasibility. Defaults to $1\times10^{50}$.
   */
  CoinKnapsackRow(size_t numCols,
                 const char *colType,
                 const double *colLB,
                 const double *colUB,
                 double primalTolerance = 1e-7,
                 double infinity = 1e50);

  ~CoinKnapsackRow();

  CoinKnapsackRow(const CoinKnapsackRow &) = delete;
  CoinKnapsackRow &operator=(const CoinKnapsackRow &) = delete;
  CoinKnapsackRow(CoinKnapsackRow &&) = delete;
  CoinKnapsackRow &operator=(CoinKnapsackRow &&) = delete;

  /**
   * @brief Analyze a constraint row and populate the knapsack representation.
   *
   * Column metadata (types and bounds) must be provided to the object at
   * construction time; this method only receives row-specific data.
   *
   * @param idx Column indexes that appear in the row.
   * @param coef Coefficients corresponding to @p rowIdx.
   * @param nz Number of non-zeros in the row (length of @p rowIdx/@p rowCoef).
   * @param sense Sense of the row (L, G, E, R, N).
   * @param multiplier Multiplier applied to the row to obtain the working
   *        knapsack form. For range/equality constraints this should be provided
   *        by the caller.
   * @param rhs Original right-hand side value.
   */
  void processRow(const int *idx,
                  const double *coef,
                  size_t nz,
                  char sense,
                  double multiplier,
                  double rhs);

  /** @brief Pointer to the processed binary/complemented entries. */
  const CoinTerm *columns() const { return columns_; }

  /** @brief Number of processed entries in #columns(). */
  size_t nzs() const { return nzs_; }

  /** @brief Return the effective right-hand side after processing.
   * If the row should be ignored, returns NaN.
   */
  double rhs() const { return rhs_; }

  /** @brief Two largest coefficients seen during processing. */
  const double *twoLargest() const { return twoLargest_; }

  /** @brief Two smallest coefficients seen during processing. */
  const double *twoSmallest() const { return twoSmallest_; }

  /** @brief Variables that should be fixed due to coefficient > RHS. */
  const int *fixedVariables() const { return fixedVariables_; }

  /** @brief Number of stored fixings. */
  size_t nFixedVariables() const { return nFixed_; }

  /**
   * @brief Copy the indexes of the stored columns into @p dest.
   * @param dest Output buffer with room for at least #nzs() entries.
   */
  void copyColumnIndices(size_t *dest) const;

  /** @brief If row is unbounded or contain variable types not handled, i.e. semicontinuous. */
  bool isUnbounded() const { return (!std::isfinite(rhs_)) || (rhs_ >= infinity_); }

  /** @brief Sort stored entries by coefficient (and index ties). */
  void sortColumns();

  /**
   * @brief Check whether the processed row is already an explicit clique.
   *
   * if the sum of the two smallest (nonnegative) coefficients
   * exceeds the adjusted right-hand side (plus the primal tolerance), then
   * every pair of variables violates the row and the whole set forms a
   * clique.
   *
   * @return `true` when the stored knapsack representation corresponds to a
   *         clique (nzs >= 2, finite RHS, nonnegative RHS, and
   *         \(a_{(1)} + a_{(2)} > b + \varepsilon\)).
   */
  bool isExplicitClique() const;

  /**
   * @brief Provide the multipliers/RHS pairs needed to materialize a row.
   *
   * Range and equality rows induce two knapsack evaluations: the original
   * row and its complemented counterpart. This helper mirrors the logic used
   * in `CoinDynamicConflictGraph` by returning the multipliers that must be
   * applied to the row (\(+1\) or \(-1\)) and the corresponding adjusted
   * right-hand sides so that callers can feed each variant into
   * #processRow.
   *
   * @param sense Row sense (L, G, E, or R). Other senses return zero.
   * @param rhs Original right-hand side of the row.
   * @param range Range value when @p sense == 'R'; ignored otherwise.
   * @param multiplier Output array (length at least 2) filled with the
   *        multipliers to apply to the row.
   * @param adjustedRHS Output array (length at least 2) with the rhs values
   *        corresponding to each multiplier.
   * @return Number of entries written to @p multiplier/@p adjustedRHS (0, 1,
   *         or 2).
   */
  static int rowIterations(
    const char sense,
    const double rhs,
    const double range,
    double multiplier[],
    double adjustedRHS[]);

private:
  void resetComputedValues();

  const char *colType_;
  const double *colLB_;
  const double *colUB_;

  CoinTerm *columns_;
  size_t nzs_;
  size_t numCols_;

  int *fixedVariables_;
  size_t nFixed_;
  double rhs_;
  double twoLargest_[2];
  double twoSmallest_[2];
  const double primalTolerance_;
  const double infinity_;
};

#endif // COIN_KNAPSACK_ROW_HPP
