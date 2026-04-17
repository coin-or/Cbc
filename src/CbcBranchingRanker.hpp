// Copyright (C) 2025, Haroldo Gambini Santos and others.
// All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

#ifndef CbcBranchingRanker_H
#define CbcBranchingRanker_H

#include <cmath>
#include <cstddef>

#include "CbcConfig.h"

/** Configurable sort-key ranker for strong branching candidate prioritization.
 *
 *  In \c CbcNode::chooseDynamicBranch, every fractional variable receives a
 *  sort key that determines its position in the strong branching loop (R1
 *  ranking).  This class augments the default pseudo-cost-based sort key with
 *  additional structural information from the conflict graph.
 *
 *  ### Conflict Graph Degree Criterion (C4)
 *
 *  For each binary variable \c x_j, the conflict graph exposes two directional
 *  degrees (both O(1) lookups):
 *  - \c d1 = degree(j)            — conflicts when x_j = 1
 *  - \c d0 = degree(j + numCols)  — conflicts when x_j = 0
 *
 *  A high conflict degree means that branching on this variable triggers many
 *  propagations, reducing the feasible region faster.  The formula_ field
 *  controls how d0 and d1 are combined into a scalar score.
 *
 *  ### Trust-Sensitive Scaling
 *
 *  The influence of the conflict score depends on how much we trust the
 *  pseudo-cost estimate for the candidate:
 *
 *  - **Trusted** (enough branching observations): pseudo-cost score is
 *    reliable; conflict info acts as a gentle tie-breaker.  The score is
 *    raised to scalingPowerTrusted_ (default 0.5 = sqrt) before weighting,
 *    which compresses the range and reduces the boost magnitude.
 *
 *  - **Untrusted** (few or no observations): pseudo-cost score is weak;
 *    conflict info is allowed stronger influence.  The score is raised to
 *    scalingPowerUntrusted_ (default 1.0 = linear).
 *
 *  ### Backward Compatibility
 *
 *  When weightConflict_ == 0.0 (the default), applyConflictBoost() returns
 *  the sort key unchanged.  Default CBC behavior is preserved exactly.
 *
 *  ### Usage
 *
 *  Attach a ranker to a model via CbcModel::setBranchingRanker().
 *  Typical CLI usage:
 *  \code
 *    -rankConflict 0.05                        # mild boost, sqrt scaling when trusted
 *    -rankConflict 0.1 -rankConflictType sum   # use d0+d1 formula
 *    -rankConflict 0.05 -rankConflictPowerTrusted 0.333   # cbrt (very mild tie-breaker)
 *  \endcode
 */
class CBCLIB_EXPORT CbcBranchingRanker {
public:
  /** How d0 and d1 are combined into the raw conflict score. */
  enum ConflictScoreFormula {
    /** min(d0, d1) — both directions must be strong (default, conservative). */
    CONFLICT_MIN = 0,
    /** d0 + d1 — total propagation power across both directions. */
    CONFLICT_SUM = 1,
    /** sqrt(d0 * d1) — product score analog; rewards balanced high degrees. */
    CONFLICT_PRODUCT = 2
  };

  // --- Construction -------------------------------------------------------

  /** Default constructor — all weights zero, no behavioral change. */
  CbcBranchingRanker();

  /** Copy constructor. */
  CbcBranchingRanker(const CbcBranchingRanker &rhs);

  /** Assignment operator. */
  CbcBranchingRanker &operator=(const CbcBranchingRanker &rhs);

  ~CbcBranchingRanker() = default;

  // --- Core computation ---------------------------------------------------

  /** Compute the raw conflict score from directional degrees.
   *
   *  \param d0  Conflict graph degree of node (col + numCols), i.e. conflicts
   *             when x_col = 0.
   *  \param d1  Conflict graph degree of node col, i.e. conflicts when
   *             x_col = 1.
   *  \return    Non-negative raw score; 0 if both degrees are zero.
   */
  double conflictScore(std::size_t d0, std::size_t d1) const;

  /** Apply the conflict boost to an existing sort key.
   *
   *  Sort keys in chooseDynamicBranch are negative (more negative = higher
   *  priority).  This method multiplies the key by (1 + w * scaledScore),
   *  making high-conflict candidates more negative (higher priority).
   *
   *  \param sortKey  Current sort key (negative).
   *  \param d0       Conflict degree when x_col = 0.
   *  \param d1       Conflict degree when x_col = 1.
   *  \param trusted  True if pseudo-cost observations are sufficient (uses
   *                  scalingPowerTrusted_, default sqrt); false uses
   *                  scalingPowerUntrusted_ (default linear).
   *  \return         Boosted sort key, or sortKey unchanged when
   *                  weightConflict_ == 0 or conflictScore == 0.
   */
  double applyConflictBoost(double sortKey, std::size_t d0, std::size_t d1,
    bool trusted) const;

  /** Apply the column non-zeros boost to an existing sort key.
   *
   *  Score = nz^scalingPower (raw count raised to a slow-growing power).
   *  No normalization is needed: the 4th-root default (scalingPowerNzTrusted_ = 0.25)
   *  grows very slowly (nz=1→1.0, nz=16→2.0, nz=100→3.16), making this a
   *  cheap O(1) tie-breaker that does not require a max-scan over all columns.
   *
   *  \param sortKey  Current sort key (negative, possibly already boosted).
   *  \param nz       Number of non-zeros in this column (constraint count).
   *  \param trusted  True if pseudo-cost observations are sufficient.
   *  \return         Boosted sort key, or sortKey unchanged when weightNonzeros_ == 0.
   */
  double applyNonzerosBoost(double sortKey, int nz, bool trusted) const;

  /** Apply the variable-range boost to an existing sort key.
   *
   *  Score = 1 / (ub - lb).  For integer variables the minimum range is 1
   *  (binary), so the score is naturally in (0, 1].  High score = tight
   *  domain = closer to being fixed = deserves higher strong-branching priority.
   *
   *  \param sortKey  Current sort key (negative, possibly already conflict-boosted).
   *  \param lb       Current lower bound of the variable (saveLower[iColumn]).
   *  \param ub       Current upper bound of the variable (saveUpper[iColumn]).
   *  \param trusted  True if pseudo-cost observations are sufficient.
   *  \return         Boosted sort key, or sortKey unchanged when weightRange_ == 0.
   */
  double applyRangeBoost(double sortKey, double lb, double ub,
    bool trusted) const;

  // --- Diagnostics --------------------------------------------------------

  /** Human-readable name of the current formula ("min", "sum", "product"). */
  const char *formulaName() const;

  /** Returns true if any criterion is active (any weight > 0). */
  bool isActive() const
  {
    return weightConflict_ > 0.0 || weightRange_ > 0.0 || weightNonzeros_ > 0.0;
  }

  /** Reset cumulative diagnostic counters (e.g. between solves). */
  void resetCounters() const;

  // --- Parameters ---------------------------------------------------------

  /** Weight for the conflict degree criterion.
   *  Default: 0.0 (disabled).  Typical useful range: [0.01, 0.5].
   *  With weight = 0, the ranker has no effect (backward compatible). */
  double weightConflict_;

  /** Formula for combining d0 and d1 into a scalar score.
   *  Default: CONFLICT_MIN. */
  ConflictScoreFormula formula_;

  /** Scaling exponent applied to the conflict score when the pseudo-cost
   *  estimate is trusted (enough branching observations).
   *  Default: 0.5 (square root) — gentle tie-breaker.
   *  Set to 0.333 for cube root (even milder), or 1.0 for full influence. */
  double scalingPowerTrusted_;

  /** Scaling exponent applied to the conflict score when the pseudo-cost
   *  estimate is untrusted (few or no observations).
   *  Default: 1.0 (linear) — stronger influence among new variables.
   *  Set to 0.5 for sqrt (moderate). */
  double scalingPowerUntrusted_;

  /** Weight for the variable range criterion (1 / (ub - lb)).
   *  Applies to all integer variables (not just binary).  A variable with
   *  a small domain is closer to being fixed and prioritized accordingly.
   *  Default: 0.0 (disabled).  Typical useful range: [0.01, 0.5]. */
  double weightRange_;

  /** Scaling exponent for the range score when pseudo-costs are trusted.
   *  Default: 0.5 (sqrt) — gentle tie-breaker, like scalingPowerTrusted_. */
  double scalingPowerRangeTrusted_;

  /** Scaling exponent for the range score when pseudo-costs are untrusted.
   *  Default: 1.0 (linear) — full influence when pseudo-costs are weak. */
  double scalingPowerRangeUntrusted_;

  /** Weight for the column non-zeros criterion (nz / maxNz).
   *  A variable appearing in many constraints propagates fixing information
   *  more broadly and is generally a more impactful branching choice.
   *  Applies to all integer variables.  Score ∈ (0,1].
   *  Default: 0.0 (disabled).  Designed as a tie-breaker; typical range: [0.01, 0.1]. */
  double weightNonzeros_;

  /** Scaling exponent for the non-zeros score when pseudo-costs are trusted.
   *  Default: 0.25 (4th root) — very slow growth, cheap tie-breaker.
   *  nz=1→1.0, nz=16→2.0, nz=100→3.16.  No max-scan needed. */
  double scalingPowerNzTrusted_;

  /** Scaling exponent for the non-zeros score when pseudo-costs are untrusted.
   *  Default: 0.5 (sqrt) — moderate growth, more influence early in the tree. */
  double scalingPowerNzUntrusted_;

  // --- Diagnostic counters (mutable — updated by const methods) -----------

  /** Total number of binary variable candidates that received a non-zero
   *  conflict boost across all chooseDynamicBranch calls. */
  mutable long long nBoostsApplied_;

  /** Binary variable candidates skipped because their conflict score was
   *  zero (both directional degrees zero). */
  mutable long long nZeroScore_;

  /** Total number of integer variable candidates that received a non-zero
   *  range boost (score = 1/(ub-lb) > 0, which is always true for integers). */
  mutable long long nRangeBoostsApplied_;

  /** Total number of integer variable candidates that received a non-zeros boost. */
  mutable long long nNzBoostsApplied_;

  /** Set to true after the one-shot startup diagnostic message is printed. */
  mutable bool headerPrinted_;
};

#endif /* CbcBranchingRanker_H */

/* vi: softtabstop=2 shiftwidth=2 expandtab tabstop=2
*/
