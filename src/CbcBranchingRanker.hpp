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
};

#endif /* CbcBranchingRanker_H */

/* vi: softtabstop=2 shiftwidth=2 expandtab tabstop=2
*/
