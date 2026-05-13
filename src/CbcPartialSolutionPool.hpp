// Copyright (C) 2026, Haroldo Gambini Santos and others.
// All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

#ifndef CbcPartialSolutionPool_H
#define CbcPartialSolutionPool_H

#include "CbcConfig.h"
#include "CoinWarmStart.hpp"
#include <memory>
#include <vector>

/** A variable fixing: column index + bounds to impose. */
struct CbcFixing {
  int col;
  double lb, ub;
};

/** A partial solution: a set of variable fixings with an LP basis
 *  and a quality score indicating how close it is to integrality.
 *
 *  Score = fraction of integer variables that are integer-feasible
 *  in the LP solution after applying the fixings. Higher is better.
 */
struct CbcPartialSolution {
  std::vector<CbcFixing> fixings;
  std::unique_ptr<CoinWarmStart> basis;
  double score; ///< fraction of integers satisfied (0..1)

  CbcPartialSolution() : score(0) {}
  CbcPartialSolution(CbcPartialSolution &&) = default;
  CbcPartialSolution &operator=(CbcPartialSolution &&) = default;

  int depth() const { return static_cast<int>(fixings.size()); }
};

/** Pool of partial solutions from diving heuristics.
 *
 *  Maintains a bounded set of the best partial solutions seen,
 *  ranked by score. When full, a new entry replaces the worst
 *  only if it scores higher.
 *
 *  Thread-safety: not thread-safe. Use external synchronization
 *  if accessed from multiple threads.
 */
class CBCLIB_EXPORT CbcPartialSolutionPool {
public:
  /** Construct with maximum pool capacity. */
  explicit CbcPartialSolutionPool(int maxSize = 32);

  /** Try to add a partial solution. Returns true if accepted.
   *  Rejected if pool is full and score <= worst entry's score. */
  bool add(std::vector<CbcFixing> &&fixings, CoinWarmStart *basis, double score);

  /** Number of entries currently in the pool. */
  int size() const { return static_cast<int>(pool_.size()); }

  /** Whether the pool is empty. */
  bool empty() const { return pool_.empty(); }

  /** Get entry by index. */
  const CbcPartialSolution &get(int i) const { return pool_[i]; }

  /** Select an entry for restart. Picks from the top half by score,
   *  using seed for deterministic selection. */
  const CbcPartialSolution &select(int seed) const;

  /** Score of the worst entry (minimum threshold for acceptance). */
  double worstScore() const;

  /** Clear all entries. */
  void clear() { pool_.clear(); }

  /** Maximum capacity. */
  int maxSize() const { return maxSize_; }

private:
  int maxSize_;
  std::vector<CbcPartialSolution> pool_;
};

#endif
