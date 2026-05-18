// Copyright (C) 2026, Haroldo Gambini Santos and others.
// All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

#ifndef ClpRacingSolver_H
#define ClpRacingSolver_H

#include "ClpConfig.h"
#include "ClpSimplex.hpp"
#include <vector>

/** Opportunistic parallel LP solver.
 *
 *  Races multiple ClpSolve configurations in parallel on clones of the
 *  model.  The first to reach optimality wins; all others are aborted
 *  via ClpAbortHandler.  The winner's solution and basis are copied back
 *  to the original model.
 *
 *  Intended for root-node LP relaxation where the best method is unknown
 *  a priori and different configurations (dual, primal+idiot, primal+sprint)
 *  can have dramatically different performance.
 *
 *  The caller should presolve the model before racing (presolve once, race
 *  on the presolved model).
 */
class CLPLIB_EXPORT ClpRacingSolver {
public:
  /** Construct with the model to solve and number of racing threads.
   *  If numThreads <= 0, uses the number of configs added. */
  ClpRacingSolver(ClpSimplex *model, int numThreads = 0);

  /// Add a configuration to race.
  void addConfig(const ClpSolve &config);

  /** Add the default racing portfolio (dual, primal+idiot, primal+sprint).
   *  Clears any previously added configs. */
  void addDefaultConfigs();

  /** Run the race.  Returns index of winning config (0-based), or -1 if
   *  all configs failed.  On success the model's solution, basis, and
   *  status are updated to the winner's result. */
  int solve();

  /// Index of the winning configuration after solve(), -1 if none.
  int winnerIndex() const { return winnerIndex_; }

  /// Wall-clock seconds taken by the winner.
  double winnerTime() const { return winnerTime_; }

  /// Iterations performed by the winner.
  int winnerIterations() const { return winnerIterations_; }

private:
  ClpSimplex *model_;
  std::vector<ClpSolve> configs_;
  int numThreads_;
  int winnerIndex_ = -1;
  double winnerTime_ = 0.0;
  int winnerIterations_ = 0;
};

#endif

/* vi: softtabstop=2 shiftwidth=2 expandtab tabstop=2
*/
