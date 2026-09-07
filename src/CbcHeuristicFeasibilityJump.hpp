// Copyright (C) 2026, Haroldo Gambini Santos and others. All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).
//
// Integrates the Feasibility Jump heuristic (MIT (c) 2022 SINTEF) into Cbc.
// Reference: Leitner, Fischetti, Toth (2023), "Feasibility Jump"
//   https://link.springer.com/article/10.1007/s12532-023-00234-8

#ifndef CbcHeuristicFeasibilityJump_H
#define CbcHeuristicFeasibilityJump_H

#include <cstdint>
#include "CbcHeuristic.hpp"

/** Feasibility Jump heuristic.
 *
 *  A primal heuristic that searches for integer-feasible solutions without
 *  requiring LP solves. It maintains a weighted score over constraints and
 *  iteratively flips integer variables toward feasibility.
 *
 *  Root-node processing (pre-processing, cut generation, heuristics) is
 *  where CBC establishes both its dual bound (from pre-processing and cuts)
 *  and its primal bound (from heuristics). Feasibility Jump is LP-free and
 *  fast, so it is attractive to run it repeatedly at several points during
 *  root processing:
 *
 *    a) on the very first fractional solution of the pre-processed problem,
 *       before any cut-generation round has run (see setRunBeforeFirstCutRound);
 *    b) after each round of root cut generation, seeded from the round's
 *       optimal LP basis (the default whereFrom=1/2 hooks used by all Cbc
 *       heuristics already provide this, once per pass after the first);
 *    c) periodically inside the tree, every N levels of depth (setMinDepth);
 *
 *  In all cases, by default the heuristic only runs while CBC still has no
 *  incumbent solution at all (setOnlyIfNoIncumbent) -- once *some* feasible
 *  solution exists, further FJ calls are of much more limited value and just
 *  add overhead, so this is skippable to save time for other cut/heuristic
 *  work. Getting *some* incumbent as early/reliably as possible is itself an
 *  important metric: without any primal bound, other parts of the search
 *  (e.g. reduced-cost fixing, best-first node selection) cannot help at all.
 *
 *  A configurable cap on the number of separate FJ invocations across the
 *  whole solve (setMaxCalls) lets experiments trade off calling FJ fewer
 *  times with a larger iteration budget each vs. calling it more often with
 *  a smaller budget each.
 */
class CBCLIB_EXPORT CbcHeuristicFeasibilityJump : public CbcHeuristic {
public:
  CbcHeuristicFeasibilityJump();
  CbcHeuristicFeasibilityJump(CbcModel &model);
  CbcHeuristicFeasibilityJump(const CbcHeuristicFeasibilityJump &);
  ~CbcHeuristicFeasibilityJump();

  virtual CbcHeuristic *clone() const override;
  CbcHeuristicFeasibilityJump &operator=(const CbcHeuristicFeasibilityJump &);

  virtual void resetModel(CbcModel *model) override;
  virtual void setModel(CbcModel *model) override;

  /// Override to enable tree execution when minDepth > 0.
  virtual bool shouldHeurRun(int whereFrom) override;

  /** Run the heuristic.  Returns 1 and fills newSolution/objectiveValue if a
   *  feasible integer solution is found; 0 otherwise. */
  using CbcHeuristic::solution;
  virtual int solution(double &objectiveValue, double *newSolution) override;

  /** Entry point used to run FJ seeded from an externally supplied point
   *  (rather than the current LP relaxation solution), e.g. Feasibility
   *  Pump's last rounded-but-infeasible attempt when it fails to find a
   *  solution (point (d) of the FJ integration plan -- see
   *  CbcHeuristicFPump::setFeasibilityJumpFallback()). Bypasses
   *  shouldHeurRun()'s throttling (this is a one-off, event-triggered call,
   *  not part of the normal per-round schedule) but still honours
   *  onlyIfNoIncumbent_/maxCalls_. Returns 1 and fills
   *  objectiveValue/newSolution on success, exactly like solution(). */
  int solveFromSeed(double &objectiveValue, double *newSolution,
    const double *seedSolution);

  /// Whether to relax continuous variables (treat them at their bounds).
  inline void setRelaxContinuous(bool val) { relaxContinuous_ = val; }
  inline bool relaxContinuous() const { return relaxContinuous_; }

  /// Random seed for the internal PRNG.
  inline void setSeed(int seed) { seed_ = seed; }
  inline int seed() const { return seed_; }

  /// Weight-update decay parameter (default 1.0 = no exponential decay).
  inline void setWeightUpdateDecay(double d) { weightUpdateDecay_ = d; }
  inline double weightUpdateDecay() const { return weightUpdateDecay_; }

  /// Maximum effort (deterministic iteration budget) for one heuristic call.
  /// If positive, used as a fixed budget. If zero (default), the budget is
  /// computed as NNZ * effortMultiplier_, scaling with problem size.
  /// Default: 0 (use NNZ-scaled effort).
  inline void setMaxEffort(int64_t e) { maxEffort_ = e; }
  inline int64_t maxEffort() const { return maxEffort_; }

  /// Multiplier for NNZ-scaled effort budget (default 1024, like HiGHS).
  /// Only used when maxEffort_ == 0.
  inline void setEffortMultiplier(int m) { effortMultiplier_ = m; }
  inline int effortMultiplier() const { return effortMultiplier_; }

  /// Stall limit: terminate when effortSinceLastImprovement exceeds
  /// NNZ * stallMultiplier_. Default 256 (like HiGHS). Set to 0 to disable.
  inline void setStallMultiplier(int m) { stallMultiplier_ = m; }
  inline int stallMultiplier() const { return stallMultiplier_; }

  /// Run FJ every N levels in the tree. Default: 0 (root only).
  /// Set to e.g. 6 to run FJ at depths 6, 12, 18... (like bound propagation).
  inline void setMinDepth(int d) { minDepth_ = d; }
  inline int minDepth() const { return minDepth_; }

  /// Stop after finding this many feasible solutions in a single call.
  /// Default: 1.
  inline void setMaxSolutions(int n) { maxSolutions_ = n; }
  inline int maxSolutions() const { return maxSolutions_; }

  /// Feasibility tolerance for constraint violation (default: use solver's).
  inline void setFeasibilityTolerance(double t) { feasibilityTolerance_ = t; }
  inline double feasibilityTolerance() const { return feasibilityTolerance_; }

  /// Integer tolerance (default: use solver's).
  inline void setIntegerTolerance(double t) { integerTolerance_ = t; }
  inline double integerTolerance() const { return integerTolerance_; }

  /// Whether to skip running FJ entirely once CBC already has at least one
  /// incumbent solution (of any origin: heuristic, MIP start, or B&B).
  /// Default: true, per the observation that FJ is most valuable for
  /// producing the very first incumbent; once one exists, repeated FJ calls
  /// mostly just add overhead relative to other root/tree work. Set false
  /// to also let FJ try to improve on an existing incumbent.
  inline void setOnlyIfNoIncumbent(bool val) { onlyIfNoIncumbent_ = val; }
  inline bool onlyIfNoIncumbent() const { return onlyIfNoIncumbent_; }

  /// Caps the total number of separate FJ invocations (across all trigger
  /// points: before-first-cut-round, root-after-cuts, and tree) for the
  /// whole solve. Default: 0 (unlimited). Each invocation is always seeded
  /// from a genuinely new fractional solution (the round's new LP optimum,
  /// or a different tree node) -- FJ is never simply re-run on an unchanged
  /// point, since a repeat run from the same seed/PRNG state would just
  /// retrace the same trajectory. Use this together with
  /// maxEffort_/effortMultiplier_ to explore the tradeoff between calling FJ
  /// fewer times (at fewer of these distinct points) with a bigger iteration
  /// budget each vs. more times (at more of these points) with a smaller
  /// budget each.
  inline void setMaxCalls(int n) { maxCalls_ = n; }
  inline int maxCalls() const { return maxCalls_; }

  /// Number of times solution() has actually run the FJ local search so far
  /// (i.e. wasn't skipped by shouldHeurRun(), onlyIfNoIncumbent_, or
  /// maxCalls_).
  inline int callsMade() const { return callsMade_; }

protected:
  /// Shared implementation for solution().
  int solveFJ(double &objectiveValue, double *newSolution, int depth,
    const double *seedSolution = nullptr);

  bool relaxContinuous_ = false;
  int seed_ = 0;
  double weightUpdateDecay_ = 1.0;
  int64_t maxEffort_ = 0; // 0 = use NNZ-scaled
  int effortMultiplier_ = 1024;
  int stallMultiplier_ = 256;
  int minDepth_ = 0; // 0 = root only
  int maxSolutions_ = 1;
  double feasibilityTolerance_ = 1.0e-6;
  double integerTolerance_ = 1.0e-6;
  bool onlyIfNoIncumbent_ = true;
  int maxCalls_ = 0; // 0 = unlimited
  int callsMade_ = 0;
};

#endif // CbcHeuristicFeasibilityJump_H
