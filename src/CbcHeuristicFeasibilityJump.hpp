// Copyright (C) 2024 MIPster contributors.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).
//
// Integrates the Feasibility Jump heuristic (MIT © 2022 SINTEF) into CBC/MIPster.
// Reference: Leitner, Fischetti, Toth (2023), "Feasibility Jump"
//   https://link.springer.com/article/10.1007/s12532-023-00234-8

#ifndef CbcHeuristicFeasibilityJump_H
#define CbcHeuristicFeasibilityJump_H

#include <cstdint>
#include "CbcHeuristic.hpp"

/** Feasibility Jump heuristic.
 *
 *  A primal heuristic that searches for integer-feasible solutions without
 *  requiring LP solves.  It maintains a weighted score over constraints and
 *  iteratively flips integer variables toward feasibility.
 *
 *  The algorithm is a port of the reference C++ implementation by SINTEF
 *  (MIT licence) into the CBC CbcHeuristic framework.
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

  /** Run the heuristic.  Returns 1 and fills newSolution/objectiveValue if a
   *  feasible integer solution is found; 0 otherwise. */
  using CbcHeuristic::solution;
  virtual int solution(double &objectiveValue, double *newSolution) override;

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

  /// Stop after finding this many feasible solutions. Default: 1.
  inline void setMaxSolutions(int n) { maxSolutions_ = n; }
  inline int maxSolutions() const { return maxSolutions_; }

  /// Feasibility tolerance for constraint violation (default: use solver's).
  inline void setFeasibilityTolerance(double t) { feasibilityTolerance_ = t; }
  inline double feasibilityTolerance() const { return feasibilityTolerance_; }

  /// Integer tolerance (default: use solver's).
  inline void setIntegerTolerance(double t) { integerTolerance_ = t; }
  inline double integerTolerance() const { return integerTolerance_; }

protected:
  bool relaxContinuous_ = false;
  int seed_ = 0;
  double weightUpdateDecay_ = 1.0;
  int64_t maxEffort_ = 0; // 0 = use NNZ-scaled
  int effortMultiplier_ = 1024;
  int stallMultiplier_ = 256;
  int maxSolutions_ = 1;
  double feasibilityTolerance_ = 1.0e-6;
  double integerTolerance_ = 1.0e-6;
};

#endif // CbcHeuristicFeasibilityJump_H
