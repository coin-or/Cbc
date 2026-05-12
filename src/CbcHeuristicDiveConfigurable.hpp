// Copyright (C) 2026, Haroldo Gambini Santos and others.
// All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

#ifndef CbcHeuristicDiveConfigurable_H
#define CbcHeuristicDiveConfigurable_H

#include "CbcHeuristicDive.hpp"

/** Configurable diving heuristic with weighted multi-criteria variable selection.
 *
 *  Instead of a fixed selection criterion (like DiveCoefficient or DiveFractional),
 *  this class scores each fractional variable using a weighted combination of:
 *
 *  - **Fractionality**: how far from integrality (higher = more fractional)
 *  - **Locks**: number of constraints blocking rounding (fewer = easier to fix)
 *  - **Conflict degree**: from the conflict graph — high-degree variables propagate
 *    more fixings when rounded
 *  - **Objective coefficient**: variables with large |c_j| impact the bound more
 *  - **Column nonzeros**: variables in many constraints propagate information broadly
 *
 *  Each criterion uses a configurable weight and scaling power (following the
 *  CbcBranchingRanker pattern from the conflict-graph-branching-ranker branch).
 *
 *  The rounding direction can optionally use the conflict graph to pick the
 *  direction with fewer conflicts against already-fixed variables.
 *
 *  A minimum fractionality threshold prevents wasting iterations on near-integer
 *  variables (the main failure mode of DiveFractional and DiveGuided at root).
 *
 *  Multiple instances with different configurations can be added to the model
 *  to form a diving schedule. Use when_=4 on expensive configurations so they
 *  are skipped if an earlier, cheaper dive already found a solution.
 */
class CBCLIB_EXPORT CbcHeuristicDiveConfigurable : public CbcHeuristicDive {
public:
  CbcHeuristicDiveConfigurable();
  CbcHeuristicDiveConfigurable(CbcModel &model);
  CbcHeuristicDiveConfigurable(const CbcHeuristicDiveConfigurable &rhs);
  ~CbcHeuristicDiveConfigurable();

  CbcHeuristicDiveConfigurable *clone() const override;
  CbcHeuristicDiveConfigurable &operator=(const CbcHeuristicDiveConfigurable &rhs);

  void generateCpp(FILE *fp) override;

  bool selectVariableToBranch(OsiSolverInterface *solver,
    const double *newSolution,
    int &bestColumn,
    int &bestRound) override;

  void initializeData() override;

  // --- Configuration setters ---

  /// Minimum fractionality to consider a variable (default 0.01)
  void setMinFractionality(double v) { minFractionality_ = v; }
  double minFractionality() const { return minFractionality_; }

  /// Weight for fractionality criterion (default 1.0)
  void setWeightFractionality(double v) { weightFractionality_ = v; }
  /// Power for fractionality (default 1.0 = linear)
  void setPowerFractionality(double v) { powerFractionality_ = v; }

  /// Weight for locks criterion (default 0.1)
  void setWeightLocks(double v) { weightLocks_ = v; }
  /// Power for locks (default 0.5 = sqrt)
  void setPowerLocks(double v) { powerLocks_ = v; }

  /// Weight for conflict degree (default 0.0 = disabled)
  void setWeightConflict(double v) { weightConflict_ = v; }
  /// Power for conflict score (default 0.5 = sqrt)
  void setPowerConflict(double v) { powerConflict_ = v; }

  /// Weight for objective coefficient (default 0.0 = disabled)
  void setWeightObjCoeff(double v) { weightObjCoeff_ = v; }
  /// Power for objective coefficient (default 0.2)
  void setPowerObjCoeff(double v) { powerObjCoeff_ = v; }

  /// Weight for column nonzeros (default 0.0 = disabled)
  void setWeightNonzeros(double v) { weightNonzeros_ = v; }
  /// Power for nonzeros (default 0.25 = 4th root)
  void setPowerNonzeros(double v) { powerNonzeros_ = v; }

  /// Use conflict graph for rounding direction (default true)
  void setUseConflictForDirection(bool v) { useConflictForDirection_ = v; }

  /// If true, prefer rounding toward objective improvement (default false)
  void setPreferObjectiveDirection(bool v) { preferObjectiveDirection_ = v; }

  /// Random noise factor (0.0 = deterministic, 0.5 = up to 50% noise on scores)
  void setRandomFactor(double v) { randomFactor_ = v; }

  /// Target fractionality (default 0.5 = most fractional; 0.6 = biased toward rounding up)
  /// Score = 1 - |frac - target| so vars closest to target score highest
  void setTargetFractionality(double v) { targetFractionality_ = v; }

  /// Fix a fixed number of at-bound vars per iteration (0 = use percentageToFix instead)
  /// When > 0, overrides percentageToFix. Useful for propagation-heavy mode.
  void setFixCount(int v) { fixCount_ = v; }

private:
  // --- Parameters ---
  double minFractionality_;
  double weightFractionality_;
  double powerFractionality_;
  double weightLocks_;
  double powerLocks_;
  double weightConflict_;
  double powerConflict_;
  double weightObjCoeff_;
  double powerObjCoeff_;
  double weightNonzeros_;
  double powerNonzeros_;
  bool useConflictForDirection_;
  bool preferObjectiveDirection_;
  double randomFactor_;
  double targetFractionality_;
  int fixCount_;
};

#endif

/* vi: softtabstop=2 shiftwidth=2 expandtab tabstop=2
*/
