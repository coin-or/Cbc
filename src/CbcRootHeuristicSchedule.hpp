// Copyright (C) 2026, Haroldo Gambini Santos and others.
// All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

#ifndef CbcRootHeuristicSchedule_H
#define CbcRootHeuristicSchedule_H

#include <atomic>
#include <vector>
#include "CbcConfig.h"

class CbcModel;
class CbcHeuristic;

/** Two-phase parallel heuristic schedule for root node.
 *
 *  Phase 1 (constructive): runs diving, FPump, rounding, greedy in parallel.
 *  Stops as soon as maxSolutionsPhase1 solutions are found.
 *
 *  Phase 2 (improvement): runs RINS, DINS, Combine, Proximity in parallel.
 *  Only executes if Phase 1 found at least one solution.
 *
 *  Uses std::atomic abort flags and ClpAbortHandler for fast cancellation.
 */
class CBCLIB_EXPORT CbcRootHeuristicSchedule {
public:
  CbcRootHeuristicSchedule(CbcModel &model);
  ~CbcRootHeuristicSchedule();

  /// Run the two-phase schedule. Returns number of solutions found.
  int run();

  /// Add the default v5 diving schedule to the model
  void addDefaultDivingConfigs();

  /// Parameters
  void setMaxSolutionsPhase1(int n) { maxSolutionsPhase1_ = n; }
  void setNumThreads(int n) { numThreads_ = n; }

  int maxSolutionsPhase1() const { return maxSolutionsPhase1_; }
  int numThreads() const { return numThreads_; }

private:
  /// Run constructive heuristics. Returns number of solutions found.
  int runPhase1(std::vector<CbcHeuristic *> &heuristics);

  /// Run improvement heuristics. Returns number of solutions found.
  int runPhase2(std::vector<CbcHeuristic *> &heuristics);

  /// Run a set of heuristics in parallel with early termination.
  /// Returns number of solutions found.
  int runParallel(std::vector<CbcHeuristic *> &heuristics,
    int maxSolutions);

  CbcModel &model_;
  std::atomic<bool> stopFlag_;
  int maxSolutionsPhase1_;
  int numThreads_;
  int solutionsFound_;
};

#endif
