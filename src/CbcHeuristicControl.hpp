// Copyright (C) COIN-OR contributors. All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

#ifndef CbcHeuristicControl_H
#define CbcHeuristicControl_H

#include "CbcConfig.h"
#include <string>

/** Node-count/depth-based scheduling policy for periodic re-invocation of a
    heuristic in the tree.
      - Legacy: the heuristic's own historical gating logic (default).
      - EveryKDepth: run only when the current depth is a multiple of K.
      - EveryKNodes: run once every K nodes, independent of depth.
      - EveryKNodesNoImprove: run every K nodes since the last incumbent
        improvement AND immediately after any improvement (which resets the
        K-node window).
    For heuristics that are triggered by a new incumbent (combine, crossover,
    proximity) any non-Legacy mode means "wait at least K nodes since the
    last run". See RINS-FIXTURES.md for the benchmarking behind the modes. */
enum class HeuristicScheduleMode {
  Legacy = 0,
  EveryKDepth = 1,
  EveryKNodes = 2,
  EveryKNodesNoImprove = 3
};

/** Everything that decides WHETHER a heuristic is attempted at a given point
    of the search, kept apart from the heuristics themselves so that the
    policy is the same for all of them, is inspectable, and can be extended
    in one place.

    The policy (what the user can set):
      - mode, everyK: the schedule (see HeuristicScheduleMode);
      - decay, window, minSuccess, maxK: back-off. After every `window` runs,
        if fewer than `minSuccess` of all runs so far improved the incumbent,
        everyK is multiplied by (1 + decay) (capped at maxK);
      - maxFailures: pause the heuristic after that many consecutive failed
        runs against the same incumbent; it is re-armed as soon as the
        incumbent changes (0 = never pause);
      - timeShare: do not start a run while the heuristic has already used
        more than this fraction of the elapsed solve time (0 = no limit).
    The state/statistics are counted for every real run (see RunGuard in
    CbcHeuristic). */
class CBCLIB_EXPORT CbcHeuristicControl {
public:
  CbcHeuristicControl() = default;

  /// \name Policy
  //@{
  HeuristicScheduleMode mode = HeuristicScheduleMode::Legacy;
  int everyK = 1;
  double decay = 0.0;
  int window = 10;
  double minSuccess = 1.0 / 3.0;
  int maxK = 1000000;
  int maxFailures = 0;
  double timeShare = 0.0;
  //@}

  /// \name Statistics (also the state of the policy)
  //@{
  int runs = 0; ///< completed real runs
  int successes = 0; ///< runs that improved the incumbent
  int failures = 0; ///< consecutive failed runs on the current incumbent
  int backoffs = 0; ///< times everyK was increased
  int skippedPaused = 0; ///< considered but paused by maxFailures
  int skippedTime = 0; ///< considered but over the time share
  int skippedSchedule = 0; ///< considered but not scheduled now
  double seconds = 0.0; ///< total time spent in real runs
  int initialK = 1; ///< everyK before any back-off
  //@}

  /** Cheap checks independent of the schedule: pause after repeated failures
      and time share. `minGapNodes` > 0 additionally requires that many nodes
      since the last run. Counts the reason it said no. */
  bool allows(int nodeCount, double elapsed, int solutionCount, bool useGap);

  /** The schedule decision for `mode` at this node. Has side effects (window
      bookkeeping): call at most once per candidate invocation. */
  bool scheduled(int nodeCount, int depth, int solutionCount);

  /// Account for a finished real run (time, outcome, back-off).
  void finishRun(bool success, double runSeconds, int nodeCount, int solutionCount);

  /// Account for time of an unfinished (skipped midway) attempt
  void addSeconds(double s) { seconds += s; }

  /// True if the heuristic is currently paused by maxFailures
  bool paused(int solutionCount) const
  {
    return maxFailures > 0 && failures >= maxFailures && failuresAt == solutionCount;
  }

  /// Short human readable policy summary
  std::string describe() const;
  static const char *modeName(HeuristicScheduleMode m);

private:
  int lastRunNode = -1000000;
  int lastSolutionCount = 0;
  int failuresAt = -1;
};

#endif
