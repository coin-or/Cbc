// Copyright (C) 2026, Haroldo Gambini Santos and others.
// All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

#ifndef CbcDiveScheduler_H
#define CbcDiveScheduler_H

#include "CbcConfig.h"
#include <atomic>
#include <vector>

class CbcModel;
class CbcHeuristicDive;

/** Orchestrates a schedule of diving heuristics at the root node.
 *
 *  Runs a sequence of diving configurations, stopping early when a feasible
 *  solution is found. Supports both sequential and parallel execution.
 *
 *  In parallel mode (the default when threads > 1), all heuristics run
 *  concurrently and are terminated as soon as any one finds a solution.
 *
 *  Usage:
 *  \code
 *    CbcDiveScheduler scheduler(model);
 *    scheduler.addDefaultSchedule();  // adds the recommended configs
 *    int status = scheduler.run();    // 1 = solution found
 *  \endcode
 */
class CBCLIB_EXPORT CbcDiveScheduler {
public:
  CbcDiveScheduler(CbcModel &model);
  ~CbcDiveScheduler();

  /// Add a diving heuristic to the schedule (takes ownership)
  void addHeuristic(CbcHeuristicDive *heuristic);

  /// Add the recommended default schedule (Frac_A, Rand_Frac_NoFix, etc.)
  void addDefaultSchedule();

  /// Set number of threads (0 = auto from model, 1 = sequential)
  void setNumThreads(int n) { numThreads_ = n; }

  /// When true (default), cancel remaining dives once a solution is found
  void setStopOnFirstSolution(bool v) { stopOnFirst_ = v; }

  /** Run the schedule.
   *  \return 1 if a feasible solution was found, 0 otherwise.
   *  The best solution is stored in the model via setBestSolution().
   */
  int run();

  /// Best objective found (valid after run() returns 1)
  double bestObjective() const { return bestObj_; }

  /// Best solution vector (valid after run() returns 1)
  const double *bestSolution() const { return bestSolution_.empty() ? nullptr : bestSolution_.data(); }

  /// Which heuristic found the solution (index in schedule, -1 if none)
  int winnerIndex() const { return winnerIdx_; }

  /// Name of the winning heuristic
  const char *winnerName() const;

  /// Time spent in run()
  double timeUsed() const { return timeUsed_; }

private:
  int runSequential();
  int runParallel();

  CbcModel &model_;
  std::vector<CbcHeuristicDive *> schedule_;
  int numThreads_;
  bool stopOnFirst_;
  double bestObj_;
  int winnerIdx_;
  double timeUsed_;
  std::vector<double> bestSolution_;

  /// Shared flag for parallel cancellation
  std::atomic<bool> solutionFound_;
};

#endif

/* vi: softtabstop=2 shiftwidth=2 expandtab tabstop=2
*/
