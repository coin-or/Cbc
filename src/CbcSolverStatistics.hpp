/**
 * Copyright COIN-OR Foundation (C) 2025
 * All Rights Reserved.
 * This file is distributed under the Eclipse Public License.
 *
 * Purpose: Store several metrics related to the performance and results
 * of the CBC MIP solver.
 *
 **/


#ifndef CBC_SOLVER_STATISTICS
#define CBC_SOLVER_STATISTICS

#include "CbcConfig.h"
#include <deque>
#include <string>
#include <vector>
#include "CbcParameters.hpp"

/** Statistics for a single cut generator collected during a MIP solve.
 *
 *  Note: in multi-thread B&B (-threads N) child models hold independent
 *  clones of each generator and their stats are NOT merged back to the
 *  master.  The fields here therefore reflect only master-thread activity
 *  in that case — the same limitation that exists for all other per-generator
 *  counters (numberCutsInTotal, etc.).  For the racing-root-LP case the
 *  stats ARE aggregated via CbcCutGenerator::addStatistics.
 */
struct CBCLIB_EXPORT CutGeneratorStats {
  std::string name;
  int nCuts = 0;        ///< Total row cuts added
  int nCalls = 0;       ///< Number of times the generator was invoked
  double time = 0.0;    ///< Time spent in this generator (seconds)
  int nColumnCuts = 0;  ///< Total column cuts added
  int minDepth = -1;    ///< Minimum B&B depth at which generator was called (-1 if never called)
  int maxDepth = -1;    ///< Maximum B&B depth at which generator was called (-1 if never called)
};

/** Statistics for a single heuristic collected during a MIP solve.
 *
 *  Same threading caveat as CutGeneratorStats: in multi-thread B&B mode
 *  the heuristics on child thread models are cloned and their stats are
 *  never merged back, so these fields reflect only master-thread activity.
 */
struct CBCLIB_EXPORT HeuristicStats {
  std::string name;
  int nExecutions = 0;  ///< Number of times the heuristic was executed
  double totalTime = 0.0; ///< Total time spent in the heuristic (seconds)
  int nSolutions = 0;   ///< Number of solutions found
  int minDepth = -1;    ///< Minimum B&B depth at which heuristic was applied (-1 if never run)
  int maxDepth = -1;    ///< Maximum B&B depth at which heuristic was applied (-1 if never run)
};

class CBCLIB_EXPORT CbcSolverStatistics {
public:
  /** Elapsed total time */
  double seconds = 0.0;

  /** Best solution found in the search */
  double obj = 0.0;

  /** True if an integer feasible solution was found during the search */
  bool integer_feasible = false;

  /** CPU time */
  double sys_seconds = 0.0;

  /** Elapsed time from solver start */
  double elapsed_seconds = 0.0;

  /** LP relaxation cost */
  double continuous = 0.0;

  /** LP relaxation time */
  double lp_seconds = 0.0;

  /** Cost after tightening LP relaxation with cuts */
  double tighter = 0.0;

  /** Total time spent generating cuts (sum across all generators) */
  double cut_time = 0.0;

  /** Nodes processed during branch-and-cut */
  int nodes = 0;

  /** Iterations processed in the linear programming algorithm */
  int iterations = 0;

  /** number of rows of original problem */
  int nrows = 0;

  /** number of columns of original problem */
  int ncols = 0;

  /** number of rows of preprocessed problem */
  int nprocessedrows = 0;

  /** number of columns of preprocessed problem */
  int nprocessedcols = 0;

  /** Solver status */
  std::string result;

  /** Time (CPU seconds) to build the conflict graph */
  double cgraph_time = 0.0;

  /** Density of the conflict graph */
  double cgraph_density = 0.0;

  /** Clique strengthening: constraints extended (before LP) */
  int clqstr_extended = 0;

  /** Clique strengthening: constraints dominated (before LP) */
  int clqstr_dominated = 0;

  /** Clique strengthening time (before LP) */
  double clqstr_time = 0.0;

  /** Per-cut-generator statistics (one entry per registered generator) */
  std::vector<CutGeneratorStats> cutStats;

  /** Per-heuristic statistics (one entry per registered heuristic) */
  std::vector<HeuristicStats> heuristicStats;

  /**
   * Canonical set of cut-generator names known to MIPster.
   * Always written to the CSV even when the generator produced no cuts,
   * so that the column set stays stable across runs with different settings.
   * Names must match what CbcSolverCutSetup.cpp passes to addCutGenerator().
   */
  static const std::vector<std::string> &knownCutGenerators();

  /**
   * Canonical set of heuristic names known to MIPster.
   * Always written to the CSV even when the heuristic was never called.
   * Names must match setHeuristicName() calls in CbcSolverHeuristics.cpp
   * and related files (case-sensitive; spaces preserved).
   */
  static const std::vector<std::string> &knownHeuristics();

  /**
   * Append the collected statistics to a CSV file.
   *
   * @param outFileName Fully qualified path to the CSV file.
   * @param inputQueue Tokens representing the original command line.
   * @return true on success, false if the file could not be opened.
   */
  bool writeCsv(
                CbcParameters &parameters,
                const std::string &outFileName,
                const std::deque<std::string> &inputQueue) const;
};

#endif
