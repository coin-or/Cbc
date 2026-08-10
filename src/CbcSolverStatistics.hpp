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
#include "CbcParameters.hpp"

class CBCLIB_EXPORT CbcSolverStatistics {
public:
  /** Elapsed total time */
  double seconds = 0.0;

  /** Best solution found in the search.
   *
   * Defaults to the same "no solution" sentinel CbcModel::bestObjective_
   * uses (see CbcModel.cpp, e.g. `std::min(bestObjective_, 1.0e50)`), NOT
   * 0.0. This field is only overwritten once CbcSolver::run() actually
   * reaches the branch-and-bound completion code (case CbcParam::BAB,
   * `statistics.obj = babModel_->getObjValue();`); if the run stops earlier
   * -- e.g. the root LP relaxation itself is proven infeasible/unbounded or
   * hits the time limit before B&B ever starts (see
   * CbcSolver::solveInitialLp()) -- this field is left untouched. A 0.0
   * default there would be silently misread by -writeStat CSV consumers as
   * a genuine (and often wrong) objective value of zero instead of "no
   * solution found"; 1e50 is unambiguous and already the sentinel every
   * other "no solution" consumer in this codebase checks for.
   */
  double obj = 1.0e50;

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

  /** Time spent generating cuts */
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

  /** Coefficient tightening: coefficients reduced (before LP) */
  int coefstr_changed = 0;

  /** Coefficient tightening: rows with at least one reduced coefficient */
  int coefstr_rows = 0;

  /** Coefficient tightening time (before LP) */
  double coefstr_time = 0.0;

  /** Row reductions: rows removed because all their columns were fixed */
  int rowred_fixed = 0;

  /** Row reductions: rows removed as exact duplicates of another row */
  int rowred_duplicate = 0;

  /** Row reductions: rows removed as scalar multiples of another row */
  int rowred_parallel = 0;

  /** Row reduction time (before LP) */
  double rowred_time = 0.0;

  /** Number of cut generators */
  int number_generators = 0;

  /** Number of cuts per cut generator */
  int *number_cuts = NULL;

  /** Time spent (seconds) per cut generator (0 if timing() not enabled) */
  double *time_generators = NULL;

  /** Cut generator name */
  const char **name_generators = NULL;

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
    const std::deque< std::string > &inputQueue) const;
};

#endif
