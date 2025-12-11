#ifndef CBC_SOLVER_STATISTICS
#define CBC_SOLVER_STATISTICS

#include <deque>
#include <string>
#include "CbcParameters.hpp"

class CbcSolverStatistics {
public:
  /** Elapsed total time */
  double seconds = 0.0;

  /** Best solution fount in the search */
  double obj = 0.0;

  /** CPU time */
  double sys_seconds = 0.0;

  /** Elapsed time from solver start */
  double elapsed_seconds = 0.0;

  /** LP relaxation cost */
  double continuous = 0.0;

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

  /** Number of cut generators */
  int number_generators = 0;

  /** Number of cuts per cut generator */
  int *number_cuts = NULL;

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
                const std::deque<std::string> &inputQueue) const;
};

#endif
