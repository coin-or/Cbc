// Copyright (C) 2026, COIN-OR Foundation and others.
// All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

#ifndef CbcHeuristicSppGrasp_H
#define CbcHeuristicSppGrasp_H

#include <string>
#include <vector>
#include "CbcHeuristic.hpp"
#include "CoinHelperFunctions.hpp"

class OsiSolverInterface;

/** GRASP for set partitioning / packing problems.

    Applies to pure binary problems whose rows are all `sum x_j = 1'
    (partitioning) or `sum x_j <= 1' (packing) with unit coefficients.
    Partitioning rows are what make feasibility hard; many real instances
    (e.g. rmp*) have, for every partitioning row, a singleton "slack"
    column that covers it alone at a penalty cost, so a trivial solution
    always exists. The heuristic builds solutions much better than the
    trivial one:

    Construction: repeatedly pick an uncovered partitioning row by a ranking
    criterion (fewest options, regret between its best and second-best
    option, ...), then a column covering it by a static score (cost per
    covered partitioning row, reduced cost, LP value bias, conflict degree),
    with randomised restricted candidate lists for both choices. Columns
    whose cheapest remaining option is a singleton are settled at once.

    Local search: insert a column, eject the columns it conflicts with and
    cover the partitioning rows left open with their slacks, whenever this
    lowers the cost (first improvement until a local optimum).

    LNS (run()): from the deterministic construction (or the incumbent),
    remove the columns conflicting with the best options of a few rows
    (mostly rows left to their slacks), rebuild and improve, accepting
    solutions within acceptRelative of the best; the neighbourhood grows on
    stagnation. Defaults were tuned on the rmp01-18 instances with
    test/spp-grasp-bench.
*/
class CBCLIB_EXPORT CbcSppGrasp {
public:
  enum RowRule {
    RowMinOptions = 0, ///< fewest alive options
    RowMaxRegret, ///< largest gap between best and second-best option
    RowOptionsOverRegret, ///< options / regret
    RowLogOptionsOverRegret, ///< log(1+options) / regret
    RowMinBest, ///< cheapest best option first (pure column greedy)
    RowMaxBest, ///< most expensive best option first
    RowRandom, ///< any uncovered row
    RowEndMarker
  };
  enum ColRule {
    ColCostPerRow = 0, ///< cost / (#partitioning rows)^beta
    ColReducedCost, ///< LP reduced cost, ties by cost per row
    ColProfit, ///< (cost - slack costs saved) / (#partitioning rows)^beta
    ColEndMarker
  };

  struct Config {
    int rowRule;
    int colRule;
    /// cost-per-row exponent (ColCostPerRow)
    double beta;
    /// relative tolerance of the column restricted candidate list
    double colAlpha;
    /// maximum size of the column restricted candidate list
    int colRclMax;
    /// relative tolerance of the row restricted candidate list
    double rowAlpha;
    /// subtract lpWeight * x*_j * (mean |score|) from column scores
    double lpWeight;
    /// add conflictWeight * (degree_j / mean degree) * (mean |score|)
    double conflictWeight;
    /// 0 none, 1 insert/eject with slack repair, 2 with greedy repair
    int localSearch;
    Config()
      : rowRule(RowMaxRegret)
      , colRule(ColProfit)
      , beta(1.0)
      , colAlpha(0.1)
      , colRclMax(COIN_INT_MAX)
      , rowAlpha(0.2)
      , lpWeight(20.0)
      , conflictWeight(0.0)
      , localSearch(2)
    {
    }
  };

  CbcSppGrasp();

  /** Extract the structure from a solver (first numberRows rows only, all
      if negative). Returns true if the problem is a set partitioning /
      packing problem this heuristic can handle; otherwise reason() says
      why not. Uses the solver's current column bounds. */
  bool load(const OsiSolverInterface &solver, int numberRows = -1);

  bool applicable() const { return applicable_; }
  const std::string &reason() const { return reason_; }
  int numberPartitioningRows() const { return numberEqual_; }
  /// partitioning rows covered by a singleton (slack) column
  int numberSlackRows() const { return numberSlackRows_; }
  /// true if every partitioning row has a singleton column
  bool trivialSolutionExists() const { return applicable_ && numberSlackRows_ == numberEqual_; }

  /** LP solution and reduced costs used by lpWeight and ColReducedCost (in
      the loaded column order); NULL clears them. */
  void setLpSolution(const double *x, const double *reducedCost);
  void setSeed(unsigned int seed) { seed_ = seed ? seed : 1; }

  /** Builds one solution (then improves it if config.localSearch);
      returns its cost or COIN_DBL_MAX if construction hit a dead end.
      solution has one 0/1 entry per column. */
  double construct(const Config &config, std::vector< char > &solution,
    const std::vector< char > *partial = NULL);

  struct RunParams {
    /// GRASP constructions (the first is deterministic)
    int graspIterations;
    /// destroy/rebuild iterations from the best GRASP solution
    int lnsIterations;
    double timeLimit;
    /// partitioning rows (left to slacks) around which to destroy
    int destroyRows;
    /// best columns of each such row whose conflicting columns are removed
    int seedColumns;
    /// probability of removing any other column
    double destroyFraction;
    /// LNS also accepts solutions within this relative distance of the best
    double acceptRelative;
    /// after this many LNS iterations without a new best, destroyRows grows
    /// (x2, up to maxDestroyRows); reset on improvement. 0 = never
    int stallIterations;
    int maxDestroyRows;
    /// seed columns drawn at random among the seedPool*seedColumns best
    int seedPool;
    /// probability that a seed row is any partitioning row, not a slack one
    double randomRowFraction;
    RunParams()
      : graspIterations(1)
      , lnsIterations(COIN_INT_MAX)
      , timeLimit(2.0)
      , destroyRows(5)
      , seedColumns(5)
      , destroyFraction(0.0)
      , acceptRelative(0.005)
      , stallIterations(200)
      , maxDestroyRows(40)
      , seedPool(1)
      , randomRowFraction(0.1)
    {
    }
  };
  struct RunStats {
    double trivial, first, time, timeToBest;
    int graspIterations, lnsIterations;
    RunStats()
      : trivial(COIN_DBL_MAX)
      , first(COIN_DBL_MAX)
      , time(0.0)
      , timeToBest(0.0)
      , graspIterations(0)
      , lnsIterations(0)
    {
    }
  };
  /** GRASP followed by LNS; best receives the best solution found (the
      trivial one at worst). Returns its cost, COIN_DBL_MAX if none. */
  double run(const Config &config, const RunParams &params,
    std::vector< char > &best, RunStats *stats = NULL,
    const std::vector< char > *initial = NULL);
  /** insert/eject local search; returns the new cost. level 1 covers the
      partitioning rows left open with slacks, level 2 greedily re-covers
      them with the cheapest (per row) columns that fit, slacks last. */
  /** Insert/eject local search; with reference, starts only from rows not
      covered by a column of reference. */
  double localSearch(std::vector< char > &solution, int level = 2,
    const std::vector< char > *reference = NULL);
  /// cost of the trivial solution (slacks only, plus forced columns)
  double trivialSolution(std::vector< char > &solution);
  /// cost of a 0/1 vector, COIN_DBL_MAX if infeasible
  double evaluate(const std::vector< char > &solution) const;

  int numberColumns() const { return numberColumns_; }
  /// number of slack (singleton) columns at 1 in solution
  int numberSlacksUsed(const std::vector< char > &solution) const;

private:
  void computeScores(const Config &config);
  int randomInt(int n);

  bool applicable_;
  std::string reason_;
  int numberColumns_;
  int numberRows_;
  int numberEqual_;
  int numberSlackRows_;
  double objOffset_;
  // column-major (rows of each column), all rows
  std::vector< int > colStart_, colRows_;
  // row-major (columns of each row), usable columns only
  std::vector< int > rowStart_, rowCols_;
  std::vector< double > cost_;
  std::vector< char > isEqual_; // per row
  std::vector< char > usable_; // per column: upper bound 1
  std::vector< char > forced_; // per column: lower bound 1
  std::vector< int > numberEqualInCol_;
  std::vector< int > degree_;
  std::vector< int > slackCol_; // cheapest singleton column per row, -1 if none
  std::vector< double > lpX_, lpDj_;
  // per partitioning row: columns by cost per partitioning row, cheapest
  std::vector< int > repairOrder_;
  std::vector< double > minCostPerRow_;
  // scores of the last configuration
  std::vector< double > score_;
  std::vector< int > sortedRowCols_; // per row, by score
  int scoreRule_;
  double scoreBeta_, scoreLpWeight_, scoreConflictWeight_;
  bool scoresValid_;
  unsigned int seed_;
  // work arrays
  std::vector< char > alive_;
  std::vector< int > count_, pointer_, owner_, mark_;
  int markStamp_;
};

/** Cbc heuristic wrapper around CbcSppGrasp; does nothing (cheaply) on
    problems that are not set partitioning / packing. */
class CBCLIB_EXPORT CbcHeuristicSppGrasp : public CbcHeuristic {
public:
  CbcHeuristicSppGrasp();
  CbcHeuristicSppGrasp(CbcModel &model);
  CbcHeuristicSppGrasp(const CbcHeuristicSppGrasp &);
  ~CbcHeuristicSppGrasp();
  CbcHeuristicSppGrasp &operator=(const CbcHeuristicSppGrasp &);
  virtual CbcHeuristic *clone() const override;
  virtual void resetModel(CbcModel *model) override;
  virtual void setModel(CbcModel *model) override;
  virtual void generateCpp(FILE *fp) override;

  using CbcHeuristic::solution;
  virtual int solution(double &objectiveValue, double *newSolution) override;

  /// True if the model's solver has the structure the heuristic needs
  static bool isApplicable(const OsiSolverInterface &solver, bool needSlacks,
    std::string *reason = NULL);

  /// GRASP constructions, LNS iterations and time limit per call
  inline void setRunParams(const CbcSppGrasp::RunParams &p) { params_ = p; }
  inline const CbcSppGrasp::RunParams &runParams() const { return params_; }
  inline void setMaxTime(double t) { params_.timeLimit = t; }
  inline double maxTime() const { return params_.timeLimit; }
  inline void setConfig(const CbcSppGrasp::Config &c) { config_ = c; }
  inline const CbcSppGrasp::Config &config() const { return config_; }
  /// 1: only run if every partitioning row has a slack column
  inline void setNeedSlacks(bool v) { needSlacks_ = v; }
  inline bool needSlacks() const { return needSlacks_; }

protected:
  CbcSppGrasp::Config config_;
  CbcSppGrasp::RunParams params_;
  bool needSlacks_;
  int numberCalls_;
};

#endif
