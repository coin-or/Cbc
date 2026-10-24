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
    /// degree exponent for cost/profit scores and reduced-cost tie-breaking
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
    /// greedy repair examines at most this many columns per freed row
    int repairScan;
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
      , repairScan(64)
    {
    }
  };

  CbcSppGrasp();

  /** Extract the structure from a solver (first numberRows rows only, all
      if negative). Returns true if the problem is a set partitioning /
      packing problem this heuristic can handle; otherwise reason() says
      why not. Uses the solver's current column bounds. */
  bool load(const OsiSolverInterface &solver, int numberRows = -1);

  /** Before load(): if factor > 0, every partitioning row without a
      singleton column gets an artificial one (internal only, columns
      numberOriginalColumns() .. numberColumns()-1) costing the row's most
      expensive option plus (factor - 1) times max(|that cost|, largest
      |cost| of any column, 1 if all are 0); run() multiplies that penalty
      by RunParams::artificialGrowth while the best solution has some.
      A solution is feasible for the real problem only if it uses none
      (numberArtificialUsed() == 0); run() returns such a solution if it
      found one. 0 (default) disables them. */
  void setArtificialSlacks(double factor) { artificialFactor_ = factor; }
  /// multiplies the artificial penalties of load() by scale (1 restores)
  void scaleArtificialPenalty(double scale);
  double artificialSlacks() const { return artificialFactor_; }
  int numberArtificial() const { return numberColumns_ - numberOriginalColumns_; }
  int numberOriginalColumns() const { return numberOriginalColumns_; }
  /** rows other than partitioning/packing ones that zero satisfies (e.g.
      x - y <= 0 left by preprocessing); columns that could violate them
      (numberRestrictedColumns) are kept at zero */
  int numberRestrictingRows() const { return numberRestrictingRows_; }
  int numberRestrictedColumns() const { return numberRestrictedColumns_; }
  /** restricting rows keep so many columns at zero (over 10%) that the
      heuristic is not worth running */
  bool tooRestricted() const { return 10 * numberRestrictedColumns_ > numberOriginalColumns_; }
  /** binary columns with only -1 coefficients, handled internally as
      1 - x (see toOriginal); objective values are unchanged */
  int numberComplemented() const { return numberComplemented_; }
  /// original column values (numberOriginalColumns()) of a solution
  void toOriginal(const std::vector< char > &solution, double *x) const;
  /// solution (numberColumns(), artificials at 0) of original values
  void fromOriginal(const double *x, std::vector< char > &solution) const;
  /// number of artificial columns at 1 in solution
  int numberArtificialUsed(const std::vector< char > &solution) const;

  struct FeatureSummary {
    int partitionRows, nativeSingletonRows, nativeSingletonMissing;
    int restrictedRows, restrictedColumns, complementedColumns, forcedColumns;
    int lpAvailable;
    double rowOptionsMin, rowOptionsQ25, rowOptionsMedian, rowOptionsQ75, rowOptionsMean;
    double columnDegreeQ25, columnDegreeMedian, columnDegreeQ75, columnDegreeMean;
    double singletonPremiumQ25, singletonPremiumMedian, singletonPremiumQ75;
    double restrictedRowFraction, restrictedColumnFraction, complementedFraction, forcedFraction;
    double lpFractionalityMean, nativeSingletonLpMass;
    FeatureSummary();
  };
  /// Expected O(nnz + n) descriptive structural/available-LP features; no LP solve.
  void featureSummary(FeatureSummary &features) const;
  /// Real-singleton coverage in retained normalized rows; artificials never count.
  /// Skipped raw rows can still touch these columns, so this is not strict raw degree.
  bool meetsNativeCoverage(double minimum, std::string *description = NULL) const;

  bool applicable() const { return applicable_; }
  const std::string &reason() const { return reason_; }
  int numberPartitioningRows() const { return numberEqual_; }
  /// partitioning rows covered by a singleton (slack) column
  int numberSlackRows() const { return numberSlackRows_; }
  /// true if every partitioning row has a singleton (maybe artificial) column
  bool trivialSolutionExists() const { return applicable_ && numberSlackRows_ == numberEqual_; }

  /** LP solution and reduced costs used by lpWeight and ColReducedCost (in
      the loaded column order, numberOriginalColumns() entries); NULL clears
      them. */
  void setLpSolution(const double *x, const double *reducedCost);
  /// Replace only the LP-value bias vector; keep reduced costs from the
  /// actual LP solve unchanged.
  void setLpGuide(const double *x);
  void setSeed(unsigned int seed) { seed_ = seed ? seed : 1; }

  /** Builds one solution (then improves it if config.localSearch);
      returns its cost or COIN_DBL_MAX if construction hit a dead end.
      solution has one 0/1 entry per column. */
  double construct(const Config &config, std::vector< char > &solution,
    const std::vector< char > *partial = NULL,
    const std::vector< char > *previous = NULL);

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
    /** each LNS iteration doubles the number of destroyed rows with this
        probability, repeatedly (geometric, up to maxDestroyRows): mostly
        small moves, occasionally larger ones. 0 = fixed size */
    double destroyGrowProbability;
    /// seed columns drawn at random among the seedPool*seedColumns best
    int seedPool;
    /// probability that a seed row is any partitioning row, not a slack one
    double randomRowFraction;
    /** probability that a further seed row is drawn among the partitioning
        rows freed by the previous seeds of the same iteration (connected
        neighbourhood) instead of independently. 0 = never */
    double relatedRowProbability;
    /// enable related rows only after finding an artificial-free solution
    bool relatedAfterFeasible;
    /** seed rows by tournament: the best of this many random candidates by
        gain (score of the row's current column minus its best option's).
        0 or 1 = uniform */
    int gainTournament;
    /// probability of supplementing greedy LNS repair with bounded search
    double backtrackProbability;
    /// skip extra search if the greedy repair is already artificial-free
    bool backtrackArtificialOnly;
    /// consider real columns before artificials when building branch options
    bool backtrackRealFirst;
    int backtrackWidth;
    int backtrackDepth;
    int backtrackNodes;
    double backtrackWork;
    /// LNS acceptance: threshold, LAHC, LAHC-better or simulated annealing
    int acceptMode;
    /// LAHC history length
    int lahcLength;
    /// simulated annealing initial temperature as a fraction of objective scale
    double saInitialTemperature;
    /// simulated annealing temperature half-life in deterministic work units
    double saCoolingWork;
    /// reheat to initial temperature after this many consecutive rejections; 0 disables
    int saReheatAfter;
    /// 0 exponential, 1 linear-to-zero with the same half-life
    int saSchedule;
    /// use threshold acceptance for this many trials to calibrate temperature; 0 disables
    int saCalibrationIterations;
    /// target probability for the selected warmup positive-delta quantile
    double saCalibrationTarget;
    /// quantile of normalized positive warmup deltas used for calibration
    double saCalibrationQuantile;
    /// stop the LNS after this many iterations without a new best (0 never)
    int maxNoImprove;
    /// stop the LNS once this much work (a deterministic count of basic
    /// operations, see CbcSppGrasp::work()) was done in run(); 0 never
    double maxWork;
    /// hard aggregate work cap across artificial-penalty rounds; 0 disabled
    double totalWorkLimit;
    /// stop after this much work without a better real-feasible incumbent;
    /// 0 disables the progress gate, and no real incumbent means no gate
    double maxRealStallWork;
    /// LNS repair overrides of Config (negative: keep config's)
    double lnsColAlpha;
    double lnsRowAlpha;
    /// -1 inherit; 0/1/2 fixed; 3: 1 then 2, 4: 0 then 2, 5: 2 then 1,
    /// switching permanently after the first genuinely feasible incumbent
    int lnsLocalSearch;
    int lnsRowRule;
    /// LNS local search starts only from rows changed by destroy/repair
    bool lnsFocus;
    /** with artificial columns: when the LNS stops (maxNoImprove/maxWork)
        without a solution free of them, multiply their penalties by
        artificialGrowth and run the LNS again from the best solution, up to
        artificialRounds more times (each with its own maxWork budget) */
    int artificialRounds;
    double artificialGrowth;
    /** if > 0: escalate the penalties after this much LNS work in a round
        without a solution free of artificials, independently of the stop
        rules; maxWork then limits the total LNS work of all rounds */
    double artificialRoundWork;
    RunParams()
      : graspIterations(1)
      , lnsIterations(COIN_INT_MAX)
      , timeLimit(60.0)
      , destroyRows(5)
      , seedColumns(2)
      , destroyFraction(0.0)
      , acceptRelative(0.005)
      , stallIterations(200)
      , maxDestroyRows(40)
      , destroyGrowProbability(0.0)
      , seedPool(1)
      , randomRowFraction(0.2)
      , relatedRowProbability(0.0)
      , relatedAfterFeasible(false)
      , gainTournament(0)
      , backtrackProbability(0.0)
      , backtrackArtificialOnly(false)
      , backtrackRealFirst(false)
      , backtrackWidth(3)
      , backtrackDepth(8)
      , backtrackNodes(200)
      , backtrackWork(1.0e6)
      , acceptMode(AcceptThreshold)
      , lahcLength(50)
      , saInitialTemperature(0.05)
      , saCoolingWork(1.0e9)
      , saReheatAfter(0)
      , saSchedule(0)
      , saCalibrationIterations(0)
      , saCalibrationTarget(0.5)
      , saCalibrationQuantile(0.5)
      , maxNoImprove(1000)
      , maxWork(1.5e9)
      , totalWorkLimit(0.0)
      , maxRealStallWork(0.0)
      , lnsColAlpha(-1.0)
      , lnsRowAlpha(-1.0)
      , lnsLocalSearch(-1)
      , lnsRowRule(RowRandom)
      , lnsFocus(true)
      , artificialRounds(2)
      , artificialGrowth(10.0)
      , artificialRoundWork(0.0)
    {
    }
  };
  enum AcceptMode {
    AcceptThreshold = 0,
    AcceptLahc = 1,
    AcceptLahcBetter = 2,
    AcceptSimulatedAnnealing = 3
  };
  /// a new best found by run(): LNS iteration (0 = before LNS), cost, time
  struct Improvement {
    int iteration;
    double value;
    double time;
    double work;
  };
  struct ProgressEvent {
    int iteration, kind, penaltyRound, artificialUsed;
    double value, time, work;
  };
  struct RepairStats {
    long long calls, skippedFeasible, nodes, leaves, realLeaves, deadEnds;
    long long workExits, nodeExits, exhausted, invalidPartial;
    long long improved, feasibilityGains, greedyArtificial, noLeaf;
    long long openRows, maxOpenRows, maxDepth;
    double work;
    RepairStats()
      : calls(0)
      , skippedFeasible(0)
      , nodes(0)
      , leaves(0)
      , realLeaves(0)
      , deadEnds(0)
      , workExits(0)
      , nodeExits(0)
      , exhausted(0)
      , invalidPartial(0)
      , improved(0)
      , feasibilityGains(0)
      , greedyArtificial(0)
      , noLeaf(0)
      , openRows(0)
      , maxOpenRows(0)
      , maxDepth(0)
      , work(0.0)
    {
    }
  };
  struct AcceptanceStats {
    long long trials, accepted, acceptedWorse, rejected, invalid, reheats;
    int calibrationTrials, calibrationPositive;
    double calibratedTemperature;
    AcceptanceStats()
      : trials(0)
      , accepted(0)
      , acceptedWorse(0)
      , rejected(0)
      , invalid(0)
      , reheats(0)
      , calibrationTrials(0)
      , calibrationPositive(0)
      , calibratedTemperature(0.0)
    {
    }
  };
  struct RunStats {
    double trivial, first, time, timeToBest, work;
    int graspIterations, lnsIterations;
    /** artificial columns in the best internal solution (-1 none found);
        trivial and first count artificial costs, improvements do not */
    int artificialUsed;
    std::vector< Improvement > improvements;
    /** Optional diagnostic events: kind 1 penalized best, 2 real best,
        3 artificial-penalty epoch reset, 4 real-best-stall stop. */
    std::vector< ProgressEvent > progress;
    RepairStats repair;
    AcceptanceStats acceptance;
    bool realStallStop;
    RunStats()
      : trivial(COIN_DBL_MAX)
      , first(COIN_DBL_MAX)
      , time(0.0)
      , timeToBest(0.0)
      , work(0.0)
      , graspIterations(0)
      , lnsIterations(0)
      , artificialUsed(0)
      , realStallStop(false)
    {
    }
  };
  /** GRASP followed by LNS; best receives the best solution found (the
      trivial one at worst). Returns its cost, COIN_DBL_MAX if none. With
      artificial columns, best is the best solution without them if one was
      found, otherwise the best one with them (check numberArtificialUsed). */
  double run(const Config &config, const RunParams &params,
    std::vector< char > &best, RunStats *stats = NULL,
    const std::vector< char > *initial = NULL);
  /** insert/eject local search; returns the new cost. level 1 covers the
      partitioning rows left open with slacks, level 2 greedily re-covers
      them with the cheapest (per row) columns that fit, slacks last. */
  /** Insert/eject local search; with reference, starts only from rows not
      covered by a column of reference; with previous (a local optimum that
      reference was cut from), only from rows whose covering column differs
      from previous's. */
  double localSearch(std::vector< char > &solution, int level = 2,
    const std::vector< char > *reference = NULL, int repairScan = 64,
    const std::vector< char > *previous = NULL);
  /// cost of the trivial solution (slacks only, plus forced columns)
  double trivialSolution(std::vector< char > &solution);
  /// cost of a 0/1 vector, COIN_DBL_MAX if infeasible
  /// accountWork=false is for validation outside the search budget
  double evaluate(const std::vector< char > &solution,
    bool accountWork = true) const;

  /// internal columns (solution vectors' size): original plus artificial
  int numberColumns() const { return numberColumns_; }
  /** Basic operations (mostly matrix entries visited) done so far; a
      deterministic, machine-independent proxy for running time. */
  double work() const { return work_; }
  /// number of slack (singleton) columns at 1 in solution
  int numberSlacksUsed(const std::vector< char > &solution) const;

private:
  double backtrackingRepair(const Config &config, const RunParams &params,
    std::vector< char > &solution, const std::vector< char > &partial,
    const std::vector< char > *previous, double bestValue, RepairStats &stats);
  /// repairOrder_ and minCostPerRow_ from cost_
  void sortRepairOrder();
  void computeScores(const Config &config);
  int randomInt(int n);

  bool applicable_;
  std::string reason_;
  int numberColumns_;
  int numberOriginalColumns_;
  double artificialFactor_;
  int numberRestrictingRows_;
  int numberRestrictedColumns_;
  int numberComplemented_;
  /// artificial k costs base + scale * (factor - 1) * unit
  std::vector< double > artificialBase_;
  std::vector< double > artificialUnit_;
  /// original column j is handled as 1 - x_j
  std::vector< char > complemented_;
  int numberRows_;
  int numberEqual_;
  int numberSlackRows_;
  double objOffset_;
  /// tolerance of local search move deltas (see sortRepairOrder)
  double costTolerance_;
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
  mutable double work_;
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
  /** 0: not applicable (reason says why), 1: every partitioning row has a
      slack column, 2: applicable with artificial slacks (artificialFactor
      > 0) on the partitioning rows without one */
  static int applicability(const OsiSolverInterface &solver, double artificialFactor,
    std::string *reason = NULL, double minimumNativeCoverage = 0.0);

  /// GRASP constructions, LNS iterations and time limit per call
  inline void setRunParams(const CbcSppGrasp::RunParams &p) { params_ = p; }
  inline const CbcSppGrasp::RunParams &runParams() const { return params_; }
  void setPumpFractionalGuide(const double *values, const double *rootLp,
    const double *rootReducedCosts, int n,
    double sourceFractionality, int sourceFractionalCount,
    double rootFractionality, int rootFractionalCount);
  inline void setMaxTime(double t) { params_.timeLimit = t; }
  inline double maxTime() const { return params_.timeLimit; }
  /// deterministic work limit per call (CbcSppGrasp::work() units), 0 none
  inline void setMaxWork(double w) { params_.maxWork = w; }
  inline double maxWork() const { return params_.maxWork; }
  inline void setTotalWorkLimit(double w) { params_.totalWorkLimit = w; }
  /// LNS iterations without a new best before a call stops, 0 none
  inline void setMaxNoImprove(int n) { params_.maxNoImprove = n; }
  inline int maxNoImprove() const { return params_.maxNoImprove; }
  inline void setConfig(const CbcSppGrasp::Config &c) { config_ = c; }
  inline const CbcSppGrasp::Config &config() const { return config_; }
  /// 1: only run if every partitioning row has a slack column
  inline void setNeedSlacks(bool v) { needSlacks_ = v; }
  inline bool needSlacks() const { return needSlacks_; }
  /** artificial slack factor (CbcSppGrasp::setArtificialSlacks), 0 none;
      lets the heuristic run when some partitioning rows lack a slack */
  inline void setArtificialSlacks(double f) { artificialFactor_ = f; }
  inline double artificialSlacks() const { return artificialFactor_; }
  inline void setMinimumNativeCoverage(double fraction) { minimumNativeCoverage_ = fraction; }
  /// Experimental parameter file applied only to forced-free active models
  /// with >= 50% real private coverage. Empty path preserves legacy behavior.
  bool setCoverageTuningProfile(const std::string &path,
    const std::vector< std::string > &explicitOverrides, std::string *reason = NULL,
    const std::string &stage = "unspecified");

protected:
  CbcSppGrasp::Config config_;
  CbcSppGrasp::RunParams params_;
  std::vector< double > pumpFractionalGuide_;
  std::vector< double > rootLpSnapshot_;
  std::vector< double > rootReducedCostsSnapshot_;
  double pumpSourceFractionality_;
  int pumpSourceFractionalCount_;
  double rootLpFractionality_;
  int rootLpFractionalCount_;
  bool needSlacks_;
  double artificialFactor_;
  double minimumNativeCoverage_;
  std::string coverageTuningProfile_;
  std::string coverageTuningStage_;
  std::vector< std::pair< std::string, double > > coverageTuningValues_;
  int numberCalls_;
};

#endif
