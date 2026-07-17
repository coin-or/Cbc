// Copyright (C) 2007, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

/*! \file CbcSolver.hpp
    \brief Defines CbcSolver, the top-level class for the cbc solver driver.

    CbcSolver wraps the full solve pipeline: parameter initialization,
    command parsing, LP solve, preprocessing, cut/heuristic setup,
    branch-and-bound, post-processing, and result reporting.

    The free functions CbcMain0 and CbcMain1 are backward-compatible
    wrappers that delegate to CbcSolver::initialize() and CbcSolver::run().
*/

#ifndef CbcSolver_H
#define CbcSolver_H

#include "CoinUtilsConfig.h"

#include <deque>
#include <string>
#include <vector>

#include "CoinMessageHandler.hpp"
#include "CoinModel.hpp"

#include "OsiClpSolverInterface.hpp"
#if CBC_OTHER_SOLVER == 1
#include "OsiCpxSolverInterface.hpp"
#endif

#include "ClpParameters.hpp"

#include "CglCutGenerator.hpp"
#include "CglPreProcess.hpp"

#include "CbcModel.hpp"
#include "CbcParameters.hpp"
#include "CbcMessage.hpp"
#include "CbcSolverStatistics.hpp"
#include "CbcSolverCutSetup.hpp"

#ifdef COINUTILS_HAS_GLPK
#include "glpk.h"
#endif

class CbcUser;
class CbcStopNow;
class CglStored;
class CbcPreprocHandler;
class CglTwomir;

//#############################################################################
//#############################################################################

/*! \brief Top-level driver class for the CBC MIP solver.

    Encapsulates the full solve pipeline: parameter initialization,
    command parsing, LP solve, preprocessing, cut/heuristic setup,
    branch-and-bound, post-processing, and result reporting.

    Usage patterns:

    \code
    // Pattern 1: One-call solve
    CbcSolver cbc;
    cbc.solve("problem.mps");
    if (cbc.hasSolution())
      printf("Obj: %g\n", cbc.objectiveValue());

    // Pattern 2: Import, configure, then solve
    CbcSolver cbc;
    cbc.importModel("problem.mps");
    cbc.parameters()[CbcParam::TIMELIMIT]->setVal(300.0);
    cbc.parameters()[CbcParam::GOMORYCUTS]->setVal("off");
    std::deque<std::string> q = {"-solve", "-quit"};
    cbc.run(q);

    // Pattern 3: argc/argv (backward compatible)
    CbcSolver cbc(solver);
    cbc.initialize();
    cbc.run(argc, argv);
    \endcode

    The free functions CbcMain0/CbcMain1 still work for backward
    compatibility and delegate to this class internally.
*/

class CBCLIB_EXPORT CbcSolver {

public:
  ///@name Constructors and destructors
  //@{
  /// Default Constructor
  CbcSolver();

  /// Constructor from solver
  explicit CbcSolver(const OsiClpSolverInterface &);

  /// Constructor from model
  explicit CbcSolver(const CbcModel &);

  /** Copy constructor. */
  CbcSolver(const CbcSolver &rhs);

  /// Assignment operator
  CbcSolver &operator=(const CbcSolver &rhs);

  /// Destructor
  ~CbcSolver();
  //@}

  ///@name Initialization (replaces CbcMain0)
  //@{
  /** Initialize default parameters, signal handlers, and LP solver state.
      Must be called before run(). Replaces the body of CbcMain0. */
  void initialize();
  //@}

  ///@name High-level solve
  //@{
  /** Import a model file and solve it with current parameters.
      Convenience method that chains initialize() + import + run(-solve -quit).
      \param filename  Path to MPS/LP/GMPL file (supports .gz/.bz2)
      \return 0 on success
  */
  int solve(const std::string &filename);

  /** Import a model from a file.
      Calls initialize() if not already done, then reads the file.
      \param filename  Path to MPS/LP/GMPL file (supports .gz/.bz2)
      \return 0 on success
  */
  int importModel(const std::string &filename);

  /** Solve the LP relaxation only (no branch-and-bound).
      Uses the current -lpMethod setting (dual/primal/barrier).
      \param method  "dual", "primal", or "barrier" (empty = use current setting)
      \return 0 on success
  */
  int solveLp(const std::string &method = "");
  //@}

  ///@name Result accessors (valid after solve/run)
  //@{
  /// Solve status: 0=optimal, 1=infeasible, 2=unbounded, 3+=stopped
  int status() const;
  /// Best objective value found (COIN_DBL_MAX if no solution)
  double objectiveValue() const;
  /// Best bound (dual bound from B&B tree)
  double bestBound() const;
  /// Best integer solution (NULL if none found)
  const double *bestSolution() const;
  /// Number of columns in the model
  int numCols() const;
  /// Number of B&B nodes explored
  int nodeCount() const;
  /// Number of simplex iterations
  int iterationCount() const;
  /// Whether a feasible solution was found
  bool hasSolution() const;
  //@}

  ///@name Run (replaces CbcMain1)
  //@{
  /** Run the solver with a command queue.
      This is the primary entry point that replaces CbcMain1.
      \param inputQueue  Command tokens (consumed during processing)
      \param callBack    User callback invoked at 6 defined points (default: no-op)
      \param info        AMPL interface data (NULL when not using AMPL)
      \return 0 on normal completion, non-zero if callback requests early exit
  */
  int run(std::deque<std::string> inputQueue,
    int callBack(CbcModel *currentSolver, int whereFrom) = nullptr,
    ampl_info *info = nullptr);

  /** Run the solver with argc/argv arguments.
      Converts to a command queue and delegates to the deque overload. */
  int run(int argc, const char *argv[],
    int callBack(CbcModel *currentSolver, int whereFrom) = nullptr,
    ampl_info *info = nullptr);
  //@}

  ///@name Heuristics and cut generator configuration
  //@{
  /** Configure/run heuristics on the given model, using this solver's
      current parameter and tuning settings (noPrinting, initialPumpTune_).
      Thin wrapper exposing the free function doHeuristics() as a
      class-level call surface, so heuristic setup can be invoked directly
      on a CbcSolver instance rather than only from run().
      \param model  the CbcModel to configure/run heuristics on
      \param type   1=register heuristics on model for the B&B tree,
                    2=run heuristics now (may set cutoff/best solution),
                    3=miplib variant (skips some heuristics),
                    10+=use lastSolution and relax a few variables first
      \return number of heuristic solutions found (type 2), or 0 for
              registration-only calls
  */
  int configureHeuristics(CbcModel *model, int type);

  /** Register all cut generators on babModel based on this solver's
      current parameter, conflict-graph, and clique-detection settings.
      Thin wrapper exposing the free function installCutGenerators() as a
      class-level call surface.
      \param babModel  the CbcModel to install cut generators on
      \param miplib    true when running the miplib test variant (some
                       cut generators are tuned differently)
      \param bkPivotingStrategy  Bron-Kerbosch pivoting strategy for
                       clique/odd-wheel cut generation
  */
  void configureCutGenerators(CbcModel &babModel, bool miplib,
    CoinBronKerbosch::PivotingStrategy bkPivotingStrategy);
  //@}

  ///@name User extensions
  //@{
  /// Add user function
  void addUserFunction(CbcUser *function);
  /// Set user call back
  void setUserCallBack(CbcStopNow *function);
  /// Add cut generator
  void addCutGenerator(CglCutGenerator *generator);
  //@}

  ///@name Accessors
  //@{
  /// Get int parameter value
  int intValue(int code);
  /// Set int parameter value
  void setIntValue(int code, int value);
  /// Get double parameter value
  double doubleValue(int code);
  /// Set double parameter value
  void setDoubleValue(int code, double value);
  /// User function (NULL if no match)
  CbcUser *userFunction(const char *name) const;
  /// Return reference CbcModel
  inline CbcModel *model() { return &model_; }
  /// Return updated CbcModel (after B&B)
  inline CbcModel *babModel() { return babModel_; }
  /// Return parameters
  inline CbcParameters &parameters() { return parameters_; }
  inline const CbcParameters &parameters() const { return parameters_; }
  /// Number of user functions
  inline int numberUserFunctions() const { return numberUserFunctions_; }
  /// User function array
  inline CbcUser **userFunctionArray() const { return userFunction_; }
  /// Original solver (will contain output solutions)
  inline OsiClpSolverInterface *originalSolver() const { return originalSolver_; }
  /// Original CoinModel
  inline CoinModel *originalCoinModel() const { return originalCoinModel_; }
  /// Set original solver
  void setOriginalSolver(OsiClpSolverInterface *originalSolver);
  /// Set original CoinModel
  void setOriginalCoinModel(CoinModel *originalCoinModel);
  /// Number of cut generators
  inline int numberCutGenerators() const { return numberCutGenerators_; }
  /// Cut generator array
  inline CglCutGenerator **cutGeneratorArray() const { return cutGenerator_; }
  /// Start time
  inline double startTime() const { return startTime_; }
  /// Whether to print to std::cout
  inline void setPrinting(bool onOff) { noPrinting_ = !onOff; }
  /// Where to start reading commands
  inline void setReadMode(int value) { readMode_ = value; }
  /// Solve statistics (populated after run())
  inline const CbcSolverStatistics &getStatistics() const { return statistics_; }
  //@}

private:
  ///@name Core state (existed in original CbcSolver class)
  //@{
  /// Reference model
  CbcModel model_;
  /// Updated model (created during B&B, deleted after)
  CbcModel *babModel_;
  /// User functions
  CbcUser **userFunction_;
  /// Status of user functions (0=not used, 1=needs load, 2=available, 3=loaded)
  int *statusUserFunction_;
  /// Original solver (will contain output solutions)
  OsiClpSolverInterface *originalSolver_;
  /// Original CoinModel
  CoinModel *originalCoinModel_;
  /// Cut generators
  CglCutGenerator **cutGenerator_;
  /// Number of user functions
  int numberUserFunctions_;
  /// Number of cut generators
  int numberCutGenerators_;
  /// Stop-now callback
  CbcStopNow *callBack_;
  /// CPU time at instantiation
  double startTime_;
  /// Parameters
  ClpParameters clpParameters_;
  CbcParameters parameters_;
  /// Whether to do miplib test
  bool doMiplib_;
  /// Whether to suppress printing
  bool noPrinting_;
  /// Where to start reading commands
  int readMode_;
  //@}

  ///@name Cross-phase state from CbcMain1 (previously local variables)
  //@{

  // --- Preprocessing state ---
  /// Preprocessing engine
  CglPreProcess process_;
  /// Pre-preprocessing solver clone (saved for postprocessing)
  OsiSolverInterface *saveSolver_;

  // --- Solve control flags ---
  /// Whether a valid model is loaded
  bool goodModel_;
  /// Whether in interactive mode
  bool interactiveMode_;
  /// False if user changed any cut/heuristic settings
  bool defaultSettings_;
  /// Presolve level (0=off, 5=default)
  int preSolve_;
  /// Preprocessing level (0=off, 4=default)
  int preProcess_;
  /// Whether to use strategy
  bool useStrategy_;
  /// Whether presolve writes to file
  bool preSolveFile_;
  /// Whether strong branching was changed by user
  bool strongChanged_;
  /// Whether FPump was changed by user
  bool pumpChanged_;
  /// Max cut passes at root (-1234567 = auto)
  int cutPass_;
  /// Max cut passes in tree (-1234567 = auto)
  int cutPassInTree_;
  /// Preprocessing tuning flags
  int tunePreProcess_;
  /// Test OSI parameters flag
  int testOsiParameters_;
  /// 0 normal, 1 from AMPL or MIQP etc (2 allows cuts)
  int complicatedInteger_;
  /// Feasibility pump tuning parameter
  int initialPumpTune_;
  /// Reduced cost fixing threshold
  double djFix_;
  /// Tighten factor
  double tightenFactor_;
  /// Normal cutoff increment
  double normalIncrement_;
  /// Return mode (0=untouched, 1=updated, 2=as babModel)
  int returnMode_;
  /// Solve result (0=optimal, 3=stopped, 6=infeasible, -1=not solved)
  int integerStatus_;
  /// Total number of valid commands processed
  int numberGoodCommands_;
  /// Node strategy
  int nodeStrategy_;
  /// Whether dominated cuts are used
  bool dominatedCuts_;
  /// SOS handling flag
  int doSOS_;
  /// Verbose level
  int verbose_;
  /// Cost-based priorities mode
  int useCosts_;
  /// Whether to use input solution (-1 = no)
  int useSolution_;
  /// Index of current best solution
  int currentBestSolution_;

  // --- LP solver control ---
  /// Idiot crash method (-1 = auto)
  int doIdiot_;
  /// Output format (2 = default)
  int outputFormat_;
  /// SLP value
  int slpValue_;
  /// C++ code generation value
  int cppValue_;
  /// Print options
  int printOptions_;
  /// Print mode
  int printMode_;
  /// Presolve options
  int presolveOptions_;
  /// Substitution level
  int substitution_;
  /// Dualize level
  int dualize_;
  /// Crash method
  int doCrash_;
  /// Vector mode
  int doVector_;
  /// Sprint method (-1 = auto)
  int doSprint_;
  /// Scaling method
  int doScaling_;

  // --- Barrier solver control ---
  /// Cholesky type
  int choleskyType_;
  /// Gamma for barrier
  int gamma_;
  /// Barrier scaling
  int scaleBarrier_;
  /// KKT method
  int doKKT_;
  /// Crossover method (2 = do unless quadratic)
  int crossover_;
  /// Whether problem is bilinear
  bool biLinearProblem_;

  // --- Cut generator modes ---
  int gomoryMode_;
  int probingMode_;
  int knapsackMode_;
  int redsplitMode_;
  int redsplit2Mode_;
  int GMIMode_;
  int cliqueMode_;
  int oldCliqueMode_;
  int oddWheelMode_;
  int mixedMode_;
  int mixedRoundStrategy_;
  int flowMode_;
  int twomirMode_;
  int landpMode_;
  int residualCapacityMode_;
  int zerohalfMode_;
  /// Conflict graph mode ("on"/"off"/"clq")
  std::string cgraphMode_;
  /// Clique strengthening mode ("before"/"after"/"off")
  std::string clqstrMode_;
  /// BK clique pivoting strategy
  int bkPivotingStrategy_;
  /// Max calls to BK
  int maxCallsBK_;
  /// BK clique extension method
  int bkClqExtMethod_;
  /// Odd wheel extension method
  int oddWExtMethod_;

  // --- Branching input (from AMPL or priority files) ---
  int *priorities_;
  int *branchDirection_;
  double *pseudoDown_;
  double *pseudoUp_;
  double *solutionIn_;
  int *prioritiesIn_;
  int numberSOS_;
  int *sosStart_;
  int *sosIndices_;
  char *sosType_;
  double *sosReference_;
  int *cut_;
  int *sosPriority_;

  // --- MIP start ---
  std::vector<std::pair<std::string, double>> mipStart_;
  std::vector<std::pair<std::string, double>> mipStartBefore_;
  std::string mipStartFile_;

  // --- Knapsack expansion ---
  int *whichColumn_;
  int *knapsackStart_;
  int *knapsackRow_;
  int numberKnapsack_;

  // --- Import control ---
  int allowImportErrors_;
  int keepImportNames_;

  // --- Names ---
  int lengthName_;
  std::vector<std::string> rowNames_;
  std::vector<std::string> columnNames_;

  // --- Debug ---
  double *debugValues_;
  int numberDebugValues_;
  int basisHasValues_;

  // --- Lot sizing ---
  struct LotStruct {
    double low;
    double high;
    int column;
  };
  LotStruct *lotsize_;
  int numberLotSizing_;

  // --- GLPK state ---
#ifdef COINUTILS_HAS_GLPK
  glp_tran *coin_glp_tran_;
  glp_prob *coin_glp_prob_;
#endif

  // --- Timing ---
  double totalTime_;
  double time0_;
  double time0Elapsed_;

  // --- Statistics ---
  CbcSolverStatistics statistics_;

  // --- Stored AMPL cuts ---
  // (managed internally during run, not exposed)

  // --- Input queue copies ---
  std::deque<std::string> saveInputQueue_;

  /// Collected parameter-change messages (printed as a section before solve)
  std::vector<std::string> paramChanges_;
  //@}

  ///@name Private helpers
  //@{
  /// Reset cross-phase state to defaults (called by initialize)
  void resetRunState();

  /// Print collected parameter changes as a formatted section, then clear.
  void printParamChanges();

  /** Handle the IMPORT action: read an MPS/LP/GMPL file into the model.
      Called from run() when the IMPORT command is encountered.
      \param inputQueue  Command queue (for reading filename and .par file)
      \param clpSolver   Current OsiClp solver interface
      \param lpSolver    Current ClpSimplex (may be updated)
      \param time1       CPU time reference (updated on success)
      \param totalTime   Accumulated time (updated on success)
  */
  void importModel(std::deque<std::string> &inputQueue,
    OsiClpSolverInterface *clpSolver, ClpSimplex *lpSolver,
    double &time1, double &totalTime);

  /** Run CglPreProcess on the model before branch-and-bound.
      Extracted from run() — the `if (preProcess && cbcParamCode == CbcParam::BAB)` block.
      \return 0=ok, 1=break, 2=continue, 3=return (returnCode set)
  */
  int preprocess(int preProcess, int cbcParamCode,
    OsiClpSolverInterface *&clpSolver, ClpSimplex *&lpSolver,
    CglPreProcess &process, CbcSolverStatistics &statistics,
    int &returnCode, ampl_info *info);

  /** Handle the WRITESOL/PRINTSOL/WRITEGMPLSOL/WRITENEXTSOL action.
      Called from run() when a solution-writing command is encountered.
      \param cbcParamCode  Which solution-write variant was requested
      \param inputQueue    Command queue (for reading filename)
      \param clpSolver     Current OsiClp solver interface
      \param lpSolver      Current ClpSimplex
  */
  void writeSolution(int cbcParamCode,
    std::deque<std::string> &inputQueue,
    OsiClpSolverInterface *clpSolver, ClpSimplex *lpSolver);

  /** Apply the LP method to solve the current LP relaxation.
      Unified LP-solve entry point called from both the BAB path
      (via solveInitialLp()) and from the SOLVECONTINUOUS (-initialSolve)
      action.  The function runs, in order:
        1. Fast MILP preprocessing (bound tightening, if enabled)
        2. Clique merging "before" (if clqstrMode_ == "before")
        3. Model-level LP settings: PSI positive-edge, objective scaling,
           vector mode — applied before racing so thread clones inherit them
        4. LP racing (if enabled; threads use intentionally varied configs)
        5. Full ClpSolve LP solve with all user-settable options:
           presolve, crash/idiot/sprint, SLP, barrier options, dualize,
           print options, time limits.
      The caller is responsible for setting up any LP-progress message
      handler on the ClpSimplex model before calling this function.
      \return -1 if infeasibility proved during preprocessing (skip solve),
               0 on success (LP solve performed or racing winner found).
  */
  int applyLpMethod();

  /** Solve the root LP relaxation.
      Called from the BAB action when !miplib.
      \return 0=success, 1=break BAB, 2=continue loop, 3=return from run()
  */
  int solveInitialLp(
    int logLevel, int cbcLogLevel,
    CbcSolverStatistics &statistics,
    int &returnCode,
    int callBack(CbcModel *currentSolver, int whereFrom));

  /** Run CglPreProcess preprocessing on the model.
      Called from the BAB action when preProcess is enabled.
      \return 0=success, 1=break BAB, 3=return from run()
  */
  int preprocess(
    CglPreProcess &process,
    OsiClpSolverInterface *&clpSolver, ClpSimplex *&lpSolver,
    OsiClpSolverInterface *originalSolver,
    CbcSolverStatistics &statistics,
    int &returnCode,
    int callBack(CbcModel *currentSolver, int whereFrom),
    ampl_info *info,
    int &truncateColumns, int &truncateRows,
    bool &redoSOS, double *&truncatedRhsLower, double *&truncatedRhsUpper,
    int *&newPriorities, bool integersOK, int numberOriginalColumns,
    CbcPreprocHandler *&preprocHandler, double &preprocStart);

  /** Run postprocessing after branch-and-bound and report results.
      Called at the end of the BAB action.
      \return 0=success (caller breaks BAB), 3=return from run()
  */
  int postprocess(
    CglPreProcess &process,
    OsiClpSolverInterface *clpSolver, ClpSimplex *lpSolver,
    OsiClpSolverInterface *originalSolver,
    CbcSolverStatistics &statistics,
    double time1, double time1Elapsed,
    int &returnCode,
    int callBack(CbcModel *currentSolver, int whereFrom),
    ampl_info *info);

  /** Handle the BOUNDPROP action: run standalone bound propagation
      (singletons / fixpoint / MILPbt) on the current model.
      Extracted from run() — the `case CbcParam::BOUNDPROP:` block.
      \return 0=success (caller falls through to break), 1=continue the
              command loop (model invalid, feature disabled, or infeasible)
  */
  int runBoundPropagation();

  /** Handle the CHECKSOLUTION action: write a solution-validation report
      (row/column diagnostics) for the current model to a file.
      Extracted from run() — the `case CbcParam::CHECKSOLUTION:` block.
      \return 0=success (caller falls through to break), 1=continue the
              command loop (model invalid or file could not be written)
  */
  int runCheckSolution(CbcParam *cbcParam, std::deque< std::string > &inputQueue);

  /** Handle the setup + root-LP-relaxation portion of the BAB action:
      cutoff sign flip, printout hints, the legacy OsiSolverLink/quadratic
      model path (COIN_HAS_LINK), and the solveInitialLp() call.
      Extracted from run() — the first part of `case CbcParam::BAB:`.
      \return 0=success (caller continues), 1=break BAB, 2=continue the
              command loop, 3=return from run() (returnCode is set)
  */
  int babSetupAndRootLp(bool miplib, int logLevel, int cbcLogLevel,
    OsiClpSolverInterface *&clpSolver, ClpSimplex *&lpSolver,
    CglStored &storedAmpl, CoinModel *&coinModel,
    CbcSolverStatistics &statistics, int &returnCode,
    int callBack(CbcModel *currentSolver, int whereFrom));

  /** Configuration for the conflict-graph branching ranker, collected from
      parameters in run() before `case CbcParam::BAB:` is reached and
      consumed by babConfigureBabModel(). Field names mirror the run()-scope
      locals they were extracted from.
  */
  struct RankerConfig {
    double rankConflictWeight;
    std::string rankConflictType;
    double rankConflictPowerTrusted;
    double rankConflictPowerUntrusted;
    double rankRangeWeight;
    double rankRangePowerTrusted;
    double rankRangePowerUntrusted;
    double rankRangeMax;
    double rankNzWeight;
    double rankNzPowerTrusted;
    double rankNzPowerUntrusted;
    double rankObjCoeffWeight;
    double rankObjCoeffPowerTrusted;
    double rankObjCoeffPowerUntrusted;
    double rankConflictMaxPercBin;
  };

  /** Handle the babModel_ (re)configuration portion of the BAB action:
      applies user settings (scaling hint, tighten-bounds-on-continuous
      infeasibility check), builds the conflict-graph branching ranker (if
      requested), (re)constructs babModel_ from model_, tunes the
      factorization frequency, adds lotsizing objects, and decides whether
      preprocessing should be switched off based on the dual objective
      limit. Extracted from run() — the second part of `case CbcParam::BAB:`,
      immediately following babSetupAndRootLp().
      \return 0=success (caller continues), 1=break BAB (tighten bounds
              found the problem infeasible)
  */
  int babConfigureBabModel(bool miplib, int logLevel,
    OsiClpSolverInterface *&clpSolver, ClpSimplex *&lpSolver,
    CglStored &storedAmpl,
    double time1, double &time2, double &totalTime,
    const RankerConfig &rankerConfig, int &numberOriginalColumns);

  /** Handle the cgraph/clique-strengthening, preprocessing phase-end
      banner cleanup, and post-preprocessing bound tightening portion of
      the BAB action. Extracted from run() — the block immediately
      following the preprocess() call in `case CbcParam::BAB:`.
      \return 0=success (caller continues), 1=break BAB (tighten bounds
              found the problem infeasible)
  */
  int babPostPreprocessCleanup(bool miplib,
    CglPreProcess &process, CglTwomir &twomirGen,
    CbcPreprocHandler *&preprocHandler, double preprocStart,
    CbcSolverStatistics &statistics);
  /** Handle the common post-preprocess BAB/MIPLIB search setup:
      linked-model knapsack expansion, cost-based priorities, heuristic/cut
      setup, repeated-use solver tuning, debug-value mapping, and the shared
      branch-and-bound option configuration that runs before the BAB-vs-MIPLIB
      split. Extracted from run() — the first remaining inline block in
      `case CbcParam::BAB:` following babPostPreprocessCleanup().
  */
  void babConfigureSearchModel(bool miplib, int cbcParamCode,
    OsiClpSolverInterface *&clpSolver, ClpSimplex *&lpSolver,
    CglPreProcess &process, CglStored &storedAmpl, CoinModel &saveCoinModel,
    CoinModel &saveTightenedModel,
    CoinBronKerbosch::PivotingStrategy bkPivotingStrategy,
    int *&newPriorities, int &testOsiOptions, bool &integersOK,
    ampl_info *info);
  /** Execute the remaining BAB/MIPLIB search phase after
      babConfigureSearchModel(): BAB-specific object/SOS setup, the actual
      branch-and-bound or MIPLIB execution path, search statistics, and the
      final postprocess() call. Extracted from run() — the remaining inline
      tail of `case CbcParam::BAB:` after babConfigureSearchModel().
      eturn 0=success (caller breaks BAB), 3=return from run()
              (returnCode is set)
  */
  int babExecuteSearchAndPostprocess(int cbcParamCode, int cbcLogLevel,
    int testOsiOptions, OsiClpSolverInterface *&clpSolver,
    ClpSimplex *&lpSolver, OsiClpSolverInterface *originalSolver,
    CglPreProcess &process, CglStored &storedAmpl,
    CbcSolverStatistics &statistics,
    double time1, double time1Elapsed, double &time2, double &totalTime,
    int truncateColumns, int truncateRows, bool redoSOS,
    double *&truncatedRhsLower, double *&truncatedRhsUpper,
    int *&newPriorities, bool integersOK, int numberOriginalColumns,
    int &returnCode, int callBack(CbcModel *currentSolver, int whereFrom),
    ampl_info *info);
  //@}
};

//#############################################################################
//#############################################################################

/*! \brief A class to allow the use of unknown user functionality

    For example, access to a modelling language (CbcAmpl).
*/
class CBCLIB_EXPORT CbcUser {

public:
  ///@name import/export methods
  //@{
  /*! \brief Import - gets full command arguments

      \return
      - -1 - no action
      -  0 - data read in without error
      -  1 - errors
    */
  virtual int importData(CbcSolver * /*model*/, int & /*argc*/, char ** /*argv[]*/)
  {
    return -1;
  }

  /*! \brief Export

      Values for mode:
      - 1 OsiClpSolver
      - 2 CbcModel
      - add 10 if infeasible from odd situation
    */
  virtual void exportSolution(CbcSolver * /*model*/,
    int /*mode*/, const char * /*message*/ = NULL) {}

  /// Export Data (i.e. at very end)
  virtual void exportData(CbcSolver * /*model*/) {}

  /// Get useful stuff
  virtual void fillInformation(CbcSolver * /*model*/,
    CbcParameters & /*info*/) {}
  //@}

  ///@name usage methods
  //@{
  /// CoinModel if valid
  inline CoinModel *coinModel() const
  {
    return coinModel_;
  }
  /// Other info - needs expanding
  virtual void *stuff()
  {
    return NULL;
  }
  /// Name
  inline std::string name() const
  {
    return userName_;
  }
  /// Solve (whatever that means)
  virtual void solve(CbcSolver *model, const char *options) = 0;
  /// Returns true if function knows about option
  virtual bool canDo(const char *options) = 0;
  //@}

  ///@name Constructors and destructors etc
  //@{
  /// Default Constructor
  CbcUser();

  /// Copy constructor
  CbcUser(const CbcUser &rhs);

  /// Assignment operator
  CbcUser &operator=(const CbcUser &rhs);

  /// Clone
  virtual CbcUser *clone() const = 0;

  /// Destructor
  virtual ~CbcUser();
  //@}

protected:
  ///@name Private member data
  //@{

  /// CoinModel
  CoinModel *coinModel_;

  /// Name of user function
  std::string userName_;

  //@}
};

//#############################################################################
//#############################################################################

/*! \brief Support the use of a call back class to decide whether to stop

  Definitely under construction.
*/

class CBCLIB_EXPORT CbcStopNow {

public:
  ///@name Decision methods
  //@{
  /*! \brief Import

      Values for whereFrom:
       - 1 after initial solve by dualsimplex etc
       - 2 after preprocessing
       - 3 just before branchAndBound (so user can override)
       - 4 just after branchAndBound (before postprocessing)
       - 5 after postprocessing
       - 6 after a user called heuristic phase

      \return 0 if good
       nonzero return code to stop
    */
  virtual int callBack(CbcModel * /*currentSolver*/, int /*whereFrom*/)
  {
    return 0;
  }
  //@}

  ///@name Constructors and destructors etc
  //@{
  /// Default Constructor
  CbcStopNow();

  /** Copy constructor. */
  CbcStopNow(const CbcStopNow &rhs);

  /// Assignment operator
  CbcStopNow &operator=(const CbcStopNow &rhs);

  /// Clone
  virtual CbcStopNow *clone() const;

  /// Destructor
  virtual ~CbcStopNow();
  //@}

private:
  ///@name Private member data
  //@{
  //@}
};

//###########################################################################
// Default no-op callback for CbcMain1 backward compatibility
//###########################################################################

static int dummyCallback(CbcModel * /*model*/, int /*whereFrom*/) { return 0; }

//#############################################################################
// Backward-compatible free functions (delegate to CbcSolver internally)
//#############################################################################

/// Initialize default parameters (delegates to CbcSolver::initialize)
CBCLIB_EXPORT
void CbcMain0(CbcModel &model, CbcParameters &parameters);

/// Run solver with command queue (delegates to CbcSolver::run)
CBCLIB_EXPORT
int CbcMain1(std::deque<std::string> inputQueue, CbcModel &model,
  CbcParameters &parameters,
  int callBack(CbcModel *currentSolver, int whereFrom) = dummyCallback,
  ampl_info *info = NULL);

void printGeneralMessage(CbcModel &model, std::string message, int type = CBC_GENERAL);
void printGeneralWarning(CbcModel &model, std::string message, int type = CBC_GENERAL_WARNING);
CBCLIB_EXPORT
int cbcReadAmpl(ampl_info *info, int argc, char **argv, CbcModel &model);

// for backward compatibility (samples)
#define CbcSolverUsefulData CbcParameters

/// Run solver with argc/argv (backward compatible)
CBCLIB_EXPORT
int CbcMain1(int argc, const char *argv[],
  CbcModel &model, CbcParameters &parameterData,
  int callBack(CbcModel *currentSolver, int whereFrom) = dummyCallback,
  ampl_info *info = NULL);

/// Run solver with argc/argv, alternate argument order (backward compatible)
CBCLIB_EXPORT
int CbcMain1(int argc, const char *argv[],
  CbcModel &model,
  int callBack(CbcModel *currentSolver, int whereFrom),
  CbcParameters &parameterData);

#endif // CbcSolver_H
