// Copyright (C) 2024, COIN-OR Foundation and others.
// All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

#ifndef CbcOutput_H
#define CbcOutput_H

/** \file CbcOutput.hpp
 *  \brief Compact, modular output for the Cbc solver header and problem summary.
 *
 *  Design principles:
 *  - Works identically for the CLI (CbcMain0/CbcMain1) and the C API.
 *  - Input-format agnostic: dimensions are queried from OsiSolverInterface.
 *  - Honors logLevel: messages print at logLevel >= 1 (the default).
 *  - Routes through the CoinMessageHandler chain, so C API callbacks intercept them.
 *  - UTF-8 mode (∈, κ, —) is auto-detected from the locale and can be overridden.
 */

#include "CbcConfig.h"
#include "CoinMessageHandler.hpp"
#include "CbcMessage.hpp"
#include "ClpEventHandler.hpp"

#include <climits>
#include <memory>
#include <string>
#include <vector>

class CbcModel;
class OsiSolverInterface;
class ClpSimplex;

// ---------------------------------------------------------------------------
// CbcPreprocHandler — intercepts CglPreProcess messages during preprocessing
// ---------------------------------------------------------------------------

/** Message handler installed on the CglPreProcess object before calling
 *  preProcessNonDefault().  Captures probing-pass statistics and renders
 *  them as a compact table, suppressing the raw Cgl0002/Cgl0003/Cgl0004
 *  message lines.
 *
 *  Usage (in CbcSolver.cpp):
 *    CbcPreprocHandler ph(fp, CbcOutput::useUtf8(), logLevel);
 *    process.passInMessageHandler(&ph);
 *    // ... call preProcessNonDefault() ...
 *    ph.printTableEnd();     // close table + processed model summary
 */
class CBCLIB_EXPORT CbcPreprocHandler : public CoinMessageHandler {
public:
  CbcPreprocHandler(FILE *fp, bool utf8, int logLevel);

  /** Destructor: auto-closes the preprocessing section if not already closed.
   *  This ensures the section banner is always printed even on early exits
   *  (infeasible detection, callback returns, etc.). */
  virtual ~CbcPreprocHandler();

  virtual int print() override;

  /** Close the probing table (bottom sep) and print the processed-model
   *  summary line.  Idempotent — safe to call multiple times. */
  void printTableEnd();

  /** Print the conflict-graph + clique-strengthening summary block.
   *  Pass cgraphBuilt=false when cgraph was not requested / not available. */
  void printCgraphSummary(bool cgraphBuilt, double cgraphTime,
    double cgraphDensity, bool clqRan, int clqExtended, int clqDominated);

  /** Print the overall preprocessing phase-end banner.
   *  Idempotent — safe to call multiple times (only prints once).
   *  totalTime is wall-clock seconds; if negative, computed from start time. */
  void printPhaseEnd(double totalTime = -1.0);

  /** Mark preprocessing as having detected infeasibility/unboundedness.
   *  The phase-end banner will reflect this when it is printed (explicitly or
   *  via the destructor). */
  void markInfeasible(const std::string &reason = "");

  // Accessors for the stored processed-model stats
  bool hasProcessedModel() const { return hasStats2_; }
  int processedRows() const { return procRows_; }
  int processedCols() const { return procCols_; }
  int processedInts() const { return procInts_; }
  int processedBinary() const { return procBinary_; }
  int processedNZ() const { return procNZ_; }
  bool tableWasPrinted() const { return headerPrinted_; }

private:
  void printTableHeader();
  void printTableRow(int pass, int fixed, int tightened, int strengthened,
    int subst);

  FILE *fp_;
  bool utf8_;
  bool compact_;
  double phaseStartTime_; // wall-clock time at construction
  bool headerPrinted_ = false;
  bool tableClosed_ = false; // set by printTableEnd()
  bool phaseClosed_ = false; // set by printPhaseEnd()
  bool infeasible_ = false;
  std::string infeasReason_;
  int passCount_ = 0;

  // Stored processed-model stats (from CGL_PROCESS_STATS2)
  bool hasStats2_ = false;
  int procRows_ = 0, procCols_ = 0, procInts_ = 0, procBinary_ = 0, procNZ_ = 0;
  int madeInteger_ = 0; // from CGL_MADE_INTEGER (ext=11)

  // SOS sets found during preprocessing (from CGL_PROCESS_SOS1, ext=5)
  struct SosInfo { int type; int count; int members; int integers; int overlaps; bool used; };
  std::vector<SosInfo> sosInfo_; // collected from CGL_PROCESS_SOS1/2 messages
};

// ---------------------------------------------------------------------------
// CbcImportHandler — captures MPS/LP parse errors during problem loading
// ---------------------------------------------------------------------------

/** Temporarily installed on the ClpSimplex model before readMps/readLp.
 *  Captures parse-error messages (detail=0, external numbers 3001-3006,
 *  6001-6005) and suppresses their immediate raw output.  The caller
 *  retrieves the captured list and displays it in the new formatted style.
 *
 *  At most MAX_STORED individual messages are stored; the total count is
 *  always tracked so callers can print "... and N more".
 *
 *  Usage:
 *    CbcImportHandler ih;
 *    lpSolver->passInMessageHandler(&ih);
 *    int status = clpSolver->readMps(fileName, keepNames, allowErrors);
 *    lpSolver->passInMessageHandler(originalHandler);
 *    CbcOutput::printImportErrors(fp, ih, status, utf8);
 */
class CBCLIB_EXPORT CbcImportHandler : public CoinMessageHandler {
public:
  static constexpr int MAX_STORED = 15;

  CbcImportHandler();
  virtual ~CbcImportHandler() = default;

  virtual int print() override;

  const std::vector<std::string> &errors() const { return errors_; }
  int totalErrors() const { return totalErrors_; }

private:
  std::vector<std::string> errors_;
  int totalErrors_ = 0;
};

// ---------------------------------------------------------------------------
// CbcRootLpEventHandler — installed on ClpSimplex during root node LP solve
// ---------------------------------------------------------------------------

/** Event handler for the root-node LP relaxation.
 *
 *  Installed on the ClpSimplex model before CbcModel::initialSolve() and
 *  removed afterwards.  Intercepts endOfIteration to print tabular progress
 *  rows, replacing Clp's built-in CLP_SIMPLEX_STATUS chatter (which is
 *  suppressed by setting the Clp log level to 0 during the root solve).
 *
 *  Clp may clone the event handler when it creates internal model copies
 *  during a single initialSolve() call (e.g. for idiot crash, barrier
 *  crossover, etc.).  All clones share the same SharedState via a
 *  std::shared_ptr, so the "header already printed" flag and timing are
 *  consistent across all internal sub-solves.
 *
 *  Frequency control:
 *   - iterFreq > 0 : print a row every iterFreq iterations
 *   - timeFreq > 0 : print a row every timeFreq seconds
 *  Whichever condition triggers first causes a row to be printed; the other
 *  counter resets at that point.  Both may be active simultaneously.
 *  Set both to 0 to suppress all per-iteration rows (still prints the header
 *  and final status line).
 */
class CBCLIB_EXPORT CbcRootLpEventHandler : public ClpEventHandler {
public:
  CbcRootLpEventHandler(CoinMessageHandler *msgHandler, int logLevel,
    int iterFreq, double timeFreq);

  virtual int event(Event whichEvent) override;
  virtual ClpEventHandler *clone() const override;

  /** Print the final status/result line.  Called once after initialSolve().
   *  numInts: total integer variables in the MIP (0 if pure LP).
   *  numFrac: number of fractional integer variables in the LP optimal solution. */
  void printFinalStatus(int numInts = 0, int numFrac = 0) const;

private:
  /** State shared across all Clp-internal clones of this handler. */
  struct SharedState {
    bool headerPrinted = false;
    double startTime = 0.0;
    double lastPrintTime = 0.0;
    int lastPrintIter = 0;
    int maxIterSeen = 0;   // tracks peak iteration count across all sub-solves
    int lastAlgo = 0;      // algorithm of last printed row (for phase label)
  };
  std::shared_ptr<SharedState> shared_;

  CoinMessageHandler *msgHandler_;
  int logLevel_;
  int iterFreq_;
  double timeFreq_;

  void printHeader() const;
  void printRow(int iter, double obj, double primalInf, double dualInf,
    double elapsed) const;
  static const char *algoName(int algo);
};

// ---------------------------------------------------------------------------
// CbcFPumpOutput — tabular output for the feasibility pump heuristic
// ---------------------------------------------------------------------------

/** Structured output handler for CbcHeuristicFPump.
 *
 *  Install on the heuristic via CbcHeuristicFPump::setFPumpOutput() before
 *  running branchAndBound(). The heuristic calls the onXxx/noteXxx methods
 *  at key events; when a handler is installed the heuristic suppresses its
 *  default CBC_FPUMP1 message output.
 *
 *  Frequency control:
 *   passFreq > 0 : print a row every passFreq passes
 *   timeFreq > 0 : print a row every timeFreq seconds
 *  If both are 0 every pass is printed. Solution events are always printed.
 */
class CBCLIB_EXPORT CbcFPumpOutput {
public:
  CbcFPumpOutput(FILE *fp, bool utf8, int logLevel,
    int passFreq = 0, double timeFreq = 5.0);
  ~CbcFPumpOutput() = default;

  bool isActive() const { return logLevel_ >= 1; }

  /** True while FPump is running (between onStart and onEnd). */
  bool isInPhase() const { return inPhase_; }

  // --- Events called by CbcHeuristicFPump::solutionInternal() ---

  /** Called once at the very start, with the initial fractional state. */
  void onStart(int numFrac, double suminf);

  /** Record a solution found during the current pass.
   *  Will show as a '*' row on the next onPass() call.
   *  type: "lp" | "rounding" | "minibab" | "cleaned" */
  void noteRowSolution(const std::string &type, double obj);

  /** Called after each FP pass.
   *  isOscillating: true when FPump's internal lastMove == 1000000
   *  (i.e. suminf went up this pass). */
  void onPass(int round, int pass, int numFrac, double suminf,
    double obj, int iters, bool isOscillating = false);

  /** Called at end of retry when rounding heuristic found a better solution. */
  void onRoundingImproved(double newObj);

  /** Called at end of retry when no solution was found. */
  void onNoSolutionInRetry();

  /** Called just before starting a new retry.
   *  newRound is the new round number (starts at 2). */
  void onRetry(int newRound, double newCutoff);

  /** Called at the very end of solutionInternal(). */
  void onEnd(double bestSol, double elapsed, int rounds, int totalPasses);

private:
  void printTableHeader(int numFrac, double suminf);
  void printTableEnd();
  void printRow(bool hasSolution, int round, int pass,
    const std::string &fracStr, const std::string &suminfStr,
    const std::string &status, double bestSol, double elapsed);
  bool shouldPrint(int pass, double elapsed) const;

  FILE *fp_;
  bool utf8_;
  bool compact_;
  int logLevel_;
  int passFreq_;
  double timeFreq_;

  bool headerPrinted_ = false;
  bool tableClosed_ = false;
  bool inPhase_ = false;
  double startTime_ = 0.0;
  double lastPrintTime_ = 0.0;
  int lastPrintPass_ = 0;
  int printedRows_ = 0;      ///< rows printed so far (first 10 are always shown)

  // Convergence tracking (updated on each printed row)
  double prevSuminf_ = 1e30;
  int prevNumFrac_ = INT_MAX;

  // Best solution seen so far
  double bestSol_ = 1e30;

  // Pending solution (set by noteRowSolution, consumed by onPass)
  bool pendingSolution_ = false;
  std::string pendingType_;
  double pendingObj_ = 1e30;

  // Last pass data (for forced prints)
  int lastRound_ = 0;
  int lastPass_ = 0;
  int lastNumFrac_val_ = 0;
  double lastSuminf_val_ = 0.0;
  bool lastPrinted_ = false;
};

// ---------------------------------------------------------------------------
// CbcCutGenOutput — tabular output for root-node cut generation
// ---------------------------------------------------------------------------

/** Formats cut generation output at the root node into:
 *  1. A per-pass progress table (pass, rows, tight cuts, objective, time).
 *  2. A per-generator summary table after all generators report (name, cuts,
 *     density, time, next-run frequency).
 *
 *  Messages are intercepted via CbcNautyHandler which calls the on*() methods.
 *  close() must be called after branchAndBound() returns to flush any pending
 *  generator table that hasn't been printed yet.
 */
class CBCLIB_EXPORT CbcCutGenOutput {
public:
  struct GenInfo {
    std::string name;
    int idx;
    int rowCuts;
    double avgDensity;
    int colCuts;
    double time;
    bool hasTime;
    int newFreq;
  };

  CbcCutGenOutput(FILE *fp, bool utf8, int logLevel);

  // Called by CbcNautyHandler for each relevant ext code
  void onStart();                                                   // ext=51
  void onPass(int pass, int rows, int tight, int frac, double suminf, double obj, double t); // ext=46
  void onSummary(int ncuts, double fromObj, double toObj, int passes); // ext=13
  void onGenerator(const GenInfo &g);                               // ext=14

  /** Flush pending generator table. Call after branchAndBound() returns.
   *  Idempotent — safe to call multiple times. */
  void close();

  /** Reset state to Idle so a second cut-generation phase can be printed.
   *  Optionally override the section title (default: "Cut generation (root node, part 2)"). */
  void resetForRestart(const char *title = nullptr);

  bool isInPhase() const { return state_ >= State::Started && state_ < State::Closed; }
  bool hasClosed()  const { return state_ == State::Closed; }
  /** True when generators have been accumulated but the table hasn't been flushed yet. */
  bool hasPendingGenerators() const { return state_ == State::Started && haveSummary_ && !genInfos_.empty() && !genTablePrinted_; }

private:
  enum class State { Idle, Started, Closed };

  void printProgressEnd();
  void printGeneratorTable();

  FILE *fp_;
  bool utf8_;
  bool compact_;
  int logLevel_;

  std::string title_ = "Cut generation (root node)";
  State state_ = State::Idle;
  bool firstPassSeen_ = false;

  // Progress table
  bool progHeaderPrinted_ = false;
  bool progTableClosed_ = false;

  // Summary from ext=13
  int sumNcuts_ = 0;
  double sumFromObj_ = 0.0;
  double sumToObj_ = 0.0;
  int sumPasses_ = 0;
  bool haveSummary_ = false;

  // Accumulated generator infos from ext=14
  std::vector<GenInfo> genInfos_;
  bool genTablePrinted_ = false;
};

// ---------------------------------------------------------------------------
// CbcOutput — static helper for header / problem-summary output
// ---------------------------------------------------------------------------

class CBCLIB_EXPORT CbcOutput {
public:
  /// Override UTF-8 output mode. Call before first print if needed.
  static void setUtf8(bool yesNo);

  /// Return the current UTF-8 mode (auto-detects from locale on first call).
  static bool useUtf8();

  /// Override compact-table mode (no box borders, thin rule separator).
  static void setCompact(bool yesNo);

  /// Return the current compact-table mode (default: false).
  static bool useCompact();

  /** Print the solver banner and optional args line:
   *    CBC devel (git:abc1234) — COIN-OR Branch and Cut
   *      args: instance.mps.gz -sec 180 -solve
   *
   *  Printed at logLevel >= 1. The args line is omitted when args is empty.
   */
  static void printSolverHeader(CoinMessageHandler *handler, int logLevel,
    const std::string &args = "", const std::string &strategyNote = "");

  /// Convenience overload using CbcModel's handler and log level.
  static void printSolverHeader(CbcModel &model, const std::string &args = "",
    const std::string &strategyNote = "");

  /** Print the multi-line problem summary:
   *    Problem: name
   *      |Rows| = N   |Cols| = N   |NZ| = N
   *      Variables: N binary, N integer, N continuous
   *    Coefficient ranges:
   *       Matrix  ∈ [lo, hi]   κ ≈ ratio
   *       Cost    ∈ [lo, hi]   κ ≈ ratio
   *       Bounds  ∈ [lo, hi]   κ ≈ ratio
   *       RHS     ∈ [lo, hi]   κ ≈ ratio
   *
   *  Printed at logLevel >= 1. No-op at logLevel 0.
   */
  static void printProblemSummary(CoinMessageHandler *handler,
    const OsiSolverInterface &solver, int logLevel,
    const CbcImportHandler *ih = nullptr);

  /// Convenience overload using CbcModel's handler and log level.
  static void printProblemSummary(CbcModel &model, const OsiSolverInterface &solver,
    const CbcImportHandler *ih = nullptr);

  /** Print captured parse errors/warnings in organized format.
   *  Call after restoring the original handler post readMps/readLp.
   *  Shows at most CbcImportHandler::MAX_STORED individual messages then
   *  "... and N more" if the total exceeds that.
   *  No-op if ih.totalErrors() == 0.
   */
  static void printImportErrors(FILE *fp, const CbcImportHandler &ih);

#ifdef CBC_HAS_NAUTY
  /** Print a clean symmetry-detection (nauty) section.
   *  Called from CbcSolver.cpp after running setupSymmetry/statsOrbits silently.
   *  errorCode == 0 means no error; non-zero prints an error message.
   */
  static void printNautySection(FILE *fp, bool utf8,
    int numUsefulOrbits, int numUsefulObjects,
    int totalOrbits, int numGenerators, double groupSize,
    double nautyTime, int errorCode = 0);
#endif

private:
  static bool utf8_;
  static bool utf8Set_;
  static bool compact_;
  static bool detectUtf8();
};

// ---------------------------------------------------------------------------
// CbcBnBOutput — tabular B&B tree progress output
// ---------------------------------------------------------------------------

/** Formats branch-and-bound progress into:
 *  1. A per-update progress table (Nodes/OnTree/BestSol/Method/BestBound/Gap%/Time).
 *  2. ★-prefixed rows whenever a new incumbent is found (always printed).
 *  3. A footer with strong-branching, depth, and orbital-branching stats.
 *
 *  Driven by CbcNautyHandler::print() which calls the on*() methods.
 */
class CBCLIB_EXPORT CbcBnBOutput {
public:
  CbcBnBOutput(FILE *fp, bool utf8, int logLevel);
  ~CbcBnBOutput();

  // ext=37 (CBC_STATUS2): periodic progress with depth
  void onProgress(long nodes, int onTree, int depth, double bestSol, double bestBound, double elapsed);

  // ext=4 (CBC_SOLUTION): incumbent found by B&B LP
  void onBnBIncumbent(double obj, long nodes, double elapsed);

  // ext=12 (CBC_ROUNDING): incumbent found by a named heuristic
  void onHeurIncumbent(double obj, const char *method, long nodes, double elapsed);

  // Queue an incumbent found before B&B starts (e.g. during root cut gen).
  // It will be emitted as a ★ row at the top of the B&B table once the
  // first progress message arrives and the bound is known.
  void queuePreProgressIncumbent(double obj, const char *method,
    long nodes, double elapsed);

  // ext=3/19/20/50: stopping reason (captured before onComplete)
  void onStopReason(const char *reason);

  // ext=1 (CBC_END_GOOD) / ext=5 (CBC_END): B&B complete
  void onComplete(bool optimal, double bestSol, double bestBound,
    long iters, long nodes, double elapsed);

  // Accumulated end stats (printed after ✔ line)
  void onStrongStats(long calls, long iters, long fathomed, int fixed);      // ext=32
  void onOtherStats(int maxDepth, double djFixed);                           // ext=35
  void onOtherStats2(int maxDepth, double djFixed,                           // ext=41
    long fathomTimes, long fathomNodes, long fathomIters);
  void onOrbitalStats(int successes, double avgExtra, int fixed,             // ext=45
    double avgFixed);

  // ext=44: reduced-cost fixing restart — show separator in B&B table
  void onRestartNote(int rows, int cols);

  bool isInPhase() const { return inPhase_; }

private:
  void startPhase();
  void openContinuation();   ///< Re-open B&B table after a restart
  void printRow(bool isIncumbent, long nodes, int onTree, int depth,
    double bestSol, const char *method, double bestBound, double wallclock);
  void closeTable();
  void printFooter();

  FILE *fp_;
  bool utf8_;
  bool compact_;
  int  logLevel_;

  bool inPhase_   = false;
  bool tableOpen_ = false;
  int  restartCount_ = 0;  ///< Number of restarts seen (drives "continued" header)

  // State from last ext=37 for ext=4/12 rows (which lack onTree/bestBound)
  int    lastOnTree_    = 0;
  int    lastDepth_     = 0;
  double lastBestBound_ = 1e50;

  // Incumbents found before the first onProgress() call (e.g. during root cut gen).
  // Emitted as ★ rows at the top of the B&B table, using the bound from
  // the first progress message.
  struct PendingHeurRow {
    double obj;
    std::string method;
    long nodes;
    double wallclock;   ///< CoinWallclockTime() at time of queuing
  };
  std::vector<PendingHeurRow> preProgressIncumbents_;

  std::string stopReason_;

  // Strong-branching stats (ext=32)
  bool hasStrong_   = false;
  long sbCalls_     = 0;
  long sbIters_     = 0;
  long sbFathomed_  = 0;
  int  sbFixed_     = 0;

  // Depth / DJ stats (ext=35 / ext=41)
  bool   hasDepth_    = false;
  int    maxDepth_    = 0;
  double djFixed_     = 0.0;
  bool   hasFathom_   = false;
  long   fathomTimes_ = 0;
  long   fathomNodes_ = 0;
  long   fathomIters_ = 0;

  // Orbital-branching stats (ext=45 "Orbital branching succeeded...")
  bool   hasOrbital_   = false;
  int    orbSuccesses_ = 0;
  double orbAvgExtra_  = 0.0;
  int    orbFixed_     = 0;
  double orbAvgFixed_  = 0.0;
};

// ---------------------------------------------------------------------------
// CbcNautyHandler — intercepts Nauty messages inside branchAndBound()
// ---------------------------------------------------------------------------

#ifdef CBC_HAS_NAUTY
/** Message handler temporarily installed on CbcModel before branchAndBound().
 *  Intercepts CBC_GENERAL messages containing "Nauty" and formats them into a
 *  clean symmetry-detection section.  All other messages are forwarded through
 *  the normal CoinMessageHandler::print() path.
 *
 *  Usage:
 *    CbcNautyHandler nh(fp, CbcOutput::useUtf8(), logLevel);
 *    CoinMessageHandler *saved = babModel_->messageHandler();
 *    babModel_->passInMessageHandler(&nh);
 *    babModel_->branchAndBound(doStatistics);
 *    babModel_->passInMessageHandler(saved);
 */
class CBCLIB_EXPORT CbcNautyHandler : public CoinMessageHandler {
public:
  CbcNautyHandler(FILE *fp, bool utf8, int logLevel);
  virtual ~CbcNautyHandler();

  /** Return a permanently-owned silent CoinMessageHandler suitable for use
   *  as an LP-solver message handler in the restart sub-model.  This keeps
   *  the CbcNautyHandler's own logLevel safe when CbcStrategyDefault's
   *  setupPrinting() calls solver->messageHandler()->setLogLevel(0). */
  CoinMessageHandler *getLpSilentHandler();
  virtual int print() override;

  /** Set the FPump output handler so we can suppress the raw Cbc0012I
   *  "found by feasibility pump" message while the FPump table is active. */
  void setFPumpOutput(CbcFPumpOutput *fp) { fpumpOut_ = fp; }

  /** Set the cut-gen output handler so root-node cut messages are reformatted. */
  void setCutGenOutput(CbcCutGenOutput *cg) { cutGenOut_ = cg; }

  /** Set the B&B output handler for tree progress reformatting. */
  void setBnBOutput(CbcBnBOutput *bnb) { bnbOut_ = bnb; }

  /** Prepare handler for a B&B restart: reset cut-gen output for a part-2
   *  run and set the restartMode_ flag so the sub-model's search-complete
   *  message is suppressed (it must not trigger bnbOut_->onComplete()). */
  void beginRestartMode();

  /** Returns true if at least one Nauty message was intercepted. */
  bool didRun() const { return sectionStarted_; }

private:
  void printSection();
  // Parse an "Integer solution of ... found by ... after ... nodes (...)"
  // message and route it to bnbOut_.
  void routeIncumbentMessage(const char *buf, int ext);

  FILE *fp_;
  bool utf8_;
  bool compact_;
  bool sectionStarted_ = false;
  bool restartMode_    = false; ///< True while intercepting a B&B restart sub-model
  CoinMessageHandler *lpSilentHandler_ = nullptr; ///< Silent handler lent to LP solver
  CbcFPumpOutput *fpumpOut_ = nullptr;   // for Cbc0012I FPump suppression
  CbcCutGenOutput *cutGenOut_ = nullptr; // for root cut generation output
  CbcBnBOutput *bnbOut_ = nullptr;       // for B&B progress table

  // Stats accumulated from the "Nauty: N orbits ..." message
  int totalOrbits_ = 0;
  int usefulOrbits_ = 0;
  int usefulVars_ = 0;
  int numGens_ = 0;
  double groupSize_ = 0.0;
  double nautyTime_ = 0.0;
  bool hasError_ = false;
  int errorCode_ = 0;
};
#endif /* CBC_HAS_NAUTY */

#endif /* CbcOutput_H */

/* vi: softtabstop=2 shiftwidth=2 expandtab tabstop=2
*/
