// Copyright (C) 2024, COIN-OR Foundation
// All Rights Reserved. This code is published under the Eclipse Public License.

#ifndef ClpOutput_HPP
#define ClpOutput_HPP

#include "ClpConfig.h"
#include "ClpSimplex.hpp"
#include "ClpEventHandler.hpp"
#include "CoinMessageHandler.hpp"
#include "CoinTime.hpp"

#include <memory>
#include <string>

// ─── ClpProgressEventHandler ──────────────────────────────────────────────────
// Event handler for LP progress output in the Clp standalone solver.
// Installed on the ClpSimplex model before solve and removed afterwards.
// Intercepts endOfIteration to print tabular progress, replacing Clp's built-in
// CLP_SIMPLEX_STATUS chatter.
//
// Presolve support: pass origRows/origCols/origNZ from the original model.
// When the first iteration fires in a smaller (presolved) model, the header
// shows the presolve reduction.
// ──────────────────────────────────────────────────────────────────────────────
class CLPLIB_EXPORT ClpProgressEventHandler : public ClpEventHandler {
public:
  ClpProgressEventHandler(CoinMessageHandler *handler, int logLevel,
    int iterFreq, double timeFreq,
    int origRows, int origCols, int origNZ);

  virtual int event(Event whichEvent) override;
  virtual ClpEventHandler *clone() const override;

  void printFinalStatus() const;

private:
  struct SharedState {
    bool headerPrinted = false;
    double startTime = 0.0;
    double lastPrintTime = 0.0;
    int lastPrintIter = 0;
    int maxIterSeen = 0;
  };
  std::shared_ptr<SharedState> shared_;

  CoinMessageHandler *handler_;
  int logLevel_;
  int iterFreq_;
  double timeFreq_;
  int origRows_;
  int origCols_;
  int origNZ_;

  void printHeader() const;
  void printRow(int iter, double obj, double pInf, double dInf, double elapsed, int algo) const;
};

// ─── ClpLpPhaseState ──────────────────────────────────────────────────────────
// Shared mutable state for ClpLpEventHandler and ClpLpMsgHandler.
// Both handlers hold a std::shared_ptr to the same instance so that data
// written by the message handler (Idiot/Sprint rows) is immediately visible
// to the event handler (LP rows), and vice-versa.
// ──────────────────────────────────────────────────────────────────────────────
struct ClpLpPhaseState {
  // ── Table bookkeeping ─────────────────────────────────────────────────────
  bool headerPrinted = false; // table header has been printed
  bool lpStarted = false;     // first LP endOfIteration has fired
  double startTime = 0.0;
  double lastPrintTime = 0.0;
  int lastPrintIter = 0;
  int maxIterSeen = 0;

  // ── Idiot: last-row tracking (force-print when LP begins or solve ends) ───
  int lastIdiotIter = -1;
  double lastIdiotInfeas = 0.0, lastIdiotObj = 0.0, lastIdiotTime = 0.0;
  bool idiotSeen = false;  // first Idiot message received
  bool sprintSeen = false; // first Sprint message received
  bool idiotPending = false; // last received Idiot row not yet printed

  // ── Sprint: accumulated LP iters + last-row tracking ─────────────────────
  int sprintCumIters = 0;
  int lastSprintCumIters = -1;
  double lastSprintObj = 0.0, lastSprintDInf = 0.0, lastSprintTime = 0.0;
  bool sprintPending = false; // last received Sprint row not yet printed

  // ── Configuration (set at construction, read-only afterwards) ─────────────
  FILE *fp = nullptr;
  bool utf8 = false;
  bool compact = true;
  int iterFreq = 0;
  double timeFreq = 5.0;
  int logLevel = 1;
  int origRows = 0, origCols = 0;
  std::string title = "LP solve";
};

// ─── ClpLpEventHandler ────────────────────────────────────────────────────────
// ClpEventHandler that drives a unified LP+Idiot+Sprint progress table.
//
// Install on a ClpSimplex model via passInEventHandler() (Clp clones it).
// Always pair with a ClpLpMsgHandler installed as the model's message handler
// so that Idiot and Sprint messages appear as rows in the same table.
//
// Usage:
//   auto state = std::make_shared<ClpLpPhaseState>();
//   // fill state->fp, utf8, logLevel, iterFreq, timeFreq, title ...
//   ClpLpEventHandler evtH(state);
//   ClpLpMsgHandler   msgH(state);
//   bool msgOldDefault;
//   CoinMessageHandler *savedMsg = model->pushMessageHandler(&msgH, msgOldDefault);
//   model->passInEventHandler(&evtH);
//   auto *evtClone = dynamic_cast<ClpLpEventHandler *>(model->eventHandler());
//   ... model->initialSolve() ...
//   if (evtClone) evtClone->printFinalStatus(numInts, numFrac);
//   model->popMessageHandler(savedMsg, msgOldDefault);
//   ClpEventHandler defEvt; model->passInEventHandler(&defEvt);
// ──────────────────────────────────────────────────────────────────────────────
class CLPLIB_EXPORT ClpLpEventHandler : public ClpEventHandler {
public:
  explicit ClpLpEventHandler(std::shared_ptr<ClpLpPhaseState> state);

  virtual int event(Event whichEvent) override;
  virtual ClpEventHandler *clone() const override;

  /** Print final status line + close the table.
   *  Call once after initialSolve() on the clone returned by eventHandler().
   *  numInts: total MIP integer variables (0 for pure LP).
   *  numFrac: fractional integer variables in LP optimal (0 for pure LP). */
  void printFinalStatus(int numInts = 0, int numFrac = 0);

  /// Whether the progress table has already been opened (Idiot/Sprint started)
  bool tableStarted() const { return s_ && (s_->idiotSeen || s_->sprintSeen || s_->lpStarted); }
  /// File pointer for output
  FILE *fp() const { return s_ ? s_->fp : nullptr; }

private:
  std::shared_ptr<ClpLpPhaseState> s_;

  void openTable();
  void printLpRow(int iter, double obj, double pInf, double dInf, double elapsed);
  /** Force-print any pending (unshown) last Idiot/Sprint row. */
  void flushPendingIdiotSprint();
};

// ─── ClpLpMsgHandler ──────────────────────────────────────────────────────────
// CoinMessageHandler that intercepts CLP_IDIOT_ITERATION (ext 30) and
// CLP_SPRINT (ext 34), emitting them as "Idiot" / "Sprint" rows into the
// unified LP table driven by ClpLpEventHandler.
//
// All other Clp messages are suppressed; detail-0 critical errors are forwarded.
//
// IMPORTANT: The Clp model's log level must NOT be set to 0 while Idiot is
// active — Idiot checks the model log level before constructing messages at all.
// ──────────────────────────────────────────────────────────────────────────────
class CLPLIB_EXPORT ClpLpMsgHandler : public CoinMessageHandler {
public:
  explicit ClpLpMsgHandler(std::shared_ptr<ClpLpPhaseState> state);

  virtual int print() override;
  virtual CoinMessageHandler *clone() const override;

private:
  std::shared_ptr<ClpLpPhaseState> s_;

  void ensureTableOpen();
  void printIdiotRow(int iter, double infeas, double obj, double elapsed);
  void printSprintRow(int cumIters, double obj, double dInf, double elapsed);
};

// ─── ClpOutput ────────────────────────────────────────────────────────────────
// Centralised clean output for the Clp standalone LP solver.
//
//  • printSolverHeader   — one-line "CLP x.y.z — COIN-OR LP Solver" + args
//  • printProblemSummary — compact problem description + coefficient-range table
//
// All output is written via the provided CoinMessageHandler's filePointer().
// ──────────────────────────────────────────────────────────────────────────────

class CLPLIB_EXPORT ClpOutput {
public:
  /// Print the one-line solver header and (if argc>1) an args line.
  static void printSolverHeader(CoinMessageHandler *handler, int logLevel,
    int argc = 0, const char *const *argv = nullptr);

  /// Print compact problem summary (dimensions + coefficient ranges).
  static void printProblemSummary(CoinMessageHandler *handler,
    const ClpSimplex &model, int logLevel);

  // ── UTF-8 toggle ──────────────────────────────────────────────────────────
  static void setUtf8(bool on) { useUtf8_ = on; }
  static bool useUtf8() { return useUtf8_; }

  // ── Compact table style toggle ────────────────────────────────────────────
  static void setCompact(bool on) { useCompact_ = on; }
  static bool useCompact() { return useCompact_; }

private:
  static bool useUtf8_;
  static bool useCompact_;

  // Helpers
  static void sendLine(CoinMessageHandler *handler, const std::string &line);
  static void sendBlank(CoinMessageHandler *handler);
  static std::string fmtNum(double v);

  struct RangeStat {
    double lo = 1e300, hi = -1e300;
    bool valid = false;
    void update(double v);
    std::string ratio() const;
  };
};

#endif // ClpOutput_HPP
