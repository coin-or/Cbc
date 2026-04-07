// Copyright (C) 2024, COIN-OR Foundation and others.
// All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

#include "CbcOutput.hpp"
#include "CbcModel.hpp"
#include "OsiSolverInterface.hpp"
#include "OsiClpSolverInterface.hpp"
#include "CoinPackedMatrix.hpp"
#include "CbcConfig.h"
#include "ClpSimplex.hpp"
#include "CoinTable.hpp"
#include "CoinTime.hpp"
#ifdef CBC_HAS_NAUTY
#include "CbcSymmetry.hpp"
#endif

#include <algorithm>
#include <climits>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <iomanip>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

// ---------------------------------------------------------------------------
// Shared output helpers
// ---------------------------------------------------------------------------

// Table column separator (│ or | or two spaces in compact mode).
static inline const char *tableBar(bool utf8, bool compact)
{
  if (compact)
    return " ";
  return utf8 ? " \xe2\x94\x82 " : " | ";
}

// UTF-8 visual width: counts Unicode code points (not bytes) for column alignment.
static int utf8VisLen(const std::string &s)
{
  int n = 0;
  for (unsigned char c : s)
    if ((c & 0xC0) != 0x80) ++n;
  return n;
}

// Print a table header block and flush.
//   Bordered: top sep → header → middle sep (3 lines)
//   Compact:  header → separator rule (2 lines, no border frame)
static void printTableOpen(FILE *fp, const CoinTable &tbl)
{
  if (tbl.compact()) {
    fprintf(fp, "%s\n", tbl.headerLine().c_str());
    fprintf(fp, "%s\n", tbl.sepLine().c_str());
  } else {
    fprintf(fp, "%s\n", tbl.sepLine(CoinTable::Top).c_str());
    fprintf(fp, "%s\n", tbl.headerLine().c_str());
    fprintf(fp, "%s\n", tbl.sepLine(CoinTable::Middle).c_str());
  }
  fflush(fp);
}

// Print the bottom-close separator and flush (no-op in compact mode).
static void printTableClose(FILE *fp, const CoinTable &tbl)
{
  if (!tbl.compact())
    fprintf(fp, "%s\n", tbl.sepLine(CoinTable::Bottom).c_str());
  fflush(fp);
}

// Format a time value (seconds) with dynamic precision:
//   t < 1   → 3 decimal places  e.g. "0.001"
//   t < 10  → 2 decimal places  e.g. "1.23"
//   t < 100 → 1 decimal place   e.g. "43.2"
//   t >= 100→ 0 decimal places  e.g. "99876"
static std::string fmtTime(double t)
{
  char buf[32];
  if (t < 1.0)        std::snprintf(buf, sizeof(buf), "%.3f", t);
  else if (t < 10.0)  std::snprintf(buf, sizeof(buf), "%.2f", t);
  else if (t < 100.0) std::snprintf(buf, sizeof(buf), "%.1f", t);
  else                std::snprintf(buf, sizeof(buf), "%.0f", t);
  return buf;
}

// Parsed fields from "Integer solution of OBJ found by METHOD after NITERS iterations and NODES nodes (ELAPSED"
struct IncumbentMsg {
  double obj = 0.0;
  std::string method;
  long nodes = 0;
  double elapsed = 0.0;
};

// Parse a CBC incumbent solution message into an IncumbentMsg.
// Returns true on success.
static bool parseIncumbentMsg(const char *buf, IncumbentMsg &out)
{
  const char *p = std::strstr(buf, "Integer solution of ");
  if (!p) return false;
  if (std::sscanf(p, "Integer solution of %lf", &out.obj) != 1) return false;
  const char *byP  = std::strstr(p, " found by ");
  const char *aftP = std::strstr(p, " after ");
  out.method = (byP && aftP && aftP > byP + 10) ? std::string(byP + 10, aftP) : "B&B";
  if (!aftP) return false;
  return std::sscanf(aftP, " after %*ld iterations and %ld nodes (%lf",
    &out.nodes, &out.elapsed) == 2;
}

// ---------------------------------------------------------------------------
// CbcPreprocHandler — probing-table message handler
// ---------------------------------------------------------------------------

// Column widths for the preprocessing probing table:
//
//   Pass │  Fixed │ Tightened │ Strengthened │  Subst │  Time
//   ─────┼────────┼───────────┼──────────────┼────────┼──────
//      1 │    480 │         0 │         7157 │     76 │ 0.12s
//
static const int PP_W_PASS  = 5;
static const int PP_W_FIXED = 7;
static const int PP_W_TIGHT = 9;
static const int PP_W_STR   = 12;
static const int PP_W_SUBST = 6;
static const int PP_W_TIME  = 7;

static CoinTable makePpTable(bool utf8, bool compact)
{
  return CoinTable({
    { "Pass",         PP_W_PASS  },
    { "Fixed",        PP_W_FIXED },
    { "Tightened",    PP_W_TIGHT },
    { "Strengthened", PP_W_STR   },
    { "Subst",        PP_W_SUBST },
    { "Time(s)",         PP_W_TIME  },
  }, utf8, /*indent=*/2, compact);
}

CbcImportHandler::CbcImportHandler()
  : CoinMessageHandler()
{
  // Suppress all immediate output — we collect and display later.
  setFilePointer(nullptr);
  setLogLevel(0);
  setPrefix(false); // no "CoinNNNNW " prefix in buffered text
}

int CbcImportHandler::print()
{
  const int ext = currentMessage().externalNumber();
  // Capture detail-0 MPS/LP parse messages by external number:
  //   3001 COIN_MPS_ILLEGAL     3002 COIN_MPS_BADIMAGE
  //   3003 COIN_MPS_DUPOBJ      3004 COIN_MPS_DUPROW
  //   3005 COIN_MPS_NOMATCHROW  3006 COIN_MPS_NOMATCHCOL
  //   6001 COIN_MPS_FILE        6002 COIN_MPS_BADFILE1
  //   6003 COIN_MPS_BADFILE2    6004 COIN_MPS_EOF
  //   6005 COIN_MPS_RETURNING
  bool isParseDiag = (ext >= 3001 && ext <= 3006)
                  || (ext >= 6001 && ext <= 6005);
  if (isParseDiag) {
    totalErrors_++;
    if ((int)errors_.size() < MAX_STORED) {
      const char *buf = messageBuffer();
      if (buf && *buf)
        errors_.emplace_back(buf);
    }
  }
  return 0; // always suppress immediate printing
}

CbcPreprocHandler::CbcPreprocHandler(FILE *fp, bool utf8, int logLevel)
  : CoinMessageHandler(fp)
  , fp_(fp)
  , utf8_(utf8)
  , compact_(CbcOutput::useCompact())
  , phaseStartTime_(CoinWallclockTime())
{
  setLogLevel(logLevel);
}

CbcPreprocHandler::~CbcPreprocHandler()
{
  // Auto-close the section on any exit path (infeasible break, early return, etc.)
  if (!phaseClosed_)
    printPhaseEnd();
}

void CbcPreprocHandler::printTableHeader()
{
  if (!fp_) return;
  printTableOpen(fp_, makePpTable(utf8_, compact_));
}

void CbcPreprocHandler::printTableRow(int pass, int fixed, int tightened,
  int strengthened, int subst)
{
  if (!fp_) return;
  const char *bar = tableBar(utf8_, compact_);
// Preprocessing time column: absolute wall-clock since program start.
  const std::string timeStr = fmtTime(CoinWallclockTime());
  fprintf(fp_, "  %*d%s%*d%s%*d%s%*d%s%*d%s%*s\n",
    PP_W_PASS,  pass,         bar,
    PP_W_FIXED, fixed,        bar,
    PP_W_TIGHT, tightened,    bar,
    PP_W_STR,   strengthened, bar,
    PP_W_SUBST, subst,        bar,
    PP_W_TIME,  timeStr.c_str());
  fflush(fp_);
}

int CbcPreprocHandler::print()
{
  const CoinOneMessage cm = currentMessage();
  std::string src = currentSource();
  const int ext = cm.externalNumber();
  const char *buf = messageBuffer();

  // CGL_MADE_INTEGER (ext=11): "N variables made integer"
  // Buffer may have prefix "Cgl0011I N variables made integer"
  if (src == "Cgl" && ext == 11) {
    const char *p = std::strstr(buf, "variables made integer");
    while (p && p > buf && *(p - 1) == ' ') --p;
    while (p && p > buf && std::isdigit(static_cast<unsigned char>(*(p - 1)))) --p;
    if (p) std::sscanf(p, "%d", &madeInteger_);
    return 0; // suppress; shown in printTableEnd()
  }

  // CGL_FIXED (ext=2): "N variables fixed" — suppress, will show in summary
  if (src == "Cgl" && ext == 2)
    return 0;

  // CGL_PROCESS_STATS (ext=3): "N fixed, N tightened bounds, N strengthened rows, N substitutions"
  if (src == "Cgl" && ext == 3) {
    int fixed = 0, tightened = 0, strengthened = 0, subst = 0;
    // messageBuffer() may contain a "Cgl0003I " prefix; find the data part
    const char *p = std::strstr(buf, " fixed,");
    if (!p) p = buf;
    else { while (p > buf && *(p-1) != ' ' && *(p-1) != '\t') --p; }
    std::sscanf(p, "%d fixed, %d tightened bounds, %d strengthened rows, %d substitutions",
      &fixed, &tightened, &strengthened, &subst);
    if (!headerPrinted_) {
      printTableHeader();
      headerPrinted_ = true;
    }
    passCount_++;
    printTableRow(passCount_, fixed, tightened, strengthened, subst);
    return 0;
  }

  // CGL_PROCESS_STATS2 (ext=4): "processed model has R rows, C columns (I integer (B of which binary)) and E elements"
  if (src == "Cgl" && ext == 4) {
    const char *p = std::strstr(buf, "processed model has ");
    if (p)
      std::sscanf(p,
        "processed model has %d rows, %d columns (%d integer (%d of which binary)) and %d elements",
        &procRows_, &procCols_, &procInts_, &procBinary_, &procNZ_);
    hasStats2_ = true;
    return 0; // suppress; caller will print via printTableEnd()
  }

  // CGL_PROCESS_SOS1 (ext=5): "%d SOS with %d members" — capture for summary line
  if (src == "Cgl" && ext == 5) {
    int count = 0, members = 0;
    // Walk back from "SOS " to the number before it
    const char *p = std::strstr(buf, "SOS ");
    if (p) {
      while (p > buf && *(p-1) == ' ') --p;
      while (p > buf && std::isdigit(static_cast<unsigned char>(*(p-1)))) --p;
      std::sscanf(p, "%d SOS with %d members", &count, &members);
    }
    if (count > 0) {
      SosInfo si; si.type=1; si.count=count; si.members=members;
      si.integers=0; si.overlaps=0; si.used=true;
      sosInfo_.push_back(si);
    }
    return 0;
  }

  // CGL_PROCESS_SOS2 (ext=6): "%d SOS (%d members out of %d) with %d overlaps - too much overlap..."
  // Fires when SOS sets are found but cannot be used (overlap/integer ratio).
  if (src == "Cgl" && ext == 6) {
    int count = 0, inSOS = 0, integers = 0, overlaps = 0;
    const char *p = std::strstr(buf, "SOS ");
    if (p) {
      while (p > buf && *(p-1) == ' ') --p;
      while (p > buf && std::isdigit(static_cast<unsigned char>(*(p-1)))) --p;
      std::sscanf(p, "%d SOS (%d members out of %d) with %d overlaps",
        &count, &inSOS, &integers, &overlaps);
    }
    if (count > 0) {
      SosInfo si; si.type=2; si.count=count; si.members=inSOS;
      si.integers=integers; si.overlaps=overlaps; si.used=false;
      sosInfo_.push_back(si);
    }
    return 0;
  }

  // All other messages: default output
  return CoinMessageHandler::print();
}

void CbcPreprocHandler::printTableEnd()
{
  if (tableClosed_) return;
  tableClosed_ = true;
  if (!fp_) return;
  if (headerPrinted_)
    printTableClose(fp_, makePpTable(utf8_, compact_));
  if (hasStats2_) {
    fprintf(fp_, "\n  Processed model: %d rows, %d cols, %d NZ\n",
      procRows_, procCols_, procNZ_);
    if (procInts_ > 0) {
      fprintf(fp_, "                   %d integer (%d binary)\n",
        procInts_, procBinary_);
    }
    if (madeInteger_ > 0) {
      fprintf(fp_, "                   %d continuous variable(s) made integer\n",
        madeInteger_);
    }
  }
  // Blank line after processed model block (before SOS / conflict graph lines)
  if (hasStats2_)
    fprintf(fp_, "\n");
  // SOS summary: one line per SOS type found
  for (const auto &s : sosInfo_) {
    if (s.used) {
      // SOS1 sets that are active
      double avg = s.count > 0 ? static_cast<double>(s.members) / s.count : 0.0;
      fprintf(fp_, "  SOS%d: %d sets, %d members (avg %.1f members/set)\n",
        s.type, s.count, s.members, avg);
    } else {
      // SOS sets found but not used due to overlap / integer ratio
      fprintf(fp_, "  SOS%d: %d sets found (%d members) — not used"
        " (%d overlaps, %d integer vars; too much overlap)\n",
        s.type, s.count, s.members, s.overlaps, s.integers);
    }
  }
  if (hasStats2_ || !sosInfo_.empty())
    fflush(fp_);
}

void CbcPreprocHandler::printCgraphSummary(bool cgraphBuilt, double cgraphTime,
  double cgraphDensity, bool clqRan, int clqExtended, int clqDominated)
{
  if (!fp_) return;
  if (cgraphBuilt) {
    fprintf(fp_, "  Conflict graph: built in %.2fs, density %.3f%%\n",
      cgraphTime, cgraphDensity * 100.0);
  }
  if (clqRan) {
    fprintf(fp_, "  Clique strengthening: extended %d, dominated %d\n",
      clqExtended, clqDominated);
  }
  fflush(fp_);
}

void CbcPreprocHandler::markInfeasible(const std::string &reason)
{
  infeasible_ = true;
  infeasReason_ = reason.empty() ? "infeasible or unbounded" : reason;
}

void CbcPreprocHandler::printPhaseEnd(double totalTime)
{
  if (phaseClosed_) return;
  phaseClosed_ = true;
  printTableEnd(); // idempotent — close table if not already done
  if (!fp_) return;
  if (totalTime < 0.0)
    totalTime = CoinWallclockTime() - phaseStartTime_;
  char summary[160];
  const std::string tStr = fmtTime(totalTime);
  if (infeasible_) {
    std::snprintf(summary, sizeof(summary),
      "Preprocessing infeasible%sTime: %ss",
      CoinTable::dashSep(utf8_), tStr.c_str());
  } else {
    std::snprintf(summary, sizeof(summary),
      "Preprocessing complete%sTime: %ss",
      CoinTable::dashSep(utf8_), tStr.c_str());
  }
  fprintf(fp_, "\n%s\n", CoinTable::phaseEnd(summary, utf8_).c_str());
  fflush(fp_);
}

// ---------------------------------------------------------------------------
// Static state
// ---------------------------------------------------------------------------

bool CbcOutput::utf8_ = false;
bool CbcOutput::utf8Set_ = false;
bool CbcOutput::compact_ = true;

void CbcOutput::setUtf8(bool yesNo)
{
  utf8_ = yesNo;
  utf8Set_ = true;
}

void CbcOutput::setCompact(bool yesNo)
{
  compact_ = yesNo;
}

bool CbcOutput::useCompact()
{
  return compact_;
}

bool CbcOutput::detectUtf8()
{
  const char *lang = std::getenv("LC_ALL");
  if (!lang || !*lang)
    lang = std::getenv("LC_CTYPE");
  if (!lang || !*lang)
    lang = std::getenv("LANG");
  if (!lang)
    return false;
  // Check for "UTF-8" or "utf8" suffix (case-insensitive via strstr)
  std::string s(lang);
  for (char &c : s)
    c = (char)std::tolower((unsigned char)c);
  return s.find("utf-8") != std::string::npos || s.find("utf8") != std::string::npos;
}

bool CbcOutput::useUtf8()
{
  if (!utf8Set_) {
    utf8_ = detectUtf8();
    utf8Set_ = true;
  }
  return utf8_;
}

// ---------------------------------------------------------------------------
// CbcRootLpEventHandler — tabular LP progress for the root node
// ---------------------------------------------------------------------------

// Column widths for the progress table:
//
//   Phase   │     Iter │       Objective │   Primal inf │    Dual inf │    Time
//   ─────────┬──────────┬─────────────────┬──────────────┬─────────────┬─────────
//   Dual     │     3123 │     3.709408e+01 │  1.7376e+02  │  2.8598e+14 │    5.0s
//
static const int LP_W_PHASE = 8;
static const int LP_W_ITER  = 8;
static const int LP_W_OBJ   = 15;
static const int LP_W_PINF  = 12;
static const int LP_W_DINF  = 12;
static const int LP_W_TIME  = 8;

static CoinTable makeLpTable(bool utf8, bool compact)
{
  return CoinTable({
    { "Phase",      LP_W_PHASE, /*leftAlign=*/true },
    { "Iter",       LP_W_ITER  },
    { "Objective",  LP_W_OBJ   },
    { "Primal inf", LP_W_PINF  },
    { "Dual inf",   LP_W_DINF  },
    { "Time(s)",       LP_W_TIME  },
  }, utf8, /*indent=*/2, compact);
}

CbcRootLpEventHandler::CbcRootLpEventHandler(CoinMessageHandler *msgHandler,
  int logLevel, int iterFreq, double timeFreq)
  : ClpEventHandler()
  , shared_(std::make_shared<SharedState>())
  , msgHandler_(msgHandler)
  , logLevel_(logLevel)
  , iterFreq_(iterFreq)
  , timeFreq_(timeFreq)
{
  shared_->startTime = CoinWallclockTime();
}

ClpEventHandler *CbcRootLpEventHandler::clone() const
{
  // std::shared_ptr copy increments ref-count; all clones share the same
  // SharedState, so headerPrinted/startTime are consistent across
  // Clp-internal model copies (e.g. idiot crash, barrier crossover).
  return new CbcRootLpEventHandler(*this);
}

const char *CbcRootLpEventHandler::algoName(int algo)
{
  switch (algo) {
  case 1: return "primal simplex";
  case -1: return "dual simplex";
  case 2: return "barrier";
  default: return "LP";
  }
}

void CbcRootLpEventHandler::printHeader() const
{
  if (logLevel_ <= 0 || !msgHandler_)
    return;
  FILE *fp = msgHandler_->filePointer();
  if (!fp)
    return;

  const bool u8 = CbcOutput::useUtf8();
  const bool compact = CbcOutput::useCompact();
  const CoinTable tbl = makeLpTable(u8, compact);
  const char *alg = model_ ? algoName(model_->algorithm()) : "LP";

  std::string title = std::string("Root LP relaxation (") + alg + ")";
  fprintf(fp, "\n%s\n\n", CoinTable::phaseStart(title, u8).c_str()); // phase start
  printTableOpen(fp, tbl);
}

void CbcRootLpEventHandler::printRow(int iter, double obj, double primalInf,
  double dualInf, double elapsed) const
{
  if (!msgHandler_)
    return;
  FILE *fp = msgHandler_->filePointer();
  if (!fp)
    return;
  const bool u8 = CbcOutput::useUtf8();
  const char *bar = tableBar(u8, CbcOutput::useCompact());
  const int algo = model_ ? model_->algorithm() : 0;
  const char *phase = nullptr;
  switch (algo) {
  case -1: phase = "Dual";    break;
  case  1: phase = "Primal";  break;
  case  2: phase = "Barrier"; break;
  default: phase = "LP";      break;
  }
  const std::string tbuf = fmtTime(elapsed);
  fprintf(fp, "  %-*s%s%*d%s%*.6e%s%*.4e%s%*.4e%s%*s\n",
    LP_W_PHASE, phase, bar,
    LP_W_ITER, iter, bar,
    LP_W_OBJ, obj, bar,
    LP_W_PINF, (primalInf > 0.0 ? primalInf : 0.0), bar,
    LP_W_DINF, (dualInf > 0.0 ? dualInf : 0.0), bar,
    LP_W_TIME, tbuf.c_str());
  fflush(fp);
}

void CbcRootLpEventHandler::printFinalStatus(int numInts, int numFrac) const
{
  if (logLevel_ <= 0 || !msgHandler_ || !model_)
    return;
  FILE *fp = msgHandler_->filePointer();
  if (!fp)
    return;

  if (!shared_->headerPrinted)
    return;

  const bool u8 = CbcOutput::useUtf8();
  const CoinTable tbl = makeLpTable(u8, CbcOutput::useCompact());
  const double elapsed = CoinWallclockTime() - shared_->startTime;
  const int iters = (model_->numberIterations() > shared_->maxIterSeen)
    ? model_->numberIterations()
    : shared_->maxIterSeen;

  // Always show the last iteration if frequency limiting skipped it.
  if (shared_->lastPrintIter < iters) {
    double now = CoinWallclockTime();
    printRow(iters, model_->objectiveValue(),
      model_->sumPrimalInfeasibilities(),
      model_->sumDualInfeasibilities(), now);
  }

  printTableClose(fp, tbl);

  // Build status string: use "LP " prefix when this is part of a MIP solve
  const char *baseStatus = "Unknown";
  int s = model_->status();
  int ss = model_->secondaryStatus();
  if (s == 0)
    baseStatus = "Optimal";
  else if (s == 1) {
    if (ss == 4)
      baseStatus = "Iteration limit";
    else if (ss == 9)
      baseStatus = "Time limit";
    else
      baseStatus = "Stopped";
  } else if (s == 2)
    baseStatus = "Infeasible";
  else if (s == 3)
    baseStatus = "Unbounded";
  else if (s == 4)
    baseStatus = "Numerical difficulties";
  else if (s == 5)
    baseStatus = "Stopped by event";

  // Prefix "LP " when inside a MIP solve
  std::string statusStr = (numInts > 0 ? std::string("LP ") : std::string()) + baseStatus;

  // Build summary line: "LP Optimal — Frac: N/M (P%)   Obj: X   Iters: N   Time: Xs"
  char summary[512];
  const std::string tStr = fmtTime(elapsed);
  if (numInts > 0 && numFrac >= 0) {
    double pct = (numInts > 0) ? 100.0 * numFrac / numInts : 0.0;
    std::snprintf(summary, sizeof(summary),
      "%s%sFrac: %d/%d (%.1f%%)   Obj: %g   Iters: %d   Time: %ss",
      statusStr.c_str(), CoinTable::dashSep(u8),
      numFrac, numInts, pct, model_->objectiveValue(), iters, tStr.c_str());
  } else {
    std::snprintf(summary, sizeof(summary),
      "%s%sObj: %g   Iters: %d   Time: %ss",
      statusStr.c_str(), CoinTable::dashSep(u8),
      model_->objectiveValue(), iters, tStr.c_str());
  }
  fprintf(fp, "\n%s\n", CoinTable::phaseEnd(summary, u8).c_str());
  fflush(fp);
}

int CbcRootLpEventHandler::event(Event whichEvent)
{
  if (whichEvent != endOfIteration || !model_)
    return -1;

  if (logLevel_ <= 0 || !msgHandler_)
    return -1;

  // Print header once across all clones (shared flag)
  if (!shared_->headerPrinted) {
    printHeader();
    shared_->headerPrinted = true;
    // Print iteration 0 row (absolute wall-clock time)
    double now = CoinWallclockTime();
    printRow(0, model_->objectiveValue(),
      model_->sumPrimalInfeasibilities(),
      model_->sumDualInfeasibilities(), now);
    shared_->lastPrintTime = now;
    shared_->lastPrintIter = 0;
    return -1;
  }

  const int iter = model_->numberIterations();
  const double now = CoinWallclockTime();

  // Track peak iteration count across all sub-solves
  if (iter > shared_->maxIterSeen)
    shared_->maxIterSeen = iter;

  // Detect restart: iteration counter dropped significantly
  // (Clp restarts after perturbation, sprint, or algorithm switch)
  const bool isRestart = (shared_->lastPrintIter > 0 && iter < shared_->lastPrintIter - 50);
  if (isRestart) {
    if (msgHandler_) {
      FILE *fp = msgHandler_->filePointer();
      if (fp) {
        const bool u8 = CbcOutput::useUtf8();
        const CoinTable tbl = makeLpTable(u8, CbcOutput::useCompact());
        const std::string label = u8 ? "\xe2\x94\x80\xe2\x94\x80 restart " : "-- restart ";
        fprintf(fp, "%s\n", tbl.markerRow(label, 2).c_str());
        fflush(fp);
      }
    }
    shared_->lastPrintIter = -1; // force next row to print
    shared_->lastPrintTime = now;
  }

  const bool doIter = (iterFreq_ > 0 && iter - shared_->lastPrintIter >= iterFreq_);
  const bool doTime = (timeFreq_ > 0.0 && now - shared_->lastPrintTime >= timeFreq_);
  const bool doForce = (shared_->lastPrintIter < 0);

  if (doIter || doTime || doForce) {
    printRow(iter, model_->objectiveValue(),
      model_->sumPrimalInfeasibilities(),
      model_->sumDualInfeasibilities(), now);
    shared_->lastPrintTime = now;
    shared_->lastPrintIter = iter;
  }

  return -1; // continue solving
}

// ---------------------------------------------------------------------------
// Number formatting helpers
// ---------------------------------------------------------------------------

/** Format a non-negative number using "natural" notation when possible.
 *  Uses %.5g which gives 5 significant figures and automatically switches
 *  to scientific notation for very large (>=1e5) or very small (<1e-4) values.
 *  Result: "1", "500", "1.5", "0.001", "1.23e+06", etc.
 */
static std::string fmtNum(double x)
{
  if (x == 0.0)
    return "0";
  char buf[32];
  std::snprintf(buf, sizeof(buf), "%.5g", std::fabs(x));
  return buf;
}

/** Right-pad string s to width w with leading spaces. */
static std::string rpad(const std::string &s, int w)
{
  if ((int)s.size() >= w)
    return s;
  return std::string(w - (int)s.size(), ' ') + s;
}

// ---------------------------------------------------------------------------
// Range statistics helper
// ---------------------------------------------------------------------------

struct RangeStat {
  double lo = std::numeric_limits<double>::max();
  double hi = 0.0;
  bool valid = false;

  void update(double v)
  {
    double av = std::fabs(v);
    if (av < 1e-300)
      return;
    valid = true;
    if (av < lo)
      lo = av;
    if (av > hi)
      hi = av;
  }

  double ratio() const
  {
    return (valid && lo > 0.0) ? hi / lo : 1.0;
  }
};

// ---------------------------------------------------------------------------
// send a single line through the handler (prefix already disabled by caller)
// ---------------------------------------------------------------------------

static void sendLine(CoinMessageHandler *handler, const std::string &line)
{
  CbcMessage msgs;
  handler->message(CBC_PROBLEM_SUMMARY, msgs) << line << CoinMessageEol;
}

static void sendBlank(CoinMessageHandler *handler)
{
  // CoinMessageHandler may suppress empty %s messages; use fprintf directly.
  // filePointer() returns stdout by default; custom handlers will still see
  // content-bearing lines via the message callback.
  FILE *fp = handler->filePointer();
  if (fp)
    fprintf(fp, "\n");
}

// ---------------------------------------------------------------------------
// printSolverHeader
// ---------------------------------------------------------------------------

void CbcOutput::printSolverHeader(CoinMessageHandler *handler, int logLevel,
  const std::string &args, const std::string &strategyNote)
{
  if (logLevel <= 0 || !handler)
    return;

  handler->setPrefix(false);
  const bool u8 = useUtf8();
  CbcMessage msgs;

  // Version line: "CBC devel (git:abc1234) — COIN-OR Branch and Cut"
  std::ostringstream oss;
  oss << "CBC ";
  if (std::strcmp(CBC_VERSION, "devel") != 0)
    oss << "v" << CBC_VERSION;
  else
    oss << "devel";

#ifdef CBC_GIT_HASH
  const char *ghash = CBC_GIT_HASH;
  if (ghash && *ghash && std::strcmp(ghash, "unknown") != 0)
    oss << " (git:" << ghash << ")";
#endif

  oss << (u8 ? " \xe2\x80\x94 " : " -- ") << "COIN-OR Branch and Cut";

  handler->message(CBC_HEADER, msgs) << oss.str() << CoinMessageEol;

  // Args line (when provided)
  if (!args.empty()) {
    std::string aline("  args:");
    aline += args;
    handler->message(CBC_HEADER, msgs) << aline << CoinMessageEol;
  }

  // Strategy note on its own line with same indentation
  if (!strategyNote.empty()) {
    std::string sline("  ");
    sline += strategyNote;
    handler->message(CBC_HEADER, msgs) << sline << CoinMessageEol;
  }
}

void CbcOutput::printSolverHeader(CbcModel &model, const std::string &args,
  const std::string &strategyNote)
{
  printSolverHeader(model.messageHandler(),
    model.messageHandler()->logLevel(), args, strategyNote);
}

// ---------------------------------------------------------------------------
// printProblemSummary
// ---------------------------------------------------------------------------

void CbcOutput::printProblemSummary(CoinMessageHandler *handler,
  const OsiSolverInterface &solver, int logLevel,
  const CbcImportHandler *ih)
{
  if (logLevel <= 0 || !handler)
    return;

  handler->setPrefix(false);
  const bool u8 = useUtf8();

  const int nRows = solver.getNumRows();
  const int nCols = solver.getNumCols();

  // Variable type counts
  int nBin = 0, nGenInt = 0, nCts = 0;
  const double *colLb = solver.getColLower();
  const double *colUb = solver.getColUpper();
  for (int j = 0; j < nCols; j++) {
    if (solver.isInteger(j)) {
      if (colLb[j] == 0.0 && colUb[j] == 1.0)
        nBin++;
      else
        nGenInt++;
    } else {
      nCts++;
    }
  }

  // Non-zero count
  const CoinPackedMatrix *mat = solver.getMatrixByRow();
  const int nNz = mat ? mat->getNumElements() : 0;

  // Problem name
  std::string probName;
  solver.getStrParam(OsiProbName, probName);

  // ---- Coefficient ranges ----
  RangeStat rsMatrix, rsCost, rsBounds, rsRhs;

  // Matrix non-zeros
  if (mat) {
    const double *elems = mat->getElements();
    int nel = mat->getNumElements();
    for (int k = 0; k < nel; k++)
      rsMatrix.update(elems[k]);
  }

  // Objective coefficients
  const double *obj = solver.getObjCoefficients();
  if (obj) {
    for (int j = 0; j < nCols; j++)
      rsCost.update(obj[j]);
  }

  // Bounds: non-zero finite lower and upper bounds
  const double inf = solver.getInfinity();
  if (colLb && colUb) {
    for (int j = 0; j < nCols; j++) {
      if (colLb[j] != 0.0 && colLb[j] > -inf && colLb[j] < inf)
        rsBounds.update(colLb[j]);
      if (colUb[j] != 0.0 && colUb[j] > -inf && colUb[j] < inf)
        rsBounds.update(colUb[j]);
    }
  }

  // RHS
  const double *rhs = solver.getRightHandSide();
  const char *sense = solver.getRowSense();
  if (rhs && sense) {
    for (int i = 0; i < nRows; i++) {
      if (sense[i] != 'N')
        rsRhs.update(rhs[i]);
    }
  }

  // ---- Format numbers for the ranges table ----
  // Columns: lo, hi, ratio  (one row per stat)
  struct RangeRow {
    const char *label;
    const RangeStat &rs;
  };
  const RangeRow rows[] = {
    { "Matrix", rsMatrix },
    { "Cost", rsCost },
    { "Bounds", rsBounds },
    { "RHS", rsRhs }
  };

  // Pre-format all numbers so we can determine column widths for alignment
  std::vector<std::string> sLo(4), sHi(4), sRatio(4);
  int wLo = 1, wHi = 1, wRatio = 1;
  for (int r = 0; r < 4; r++) {
    if (rows[r].rs.valid) {
      sLo[r] = fmtNum(rows[r].rs.lo);
      sHi[r] = fmtNum(rows[r].rs.hi);
      sRatio[r] = fmtNum(rows[r].rs.ratio());
    } else {
      sLo[r] = sHi[r] = sRatio[r] = "-";
    }
    wLo = std::max(wLo, (int)sLo[r].size());
    wHi = std::max(wHi, (int)sHi[r].size());
    wRatio = std::max(wRatio, (int)sRatio[r].size());
  }

  // ---- Emit lines ----

  FILE *fp = handler->filePointer();

  // Blank line + phase start
  sendBlank(handler);
  if (fp)
    fprintf(fp, "%s\n\n", CoinTable::phaseStart("Problem loading", u8).c_str());

  // Problem name line
  if (!probName.empty() && probName != "BLANK") {
    sendLine(handler, "  Problem: " + probName);
  }

  // Dimensions line
  {
    std::ostringstream d;
    d << "  |Rows| = " << nRows
      << "   |Cols| = " << nCols
      << "   |NZ| = " << nNz;
    sendLine(handler, d.str());
  }

  // Variables line
  {
    std::ostringstream v;
    v << "  Variables: ";
    const int nInt = nBin + nGenInt;
    if (nInt > 0)
      v << nBin << " binary, " << nGenInt << " integer, " << nCts << " continuous";
    else
      v << nCts << " continuous";
    sendLine(handler, v.str());
  }

  // SOS sets declared in the instance (from MPS SOS section or equivalent)
  {
    const OsiClpSolverInterface *clpSolver
      = dynamic_cast< const OsiClpSolverInterface * >(&solver);
    if (clpSolver) {
      int nSOS = clpSolver->numberSOS();
      if (nSOS > 0) {
        const CoinSet *sets = clpSolver->setInfo();
        int nSOS1 = 0, nSOS2 = 0, mem1 = 0, mem2 = 0;
        for (int i = 0; i < nSOS; i++) {
          int n = sets[i].numberEntries();
          if (sets[i].setType() == 1) { nSOS1++; mem1 += n; }
          else                         { nSOS2++; mem2 += n; }
        }
        if (nSOS1 > 0) {
          std::ostringstream s;
          s << "  SOS type 1: " << nSOS1 << " sets, " << mem1 << " members"
            << " (avg " << std::fixed << std::setprecision(1)
            << (double)mem1 / nSOS1 << " members/set)";
          sendLine(handler, s.str());
        }
        if (nSOS2 > 0) {
          std::ostringstream s;
          s << "  SOS type 2: " << nSOS2 << " sets, " << mem2 << " members"
            << " (avg " << std::fixed << std::setprecision(1)
            << (double)mem2 / nSOS2 << " members/set)";
          sendLine(handler, s.str());
        }
      }
    }
  }

  // Coefficient ranges section (only when matrix has entries)
  if (rsMatrix.valid) {
    sendBlank(handler);
    sendLine(handler, "  Coefficient ranges:");

    // UTF-8 symbols
    const char *sym_in = u8 ? "\xe2\x88\x88" : "in";      // ∈
    const char *sym_kap = u8 ? "\xce\xba \xe2\x89\x88" : "ratio ="; // κ ≈

    for (int r = 0; r < 4; r++) {
      std::ostringstream row;
      // label: left-aligned in 8 chars
      char lbuf[16];
      std::snprintf(lbuf, sizeof(lbuf), "    %-8s", rows[r].label);
      row << lbuf;
      row << sym_in << " ["
          << rpad(sLo[r], wLo) << ", "
          << rpad(sHi[r], wHi) << "]"
          << "   " << sym_kap << " "
          << rpad(sRatio[r], wRatio);
      sendLine(handler, row.str());
    }
  }

  // Phase-end line — replaces the old trailing blank line.
  // One blank line before the next phase comes from its own "\n▶ ..." prefix.
  if (fp) {
    // Print any parse warnings before closing the section.
    if (ih && ih->totalErrors() > 0)
      printImportErrors(fp, *ih);
    std::ostringstream summary;
    summary << "Problem loaded";
    if (ih && ih->totalErrors() > 0)
      summary << " with " << ih->totalErrors() << " parse warning(s)";
    summary << CoinTable::dashSep(u8) << nRows << " rows"
            << ", " << nCols << " cols"
            << ", " << nNz << " NZ";
    if (nBin + nGenInt > 0)
      summary << ", " << (nBin + nGenInt) << " integer";
    fprintf(fp, "\n%s\n", CoinTable::phaseEnd(summary.str(), u8).c_str());
    fflush(fp);
  }
}

void CbcOutput::printProblemSummary(CbcModel &model,
  const OsiSolverInterface &solver, const CbcImportHandler *ih)
{
  printProblemSummary(model.messageHandler(), solver,
    model.messageHandler()->logLevel(), ih);
}

void CbcOutput::printImportErrors(FILE *fp, const CbcImportHandler &ih)
{
  if (ih.totalErrors() == 0)
    return;
  if (!fp)
    fp = stdout;
  const bool u8 = useUtf8();
  const char *warn = u8 ? "⚠" : "!";
  const int total = ih.totalErrors();
  const auto &msgs = ih.errors();
  const int shown = (int)msgs.size();
  fprintf(fp, "\n  %s %d parse warning(s) in input file:\n", warn, total);
  for (const auto &m : msgs)
    fprintf(fp, "    %s\n", m.c_str());
  if (total > shown)
    fprintf(fp, "    ... and %d more\n", total - shown);
  fflush(fp);
}

#ifdef CBC_HAS_NAUTY
// ===========================================================================
// CbcOutput::printNautySection — formatted symmetry detection output (static)
// ===========================================================================

void CbcOutput::printNautySection(FILE *fp, bool utf8,
  int numUsefulOrbits, int numUsefulObjects,
  int totalOrbits, int numGenerators, double groupSize,
  double nautyTime, int errorCode)
{
  if (!fp) fp = stdout;

  fprintf(fp, "\n%s\n\n", CoinTable::phaseStart("Symmetry detection (nauty)", utf8).c_str());

  char endSummary[256];
  const std::string nStr = fmtTime(nautyTime);
  if (errorCode != 0) {
    fprintf(fp, "  Error (code %d)   Time: %ss\n", errorCode, nStr.c_str());
    std::snprintf(endSummary, sizeof(endSummary),
      "Symmetry detection complete%sError (code %d)   Time: %ss",
      CoinTable::dashSep(utf8), errorCode, nStr.c_str());
  } else if (numUsefulOrbits > 0) {
    char groupStr[32];
    if (groupSize >= 1e15)
      std::snprintf(groupStr, sizeof(groupStr), "%s", utf8 ? "\xe2\x88\x9e" : "inf");
    else if (groupSize >= 1e9)
      std::snprintf(groupStr, sizeof(groupStr), "%.3g", groupSize);
    else
      std::snprintf(groupStr, sizeof(groupStr), "%.6g", groupSize);
    fprintf(fp, "  Orbits: %d (%d useful, %d vars)   Generators: %d   Group size: %s   Time: %ss\n",
      totalOrbits, numUsefulOrbits, numUsefulObjects, numGenerators, groupStr, nStr.c_str());
    std::snprintf(endSummary, sizeof(endSummary),
      "Symmetry detection complete%s%d orbits (%d useful)   Time: %ss",
      CoinTable::dashSep(utf8), totalOrbits, numUsefulOrbits, nStr.c_str());
  } else {
    fprintf(fp, "  No useful orbits found   Time: %ss\n", nStr.c_str());
    std::snprintf(endSummary, sizeof(endSummary),
      "Symmetry detection complete%sno useful orbits   Time: %ss",
      CoinTable::dashSep(utf8), nStr.c_str());
  }
  fprintf(fp, "\n%s\n", CoinTable::phaseEnd(endSummary, utf8).c_str());
  fflush(fp);
}

// ===========================================================================
// CbcNautyHandler — intercepts Nauty messages during branchAndBound()
// ===========================================================================

CbcNautyHandler::CbcNautyHandler(FILE *fp, bool utf8, int logLevel)
  : fp_(fp ? fp : stdout)
  , utf8_(utf8)
  , compact_(CbcOutput::useCompact())
{
  CoinMessageHandler::setLogLevel(logLevel);
  setFilePointer(fp_);
}

CbcNautyHandler::~CbcNautyHandler()
{
  delete lpSilentHandler_;
}

CoinMessageHandler *CbcNautyHandler::getLpSilentHandler()
{
  if (!lpSilentHandler_) {
    lpSilentHandler_ = new CoinMessageHandler();
    lpSilentHandler_->setLogLevel(0);
  }
  return lpSilentHandler_;
}

int CbcNautyHandler::print()
{
  const char *buf = messageBuffer();
  const int ext   = currentMessage().externalNumber();

  if (currentSource() == "Cbc") {
    // Suppress internal priority/branching setup CBC_GENERAL messages.
    if (ext == 45 && (std::strstr(buf, "have cost of") || std::strstr(buf, "branch on satisfied")
          || std::strstr(buf, "have costs") || std::strstr(buf, "have cost")))
      return 0;

    // CBC_INFEAS (ext=6): "The LP relaxation is infeasible or too expensive"
    // Show as a structured note; the calling code will print its own result banner.
    if (ext == 6) {
      FILE *fp = filePointer();
      if (fp) {
        fprintf(fp, "\n  ✘ LP relaxation infeasible or too expensive\n\n");
        fflush(fp);
      }
      return 0;
    }

    // CBC_GENERAL_WARNING (ext=3009): strip the raw prefix and print as indented warning.
    if (ext == 3009) {
      FILE *fp = filePointer();
      if (fp) {
        // Skip past any "CbcNNNNW " prefix
        const char *p = buf;
        if (std::strncmp(p, "Cbc", 3) == 0) {
          while (*p && *p != ' ') ++p;
          if (*p == ' ') ++p;
          // Also skip a possible "ClpNNNNW " prefix that CBC may embed
          if (std::strncmp(p, "Clp", 3) == 0) {
            while (*p && *p != ' ') ++p;
            if (*p == ' ') ++p;
          }
        }
        fprintf(fp, "  ⚠ %s\n", p);
        fflush(fp);
      }
      return 0;
    }

    // Suppress "Full problem N rows M columns, reduced to..." (CBC_FPUMP1=38)
    // This is a heuristic sub-model reduction notice not needed in tabular output.
    if (ext == 38)
      return 0;

    // ── B&B tree progress ────────────────────────────────────────────────
    if (bnbOut_) {
      // Any first B&B message must flush pending cut-gen generator table first.
      auto flushCutGen = [this]() {
        if (cutGenOut_ && cutGenOut_->hasPendingGenerators())
          cutGenOut_->close();
      };

      // CBC_GAP (ext=11): "Exiting as integer gap of X less than Y or Z%"
      // Suppress — the ✔ Optimal banner from ext=1 already conveys this.
      if (ext == 11)
        return 0;

      // ext=37 (CBC_STATUS2): periodic progress with depth
      // Format: "NNN nodes, N on tree, best G - possible G depth N unsat N value G its N (G seconds)"
      // Note: messageBuffer() includes the "Cbc0037I " prefix; use %*s to skip it.
      if (ext == 37) {
        flushCutGen();
        long nodes, iters; int onTree, depth, unsat; double bestSol, bestBound, value, elapsed;
        if (std::sscanf(buf, "%*s %ld nodes, %d on tree, best %lf - possible %lf depth %d unsat %d value %lf its %ld (%lf",
              &nodes, &onTree, &bestSol, &bestBound, &depth, &unsat, &value, &iters, &elapsed) == 9)
          bnbOut_->onProgress(nodes, onTree, depth, bestSol, bestBound, elapsed);
        return 0;
      }
      // ext=4: B&B LP found a solution (no heuristic name in message)
      if (ext == 4) {
        // Root-node heuristics fire before/during cut gen. Only treat as a B&B
        // event once cut generation has fully completed.
        if (cutGenOut_ && !cutGenOut_->hasClosed()) return 0;
        const char *p = std::strstr(buf, "Integer solution of ");
        if (p) {
          double obj, elapsed; long nodes;
          if (std::sscanf(p, "Integer solution of %lf found after %*ld iterations and %ld nodes (%lf",
                &obj, &nodes, &elapsed) == 3)
            bnbOut_->onBnBIncumbent(obj, nodes, elapsed);
        }
        return 0;
      }
      // ext=12: heuristic found a solution (suppress during FPump phase; show as ★ during B&B)
      // ext=16: strong branching found a solution (show as ★ during B&B)
      if (ext == 12 || ext == 16) {
        if (ext == 12 && fpumpOut_ && fpumpOut_->isInPhase()) return 0;
        IncumbentMsg im;
        if (cutGenOut_ && !cutGenOut_->hasClosed()) {
          // Heuristic found before/during root cut gen — queue for the B&B table.
          // Will be emitted as ★ rows once the first progress message arrives.
          if (parseIncumbentMsg(buf, im))
            bnbOut_->queuePreProgressIncumbent(im.obj, im.method.c_str(), im.nodes, im.elapsed);
        } else {
          // Cut gen is done (or was never started): emit directly as ★ row.
          routeIncumbentMessage(buf, ext);
        }
        return 0;
      }
      // ext=1: optimal (CBC_END_GOOD)
      if (ext == 1) {
        const char *p = std::strstr(buf, "Search completed");
        if (p) {
          double bestSol, elapsed; long iters, nodes;
          if (restartMode_) {
            // This is the sub-model's completion after a B&B restart.
            // Close cut-gen part 2 if still open; do NOT trigger onComplete().
            if (cutGenOut_) cutGenOut_->close();
            restartMode_ = false;
            return 0;
          }
          if (std::sscanf(p, "Search completed - best objective %lf , took %ld iterations and %ld nodes (%lf",
                &bestSol, &iters, &nodes, &elapsed) == 4 ||
              std::sscanf(p, "Search completed - best objective %lf, took %ld iterations and %ld nodes (%lf",
                &bestSol, &iters, &nodes, &elapsed) == 4)
            bnbOut_->onComplete(true, bestSol, bestSol, iters, nodes, elapsed);
        }
        return 0;
      }
      // ext=5: partial search (CBC_END)
      if (ext == 5) {
        const char *p = std::strstr(buf, "Partial search");
        if (p) {
          double bestSol, bestBound, elapsed; long iters, nodes;
          if (restartMode_) {
            // Sub-model ended without finding optimal — still not full B&B done.
            if (cutGenOut_) cutGenOut_->close();
            restartMode_ = false;
            return 0;
          }
          if (std::sscanf(p,
                "Partial search - best objective %lf (best possible %lf), took %ld iterations and %ld nodes (%lf",
                &bestSol, &bestBound, &iters, &nodes, &elapsed) == 5)
            bnbOut_->onComplete(false, bestSol, bestBound, iters, nodes, elapsed);
          else
            bnbOut_->onComplete(false, 1e50, 1e50, 0, 0, 0.0);
        }
        return 0;
      }
      // Stopping reasons
      if (ext == 3)  { bnbOut_->onStopReason("node limit");      return 0; }
      if (ext == 19) { bnbOut_->onStopReason("solution limit");  return 0; }
      if (ext == 20) { bnbOut_->onStopReason("time limit");      return 0; }
      if (ext == 50) { bnbOut_->onStopReason("iteration limit"); return 0; }
      // ext=32: strong branching stats (CBC_STRONG_STATS)
      if (ext == 32) {
        const char *p = std::strstr(buf, "Strong branching done ");
        if (p) {
          long calls, iters, fathomed; int fixed;
          if (std::sscanf(p,
                "Strong branching done %ld times (%ld iterations), fathomed %ld nodes and fixed %d variables",
                &calls, &iters, &fathomed, &fixed) == 4)
            bnbOut_->onStrongStats(calls, iters, fathomed, fixed);
        }
        return 0;
      }
      // ext=35: depth/DJ stats (CBC_OTHER_STATS)
      if (ext == 35) {
        const char *p = std::strstr(buf, "Maximum depth ");
        if (p) {
          int maxDepth; double djFixed;
          if (std::sscanf(p, "Maximum depth %d, %lf variables fixed on reduced cost",
                &maxDepth, &djFixed) == 2)
            bnbOut_->onOtherStats(maxDepth, djFixed);
        }
        return 0;
      }
      // ext=41: depth/fathoming stats (CBC_OTHER_STATS2)
      if (ext == 41) {
        const char *p = std::strstr(buf, "Maximum depth ");
        if (p) {
          int maxDepth; double djFixed; long fTimes, fNodes, fIters;
          if (std::sscanf(p,
                "Maximum depth %d, %lf variables fixed on reduced cost (complete fathoming %ld times, %ld nodes taking %ld iterations)",
                &maxDepth, &djFixed, &fTimes, &fNodes, &fIters) == 5)
            bnbOut_->onOtherStats2(maxDepth, djFixed, fTimes, fNodes, fIters);
          else {
            if (std::sscanf(p, "Maximum depth %d, %lf variables fixed on reduced cost",
                  &maxDepth, &djFixed) == 2)
              bnbOut_->onOtherStats(maxDepth, djFixed);
          }
        }
        return 0;
      }
      // ext=45: orbital branching stats
      if (ext == 45 && std::strstr(buf, "Orbital branching succeeded")) {
        const char *p = std::strstr(buf, "Orbital branching succeeded");
        if (p) {
          int successes, fixed; double avgExtra, avgFixed;
          if (std::sscanf(p, "Orbital branching succeeded %d times - average extra %lf , fixing (%d, %lf)",
                &successes, &avgExtra, &fixed, &avgFixed) == 4 ||
              std::sscanf(p, "Orbital branching succeeded %d times - average extra %lf, fixing (%d, %lf)",
                &successes, &avgExtra, &fixed, &avgFixed) == 4)
            bnbOut_->onOrbitalStats(successes, avgExtra, fixed, avgFixed);
        }
        return 0;
      }
      // ext=44: CBC_RESTART — reduced cost fixing, restarting root node
      // "Reduced cost fixing - N rows, M columns - restarting search"
      if (ext == 44) {
        int rows = 0, cols = 0;
        const char *p = std::strstr(buf, "Reduced cost fixing - ");
        if (p)
          std::sscanf(p, "Reduced cost fixing - %d rows, %d columns", &rows, &cols);
        bnbOut_->onRestartNote(rows, cols);
        return 0;
      }
    } // end bnbOut_ block

    // ext=12 without bnbOut_ but with fpumpOut_: suppress between phases
    if (ext == 12 && fpumpOut_) return 0;
  }

  // ── Cut generation at root node ─────────────────────────────────────────
  if (cutGenOut_ && currentSource() == "Cbc") {

    // ext=51: "Starting cut generation at root node"
    if (ext == 51) {
      cutGenOut_->onStart();
      return 0;
    }
    // ext=46: "Root node pass %d, %d rows, %d tight cuts, %d frac, %g suminf - objective %g (%.2f seconds)"
    if (ext == 46) {
      const char *p = std::strstr(buf, "Root node pass ");
      if (p) {
        int pass, rows, tight, frac;
        double suminf, obj, t;
        if (std::sscanf(p,
              "Root node pass %d, %d rows, %d tight cuts, %d frac, %lf suminf - objective %lf (%lf",
              &pass, &rows, &tight, &frac, &suminf, &obj, &t) == 7)
          cutGenOut_->onPass(pass, rows, tight, frac, suminf, obj, t);
      }
      return 0;
    }
    // ext=31: "N added rows had average density of X" — suppress (not per-pass detail)
    if (ext == 31)
      return 0;
    // ext=13: "At root node, %d cuts changed objective from %g to %g in %d passes"
    if (ext == 13) {
      const char *p = std::strstr(buf, "At root node, ");
      if (p) {
        int ncuts, passes;
        double fromObj, toObj;
        if (std::sscanf(p,
              "At root node, %d cuts changed objective from %lf to %lf in %d passes",
              &ncuts, &fromObj, &toObj, &passes) == 4)
          cutGenOut_->onSummary(ncuts, fromObj, toObj, passes);
      }
      return 0;
    }
    // ext=14: "Cut generator %d (%s) - %d row cuts average %.1f elements, ..."
    if (ext == 14) {
      const char *p = std::strstr(buf, "Cut generator ");
      if (p) {
        CbcCutGenOutput::GenInfo g;
        char name[64] = "";
        if (std::sscanf(p,
              "Cut generator %d (%63[^)]) - %d row cuts average %lf elements, %d column cuts",
              &g.idx, name, &g.rowCuts, &g.avgDensity, &g.colCuts) >= 5) {
          g.name = name;
          // Optional timing: "in %.3f seconds"
          const char *freqPtr = std::strstr(p, "new frequency is ");
          const char *inPtr   = std::strstr(p, " in ");
          g.hasTime = (inPtr && freqPtr && inPtr < freqPtr);
          g.time = 0.0;
          if (g.hasTime)
            std::sscanf(inPtr, " in %lf seconds", &g.time);
          g.newFreq = -100;
          if (freqPtr)
            std::sscanf(freqPtr, "new frequency is %d", &g.newFreq);
          cutGenOut_->onGenerator(g);
        }
      }
      return 0;
    }
  }

  // ── Nauty messages ───────────────────────────────────────────────────────
  // We only intercept CBC_GENERAL messages from "Cbc" that contain "Nauty"
  if (currentSource() == "Cbc" && std::strstr(buf, "Nauty")) {
    sectionStarted_ = true;

    // Find the actual "Nauty " start (skip any "CbcNNNNI " prefix)
    const char *p = std::strstr(buf, "Nauty");
    if (!p) p = buf;

    // "Nauty: N orbits (U useful covering V variables), G generators, group size: S - sparse size X - took T seconds"
    if (std::sscanf(p,
          "Nauty: %d orbits (%d useful covering %d variables), %d generators, group size: %lf",
          &totalOrbits_, &usefulOrbits_, &usefulVars_, &numGens_, &groupSize_) >= 5) {
      const char *tookPtr = std::strstr(p, "took ");
      if (tookPtr) std::sscanf(tookPtr, "took %lf", &nautyTime_);
      printSection();
      return 0;
    }

    // "Nauty did not find any useful orbits in time T"
    if (std::strstr(p, "did not find any useful orbits")) {
      const char *tp = std::strstr(p, "in time ");
      if (tp) std::sscanf(tp, "in time %lf", &nautyTime_);
      printSection();
      return 0;
    }

    // "Nauty failed with error code N (T seconds)"
    if (std::strstr(p, "failed with error code")) {
      hasError_ = true;
      std::sscanf(p, "Nauty failed with error code %d (%lf", &errorCode_, &nautyTime_);
      printSection();
      return 0;
    }

    // "Nauty sparseSpace ..." or "Nauty - initial level N - will probably take too long"
    // or any other Nauty informational prefix — suppress silently
    if (std::strstr(p, "sparseSpace") || std::strstr(p, "initial level"))
      return 0;

    // Unknown Nauty message — suppress to avoid raw prefix output
    return 0;
  } // end Nauty block

  // If cut gen generators are accumulated and a non-cut-gen message has arrived,
  // flush the generator table now (before B&B tree output begins).
  if (cutGenOut_ && cutGenOut_->hasPendingGenerators())
    cutGenOut_->close();

  // CLP_GENERAL_WARNING (ext=3006) from Clp: accuracy/infeasibility warnings during B&B.
  // Strip the raw "ClpNNNNW" prefix and print as an indented warning.
  if (currentSource() == "Clp" && ext == 3006) {
    FILE *fp = filePointer();
    if (fp) {
      const char *p = buf;
      if (std::strncmp(p, "Clp", 3) == 0) {
        while (*p && *p != ' ') ++p;
        if (*p == ' ') ++p;
      }
      fprintf(fp, "  ⚠ %s\n", p);
      fflush(fp);
    }
    return 0;
  }

  return CoinMessageHandler::print();
}

void CbcNautyHandler::printSection()
{
  CbcOutput::printNautySection(fp_, utf8_,
    usefulOrbits_, usefulVars_, totalOrbits_,
    numGens_, groupSize_, nautyTime_,
    hasError_ ? errorCode_ : 0);
}

void CbcNautyHandler::beginRestartMode()
{
  // Reset cut-gen output for the part-2 run and flag that the next
  // ext==1 / ext==5 from a sub-model must not close the B&B section.
  if (cutGenOut_) cutGenOut_->resetForRestart();
  restartMode_ = true;
}

void CbcNautyHandler::routeIncumbentMessage(const char *buf, int /*ext*/)
{
  if (!bnbOut_) return;
  IncumbentMsg im;
  if (parseIncumbentMsg(buf, im))
    bnbOut_->onHeurIncumbent(im.obj, im.method.c_str(), im.nodes, im.elapsed);
}

#endif /* CBC_HAS_NAUTY */

// ===========================================================================
// CbcFPumpOutput — tabular output for the feasibility pump heuristic
// ===========================================================================

static const int FP_W_ROUND  = 5;
static const int FP_W_PASS   = 6;
static const int FP_W_FRAC   = 12;
static const int FP_W_SUMINF = 12;
static const int FP_W_STATUS = 14;
static const int FP_W_BEST   = 14;
static const int FP_W_TIME   = 8;

static CoinTable makeFpTable(bool utf8, bool compact)
{
  return CoinTable({
    { "Round",      FP_W_ROUND  },
    { "Pass",       FP_W_PASS   },
    { "Fractional", FP_W_FRAC   },
    { "Suminf",     FP_W_SUMINF },
    { "Status",     FP_W_STATUS, /*leftAlign=*/true },
    { "BestSol",    FP_W_BEST   },
    { "Time(s)",       FP_W_TIME   },
  }, utf8, /*indent=*/2, compact);
}

CbcFPumpOutput::CbcFPumpOutput(FILE *fp, bool utf8, int logLevel,
  int passFreq, double timeFreq)
  : fp_(fp)
  , utf8_(utf8)
  , compact_(CbcOutput::useCompact())
  , logLevel_(logLevel)
  , passFreq_(passFreq)
  , timeFreq_(timeFreq)
{
}

void CbcFPumpOutput::onStart(int numFrac, double suminf)
{
  if (!isActive()) return;
  startTime_ = CoinWallclockTime();
  lastPrintTime_ = startTime_;
  inPhase_ = true;
  prevSuminf_ = suminf;
  prevNumFrac_ = numFrac;
  // Print header now with the true initial LP state, so the title correctly
  // reflects the number of fractional variables before any rounding takes place.
  printTableHeader(numFrac, suminf);
}

void CbcFPumpOutput::noteRowSolution(const std::string &type, double obj)
{
  if (!isActive()) return;
  pendingSolution_ = true;
  pendingType_ = type;
  pendingObj_ = obj;
  if (obj < bestSol_)
    bestSol_ = obj;
}

bool CbcFPumpOutput::shouldPrint(int pass, double elapsed) const
{
  if (passFreq_ <= 0 && timeFreq_ <= 0.0)
    return true;
  if (printedRows_ < 10)   // always show the first 10 rows
    return true;
  if (passFreq_ > 0 && (pass - lastPrintPass_) >= passFreq_)
    return true;
  if (timeFreq_ > 0.0 && (elapsed - lastPrintTime_) >= timeFreq_)
    return true;
  return false;
}

void CbcFPumpOutput::printTableHeader(int numFrac, double suminf)
{
  if (headerPrinted_) return;
  headerPrinted_ = true;

  char title[80];
  std::snprintf(title, sizeof(title),
    "Feasibility pump (%d fractional, suminf %.4g)", numFrac, suminf);
  const CoinTable tbl = makeFpTable(utf8_, compact_);
  fprintf(fp_, "\n%s\n\n", CoinTable::phaseStart(title, utf8_).c_str());
  printTableOpen(fp_, tbl);
}

void CbcFPumpOutput::printTableEnd()
{
  if (!headerPrinted_ || tableClosed_) return;
  tableClosed_ = true;
  printTableClose(fp_, makeFpTable(utf8_, compact_));
}

void CbcFPumpOutput::printRow(bool hasSolution, int round, int pass,
  const std::string &fracStr, const std::string &suminfStr,
  const std::string &status, double bestSol, double elapsed)
{
  const char *bar = tableBar(utf8_, compact_);
  const std::string timeStr = fmtTime(elapsed);
  char bestStr[24] = "";
  if (bestSol < 1e30)
    std::snprintf(bestStr, sizeof(bestStr), "%.6g", bestSol);
  const char *pfx = hasSolution
    ? (utf8_ ? " \xe2\x98\x85" : " *")  // ★ or *
    : "  ";

  int fracVis = utf8VisLen(fracStr);
  int sumVis  = utf8VisLen(suminfStr);
  int fracPad = FP_W_FRAC   - fracVis;
  int sumPad  = FP_W_SUMINF - sumVis;
  if (fracPad < 0) fracPad = 0;
  if (sumPad  < 0) sumPad  = 0;

  fprintf(fp_, "%s%*d%s%*d%s%*s%s%*s%s%-*s%s%*s%s%*s\n",
    pfx,
    FP_W_ROUND, round, bar,
    FP_W_PASS,  pass,  bar,
    fracPad  + (int)fracStr.size(),  fracStr.c_str(),  bar,
    sumPad   + (int)suminfStr.size(), suminfStr.c_str(), bar,
    FP_W_STATUS, status.c_str(), bar,
    FP_W_BEST, bestStr, bar,
    FP_W_TIME, timeStr.c_str());
  fflush(fp_);
}

void CbcFPumpOutput::onPass(int round, int pass, int numFrac, double suminf,
  double obj, int /*iters*/, bool isOscillating)
{
  if (!isActive()) return;

  lastRound_       = round;
  lastPass_        = pass;
  lastNumFrac_val_ = numFrac;
  lastSuminf_val_  = suminf;

  double elapsed = CoinWallclockTime();
  bool hasPendingSol = pendingSolution_;

  if (!headerPrinted_)
    printTableHeader(numFrac, suminf);

  if (!hasPendingSol && !shouldPrint(pass, elapsed)) {
    lastPrinted_ = false;
    return;
  }

  // Determine status string
  std::string status;
  if (hasPendingSol) {
    if (pendingType_ == "rounding")
      status = "rounding solution";
    else if (pendingType_ == "minibab")
      status = "b&b solution";
    else if (pendingType_ == "cont")
      status = "cont solution";
    else
      status = "solution";
  } else if (isOscillating) {
    status = "oscillating";
  } else if (suminf < prevSuminf_ - 1e-6 || numFrac < prevNumFrac_) {
    status = "converging";
  } else {
    status = "not converging";
  }

  prevSuminf_  = suminf;
  prevNumFrac_ = numFrac;

  // Format fractional and suminf columns
  char fracBuf[24], sumBuf[24];
  if (hasPendingSol && pendingType_ == "rounding") {
    std::snprintf(fracBuf, sizeof(fracBuf), "%s", utf8_ ? "\xe2\x80\x94" : "---");
    std::snprintf(sumBuf,  sizeof(sumBuf),  "rounding");
  } else {
    std::snprintf(fracBuf, sizeof(fracBuf), "%d",    numFrac);
    std::snprintf(sumBuf,  sizeof(sumBuf),  "%.4f",  suminf);
  }

  double bestSolToShow = hasPendingSol ? pendingObj_ : 1e30;
  printRow(hasPendingSol, round, pass, fracBuf, sumBuf, status, bestSolToShow, elapsed);

  lastPrintTime_ = elapsed;
  lastPrintPass_ = pass;
  lastPrinted_   = true;
  printedRows_++;

  pendingSolution_ = false;
  pendingType_.clear();
  pendingObj_ = 1e30;
}

void CbcFPumpOutput::onRoundingImproved(double newObj)
{
  if (!isActive()) return;
  if (newObj < bestSol_)
    bestSol_ = newObj;
}

void CbcFPumpOutput::onNoSolutionInRetry()
{
  if (!isActive()) return;
  if (!headerPrinted_ || lastPass_ == 0) return;
  if (!lastPrinted_) {
    double elapsed = CoinWallclockTime();
    char fracBuf[24], sumBuf[24];
    std::snprintf(fracBuf, sizeof(fracBuf), "%d",   lastNumFrac_val_);
    std::snprintf(sumBuf,  sizeof(sumBuf),  "%.4f", lastSuminf_val_);
    printRow(false, lastRound_, lastPass_, fracBuf, sumBuf, "no solution", 1e30, elapsed);
    lastPrinted_ = true;
    printedRows_++;
  }
}

void CbcFPumpOutput::onRetry(int /*newRound*/, double /*newCutoff*/)
{
  if (!isActive() || !headerPrinted_) return;
  lastPrintTime_ = CoinWallclockTime();
  lastPrintPass_ = 0;
  lastPrinted_   = false;
  prevSuminf_    = 1e30;
  prevNumFrac_   = INT_MAX;
}

void CbcFPumpOutput::onEnd(double bestSol, double /*elapsed_cpu*/, int rounds, int totalPasses)
{
  if (!isActive()) return;
  // Always show the last pass if it wasn't printed (frequency limiting may have skipped it).
  if (headerPrinted_ && !lastPrinted_ && lastPass_ > 0) {
    double elapsed = CoinWallclockTime();
    char fracBuf[24], sumBuf[24];
    std::snprintf(fracBuf, sizeof(fracBuf), "%d",    lastNumFrac_val_);
    std::snprintf(sumBuf,  sizeof(sumBuf),  "%.4f",  lastSuminf_val_);
    printRow(false, lastRound_, lastPass_, fracBuf, sumBuf, "", 1e30, elapsed);
    lastPrinted_ = true;
  }
  if (headerPrinted_)
    printTableEnd();
  // Use wall-clock duration for the summary line (not the CPU time passed from the heuristic).
  const double duration = CoinWallclockTime() - startTime_;
  const std::string dStr = fmtTime(duration);
  char detail[256];
  if (bestSol < 1e30) {
    std::snprintf(detail, sizeof(detail),
      "Feasibility pump%sbest %.6g in %ss (%d round%s, %d pass%s)",
      CoinTable::dashSep(utf8_),
      bestSol, dStr.c_str(),
      rounds, rounds == 1 ? "" : "s",
      totalPasses, totalPasses == 1 ? "" : "es");
  } else {
    std::snprintf(detail, sizeof(detail),
      "Feasibility pump%sno solution in %ss (%d round%s, %d pass%s)",
      CoinTable::dashSep(utf8_),
      dStr.c_str(),
      rounds, rounds == 1 ? "" : "s",
      totalPasses, totalPasses == 1 ? "" : "es");
  }
  fprintf(fp_, "\n%s\n", CoinTable::phaseEnd(detail, utf8_).c_str());
  fflush(fp_);
  inPhase_ = false;
}

// ===========================================================================
// CbcCutGenOutput — root-node cut generation tabular output
// ===========================================================================

// Progress table column widths
static const int CG_W_PASS   = 4;
static const int CG_W_ROWS   = 8;
static const int CG_W_TIGHT  = 8;
static const int CG_W_FRAC   = 6;
static const int CG_W_SUMINF = 10;
static const int CG_W_OBJ    = 16;
static const int CG_W_TIME   = 8;

static CoinTable makeCgProgTable(bool utf8, bool compact)
{
  return CoinTable({
    { "Pass",      CG_W_PASS   },
    { "Rows",      CG_W_ROWS   },
    { "Tight",     CG_W_TIGHT  },
    { "Frac",      CG_W_FRAC   },
    { "Suminf",    CG_W_SUMINF },
    { "Objective", CG_W_OBJ    },
    { "Time(s)",      CG_W_TIME   },
  }, utf8, /*indent=*/2, compact);
}

// Generator summary table column widths
static const int CG_W_GENNAME  = 22;
static const int CG_W_ROWCUTS  = 8;
static const int CG_W_DENSITY  = 11;
static const int CG_W_COLCUTS  = 8;
static const int CG_W_GENTIME  = 9;
static const int CG_W_NEXTRUN  = 12;

static CoinTable makeCgGenTable(bool utf8, bool compact)
{
  return CoinTable({
    { "Generator",   CG_W_GENNAME, /*leftAlign=*/true },
    { "Row cuts",    CG_W_ROWCUTS },
    { "Avg density", CG_W_DENSITY },
    { "Col cuts",    CG_W_COLCUTS },
    { "Time(s)",        CG_W_GENTIME },
    { "Next run",    CG_W_NEXTRUN, /*leftAlign=*/true },
  }, utf8, /*indent=*/2, compact);
}

CbcCutGenOutput::CbcCutGenOutput(FILE *fp, bool utf8, int logLevel)
  : fp_(fp), utf8_(utf8), compact_(CbcOutput::useCompact()), logLevel_(logLevel)
{
}

void CbcCutGenOutput::onStart()
{
  if (state_ != State::Idle) return;
  state_ = State::Started;
  fprintf(fp_, "\n%s\n\n", CoinTable::phaseStart(title_.c_str(), utf8_).c_str());
  fflush(fp_);
}

void CbcCutGenOutput::resetForRestart(const char *title)
{
  state_             = State::Idle;
  progHeaderPrinted_ = false;
  progTableClosed_   = false;
  genTablePrinted_   = false;
  haveSummary_       = false;
  genInfos_.clear();
  sumNcuts_   = 0;
  sumFromObj_ = 0.0;
  sumToObj_   = 0.0;
  sumPasses_  = 0;
  title_ = title ? title : "Cut generation (root node, part 2)";
}

void CbcCutGenOutput::onPass(int pass, int rows, int tight, int frac, double suminf, double obj, double t)
{
  if (state_ != State::Started) return;

  if (!progHeaderPrinted_) {
    progHeaderPrinted_ = true;
    printTableOpen(fp_, makeCgProgTable(utf8_, compact_));
  }

  const char *bar = tableBar(utf8_, compact_);
  char objBuf[24], suminfBuf[16];
  std::snprintf(objBuf,    sizeof(objBuf),    "%.6g", obj);
  std::snprintf(suminfBuf, sizeof(suminfBuf), "%.4g", suminf);
  const std::string timeBuf = fmtTime(CoinWallclockTime());
  fprintf(fp_, "  %*d%s%*d%s%*d%s%*d%s%*s%s%*s%s%*s\n",
    CG_W_PASS,   pass,             bar,
    CG_W_ROWS,   rows,             bar,
    CG_W_TIGHT,  tight,            bar,
    CG_W_FRAC,   frac,             bar,
    CG_W_SUMINF, suminfBuf,        bar,
    CG_W_OBJ,    objBuf,           bar,
    CG_W_TIME,   timeBuf.c_str());
  fflush(fp_);
}

void CbcCutGenOutput::printProgressEnd()
{
  if (!progHeaderPrinted_ || progTableClosed_) return;
  progTableClosed_ = true;
  printTableClose(fp_, makeCgProgTable(utf8_, compact_));
}

void CbcCutGenOutput::onSummary(int ncuts, double fromObj, double toObj, int passes)
{
  if (state_ != State::Started) return;
  sumNcuts_   = ncuts;
  sumFromObj_ = fromObj;
  sumToObj_   = toObj;
  sumPasses_  = passes;
  haveSummary_ = true;
  printProgressEnd();
}

void CbcCutGenOutput::onGenerator(const GenInfo &g)
{
  if (state_ == State::Idle) return;
  genInfos_.push_back(g);
}

static std::string nextRunStr(int freq)
{
  if (freq <= -90)   return "disabled";
  if (freq == 1)     return "every node";
  if (freq > 1) {
    char buf[24];
    std::snprintf(buf, sizeof(buf), "every %d nodes", freq);
    return buf;
  }
  // freq == 0 or other negative: treat as disabled
  return "disabled";
}

void CbcCutGenOutput::printGeneratorTable()
{
  if (genTablePrinted_ || genInfos_.empty()) return;
  genTablePrinted_ = true;

  const CoinTable tbl = makeCgGenTable(utf8_, compact_);
  const char *bar = tableBar(utf8_, compact_);

  fprintf(fp_, "\n  Cut generator summary:\n");
  printTableOpen(fp_, tbl);

  for (const auto &g : genInfos_) {
    char densityBuf[16];
    std::snprintf(densityBuf, sizeof(densityBuf), "%.1f", g.avgDensity);
    const std::string timeBuf = g.hasTime ? fmtTime(g.time) : std::string();
    std::string nextRun = nextRunStr(g.newFreq);

    // Generator name: left-aligned, truncate if needed
    char nameBuf[CG_W_GENNAME + 1];
    std::snprintf(nameBuf, sizeof(nameBuf), "%-*s", CG_W_GENNAME, g.name.c_str());

    fprintf(fp_, "  %s%s%*d%s%*s%s%*d%s%*s%s%-*s\n",
      nameBuf, bar,
      CG_W_ROWCUTS, g.rowCuts, bar,
      CG_W_DENSITY, densityBuf, bar,
      CG_W_COLCUTS, g.colCuts, bar,
      CG_W_GENTIME, timeBuf.c_str(), bar,
      CG_W_NEXTRUN, nextRun.c_str());
  }

  printTableClose(fp_, tbl);
}

void CbcCutGenOutput::close()
{
  if (state_ == State::Closed || state_ == State::Idle) return;
  state_ = State::Closed;

  printProgressEnd(); // idempotent — close progress table if not done
  printGeneratorTable();

  // Phase end line
  char detail[256];
  if (haveSummary_) {
    std::snprintf(detail, sizeof(detail),
      "Cut generation complete%s%d cuts, obj %g%s%g in %d pass%s",
      CoinTable::dashSep(utf8_),
      sumNcuts_, sumFromObj_, utf8_ ? " \xe2\x86\x92 " : " -> ",
      sumToObj_, sumPasses_, sumPasses_ == 1 ? "" : "es");
  } else {
    std::snprintf(detail, sizeof(detail), "Cut generation complete");
  }
  fprintf(fp_, "\n%s\n", CoinTable::phaseEnd(detail, utf8_).c_str());
  fflush(fp_);
}

// ============================================================================
// CbcBnBOutput — Branch-and-bound tree tabular output
// ============================================================================

static const int BB_W_NODES   = 8;
static const int BB_W_ONTREE  = 8;
static const int BB_W_DEPTH   = 6;
static const int BB_W_BESTSOL = 15;
static const int BB_W_METHOD  = 16;
static const int BB_W_BOUND   = 15;
static const int BB_W_GAP     = 8;
static const int BB_W_TIME    = 9;

static CoinTable makeBnBTable(bool utf8, bool compact)
{
  return CoinTable({
    { "Nodes",     BB_W_NODES,   false },
    { "OnTree",    BB_W_ONTREE,  false },
    { "Depth",     BB_W_DEPTH,   false },
    { "BestSol",   BB_W_BESTSOL, false },
    { "Method",    BB_W_METHOD,  true  },
    { "BestBound", BB_W_BOUND,   false },
    { "Gap%",      BB_W_GAP,     false },
    { "Time(s)",      BB_W_TIME,    false },
  }, utf8, /*indent=*/2, compact);
}

// Format an objective value: "—" for 1e49+, otherwise %.6g
static std::string bnbFmtObj(double v, bool utf8)
{
  if (v >= 1e49) return utf8 ? "\xe2\x80\x94" : "-";
  char buf[32];
  std::snprintf(buf, sizeof(buf), "%.6g", v);
  return buf;
}

// Format the MIP optimality gap as a percentage string
static std::string bnbFmtGap(double bestSol, double bestBound, bool utf8)
{
  if (bestSol >= 1e49 || bestBound >= 1e49)
    return utf8 ? "\xe2\x80\x94" : "-";
  double ref = std::max(std::fabs(bestSol), 1e-10);
  double gap = 100.0 * std::fabs(bestSol - bestBound) / ref;
  char buf[24];
  if (gap >= 9999.5)
    std::snprintf(buf, sizeof(buf), ">9999%%");
  else
    std::snprintf(buf, sizeof(buf), "%.2f%%", gap);
  return buf;
}

// Format a large integer count compactly: 12345 → "12.3K", 1234567 → "1.23M"
static std::string bnbFmtCount(long n)
{
  char buf[24];
  if (n < 10000L)
    std::snprintf(buf, sizeof(buf), "%ld", n);
  else if (n < 1000000L)
    std::snprintf(buf, sizeof(buf), "%.1fK", n / 1000.0);
  else
    std::snprintf(buf, sizeof(buf), "%.2fM", n / 1e6);
  return buf;
}

// Shorten or truncate a heuristic method name to fit BB_W_METHOD chars
static std::string bnbMethodName(const char *raw)
{
  if (!raw || !*raw) return "B&B";
  if (std::strstr(raw, "rounding in feaspump")) return "FP rounding";
  if (std::strstr(raw, "feasibility pump"))     return "FP";
  std::string s(raw);
  if ((int)s.size() > BB_W_METHOD) s.resize(BB_W_METHOD);
  return s;
}

// ── CbcBnBOutput methods ────────────────────────────────────────────────────

CbcBnBOutput::CbcBnBOutput(FILE *fp, bool utf8, int logLevel)
  : fp_(fp), utf8_(utf8), compact_(CbcOutput::useCompact()), logLevel_(logLevel)
{}

CbcBnBOutput::~CbcBnBOutput()
{
  if (tableOpen_) closeTable();
}

void CbcBnBOutput::startPhase()
{
  inPhase_  = true;
  tableOpen_ = true;
  const CoinTable tbl = makeBnBTable(utf8_, compact_);
  fprintf(fp_, "\n%s\n\n", CoinTable::phaseStart("Branch and bound", utf8_).c_str());
  printTableOpen(fp_, tbl);
}

void CbcBnBOutput::printRow(bool isIncumbent, long nodes, int onTree, int depth,
  double bestSol, const char *method, double bestBound, double wallclock)
{
  if (!inPhase_)
    startPhase();
  else if (!tableOpen_)
    openContinuation();

  const char *bar = tableBar(utf8_, compact_);
  const char *star = isIncumbent
    ? (utf8_ ? " \xe2\x98\x85" : " *")
    : "  ";

  char nodeStr[24], onTreeStr[24], depthStr[16];
  std::snprintf(nodeStr,   sizeof(nodeStr),   "%ld", nodes);
  std::snprintf(onTreeStr, sizeof(onTreeStr), "%d",  onTree);
  std::snprintf(depthStr,  sizeof(depthStr),  "%d",  depth);
  const std::string timeStr = fmtTime(wallclock);

  const std::string bsStr  = bnbFmtObj(bestSol,   utf8_);
  const std::string bbStr  = bnbFmtObj(bestBound,  utf8_);
  const std::string gapStr = bnbFmtGap(bestSol, bestBound, utf8_);
  const std::string mStr   = method ? bnbMethodName(method) : "";

  fprintf(fp_, "%s", star);
  fprintf(fp_, "%*s%s", BB_W_NODES,  nodeStr,   bar);
  fprintf(fp_, "%*s%s", BB_W_ONTREE, onTreeStr, bar);
  fprintf(fp_, "%*s%s", BB_W_DEPTH,  depthStr,  bar);
  // BestSol: right-align with visual-width correction for UTF-8 dashes
  { int pad = BB_W_BESTSOL - utf8VisLen(bsStr); if (pad > 0) fprintf(fp_, "%*s", pad, ""); }
  fprintf(fp_, "%s%s", bsStr.c_str(), bar);
  // Method: left-align
  fprintf(fp_, "%s", mStr.c_str());
  { int vis = utf8VisLen(mStr); if (vis < BB_W_METHOD) fprintf(fp_, "%*s", BB_W_METHOD - vis, ""); }
  fprintf(fp_, "%s", bar);
  // BestBound: right-align
  { int pad = BB_W_BOUND - utf8VisLen(bbStr);  if (pad > 0) fprintf(fp_, "%*s", pad, ""); }
  fprintf(fp_, "%s%s", bbStr.c_str(), bar);
  // Gap%: right-align
  { int pad = BB_W_GAP - utf8VisLen(gapStr); if (pad > 0) fprintf(fp_, "%*s", pad, ""); }
  fprintf(fp_, "%s%s", gapStr.c_str(), bar);
  // Time: right-align
  fprintf(fp_, "%*s\n", BB_W_TIME, timeStr.c_str());
  fflush(fp_);
}

void CbcBnBOutput::onProgress(long nodes, int onTree, int depth, double bestSol,
  double bestBound, double /*elapsed*/)
{
  lastOnTree_    = onTree;
  lastDepth_     = depth;
  lastBestBound_ = bestBound;
  // Flush incumbents found during root cut gen (before first status row).
  // We use the bound from this first progress call so gap% is meaningful.
  for (auto &pi : preProgressIncumbents_)
    printRow(true, pi.nodes, onTree, depth, pi.obj, pi.method.c_str(), bestBound, pi.wallclock);
  preProgressIncumbents_.clear();
  printRow(false, nodes, onTree, depth, bestSol, nullptr, bestBound, CoinWallclockTime());
}

void CbcBnBOutput::queuePreProgressIncumbent(double obj, const char *method,
  long nodes, double /*elapsed*/)
{
  // Record wall-clock time now (when the event is intercepted) so the time
  // shown in the table reflects when the incumbent was actually found.
  preProgressIncumbents_.push_back({ obj, method ? method : "", nodes, CoinWallclockTime() });
}

void CbcBnBOutput::onBnBIncumbent(double obj, long nodes, double /*elapsed*/)
{
  printRow(true, nodes, lastOnTree_, lastDepth_, obj, "B&B", lastBestBound_, CoinWallclockTime());
}

void CbcBnBOutput::onHeurIncumbent(double obj, const char *method,
  long nodes, double /*elapsed*/)
{
  printRow(true, nodes, lastOnTree_, lastDepth_, obj, method, lastBestBound_, CoinWallclockTime());
}

void CbcBnBOutput::onStopReason(const char *reason)
{
  stopReason_ = reason ? reason : "";
}

void CbcBnBOutput::onRestartNote(int rows, int cols)
{
  if (!inPhase_) return;
  const char *dash = utf8_ ? "──" : "--";
  if (rows > 0 && cols > 0)
    fprintf(fp_, "  %s Reduced cost fixing: restarting root node (%d rows, %d cols fixed) %s\n",
      dash, rows, cols, dash);
  else
    fprintf(fp_, "  %s Reduced cost fixing: restarting root node %s\n", dash, dash);
  fflush(fp_);
  // Close the progress table so the restart phases can be printed cleanly.
  // inPhase_ stays true — openContinuation() will re-open the table when B&B resumes.
  closeTable();
}

void CbcBnBOutput::openContinuation()
{
  ++restartCount_;
  const CoinTable tbl = makeBnBTable(utf8_, compact_);
  fprintf(fp_, "\n%s\n\n", CoinTable::phaseStart("Branch and bound (continued)", utf8_).c_str());
  printTableOpen(fp_, tbl);
  tableOpen_ = true;
}

void CbcBnBOutput::closeTable()
{
  if (!tableOpen_) return;
  printTableClose(fp_, makeBnBTable(utf8_, compact_));
  tableOpen_ = false;
}

void CbcBnBOutput::onComplete(bool optimal, double bestSol, double bestBound,
  long iters, long nodes, double /*elapsed*/)
{
  // Flush pre-B&B incumbents (e.g. Proximity Search) if the time limit fired
  // before the first onProgress() call — they would otherwise be silently dropped.
  if (!preProgressIncumbents_.empty()) {
    for (auto &pi : preProgressIncumbents_)
      printRow(true, pi.nodes, lastOnTree_, lastDepth_, pi.obj, pi.method.c_str(), bestBound, pi.wallclock);
    preProgressIncumbents_.clear();
  }
  closeTable();
  inPhase_ = false;

  const double now = CoinWallclockTime();

  // Build summary line for ✔
  std::string reasonPart = stopReason_.empty() ? "" : " (" + stopReason_ + ")";
  const char *dash = CoinTable::dashSep(utf8_);
  const std::string nowStr = fmtTime(now);
  char detail[512];
  if (optimal) {
    std::snprintf(detail, sizeof(detail),
      "Optimal%s%sObj: %s   Bound: %s   Gap: %s   Nodes: %s   Iters: %s   Time: %ss",
      reasonPart.c_str(), dash,
      bnbFmtObj(bestSol,   utf8_).c_str(),
      bnbFmtObj(bestBound, utf8_).c_str(),
      bnbFmtGap(bestSol, bestBound, utf8_).c_str(),
      bnbFmtCount(nodes).c_str(),
      bnbFmtCount(iters).c_str(),
      nowStr.c_str());
  } else if (bestSol >= 1e49) {
    std::snprintf(detail, sizeof(detail),
      "Stopped%s%sNo integer solution   Bound: %s   Nodes: %s   Iters: %s   Time: %ss",
      reasonPart.c_str(), dash,
      bnbFmtObj(bestBound, utf8_).c_str(),
      bnbFmtCount(nodes).c_str(),
      bnbFmtCount(iters).c_str(),
      nowStr.c_str());
  } else {
    std::snprintf(detail, sizeof(detail),
      "Stopped%s%sBestSol: %s   Bound: %s   Gap: %s   Nodes: %s   Iters: %s   Time: %ss",
      reasonPart.c_str(), dash,
      bnbFmtObj(bestSol,   utf8_).c_str(),
      bnbFmtObj(bestBound, utf8_).c_str(),
      bnbFmtGap(bestSol, bestBound, utf8_).c_str(),
      bnbFmtCount(nodes).c_str(),
      bnbFmtCount(iters).c_str(),
      nowStr.c_str());
  }
  fprintf(fp_, "\n%s\n", CoinTable::phaseEnd(detail, utf8_).c_str());
  printFooter();
  fflush(fp_);
}

void CbcBnBOutput::onStrongStats(long calls, long iters, long fathomed, int fixed)
{
  hasStrong_  = true;
  sbCalls_    = calls;
  sbIters_    = iters;
  sbFathomed_ = fathomed;
  sbFixed_    = fixed;
}

void CbcBnBOutput::onOtherStats(int maxDepth, double djFixed)
{
  hasDepth_ = true;
  maxDepth_ = maxDepth;
  djFixed_  = djFixed;
}

void CbcBnBOutput::onOtherStats2(int maxDepth, double djFixed,
  long fathomTimes, long fathomNodes, long fathomIters)
{
  hasDepth_    = true;
  maxDepth_    = maxDepth;
  djFixed_     = djFixed;
  hasFathom_   = true;
  fathomTimes_ = fathomTimes;
  fathomNodes_ = fathomNodes;
  fathomIters_ = fathomIters;
}

void CbcBnBOutput::onOrbitalStats(int successes, double avgExtra,
  int fixed, double avgFixed)
{
  hasOrbital_   = true;
  orbSuccesses_ = successes;
  orbAvgExtra_  = avgExtra;
  orbFixed_     = fixed;
  orbAvgFixed_  = avgFixed;
}

void CbcBnBOutput::printFooter()
{
  if (hasStrong_) {
    fprintf(fp_, "  Strong branching: %s calls, %s iters, %ld nodes fathomed, %d vars fixed\n",
      bnbFmtCount(sbCalls_).c_str(), bnbFmtCount(sbIters_).c_str(),
      sbFathomed_, sbFixed_);
  }
  if (hasDepth_) {
    if (hasFathom_)
      fprintf(fp_,
        "  Max depth: %d   Vars fixed on cost: %.0f"
        "   Fathoming: %ld times, %ld nodes, %s iters\n",
        maxDepth_, djFixed_, fathomTimes_, fathomNodes_,
        bnbFmtCount(fathomIters_).c_str());
    else
      fprintf(fp_, "  Max depth: %d   Vars fixed on cost: %.0f\n",
        maxDepth_, djFixed_);
  }
  if (hasOrbital_ && orbSuccesses_ > 0)
    fprintf(fp_,
      "  Orbital branching: %d successes   avg extra: %.3f"
      "   vars fixed: %d (avg %.3f)\n",
      orbSuccesses_, orbAvgExtra_, orbFixed_, orbAvgFixed_);
}
