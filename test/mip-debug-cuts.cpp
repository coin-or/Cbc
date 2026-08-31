/*
 * mip-debug-cuts — activates Osi's row-cut debugger with a known-correct
 * reference solution (e.g. a certified optimum from mip-sanity-data) and
 * solves the instance, so that every cut, bound fixing, or branching
 * decision that would exclude that reference solution is flagged during
 * the search.
 *
 * This is meant for diagnosing WRONG_OBJ / INFEASIBLE_WRONG regressions:
 * Cbc claims optimality (or infeasibility) at the wrong value, meaning some
 * cut generator, bound-propagation pass, or branching decision incorrectly
 * excluded the true optimum. Feeding that true optimum (from a .sol file)
 * to the row-cut debugger makes CbcModel itself check, at every step, that
 * the known solution remains reachable — see CHECK_KNOWN_SOLUTION-guarded
 * code throughout CbcModel.cpp and OsiRowCutDebugger. Any violation prints
 * a diagnostic; in a --debug (non-NDEBUG) build the corresponding CoinAssert
 * additionally throws a CoinError, whose message/file/line pinpoints exactly
 * which check tripped.
 *
 * IMPORTANT: build Cbc in --debug mode (see ./config --debug --install at
 * the workspace root) for the CoinAssert checks to be active; in a
 * default --opt (NDEBUG) build they are compiled out and only the
 * printf-style diagnostics (e.g. "Conflict cut...", writeMps("bad")) fire.
 *
 * Usage:
 *   mip-debug-cuts <instance-name> [options]
 *   mip-debug-cuts <instanceFileName> <solFileName> <timeLimit> <nodesLimit> [options]
 *
 * The first form looks up <instance-name> in the mip-sanity-data checkout
 * (mips/<name>.mps.gz + sols/<name>.sol), applying limits.tsv automatically.
 * The second, generic form takes any MPS/MPS.gz/LP file plus an explicit
 * reference .sol file and time/node limits (0 = unlimited nodes).
 *
 * The .sol file is read by CbcMipStart::read(), which accepts both Cbc's
 * own "<index> <colName> <value> [<extra column, ignored>]" format and the
 * plain MIPLIB 2017 "<colName> <value>" format, one line per variable (or
 * per nonzero variable -- columns not listed default to 0). The reference
 * solution is loaded via the same "-debugCuts" parameter CbcSolver's CLI
 * uses (CbcParam::DEBUGCUTS -> Cbc_setParameter(model, "debugCuts",
 * solFile)) rather than a direct Cbc_activateRowCutDebugger() call, so
 * that when CglPreProcess renumbers/substitutes columns the reference
 * solution is remapped via CglPreProcess::originalColumns() and its
 * continuous values are re-solved for consistently before the debugger is
 * actually activated on the *processed* model (see CbcSolver.cpp's
 * DEBUGCUTS case and the babPostPreprocessCleanup()-adjacent remap block).
 * A raw, unmapped activation would otherwise report false violations from
 * essentially every generator once preprocessing changes the column
 * space. --raw-continuous is a no-op with this mechanism (kept only for
 * back-compat) since the remap step's continuous re-solve is unconditional.
 *
 * Exit codes:
 *   0  solve completed, no violation of the reference solution detected
 *   1  a cut/bound-fixing invalidating the reference solution was detected
 *      (CoinError caught) or Cbc's final answer still mismatches bks
 *   2  usage / file / instance-lookup error
 */

#include "Cbc_C_Interface.h"
#include "CoinError.hpp"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <map>
#include <sstream>
#include <string>
#include <vector>
#ifdef _WIN32
#include <windows.h>
#else
#include <unistd.h>
#endif

static std::vector<std::string> splitTab(const std::string &line)
{
  std::vector<std::string> fields;
  std::stringstream ss(line);
  std::string field;
  while (std::getline(ss, field, '\t'))
    fields.push_back(field);
  return fields;
}

static std::vector<std::string> lookupRow(const std::string &tsvPath, const std::string &key)
{
  std::ifstream in(tsvPath);
  if (!in.is_open())
    return {};
  std::string line;
  bool first = true;
  while (std::getline(in, line)) {
    if (first) { first = false; continue; }
    if (line.empty())
      continue;
    std::vector<std::string> fields = splitTab(line);
    if (!fields.empty() && fields[0] == key)
      return fields;
  }
  return {};
}

static bool fileExists(const std::string &path)
{
  std::ifstream f(path);
  return f.good();
}

static std::string exeDir()
{
  char buf[4096];
#ifdef _WIN32
  DWORD n = GetModuleFileNameA(nullptr, buf, sizeof(buf) - 1);
  if (n == 0 || n >= sizeof(buf) - 1)
    return ".";
  buf[n] = '\0';
  std::string path(buf);
  size_t slash = path.find_last_of("/\\");
#else
  ssize_t n = readlink("/proc/self/exe", buf, sizeof(buf) - 1);
  if (n <= 0)
    return ".";
  buf[n] = '\0';
  std::string path(buf);
  size_t slash = path.find_last_of('/');
#endif
  return (slash == std::string::npos) ? "." : path.substr(0, slash);
}

/* Parse a whitespace-delimited .sol file into a name -> value map. Lines
 * whose first token is not numeric (header lines, e.g. "Optimal - objective
 * value ...") are skipped. */
static bool readSolFile(const std::string &path, std::map<std::string, double> &values)
{
  std::ifstream in(path);
  if (!in.is_open())
    return false;
  std::string line;
  while (std::getline(in, line)) {
    std::istringstream ls(line);
    std::string tok0, tok1, tok2;
    if (!(ls >> tok0))
      continue;
    if (tok0.empty() || !(isdigit((unsigned char)tok0[0]) || tok0[0] == '-'))
      continue; /* not an "<index> <name> <value> ..." line */
    if (!(ls >> tok1 >> tok2))
      continue;
    char *endp = nullptr;
    double value = strtod(tok2.c_str(), &endp);
    if (endp == tok2.c_str())
      continue; /* not numeric */
    values[tok1] = value;
  }
  return true;
}

static void printUsage(const char *progName)
{
  fprintf(stderr,
    "Usage: %s <instance-name> [options]\n"
    "       %s <instanceFileName> <solFileName> <timeLimit> <nodesLimit> [options]\n"
    "\n"
    "Activates the OsiRowCutDebugger with the reference solution and solves\n"
    "the instance; any cut/bound-fixing/branch that excludes the reference\n"
    "solution is flagged (requires a --debug/non-NDEBUG Cbc build for the\n"
    "CoinAssert checks to actually raise -- see ./config --debug --install).\n"
    "\n"
    "First form: looks up <instance-name> in the mip-sanity-data checkout\n"
    "(mips/<name>.mps.gz, sols/<name>.sol, limits.tsv, bks.tsv).\n"
    "\n"
    "Second form: any MPS/MPS.gz/LP file plus an explicit reference .sol\n"
    "file and time/node limits (0 nodes = unlimited) -- no mip-sanity-data\n"
    "lookup performed.\n"
    "\n"
    "Options:\n"
    "  --data-dir=PATH        mip-sanity-data checkout (default: <exe-dir>/mip-sanity-data)\n"
    "  --threads=N            Threads passed to Cbc (default: 1)\n"
    "  --node-limit=N         Override node limit (default: from limits.tsv; 0 = unlimited)\n"
    "  --time-limit=SEC       Override time limit in seconds (default: from limits.tsv, or 120)\n"
    "  --log-level=N          Cbc_setLogLevel (default: 1, so debugger diagnostics are visible)\n"
    "  --raw-continuous       No-op with this tool's -debugCuts-based activation (kept\n"
    "                         for back-compat); continuous values are always recomputed\n"
    "                         by CbcSolver's own remap step when preprocessing changes\n"
    "                         the column count.\n"
    "  --no-preprocess        Disable CglPreProcess (-preprocess off). Reference solution\n"
    "                         remapping across preprocessing is handled correctly (via\n"
    "                         -debugCuts), so this is now purely a bisection tool: use it\n"
    "                         to tell apart a CglPreProcess-internal bug (violation only\n"
    "                         with preprocessing on) from a cut/bound-fixing bug in the\n"
    "                         main search (violation persists with preprocessing off).\n"
    "  --no-heuristics        Disable all heuristics (-heuristicsOnOff off). Heuristics run\n"
    "                         private sub-MIPs via smallBranchAndBound() with no\n"
    "                         parentModel_ set, so their onpath-trace/incumbent diagnostics\n"
    "                         are indistinguishable from a real top-level violation by\n"
    "                         parentModel_ alone. Use this to eliminate that confounder.\n"
    "  -h, --help             Show this help\n"
    "\n"
    "Exit code: 0 = OK,  1 = violation detected / obj mismatch,  2 = usage/file error\n",
    progName, progName);
}

int main(int argc, char *argv[])
{
  std::string dataDir;
  int threads = 1;
  long nodeLimitOverride = -1;
  double timeLimitOverride = -1.0;
  int logLevel = 1;
  bool rawContinuous = false;
  bool noPreprocess = false;
  bool noHeuristics = false;

  std::vector<std::string> positional;
  for (int i = 1; i < argc; ++i) {
    std::string arg = argv[i];
    auto valueOf = [&](const std::string &prefix) -> std::string {
      return arg.substr(prefix.size());
    };
    if (arg == "-h" || arg == "--help") {
      printUsage(argv[0]);
      return 0;
    } else if (arg.rfind("--data-dir=", 0) == 0) {
      dataDir = valueOf("--data-dir=");
    } else if (arg.rfind("--threads=", 0) == 0) {
      threads = atoi(valueOf("--threads=").c_str());
    } else if (arg.rfind("--node-limit=", 0) == 0) {
      nodeLimitOverride = atol(valueOf("--node-limit=").c_str());
    } else if (arg.rfind("--time-limit=", 0) == 0) {
      timeLimitOverride = atof(valueOf("--time-limit=").c_str());
    } else if (arg.rfind("--log-level=", 0) == 0) {
      logLevel = atoi(valueOf("--log-level=").c_str());
    } else if (arg == "--raw-continuous") {
      rawContinuous = true;
    } else if (arg == "--no-preprocess") {
      noPreprocess = true;
    } else if (arg == "--no-heuristics") {
      noHeuristics = true;
    } else if (!arg.empty() && arg[0] == '-') {
      fprintf(stderr, "Unknown option: %s\n", arg.c_str());
      printUsage(argv[0]);
      return 2;
    } else {
      positional.push_back(arg);
    }
  }

  if (positional.empty()) {
    printUsage(argv[0]);
    return 2;
  }

  std::string instanceFile, solFile;
  double timeLimit = 120.0;
  long nodeLimit = -1; /* -1 == unlimited */
  double expectedObj = 0.0;
  bool haveExpectedObj = false;

  if (positional.size() >= 4) {
    /* Generic form: file solFile timeLimit nodesLimit */
    instanceFile = positional[0];
    solFile = positional[1];
    timeLimit = atof(positional[2].c_str());
    long nl = atol(positional[3].c_str());
    nodeLimit = (nl == 0) ? -1 : nl;
  } else if (positional.size() == 1) {
    /* mip-sanity-data lookup form */
    std::string name = positional[0];
    if (dataDir.empty())
      dataDir = exeDir() + "/mip-sanity-data";

    instanceFile = dataDir + "/mips/" + name + ".mps.gz";
    if (!fileExists(instanceFile)) {
      fprintf(stderr, "Instance file not found: %s\n", instanceFile.c_str());
      return 2;
    }
    solFile = dataDir + "/sols/" + name + ".sol";
    if (!fileExists(solFile)) {
      fprintf(stderr, "Reference solution file not found: %s\n", solFile.c_str());
      return 2;
    }

    std::vector<std::string> limitsRow = lookupRow(dataDir + "/limits.tsv", name);
    if (limitsRow.size() >= 5) {
      timeLimit = atof(limitsRow[4].c_str());
    }
    std::vector<std::string> bksRow = lookupRow(dataDir + "/bks.tsv", name);
    if (bksRow.size() >= 3) {
      expectedObj = atof(bksRow[2].c_str());
      haveExpectedObj = true;
    }
  } else {
    fprintf(stderr, "Wrong number of positional arguments.\n");
    printUsage(argv[0]);
    return 2;
  }

  if (nodeLimitOverride >= 0)
    nodeLimit = (nodeLimitOverride == 0) ? -1 : nodeLimitOverride;
  if (timeLimitOverride >= 0.0)
    timeLimit = timeLimitOverride;

  std::map<std::string, double> refValues;
  if (!readSolFile(solFile, refValues)) {
    fprintf(stderr, "Unable to open reference solution file: %s\n", solFile.c_str());
    return 2;
  }
  if (refValues.empty()) {
    fprintf(stderr, "No variable values parsed from reference solution file: %s\n", solFile.c_str());
    return 2;
  }

  Cbc_Model *model = Cbc_newModel();
  Cbc_setLogLevel(model, logLevel);

  int mpsStatus = Cbc_readMps(model, instanceFile.c_str());
  if (mpsStatus != 0) {
    fprintf(stderr, "Error reading MPS file: %s (status=%d)\n", instanceFile.c_str(), mpsStatus);
    Cbc_deleteModel(model);
    return 2;
  }

  int nCols = Cbc_getNumCols(model);
  int matched = 0;
  for (int j = 0; j < nCols; ++j) {
    char nameBuf[512];
    Cbc_getColName(model, j, nameBuf, sizeof(nameBuf));
    if (refValues.find(nameBuf) != refValues.end())
      ++matched;
  }
  printf("[mip-debug-cuts] matched %d of %d reference values to model columns "
         "(%d columns default to 0)\n",
    matched, (int)refValues.size(), nCols - matched);
  if (matched < (int)refValues.size()) {
    printf("[mip-debug-cuts] WARNING: %d reference value(s) did not match any "
           "column name -- check the .sol file matches this instance.\n",
      (int)refValues.size() - matched);
  }

  // Route the reference solution through the same "-debugCuts" parameter
  // CbcSolver's CLI uses (CbcParam::DEBUGCUTS), rather than calling
  // Cbc_activateRowCutDebugger() directly here. CglPreProcess renumbers,
  // substitutes and drops columns, so a debugger activated on the
  // as-loaded solver (indexed in the ORIGINAL column space) is checked
  // against the WRONG columns once preprocessing has run -- every
  // generator that then runs on the renumbered model appears to violate
  // it, drowning out genuine violations in false positives. "-debugCuts"
  // avoids this: it name-matches the reference solution against the
  // original solver once (see CbcSolver.cpp's DEBUGCUTS case), then, after
  // preprocessing, remaps it via CglPreProcess::originalColumns() and
  // re-solves the LP with integers fixed to get consistent continuous
  // values before activating the debugger on the *processed* model (see
  // the "if (debugValues)" block feeding babModel_->solver()->
  // activateRowCutDebugger() in CbcSolver.cpp) -- CbcMipStart's file
  // reader (used by both -mipStart and -debugCuts) already accepts this
  // tool's "<idx> <name> <value>" / MIPLIB "<name> <value>" formats
  // directly, so the same solFile can be handed to it unchanged.
  Cbc_setParameter(model, "debugCuts", solFile.c_str());
  printf("[mip-debug-cuts] Row cut debugger will be activated (via -debugCuts, "
         "remapped across CglPreProcess if it renumbers columns) with the "
         "%d-column reference solution. Every cut/bound-fixing that excludes "
         "it will be flagged below.\n\n",
    nCols);
  if (rawContinuous) {
    printf("[mip-debug-cuts] NOTE: --raw-continuous has no effect with the "
           "-debugCuts-based activation used here -- continuous values are "
           "always recomputed by CbcSolver's own remap step when "
           "preprocessing changes the column count (a no-op for instances "
           "with 0 continuous variables, as here).\n");
  }


  Cbc_setDblParam(model, DBL_PARAM_TIME_LIMIT, timeLimit);
  Cbc_setParameter(model, "threads", std::to_string(threads).c_str());
  if (nodeLimit >= 0)
    Cbc_setParameter(model, "maxNodes", std::to_string(nodeLimit).c_str());
  if (noPreprocess) {
    // Reference-solution remapping across preprocessing is now handled
    // correctly (via -debugCuts, see above), so this flag is purely a
    // bisection aid: if a violation only appears with preprocessing on,
    // it points at CglPreProcess itself; if it persists with preprocessing
    // off, the bug is in the main-search cut/bound-fixing machinery.
    Cbc_setParameter(model, "preprocess", "off");
    printf("[mip-debug-cuts] --no-preprocess: preprocessing disabled -- any\n"
           "violation that only occurred with preprocessing on points at\n"
           "CglPreProcess; one that persists here is a main-search bug.\n");
  }
  if (noHeuristics) {
    // Heuristics (FPump, RINS, Dive*, Local Search, ...) run their own
    // private CbcModel/solver clones via CbcHeuristic::smallBranchAndBound(),
    // which do NOT set parentModel_ -- so an [onpath-trace]/[incumbent]
    // diagnostic from one of these sub-MIPs is indistinguishable from a
    // real top-level violation by parentModel_ alone (both show up as
    // "parentModel_=(nil)"). Disabling all heuristics removes that
    // confounder entirely: any remaining onpath-trace hit is then
    // guaranteed to be from the real top-level search tree.
    Cbc_setParameter(model, "heuristicsOnOff", "off");
    printf("[mip-debug-cuts] --no-heuristics: all heuristics disabled -- any\n"
           "onpath-trace/incumbent diagnostic remaining is from the real\n"
           "top-level search tree, not a heuristic's private sub-MIP.\n");
  }

  bool violationDetected = false;
  try {
    Cbc_solve(model);
  } catch (CoinError &e) {
    violationDetected = true;
    printf("\n[mip-debug-cuts] *** CoinError caught during solve ***\n");
    printf("[mip-debug-cuts] %s:%d method %s: assertion '%s' failed.\n",
      e.fileName().c_str(), e.lineNumber(), e.methodName().c_str(),
      e.message().c_str());
    if (!e.className().empty())
      printf("[mip-debug-cuts] Possible reason: %s\n", e.className().c_str());
  }

  bool isProven = Cbc_isProvenOptimal(model) != 0;
  double obj = Cbc_getObjValue(model);

  printf("\n[mip-debug-cuts] result: %s  obj=%.10g",
    isProven ? "proven optimal" : "not proven", obj);
  if (haveExpectedObj)
    printf("  (certified=%.10g)", expectedObj);
  printf("\n");

  bool objMismatch = false;
  if (isProven && haveExpectedObj) {
    double tol = std::max(1e-4, 1e-4 * fabs(expectedObj));
    if (fabs(obj - expectedObj) > tol) {
      objMismatch = true;
      printf("[mip-debug-cuts] OBJECTIVE MISMATCH vs certified best-known "
             "(delta=%.6g, tol=%.4g) -- yet no cut/bound-fixing violation was "
             "detected against the reference solution. This suggests the bug "
             "is NOT a classically-invalid cut, but rather LP-bound inflation "
             "(e.g. floating-point rounding pushing the relaxation bound "
             "above the fathoming threshold) or a similar non-cut-level "
             "issue.\n",
        obj - expectedObj, tol);
    }
  }

  Cbc_deleteModel(model);

  if (violationDetected || objMismatch)
    return 1;
  return 0;
}
