/*
 * mip-c-solve — manual driver to solve one mip-sanity-data instance through
 * Cbc's C interface (Cbc_C_Interface.h) and exercise/validate the API that
 * matters for that path: loading the instance, applying the suggested
 * node/time limits from limits.tsv, calling Cbc_solve(), then querying the
 * best solution and the whole solution pool (Cbc_bestSolution /
 * Cbc_numberSavedSolutions / Cbc_savedSolution / Cbc_savedSolutionObj) and
 * the final status accessors (Cbc_status, Cbc_secondaryStatus,
 * Cbc_isProvenOptimal, Cbc_isProvenInfeasible, ...).
 *
 * Unlike cbc_validate_sol (which validates a .sol file already produced by
 * the `cbc` command-line executable), this tool solves the instance itself,
 * in-process, purely through the C API — useful to manually sanity-check
 * that interface end to end. Cbc's own solve messages (from Cbc_readMps and
 * Cbc_solve) are deliberately left unsuppressed, so this doubles as a way to
 * eyeball exactly what Cbc prints during a normal solve.
 *
 * After Cbc_solve() returns, this tool:
 *   a) checks the feasibility of every produced solution — the best solution
 *      and every solution saved in the pool — using Cbc_checkFeasibility
 *      (bounds, integrality, row activity), reporting any violation;
 *   b) recomputes each solution's objective from its vector and the model's
 *      objective coefficients, and compares it against the value reported by
 *      the C interface (Cbc_getObjValue / Cbc_savedSolutionObj);
 *   c) if optimality is claimed (Cbc_isProvenOptimal), cross-checks the
 *      reported objective against mip-sanity-data's bks.tsv best-known value
 *      (and, if infeasibility is claimed, cross-checks against bks.tsv's
 *      expected status too).
 *
 * Usage:
 *   mip-c-solve <instance-name> [options]
 *   mip-c-solve <instanceFileName> <timeLimit> <nodesLimit> [options]
 *
 * The first form looks up <instance-name> in the mip-sanity-data checkout
 * (mips/<name>.mps.gz), applying limits.tsv/bks.tsv automatically. The
 * second, generic form takes any MPS/MPS.gz/LP file directly plus an
 * explicit time limit (seconds) and node limit (0 = unlimited) — handy for
 * quickly checking, against any instance, that the C interface solves it
 * the same way the `cbc` command line would (e.g.
 * `cbc <file> -sec <timeLimit> -maxNodes <nodesLimit> -solve`).
 *
 * Exit codes:
 *   0  solve completed and every check above passed
 *   1  one or more violations/mismatches found
 *   2  usage / file / instance-lookup error
 */

#include "Cbc_C_Interface.h"
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

/* ── Tolerances (same defaults as cbc_validate_sol) ────────────────────── */

static double OBJ_TOL = 1e-4;     /* relative: C-interface obj vs recomputed */
static double BKS_ABS_TOL = 1e-4; /* absolute tolerance vs best-known objective */
static double BKS_PCT_TOL = 0.01; /* percent tolerance vs best-known objective */

/* ── Small TSV lookup helper ────────────────────────────────────────────── */

static std::vector<std::string> splitTab(const std::string &line)
{
  std::vector<std::string> fields;
  std::stringstream ss(line);
  std::string field;
  while (std::getline(ss, field, '\t'))
    fields.push_back(field);
  return fields;
}

/* Returns the row (split on tab) whose first column equals key, or an empty
 * vector if not found / file missing. */
static std::vector<std::string> lookupRow(const std::string &tsvPath, const std::string &key)
{
  std::ifstream in(tsvPath);
  if (!in.is_open())
    return {};
  std::string line;
  bool first = true;
  while (std::getline(in, line)) {
    if (first) { first = false; continue; } /* skip header */
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

/* Parses mip-sanity-data's sols/<instance>.sol reference-solution format:
 *   # Solution for model ...
 *   # Objective value = ...
 *   <colName> <value>
 *   ...
 * (comment lines start with '#'; one "name value" pair per remaining line).
 * Returns a map from column name to value, empty if the file is missing,
 * empty, or unreadable. */
static std::map<std::string, double> parseReferenceSol(const std::string &path)
{
  std::map<std::string, double> values;
  std::ifstream in(path);
  if (!in.is_open())
    return values;
  std::string line;
  while (std::getline(in, line)) {
    if (line.empty() || line[0] == '#')
      continue;
    std::istringstream ss(line);
    std::string name;
    double value;
    if (ss >> name >> value)
      values[name] = value;
  }
  return values;
}

/* True if path's extension (ignoring a trailing .gz) is ".lp" (case-insensitive). */
static bool looksLikeLpFile(const std::string &path)
{
  std::string p = path;
  std::transform(p.begin(), p.end(), p.begin(), ::tolower);
  if (p.size() >= 3 && p.substr(p.size() - 3) == ".gz")
    p = p.substr(0, p.size() - 3);
  return p.size() >= 3 && p.substr(p.size() - 3) == ".lp";
}


/* Directory containing this executable (used to locate the default
 * mip-sanity-data checkout next to it, mirroring run-mip-sanity-tests). */
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

static void printUsage(const char *progName)
{
  fprintf(stderr,
    "Usage: %s <instance-name> [options]\n"
    "       %s <instanceFileName> <timeLimit> <nodesLimit> [options]\n"
    "\n"
    "First form: solves one mip-sanity-data instance through Cbc's C\n"
    "interface (Cbc_C_Interface.h), applying the suggested limits.tsv\n"
    "node/time limits, then validates the best solution and the whole\n"
    "solution pool (feasibility, objective consistency, optimality/\n"
    "infeasibility claims vs bks.tsv).\n"
    "\n"
    "Second form: solves any MPS/MPS.gz/LP file directly, with an explicit\n"
    "time limit (seconds) and node limit (0 = unlimited) — no mip-sanity-data\n"
    "lookup is performed. Useful to quickly check that the C interface solves\n"
    "an arbitrary instance the same way the `cbc` command line would (e.g.\n"
    "`cbc <file> -sec <timeLimit> -maxNodes <nodesLimit> -solve`).\n"
    "\n"
    "Options:\n"
    "  --stage=STAGE          Which pipeline stage to exercise instead of a\n"
    "                         full solve (default: solve). One of:\n"
    "                           solve              Full Cbc_solve() (default)\n"
    "                           bound-propagation  Cbc_propagateBounds() only --\n"
    "                             reports how many column bounds were tightened\n"
    "                             and whether infeasibility was proved, then exits\n"
    "                             (no branch-and-bound is run). If bks.tsv marks\n"
    "                             this instance optimal and a reference solution\n"
    "                             exists (sols/<instance>.sol), also checks that\n"
    "                             the tightened bounds do not exclude that\n"
    "                             solution's values.\n"
    "  --data-dir=PATH        mip-sanity-data checkout (default: <exe-dir>/mip-sanity-data)\n"
    "                         (ignored in the <instanceFileName> <timeLimit> <nodesLimit> form)\n"
    "  --threads=N            Threads passed to Cbc (default: 1)\n"
    "  --node-limit=N         Override node limit (default: from limits.tsv; 0 = unlimited)\n"
    "  --time-limit=SEC       Override time limit in seconds (default: from limits.tsv, or 120)\n"
    "  --log-level=N          Cbc_setLogLevel (default: leave Cbc's own default)\n"
    "  --expected-status=S    Override bks.tsv expected status (optimal|infeasible)\n"
    "  --expected-obj=VAL     Override bks.tsv expected objective\n"
    "  --verbose              Also print nonzero variable values for every solution\n"
    "  -h, --help             Show this help\n"
    "\n"
    "Exit code: 0 = OK,  1 = violations/mismatches found,  2 = usage/file error\n",
    progName, progName);
}

/* ── Feasibility check + report for one solution vector ────────────────── */

static bool checkOneSolution(Cbc_Model *m, const double *x, const char *label,
    double reportedObj, bool haveReportedObj, bool verbose,
    std::vector<std::string> &violations)
{
  int ncols = Cbc_getNumCols(m);
  double maxViolRow, maxViolCol;
  int rowIdx, colIdx;
  char feasible = Cbc_checkFeasibility(m, x, &maxViolRow, &rowIdx, &maxViolCol, &colIdx);

  const double *objCoef = Cbc_getObjCoefficients(m);
  double computedObj = 0.0;
  for (int j = 0; j < ncols; ++j)
    computedObj += objCoef[j] * x[j];

  printf("  [%s] computed obj = %.10g", label, computedObj);
  if (haveReportedObj)
    printf("   (C-interface reported = %.10g)", reportedObj);
  printf("\n");

  if (verbose) {
    printf("    nonzero variables:");
    int printed = 0;
    for (int j = 0; j < ncols && printed < 40; ++j) {
      if (fabs(x[j]) > 1e-9) {
        char nameBuf[512];
        Cbc_getColName(m, j, nameBuf, sizeof(nameBuf));
        printf(" %s=%.6g", nameBuf, x[j]);
        ++printed;
      }
    }
    if (printed == 40)
      printf(" ...");
    printf("\n");
  }

  bool ok = true;

  if (!feasible) {
    ok = false;
    char buf[512];
    if (rowIdx >= 0) {
      char rowName[512];
      Cbc_getRowName(m, rowIdx, rowName, sizeof(rowName));
      snprintf(buf, sizeof(buf),
        "  %s: row feasibility violated, row %d (%s), max violation=%.3e",
        label, rowIdx, rowName, maxViolRow);
      violations.push_back(buf);
    }
    if (colIdx >= 0) {
      char colName[512];
      Cbc_getColName(m, colIdx, colName, sizeof(colName));
      snprintf(buf, sizeof(buf),
        "  %s: column bound/integrality violated, col %d (%s), max violation=%.3e",
        label, colIdx, colName, maxViolCol);
      violations.push_back(buf);
    }
  } else {
    printf("    feasibility: OK (max row viol=%.3e, max col viol=%.3e)\n", maxViolRow, maxViolCol);
  }

  if (haveReportedObj) {
    double denom = std::max(1.0, fabs(reportedObj));
    double objDelta = fabs(computedObj - reportedObj) / denom;
    if (objDelta > OBJ_TOL) {
      ok = false;
      char buf[256];
      snprintf(buf, sizeof(buf),
        "  %s: OBJ MISMATCH computed=%.10g reported=%.10g rel_delta=%.3e (exceeds tol %.1e)",
        label, computedObj, reportedObj, objDelta, OBJ_TOL);
      violations.push_back(buf);
    }
  }

  return ok;
}

int main(int argc, char *argv[])
{
  std::string dataDir;
  std::string stage = "solve";
  int threads = 1;
  long nodeLimitOverride = -1; /* -1 = use limits.tsv */
  double timeLimitOverride = -1.0;
  int logLevelOverride = -1;
  std::string expectedStatusOverride;
  double expectedObjOverride = 0.0;
  bool haveExpectedObjOverride = false;
  bool verbose = false;

  std::vector<std::string> positional;
  for (int i = 1; i < argc; ++i) {
    std::string arg = argv[i];
    auto valueOf = [&](const std::string &prefix) -> std::string {
      return arg.substr(prefix.size());
    };
    if (arg == "-h" || arg == "--help") {
      printUsage(argv[0]);
      return 2;
    } else if (arg.rfind("--data-dir=", 0) == 0) {
      dataDir = valueOf("--data-dir=");
    } else if (arg.rfind("--stage=", 0) == 0) {
      stage = valueOf("--stage=");
      if (stage != "solve" && stage != "bound-propagation") {
        fprintf(stderr, "Error: unknown --stage value '%s' (expected 'solve' or 'bound-propagation')\n\n",
          stage.c_str());
        printUsage(argv[0]);
        return 2;
      }
    } else if (arg.rfind("--threads=", 0) == 0) {
      threads = atoi(valueOf("--threads=").c_str());
    } else if (arg.rfind("--node-limit=", 0) == 0) {
      nodeLimitOverride = atol(valueOf("--node-limit=").c_str());
    } else if (arg.rfind("--time-limit=", 0) == 0) {
      timeLimitOverride = atof(valueOf("--time-limit=").c_str());
    } else if (arg.rfind("--log-level=", 0) == 0) {
      logLevelOverride = atoi(valueOf("--log-level=").c_str());
    } else if (arg.rfind("--expected-status=", 0) == 0) {
      expectedStatusOverride = valueOf("--expected-status=");
      std::transform(expectedStatusOverride.begin(), expectedStatusOverride.end(),
        expectedStatusOverride.begin(), ::tolower);
    } else if (arg.rfind("--expected-obj=", 0) == 0) {
      expectedObjOverride = atof(valueOf("--expected-obj=").c_str());
      haveExpectedObjOverride = true;
    } else if (arg == "--verbose") {
      verbose = true;
    } else if (!arg.empty() && arg[0] == '-') {
      fprintf(stderr, "Error: unknown option '%s'\n\n", arg.c_str());
      printUsage(argv[0]);
      return 2;
    } else {
      positional.push_back(arg);
    }
  }

  const bool directMode = (positional.size() == 3);
  if (positional.size() != 1 && !directMode) {
    fprintf(stderr, "Error: expected either one instance name, or "
      "<instanceFileName> <timeLimit> <nodesLimit>\n\n");
    printUsage(argv[0]);
    return 2;
  }

  std::string instance;   /* label used only for printing/errors */
  std::string mpsPath;    /* file actually passed to Cbc_readMps/Cbc_readLp */
  long nodeLimit = 0;     /* 0 = unlimited, matches run-mip-sanity-tests default */
  double timeLimit = 120.0;
  double hardKill = -1.0; /* -1 = not applicable (direct mode) */
  std::string expectedStatus;
  bool haveExpectedObj = false;
  double expectedObj = 0.0;

  if (directMode) {
    /* Generic form: mip-c-solve <instanceFileName> <timeLimit> <nodesLimit> */
    mpsPath = positional[0];
    instance = mpsPath;
    if (!fileExists(mpsPath)) {
      fprintf(stderr, "Error: instance file not found: %s\n", mpsPath.c_str());
      return 2;
    }
    timeLimit = atof(positional[1].c_str());
    nodeLimit = atol(positional[2].c_str());
    if (nodeLimitOverride >= 0)
      nodeLimit = nodeLimitOverride;
    if (timeLimitOverride >= 0.0)
      timeLimit = timeLimitOverride;
    if (!expectedStatusOverride.empty())
      expectedStatus = expectedStatusOverride;
    if (haveExpectedObjOverride) {
      expectedObj = expectedObjOverride;
      haveExpectedObj = true;
    }
  } else {
    /* mip-sanity-data instance-name lookup form */
    instance = positional[0];

    if (dataDir.empty())
      dataDir = exeDir() + "/mip-sanity-data";

    mpsPath = dataDir + "/mips/" + instance + ".mps.gz";
    if (!fileExists(mpsPath)) {
      fprintf(stderr, "Error: instance file not found: %s\n", mpsPath.c_str());
      fprintf(stderr, "       (looked under data-dir '%s' — pass --data-dir=PATH to override)\n",
        dataDir.c_str());
      return 2;
    }

    /* Suggested limits (limits.tsv), with CLI overrides */
    hardKill = 180.0;
    std::vector<std::string> limRow = lookupRow(dataDir + "/limits.tsv", instance);
    if (limRow.size() >= 5) {
      nodeLimit = atol(limRow[2].c_str());
      timeLimit = atof(limRow[3].c_str());
      hardKill = atof(limRow[4].c_str());
    }
    if (nodeLimitOverride >= 0)
      nodeLimit = nodeLimitOverride;
    if (timeLimitOverride >= 0.0)
      timeLimit = timeLimitOverride;

    /* Best-known reference (bks.tsv), with CLI overrides */
    std::vector<std::string> bksRow = lookupRow(dataDir + "/bks.tsv", instance);
    if (bksRow.size() >= 3) {
      expectedStatus = bksRow[1];
      std::transform(expectedStatus.begin(), expectedStatus.end(), expectedStatus.begin(), ::tolower);
      expectedObj = atof(bksRow[2].c_str());
      haveExpectedObj = true;
    }
    if (!expectedStatusOverride.empty())
      expectedStatus = expectedStatusOverride;
    if (haveExpectedObjOverride) {
      expectedObj = expectedObjOverride;
      haveExpectedObj = true;
    }
  }

  printf("=== mip-c-solve: %s ===\n", instance.c_str());
  printf("  instance file : %s\n", mpsPath.c_str());
  printf("  threads       : %d\n", threads);
  printf("  node limit    : %s\n", nodeLimit > 0 ? std::to_string(nodeLimit).c_str() : "unlimited");
  if (hardKill >= 0.0)
    printf("  time limit    : %.0f s  (hard-kill reference: %.0f s, not enforced here)\n", timeLimit, hardKill);
  else
    printf("  time limit    : %.0f s\n", timeLimit);
  if (!expectedStatus.empty())
    printf("  bks.tsv       : status=%s%s\n", expectedStatus.c_str(),
      haveExpectedObj ? (std::string("  objective=") + std::to_string(expectedObj)).c_str() : "");
  else if (!directMode)
    printf("  bks.tsv       : no reference found for this instance\n");
  printf("\n--- Cbc output (unmodified) ------------------------------------------\n");

  /* ── Load & solve, deliberately NOT silencing Cbc's own messages ───────── */
  Cbc_Model *m = Cbc_newModel();
  Cbc_storeNameIndexes(m, 1); /* must precede readMps/readLp so column/row
                                 name lookups (Cbc_getColNameIndex(), used
                                 below by the bound-propagation stage's
                                 reference-solution check) are populated */
  const bool isLp = looksLikeLpFile(mpsPath);
  int readErr = isLp ? Cbc_readLp(m, mpsPath.c_str()) : Cbc_readMps(m, mpsPath.c_str());
  if (readErr) {
    fprintf(stderr, "Error: %s failed to read '%s' (code %d)\n",
      isLp ? "Cbc_readLp" : "Cbc_readMps", mpsPath.c_str(), readErr);
    Cbc_deleteModel(m);
    return 2;
  }

  if (logLevelOverride >= 0)
    Cbc_setLogLevel(m, logLevelOverride);

  /* ── --stage=bound-propagation: exercise Cbc_propagateBounds() only,   ──
   * skipping the full Cbc_solve() path entirely (no branch-and-bound). */
  if (stage == "bound-propagation") {
    int ncols = Cbc_getNumCols(m);
    std::vector<double> lbBefore(ncols), ubBefore(ncols);
    {
      const double *lb0 = Cbc_getColLower(m);
      const double *ub0 = Cbc_getColUpper(m);
      for (int j = 0; j < ncols; ++j) {
        lbBefore[j] = lb0[j];
        ubBefore[j] = ub0[j];
      }
    }

    int rc = Cbc_propagateBounds(m);

    int nTightened = 0;
    {
      const double *lb1 = Cbc_getColLower(m);
      const double *ub1 = Cbc_getColUpper(m);
      for (int j = 0; j < ncols; ++j) {
        if (lb1[j] > lbBefore[j] + 1e-9 || ub1[j] < ubBefore[j] - 1e-9)
          ++nTightened;
      }
    }

    printf("------------------------------------------------------------------------\n\n");
    printf("=== Bound propagation stage ===\n");
    printf("  Cbc_propagateBounds()  = %d (%s)\n", rc,
      rc == 0 ? "ok, no infeasibility detected" : "infeasibility proved");
    printf("  columns tightened     = %d / %d\n", nTightened, ncols);

    std::vector<std::string> violations;
    if (rc != 0 && expectedStatus == "optimal") {
      violations.push_back(
        "MISMATCH: Cbc_propagateBounds() claims infeasibility, but bks.tsv "
        "says this instance has an optimal solution. This is a solver bug.");
    }

    /* If this instance's bks.tsv status is a confirmed "optimal" and a
     * reference solution is available (mip-sanity-data/sols/<instance>.sol),
     * make sure bound propagation did not tighten any variable's bounds past
     * that solution's value: the updated bounds must not exclude the known
     * optimal solution's values -- if they do, bound propagation is
     * incorrectly cutting off the true optimum, a serious correctness bug
     * (not just a performance one, since it can make CBC miss/misreport the
     * optimal solution entirely on this and any other instance triggering
     * the same faulty propagation rule). Only meaningful in instance-lookup
     * mode (direct mode has no known instance name to look a reference
     * solution up by).
     */
    if (!directMode && expectedStatus == "optimal") {
      std::string solPath = dataDir + "/sols/" + instance + ".sol";
      std::map<std::string, double> refSol = parseReferenceSol(solPath);
      if (!refSol.empty()) {
        const double *lb1 = Cbc_getColLower(m);
        const double *ub1 = Cbc_getColUpper(m);
        const double tol = 1e-6;
        int nChecked = 0, nViolated = 0;
        for (const auto &kv : refSol) {
          int col = Cbc_getColNameIndex(m, kv.first.c_str());
          if (col < 0)
            continue;
          ++nChecked;
          double v = kv.second;
          if (v < lb1[col] - tol || v > ub1[col] + tol) {
            ++nViolated;
            char buf[512];
            snprintf(buf, sizeof(buf),
              "MISMATCH: bound propagation excludes the reference optimal "
              "solution -- column '%s' = %.10g falls outside its tightened "
              "bounds [%.10g, %.10g] (sols/%s.sol). This is a solver bug.",
              kv.first.c_str(), v, lb1[col], ub1[col], instance.c_str());
            violations.push_back(buf);
          }
        }
        printf("  reference solution     = sols/%s.sol (%d/%zu variables checked, %d violation(s))\n",
          instance.c_str(), nChecked, refSol.size(), nViolated);
      } else {
        printf("  reference solution     = sols/%s.sol not found — skipping optimum-exclusion check\n",
          instance.c_str());
      }
    }
    printf("\n");

    bool ok = violations.empty();
    if (ok) {
      printf("RESULT: OK\n");
    } else {
      printf("RESULT: FAILED\n\n");
      printf("--- VIOLATIONS FOUND (%zu) -----------------------------------------\n", violations.size());
      for (const auto &v : violations)
        printf("%s\n", v.c_str());
      printf("\n");
    }

    Cbc_deleteModel(m);
    return ok ? 0 : 1;
  }

  char threadsStr[32];
  snprintf(threadsStr, sizeof(threadsStr), "%d", threads);
  Cbc_setParameter(m, "threads", threadsStr);
  if (nodeLimit > 0)
    Cbc_setMaximumNodes(m, (int)nodeLimit);
  Cbc_setMaximumSeconds(m, timeLimit);

  Cbc_solve(m);

  printf("------------------------------------------------------------------------\n\n");

  /* ── Status accessors ───────────────────────────────────────────────────── */
  int status = Cbc_status(m);
  int secStatus = Cbc_secondaryStatus(m);
  int isOptimal = Cbc_isProvenOptimal(m);
  int isInfeasible = Cbc_isProvenInfeasible(m);
  int isAbandoned = Cbc_isAbandoned(m);
  int nodeCount = Cbc_getNodeCount(m);
  int numSaved = Cbc_numberSavedSolutions(m);
  const double *best = Cbc_bestSolution(m);

  printf("=== Status (queried via the C interface) ===\n");
  printf("  Cbc_status()            = %d\n", status);
  printf("  Cbc_secondaryStatus()   = %d\n", secStatus);
  printf("  Cbc_isProvenOptimal()   = %d\n", isOptimal);
  printf("  Cbc_isProvenInfeasible()= %d\n", isInfeasible);
  printf("  Cbc_isAbandoned()       = %d\n", isAbandoned);
  printf("  Cbc_getNodeCount()      = %d\n", nodeCount);
  printf("  Cbc_numberSavedSolutions() = %d\n", numSaved);
  if (best)
    printf("  Cbc_getObjValue()       = %.10g\n", Cbc_getObjValue(m));
  printf("  Cbc_getBestPossibleObjValue() = %.10g\n", Cbc_getBestPossibleObjValue(m));
  printf("\n");

  std::vector<std::string> violations;

  /* ── a) + b) feasibility & objective consistency for every solution ────── */
  if (best) {
    printf("=== Best solution ===\n");
    checkOneSolution(m, best, "best", Cbc_getObjValue(m), true, verbose, violations);
    printf("\n");
  } else if (isInfeasible) {
    printf("No integer-feasible solution (infeasibility proven) — nothing to check.\n\n");
  } else {
    printf("No integer-feasible solution was found — nothing to check.\n\n");
  }

  if (numSaved > 0) {
    printf("=== Solution pool (%d saved solution(s)) ===\n", numSaved);
    for (int i = 0; i < numSaved; ++i) {
      const double *sol = Cbc_savedSolution(m, i);
      double obj = Cbc_savedSolutionObj(m, i);
      char label[64];
      snprintf(label, sizeof(label), "pool #%d", i);
      if (sol)
        checkOneSolution(m, sol, label, obj, true, verbose, violations);
      else
        printf("  [%s] Cbc_savedSolution returned NULL\n", label);
    }
    printf("\n");
  }

  /* ── c) optimality / infeasibility claims vs bks.tsv ───────────────────── */
  printf("=== Cross-check vs bks.tsv ===\n");
  if (isOptimal) {
    if (expectedStatus == "infeasible") {
      violations.push_back(
        "MISMATCH: Cbc_isProvenOptimal() is true, but bks.tsv says this instance "
        "is infeasible. This is a solver bug.");
    } else if (expectedStatus == "optimal" && haveExpectedObj) {
      double reportedObj = Cbc_getObjValue(m);
      double tolAbs = std::max(BKS_ABS_TOL, fabs(expectedObj) * BKS_PCT_TOL / 100.0);
      double delta = fabs(reportedObj - expectedObj);
      if (delta > tolAbs) {
        char buf[512];
        snprintf(buf, sizeof(buf),
          "MISMATCH: claimed-optimal objective %.10g does not match best-known "
          "value %.10g (delta=%.3e exceeds tol=%.3e = max(abs=%.1e, pct=%.4g%%))",
          reportedObj, expectedObj, delta, tolAbs, BKS_ABS_TOL, BKS_PCT_TOL);
        violations.push_back(buf);
      } else {
        printf("  Optimal objective matches best-known value (delta=%.3e <= tol=%.3e): OK\n",
          delta, tolAbs);
      }
    } else {
      printf("  Cbc claims optimality; no comparable bks.tsv reference available.\n");
    }
  } else if (isInfeasible) {
    if (expectedStatus == "optimal") {
      violations.push_back(
        "MISMATCH: Cbc_isProvenInfeasible() is true, but bks.tsv says this "
        "instance has an optimal solution. This is a solver bug.");
    } else if (expectedStatus == "infeasible") {
      printf("  Infeasibility claim matches bks.tsv: OK\n");
    } else {
      printf("  Cbc claims infeasibility; no comparable bks.tsv reference available.\n");
    }
  } else {
    printf("  Run stopped without an optimality/infeasibility proof (status=%d, secondary=%d) — "
      "no bks.tsv cross-check applies.\n", status, secStatus);
  }
  printf("\n");

  /* ── Results ────────────────────────────────────────────────────────────── */
  bool ok = violations.empty();
  if (ok) {
    printf("RESULT: OK\n");
  } else {
    printf("RESULT: FAILED\n\n");
    printf("--- VIOLATIONS FOUND (%zu) -----------------------------------------\n", violations.size());
    for (const auto &v : violations)
      printf("%s\n", v.c_str());
    printf("\n");
  }

  Cbc_deleteModel(m);
  return ok ? 0 : 1;
}
