/*
 * cbc_validate_sol — solution feasibility/optimality checker for Cbc's
 * mip-sanity-data regression suite.
 *
 * Ported from MIPster's `mipster_validate_sol` (h-g-s/mipster) and extended
 * with status-aware classification and best-known-solution (bks) cross
 * checking, driven by this repo's `test/run-mip-sanity-tests` harness.
 *
 * Given a problem file and the `.sol` file Cbc wrote for it (via
 * `-solve -solu <file>`), this tool:
 *   - Parses the standard Cbc solution header (see CbcSolver.cpp's
 *     `saveSolution`) to determine the claimed status. Importantly, when no
 *     integer-feasible solution was found (e.g. status "Stopped on time
 *     (no integer solution - continuous used)", or "Integer infeasible"),
 *     the listed variable values are the *fractional* LP relaxation, not an
 *     integer-feasible solution — such runs are NOT checked for
 *     integrality/row feasibility, only reported.
 *   - For genuine integer-feasible solutions: checks variable bounds,
 *     integrality, row/constraint activities, and recomputes the objective.
 *   - If the run claims optimality and an expected best-known objective was
 *     supplied (--expected-status optimal --expected-obj <v>), verifies the
 *     claimed value matches within tolerance.
 *   - If the run claims infeasibility, cross-checks against the expected
 *     status (--expected-status infeasible) coming from the trusted
 *     mip-sanity-data `bks.tsv` reference.
 *
 * Usage:
 *   cbc_validate_sol [options] <problem.{mps,lp}[.gz]> <solution.sol>
 *
 * Exit codes:
 *   0  solution/status is consistent with expectations (within tolerances)
 *   1  one or more violations/mismatches found
 *   2  usage / file error
 */

#include "Cbc_C_Interface.h"
#include <cctype>
#include <cerrno>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fcntl.h>
#include <map>
#include <string>
#include <vector>
#include <algorithm>
#include <unistd.h>

/* ── Tolerances ─────────────────────────────────────────────────────── */

/* Relaxed defaults account for limited decimal places / rounding in .sol
 * files, and for cut/presolve post-processing/un-scaling error.           */
static double PRIMAL_TOL = 1e-4; /* bound / constraint feasibility */
static double INT_TOL = 1e-4;    /* integrality gap                 */
static double OBJ_TOL = 1e-4;    /* relative objective discrepancy (sol vs recomputed) */
static double BKS_ABS_TOL = 1e-4;   /* absolute tolerance vs best-known objective */
static double BKS_PCT_TOL = 0.01;   /* percent tolerance vs best-known objective (0.01%) */

/* ── Utilities ──────────────────────────────────────────────────────── */

static bool isNumericStr(const char *s)
{
  if (!s || !*s)
    return false;
  for (const char *p = s; *p; ++p)
    if (!isdigit(*p) && *p != '.' && *p != '-' && *p != '+' && *p != 'e' && *p != 'E')
      return false;
  return true;
}

/* Silence all output from the model-loading call (no need to show Cbc banner). */
static int devNull_ = -1;
static int savedOut_ = -1;
static int savedErr_ = -1;

static void silenceBegin()
{
  devNull_ = open("/dev/null", O_WRONLY);
  if (devNull_ < 0)
    return;
  savedOut_ = dup(STDOUT_FILENO);
  savedErr_ = dup(STDERR_FILENO);
  dup2(devNull_, STDOUT_FILENO);
  dup2(devNull_, STDERR_FILENO);
}

static void silenceEnd()
{
  fflush(stdout);
  fflush(stderr);
  if (savedOut_ >= 0) {
    dup2(savedOut_, STDOUT_FILENO);
    close(savedOut_);
    savedOut_ = -1;
  }
  if (savedErr_ >= 0) {
    dup2(savedErr_, STDERR_FILENO);
    close(savedErr_);
    savedErr_ = -1;
  }
  if (devNull_ >= 0) {
    close(devNull_);
    devNull_ = -1;
  }
}

/* ── Sol file parser ────────────────────────────────────────────────── */

struct SolFile {
  std::string statusLine;   /* full header line, minus trailing " - objective value X" */
  double claimedObj;        /* objective value from header line            */
  bool hasClaimedObj;       /* false if header could not be parsed         */
  /* variable values: col_name -> value (only listed entries stored)      */
  std::map< std::string, double > values;
};

static bool parseSolFile(const char *path, SolFile &out)
{
  FILE *f = fopen(path, "r");
  if (!f) {
    fprintf(stderr, "Cannot open solution file: %s\n", path);
    return false;
  }

  out.claimedObj = 0.0;
  out.hasClaimedObj = false;
  bool firstLine = true;

  char line[1024];
  while (fgets(line, sizeof(line), f)) {
    /* Strip trailing newline */
    size_t len = strlen(line);
    while (len > 0 && (line[len - 1] == '\n' || line[len - 1] == '\r'))
      line[--len] = '\0';

    if (firstLine) {
      firstLine = false;
      /* Cbc's format (CbcSolver.cpp saveSolution):
       *   "<Status text> - objective value <v>"
       * e.g. "Optimal - objective value 3.00000000"
       *      "Integer infeasible - objective value 20470.08333333"
       *      "Stopped on time (no integer solution - continuous used) - objective value 123.0"
       * The status text itself may contain parentheses with more text, so
       * split on the *last* " - objective value" marker.                   */
      const char *marker = strstr(line, " - objective value");
      if (marker) {
        out.statusLine.assign(line, marker - line);
        const char *p = marker + strlen(" - objective value");
        char *end;
        double v = strtod(p, &end);
        if (end != p) {
          out.claimedObj = v;
          out.hasClaimedObj = true;
        }
      } else {
        out.statusLine = line;
      }
      continue;
    }

    /* Variable lines: "idx name value [extra]" or "name value", optionally
     * prefixed with a "**" marker that CbcSolver.cpp writes ahead of the
     * index for rows/columns whose value lies outside its bounds beyond
     * tolerance (see CbcSolver.cpp's use of "** " before printFormat). That
     * marker must be skipped, not parsed as the variable name, or the real
     * name/value one token to the right is silently dropped.               */
    char col[5][256];
    int nread = sscanf(line, "%255s %255s %255s %255s %255s",
      col[0], col[1], col[2], col[3], col[4]);
    if (nread <= 0 || col[0][0] == '#')
      continue;

    int base = 0;
    if (strcmp(col[0], "**") == 0) {
      base = 1;
      --nread;
    }

    const char *name;
    const char *valStr;
    if (isNumericStr(col[base]) && nread >= 3) {
      /* CBC .sol format: index name value */
      name = col[base + 1];
      valStr = col[base + 2];
    } else if (nread >= 2) {
      /* name value format */
      name = col[base];
      valStr = col[base + 1];
    } else {
      continue;
    }

    if (!isNumericStr(valStr))
      continue;

    double val = atof(valStr);
    out.values[name] = val;
  }

  fclose(f);
  return true;
}

/* ── Status classification ─────────────────────────────────────────── */

enum class RunStatus {
  Optimal,          /* proven optimal (possibly within gap tolerance)      */
  Infeasible,       /* proven infeasible (LP-infeasible or bab-infeasible) */
  Unbounded,        /* proven unbounded                                    */
  StoppedWithSol,   /* hit a limit, but an integer-feasible sol is listed  */
  StoppedNoSol,     /* hit a limit, only the fractional LP relaxation is listed */
  Unknown           /* unparseable / unexpected status text                */
};

static RunStatus classifyStatus(const std::string &statusLine)
{
  if (statusLine.rfind("Optimal", 0) == 0)
    return RunStatus::Optimal;
  if (statusLine == "Integer infeasible" || statusLine == "Infeasible")
    return RunStatus::Infeasible;
  if (statusLine == "Unbounded")
    return RunStatus::Unbounded;
  if (statusLine.rfind("Stopped", 0) == 0) {
    bool noIntSol = statusLine.find("no integer solution") != std::string::npos;
    return noIntSol ? RunStatus::StoppedNoSol : RunStatus::StoppedWithSol;
  }
  return RunStatus::Unknown;
}

static const char *statusName(RunStatus s)
{
  switch (s) {
  case RunStatus::Optimal: return "OPTIMAL";
  case RunStatus::Infeasible: return "INFEASIBLE";
  case RunStatus::Unbounded: return "UNBOUNDED";
  case RunStatus::StoppedWithSol: return "STOPPED (integer solution found)";
  case RunStatus::StoppedNoSol: return "STOPPED (no integer solution - fractional relaxation only)";
  default: return "UNKNOWN";
  }
}

/* ── Main ────────────────────────────────────────────────────────────── */

static void printUsage(const char *progName)
{
  fprintf(stderr,
    "Usage: %s [options] <problem.{mps,lp}[.gz]> <solution.sol>\n"
    "\n"
    "Options:\n"
    "  -p, --primal-tol <val>     Primal feasibility tolerance (default: %.1e)\n"
    "  -i, --int-tol <val>        Integrality tolerance (default: %.1e)\n"
    "  -o, --obj-tol <val>        Relative tolerance: recomputed vs claimed obj (default: %.1e)\n"
    "  --expected-status <s>      Cross-check vs mip-sanity-data bks.tsv: 'optimal' or 'infeasible'\n"
    "  --expected-obj <val>       Best-known objective value (used when --expected-status optimal)\n"
    "  --bks-abs-tol <val>        Absolute tolerance vs expected-obj (default: %.1e)\n"
    "  --bks-pct-tol <val>        Percent tolerance vs expected-obj (default: %.4g%%)\n"
    "  -h, --help                 Show this help message\n"
    "\n"
    "Validates feasibility/optimality of a Cbc solution against the problem model\n"
    "and (optionally) a trusted best-known-solution reference.\n"
    "\n"
    "Exit code: 0 = OK,  1 = violations/mismatches found,  2 = file/usage error\n",
    progName, PRIMAL_TOL, INT_TOL, OBJ_TOL, BKS_ABS_TOL, BKS_PCT_TOL);
}

int main(int argc, char *argv[])
{
  std::string expectedStatus; /* "" (none), "optimal", "infeasible" */
  bool haveExpectedObj = false;
  double expectedObj = 0.0;

  int argIdx = 1;
  while (argIdx < argc && argv[argIdx][0] == '-') {
    std::string arg = argv[argIdx];
    if (arg == "-h" || arg == "--help") {
      printUsage(argv[0]);
      return 2;
    } else if (arg == "-p" || arg == "--primal-tol") {
      if (argIdx + 1 >= argc) { fprintf(stderr, "Error: missing value for %s\n", arg.c_str()); return 2; }
      PRIMAL_TOL = atof(argv[++argIdx]);
    } else if (arg == "-i" || arg == "--int-tol") {
      if (argIdx + 1 >= argc) { fprintf(stderr, "Error: missing value for %s\n", arg.c_str()); return 2; }
      INT_TOL = atof(argv[++argIdx]);
    } else if (arg == "-o" || arg == "--obj-tol") {
      if (argIdx + 1 >= argc) { fprintf(stderr, "Error: missing value for %s\n", arg.c_str()); return 2; }
      OBJ_TOL = atof(argv[++argIdx]);
    } else if (arg == "--expected-status") {
      if (argIdx + 1 >= argc) { fprintf(stderr, "Error: missing value for %s\n", arg.c_str()); return 2; }
      expectedStatus = argv[++argIdx];
      std::transform(expectedStatus.begin(), expectedStatus.end(), expectedStatus.begin(), ::tolower);
    } else if (arg == "--expected-obj") {
      if (argIdx + 1 >= argc) { fprintf(stderr, "Error: missing value for %s\n", arg.c_str()); return 2; }
      expectedObj = atof(argv[++argIdx]);
      haveExpectedObj = true;
    } else if (arg == "--bks-abs-tol") {
      if (argIdx + 1 >= argc) { fprintf(stderr, "Error: missing value for %s\n", arg.c_str()); return 2; }
      BKS_ABS_TOL = atof(argv[++argIdx]);
    } else if (arg == "--bks-pct-tol") {
      if (argIdx + 1 >= argc) { fprintf(stderr, "Error: missing value for %s\n", arg.c_str()); return 2; }
      BKS_PCT_TOL = atof(argv[++argIdx]);
    } else {
      fprintf(stderr, "Error: unknown option '%s'\n\n", arg.c_str());
      printUsage(argv[0]);
      return 2;
    }
    ++argIdx;
  }

  if (argc - argIdx < 2) {
    fprintf(stderr, "Error: missing problem file and/or solution file\n\n");
    printUsage(argv[0]);
    return 2;
  }

  const char *problemFile = argv[argIdx];
  const char *solFile = argv[argIdx + 1];

  /* ── Load problem ──────────────────────────────────────────────────── */

  silenceBegin();
  Cbc_Model *m = Cbc_newModel();
  int readErr = Cbc_readMps(m, problemFile);
  if (readErr) {
    silenceEnd();
    silenceBegin();
    readErr = Cbc_readLp(m, problemFile);
  }
  silenceEnd();

  if (readErr) {
    fprintf(stderr, "Error: could not read problem file '%s'\n", problemFile);
    Cbc_deleteModel(m);
    return 2;
  }

  int ncols = Cbc_getNumCols(m);
  int nrows = Cbc_getNumRows(m);

  /* ── Parse solution file ───────────────────────────────────────────── */

  SolFile sol;
  if (!parseSolFile(solFile, sol)) {
    Cbc_deleteModel(m);
    return 2;
  }

  RunStatus status = classifyStatus(sol.statusLine);

  printf("cbc_validate_sol\n");
  printf("  Problem:  %s\n", problemFile);
  printf("  Solution: %s\n\n", solFile);
  printf("  Claimed status: \"%s\"  -> %s\n", sol.statusLine.c_str(), statusName(status));
  if (sol.hasClaimedObj)
    printf("  Claimed objective: %.10g\n", sol.claimedObj);
  printf("\n");

  std::vector< std::string > violations;

  /* ── Infeasible / Unbounded / Unknown claims ───────────────────────── */

  if (status == RunStatus::Infeasible) {
    if (expectedStatus == "optimal") {
      char buf[512];
      snprintf(buf, sizeof(buf),
        "MISMATCH: Cbc claims the problem is infeasible, but the trusted "
        "reference (bks.tsv) says it is optimal (expected obj=%.10g). "
        "This is a solver bug (false infeasibility).", expectedObj);
      violations.push_back(buf);
    } else if (expectedStatus == "infeasible") {
      printf("  Infeasibility claim matches expected status (bks.tsv): OK\n\n");
    } else {
      printf("  NOTE: infeasibility claimed, but no --expected-status was supplied "
             "to cross-check against a trusted reference.\n\n");
    }
  } else if (status == RunStatus::Unbounded) {
    violations.push_back("UNBOUNDED: solver claims the problem is unbounded — "
                          "unexpected for a mip-sanity-data instance (check for a modeling/limit bug).");
  } else if (status == RunStatus::Unknown) {
    char buf[512];
    snprintf(buf, sizeof(buf), "UNKNOWN STATUS: could not classify solution header \"%s\".",
      sol.statusLine.c_str());
    violations.push_back(buf);
  }

  /* ── Feasibility / objective checks (only meaningful when an actual
   *    integer-feasible solution is listed)                             */

  bool feasible = true;
  double computedObj = 0.0;

  if (status == RunStatus::Optimal || status == RunStatus::StoppedWithSol) {
    /* Build full solution vector (unmentioned variables = 0) */
    std::vector< double > x(ncols, 0.0);
    int matched = 0;
    char nameBuf[512];
    for (int i = 0; i < ncols; ++i) {
      Cbc_getColName(m, i, nameBuf, sizeof(nameBuf));
      auto it = sol.values.find(nameBuf);
      if (it != sol.values.end()) {
        x[i] = it->second;
        ++matched;
      }
    }
    int solNonzero = static_cast< int >(sol.values.size());
    int unmatched = solNonzero - matched;

    const double *colLB = Cbc_getColLower(m);
    const double *colUB = Cbc_getColUpper(m);
    const double *objCoef = Cbc_getObjCoefficients(m);
    const double *rowLB = Cbc_getRowLower(m);
    const double *rowUB = Cbc_getRowUpper(m);

    printf("  Model: %d variables (%d integer),  %d constraints\n", ncols, Cbc_getNumIntegers(m), nrows);
    printf("  Solution file: %d listed variable(s)", solNonzero);
    if (unmatched > 0)
      printf("  (WARNING: %d not found in model)", unmatched);
    printf("\n\n");

    /* 1. Column bounds & integrality */
    int nBoundViol = 0, nIntViol = 0;
    for (int i = 0; i < ncols; ++i) {
      double v = x[i];
      double lb = colLB[i];
      double ub = colUB[i];

      double viol = 0.0;
      const char *dir = "";
      double lbTol = (lb > -1e20) ? PRIMAL_TOL * (1.0 + fabs(lb)) : 0.0;
      double ubTol = (ub < 1e20) ? PRIMAL_TOL * (1.0 + fabs(ub)) : 0.0;
      if (lb > -1e20 && v < lb - lbTol) { viol = lb - v; dir = "BELOW LB"; }
      else if (ub < 1e20 && v > ub + ubTol) { viol = v - ub; dir = "ABOVE UB"; }
      if (viol > 0.0) {
        ++nBoundViol;
        Cbc_getColName(m, i, nameBuf, sizeof(nameBuf));
        char buf[512];
        snprintf(buf, sizeof(buf),
          "  BOUND: col %d (%s)  val=%.10g  lb=%.10g  ub=%.10g  viol=%.3e [%s]",
          i, nameBuf, v, lb, ub, viol, dir);
        violations.push_back(buf);
      }

      if (Cbc_isInteger(m, i)) {
        double nearest = round(v);
        double intErr = fabs(v - nearest);
        if (intErr > INT_TOL) {
          ++nIntViol;
          Cbc_getColName(m, i, nameBuf, sizeof(nameBuf));
          char buf[512];
          snprintf(buf, sizeof(buf),
            "  INTEG: col %d (%s)  val=%.10g  nearest_int=%.0f  error=%.3e",
            i, nameBuf, v, nearest, intErr);
          violations.push_back(buf);
        }
      }
    }

    /* 2. Row activities */
    int nRowViol = 0;
    for (int r = 0; r < nrows; ++r) {
      int nz = Cbc_getRowNz(m, r);
      const int *idx = Cbc_getRowIndices(m, r);
      const double *coef = Cbc_getRowCoeffs(m, r);

      double activity = 0.0;
      for (int k = 0; k < nz; ++k)
        activity += coef[k] * x[idx[k]];

      double lb = rowLB[r];
      double ub = rowUB[r];

      double viol = 0.0;
      const char *dir = "";
      double lbTol = (lb > -1e20) ? PRIMAL_TOL * (1.0 + fabs(lb)) : 0.0;
      double ubTol = (ub < 1e20) ? PRIMAL_TOL * (1.0 + fabs(ub)) : 0.0;
      if (lb > -1e20 && activity < lb - lbTol) { viol = lb - activity; dir = "BELOW LB"; }
      else if (ub < 1e20 && activity > ub + ubTol) { viol = activity - ub; dir = "ABOVE UB"; }
      if (viol > 0.0) {
        ++nRowViol;
        char rowName[512];
        Cbc_getRowName(m, r, rowName, sizeof(rowName));
        char buf[512];
        snprintf(buf, sizeof(buf),
          "  ROW:   row %d (%s)  activity=%.10g  bounds=[%.10g, %.10g]  viol=%.3e [%s]",
          r, rowName, activity, lb, ub, viol, dir);
        violations.push_back(buf);
      }
    }

    /* 3. Objective (recomputed vs claimed in .sol) */
    computedObj = 0.0;
    for (int i = 0; i < ncols; ++i)
      computedObj += objCoef[i] * x[i];

    if (sol.hasClaimedObj) {
      double denom = std::max(1.0, fabs(sol.claimedObj));
      double objDelta = fabs(computedObj - sol.claimedObj) / denom;
      if (objDelta > OBJ_TOL) {
        char buf[256];
        snprintf(buf, sizeof(buf),
          "  OBJ:   computed=%.10g  claimed(.sol)=%.10g  rel_delta=%.3e (exceeds tol %.1e)",
          computedObj, sol.claimedObj, objDelta, OBJ_TOL);
        violations.push_back(buf);
      }
    }

    /* 4. Cross-check against expected best-known status/objective */
    if (status == RunStatus::Optimal) {
      if (expectedStatus == "infeasible") {
        violations.push_back(
          "MISMATCH: Cbc proved optimality, but the trusted reference (bks.tsv) "
          "says this instance is infeasible. This is a solver bug.");
      } else if (expectedStatus == "optimal" && haveExpectedObj) {
        double tolAbs = std::max(BKS_ABS_TOL, fabs(expectedObj) * BKS_PCT_TOL / 100.0);
        double delta = fabs(computedObj - expectedObj);
        if (delta > tolAbs) {
          char buf[512];
          snprintf(buf, sizeof(buf),
            "MISMATCH: claimed-optimal objective %.10g does not match best-known "
            "value %.10g (delta=%.3e exceeds tol=%.3e = max(abs=%.1e, pct=%.4g%%))",
            computedObj, expectedObj, delta, tolAbs, BKS_ABS_TOL, BKS_PCT_TOL);
          violations.push_back(buf);
        } else {
          printf("  Optimal objective matches best-known value (delta=%.3e <= tol=%.3e): OK\n\n",
            delta, tolAbs);
        }
      }
    }
  } else if (status == RunStatus::StoppedNoSol) {
    printf("  No integer-feasible solution was found (status indicates the listed\n"
           "  values are the fractional LP relaxation) — skipping integrality/row\n"
           "  feasibility checks for this run, as expected.\n\n");
  }

  feasible = violations.empty();

  /* ── Results ────────────────────────────────────────────────────────── */

  if (feasible) {
    printf("RESULT: OK\n");
  } else {
    printf("RESULT: FAILED\n\n");
    printf("─── VIOLATIONS FOUND (%zu) ─────────────────────────────────────\n", violations.size());
    for (const auto &v : violations)
      printf("%s\n", v.c_str());
    printf("\n");
  }

  Cbc_deleteModel(m);
  return feasible ? 0 : 1;
}
