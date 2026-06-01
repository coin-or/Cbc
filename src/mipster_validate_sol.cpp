/*
 * mipster_validate_sol — solution feasibility checker for MIPster
 *
 * Validates whether a solution file (.sol) is feasible for a given MIP problem.
 * Checks:
 *   - Variable bounds (lb <= x[i] <= ub)
 *   - Integrality (for integer/binary variables)
 *   - Row/constraint activities (rowLB <= Ax <= rowUB)
 *   - Objective value (computed vs. claimed in .sol header)
 *
 * Usage:
 *   mipster_validate_sol <problem.{mps,lp}[.gz]> <solution.sol>
 *
 * Exit codes:
 *   0  solution is feasible (within tolerances)
 *   1  one or more violations found
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
#include <unistd.h>
#include <vector>

/* ── Tolerances ─────────────────────────────────────────────────────── */

static const double PRIMAL_TOL = 1e-7; /* bound / constraint feasibility */
static const double INT_TOL = 1e-5;    /* integrality gap                 */
static const double OBJ_TOL = 1e-6;    /* relative objective discrepancy  */

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
  std::string statusStr;    /* first token on header line, e.g. "Optimal" */
  double claimedObj;        /* objective value from header line            */
  bool hasClaimedObj;       /* false if header could not be parsed         */
  /* variable values: col_name -> value (only nonzero entries stored)     */
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
      /* Expected format: "Optimal - objective value 37.00000000"
       * or just         "Infeasible - objective value ..."
       * The last token is the objective value.                            */
      out.statusStr = "";
      char tok[256];
      if (sscanf(line, "%255s", tok) == 1)
        out.statusStr = tok;

      /* Scan for the objective value — last numeric token */
      const char *p = line;
      double lastVal = 0.0;
      bool found = false;
      while (*p) {
        while (*p && !isdigit(*p) && *p != '-' && *p != '+')
          ++p;
        if (!*p)
          break;
        char *end;
        double v = strtod(p, &end);
        if (end != p) {
          lastVal = v;
          found = true;
          p = end;
        } else {
          ++p;
        }
      }
      if (found) {
        out.claimedObj = lastVal;
        out.hasClaimedObj = true;
      }
      continue;
    }

    /* Variable lines: either "idx name value [extra]" or "name value" */
    char col[4][256];
    int nread = sscanf(line, "%255s %255s %255s %255s",
      col[0], col[1], col[2], col[3]);
    if (nread <= 0 || col[0][0] == '#')
      continue;

    const char *name;
    const char *valStr;
    if (isNumericStr(col[0]) && nread >= 3) {
      /* CBC .sol format: index name value */
      name = col[1];
      valStr = col[2];
    } else if (nread >= 2) {
      /* name value format */
      name = col[0];
      valStr = col[1];
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

/* ── Main ────────────────────────────────────────────────────────────── */

int main(int argc, char *argv[])
{
  if (argc < 3) {
    fprintf(stderr,
      "Usage: mipster_validate_sol <problem.{mps,lp}[.gz]> <solution.sol>\n"
      "\n"
      "Validates feasibility of a MIP solution against the problem model.\n"
      "Checks variable bounds, integrality, constraints, and objective value.\n"
      "\n"
      "Exit code: 0 = feasible,  1 = violations found,  2 = file/usage error\n");
    return 2;
  }

  const char *problemFile = argv[1];
  const char *solFile = argv[2];

  /* ── Load problem ──────────────────────────────────────────────────── */

  silenceBegin();
  Cbc_Model *m = Cbc_newModel();
  int readErr = Cbc_readMps(m, problemFile);
  if (readErr) {
    silenceEnd();
    /* Try LP format */
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

  /* Build full solution vector (unmentioned variables = 0) */
  std::vector< double > x(ncols, 0.0);
  int matched = 0;
  int unmatched = 0;
  char nameBuf[512];
  for (int i = 0; i < ncols; ++i) {
    Cbc_getColName(m, i, nameBuf, sizeof(nameBuf));
    auto it = sol.values.find(nameBuf);
    if (it != sol.values.end()) {
      x[i] = it->second;
      ++matched;
    }
  }
  for (auto &kv : sol.values)
    if (kv.second != 0.0)
      ++unmatched; /* will subtract matched below */
  /* unmatched = variables in sol file that were not found in the model */
  int solNonzero = static_cast< int >(sol.values.size());
  unmatched = solNonzero - matched;

  /* ── Retrieve model data ───────────────────────────────────────────── */

  const double *colLB = Cbc_getColLower(m);
  const double *colUB = Cbc_getColUpper(m);
  const double *objCoef = Cbc_getObjCoefficients(m);
  const double *rowLB = Cbc_getRowLower(m);
  const double *rowUB = Cbc_getRowUpper(m);
  double objSense = Cbc_getObjSense(m); /* +1 = min, -1 = max */

  /* ── Header output ─────────────────────────────────────────────────── */

  printf("mipster_validate_sol\n");
  printf("  Problem:  %s\n", problemFile);
  printf("  Solution: %s\n\n", solFile);
  printf("  Tolerances: primal=%.0e  integer=%.0e  objective=%.0e\n\n",
    PRIMAL_TOL, INT_TOL, OBJ_TOL);
  printf("  Model: %d variables (%d integer),  %d constraints\n",
    ncols, Cbc_getNumIntegers(m), nrows);
  printf("  Solution file: %d non-zero variables listed",
    solNonzero);
  if (unmatched > 0)
    printf("  (WARNING: %d variable(s) in .sol file not found in model)", unmatched);
  printf("\n\n");

  if (sol.hasClaimedObj) {
    printf("  Claimed status: %s   claimed obj = %.10g\n\n",
      sol.statusStr.c_str(), sol.claimedObj);
  }

  /* ── Validation ────────────────────────────────────────────────────── */

  struct Violation {
    std::string what;
  };
  std::vector< Violation > violations;

  /* Track worst-case metrics (even if within tolerance) */
  double worstBoundViol = 0.0;
  int worstBoundCol = -1;
  char worstBoundDir[16] = "";

  double worstIntViol = 0.0;
  int worstIntCol = -1;

  double worstRowViol = 0.0;
  int worstRowIdx = -1;
  char worstRowDir[16] = "";

  /* 1. Column bounds & integrality ──────────────────────────────────── */

  int nBoundViol = 0;
  int nIntViol = 0;

  for (int i = 0; i < ncols; ++i) {
    double v = x[i];
    double lb = colLB[i];
    double ub = colUB[i];

    /* Bound check */
    double viol = 0.0;
    const char *dir = "";
    if (v < lb - PRIMAL_TOL) {
      viol = lb - v;
      dir = "BELOW LB";
    } else if (v > ub + PRIMAL_TOL) {
      viol = v - ub;
      dir = "ABOVE UB";
    }
    if (viol > 0.0) {
      ++nBoundViol;
      Cbc_getColName(m, i, nameBuf, sizeof(nameBuf));
      char buf[512];
      snprintf(buf, sizeof(buf),
        "  BOUND: col %d (%s)  val=%.10g  lb=%.10g  ub=%.10g  viol=%.3e [%s]",
        i, nameBuf, v, lb, ub, viol, dir);
      violations.push_back({ buf });
    }
    /* Track worst (including within-tolerance) */
    double boundErr = std::max(lb - v, v - ub);
    if (boundErr > worstBoundViol) {
      worstBoundViol = boundErr;
      worstBoundCol = i;
      strncpy(worstBoundDir, (v < lb) ? "below lb" : "above ub",
        sizeof(worstBoundDir) - 1);
    }

    /* Integrality check */
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
        violations.push_back({ buf });
      }
      if (intErr > worstIntViol) {
        worstIntViol = intErr;
        worstIntCol = i;
      }
    }
  }

  /* 2. Row activities ─────────────────────────────────────────────────  */

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
    if (activity < lb - PRIMAL_TOL) {
      viol = lb - activity;
      dir = "BELOW LB";
    } else if (activity > ub + PRIMAL_TOL) {
      viol = activity - ub;
      dir = "ABOVE UB";
    }
    if (viol > 0.0) {
      ++nRowViol;
      char rowName[512];
      Cbc_getRowName(m, r, rowName, sizeof(rowName));
      char buf[512];
      snprintf(buf, sizeof(buf),
        "  ROW:   row %d (%s)  activity=%.10g  bounds=[%.10g, %.10g]  viol=%.3e [%s]",
        r, rowName, activity, lb, ub, viol, dir);
      violations.push_back({ buf });
    }
    /* Track worst */
    double rowErr = std::max(lb - activity, activity - ub);
    if (rowErr > worstRowViol) {
      worstRowViol = rowErr;
      worstRowIdx = r;
      strncpy(worstRowDir, (activity < lb) ? "below lb" : "above ub",
        sizeof(worstRowDir) - 1);
    }
  }

  /* 3. Objective ──────────────────────────────────────────────────────── */

  double computedObj = 0.0;
  for (int i = 0; i < ncols; ++i)
    computedObj += objCoef[i] * x[i];

  bool objMismatch = false;
  double objDelta = 0.0;
  if (sol.hasClaimedObj) {
    double denom = std::max(1.0, fabs(sol.claimedObj));
    objDelta = fabs(computedObj - sol.claimedObj) / denom;
    if (objDelta > OBJ_TOL) {
      objMismatch = true;
      char buf[256];
      snprintf(buf, sizeof(buf),
        "  OBJ:   computed=%.10g  claimed=%.10g  rel_delta=%.3e",
        computedObj, sol.claimedObj, objDelta);
      violations.push_back({ buf });
    }
  }

  /* ── Results ────────────────────────────────────────────────────────── */

  bool feasible = violations.empty();

  if (feasible) {
    printf("RESULT: FEASIBLE ✓\n\n");
    printf("  Bounds:       %d variables — no violations\n", ncols);
    printf("  Integrality:  %d integer variables — no violations\n",
      Cbc_getNumIntegers(m));
    printf("  Constraints:  %d rows — no violations\n", nrows);
    if (sol.hasClaimedObj)
      printf("  Objective:    computed=%.10g  claimed=%.10g  rel_delta=%.3e\n\n",
        computedObj, sol.claimedObj, objDelta);
    else
      printf("  Objective:    computed=%.10g\n\n", computedObj);

    printf("─── CLOSEST TO TOLERANCE ────────────────────────────────────────\n");
    if (worstBoundCol >= 0) {
      Cbc_getColName(m, worstBoundCol, nameBuf, sizeof(nameBuf));
      printf("  Largest bound error:        col %d (%s)  error=%.3e [%s]\n",
        worstBoundCol, nameBuf, worstBoundViol, worstBoundDir);
    }
    if (worstIntCol >= 0) {
      Cbc_getColName(m, worstIntCol, nameBuf, sizeof(nameBuf));
      double v = x[worstIntCol];
      printf("  Largest integrality error:  col %d (%s)  val=%.10g  error=%.3e\n",
        worstIntCol, nameBuf, v, worstIntViol);
    }
    if (worstRowIdx >= 0) {
      char rowName[512];
      Cbc_getRowName(m, worstRowIdx, rowName, sizeof(rowName));
      printf("  Largest constraint error:   row %d (%s)  error=%.3e [%s]\n",
        worstRowIdx, rowName, worstRowViol, worstRowDir);
    }
    if (sol.hasClaimedObj)
      printf("  Objective discrepancy:      rel_delta=%.3e\n", objDelta);
    printf("\n");
  } else {
    printf("RESULT: INFEASIBLE ✗\n\n");
    printf("─── VIOLATIONS FOUND (%zu) ─────────────────────────────────────\n",
      violations.size());
    for (const auto &v : violations)
      printf("%s\n", v.what.c_str());
    printf("\n");
    printf("Summary:\n");
    if (nBoundViol > 0)
      printf("  Variable bound violations:   %d\n", nBoundViol);
    if (nIntViol > 0)
      printf("  Integrality violations:      %d\n", nIntViol);
    if (nRowViol > 0)
      printf("  Constraint violations:       %d\n", nRowViol);
    if (objMismatch)
      printf("  Objective mismatch:          computed=%.10g  claimed=%.10g\n",
        computedObj, sol.claimedObj);
    printf("\n");
  }

  Cbc_deleteModel(m);
  return feasible ? 0 : 1;
}
