/*
 * lp-relaxation-test — targeted tests for Cbc_solveLinearProgram()/
 * Cbc_resolve() (the C interface's LP-relaxation-only solve/reoptimize
 * entry points), covering:
 *
 *   1. Root LP relaxation solve with different settings (LPMethod,
 *      DualPivot, LPReductions) on a couple of small MIPLIB instances --
 *      every configuration must agree on the LP relaxation objective
 *      (cross-checked against an independent solver's LP relaxation bound,
 *      hardcoded below as `expectedLpObj`).
 *   2. Reoptimization (Cbc_resolve()) reusing the warm-start basis:
 *        a) fixing a few variables to their optimal LP relaxation values
 *           and reoptimizing must not change the objective (a no-op
 *           tightening);
 *        b) removing a binding constraint, resolving (bound must relax to
 *           a known, independently-computed value), then adding it back
 *           and resolving again (bound must return to the original root
 *           LP value).
 *   3. Solution validation: after solving, the reported primal solution,
 *      row activity/slack, and reduced costs must all be internally
 *      consistent (row activity within row bounds, reduced-cost
 *      complementary slackness sign at variable bounds) -- since a
 *      different LP method may legitimately return a different optimal
 *      basis/vertex, these checks validate feasibility & LP optimality
 *      conditions rather than expecting one specific basis.
 *   4. LP racing (LPM_Racing) and its thread precondition: racing must
 *      reach the same LP optimum as the sequential methods when >= 2
 *      threads are available, and must degrade gracefully (not abort the
 *      process) when asked for with a single thread.
 *   5. Perturbation (INT_PARAM_PERT_VALUE) reaching Clp on *every* solve
 *      path. Clp's own default is 100 ("off") whereas Cbc's tuned default
 *      is the mild 50 that fast warm-started reoptimization relies on, so a
 *      path that forgets to forward this silently reoptimizes with
 *      perturbation off -- a pure performance regression that no
 *      objective-value check can catch. These tests therefore read the
 *      perturbation actually installed on the underlying ClpSimplex.
 *
 * Usage: lp-relaxation-test [--data-dir=PATH]
 * Exit code: 0 = all checks passed, 1 = one or more checks failed.
 */

#include "Cbc_C_Interface.h"
#include "ClpSimplex.hpp"
#include "ClpPEDualRowDantzig.hpp"
#include "ClpPEDualRowSteepest.hpp"
#include "OsiClpSolverInterface.hpp"
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <typeinfo>
#include <unistd.h>
#include <vector>

namespace {

int g_failures = 0;

void check(bool cond, const std::string &msg)
{
  if (!cond) {
    ++g_failures;
    printf("  FAIL: %s\n", msg.c_str());
  }
}

void checkClose(double actual, double expected, double tol, const std::string &msg)
{
  double delta = fabs(actual - expected);
  if (delta > tol) {
    ++g_failures;
    printf("  FAIL: %s (actual=%.10g expected=%.10g delta=%.3e > tol=%.3e)\n",
      msg.c_str(), actual, expected, delta, tol);
  } else {
    printf("  ok:   %s (actual=%.10g expected=%.10g delta=%.3e <= tol=%.3e)\n",
      msg.c_str(), actual, expected, delta, tol);
  }
}

void checkEqInt(int actual, int expected, const std::string &msg)
{
  if (actual != expected) {
    ++g_failures;
    printf("  FAIL: %s (actual=%d expected=%d)\n", msg.c_str(), actual, expected);
  } else {
    printf("  ok:   %s (=%d)\n", msg.c_str(), actual);
  }
}

/* Perturbation actually installed on the ClpSimplex behind the C model --
 * i.e. what the next simplex call will really use, as opposed to what was
 * merely requested via Cbc_setIntParam(). */
int clpPerturbation(Cbc_Model *m)
{
  OsiClpSolverInterface *osi =
    static_cast<OsiClpSolverInterface *>(Cbc_getSolverPtr(m));
  return osi->getModelPtr()->perturbation();
}

/* The two time-limit deadlines actually armed on the ClpSimplex. Each is an
 * absolute deadline (Clp stores limit + now()) with -1.0 as the sole "no limit"
 * sentinel, and hitMaximumIterations() tests both -- so which one is set, not
 * just its value, decides whether INT_PARAM_ELAPSED_TIME is honored. There is
 * no maximumWallSeconds() getter, hence the getDblParam() for that one. */
void clpTimeLimits(Cbc_Model *m, double &cpuDeadline, double &wallDeadline)
{
  ClpSimplex *clps =
    static_cast<OsiClpSolverInterface *>(Cbc_getSolverPtr(m))->getModelPtr();
  cpuDeadline = clps->maximumSeconds();
  wallDeadline = -1.0;
  clps->getDblParam(ClpMaxWallSeconds, wallDeadline);
}

/* Whether the dual pivot rule currently installed on the ClpSimplex is a
 * Positive Edge variant. Note the test is only meaningful in this direction:
 * the PE classes *derive* from the plain ones, so casting to the plain type
 * would match a PE object too and could never distinguish them. */
bool clpDualPivotIsPE(Cbc_Model *m)
{
  OsiClpSolverInterface *osi =
    static_cast<OsiClpSolverInterface *>(Cbc_getSolverPtr(m));
  ClpDualRowPivot *d = osi->getModelPtr()->dualRowPivot();
  return dynamic_cast<ClpPEDualRowSteepest *>(d) != nullptr
      || dynamic_cast<ClpPEDualRowDantzig *>(d) != nullptr;
}

/* Identity of the installed dual pivot rule: its exact dynamic type plus, for
 * the steepest-edge family, its mode. Both matter -- re-asserting psi must
 * neither change the rule's class nor silently drift its mode (e.g. lose
 * DP_Partial's mode 2 to ClpPEDualRowSteepest's own default of 3), and the
 * type name alone would not catch the latter. */
std::string clpDualPivotIdentity(Cbc_Model *m)
{
  OsiClpSolverInterface *osi =
    static_cast<OsiClpSolverInterface *>(Cbc_getSolverPtr(m));
  ClpDualRowPivot *d = osi->getModelPtr()->dualRowPivot();
  std::string id = typeid(*d).name();
  if (ClpDualRowSteepest *s = dynamic_cast<ClpDualRowSteepest *>(d))
    id += " mode=" + std::to_string(s->mode());
  return id;
}

std::string exeDir()
{
  char buf[4096];
  ssize_t n = readlink("/proc/self/exe", buf, sizeof(buf) - 1);
  if (n <= 0)
    return ".";
  buf[n] = '\0';
  std::string path(buf);
  size_t slash = path.find_last_of('/');
  return (slash == std::string::npos) ? "." : path.substr(0, slash);
}

/* ── Basic LP-optimality validation for the current solution ──────────── */
/* Checks feasibility of the primal solution against row/col bounds and
 * complementary-slackness sign consistency of reduced costs -- not tied to
 * any specific basis, so it is valid regardless of which LP method solved
 * the problem. */
void validateLpSolution(Cbc_Model *m, const std::string &label)
{
  int ncols = Cbc_getNumCols(m);
  int nrows = Cbc_getNumRows(m);
  const double *x = Cbc_getColSolution(m);
  const double *colLB = Cbc_getColLower(m);
  const double *colUB = Cbc_getColUpper(m);
  const double *rowAct = Cbc_getRowActivity(m);
  const double *rowSlack = Cbc_getRowSlack(m);
  const double *rc = Cbc_getReducedCost(m);
  const double feasTol = 1e-6;

  /* Column bound feasibility */
  double maxColViol = 0.0;
  for (int j = 0; j < ncols; ++j) {
    double viol = std::max(colLB[j] - x[j], x[j] - colUB[j]);
    maxColViol = std::max(maxColViol, viol);
  }
  check(maxColViol <= 1e-6, label + ": column bounds respected (max viol="
    + std::to_string(maxColViol) + ")");

  /* Row bound feasibility (rowAct/rowLB/rowUB) */
  const double *rowLB = Cbc_getRowLower(m);
  const double *rowUB = Cbc_getRowUpper(m);
  double maxRowViol = 0.0;
  for (int i = 0; i < nrows; ++i) {
    double viol = std::max(rowLB[i] - rowAct[i], rowAct[i] - rowUB[i]);
    maxRowViol = std::max(maxRowViol, viol);
    /* row slack should be consistent with activity vs bounds (Cbc reports
     * slack = rhs - activity for the binding side) */
    (void)rowSlack;
  }
  check(maxRowViol <= 1e-6, label + ": row activity within bounds (max viol="
    + std::to_string(maxRowViol) + ")");

  /* Reduced-cost complementary slackness: for a variable strictly between
   * its bounds, the reduced cost must be ~0 (within tolerance); at a bound,
   * the reduced cost's sign must be consistent with optimality (assuming a
   * minimization sense, which is what mip-sanity-data instances use). */
  int objSense = 1; /* Cbc_getObjSense not queried here; MIPLIB set defaults to min == +1 */
  int rcViolations = 0;
  for (int j = 0; j < ncols; ++j) {
    bool atLower = (x[j] - colLB[j]) <= 1e-6;
    bool atUpper = (colUB[j] - x[j]) <= 1e-6;
    if (!atLower && !atUpper) {
      if (fabs(rc[j]) > 1e-4)
        ++rcViolations;
    } else if (atLower && !atUpper) {
      if (objSense * rc[j] < -1e-4)
        ++rcViolations;
    } else if (atUpper && !atLower) {
      if (objSense * rc[j] > 1e-4)
        ++rcViolations;
    }
  }
  check(rcViolations == 0, label + ": reduced-cost signs consistent with optimality ("
    + std::to_string(rcViolations) + " violation(s))");
}

/* ── Test 1: root LP relaxation with different settings ────────────────── */

void testRootLpMethods(const std::string &dataDir, const std::string &inst, double expectedLpObj)
{
  printf("\n=== Test 1: root LP relaxation, multiple settings — %s (expected obj=%.10g) ===\n",
    inst.c_str(), expectedLpObj);

  struct Cfg {
    const char *label;
    LPMethod method;
    DualPivot pivot;
  };
  const Cfg configs[] = {
    { "auto/auto",       LPM_Auto,      DP_Auto },
    { "dual/dantzig",    LPM_Dual,      DP_Dantzig },
    { "dual/steepest",   LPM_Dual,      DP_Steepest },
    { "dual/partial",    LPM_Dual,      DP_Partial },
    { "dual/PEsteepest", LPM_Dual,      DP_PESteepest },
    { "primal/auto",     LPM_Primal,    DP_Auto },
    { "barrier/auto",    LPM_Barrier,   DP_Auto },
    { "recommend/auto",  LPM_Recommend, DP_Auto },
  };

  const std::string mpsPath = dataDir + "/mips/" + inst + ".mps.gz";
  for (const Cfg &cfg : configs) {
    Cbc_Model *m = Cbc_newModel();
    if (Cbc_readMps(m, mpsPath.c_str()) != 0) {
      check(false, std::string("could not read instance ") + mpsPath);
      Cbc_deleteModel(m);
      continue;
    }
    Cbc_setLogLevel(m, 0);
    Cbc_setLPmethod(m, cfg.method);
    Cbc_setDualPivot(m, cfg.pivot);

    int rc = Cbc_solveLinearProgram(m);
    check(rc == 0, std::string(cfg.label) + ": Cbc_solveLinearProgram returned optimal (rc=0)");
    if (rc == 0) {
      checkClose(Cbc_getObjValue(m), expectedLpObj, 1e-6,
        std::string(cfg.label) + ": LP relaxation objective matches reference");
      validateLpSolution(m, cfg.label);
    }
    Cbc_deleteModel(m);
  }

  /* LPR_NoDualReds / LPM_BarrierNoCross exercise the legacy cold-solve path */
  {
    Cbc_Model *m = Cbc_newModel();
    Cbc_readMps(m, mpsPath.c_str());
    Cbc_setLogLevel(m, 0);
    Cbc_setDualReductionsType(m, LPR_NoDualReds);
    int rc = Cbc_solveLinearProgram(m);
    check(rc == 0, "LPR_NoDualReds (legacy path): Cbc_solveLinearProgram returned optimal");
    if (rc == 0)
      checkClose(Cbc_getObjValue(m), expectedLpObj, 1e-6,
        "LPR_NoDualReds (legacy path): LP relaxation objective matches reference");
    Cbc_deleteModel(m);
  }
  {
    Cbc_Model *m = Cbc_newModel();
    Cbc_readMps(m, mpsPath.c_str());
    Cbc_setLogLevel(m, 0);
    Cbc_setLPmethod(m, LPM_BarrierNoCross);
    int rc = Cbc_solveLinearProgram(m);
    check(rc == 0, "LPM_BarrierNoCross (legacy path): Cbc_solveLinearProgram returned optimal");
    if (rc == 0)
      checkClose(Cbc_getObjValue(m), expectedLpObj, 1e-6,
        "LPM_BarrierNoCross (legacy path): LP relaxation objective matches reference");
    Cbc_deleteModel(m);
  }
}

/* ── Test 2a: reoptimization via variable fixing (warm start) ──────────── */

void testResolveAfterFixing(const std::string &dataDir, const std::string &inst, double expectedLpObj)
{
  printf("\n=== Test 2a: reoptimization after fixing variables — %s ===\n", inst.c_str());

  const std::string mpsPath = dataDir + "/mips/" + inst + ".mps.gz";
  Cbc_Model *m = Cbc_newModel();
  if (Cbc_readMps(m, mpsPath.c_str()) != 0) {
    check(false, std::string("could not read instance ") + mpsPath);
    Cbc_deleteModel(m);
    return;
  }
  Cbc_setLogLevel(m, 0);

  int rc = Cbc_solveLinearProgram(m);
  check(rc == 0, "initial root LP solve is optimal");
  if (rc != 0) {
    Cbc_deleteModel(m);
    return;
  }
  checkClose(Cbc_getObjValue(m), expectedLpObj, 1e-6, "root LP objective matches reference");

  /* Fixing a handful of variables to their current (optimal) LP relaxation
   * values is a no-op tightening: it must not change the objective, and
   * the reoptimization must reuse the existing basis (Cbc_resolve(), not a
   * cold solve). */
  int ncols = Cbc_getNumCols(m);
  std::vector<double> xBefore(Cbc_getColSolution(m), Cbc_getColSolution(m) + ncols);
  int nFix = std::min(ncols, 5);
  for (int j = 0; j < nFix; ++j) {
    Cbc_setColLower(m, j, xBefore[j]);
    Cbc_setColUpper(m, j, xBefore[j]);
  }

  int rc2 = Cbc_resolve(m);
  check(rc2 == 0, "Cbc_resolve() after fixing variables to their optimal values is optimal");
  if (rc2 == 0) {
    checkClose(Cbc_getObjValue(m), expectedLpObj, 1e-6,
      "objective unchanged after fixing variables to their already-optimal values");
    validateLpSolution(m, "after fixing variables");
  }

  Cbc_deleteModel(m);
}

/* ── Test 2b: reoptimization via removing/re-adding a binding constraint ── */

void testResolveAfterConstraintChange(const std::string &dataDir, const std::string &inst,
  double expectedLpObj, const std::string &rowName, double expectedLpObjWithoutRow)
{
  printf("\n=== Test 2b: reoptimization after removing/re-adding a binding row — %s (row %s) ===\n",
    inst.c_str(), rowName.c_str());

  const std::string mpsPath = dataDir + "/mips/" + inst + ".mps.gz";
  Cbc_Model *m = Cbc_newModel();
  Cbc_storeNameIndexes(m, 1);
  if (Cbc_readMps(m, mpsPath.c_str()) != 0) {
    check(false, std::string("could not read instance ") + mpsPath);
    Cbc_deleteModel(m);
    return;
  }
  Cbc_setLogLevel(m, 0);

  int rc = Cbc_solveLinearProgram(m);
  check(rc == 0, "initial root LP solve (with the binding row present) is optimal");
  if (rc != 0) {
    Cbc_deleteModel(m);
    return;
  }
  checkClose(Cbc_getObjValue(m), expectedLpObj, 1e-6,
    "root LP objective (row present) matches reference");

  int rowIdx = Cbc_getRowNameIndex(m, rowName.c_str());
  check(rowIdx >= 0, "target row found by name");
  if (rowIdx < 0) {
    Cbc_deleteModel(m);
    return;
  }

  /* Save the row's data so it can be re-added identically afterwards. */
  int nz = Cbc_getRowNz(m, rowIdx);
  std::vector<int> cols(Cbc_getRowIndices(m, rowIdx), Cbc_getRowIndices(m, rowIdx) + nz);
  std::vector<double> coefs(Cbc_getRowCoeffs(m, rowIdx), Cbc_getRowCoeffs(m, rowIdx) + nz);
  char sense = Cbc_getRowSense(m, rowIdx);
  double rhs = Cbc_getRowRHS(m, rowIdx);
  char rowNameBuf[512];
  Cbc_getRowName(m, rowIdx, rowNameBuf, sizeof(rowNameBuf));

  /* Remove the row and reoptimize (warm start): the bound must relax to
   * (weaken to) the independently-computed value for "row removed". */
  int rowsToDelete[1] = { rowIdx };
  Cbc_deleteRows(m, 1, rowsToDelete);
  int rc2 = Cbc_resolve(m);
  check(rc2 == 0, "Cbc_resolve() after removing the binding row is optimal");
  if (rc2 == 0) {
    checkClose(Cbc_getObjValue(m), expectedLpObjWithoutRow, 1e-6,
      "LP objective after removing the row matches independently-computed reference");
    validateLpSolution(m, "after removing row");
  }

  /* Add the row back and reoptimize: the bound must return to the original
   * root LP value. */
  Cbc_addRow(m, rowNameBuf, nz, cols.data(), coefs.data(), sense, rhs);
  int rc3 = Cbc_resolve(m);
  check(rc3 == 0, "Cbc_resolve() after re-adding the row is optimal");
  if (rc3 == 0) {
    checkClose(Cbc_getObjValue(m), expectedLpObj, 1e-6,
      "LP objective after re-adding the row returns to the original reference");
    validateLpSolution(m, "after re-adding row");
  }

  Cbc_deleteModel(m);
}

/* ── Test 3: LP racing and its thread precondition ─────────────────────── */

/* LPM_Racing runs several LP configurations (dual, primal+Idiot,
 * primal+Sprint) in parallel threads and takes the first to reach
 * optimality. Two things are checked:
 *   a) with >= 2 threads it must reach the same LP optimum as every
 *      sequential method (the winner is a genuine LP optimum, not a
 *      truncated/loosened one);
 *   b) with a single thread it must NOT abort. Racing with one thread is a
 *      hard error deep inside CbcSolver::applyLpMethod() (CoinError), and an
 *      exception escaping through the extern "C" boundary would terminate
 *      the calling process -- for an embedding interpreter (python-mip) that
 *      turns a mis-set parameter into a crash. The C interface resolves the
 *      combination up front instead, falling back to the sequential
 *      recommended method.
 * LPM_Auto is covered for both thread counts too, since it dispatches to
 * racing / the sequential recommendation based on the same thread count. */
void testRacingLp(const std::string &dataDir, const std::string &inst, double expectedLpObj)
{
  printf("\n=== Test 3: LP racing / thread precondition — %s ===\n", inst.c_str());

  struct Cfg {
    const char *label;
    LPMethod method;
    int threads;
  };
  const Cfg configs[] = {
    { "racing, 4 threads",             LPM_Racing, 4 },
    { "racing, 2 threads",             LPM_Racing, 2 },
    { "racing, 1 thread (must not abort)", LPM_Racing, 1 },
    { "racing, threads unset",         LPM_Racing, -1 },
    { "auto, 4 threads",               LPM_Auto,   4 },
    { "auto, 1 thread",                LPM_Auto,   1 },
  };

  const std::string mpsPath = dataDir + "/mips/" + inst + ".mps.gz";
  for (const Cfg &cfg : configs) {
    Cbc_Model *m = Cbc_newModel();
    if (Cbc_readMps(m, mpsPath.c_str()) != 0) {
      check(false, std::string("could not read instance ") + mpsPath);
      Cbc_deleteModel(m);
      continue;
    }
    Cbc_setLogLevel(m, 0);
    Cbc_setLPmethod(m, cfg.method);
    if (cfg.threads > 0)
      Cbc_setIntParam(m, INT_PARAM_THREADS, cfg.threads);

    int rc = Cbc_solveLinearProgram(m);
    check(rc == 0, std::string(cfg.label) + ": Cbc_solveLinearProgram returned optimal (rc=0)");
    if (rc == 0) {
      checkClose(Cbc_getObjValue(m), expectedLpObj, 1e-6,
        std::string(cfg.label) + ": LP relaxation objective matches reference");
      validateLpSolution(m, cfg.label);
    }
    Cbc_deleteModel(m);
  }

  /* Racing leaves the winning thread's basis in the caller's solver, so a
   * subsequent warm-started reoptimization must work off it as usual. */
  {
    Cbc_Model *m = Cbc_newModel();
    Cbc_readMps(m, mpsPath.c_str());
    Cbc_setLogLevel(m, 0);
    Cbc_setLPmethod(m, LPM_Racing);
    Cbc_setIntParam(m, INT_PARAM_THREADS, 4);
    if (Cbc_solveLinearProgram(m) == 0) {
      int ncols = Cbc_getNumCols(m);
      std::vector<double> x(Cbc_getColSolution(m), Cbc_getColSolution(m) + ncols);
      int nFix = std::min(ncols, 5);
      for (int j = 0; j < nFix; ++j) {
        Cbc_setColLower(m, j, x[j]);
        Cbc_setColUpper(m, j, x[j]);
      }
      int rc = Cbc_resolve(m);
      check(rc == 0, "racing: Cbc_resolve() after racing solve is optimal");
      if (rc == 0)
        checkClose(Cbc_getObjValue(m), expectedLpObj, 1e-6,
          "racing: objective unchanged when reoptimizing after racing (basis reused)");
    }
    Cbc_deleteModel(m);
  }
}

/* ── Test 4: perturbation reaches Clp on every solve path ──────────────── */

/* Regression test for a silent LP-performance bug: INT_PARAM_PERT_VALUE was
 * only forwarded to Clp on one of the three solve paths (as a "-pertV" token
 * on the CbcSolver-driven cold path), so
 *   - a reoptimize following a *legacy* cold solve (LPR_NoDualReds /
 *     LPM_BarrierNoCross, neither of which touches perturbation) ran with
 *     perturbation at Clp's default 100 = OFF, never getting Cbc's tuned 50;
 *   - Cbc_setIntParam(INT_PARAM_PERT_VALUE) issued after the first solve was
 *     silently ignored by every later Cbc_resolve().
 * Both are invisible to objective checks -- the LP still solves to the same
 * optimum, just slower -- so this asserts on the perturbation value actually
 * installed on the ClpSimplex. */
void testPerturbationPaths(const std::string &dataDir, const std::string &inst)
{
  printf("\n=== Test 4: INT_PARAM_PERT_VALUE honored on all solve paths — %s ===\n",
    inst.c_str());

  const std::string mpsPath = dataDir + "/mips/" + inst + ".mps.gz";

  /* Default: Cbc's tuned 50 (mild perturbation), not Clp's 100 (off). */
  {
    Cbc_Model *m = Cbc_newModel();
    Cbc_readMps(m, mpsPath.c_str());
    Cbc_setLogLevel(m, 0);
    checkEqInt(Cbc_getIntParam(m, INT_PARAM_PERT_VALUE), 50,
      "default INT_PARAM_PERT_VALUE is Cbc's tuned 50 (not Clp's 100 = off)");
    Cbc_solveLinearProgram(m);
    checkEqInt(clpPerturbation(m), 50,
      "default: perturbation installed on ClpSimplex after CbcSolver cold path");
    Cbc_deleteModel(m);
  }

  /* Explicit value, on each of the three solve paths. 61 is arbitrary but
   * distinct from both 50 and 100, so neither default can mask a failure. */
  struct Path {
    const char *label;
    LPMethod method;
    int noDualReds;
  };
  const Path paths[] = {
    { "CbcSolver cold path (auto)",         LPM_Auto,             0 },
    { "legacy cold path (LPR_NoDualReds)",  LPM_Auto,             1 },
    { "legacy cold path (BarrierNoCross)",  LPM_BarrierNoCross,   0 },
  };
  for (const Path &p : paths) {
    Cbc_Model *m = Cbc_newModel();
    Cbc_readMps(m, mpsPath.c_str());
    Cbc_setLogLevel(m, 0);
    Cbc_setIntParam(m, INT_PARAM_PERT_VALUE, 61);
    Cbc_setLPmethod(m, p.method);
    if (p.noDualReds)
      Cbc_setDualReductionsType(m, LPR_NoDualReds);
    Cbc_solveLinearProgram(m);
    checkEqInt(clpPerturbation(m), 61,
      std::string(p.label) + ": requested pertV=61 reached ClpSimplex");
    Cbc_deleteModel(m);
  }

  /* After a legacy cold solve, a warm reoptimize must still be running with
   * the configured perturbation rather than Clp's "off" default. */
  {
    Cbc_Model *m = Cbc_newModel();
    Cbc_readMps(m, mpsPath.c_str());
    Cbc_setLogLevel(m, 0);
    Cbc_setDualReductionsType(m, LPR_NoDualReds);
    if (Cbc_solveLinearProgram(m) == 0) {
      int ncols = Cbc_getNumCols(m);
      std::vector<double> x(Cbc_getColSolution(m), Cbc_getColSolution(m) + ncols);
      int nFix = std::min(ncols, 5);
      for (int j = 0; j < nFix; ++j) {
        Cbc_setColLower(m, j, x[j]);
        Cbc_setColUpper(m, j, x[j]);
      }
      Cbc_resolve(m);
      checkEqInt(clpPerturbation(m), 50,
        "reoptimize after legacy cold solve uses Cbc's tuned 50, not Clp's 100 (off)");
    }
    Cbc_deleteModel(m);
  }

  /* A perturbation change made *after* the first solve must take effect on
   * the next Cbc_resolve(). */
  {
    Cbc_Model *m = Cbc_newModel();
    Cbc_readMps(m, mpsPath.c_str());
    Cbc_setLogLevel(m, 0);
    if (Cbc_solveLinearProgram(m) == 0) {
      Cbc_setIntParam(m, INT_PARAM_PERT_VALUE, 97);
      int ncols = Cbc_getNumCols(m);
      std::vector<double> x(Cbc_getColSolution(m), Cbc_getColSolution(m) + ncols);
      int nFix = std::min(ncols, 5);
      for (int j = 0; j < nFix; ++j) {
        Cbc_setColLower(m, j, x[j]);
        Cbc_setColUpper(m, j, x[j]);
      }
      Cbc_resolve(m);
      checkEqInt(clpPerturbation(m), 97,
        "INT_PARAM_PERT_VALUE set after the first solve takes effect on Cbc_resolve()");
    }
    Cbc_deleteModel(m);
  }
}

/* ── Test 5: Positive Edge pivoting survives reoptimization ─────────────── */

/* Companion to Test 4, for the other persistent LP setting that has to reach
 * every solve rather than just the root one: DBL_PARAM_PSI.
 *
 * A positive psi asks for Positive Edge pivoting, which Cbc implements by
 * *wrapping* the installed pivot rules in their PE variants. The cold path
 * requests it via a "-psi" token, but the pivot rule is (re)installed before
 * every solve, including each Cbc_resolve() -- so re-installing a plain
 * ClpDualRowSteepest there used to strip the PE wrapper off the dual side on
 * the first reoptimization while leaving the primal side wrapped. The
 * requested pivot rule then applied to the root LP only, asymmetrically, and
 * to none of the reoptimizations.
 *
 * Like Test 4 this is invisible to objective checks -- the LP still reaches
 * the same optimum, just with more degenerate pivots -- so it asserts on the
 * pivot rule actually installed on the ClpSimplex.
 *
 * Note psi is deliberately *not* treated like the one-shot overrides CbcSolver
 * rolls back via ClpProtectedLpSettings (auto-mode pivot/scaling/perturbation
 * choices that exist only to speed up a single root solve): it is a persistent
 * user-requested setting, so the correct behavior is to re-assert it. */
void testPositiveEdgePersistence(const std::string &dataDir,
  const std::string &inst, double expectedLpObj)
{
  printf("\n=== Test 5: DBL_PARAM_PSI (Positive Edge) survives reoptimization — %s ===\n",
    inst.c_str());

  const std::string mpsPath = dataDir + "/mips/" + inst + ".mps.gz";

  /* Without psi, no PE wrapper should appear anywhere -- guards against the
   * fix over-applying and turning PE on when it was never requested. */
  {
    Cbc_Model *m = Cbc_newModel();
    Cbc_readMps(m, mpsPath.c_str());
    Cbc_setLogLevel(m, 0);
    if (Cbc_solveLinearProgram(m) == 0) {
      check(!clpDualPivotIsPE(m), "psi unset: cold solve installs a plain (non-PE) dual pivot");
      Cbc_resolve(m);
      check(!clpDualPivotIsPE(m), "psi unset: reoptimize keeps a plain (non-PE) dual pivot");
    }
    Cbc_deleteModel(m);
  }

  /* With psi > 0, the PE dual pivot must be in place after the cold solve and
   * still in place after reoptimizing -- on both the CbcSolver-driven cold
   * path and the legacy one, and with DP_Auto (which does not itself name a
   * PE rule) as well as an explicitly selected DP_PESteepest. */
  struct Case {
    const char *label;
    LPMethod method;
    DualPivot pivot;
  };
  const Case cases[] = {
    { "CbcSolver cold path, DP_Auto",       LPM_Auto,           DP_Auto },
    { "CbcSolver cold path, DP_Dual",       LPM_Dual,           DP_Auto },
    { "CbcSolver cold path, DP_Steepest",   LPM_Dual,           DP_Steepest },
    { "CbcSolver cold path, DP_Partial",    LPM_Dual,           DP_Partial },
    { "CbcSolver cold path, DP_Dantzig",    LPM_Dual,           DP_Dantzig },
    { "CbcSolver cold path, DP_PESteepest", LPM_Dual,           DP_PESteepest },
    { "legacy cold path (BarrierNoCross)",  LPM_BarrierNoCross, DP_Auto },
  };
  for (const Case &c : cases) {
    Cbc_Model *m = Cbc_newModel();
    Cbc_readMps(m, mpsPath.c_str());
    Cbc_setLogLevel(m, 0);
    Cbc_setDblParam(m, DBL_PARAM_PSI, 0.5);
    Cbc_setLPmethod(m, c.method);
    Cbc_setDualPivot(m, c.pivot);
    if (Cbc_solveLinearProgram(m) != 0) {
      Cbc_deleteModel(m);
      continue;
    }
    check(clpDualPivotIsPE(m),
      std::string(c.label) + ": psi=0.5 installs a PE dual pivot for the cold solve");
    const std::string coldPivot = clpDualPivotIdentity(m);

    /* Reoptimize repeatedly, after fixing a few variables at their LP values
     * so each resolve has real work to do. The pivot rule must stay PE, and
     * stay the *same* PE variant: re-asserting psi must not drift the rule to
     * a different mode (e.g. losing DP_Partial's mode 2 to the PE class's own
     * default of 3) on each pass through. */
    int ncols = Cbc_getNumCols(m);
    const double *sol = Cbc_getColSolution(m);
    std::vector<double> x(sol, sol + ncols);
    for (int k = 0; k < 3; ++k) {
      int lo = k * 2, hi = std::min(ncols, lo + 2);
      for (int j = lo; j < hi; ++j) {
        Cbc_setColLower(m, j, x[j]);
        Cbc_setColUpper(m, j, x[j]);
      }
      Cbc_resolve(m);
      check(clpDualPivotIsPE(m),
        std::string(c.label) + ": PE dual pivot still installed after reoptimize");
      check(clpDualPivotIdentity(m) == coldPivot,
        std::string(c.label) + ": PE dual pivot rule and mode unchanged by reoptimize");
    }
    Cbc_deleteModel(m);
  }

  /* PE pivoting must not change the optimum it converges to -- it only avoids
   * degenerate moves. Solved fresh (no fixed bounds) so the reference
   * objective applies. */
  {
    Cbc_Model *m = Cbc_newModel();
    Cbc_readMps(m, mpsPath.c_str());
    Cbc_setLogLevel(m, 0);
    Cbc_setDblParam(m, DBL_PARAM_PSI, 0.5);
    if (Cbc_solveLinearProgram(m) == 0) {
      checkClose(Cbc_getObjValue(m), expectedLpObj, 1e-6,
        "psi=0.5 reaches the same LP optimum");
      Cbc_resolve(m);
      checkClose(Cbc_getObjValue(m), expectedLpObj, 1e-6,
        "psi=0.5 reoptimize reaches the same LP optimum");
      validateLpSolution(m, "psi=0.5 reoptimize");
    }
    Cbc_deleteModel(m);
  }

  /* A psi set *after* the first solve must take effect on the next resolve,
   * mirroring the same guarantee Test 4 checks for INT_PARAM_PERT_VALUE. */
  {
    Cbc_Model *m = Cbc_newModel();
    Cbc_readMps(m, mpsPath.c_str());
    Cbc_setLogLevel(m, 0);
    if (Cbc_solveLinearProgram(m) == 0) {
      check(!clpDualPivotIsPE(m), "psi unset before first solve: no PE dual pivot");
      Cbc_setDblParam(m, DBL_PARAM_PSI, 0.5);
      Cbc_resolve(m);
      check(clpDualPivotIsPE(m),
        "DBL_PARAM_PSI set after the first solve takes effect on Cbc_resolve()");
    }
    Cbc_deleteModel(m);
  }
}

/* ── Test 6: LP time limit uses the clock INT_PARAM_ELAPSED_TIME selects ── */

/* Two failure modes that no objective check can catch, since a correctly
 * solved LP looks identical either way:
 *
 *  - Clp keeps CPU and wall-clock deadlines as two independent parameters and
 *    checks both, so arming only the CPU one makes INT_PARAM_ELAPSED_TIME (which
 *    *defaults* to 1 = wall clock) a no-op for LP solves. Under LP racing that
 *    matters materially: CPU time accrues once per racing thread, so an N-second
 *    budget would be enforced after roughly N/threads seconds of wall time.
 *
 *  - COIN_DBL_MAX is the C interface's "no limit", but it is not Clp's:
 *    setMaximumSeconds() stores value + now() for any value >= 0, so passing it
 *    through arms a finite deadline rather than disarming the limit. -1.0 is
 *    Clp's only "no limit" sentinel. */
void testLpTimeLimitClock(const std::string &dataDir, const std::string &inst)
{
  printf("\n=== Test 6: LP time limit honors INT_PARAM_ELAPSED_TIME — %s ===\n",
    inst.c_str());

  const std::string mpsPath = dataDir + "/mips/" + inst + ".mps.gz";
  double cpu = 0.0, wall = 0.0;

  struct Path {
    const char *label;
    LPMethod method;
    int noDualReds;
  };
  const Path paths[] = {
    { "CbcSolver cold path (auto)",         LPM_Auto,             0 },
    { "legacy cold path (LPR_NoDualReds)",  LPM_Auto,             1 },
    { "legacy cold path (BarrierNoCross)",  LPM_BarrierNoCross,   0 },
  };

  /* Default (no limit set) must leave *both* deadlines disarmed, after the cold
   * solve and after a reoptimize. */
  for (const Path &p : paths) {
    Cbc_Model *m = Cbc_newModel();
    Cbc_readMps(m, mpsPath.c_str());
    Cbc_setLogLevel(m, 0);
    Cbc_setLPmethod(m, p.method);
    if (p.noDualReds)
      Cbc_setDualReductionsType(m, LPR_NoDualReds);
    const std::string lbl = std::string(p.label) + ": default time limit";
    if (Cbc_solveLinearProgram(m) == 0) {
      clpTimeLimits(m, cpu, wall);
      check(cpu < 0.0 && wall < 0.0, lbl + " leaves both LP deadlines disarmed after cold solve");
      Cbc_resolve(m);
      clpTimeLimits(m, cpu, wall);
      check(cpu < 0.0 && wall < 0.0, lbl + " leaves both LP deadlines disarmed after reoptimize");
    }
    Cbc_deleteModel(m);
  }

  /* An explicit limit must arm exactly the selected clock and disarm the other,
   * on every path and on the reoptimize too. The deadlines are absolute, so the
   * check is "armed / not armed" rather than an equality against 30. */
  for (int elapsed = 1; elapsed >= 0; --elapsed) {
    for (const Path &p : paths) {
      Cbc_Model *m = Cbc_newModel();
      Cbc_readMps(m, mpsPath.c_str());
      Cbc_setLogLevel(m, 0);
      Cbc_setLPmethod(m, p.method);
      if (p.noDualReds)
        Cbc_setDualReductionsType(m, LPR_NoDualReds);
      Cbc_setIntParam(m, INT_PARAM_ELAPSED_TIME, elapsed);
      Cbc_setDblParam(m, DBL_PARAM_TIME_LIMIT, 30.0);
      const std::string lbl = std::string(p.label) + ", ELAPSED_TIME=" +
        std::to_string(elapsed) + " (" + (elapsed ? "wall" : "CPU") + ")";
      if (Cbc_solveLinearProgram(m) == 0) {
        /* Only the *unselected* clock can be asserted about after the cold
         * solve: the CbcSolver path ends with clearClpTimeLimits(), which
         * deliberately disarms both (its own -sec token carried the limit
         * through the solve, routed via applyClpTimeLimit() and so likewise
         * honoring the selected clock), whereas the legacy paths leave the
         * deadline armed. Both agree that the wrong clock is never armed. */
        clpTimeLimits(m, cpu, wall);
        check(elapsed ? cpu < 0.0 : wall < 0.0,
          lbl + ": cold solve never arms the unselected clock");
        Cbc_resolve(m);
        clpTimeLimits(m, cpu, wall);
        check(elapsed ? (wall >= 0.0 && cpu < 0.0) : (cpu >= 0.0 && wall < 0.0),
          lbl + ": reoptimize arms only the selected clock");
      }
      Cbc_deleteModel(m);
    }
  }

  /* Switching clocks between two solves must not leave the previous solve's
   * deadline armed -- otherwise a stale limit keeps firing on the wrong clock. */
  {
    Cbc_Model *m = Cbc_newModel();
    Cbc_readMps(m, mpsPath.c_str());
    Cbc_setLogLevel(m, 0);
    Cbc_setIntParam(m, INT_PARAM_ELAPSED_TIME, 0);
    Cbc_setDblParam(m, DBL_PARAM_TIME_LIMIT, 30.0);
    if (Cbc_solveLinearProgram(m) == 0) {
      Cbc_setIntParam(m, INT_PARAM_ELAPSED_TIME, 1);
      Cbc_resolve(m);
      clpTimeLimits(m, cpu, wall);
      check(wall >= 0.0 && cpu < 0.0,
        "switching ELAPSED_TIME between solves disarms the previously armed clock");
    }
    Cbc_deleteModel(m);
  }

  /* Clearing the limit back to "no limit" must disarm, not re-arm. */
  {
    Cbc_Model *m = Cbc_newModel();
    Cbc_readMps(m, mpsPath.c_str());
    Cbc_setLogLevel(m, 0);
    Cbc_setDblParam(m, DBL_PARAM_TIME_LIMIT, 30.0);
    if (Cbc_solveLinearProgram(m) == 0) {
      Cbc_setDblParam(m, DBL_PARAM_TIME_LIMIT, COIN_DBL_MAX);
      Cbc_resolve(m);
      clpTimeLimits(m, cpu, wall);
      check(cpu < 0.0 && wall < 0.0,
        "resetting DBL_PARAM_TIME_LIMIT to COIN_DBL_MAX disarms both LP deadlines");
    }
    Cbc_deleteModel(m);
  }

  /* An iteration limit must be reported as truncated (1), not optimal (0), and
   * must be visible through Osi's own predicate. */
  {
    Cbc_Model *m = Cbc_newModel();
    Cbc_readMps(m, mpsPath.c_str());
    Cbc_setLogLevel(m, 0);
    Cbc_setIntParam(m, INT_PARAM_MAX_ITER, 3);
    checkEqInt(Cbc_solveLinearProgram(m), 1,
      "INT_PARAM_MAX_ITER=3 reports status 1 (truncated), not 0 (optimal)");
    OsiClpSolverInterface *osi =
      static_cast<OsiClpSolverInterface *>(Cbc_getSolverPtr(m));
    check(osi->isIterationLimitReached(),
      "INT_PARAM_MAX_ITER=3: Osi reports the iteration limit as reached");
    check(!osi->isProvenOptimal(),
      "INT_PARAM_MAX_ITER=3: LP is not claimed proven optimal");
    Cbc_deleteModel(m);
  }
}

} // namespace

int main(int argc, char *argv[])
{
  std::string dataDir;
  for (int i = 1; i < argc; ++i) {
    std::string arg = argv[i];
    if (arg == "--help" || arg == "-h") {
      printf(
        "Usage: lp-relaxation-test [--data-dir=PATH]\n"
        "\n"
        "Targeted tests for the C interface's LP-relaxation-only entry points,\n"
        "Cbc_solveLinearProgram() and Cbc_resolve(). Covers:\n"
        "  1. root LP solves under every LPMethod/DualPivot/LPReductions setting,\n"
        "     each checked against an independently computed reference objective;\n"
        "  2. warm-start reoptimization (fixing variables; removing and re-adding\n"
        "     a binding row);\n"
        "  3. LP solution validation (bound feasibility, reduced-cost\n"
        "     complementary slackness);\n"
        "  4. LP racing (LPM_Racing) and its >= 2 threads precondition;\n"
        "  5. INT_PARAM_PERT_VALUE reaching Clp on every solve path, including\n"
        "     warm-started reoptimization;\n"
        "  6. DBL_PARAM_PSI (Positive Edge pivoting) likewise surviving into\n"
        "     reoptimization, without the wrapper nesting on repeated solves;\n"
        "  7. the LP time limit being armed on the clock INT_PARAM_ELAPSED_TIME\n"
        "     selects (and only that one), \"no limit\" really disarming it, and\n"
        "     INT_PARAM_MAX_ITER reporting truncation rather than optimality.\n"
        "\n"
        "  --data-dir=PATH  mip-sanity-data checkout to read instances from\n"
        "                   (default: <dir of this executable>/mip-sanity-data)\n"
        "\n"
        "Exit code: 0 = all checks passed, 1 = one or more checks failed.\n");
      return 0;
    }
    if (arg.rfind("--data-dir=", 0) == 0) {
      dataDir = arg.substr(strlen("--data-dir="));
      continue;
    }
    fprintf(stderr, "lp-relaxation-test: unrecognized argument '%s' "
      "(try --help)\n", arg.c_str());
    return 1;
  }
  if (dataDir.empty())
    dataDir = exeDir() + "/mip-sanity-data";

  /* Reference LP relaxation objectives below were computed independently
   * (two different external LP/MIP solvers agreed exactly) and cross-checked
   * against Cbc's own LP relaxation solve before writing this test. All four
   * were re-confirmed with HiGHS 1.13.1 (`solve_relaxation = true`,
   * `presolve = off`; the "row removed" values against a copy of the MPS with
   * that row and all its coefficients stripped):
   *   p0033 2520.5717391304  p0033 less R127 2035.9932608696
   *   lseu   834.6823529412  lseu  less R102  662.9019607843 */
  testRootLpMethods(dataDir, "p0033", 2520.5717391304347);
  testRootLpMethods(dataDir, "lseu", 834.6823529411765);

  testResolveAfterFixing(dataDir, "p0033", 2520.5717391304347);
  testResolveAfterFixing(dataDir, "lseu", 834.6823529411765);

  testResolveAfterConstraintChange(dataDir, "p0033", 2520.5717391304347,
    "R127", 2035.9932608695651);
  testResolveAfterConstraintChange(dataDir, "lseu", 834.6823529411765,
    "R102", 662.9019607843137);

  testRacingLp(dataDir, "p0033", 2520.5717391304347);
  testRacingLp(dataDir, "lseu", 834.6823529411765);

  testPerturbationPaths(dataDir, "p0033");
  testPerturbationPaths(dataDir, "lseu");

  testPositiveEdgePersistence(dataDir, "p0033", 2520.5717391304347);
  testPositiveEdgePersistence(dataDir, "lseu", 834.6823529411765);

  testLpTimeLimitClock(dataDir, "p0033");
  testLpTimeLimitClock(dataDir, "lseu");

  printf("\n=== Summary ===\n");
  if (g_failures == 0) {
    printf("RESULT: OK (all checks passed)\n");
    return 0;
  }
  printf("RESULT: FAILED (%d check(s) failed)\n", g_failures);
  return 1;
}
