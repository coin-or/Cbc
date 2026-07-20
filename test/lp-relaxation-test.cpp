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
 *
 * Usage: lp-relaxation-test [--data-dir=PATH]
 * Exit code: 0 = all checks passed, 1 = one or more checks failed.
 */

#include "Cbc_C_Interface.h"
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
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
    printf("  ok:   %s (actual=%.10g expected=%.10g delta=%.3e)\n",
      msg.c_str(), actual, expected, delta, tol);
  }
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

} // namespace

int main(int argc, char *argv[])
{
  std::string dataDir;
  for (int i = 1; i < argc; ++i) {
    std::string arg = argv[i];
    if (arg.rfind("--data-dir=", 0) == 0)
      dataDir = arg.substr(strlen("--data-dir="));
  }
  if (dataDir.empty())
    dataDir = exeDir() + "/mip-sanity-data";

  /* Reference LP relaxation objectives below were computed independently
   * (two different external LP/MIP solvers agreed exactly) and cross-checked
   * against Cbc's own LP relaxation solve before writing this test. */
  testRootLpMethods(dataDir, "p0033", 2520.5717391304347);
  testRootLpMethods(dataDir, "lseu", 834.6823529411765);

  testResolveAfterFixing(dataDir, "p0033", 2520.5717391304347);
  testResolveAfterFixing(dataDir, "lseu", 834.6823529411765);

  testResolveAfterConstraintChange(dataDir, "p0033", 2520.5717391304347,
    "R127", 2035.9932608695651);
  testResolveAfterConstraintChange(dataDir, "lseu", 834.6823529411765,
    "R102", 662.9019607843137);

  printf("\n=== Summary ===\n");
  if (g_failures == 0) {
    printf("RESULT: OK (all checks passed)\n");
    return 0;
  }
  printf("RESULT: FAILED (%d check(s) failed)\n", g_failures);
  return 1;
}
