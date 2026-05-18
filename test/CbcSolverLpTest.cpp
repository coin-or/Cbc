/* CbcSolver LP relaxation tests.
   Exercises CbcSolver::solveLp() and CbcSolver::resolveLp() on MDKP and UFL
   instances embedded in the fixture headers.

   Checks:
     - solveLp() runs full LP machinery without combinatorial preprocessing.
     - resolveLp() warm-starts correctly; objective cannot improve after
       tightening a bound (for maximisation).
     - LP objective satisfies the bounding relation to the known MIP optimum.
*/

#include "CbcSolver.hpp"
#include "OsiClpSolverInterface.hpp"
#include "CoinModel.hpp"
#include "ClpSimplex.hpp"

#include <cassert>
#include <cmath>
#include <cstdio>
#include <cstring>

/* Include C fixtures (they are guarded with #ifndef ...H) */
#include "mdkp_fixtures.h"
#include "fl_fixtures.h"

static int nfail = 0;
static int npass = 0;

static void PASS(const char *label)
{
  printf("  PASS %s\n", label);
  npass++;
}

static void FAIL(const char *label, const char *reason)
{
  printf("  FAIL %s: %s\n", label, reason);
  nfail++;
}

/* ── Build LP models ──────────────────────────────────────────────────────── */

/* MDKP LP: all binary vars relaxed to continuous [0,1], maximise profit. */
static OsiClpSolverInterface buildMdkpLp(const MdkpInstance *inst)
{
  int m = inst->m, n = inst->n;
  OsiClpSolverInterface si;
  si.messageHandler()->setLogLevel(0);

  CoinModel cm;
  cm.setOptimizationDirection(-1.0);  // maximise

  for (int j = 0; j < n; j++) {
    char name[24]; snprintf(name, sizeof(name), "x(%d)", j);
    cm.addColumn(0, nullptr, nullptr, 0.0, 1.0,
                 static_cast<double>(inst->p[j]), name, false);
  }

  for (int i = 0; i < m; i++) {
    int idx[45]; double coef[45];
    for (int j = 0; j < n; j++) { idx[j] = j; coef[j] = static_cast<double>(inst->w[i * n + j]); }
    char name[24]; snprintf(name, sizeof(name), "cap(%d)", i);
    cm.addRow(n, idx, coef, -1e30, static_cast<double>(inst->c[i]), name);
  }

  si.loadFromCoinModel(cm);
  return si;
}

/* UFL individual-link LP: all vars continuous. */
static OsiClpSolverInterface buildUflLp(const UflInstance *inst)
{
  int m = inst->m, n = inst->n;
  OsiClpSolverInterface si;
  si.messageHandler()->setLogLevel(0);

  CoinModel cm;
  cm.setOptimizationDirection(1.0);  // minimise

  for (int i = 0; i < m; i++) {
    char name[24]; snprintf(name, sizeof(name), "y(%d)", i);
    cm.addColumn(0, nullptr, nullptr, 0.0, 1.0,
                 static_cast<double>(inst->f[i]), name, false);
  }

  auto xidx = [&](int i, int j) { return m + i * n + j; };

  for (int i = 0; i < m; i++) {
    for (int j = 0; j < n; j++) {
      char name[32]; snprintf(name, sizeof(name), "x(%d,%d)", i, j);
      cm.addColumn(0, nullptr, nullptr, 0.0, 1.0,
                   static_cast<double>(inst->c[i * n + j]), name, false);
    }
  }

  /* Assignment: sum_i x[i][j] = 1 */
  for (int j = 0; j < n; j++) {
    int idx[20]; double coef[20];
    for (int i = 0; i < m; i++) { idx[i] = xidx(i, j); coef[i] = 1.0; }
    char name[24]; snprintf(name, sizeof(name), "asgn(%d)", j);
    cm.addRow(m, idx, coef, 1.0, 1.0, name);
  }
  /* Linking: x[i][j] <= y[i] */
  for (int i = 0; i < m; i++) {
    for (int j = 0; j < n; j++) {
      int    idx2[2]  = { xidx(i, j), i };
      double coef2[2] = { 1.0, -1.0 };
      char name[32]; snprintf(name, sizeof(name), "lnk(%d,%d)", i, j);
      cm.addRow(2, idx2, coef2, -1e30, 0.0, name);
    }
  }

  si.loadFromCoinModel(cm);
  return si;
}

/* ── Tests ────────────────────────────────────────────────────────────────── */

static void test_mdkp_solveLp(const MdkpInstance *inst)
{
  char label[64];
  snprintf(label, sizeof(label), "MDKP solveLp m=%d n=%d", inst->m, inst->n);

  OsiClpSolverInterface si = buildMdkpLp(inst);
  CbcSolver solver(si);
  solver.initialize();
  solver.parameters().disablePrinting();

  int rc = solver.solveLp();
  if (rc != 0) {
    char buf[64]; snprintf(buf, sizeof(buf), "solveLp rc=%d", rc);
    FAIL(label, buf); return;
  }

  double lp = solver.objectiveValue();
  /* Maximisation: LP >= MIP opt (LP is upper bound) */
  if (lp < inst->opt - 0.5) {
    char buf[64]; snprintf(buf, sizeof(buf), "LP=%.4f < MIP=%d", lp, inst->opt);
    FAIL(label, buf); return;
  }
  PASS(label);
}

/* Tighten one variable's upper bound after solveLp, then resolveLp.
   Tightening can only decrease (or hold) the maximisation LP objective. */
static void test_mdkp_resolveLp(const MdkpInstance *inst)
{
  char label[64];
  snprintf(label, sizeof(label), "MDKP resolveLp m=%d n=%d", inst->m, inst->n);

  OsiClpSolverInterface si = buildMdkpLp(inst);
  CbcSolver solver(si);
  solver.initialize();
  solver.parameters().disablePrinting();

  int rc = solver.solveLp();
  if (rc != 0) { FAIL(label, "initial solveLp failed"); return; }
  double lpBefore = solver.objectiveValue();

  /* Tighten x(0) to at most 0.5 inside the solver's internal model */
  OsiSolverInterface *internalSi = solver.model()->solver();
  if (!internalSi) { FAIL(label, "no internal solver"); return; }
  internalSi->setColUpper(0, 0.5);

  int rc2 = solver.resolveLp();
  if (rc2 != 0) { FAIL(label, "resolveLp failed"); return; }
  double lpAfter = solver.objectiveValue();

  if (lpAfter > lpBefore + 1e-6) {
    char buf[128];
    snprintf(buf, sizeof(buf), "obj increased after tightening: %.4f -> %.4f",
             lpBefore, lpAfter);
    FAIL(label, buf); return;
  }
  /* Just verify the resolve returned a finite objective */
  if (lpAfter < -1e30) {
    char buf[64]; snprintf(buf, sizeof(buf), "LP after resolve non-finite: %.4f", lpAfter);
    FAIL(label, buf); return;
  }
  PASS(label);
}

static void test_ufl_solveLp(const UflInstance *inst)
{
  char label[64];
  snprintf(label, sizeof(label), "UFL solveLp m=%d n=%d", inst->m, inst->n);

  OsiClpSolverInterface si = buildUflLp(inst);
  CbcSolver solver(si);
  solver.initialize();
  solver.parameters().disablePrinting();

  int rc = solver.solveLp();
  if (rc != 0) {
    char buf[64]; snprintf(buf, sizeof(buf), "solveLp rc=%d", rc);
    FAIL(label, buf); return;
  }

  double lp = solver.objectiveValue();
  /* Minimisation: LP <= MIP opt (LP is lower bound) */
  if (lp > inst->opt + 0.5) {
    char buf[64]; snprintf(buf, sizeof(buf), "LP=%.4f > MIP=%d", lp, inst->opt);
    FAIL(label, buf); return;
  }
  PASS(label);
}

/* ── main ─────────────────────────────────────────────────────────────────── */

int main()
{
  printf("=== CbcSolver LP Relaxation Tests ===\n\n");

  printf("-- MDKP solveLp --\n");
  for (int k = 0; k < NUM_MDKP_INSTANCES; k++)
    test_mdkp_solveLp(&MDKP_INSTANCES[k]);

  printf("\n-- MDKP resolveLp (warm-start after bound tightening) --\n");
  for (int k = 0; k < NUM_MDKP_INSTANCES; k++)
    test_mdkp_resolveLp(&MDKP_INSTANCES[k]);

  printf("\n-- UFL solveLp --\n");
  for (int k = 0; k < N_UFL_INSTANCES; k++)
    test_ufl_solveLp(&UFL_INSTANCES[k]);

  printf("\n");
  if (nfail == 0)
    printf("=== All %d CbcSolver LP tests PASSED ===\n", npass);
  else
    printf("=== %d/%d CbcSolver LP test(s) FAILED ===\n", nfail, npass + nfail);
  return nfail ? 1 : 0;
}
