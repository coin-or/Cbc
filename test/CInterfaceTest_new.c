/* Test for new C interface functions: Cbc_setParam, Cbc_getParam, Cbc_resolve,
   LP reoptimization, and the knapsack problem with parameter control.

   All data is embedded — no file I/O needed. */

#include "Cbc_C_Interface.h"
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

/* ── Test 1: Knapsack with Cbc_setParam ────────────────────────────── */
void testKnapsackWithParams()
{
  printf("Test 1: Knapsack with Cbc_setParam\n");

  /* max 10x0 + 13x1 + 18x2 + 31x3 + 7x4 + 15x5
     s.t. 11x0 + 15x1 + 20x2 + 35x3 + 10x4 + 33x5 <= 47
          xi in {0,1} */
  double p[] = { 10, 13, 18, 31, 7, 15 };
  double w[] = { 11, 15, 20, 35, 10, 33 };
  int n = 6;
  double cap = 47;

  Cbc_Model *m = Cbc_newModel();
  Cbc_setObjSense(m, -1.0); /* maximize */

  /* add binary variables */
  for (int i = 0; i < n; i++) {
    char name[8];
    snprintf(name, sizeof(name), "x%d", i);
    Cbc_addCol(m, name, 0.0, 1.0, p[i], 1, 0, NULL, NULL);
  }

  /* capacity constraint */
  int idx[] = { 0, 1, 2, 3, 4, 5 };
  Cbc_addRow(m, "cap", n, idx, w, 'L', cap);

  /* Set parameters using the new Cbc_setParam */
  int rc;
  rc = Cbc_setParam(m, "log", "0");
  assert(rc == 0);
  rc = Cbc_setParam(m, "sec", "60");
  assert(rc == 0);

  /* Verify Cbc_getParam */
  char buf[64];
  rc = Cbc_getParam(m, "log", buf, sizeof(buf));
  assert(rc == 0);
  assert(strcmp(buf, "0") == 0);

  /* Solve */
  Cbc_solve(m);
  assert(Cbc_isProvenOptimal(m));
  assert(fabs(Cbc_getObjValue(m) - 41.0) < 1e-4);

  printf("  Optimal obj = %g (expected 41)\n", Cbc_getObjValue(m));
  Cbc_deleteModel(m);
  printf("  PASSED\n\n");
}

/* ── Test 2: LP reoptimization with Cbc_resolve ───────────────────── */
void testLPReoptimization()
{
  printf("Test 2: LP reoptimization with Cbc_resolve\n");

  /* Simple LP: max 5x + 4y
                s.t. 6x + 4y <= 24
                     x + 2y <= 6
                     x, y >= 0 */
  Cbc_Model *m = Cbc_newModel();
  Cbc_setObjSense(m, -1.0);

  Cbc_addCol(m, "x", 0.0, 1e30, 5.0, 0, 0, NULL, NULL);
  Cbc_addCol(m, "y", 0.0, 1e30, 4.0, 0, 0, NULL, NULL);

  int idx0[] = { 0, 1 };
  double coef0[] = { 6.0, 4.0 };
  Cbc_addRow(m, "c1", 2, idx0, coef0, 'L', 24.0);

  double coef1[] = { 1.0, 2.0 };
  Cbc_addRow(m, "c2", 2, idx0, coef1, 'L', 6.0);

  Cbc_setParam(m, "log", "0");

  /* Initial solve */
  int rc = Cbc_solveLinearProgram(m);
  assert(rc == 0);
  double obj1 = Cbc_getObjValue(m);
  printf("  Initial LP obj = %g\n", obj1);
  assert(fabs(obj1 - 21.0) < 1e-4); /* optimal: x=3, y=1.5 → 5*3+4*1.5=21 */

  /* Tighten bound: x <= 2 */
  Cbc_setColUpper(m, 0, 2.0);
  rc = Cbc_resolve(m);
  assert(rc == 0);
  double obj2 = Cbc_getObjValue(m);
  printf("  After x<=2: obj = %g\n", obj2);
  assert(obj2 < obj1 - 1e-6); /* should be worse */

  /* Add constraint: x + y <= 3 */
  double coef2[] = { 1.0, 1.0 };
  Cbc_addRow(m, "c3", 2, idx0, coef2, 'L', 3.0);
  rc = Cbc_resolve(m);
  assert(rc == 0);
  double obj3 = Cbc_getObjValue(m);
  printf("  After x+y<=3: obj = %g\n", obj3);
  assert(obj3 <= obj2 + 1e-6); /* should be same or worse */

  /* Relax bound: x <= 100 */
  Cbc_setColUpper(m, 0, 100.0);
  rc = Cbc_resolve(m);
  assert(rc == 0);
  double obj4 = Cbc_getObjValue(m);
  printf("  After relaxing x: obj = %g\n", obj4);

  /* Make infeasible: x >= 10 */
  Cbc_setColLower(m, 0, 10.0);
  rc = Cbc_resolve(m);
  printf("  After x>=10: rc = %d (expect 2=infeasible)\n", rc);
  assert(rc == 2);

  Cbc_deleteModel(m);
  printf("  PASSED\n\n");
}

/* ── Test 3: MIP with solution pool ──────────────────────────────── */
void testSolutionPool()
{
  printf("Test 3: MIP with solution pool\n");

  /* Small knapsack with multiple optimal solutions */
  /* max x0 + x1 + x2, s.t. x0 + x1 + x2 <= 2, xi binary */
  Cbc_Model *m = Cbc_newModel();
  Cbc_setObjSense(m, -1.0);

  for (int i = 0; i < 3; i++) {
    char name[4];
    snprintf(name, sizeof(name), "x%d", i);
    Cbc_addCol(m, name, 0.0, 1.0, 1.0, 1, 0, NULL, NULL);
  }
  int idx[] = { 0, 1, 2 };
  double coef[] = { 1.0, 1.0, 1.0 };
  Cbc_addRow(m, "cap", 3, idx, coef, 'L', 2.0);

  Cbc_setParam(m, "log", "0");
  Cbc_setParam(m, "maxSaved", "10");

  Cbc_solve(m);
  assert(Cbc_isProvenOptimal(m));
  assert(fabs(Cbc_getObjValue(m) - 2.0) < 1e-4);

  int nSols = Cbc_numberSavedSolutions(m);
  printf("  Optimal obj = %g, saved solutions = %d\n",
    Cbc_getObjValue(m), nSols);
  assert(nSols >= 1);

  Cbc_deleteModel(m);
  printf("  PASSED\n\n");
}

/* ── Test 4: MIP start ───────────────────────────────────────────── */
void testMIPStart()
{
  printf("Test 4: MIP start\n");

  /* Knapsack: provide a feasible starting solution */
  double p[] = { 10, 13, 18, 31, 7, 15 };
  double w[] = { 11, 15, 20, 35, 10, 33 };
  int n = 6;

  Cbc_Model *m = Cbc_newModel();
  Cbc_setObjSense(m, -1.0);

  for (int i = 0; i < n; i++) {
    char name[8];
    snprintf(name, sizeof(name), "x%d", i);
    Cbc_addCol(m, name, 0.0, 1.0, p[i], 1, 0, NULL, NULL);
  }
  int idx[] = { 0, 1, 2, 3, 4, 5 };
  Cbc_addRow(m, "cap", n, idx, w, 'L', 47.0);

  /* Provide MIP start: x0=1, x2=1 (weight=31, profit=28) */
  int startIdx[] = { 0, 2 };
  double startVal[] = { 1.0, 1.0 };
  Cbc_setMIPStartI(m, 2, startIdx, startVal);

  Cbc_setParam(m, "log", "0");
  Cbc_solve(m);

  assert(Cbc_isProvenOptimal(m));
  assert(fabs(Cbc_getObjValue(m) - 41.0) < 1e-4);
  printf("  Optimal obj = %g with MIP start\n", Cbc_getObjValue(m));

  Cbc_deleteModel(m);
  printf("  PASSED\n\n");
}

static int heuristic_stats_printed = 0;

static void CBC_LINKAGE_CB stats_callback(Cbc_Model *model, int msgno, int ndouble,
  const double *dvec, int nint, const int *ivec,
  int nchar, char **cvec)
{
  if (msgno == 45 && nchar > 0 && cvec && cvec[0]) {
    if (strstr(cvec[0], "Heuristic") != NULL) {
      heuristic_stats_printed = 1;
    }
  }
}

/* ── Test 5: Heuristic stats flag (-heuristicStats) ─────────────────── */
void testHeuristicStats()
{
  printf("Test 5: Heuristic stats flag (-heuristicStats)\n");

  /* Small knapsack problem */
  double p[] = { 10, 13, 18, 31, 7, 15 };
  double w[] = { 11, 15, 20, 35, 10, 33 };
  int n = 6;
  double cap = 47;

  // Case 1: heuristicStats = off (default)
  {
    Cbc_Model *m = Cbc_newModel();
    Cbc_setObjSense(m, -1.0); /* maximize */

    for (int i = 0; i < n; i++) {
      char name[8];
      snprintf(name, sizeof(name), "x%d", i);
      Cbc_addCol(m, name, 0.0, 1.0, p[i], 1, 0, NULL, NULL);
    }
    int idx[] = { 0, 1, 2, 3, 4, 5 };
    Cbc_addRow(m, "cap", n, idx, w, 'L', cap);

    heuristic_stats_printed = 0;
    Cbc_registerCallBack(m, stats_callback);

    Cbc_setParam(m, "log", "1");
    Cbc_setParam(m, "heuristicStats", "off");

    Cbc_solve(m);

    assert(Cbc_isProvenOptimal(m));
    assert(heuristic_stats_printed == 0);

    Cbc_deleteModel(m);
  }

  // Case 2: heuristicStats = on
  {
    Cbc_Model *m = Cbc_newModel();
    Cbc_setObjSense(m, -1.0); /* maximize */

    for (int i = 0; i < n; i++) {
      char name[8];
      snprintf(name, sizeof(name), "x%d", i);
      Cbc_addCol(m, name, 0.0, 1.0, p[i], 1, 0, NULL, NULL);
    }
    int idx[] = { 0, 1, 2, 3, 4, 5 };
    Cbc_addRow(m, "cap", n, idx, w, 'L', cap);

    heuristic_stats_printed = 0;
    Cbc_registerCallBack(m, stats_callback);

    Cbc_setParam(m, "log", "1");
    Cbc_setParam(m, "heuristicStats", "on");

    Cbc_solve(m);

    assert(Cbc_isProvenOptimal(m));
    assert(heuristic_stats_printed == 1);

    Cbc_deleteModel(m);
  }

  printf("  PASSED\n\n");
}

int main()
{
  printf("=== CBC C Interface Tests ===\n\n");
  testKnapsackWithParams();
  testLPReoptimization();
  testSolutionPool();
  testMIPStart();
  testHeuristicStats();
  printf("All tests PASSED!\n");
  return 0;
}
