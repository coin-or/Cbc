/* Tests for SOS constraints, lazy constraints, and row/col deletion.
   All data embedded. */

#include "Cbc_C_Interface.h"
#include <assert.h>
#include <math.h>
#include <stdio.h>

/* ── Test 1: SOS Type 1 ──────────────────────────────────────────── */
/* At most one of x0, x1, x2 can be nonzero (SOS1).
   max 3*x0 + 2*x1 + 5*x2, x0+x1+x2 <= 10, xi >= 0 */
void testSOS1(void) {
  printf("Test SOS1: ");
  Cbc_Model *m = Cbc_newModel();
  Cbc_setObjSense(m, -1.0);
  Cbc_setParam(m, "log", "0");

  Cbc_addCol(m, "x0", 0.0, 10.0, 3.0, 0, 0, NULL, NULL);
  Cbc_addCol(m, "x1", 0.0, 10.0, 2.0, 0, 0, NULL, NULL);
  Cbc_addCol(m, "x2", 0.0, 10.0, 5.0, 0, 0, NULL, NULL);

  int idx[] = {0, 1, 2};
  double coef[] = {1.0, 1.0, 1.0};
  Cbc_addRow(m, "cap", 3, idx, coef, 'L', 10.0);

  /* SOS1: at most one of x0, x1, x2 nonzero */
  int sosStart[] = {0, 3};
  int sosIdx[] = {0, 1, 2};
  double sosWt[] = {1.0, 2.0, 3.0};
  Cbc_addSOS(m, 1, sosStart, sosIdx, sosWt, 1);

  Cbc_solve(m);
  assert(Cbc_isProvenOptimal(m));

  double obj = Cbc_getObjValue(m);
  const double *sol = Cbc_getColSolution(m);

  /* With SOS1, only one variable can be nonzero.
     Best: x2=10, obj=50 */
  assert(fabs(obj - 50.0) < 1e-4);

  /* Verify SOS1: at most one nonzero */
  int nNonzero = 0;
  for (int i = 0; i < 3; i++)
    if (sol[i] > 1e-6) nNonzero++;
  assert(nNonzero <= 1);

  Cbc_deleteModel(m);
  printf("obj=%.0f PASSED\n", obj);
}

/* ── Test 2: Lazy constraints (TSP subtour as lazy) ──────────────── */
/* Solve small TSP with subtour elimination as lazy constraint */
void testLazy(void) {
  printf("Test Lazy: ");
  Cbc_Model *m = Cbc_newModel();
  Cbc_setObjSense(m, 1.0);
  Cbc_setParam(m, "log", "0");

  /* 4-city TSP: 0-1-2-3 */
  int dist[4][4] = {{0,10,15,20},{12,0,35,25},{14,32,0,30},{19,24,28,0}};
  int nArcs = 0;
  int from[12], to[12];
  for (int i = 0; i < 4; i++)
    for (int j = 0; j < 4; j++) {
      if (i == j) continue;
      char name[8];
      snprintf(name, sizeof(name), "x%d%d", i, j);
      Cbc_addCol(m, name, 0.0, 1.0, (double)dist[i][j], 1, 0, NULL, NULL);
      from[nArcs] = i; to[nArcs] = j; nArcs++;
    }

  /* Degree constraints */
  for (int node = 0; node < 4; node++) {
    int outI[3], inI[3]; double c[3] = {1,1,1}; int cnt = 0;
    for (int a = 0; a < nArcs; a++) {
      if (from[a] == node) outI[cnt] = a;
      if (to[a] == node) inI[cnt] = a;
      if (from[a] == node) cnt++;
    }
    cnt = 0;
    for (int a = 0; a < nArcs; a++) if (from[a] == node) { outI[cnt] = a; cnt++; }
    Cbc_addRow(m, "", 3, outI, c, 'E', 1.0);
    cnt = 0;
    for (int a = 0; a < nArcs; a++) if (to[a] == node) { inI[cnt] = a; cnt++; }
    Cbc_addRow(m, "", 3, inI, c, 'E', 1.0);
  }

  /* Add subtour elimination for {0,1} as LAZY constraint:
     x[0][1] + x[1][0] <= 1 */
  {
    int idx[2]; double coef[2] = {1.0, 1.0};
    for (int a = 0; a < nArcs; a++) {
      if (from[a]==0 && to[a]==1) idx[0] = a;
      if (from[a]==1 && to[a]==0) idx[1] = a;
    }
    Cbc_addLazyConstraint(m, 2, idx, coef, 'L', 1.0);
  }
  /* SEC for {2,3} */
  {
    int idx[2]; double coef[2] = {1.0, 1.0};
    for (int a = 0; a < nArcs; a++) {
      if (from[a]==2 && to[a]==3) idx[0] = a;
      if (from[a]==3 && to[a]==2) idx[1] = a;
    }
    Cbc_addLazyConstraint(m, 2, idx, coef, 'L', 1.0);
  }

  Cbc_solve(m);
  assert(Cbc_isProvenOptimal(m));
  double obj = Cbc_getObjValue(m);
  /* Optimal 4-city tour with these distances: 0->1->3->2->0 = 10+25+30+14 = 79 */
  printf("obj=%.0f ", obj);
  assert(obj >= 70 && obj <= 100);

  Cbc_deleteModel(m);
  printf("PASSED\n");
}

/* ── Test 3: Delete rows and columns ─────────────────────────────── */
void testDelete(void) {
  printf("Test Delete: ");
  Cbc_Model *m = Cbc_newModel();
  Cbc_setObjSense(m, -1.0);
  Cbc_setParam(m, "log", "0");

  /* max 5x + 4y + 3z, x+y+z <= 10 */
  Cbc_addCol(m, "x", 0.0, 10.0, 5.0, 0, 0, NULL, NULL);
  Cbc_addCol(m, "y", 0.0, 10.0, 4.0, 0, 0, NULL, NULL);
  Cbc_addCol(m, "z", 0.0, 10.0, 3.0, 0, 0, NULL, NULL);
  int idx[] = {0, 1, 2};
  double coef[] = {1.0, 1.0, 1.0};
  Cbc_addRow(m, "cap", 3, idx, coef, 'L', 10.0);

  /* Add a tight constraint: x <= 3 */
  int idx2[] = {0};
  double coef2[] = {1.0};
  Cbc_addRow(m, "tight", 1, idx2, coef2, 'L', 3.0);

  Cbc_solveLinearProgram(m);
  double obj1 = Cbc_getObjValue(m);
  /* With x<=3: x=3, y=7, z=0 → obj = 15+28 = 43 */
  assert(fabs(obj1 - 43.0) < 1e-4);

  /* Delete the tight constraint (row 1) */
  int delRow[] = {1};
  Cbc_deleteRows(m, 1, delRow);

  Cbc_resolve(m);
  double obj2 = Cbc_getObjValue(m);
  /* Without x<=3: x=10, y=0, z=0 → obj = 50 */
  assert(fabs(obj2 - 50.0) < 1e-4);

  /* Delete column z (col 2) */
  int delCol[] = {2};
  Cbc_deleteCols(m, 1, delCol);
  assert(Cbc_getNumCols(m) == 2);

  Cbc_resolve(m);
  double obj3 = Cbc_getObjValue(m);
  /* Only x,y with x+y<=10: x=10, y=0 → obj=50 (same) */
  assert(fabs(obj3 - 50.0) < 1e-4);

  Cbc_deleteModel(m);
  printf("obj1=%.0f obj2=%.0f obj3=%.0f PASSED\n", obj1, obj2, obj3);
}

int main() {
  printf("=== SOS, Lazy, Delete Tests ===\n");
  testSOS1();
  testLazy();
  testDelete();
  printf("All tests PASSED!\n");
  return 0;
}
