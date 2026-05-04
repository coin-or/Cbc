/* MIP reoptimization: solve TSP, forbid a 3-node subsequence, re-solve.
   
   Tests: Cbc_solve, add constraint, Cbc_solve again (MIP reoptimization).
   Verifies second solution is worse and doesn't contain the forbidden path. */

#include "Cbc_C_Interface.h"
#include <assert.h>
#include <math.h>
#include <stdio.h>

#define N 5

static int dist[N][N] = {
  { 0, 10, 15, 20, 25},
  {12,  0, 35, 25, 18},
  {14, 32,  0, 30, 20},
  {19, 24, 28,  0, 15},
  {22, 17, 21, 14,  0}
};

/* arc variable index for (i,j) */
static int arcIdx(int i, int j) {
  /* skip diagonal: index = i*(N-1) + (j > i ? j-1 : j) */
  return i * (N - 1) + (j > i ? j - 1 : j);
}

static Cbc_Model *buildTSP(void) {
  Cbc_Model *m = Cbc_newModel();
  Cbc_setObjSense(m, 1.0);
  Cbc_setParam(m, "log", "0");

  /* Binary arc variables */
  for (int i = 0; i < N; i++)
    for (int j = 0; j < N; j++) {
      if (i == j) continue;
      char name[12];
      snprintf(name, sizeof(name), "x_%d_%d", i, j);
      Cbc_addCol(m, name, 0.0, 1.0, (double)dist[i][j], 1, 0, NULL, NULL);
    }

  /* Degree constraints */
  for (int node = 0; node < N; node++) {
    int outIdx[N - 1], inIdx[N - 1];
    double coef[N - 1];
    int cnt = 0;
    for (int j = 0; j < N; j++) {
      if (j == node) continue;
      outIdx[cnt] = arcIdx(node, j);
      inIdx[cnt] = arcIdx(j, node);
      coef[cnt] = 1.0;
      cnt++;
    }
    char name[12];
    snprintf(name, sizeof(name), "out_%d", node);
    Cbc_addRow(m, name, cnt, outIdx, coef, 'E', 1.0);
    snprintf(name, sizeof(name), "in_%d", node);
    Cbc_addRow(m, name, cnt, inIdx, coef, 'E', 1.0);
  }

  /* Subtour elimination: for all subsets S with 2 <= |S| <= N-2,
     sum of arcs within S <= |S| - 1.
     For this small instance, add all 2-node SECs (equivalent to
     no 2-cycles: x[i][j] + x[j][i] <= 1) */
  for (int i = 0; i < N; i++)
    for (int j = i + 1; j < N; j++) {
      int idx[2] = {arcIdx(i, j), arcIdx(j, i)};
      double coef[2] = {1.0, 1.0};
      char name[16];
      snprintf(name, sizeof(name), "sec2_%d_%d", i, j);
      Cbc_addRow(m, name, 2, idx, coef, 'L', 1.0);
    }

  /* 3-node SECs */
  for (int i = 0; i < N; i++)
    for (int j = 0; j < N; j++) {
      if (j == i) continue;
      for (int k = 0; k < N; k++) {
        if (k == i || k == j) continue;
        int idx[6] = {arcIdx(i,j), arcIdx(j,i), arcIdx(i,k),
                      arcIdx(k,i), arcIdx(j,k), arcIdx(k,j)};
        double coef[6] = {1,1,1,1,1,1};
        char name[20];
        snprintf(name, sizeof(name), "sec3_%d_%d_%d", i, j, k);
        Cbc_addRow(m, name, 6, idx, coef, 'L', 2.0);
      }
    }

  return m;
}

int main() {
  printf("=== MIP Reoptimization Test ===\n");

  Cbc_Model *m = buildTSP();

  /* First solve */
  Cbc_solve(m);
  assert(Cbc_isProvenOptimal(m));
  double obj1 = Cbc_getObjValue(m);
  printf("First solve: obj = %.0f\n", obj1);
  assert(fabs(obj1 - 84.0) < 1e-4);

  /* Extract the tour */
  const double *sol = Cbc_getColSolution(m);
  printf("Tour 1: ");
  int tour1[N + 1];
  int cur = 0;
  tour1[0] = 0;
  for (int step = 0; step < N - 1; step++) {
    for (int j = 0; j < N; j++) {
      if (j == cur) continue;
      if (sol[arcIdx(cur, j)] > 0.5) {
        cur = j;
        tour1[step + 1] = cur;
        break;
      }
    }
  }
  tour1[N] = 0;
  for (int i = 0; i <= N; i++)
    printf("%d%s", tour1[i], i < N ? " -> " : "\n");

  /* Forbid the 3-node path that appears in the optimal tour.
     If tour is 0->1->4->3->2->0, forbid path 1->4->3:
     x[1][4] + x[4][3] <= 1 (can't use both arcs) */
  int a = tour1[1], b = tour1[2], c = tour1[3];
  printf("Forbidding path %d -> %d -> %d\n", a, b, c);
  {
    int idx[2] = {arcIdx(a, b), arcIdx(b, c)};
    double coef[2] = {1.0, 1.0};
    Cbc_addRow(m, "forbid", 2, idx, coef, 'L', 1.0);
  }

  /* Re-solve MIP */
  Cbc_solve(m);
  assert(Cbc_isProvenOptimal(m));
  double obj2 = Cbc_getObjValue(m);
  printf("Second solve: obj = %.0f\n", obj2);

  /* Second solution must be worse (or equal, but not better) */
  assert(obj2 >= obj1 - 1e-6);
  /* And strictly worse for this instance */
  assert(obj2 > obj1 + 0.5);

  /* Verify the forbidden path is not in the new solution */
  const double *sol2 = Cbc_getColSolution(m);
  double pathUsage = sol2[arcIdx(a, b)] + sol2[arcIdx(b, c)];
  printf("Forbidden path usage: %.1f (must be <= 1)\n", pathUsage);
  assert(pathUsage <= 1.0 + 1e-6);

  /* Print new tour */
  printf("Tour 2: ");
  cur = 0;
  printf("0");
  for (int step = 0; step < N - 1; step++) {
    for (int j = 0; j < N; j++) {
      if (j == cur) continue;
      if (sol2[arcIdx(cur, j)] > 0.5) {
        cur = j;
        printf(" -> %d", cur);
        break;
      }
    }
  }
  printf(" -> 0\n");

  Cbc_deleteModel(m);
  printf("MIP reoptimization test PASSED!\n");
  return 0;
}
