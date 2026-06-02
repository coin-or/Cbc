/* TSP solved by cutting planes (subtour elimination).
   Simplified port of python-mip/examples/cutting_planes.py.
   Uses manual BFS for subtour detection (no external dependencies).
   
   Tests: LP relaxation, Cbc_resolve, Cbc_addRow (cuts), iterative solving.
   All data embedded — 5-city asymmetric TSP. */

#include "Cbc_C_Interface.h"
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

#define NCITIES 5
#define MAXARCS (NCITIES * (NCITIES - 1))

/* Distance matrix (asymmetric) */
static int dist[NCITIES][NCITIES] = {
  { 0, 10, 15, 20, 25},
  {12,  0, 35, 25, 18},
  {14, 32,  0, 30, 20},
  {19, 24, 28,  0, 15},
  {22, 17, 21, 14,  0}
};

/* Find subtour containing node 0 using BFS on fractional solution */
static int findSubtour(int n, const double *xval, int *arcFrom, int *arcTo,
  int nArcs, int *subtour)
{
  char visited[NCITIES] = {0};
  int queue[NCITIES], qHead = 0, qTail = 0;
  queue[qTail++] = 0;
  visited[0] = 1;
  int count = 0;

  while (qHead < qTail) {
    int node = queue[qHead++];
    subtour[count++] = node;
    for (int a = 0; a < nArcs; a++) {
      if (arcFrom[a] == node && xval[a] > 0.5 && !visited[arcTo[a]]) {
        visited[arcTo[a]] = 1;
        queue[qTail++] = arcTo[a];
      }
    }
  }
  return count; /* size of component containing node 0 */
}

int main() {
  printf("=== TSP Cutting Planes (C interface) ===\n");

  Cbc_Model *model = Cbc_newModel();
  Cbc_setObjSense(model, 1.0); /* minimize */
  Cbc_setParam(model, "log", "0");

  /* Arc variables: x[i][j] binary, index = i*(n-1) + (j > i ? j-1 : j) */
  int nArcs = 0;
  int arcFrom[MAXARCS], arcTo[MAXARCS];
  for (int i = 0; i < NCITIES; i++) {
    for (int j = 0; j < NCITIES; j++) {
      if (i == j) continue;
      char name[16];
      snprintf(name, sizeof(name), "x_%d_%d", i, j);
      Cbc_addCol(model, name, 0.0, 1.0, (double)dist[i][j], 1, 0, NULL, NULL);
      arcFrom[nArcs] = i;
      arcTo[nArcs] = j;
      nArcs++;
    }
  }

  /* Degree constraints: out-degree = 1, in-degree = 1 */
  for (int node = 0; node < NCITIES; node++) {
    /* out-degree */
    int idx[NCITIES - 1];
    double coef[NCITIES - 1];
    int cnt = 0;
    for (int a = 0; a < nArcs; a++) {
      if (arcFrom[a] == node) {
        idx[cnt] = a;
        coef[cnt] = 1.0;
        cnt++;
      }
    }
    char name[16];
    snprintf(name, sizeof(name), "out_%d", node);
    Cbc_addRow(model, name, cnt, idx, coef, 'E', 1.0);

    /* in-degree */
    cnt = 0;
    for (int a = 0; a < nArcs; a++) {
      if (arcTo[a] == node) {
        idx[cnt] = a;
        coef[cnt] = 1.0;
        cnt++;
      }
    }
    snprintf(name, sizeof(name), "in_%d", node);
    Cbc_addRow(model, name, cnt, idx, coef, 'E', 1.0);
  }

  /* Cutting plane loop */
  int iter = 0;
  while (1) {
    iter++;

    /* Solve LP relaxation */
    int rc;
    if (iter == 1)
      rc = Cbc_solveLinearProgram(model);
    else
      rc = Cbc_resolve(model);
    if (rc != 0) break;

    double obj = Cbc_getObjValue(model);
    printf("Iter %d: LP obj = %.2f\n", iter, obj);

    /* Check for subtours */
    const double *xval = Cbc_getColSolution(model);
    int subtour[NCITIES];
    int subSize = findSubtour(NCITIES, xval, arcFrom, arcTo, nArcs, subtour);

    if (subSize == NCITIES) {
      /* No subtour — we have a valid tour */
      printf("  No subtour found — optimal tour!\n");
      break;
    }

    /* Add subtour elimination constraint: sum of arcs within S <= |S|-1 */
    printf("  Subtour of size %d, adding SEC\n", subSize);
    char inS[NCITIES] = {0};
    for (int k = 0; k < subSize; k++)
      inS[subtour[k]] = 1;

    int cutIdx[MAXARCS];
    double cutCoef[MAXARCS];
    int cutNz = 0;
    for (int a = 0; a < nArcs; a++) {
      if (inS[arcFrom[a]] && inS[arcTo[a]]) {
        cutIdx[cutNz] = a;
        cutCoef[cutNz] = 1.0;
        cutNz++;
      }
    }
    char cname[32];
    snprintf(cname, sizeof(cname), "sec_%d", iter);
    Cbc_addRow(model, cname, cutNz, cutIdx, cutCoef, 'L', (double)(subSize - 1));

    if (iter > 20) {
      printf("  Too many iterations, stopping.\n");
      break;
    }
  }

  /* Now solve as MIP */
  printf("\nSolving MIP...\n");
  Cbc_solve(model);
  assert(Cbc_isProvenOptimal(model));

  double obj = Cbc_getObjValue(model);
  printf("Optimal tour cost: %.0f\n", obj);

  /* Validate best solution feasibility */
  const double *bestSol = Cbc_getColSolution(model);
  double maxViolRow = 0.0, maxViolCol = 0.0;
  int rowIdx = -1, colIdx = -1;
  if (!Cbc_checkFeasibility(model, bestSol, &maxViolRow, &rowIdx, &maxViolCol, &colIdx)) {
    printf("FAIL: best solution infeasible — rowViol=%.2e (row %d) colViol=%.2e (col %d)\n",
           maxViolRow, rowIdx, maxViolCol, colIdx);
    assert(0);
  }

  /* Validate ALL solutions in the pool */
  int nSol = Cbc_numberSavedSolutions(model);
  for (int s = 0; s < nSol; s++) {
    const double *sol = Cbc_savedSolution(model, s);
    double solObj = Cbc_savedSolutionObj(model, s);
    maxViolRow = 0.0; maxViolCol = 0.0;
    rowIdx = -1; colIdx = -1;

    if (!Cbc_checkFeasibility(model, sol, &maxViolRow, &rowIdx, &maxViolCol, &colIdx)) {
      printf("FAIL: pool solution %d/%d infeasible (obj=%.0f) — "
             "rowViol=%.2e (row %d) colViol=%.2e (col %d)\n",
             s+1, nSol, solObj, maxViolRow, rowIdx, maxViolCol, colIdx);
      assert(0);
    }
  }

  /* Print tour */
  const double *sol = Cbc_getColSolution(model);
  printf("Tour: 0");
  int current = 0;
  for (int step = 0; step < NCITIES - 1; step++) {
    for (int a = 0; a < nArcs; a++) {
      if (arcFrom[a] == current && sol[a] > 0.5) {
        current = arcTo[a];
        printf(" -> %d", current);
        break;
      }
    }
  }
  printf(" -> 0\n");

  /* Verify: optimal tour cost is 84 for this 5-city instance
     (tour: 0->1->4->3->2->0 = 10+18+14+28+14 = 84) */
  assert(fabs(obj - 84.0) < 1e-4);

  Cbc_deleteModel(model);
  printf("TSP cutting planes test PASSED!\n");
  return 0;
}
