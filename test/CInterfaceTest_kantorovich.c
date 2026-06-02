/* Cutting Stock - Kantorovich formulation (compact MIP).
   Ported from python-mip/examples/cuttingstock_kantorovich.py.
   
   Tests: model building, MIP solve, parameter setting.
   All data embedded. */

#include "Cbc_C_Interface.h"
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

#define N_BARS 10   /* maximum number of bars */
#define L 250       /* bar length */
#define M 4         /* number of item types */

static int w[M] = {187, 119, 74, 90};
static int b[M] = {1, 2, 2, 1};

int main() {
  printf("=== Cutting Stock - Kantorovich (C interface) ===\n");

  Cbc_Model *model = Cbc_newModel();
  Cbc_setObjSense(model, 1.0); /* minimize */
  Cbc_setParam(model, "log", "0");
  Cbc_setParam(model, "sec", "60");

  /* Variables: x[i][j] = how many of item i cut from bar j (integer)
                y[j] = 1 if bar j is used (binary) */
  /* x[i][j]: index = i*N_BARS + j */
  int xStart = 0;
  for (int i = 0; i < M; i++) {
    for (int j = 0; j < N_BARS; j++) {
      char name[32];
      snprintf(name, sizeof(name), "x_%d_%d", i, j);
      Cbc_addCol(model, name, 0.0, (double)(L / w[i]), 0.0, 1, 0, NULL, NULL);
    }
  }
  /* y[j]: index = M*N_BARS + j */
  int yStart = M * N_BARS;
  for (int j = 0; j < N_BARS; j++) {
    char name[16];
    snprintf(name, sizeof(name), "y_%d", j);
    Cbc_addCol(model, name, 0.0, 1.0, 1.0, 1, 0, NULL, NULL);
  }

  /* Demand constraints: sum_j x[i][j] >= b[i] */
  for (int i = 0; i < M; i++) {
    int idx[N_BARS];
    double coef[N_BARS];
    for (int j = 0; j < N_BARS; j++) {
      idx[j] = xStart + i * N_BARS + j;
      coef[j] = 1.0;
    }
    char name[32];
    snprintf(name, sizeof(name), "demand_%d", i);
    Cbc_addRow(model, name, N_BARS, idx, coef, 'G', (double)b[i]);
  }

  /* Capacity constraints: sum_i w[i]*x[i][j] <= L*y[j] */
  for (int j = 0; j < N_BARS; j++) {
    int idx[M + 1];
    double coef[M + 1];
    for (int i = 0; i < M; i++) {
      idx[i] = xStart + i * N_BARS + j;
      coef[i] = (double)w[i];
    }
    idx[M] = yStart + j;
    coef[M] = -(double)L;
    char name[32];
    snprintf(name, sizeof(name), "cap_%d", j);
    Cbc_addRow(model, name, M + 1, idx, coef, 'L', 0.0);
  }

  /* Symmetry breaking: y[j-1] >= y[j] */
  for (int j = 1; j < N_BARS; j++) {
    int idx[2] = {yStart + j - 1, yStart + j};
    double coef[2] = {1.0, -1.0};
    char name[32];
    snprintf(name, sizeof(name), "sym_%d", j);
    Cbc_addRow(model, name, 2, idx, coef, 'G', 0.0);
  }

  /* Solve */
  Cbc_solve(model);
  assert(Cbc_isProvenOptimal(model));

  double obj = Cbc_getObjValue(model);
  printf("Optimal: %.0f bars\n", obj);
  assert(fabs(obj - 3.0) < 1e-4);

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

  /* Print solution */
  const double *sol = Cbc_getColSolution(model);
  for (int j = 0; j < N_BARS; j++) {
    if (sol[yStart + j] > 0.5) {
      printf("  Bar %d:", j);
      for (int i = 0; i < M; i++) {
        int cnt = (int)(sol[xStart + i * N_BARS + j] + 0.5);
        if (cnt > 0)
          printf(" %dx item%d(%d)", cnt, i, w[i]);
      }
      printf("\n");
    }
  }

  Cbc_deleteModel(model);
  printf("Kantorovich test PASSED!\n");
  return 0;
}
