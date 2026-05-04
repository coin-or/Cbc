/* Cutting Stock Problem solved by Column Generation using the CBC C interface.
   
   Ported from python-mip/examples/cuttingstock_cg.py.
   All data is embedded — no file I/O needed.
   
   Tests: Cbc_addCol (with column coefficients), Cbc_solveLinearProgram,
   Cbc_resolve, Cbc_getRowPrice (duals), Cbc_solve (pricing MIP).
*/

#include "Cbc_C_Interface.h"
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

#define L 250       /* bar length */
#define M 4         /* number of item types */

static int w[M] = {187, 119, 74, 90};  /* item widths */
static int b[M] = {1, 2, 2, 1};        /* demand for each item */

int main() {
  printf("=== Cutting Stock Column Generation (C interface) ===\n\n");

  /* ── Build master problem ── */
  Cbc_Model *master = Cbc_newModel();
  Cbc_setObjSense(master, 1.0); /* minimize */
  Cbc_setParam(master, "log", "0");

  /* Initial patterns: one item per bar (identity columns) */
  /* Add columns first (no row references yet) */
  for (int j = 0; j < M; j++) {
    char name[32];
    snprintf(name, sizeof(name), "lambda_%d", j + 1);
    Cbc_addCol(master, name, 0.0, 1e30, 1.0, 0, 0, NULL, NULL);
  }

  /* Demand constraints: lambda_i >= b[i] */
  for (int i = 0; i < M; i++) {
    char name[32];
    snprintf(name, sizeof(name), "demand_%d", i + 1);
    int colIdx = i;
    double coef = 1.0;
    Cbc_addRow(master, name, 1, &colIdx, &coef, 'G', (double)b[i]);
  }

  int nPatterns = M;
  int iter = 0;

  while (1) {
    iter++;

    /* ── Step 1: Solve master LP relaxation ── */
    int rc;
    if (iter == 1)
      rc = Cbc_solveLinearProgram(master);
    else
      rc = Cbc_resolve(master);
    if (rc != 0) {
      printf("  LP solve failed (rc=%d)\n", rc);
      break;
    }

    double masterObj = Cbc_getObjValue(master);
    printf("Iter %d: master obj = %.4f (%d patterns)\n", iter, masterObj, nPatterns);

    /* ── Step 2: Get duals and solve pricing problem ── */
    const double *duals = Cbc_getRowPrice(master);
    assert(duals != NULL);

    /* Build pricing problem: max sum(dual[i] * a[i]) s.t. sum(w[i]*a[i]) <= L */
    Cbc_Model *pricing = Cbc_newModel();
    Cbc_setObjSense(pricing, -1.0); /* maximize */
    Cbc_setParam(pricing, "log", "0");

    for (int i = 0; i < M; i++) {
      char name[32];
      snprintf(name, sizeof(name), "a_%d", i + 1);
      int maxCopies = L / w[i];
      /* objective = dual value for this row */
      Cbc_addCol(pricing, name, 0.0, (double)maxCopies, duals[i], 1, 0, NULL, NULL);
    }

    /* Knapsack constraint: sum(w[i] * a[i]) <= L */
    int idx[M];
    double wd[M];
    for (int i = 0; i < M; i++) {
      idx[i] = i;
      wd[i] = (double)w[i];
    }
    Cbc_addRow(pricing, "bar_length", M, idx, wd, 'L', (double)L);

    /* Solve pricing as MIP */
    Cbc_solve(pricing);
    assert(Cbc_isProvenOptimal(pricing));

    double pricingObj = Cbc_getObjValue(pricing);
    printf("  Pricing obj = %.4f (reduced cost = %.4f)\n",
      pricingObj, 1.0 - pricingObj);

    /* ── Step 3: Check reduced cost and add column ── */
    if (1.0 - pricingObj >= -1e-5) {
      /* No improving column found — done */
      Cbc_deleteModel(pricing);
      break;
    }

    /* Get the new pattern */
    const double *pattern = Cbc_getColSolution(pricing);
    printf("  New pattern: [");
    for (int i = 0; i < M; i++)
      printf("%.0f%s", pattern[i], i < M - 1 ? ", " : "");
    printf("]\n");

    /* Add new column to master */
    char name[32];
    snprintf(name, sizeof(name), "lambda_%d", nPatterns + 1);
    double colCoefs[M];
    int rowIdx[M];
    for (int i = 0; i < M; i++) {
      rowIdx[i] = i;
      colCoefs[i] = pattern[i];
    }
    Cbc_addCol(master, name, 0.0, 1e30, 1.0, 0, M, rowIdx, colCoefs);
    nPatterns++;

    Cbc_deleteModel(pricing);
  }

  /* Final result */
  double finalObj = Cbc_getObjValue(master);
  printf("\nOptimal LP value: %.4f (need at least %.0f bars)\n",
    finalObj, ceil(finalObj));

  /* Verify: LP relaxation optimal is exactly 3 bars for this instance
     (same data as python-mip cuttingstock_cg.py) */
  assert(fabs(finalObj - 3.0) < 1e-4);

  const double *sol = Cbc_getColSolution(master);
  printf("Pattern usage:\n");
  for (int j = 0; j < nPatterns; j++) {
    if (sol[j] > 1e-6)
      printf("  lambda_%d = %.4f\n", j + 1, sol[j]);
  }

  Cbc_deleteModel(master);
  printf("\nColumn generation test PASSED!\n");
  return 0;
}
