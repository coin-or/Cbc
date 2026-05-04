/* Test model construction in different orders:
   1. Row-wise: add empty columns, then add rows with column references
   2. Column-wise: add empty rows, then add columns with row references
   3. Mixed: interleave column and row additions
   
   All three should produce the same optimal solution.
   Uses a small transportation problem (3 sources, 4 destinations). */

#include "Cbc_C_Interface.h"
#include <assert.h>
#include <math.h>
#include <stdio.h>

#define NS 3  /* sources */
#define ND 4  /* destinations */

/* Supply at each source */
static double supply[NS] = {20, 30, 25};
/* Demand at each destination */
static double demand[ND] = {10, 15, 20, 30};
/* Transport cost[source][dest] */
static double cost[NS][ND] = {
  { 8, 6, 10, 9},
  {9, 12,  7, 5},
  {14, 9, 16, 5}
};
/* Expected optimal (LP relaxation of transportation problem) */
#define EXPECTED_OPT 465.0

/* Build model row-wise: add all columns first, then all rows */
static double solveRowWise(void) {
  Cbc_Model *m = Cbc_newModel();
  Cbc_setObjSense(m, 1.0);
  Cbc_setParam(m, "log", "0");

  /* Add all x[i][j] columns (no row references) */
  for (int i = 0; i < NS; i++)
    for (int j = 0; j < ND; j++) {
      char name[12];
      snprintf(name, sizeof(name), "x_%d_%d", i, j);
      Cbc_addCol(m, name, 0.0, 1e30, cost[i][j], 0, 0, NULL, NULL);
    }

  /* Supply constraints: sum_j x[i][j] <= supply[i] */
  for (int i = 0; i < NS; i++) {
    int idx[ND];
    double coef[ND];
    for (int j = 0; j < ND; j++) {
      idx[j] = i * ND + j;
      coef[j] = 1.0;
    }
    char name[12];
    snprintf(name, sizeof(name), "sup_%d", i);
    Cbc_addRow(m, name, ND, idx, coef, 'L', supply[i]);
  }

  /* Demand constraints: sum_i x[i][j] >= demand[j] */
  for (int j = 0; j < ND; j++) {
    int idx[NS];
    double coef[NS];
    for (int i = 0; i < NS; i++) {
      idx[i] = i * ND + j;
      coef[i] = 1.0;
    }
    char name[12];
    snprintf(name, sizeof(name), "dem_%d", j);
    Cbc_addRow(m, name, NS, idx, coef, 'G', demand[j]);
  }

  Cbc_solveLinearProgram(m);
  double obj = Cbc_getObjValue(m);
  Cbc_deleteModel(m);
  return obj;
}

/* Build model column-wise: use Cbc_addCol with row refs after rows exist in solver.
   We add a dummy column first to force rows to be flushed, then delete it. */
static double solveColumnWise(void) {
  Cbc_Model *m = Cbc_newModel();
  Cbc_setObjSense(m, 1.0);
  Cbc_setParam(m, "log", "0");

  /* Add supply rows with a dummy column to force flush */
  /* First add one dummy column so rows can reference it */
  Cbc_addCol(m, "dummy", 0.0, 0.0, 0.0, 0, 0, NULL, NULL);

  for (int i = 0; i < NS; i++) {
    char name[12];
    snprintf(name, sizeof(name), "sup_%d", i);
    int idx = 0; double coef = 0.0;
    Cbc_addRow(m, name, 1, &idx, &coef, 'L', supply[i]);
  }
  for (int j = 0; j < ND; j++) {
    char name[12];
    snprintf(name, sizeof(name), "dem_%d", j);
    int idx = 0; double coef = 0.0;
    Cbc_addRow(m, name, 1, &idx, &coef, 'G', demand[j]);
  }

  /* Now add real columns with row references (rows are flushed) */
  for (int i = 0; i < NS; i++)
    for (int j = 0; j < ND; j++) {
      char name[12];
      snprintf(name, sizeof(name), "x_%d_%d", i, j);
      int rowIdx[2] = {i, NS + j};
      double coef[2] = {1.0, 1.0};
      Cbc_addCol(m, name, 0.0, 1e30, cost[i][j], 0, 2, rowIdx, coef);
    }

  Cbc_solveLinearProgram(m);
  double obj = Cbc_getObjValue(m);
  Cbc_deleteModel(m);
  return obj;
}

/* Build model mixed: add columns and rows interleaved */
static double solveMixed(void) {
  Cbc_Model *m = Cbc_newModel();
  Cbc_setObjSense(m, 1.0);
  Cbc_setParam(m, "log", "0");

  /* Add first source's columns, then its supply row */
  for (int j = 0; j < ND; j++) {
    char name[12];
    snprintf(name, sizeof(name), "x_0_%d", j);
    Cbc_addCol(m, name, 0.0, 1e30, cost[0][j], 0, 0, NULL, NULL);
  }
  {
    int idx[ND] = {0, 1, 2, 3};
    double coef[ND] = {1.0, 1.0, 1.0, 1.0};
    Cbc_addRow(m, "sup_0", ND, idx, coef, 'L', supply[0]);
  }

  /* Add second source's columns, then its supply row */
  for (int j = 0; j < ND; j++) {
    char name[12];
    snprintf(name, sizeof(name), "x_1_%d", j);
    Cbc_addCol(m, name, 0.0, 1e30, cost[1][j], 0, 0, NULL, NULL);
  }
  {
    int idx[ND] = {4, 5, 6, 7};
    double coef[ND] = {1.0, 1.0, 1.0, 1.0};
    Cbc_addRow(m, "sup_1", ND, idx, coef, 'L', supply[1]);
  }

  /* Add third source's columns, then its supply row */
  for (int j = 0; j < ND; j++) {
    char name[12];
    snprintf(name, sizeof(name), "x_2_%d", j);
    Cbc_addCol(m, name, 0.0, 1e30, cost[2][j], 0, 0, NULL, NULL);
  }
  {
    int idx[ND] = {8, 9, 10, 11};
    double coef[ND] = {1.0, 1.0, 1.0, 1.0};
    Cbc_addRow(m, "sup_2", ND, idx, coef, 'L', supply[2]);
  }

  /* Now add demand rows referencing all columns */
  for (int j = 0; j < ND; j++) {
    int idx[NS] = {0*ND + j, 1*ND + j, 2*ND + j};
    double coef[NS] = {1.0, 1.0, 1.0};
    char name[12];
    snprintf(name, sizeof(name), "dem_%d", j);
    Cbc_addRow(m, name, NS, idx, coef, 'G', demand[j]);
  }

  Cbc_solveLinearProgram(m);
  double obj = Cbc_getObjValue(m);
  Cbc_deleteModel(m);
  return obj;
}

int main() {
  printf("=== Model Construction Order Test ===\n");

  double objRow = solveRowWise();
  printf("Row-wise:    obj = %.2f\n", objRow); fflush(stdout);

  double objCol = solveColumnWise();
  printf("Column-wise: obj = %.2f\n", objCol); fflush(stdout);

  double objMix = solveMixed();
  printf("Mixed:       obj = %.2f\n", objMix); fflush(stdout);

  /* All three must give the same optimal */
  assert(fabs(objRow - objCol) < 1e-4);
  assert(fabs(objRow - objMix) < 1e-4);
  printf("All construction orders produce obj = %.2f\n", objRow);
  printf("Model construction order test PASSED!\n");
  return 0;
}
