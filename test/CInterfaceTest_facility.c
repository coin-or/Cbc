/* Capacitated Facility Location Problem.
   
   5 potential facility locations, 8 customers.
   Decide which facilities to open and how to assign customers.
   Tests: binary + continuous variables, linking constraints. */

#include "Cbc_C_Interface.h"
#include <assert.h>
#include <math.h>
#include <stdio.h>

#define NF 5  /* facilities */
#define NC 8  /* customers */

/* Fixed cost of opening each facility */
static double fixCost[NF] = {100, 120, 110, 90, 130};
/* Capacity of each facility */
static double cap[NF] = {50, 60, 55, 45, 70};
/* Customer demands */
static double dem[NC] = {10, 12, 8, 15, 9, 11, 14, 7};
/* Transport cost[facility][customer] */
static double tcost[NF][NC] = {
  {10, 15, 20, 25, 12, 18, 22, 30},
  {20, 10, 15, 12, 25, 14, 16, 22},
  {15, 20, 10, 18, 15, 20, 12, 14},
  {25, 12, 18, 10, 20, 15, 25, 18},
  {18, 22, 14, 20, 10, 12, 18, 10}
};

int main() {
  printf("=== Facility Location (C interface) ===\n");

  Cbc_Model *m = Cbc_newModel();
  Cbc_setObjSense(m, 1.0);
  Cbc_setParam(m, "log", "0");

  /* y[f]: binary, open facility f. Index: 0..NF-1 */
  for (int f = 0; f < NF; f++) {
    char name[8];
    snprintf(name, sizeof(name), "y_%d", f);
    Cbc_addCol(m, name, 0.0, 1.0, fixCost[f], 1, 0, NULL, NULL);
  }

  /* x[f][c]: fraction of customer c served by facility f. Index: NF + f*NC + c */
  for (int f = 0; f < NF; f++)
    for (int c = 0; c < NC; c++) {
      char name[12];
      snprintf(name, sizeof(name), "x_%d_%d", f, c);
      Cbc_addCol(m, name, 0.0, 1.0, tcost[f][c] * dem[c], 0, 0, NULL, NULL);
    }

  /* Demand: sum_f x[f][c] = 1 for each customer */
  for (int c = 0; c < NC; c++) {
    int idx[NF];
    double coef[NF];
    for (int f = 0; f < NF; f++) {
      idx[f] = NF + f * NC + c;
      coef[f] = 1.0;
    }
    char name[16];
    snprintf(name, sizeof(name), "dem_%d", c);
    Cbc_addRow(m, name, NF, idx, coef, 'E', 1.0);
  }

  /* Capacity: sum_c dem[c]*x[f][c] <= cap[f]*y[f] */
  for (int f = 0; f < NF; f++) {
    int idx[NC + 1];
    double coef[NC + 1];
    for (int c = 0; c < NC; c++) {
      idx[c] = NF + f * NC + c;
      coef[c] = dem[c];
    }
    idx[NC] = f;
    coef[NC] = -cap[f];
    char name[16];
    snprintf(name, sizeof(name), "cap_%d", f);
    Cbc_addRow(m, name, NC + 1, idx, coef, 'L', 0.0);
  }

  Cbc_solve(m);
  assert(Cbc_isProvenOptimal(m));

  double obj = Cbc_getObjValue(m);
  printf("Optimal cost: %.2f\n", obj);

  /* Print open facilities */
  const double *sol = Cbc_getColSolution(m);
  int nOpen = 0;
  for (int f = 0; f < NF; f++) {
    if (sol[f] > 0.5) {
      printf("  Facility %d open (cost=%.0f, cap=%.0f)\n", f, fixCost[f], cap[f]);
      nOpen++;
    }
  }
  printf("  %d facilities opened\n", nOpen);

  /* Verify: must open enough capacity for total demand (86) */
  assert(nOpen >= 2 && nOpen <= 4);
  /* Verify: optimal cost is 1290 for this instance */
  assert(fabs(obj - 1290.0) < 1e-4);

  Cbc_deleteModel(m);
  printf("Facility location test PASSED!\n");
  return 0;
}
