/* Job Shop Scheduling Problem (3x3 instance).
   Ported from python-mip/examples/jssp.py.
   
   3 jobs, 3 machines. Optimal makespan = 7.
   Tests: big-M constraints, binary variables, minimization. */

#include "Cbc_C_Interface.h"
#include <assert.h>
#include <math.h>
#include <stdio.h>

#define N 3  /* jobs */
#define MM 3  /* machines */

static int times[N][MM] = {{1, 2, 2}, {2, 1, 2}, {1, 2, 1}};
static int machines[N][MM] = {{2, 0, 1}, {1, 2, 0}, {2, 1, 0}};

int main() {
  printf("=== Job Shop Scheduling 3x3 (C interface) ===\n");

  int bigM = 0;
  for (int i = 0; i < N; i++)
    for (int j = 0; j < MM; j++)
      bigM += times[i][j];

  Cbc_Model *m = Cbc_newModel();
  Cbc_setObjSense(m, 1.0); /* minimize */
  Cbc_setParam(m, "log", "0");

  /* Variable C (makespan): index 0 */
  Cbc_addCol(m, "C", 0.0, 1e30, 1.0, 0, 0, NULL, NULL);

  /* Variables x[j][i] (start time of job j on machine i): continuous
     index = 1 + j*MM + i */
  for (int j = 0; j < N; j++)
    for (int i = 0; i < MM; i++) {
      char name[16];
      snprintf(name, sizeof(name), "x_%d_%d", j + 1, i + 1);
      Cbc_addCol(m, name, 0.0, 1e30, 0.0, 0, 0, NULL, NULL);
    }
  /* x[j][i] index macro */
  #define X(j,i) (1 + (j)*MM + (i))

  /* Variables y[j][k][i] (binary precedence): 
     index = 1 + N*MM + j*N*MM + k*MM + i */
  int yBase = 1 + N * MM;
  for (int j = 0; j < N; j++)
    for (int k = 0; k < N; k++)
      for (int i = 0; i < MM; i++) {
        char name[16];
        snprintf(name, sizeof(name), "y_%d_%d_%d", j+1, k+1, i+1);
        Cbc_addCol(m, name, 0.0, 1.0, 0.0, 1, 0, NULL, NULL);
      }

  /* Precedence within job: x[j][machines[j][i]] - x[j][machines[j][i-1]] >= times[j][machines[j][i-1]] */
  for (int j = 0; j < N; j++) {
    for (int i = 1; i < MM; i++) {
      int idx[2] = {1 + j*MM + machines[j][i], 1 + j*MM + machines[j][i-1]};
      double coef[2] = {1.0, -1.0};
      char name[32];
      snprintf(name, sizeof(name), "prec_%d_%d", j+1, i);
      Cbc_addRow(m, name, 2, idx, coef, 'G', (double)times[j][machines[j][i-1]]);
    }
  }

  /* No-overlap: for j!=k, machine i:
     x[j][i] - x[k][i] + M*y[j][k][i] >= times[k][i]
     -x[j][i] + x[k][i] - M*y[j][k][i] >= times[j][i] - M */
  for (int j = 0; j < N; j++) {
    for (int k = 0; k < N; k++) {
      if (k == j) continue;
      for (int i = 0; i < MM; i++) {
        int yVar = yBase + j*N*MM + k*MM + i;
        /* x[j][i] - x[k][i] + M*y >= times[k][i] */
        {
          int idx[3] = {1 + j*MM + i, 1 + k*MM + i, yVar};
          double coef[3] = {1.0, -1.0, (double)bigM};
          char name[32];
          snprintf(name, sizeof(name), "ov1_%d_%d_%d", j+1, k+1, i+1);
          Cbc_addRow(m, name, 3, idx, coef, 'G', (double)times[k][i]);
        }
        /* -x[j][i] + x[k][i] - M*y >= times[j][i] - M */
        {
          int idx[3] = {1 + j*MM + i, 1 + k*MM + i, yVar};
          double coef[3] = {-1.0, 1.0, -(double)bigM};
          char name[32];
          snprintf(name, sizeof(name), "ov2_%d_%d_%d", j+1, k+1, i+1);
          Cbc_addRow(m, name, 3, idx, coef, 'G', (double)(times[j][i] - bigM));
        }
      }
    }
  }

  /* Makespan: C - x[j][machines[j][m-1]] >= times[j][machines[j][m-1]] */
  for (int j = 0; j < N; j++) {
    int lastMach = machines[j][MM - 1];
    int idx[2] = {0, 1 + j*MM + lastMach};
    double coef[2] = {1.0, -1.0};
    char name[16];
    snprintf(name, sizeof(name), "mksp_%d", j+1);
    Cbc_addRow(m, name, 2, idx, coef, 'G', (double)times[j][lastMach]);
  }

  Cbc_solve(m);
  assert(Cbc_isProvenOptimal(m));

  /* Check feasibility of best solution */
  const double *bestSol = Cbc_getColSolution(m);
  double maxViolRow = 0.0; int rowIdx = -1;
  double maxViolCol = 0.0; int colIdx = -1;

  if (!Cbc_checkFeasibility(m, bestSol, &maxViolRow, &rowIdx, &maxViolCol, &colIdx)) {
    printf("FAIL: best solution infeasible  "
           "maxViolRow=%.2e (row %d)  maxViolCol=%.2e (col %d)\n",
           maxViolRow, rowIdx, maxViolCol, colIdx);
    return 1;
  }

  /* Check feasibility of ALL solutions in the pool */
  int nSol = Cbc_numberSavedSolutions(m);
  for (int s = 0; s < nSol; s++) {
    const double *sol = Cbc_savedSolution(m, s);
    double solObj = Cbc_savedSolutionObj(m, s);
    maxViolRow = 0.0; rowIdx = -1;
    maxViolCol = 0.0; colIdx = -1;

    if (!Cbc_checkFeasibility(m, sol, &maxViolRow, &rowIdx, &maxViolCol, &colIdx)) {
      printf("FAIL: pool solution %d/%d infeasible (obj=%.0f)  "
             "maxViolRow=%.2e (row %d)  maxViolCol=%.2e (col %d)\n",
             s+1, nSol, solObj, maxViolRow, rowIdx, maxViolCol, colIdx);
      return 1;
    }
  }

  double makespan = Cbc_getObjValue(m);
  printf("Optimal makespan: %.0f (expected 7)\n", makespan);
  assert(fabs(makespan - 7.0) < 1e-4);

  Cbc_deleteModel(m);
  printf("JSSP test PASSED!\n");
  return 0;
}
