#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <assert.h>

#include "Cbc_C_Interface.h"

int main(int argc, char *argv[]) {
  printf("=== GF(2) Parity Presolve Test ===\n");

  char lp_path_buf[4096];
  const char *lp_path = lp_path_buf;
  if (argc > 1) {
    snprintf(lp_path_buf, sizeof(lp_path_buf), "%s/parity_test.lp", argv[1]);
  } else {
    lp_path = "fixtures/parity_test.lp";
  }

  // 1. Solve with preprocessing OFF (baseline)
  printf("Solving with preprocessing OFF...\n");
  Cbc_Model *m_off = Cbc_newModel();
  Cbc_setLogLevel(m_off, 1);
  int rc_off = Cbc_readLp(m_off, lp_path);
  if (rc_off != 0) {
    fprintf(stderr, "ERROR: Failed to read %s\n", lp_path);
    Cbc_deleteModel(m_off);
    return 1;
  }
  Cbc_setParam(m_off, "preprocess", "off");
  Cbc_solve(m_off);

  int nsols_off = Cbc_numberSavedSolutions(m_off);
  printf("Baseline solutions: %d\n", nsols_off);
  if (nsols_off <= 0) {
    fprintf(stderr, "FAIL: No solution found with preprocessing OFF\n");
    Cbc_deleteModel(m_off);
    return 1;
  }
  double obj_off = Cbc_getObjValue(m_off);
  printf("Baseline objective: %g\n", obj_off);

  // Validate feasibility of baseline solution
  const double *sol_off = Cbc_getColSolution(m_off);
  double maxViolRow_off = 0.0; int rowIdx_off = -1;
  double maxViolCol_off = 0.0; int colIdx_off = -1;
  if (!Cbc_checkFeasibility(m_off, sol_off, &maxViolRow_off, &rowIdx_off, &maxViolCol_off, &colIdx_off)) {
    fprintf(stderr, "FAIL: Baseline solution is infeasible\n");
    Cbc_deleteModel(m_off);
    return 1;
  }

  // 2. Solve with preprocessing ON (enables parity presolve)
  printf("Solving with preprocessing ON...\n");
  Cbc_Model *m_on = Cbc_newModel();
  Cbc_setLogLevel(m_on, 1);
  int rc_on = Cbc_readLp(m_on, lp_path);
  if (rc_on != 0) {
    fprintf(stderr, "ERROR: Failed to read %s\n", lp_path);
    Cbc_deleteModel(m_off);
    Cbc_deleteModel(m_on);
    return 1;
  }
  Cbc_setParam(m_on, "preprocess", "on");
  Cbc_solve(m_on);

  int nsols_on = Cbc_numberSavedSolutions(m_on);
  printf("Presolve solutions: %d\n", nsols_on);
  if (nsols_on <= 0) {
    fprintf(stderr, "FAIL: No solution found with preprocessing ON\n");
    Cbc_deleteModel(m_off);
    Cbc_deleteModel(m_on);
    return 1;
  }
  double obj_on = Cbc_getObjValue(m_on);
  printf("Presolve objective: %g\n", obj_on);

  // Validate feasibility of presolved/postsolved solution
  const double *sol_on = Cbc_getColSolution(m_on);
  double maxViolRow_on = 0.0; int rowIdx_on = -1;
  double maxViolCol_on = 0.0; int colIdx_on = -1;
  if (!Cbc_checkFeasibility(m_on, sol_on, &maxViolRow_on, &rowIdx_on, &maxViolCol_on, &colIdx_on)) {
    fprintf(stderr, "FAIL: Postsolved solution is infeasible (maxViolRow=%g, rowIdx=%d)\n", maxViolRow_on, rowIdx_on);
    Cbc_deleteModel(m_off);
    Cbc_deleteModel(m_on);
    return 1;
  }

  // Compare objective values
  if (fabs(obj_off - obj_on) > 1e-4) {
    fprintf(stderr, "FAIL: Objectives do not match! Baseline: %g, Presolved: %g\n", obj_off, obj_on);
    Cbc_deleteModel(m_off);
    Cbc_deleteModel(m_on);
    return 1;
  }

  printf("Objective values match exactly: %g\n", obj_on);
  printf("GF(2) Parity Presolve Test PASSED!\n");

  Cbc_deleteModel(m_off);
  Cbc_deleteModel(m_on);
  return 0;
}
