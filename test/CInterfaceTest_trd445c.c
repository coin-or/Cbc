/* CInterfaceTest_trd445c.c
 *
 * Tests for the trd445c instance (MIP-feasible, solved to optimality).
 *
 * Known values (cross-validated with an independent LP solver):
 *   LP relaxation (pure, all rows):  -158420.57871552
 *   LP without row 0 (R0000001):     -213549.02397840
 *   LP without row 2 (R0000003):     -370800.33438784
 *   LP without rows 0+2:             -370800.33438784
 *   MIP optimal:                     -153419.07883600
 *
 * Instance structure: 96133 rows, 1431 cols (all binary)
 *   96132 packing constraints (R000xxxx, sense ≤ 1)
 *   1 coverage constraint  (atendMin, sense ≥ 3068)
 *
 * Reoptimisation test:
 *   S1 = {row 2}  (R0000003)
 *   S2 = {row 0}  (R0000001)
 *   Remove S1+S2  → LP = -370800.33438784
 *   Add S1        → LP = -213549.02397840
 *   Add S2        → LP = -158420.57871552  (= full LP)
 */

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "Cbc_C_Interface.h"

#ifndef FIXTURE_DIR
#  define FIXTURE_DIR "."
#endif
#define TRD445C_MPS FIXTURE_DIR "/fixtures/trd445c.mps.gz"

/* Reference values (cross-validated) */
#define LP_FULL        -158420.57871552
#define LP_NO_S2       -213549.02397840   /* full LP minus S2={row 0} */
#define LP_NO_S1_S2    -370800.33438784   /* full LP minus S1={row 2} and S2={row 0} */
#define MIP_OPT        -153419.07883600
#define LP_TOL          1e-4
#define MIP_TOL         1e-2

#define S1_ROW  2   /* R0000003 */
#define S2_ROW  0   /* R0000001 */

static int tests_run    = 0;
static int tests_passed = 0;

#define CHECK(cond, msg)                                                \
  do {                                                                  \
    ++tests_run;                                                        \
    if (cond) {                                                         \
      ++tests_passed;                                                   \
    } else {                                                            \
      fprintf(stderr, "FAIL [%s:%d]: %s\n", __FILE__, __LINE__, msg);  \
    }                                                                   \
  } while (0)

static Cbc_Model *load_model(void)
{
  Cbc_Model *m = Cbc_newModel();
  Cbc_setLogLevel(m, 0);
  int rc = Cbc_readMps(m, TRD445C_MPS);
  assert(rc == 0);
  return m;
}

/* Helpers to save/restore a single row -------------------------------- */
typedef struct {
  char   name[128];
  char   sense;
  double rhs;
  int    nz;
  int   *cols;
  double *coefs;
} RowData;

static void save_row(Cbc_Model *m, int r, RowData *d)
{
  Cbc_getRowName(m, r, d->name, sizeof(d->name));
  d->sense = Cbc_getRowSense(m, r);
  d->rhs   = Cbc_getRowRHS(m, r);
  d->nz    = Cbc_getRowNz(m, r);
  d->cols  = malloc(d->nz * sizeof(int));
  d->coefs = malloc(d->nz * sizeof(double));
  const int    *src_cols  = Cbc_getRowIndices(m, r);
  const double *src_coefs = Cbc_getRowCoeffs(m, r);
  memcpy(d->cols,  src_cols,  d->nz * sizeof(int));
  memcpy(d->coefs, src_coefs, d->nz * sizeof(double));
}

static void free_row(RowData *d)
{
  free(d->cols);
  free(d->coefs);
}

/* ------------------------------------------------------------------ */
static void test_model_dimensions(void)
{
  printf("test_model_dimensions\n");
  Cbc_Model *m = load_model();

  CHECK(Cbc_getNumRows(m) == 96133, "expected 96133 rows");
  CHECK(Cbc_getNumCols(m) == 1431,  "expected 1431 cols");

  int nbin = 0, ncont = 0;
  const double *lb = Cbc_getColLower(m);
  const double *ub = Cbc_getColUpper(m);
  for (int j = 0; j < Cbc_getNumCols(m); ++j) {
    if (Cbc_isInteger(m, j) && lb[j] == 0.0 && ub[j] == 1.0) nbin++;
    else ncont++;
  }
  CHECK(nbin  == 1431, "expected 1431 binary vars");
  CHECK(ncont ==    0, "expected 0 non-binary vars");

  Cbc_deleteModel(m);
}

/* ------------------------------------------------------------------ */
static void test_pure_lp(void)
{
  printf("test_pure_lp\n");
  Cbc_Model *m = load_model();

  Cbc_setIntParam(m, INT_PARAM_LP_FAST_PREPROCESS, 0);
  Cbc_setIntParam(m, INT_PARAM_CGRAPH,              0);
  Cbc_setIntParam(m, INT_PARAM_CLIQUE_MERGING,      0);

  int rc = Cbc_solveLinearProgram(m);
  CHECK(rc == 0, "LP solve must succeed");
  CHECK(Cbc_isProvenOptimal(m), "LP must be proven optimal");

  double obj = Cbc_getObjValue(m);
  CHECK(fabs(obj - LP_FULL) < LP_TOL,
        "LP obj must match reference value -158420.579");

  Cbc_deleteModel(m);
}

/* ------------------------------------------------------------------ */
static void test_lp_reoptimization(void)
{
  printf("test_lp_reoptimization\n");
  Cbc_Model *m = load_model();

  Cbc_setIntParam(m, INT_PARAM_LP_FAST_PREPROCESS, 0);
  Cbc_setIntParam(m, INT_PARAM_CGRAPH,              0);
  Cbc_setIntParam(m, INT_PARAM_CLIQUE_MERGING,      0);

  /* Save rows before deletion (indices shift after each delete) */
  RowData s1, s2;
  /* S2_ROW < S1_ROW so save S1 first (higher index unchanged) */
  save_row(m, S1_ROW, &s1);
  save_row(m, S2_ROW, &s2);

  /* Delete S1 and S2 (descending order to keep indices stable) */
  int idx;
  idx = S1_ROW; Cbc_deleteRows(m, 1, &idx);
  idx = S2_ROW; Cbc_deleteRows(m, 1, &idx);

  /* Step 0 — LP without S1 and S2 */
  int rc = Cbc_solveLinearProgram(m);
  CHECK(rc == 0, "step-0 LP solve must succeed");
  CHECK(Cbc_isProvenOptimal(m), "step-0 LP must be optimal");
  double obj0 = Cbc_getObjValue(m);
  CHECK(fabs(obj0 - LP_NO_S1_S2) < LP_TOL,
        "step-0 LP obj must match reference -370800.334");

  /* Step 1 — add S1 back (appended as the last row) */
  Cbc_addRow(m, s1.name, s1.nz, s1.cols, s1.coefs, s1.sense, s1.rhs);
  rc = Cbc_resolve(m);
  CHECK(rc == 0, "step-1 resolve must succeed");
  CHECK(Cbc_isProvenOptimal(m), "step-1 LP must be optimal");
  double obj1 = Cbc_getObjValue(m);
  CHECK(fabs(obj1 - LP_NO_S2) < LP_TOL,
        "step-1 LP obj must match reference -213549.024");

  /* Step 2 — add S2 back */
  Cbc_addRow(m, s2.name, s2.nz, s2.cols, s2.coefs, s2.sense, s2.rhs);
  rc = Cbc_resolve(m);
  CHECK(rc == 0, "step-2 resolve must succeed");
  CHECK(Cbc_isProvenOptimal(m), "step-2 LP must be optimal");
  double obj2 = Cbc_getObjValue(m);
  CHECK(fabs(obj2 - LP_FULL) < LP_TOL,
        "step-2 LP obj must match full LP reference -158420.579");

  free_row(&s1);
  free_row(&s2);
  Cbc_deleteModel(m);
}

/* ------------------------------------------------------------------ */
static void test_mip_optimal(void)
{
  printf("test_mip_optimal\n");
  Cbc_Model *m = load_model();

  Cbc_setDblParam(m, DBL_PARAM_TIME_LIMIT, 120.0);

  Cbc_solve(m);

  CHECK(Cbc_isProvenOptimal(m), "MIP must be proven optimal");

  double obj = Cbc_getObjValue(m);
  CHECK(fabs(obj - MIP_OPT) < MIP_TOL,
        "MIP obj must match reference -153419.079");

  Cbc_deleteModel(m);
}

/* ------------------------------------------------------------------ */
int main(void)
{
  printf("=== trd445c C-interface tests ===\n");

  test_model_dimensions();
  test_pure_lp();
  test_lp_reoptimization();
  test_mip_optimal();

  printf("=== %d / %d tests passed ===\n", tests_passed, tests_run);
  return (tests_passed == tests_run) ? 0 : 1;
}
