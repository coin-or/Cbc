/* CInterfaceTest_yue20013.c
 *
 * Tests for the yue20013.1.150 instance (MIP-feasible, solved to optimality).
 *
 * Known values (cross-validated with an independent LP solver):
 *   LP relaxation (pure, all rows):   0.00000000000000 (trivially zero)
 *   LP without rows 100+200:          0.00000000000000
 *   MIP optimal:                     29.0 (integer)
 *
 * Instance structure: 3684 rows, 874 cols
 *   824 binary variables, 50 general integer, 0 continuous
 *
 * Reoptimisation test:
 *   S1 = {row 200} (R0000201)
 *   S2 = {row 100} (R0000101)
 *   The LP relaxation is zero throughout; the test validates that the
 *   delete/add-row/resolve API works correctly even when the objective
 *   does not change.
 */

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "Cbc_C_Interface.h"

#include "test_utils.h"
#define YUE20013_MPS fixture_path("yue20013.1.150.mps.gz")

/* Reference values (cross-validated) */
#define LP_FULL       0.0
#define LP_NO_S1_S2   0.0
#define LP_NO_S2      0.0
#define MIP_OPT       29.0
#define LP_TOL        1e-4
#define MIP_TOL       1e-2

#define S1_ROW  200   /* R0000201 */
#define S2_ROW  100   /* R0000101 */

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
  int rc = Cbc_readMps(m, YUE20013_MPS);
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

  CHECK(Cbc_getNumRows(m) == 3684, "expected 3684 rows");
  CHECK(Cbc_getNumCols(m) ==  874, "expected 874 cols");

  int nbin = 0, nint = 0, ncont = 0;
  const double *lb = Cbc_getColLower(m);
  const double *ub = Cbc_getColUpper(m);
  for (int j = 0; j < Cbc_getNumCols(m); ++j) {
    if (Cbc_isInteger(m, j)) {
      if (lb[j] == 0.0 && ub[j] == 1.0) nbin++;
      else nint++;
    } else {
      ncont++;
    }
  }
  CHECK(nbin  == 824, "expected 824 binary vars");
  CHECK(nint  ==  50, "expected 50 general integer vars");
  CHECK(ncont ==   0, "expected 0 continuous vars");

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
        "LP obj must match reference value 0.0");

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

  /* Save rows before deletion (S1_ROW > S2_ROW so save S1 first) */
  RowData s1, s2;
  save_row(m, S1_ROW, &s1);
  save_row(m, S2_ROW, &s2);

  /* Delete in descending order to keep indices stable */
  int idx;
  idx = S1_ROW; Cbc_deleteRows(m, 1, &idx);
  idx = S2_ROW; Cbc_deleteRows(m, 1, &idx);

  /* Step 0 — LP without S1 and S2 */
  int rc = Cbc_solveLinearProgram(m);
  CHECK(rc == 0, "step-0 LP solve must succeed");
  CHECK(Cbc_isProvenOptimal(m), "step-0 LP must be optimal");
  double obj0 = Cbc_getObjValue(m);
  CHECK(fabs(obj0 - LP_NO_S1_S2) < LP_TOL,
        "step-0 LP obj must match reference 0.0");

  /* Step 1 — add S1 back (appended as last row) */
  Cbc_addRow(m, s1.name, s1.nz, s1.cols, s1.coefs, s1.sense, s1.rhs);
  rc = Cbc_resolve(m);
  CHECK(rc == 0, "step-1 resolve must succeed");
  CHECK(Cbc_isProvenOptimal(m), "step-1 LP must be optimal");
  double obj1 = Cbc_getObjValue(m);
  CHECK(fabs(obj1 - LP_NO_S2) < LP_TOL,
        "step-1 LP obj must match reference 0.0");

  /* Step 2 — add S2 back */
  Cbc_addRow(m, s2.name, s2.nz, s2.cols, s2.coefs, s2.sense, s2.rhs);
  rc = Cbc_resolve(m);
  CHECK(rc == 0, "step-2 resolve must succeed");
  CHECK(Cbc_isProvenOptimal(m), "step-2 LP must be optimal");
  double obj2 = Cbc_getObjValue(m);
  CHECK(fabs(obj2 - LP_FULL) < LP_TOL,
        "step-2 LP obj must match full LP reference 0.0");

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
  /* Disable ratio gap so the solver finds the exact integer optimal,
     regardless of platform-specific B&B ordering differences. */
  Cbc_setDblParam(m, DBL_PARAM_GAP_RATIO, 0.0);

  Cbc_solve(m);

  CHECK(Cbc_isProvenOptimal(m), "MIP must be proven optimal");

  double obj = Cbc_getObjValue(m);
  CHECK(fabs(obj - MIP_OPT) < MIP_TOL,
        "MIP obj must match reference 29.0");

  /* Validate every solution in the pool: feasibility + obj on correct side */
  int fails = validate_all_saved_solutions(m, MIP_OPT, MIP_TOL, "yue20013");
  CHECK(fails == 0, "all saved solutions are feasible with obj on correct side of MIP_OPT");

  Cbc_deleteModel(m);
}

/* ------------------------------------------------------------------ */
int main(void)
{
  printf("=== yue20013.1.150 C-interface tests ===\n");

  test_model_dimensions();
  test_pure_lp();
  test_lp_reoptimization();
  test_mip_optimal();

  printf("=== %d / %d tests passed ===\n", tests_passed, tests_run);
  return (tests_passed == tests_run) ? 0 : 1;
}
