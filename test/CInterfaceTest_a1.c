/* CInterfaceTest_a1.c
 *
 * Tests for the project-scheduling instance A-1 (MIPLIB 2017).
 *
 * Known values (cross-validated with an independent LP solver):
 *   LP relaxation (full, 1015 rows):       3,100,022
 *   LP without lastP[0] and lastP[1]:      3,100,000   (1013 rows)
 *   LP + lastP[0] re-added (1014 rows):    3,100,009   (+9 from base)
 *   LP + lastP[0]+lastP[1] re-added:       3,100,022   (+13 further)
 *   MIP optimum:                           3,200,022
 *
 * The two lastP constraints (rows 1013 and 1014) are used as subsets
 * S1 = {lastP(0)} and S2 = {lastP(1)}.  Adding S1 alone tightens the
 * LP bound by 9; adding S2 further tightens it by 13.
 */

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "Cbc_C_Interface.h"

#include "test_utils.h"
#define A1_MPS fixture_path("A-1.mps.gz")

/* Known solution values */
#define A1_LP_FULL         3100022.0
#define A1_LP_NO_LASTP     3100000.0
#define A1_LP_WITH_S1      3100009.0
#define A1_MIP_OPT         3200022.0
#define LP_TOL             1e-3
#define MIP_TOL            1e-5

/* Row indices of the two lastP constraints in the original model */
#define S1_ROW  1013   /* lastP(0) */
#define S2_ROW  1014   /* lastP(1) */

static int tests_run    = 0;
static int tests_passed = 0;

#define CHECK(cond, msg)                                              \
  do {                                                               \
    ++tests_run;                                                     \
    if (cond) {                                                      \
      ++tests_passed;                                                \
    } else {                                                         \
      fprintf(stderr, "FAIL [%s:%d] %s\n", __FILE__, __LINE__, msg);\
    }                                                                \
  } while (0)

/* ------------------------------------------------------------------ */
/* helpers                                                             */
/* ------------------------------------------------------------------ */

static Cbc_Model *load_model(int log_level)
{
  Cbc_Model *m = Cbc_newModel();
  Cbc_setLogLevel(m, log_level);
  int rc = Cbc_readMps(m, A1_MPS);
  if (rc != 0) {
    fprintf(stderr, "ERROR: could not read %s\n", A1_MPS);
    Cbc_deleteModel(m);
    return NULL;
  }
  return m;
}

/* Saved row data for later re-insertion */
typedef struct {
  char   name[64];
  int    nz;
  int   *idx;
  double *coef;
  char   sense;
  double rhs;
} RowData;

static RowData save_row(Cbc_Model *m, int r)
{
  RowData rd;
  Cbc_getRowName(m, r, rd.name, sizeof(rd.name));
  rd.nz    = Cbc_getRowNz(m, r);
  rd.sense = Cbc_getRowSense(m, r);
  rd.rhs   = Cbc_getRowRHS(m, r);
  rd.idx   = (int*)malloc(rd.nz * sizeof(int));
  rd.coef  = (double*)malloc(rd.nz * sizeof(double));
  const int    *ri = Cbc_getRowIndices(m, r);
  const double *rc = Cbc_getRowCoeffs(m, r);
  for (int j = 0; j < rd.nz; j++) {
    rd.idx[j]  = ri[j];
    rd.coef[j] = rc[j];
  }
  return rd;
}

static void free_row(RowData *rd)
{
  free(rd->idx);
  free(rd->coef);
}

/* ------------------------------------------------------------------ */
/* test_model_dimensions                                               */
/* ------------------------------------------------------------------ */

static void test_model_dimensions(void)
{
  Cbc_Model *m = load_model(0);
  if (!m) return;

  CHECK(Cbc_getNumRows(m) == 1015,  "A-1: 1015 rows");
  CHECK(Cbc_getNumCols(m) == 18598, "A-1: 18598 columns");

  Cbc_deleteModel(m);
}

/* ------------------------------------------------------------------ */
/* test_pure_lp                                                        */
/* ------------------------------------------------------------------ */

static void test_pure_lp(void)
{
  Cbc_Model *m = load_model(0);
  if (!m) return;

  Cbc_setIntParam(m, INT_PARAM_LP_FAST_PREPROCESS, 0);
  Cbc_setIntParam(m, INT_PARAM_CGRAPH, 0);
  Cbc_setIntParam(m, INT_PARAM_CLIQUE_MERGING, 0);

  int rc = Cbc_solveLinearProgram(m);
  CHECK(rc == 0, "pure LP: solve status optimal");

  double obj = Cbc_getObjValue(m);
  CHECK(fabs(obj - A1_LP_FULL) < LP_TOL,
        "pure LP: objective matches known LP relaxation");
  CHECK(obj <= A1_MIP_OPT + MIP_TOL,
        "pure LP: LP obj <= MIP opt (valid lower bound)");

  Cbc_deleteModel(m);
}

/* ------------------------------------------------------------------ */
/* test_lp_reoptimization                                              */
/*                                                                     */
/* Remove S1={lastP(0)} and S2={lastP(1)}, solve LP, then add S1      */
/* back and resolve, then add S2 back and resolve.  The LP bound       */
/* should increase at each step.                                       */
/* ------------------------------------------------------------------ */

static void test_lp_reoptimization(void)
{
  Cbc_Model *m = load_model(0);
  if (!m) return;

  Cbc_setIntParam(m, INT_PARAM_LP_FAST_PREPROCESS, 0);
  Cbc_setIntParam(m, INT_PARAM_CGRAPH, 0);
  Cbc_setIntParam(m, INT_PARAM_CLIQUE_MERGING, 0);

  /* Save S1 and S2 before deletion (copy internal pointers now). */
  RowData s1 = save_row(m, S1_ROW);
  RowData s2 = save_row(m, S2_ROW);

  /* Delete S2 first (higher index), then S1 to preserve index stability. */
  int del_s2[1] = { S2_ROW };
  int del_s1[1] = { S1_ROW };
  Cbc_deleteRows(m, 1, del_s2);
  Cbc_deleteRows(m, 1, del_s1);

  CHECK(Cbc_getNumRows(m) == 1013, "reopt: 1013 rows after removing S1+S2");

  int rc = Cbc_solveLinearProgram(m);
  CHECK(rc == 0, "reopt: base LP optimal");
  double obj_base = Cbc_getObjValue(m);
  CHECK(fabs(obj_base - A1_LP_NO_LASTP) < LP_TOL,
        "reopt: base LP bound matches reference (no lastP)");

  /* Add S1 (lastP(0)) back and warm-start resolve. */
  Cbc_addRow(m, s1.name, s1.nz, s1.idx, s1.coef, s1.sense, s1.rhs);
  CHECK(Cbc_getNumRows(m) == 1014, "reopt: 1014 rows after adding S1");

  rc = Cbc_resolve(m);
  CHECK(rc == 0, "reopt: LP with S1 optimal after resolve");
  double obj_s1 = Cbc_getObjValue(m);
  CHECK(fabs(obj_s1 - A1_LP_WITH_S1) < LP_TOL,
        "reopt: LP bound with S1 matches reference (+9)");
  CHECK(obj_s1 >= obj_base - LP_TOL,
        "reopt: bound non-decreasing after adding S1");

  /* Add S2 (lastP(1)) back and warm-start resolve. */
  Cbc_addRow(m, s2.name, s2.nz, s2.idx, s2.coef, s2.sense, s2.rhs);
  CHECK(Cbc_getNumRows(m) == 1015, "reopt: 1015 rows after adding S2");

  rc = Cbc_resolve(m);
  CHECK(rc == 0, "reopt: LP with S1+S2 optimal after resolve");
  double obj_s1s2 = Cbc_getObjValue(m);
  CHECK(fabs(obj_s1s2 - A1_LP_FULL) < LP_TOL,
        "reopt: LP bound with S1+S2 matches full LP reference");
  CHECK(obj_s1s2 >= obj_s1 - LP_TOL,
        "reopt: bound non-decreasing after adding S2");

  free_row(&s1);
  free_row(&s2);
  Cbc_deleteModel(m);
}

/* ------------------------------------------------------------------ */
/* test_mip                                                            */
/* ------------------------------------------------------------------ */

static void test_mip(void)
{
  Cbc_Model *m = load_model(0);
  if (!m) return;

  /* Node limit instead of time limit for determinism (3444 nodes at opt with seed=1) */
  Cbc_setIntParam(m, INT_PARAM_MAX_NODES, 5000);
  int rc = Cbc_solve(m);
  (void)rc;

  int is_proven_optimal = Cbc_isProvenOptimal(m);
  double obj = Cbc_getObjValue(m);

  /* Feasibility is mandatory — always check it */
  CHECK(Cbc_bestSolution(m) != NULL, "MIP: best solution exists");

  /* Only require correct objective if solver claims optimality */
  if (is_proven_optimal) {
    CHECK(fabs(obj - A1_MIP_OPT) < MIP_TOL,
          "MIP: optimal value matches known optimum");
  }

  int fails = validate_all_saved_solutions(m, A1_MIP_OPT, MIP_TOL, "A-1");
  CHECK(fails == 0, "MIP: all saved solutions are feasible with correct objective");

  Cbc_deleteModel(m);
}

/* ------------------------------------------------------------------ */
/* main                                                                */
/* ------------------------------------------------------------------ */

int main(void)
{
  printf("=== A-1 project-scheduling C interface tests ===\n");

  test_model_dimensions();
  test_pure_lp();
  test_lp_reoptimization();
  test_mip();

  printf("\n%d/%d tests passed.\n", tests_passed, tests_run);
  return (tests_passed == tests_run) ? 0 : 1;
}
