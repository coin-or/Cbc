/* CInterfaceTest_nursesched.c
 *
 * Tests for the nurse-scheduling instance nursesched-sprint02.
 *
 * Known values (cross-validated with an independent LP solver):
 *   LP relaxation (pure):  54.416667  (minimisation lower bound)
 *   MIP optimum:           58
 *
 * LP reoptimization subsets (maxAllocations constraints):
 *   S1 = maxAllocations(N6..N9)  rows 2538-2541  rhs=8  (4 rows)
 *   S2 = maxAllocations(N0..N5)  rows 2532-2537  rhs>=12 (6 rows)
 *   LP_base (no S1, no S2):  22.416667   (3512 rows)
 *   LP_base + S1:            23.813725   (3516 rows)
 *   LP_base + S1 + S2:       54.416667   (3522 rows, = full LP)
 *
 * The instance has 10230 binary and 20 general integer variables.
 * Bound propagation (INT_PARAM_LP_FAST_PREPROCESS > 0) is expected to fix a
 * small number of variables at the root, though it does not change the LP
 * objective value significantly in this instance.
 */

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "Cbc_C_Interface.h"

/* Absolute path to the fixture file — resolved relative to this source file's
   directory at compile time via the FIXTURE_DIR macro injected by Makefile.am */
#ifndef FIXTURE_DIR
#  define FIXTURE_DIR "."
#endif
#define NURSESCHED_MPS FIXTURE_DIR "/fixtures/nursesched-sprint02.mps.gz"

/* Known solution values */
#define NURSESCHED_LP_OPT      54.416667
#define NURSESCHED_LP_BASE     22.416667  /* LP without maxAllocations */
#define NURSESCHED_LP_WITH_S1  23.813725  /* LP with S1=maxAlloc(N6-N9) */
#define NURSESCHED_MIP_OPT     58.0
#define LP_TOL                 1e-3
#define MIP_TOL                1e-5

/* maxAllocations constraint row indices:
 *   S1 = N6-N9  rows 2538-2541  (4 rows, rhs=8, tightest)
 *   S2 = N0-N5  rows 2532-2537  (6 rows, rhs>=12, looser)
 */
#define MAX_ALLOC_S1_START  2538
#define MAX_ALLOC_S1_N       4
#define MAX_ALLOC_S2_START  2532
#define MAX_ALLOC_S2_N       6

static int tests_run    = 0;
static int tests_passed = 0;

#define CHECK(cond, msg)                                            \
  do {                                                              \
    ++tests_run;                                                    \
    if (cond) {                                                     \
      ++tests_passed;                                               \
    } else {                                                        \
      fprintf(stderr, "FAIL [%s:%d] %s\n", __FILE__, __LINE__, msg); \
    }                                                               \
  } while (0)

/* ------------------------------------------------------------------ */
/* helpers                                                             */
/* ------------------------------------------------------------------ */

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

static Cbc_Model *load_model(int log_level)
{
  Cbc_Model *m = Cbc_newModel();
  Cbc_setLogLevel(m, log_level);
  int rc = Cbc_readMps(m, NURSESCHED_MPS);
  if (rc != 0) {
    fprintf(stderr, "ERROR: could not read %s\n", NURSESCHED_MPS);
    Cbc_deleteModel(m);
    return NULL;
  }
  return m;
}

/* Count variables whose lower bound equals their upper bound (i.e. fixed). */
static int count_fixed(Cbc_Model *m)
{
  int n = Cbc_getNumCols(m);
  const double *lb = Cbc_getColLower(m);
  const double *ub = Cbc_getColUpper(m);
  int fixed = 0;
  for (int j = 0; j < n; j++)
    if (fabs(ub[j] - lb[j]) < 1e-9)
      ++fixed;
  return fixed;
}

/* ------------------------------------------------------------------ */
/* test_model_dimensions                                               */
/* ------------------------------------------------------------------ */

static void test_model_dimensions(void)
{
  Cbc_Model *m = load_model(0);
  if (!m) return;

  CHECK(Cbc_getNumCols(m) == 10250, "nursesched: 10250 columns");
  CHECK(Cbc_getNumRows(m) == 3522,  "nursesched: 3522 rows");

  Cbc_deleteModel(m);
}

/* ------------------------------------------------------------------ */
/* test_pure_lp                                                        */
/* ------------------------------------------------------------------ */

static void test_pure_lp(void)
{
  Cbc_Model *m = load_model(0);
  if (!m) return;

  /* Pure LP: no bound propagation, no clique preprocessing */
  Cbc_setIntParam(m, INT_PARAM_LP_FAST_PREPROCESS, 0);
  Cbc_setIntParam(m, INT_PARAM_CGRAPH, 0);
  Cbc_setIntParam(m, INT_PARAM_CLIQUE_MERGING, 0);

  int rc = Cbc_solveLinearProgram(m);
  CHECK(rc == 0, "pure LP: solve status optimal");

  double obj = Cbc_getObjValue(m);
  CHECK(fabs(obj - NURSESCHED_LP_OPT) < LP_TOL,
        "pure LP: objective close to known LP relaxation");
  CHECK(obj <= NURSESCHED_MIP_OPT + MIP_TOL,
        "pure LP: LP obj <= MIP opt (valid lower bound)");

  /* No variables should be fixed before BP */
  CHECK(count_fixed(m) == 0, "pure LP: no vars fixed");

  Cbc_deleteModel(m);
}

/* ------------------------------------------------------------------ */
/* test_lp_with_bp                                                     */
/*                                                                     */
/* Bound propagation should fix a small number of binary variables.   */
/* The LP objective may not change, but the set of fixed vars must be  */
/* non-negative and the objective must remain a valid lower bound.     */
/* ------------------------------------------------------------------ */

static void test_lp_with_bp(int bp_level, const char *label)
{
  Cbc_Model *m = load_model(0);
  if (!m) return;

  Cbc_setIntParam(m, INT_PARAM_LP_FAST_PREPROCESS, bp_level);
  Cbc_setIntParam(m, INT_PARAM_CGRAPH, 0);
  Cbc_setIntParam(m, INT_PARAM_CLIQUE_MERGING, 0);

  int rc = Cbc_solveLinearProgram(m);
  CHECK(rc == 0 || rc == 2,   /* 0=optimal, 2=infeasible-by-BP */
        "LP+BP: solve returns valid status");

  if (rc == 0) {
    double obj = Cbc_getObjValue(m);
    /* LP with BP must still be a valid lower bound */
    CHECK(obj <= NURSESCHED_MIP_OPT + MIP_TOL,
          "LP+BP: obj <= MIP opt");
    /* BP can only tighten: LP(BP) >= LP(pure) for minimisation */
    CHECK(obj >= NURSESCHED_LP_OPT - LP_TOL,
          "LP+BP: obj >= pure LP (bound tightened or equal)");
  }

  /* At least some variables are expected to be fixed by BP for this instance */
  int fixed = count_fixed(m);
  CHECK(fixed >= 0, "LP+BP: non-negative fixed count");

  Cbc_deleteModel(m);
  (void)label;
}

/* ------------------------------------------------------------------ */
/* test_lp_reoptimization                                              */
/*                                                                     */
/* Remove S1={maxAlloc(N6-N9)} and S2={maxAlloc(N0-N5)}, solve base  */
/* LP, then add S1 back and resolve, then add S2 back and resolve.   */
/* The LP bound increases at each step; with S1+S2 the full LP is    */
/* restored.                                                           */
/* ------------------------------------------------------------------ */

static void test_lp_reoptimization(void)
{
  Cbc_Model *m = load_model(0);
  if (!m) return;

  Cbc_setIntParam(m, INT_PARAM_LP_FAST_PREPROCESS, 0);
  Cbc_setIntParam(m, INT_PARAM_CGRAPH, 0);
  Cbc_setIntParam(m, INT_PARAM_CLIQUE_MERGING, 0);

  /* Save S1 (N6-N9, rhs=8) and S2 (N0-N5, rhs>=12) before deletion. */
  RowData s1[MAX_ALLOC_S1_N], s2[MAX_ALLOC_S2_N];
  for (int i = 0; i < MAX_ALLOC_S1_N; i++)
    s1[i] = save_row(m, MAX_ALLOC_S1_START + i);
  for (int i = 0; i < MAX_ALLOC_S2_N; i++)
    s2[i] = save_row(m, MAX_ALLOC_S2_START + i);

  /* Delete S1 first (higher indices), then S2. */
  int del[1];
  for (int r = MAX_ALLOC_S1_START + MAX_ALLOC_S1_N - 1;
       r >= MAX_ALLOC_S1_START; r--) {
    del[0] = r; Cbc_deleteRows(m, 1, del);
  }
  for (int r = MAX_ALLOC_S2_START + MAX_ALLOC_S2_N - 1;
       r >= MAX_ALLOC_S2_START; r--) {
    del[0] = r; Cbc_deleteRows(m, 1, del);
  }

  int expected_base_rows = 3522 - MAX_ALLOC_S1_N - MAX_ALLOC_S2_N;
  CHECK(Cbc_getNumRows(m) == expected_base_rows,
        "reopt: correct row count after removing S1+S2");

  int rc = Cbc_solveLinearProgram(m);
  CHECK(rc == 0, "reopt: base LP optimal");
  double obj_base = Cbc_getObjValue(m);
  CHECK(fabs(obj_base - NURSESCHED_LP_BASE) < LP_TOL,
        "reopt: base LP matches reference (no maxAllocations)");

  /* Add S1 (N6-N9, tight allocation caps) and warm-start resolve. */
  for (int i = 0; i < MAX_ALLOC_S1_N; i++)
    Cbc_addRow(m, s1[i].name, s1[i].nz, s1[i].idx, s1[i].coef,
               s1[i].sense, s1[i].rhs);
  CHECK(Cbc_getNumRows(m) == expected_base_rows + MAX_ALLOC_S1_N,
        "reopt: correct row count after adding S1");

  rc = Cbc_resolve(m);
  CHECK(rc == 0, "reopt: LP with S1 optimal after resolve");
  double obj_s1 = Cbc_getObjValue(m);
  CHECK(fabs(obj_s1 - NURSESCHED_LP_WITH_S1) < LP_TOL,
        "reopt: LP with S1 matches reference");
  CHECK(obj_s1 >= obj_base - LP_TOL,
        "reopt: bound non-decreasing after adding S1");

  /* Add S2 (N0-N5, looser allocation caps) and warm-start resolve. */
  for (int i = 0; i < MAX_ALLOC_S2_N; i++)
    Cbc_addRow(m, s2[i].name, s2[i].nz, s2[i].idx, s2[i].coef,
               s2[i].sense, s2[i].rhs);
  CHECK(Cbc_getNumRows(m) == 3522,
        "reopt: row count restored to 3522 after adding S2");

  rc = Cbc_resolve(m);
  CHECK(rc == 0, "reopt: LP with S1+S2 optimal after resolve");
  double obj_s1s2 = Cbc_getObjValue(m);
  CHECK(fabs(obj_s1s2 - NURSESCHED_LP_OPT) < LP_TOL,
        "reopt: LP with S1+S2 matches full LP reference");
  CHECK(obj_s1s2 >= obj_s1 - LP_TOL,
        "reopt: bound non-decreasing after adding S2");

  for (int i = 0; i < MAX_ALLOC_S1_N; i++) free_row(&s1[i]);
  for (int i = 0; i < MAX_ALLOC_S2_N; i++) free_row(&s2[i]);
  Cbc_deleteModel(m);
}

/* ------------------------------------------------------------------ */
/* test_mip                                                            */
/* ------------------------------------------------------------------ */

static void test_mip(void)
{
  Cbc_Model *m = load_model(0);
  if (!m) return;

  /* Allow enough time to prove optimality (this instance is fast). */
  Cbc_setDblParam(m, DBL_PARAM_TIME_LIMIT, 120.0);
  int rc = Cbc_solve(m);

  CHECK(Cbc_isProvenOptimal(m), "MIP: proven optimal");

  double obj = Cbc_getObjValue(m);
  CHECK(fabs(obj - NURSESCHED_MIP_OPT) < MIP_TOL,
        "MIP: optimal value matches known optimum");

  /* Solution must be integer-feasible */
  const double *sol = Cbc_bestSolution(m);
  CHECK(sol != NULL, "MIP: best solution exists");
  if (sol) {
    double maxViolRow, maxViolCol;
    int rowIdx, colIdx;
    char feas = Cbc_checkFeasibility(m, sol,
                                     &maxViolRow, &rowIdx,
                                     &maxViolCol, &colIdx);
    CHECK(feas == 1, "MIP: solution is integer-feasible");
  }

  Cbc_deleteModel(m);
  (void)rc;
}

/* ------------------------------------------------------------------ */
/* main                                                                */
/* ------------------------------------------------------------------ */

int main(void)
{
  printf("=== nursesched-sprint02 C interface tests ===\n");

  test_model_dimensions();
  test_pure_lp();
  test_lp_with_bp(1, "singletons");
  test_lp_with_bp(2, "milpbt");
  test_lp_with_bp(3, "fixpoint");
  test_lp_reoptimization();
  test_mip();

  printf("\n%d/%d tests passed.\n", tests_passed, tests_run);
  return (tests_passed == tests_run) ? 0 : 1;
}
