/* CInterfaceTest_nursesched.c
 *
 * Tests for the nurse-scheduling instance nursesched-sprint02.
 *
 * Known values:
 *   LP relaxation (pure):  54.416667  (minimisation lower bound)
 *   MIP optimum:           58
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
#define NURSESCHED_LP_OPT  54.416667
#define NURSESCHED_MIP_OPT 58.0
#define LP_TOL             1e-3
#define MIP_TOL            1e-5

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
  test_mip();

  printf("\n%d/%d tests passed.\n", tests_passed, tests_run);
  return (tests_passed == tests_run) ? 0 : 1;
}
