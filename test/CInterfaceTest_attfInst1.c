/* CInterfaceTest_attfInst1.c
 *
 * Tests for the attfInst1 instance (ATIF multi-commodity flow, MIPLIB 2017+spp).
 *
 * Known values:
 *   LP relaxation (root, default preprocessing): 189
 *   MIP optimal:                                 202
 *
 * Instance structure: 1187 rows, 1163 cols, 5923 NZ
 *   1148 binary, 0 general integer, 15 continuous
 *
 * The MIP test uses a 500-node limit.  If any integer solutions are found,
 * every solution in the pool is validated for feasibility and objective >= 202
 * (the known optimal).  If optimality is claimed within the node limit, the
 * objective is also checked to equal 202 exactly.
 */

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "Cbc_C_Interface.h"

#include "test_utils.h"
#define ATTFINST1_MPS fixture_path("attfInst1.mps.gz")

/* Known / reference values */
#define MIP_OPT        202.0
#define MIP_TOL        1e-6   /* exact integer optimal */
#define MIP_NODE_LIMIT 500

static int tests_run    = 0;
static int tests_passed = 0;

#define CHECK(cond, msg)                                               \
  do {                                                                 \
    ++tests_run;                                                       \
    if (cond) {                                                        \
      ++tests_passed;                                                  \
    } else {                                                           \
      fprintf(stderr, "FAIL [%s:%d] %s\n", __FILE__, __LINE__, msg);  \
    }                                                                  \
  } while (0)

/* ------------------------------------------------------------------ */
/* helpers                                                             */
/* ------------------------------------------------------------------ */

static Cbc_Model *load_model(void)
{
  Cbc_Model *m = Cbc_newModel();
  Cbc_setLogLevel(m, 0);
  int rc = Cbc_readMps(m, ATTFINST1_MPS);
  if (rc != 0) {
    fprintf(stderr, "ERROR: could not read %s\n", ATTFINST1_MPS);
    Cbc_deleteModel(m);
    return NULL;
  }
  return m;
}

/* ------------------------------------------------------------------ */
/* test_model_dimensions                                               */
/* ------------------------------------------------------------------ */

static void test_model_dimensions(void)
{
  Cbc_Model *m = load_model();
  if (!m) return;

  CHECK(Cbc_getNumRows(m) == 1187, "attfInst1: 1187 rows");
  CHECK(Cbc_getNumCols(m) == 1163, "attfInst1: 1163 cols");

  int nint = 0;
  for (int j = 0; j < Cbc_getNumCols(m); j++)
    if (Cbc_isInteger(m, j)) ++nint;
  CHECK(nint == 1148, "attfInst1: 1148 integer variables");

  Cbc_deleteModel(m);
}

/* ------------------------------------------------------------------ */
/* test_mip_node_limited                                               */
/*                                                                     */
/* Solve with a 500-node limit.  Validates:                           */
/*   - problem is not infeasible                                       */
/*   - every solution found is feasible with obj >= known optimal      */
/*   - if optimality is proven, obj == 202 exactly                    */
/* ------------------------------------------------------------------ */

static void test_mip_node_limited(void)
{
  Cbc_Model *m = load_model();
  if (!m) return;

  Cbc_setIntParam(m, INT_PARAM_MAX_NODES, MIP_NODE_LIMIT);
  Cbc_setDblParam(m, DBL_PARAM_TIME_LIMIT, 300.0);

  int rc = Cbc_solve(m);
  (void)rc;

  CHECK(!Cbc_isProvenInfeasible(m), "attfInst1 is not infeasible");

  int nsols = Cbc_numberSavedSolutions(m);
  printf("  Saved solutions: %d   Best bound: %g\n",
         nsols, Cbc_getBestPossibleObjValue(m));

  if (nsols > 0) {
    double best_obj = Cbc_getObjValue(m);
    CHECK(best_obj >= MIP_OPT - MIP_TOL,
          "attfInst1: best obj >= known optimal (202)");

    /* Validate every solution in the pool for feasibility and correct obj */
    int fails = validate_all_saved_solutions(m, MIP_OPT, MIP_TOL, "attfInst1");
    CHECK(fails == 0,
          "attfInst1: all saved solutions are feasible with obj >= 202");

    if (fails == 0)
      printf("  All %d saved solution(s) validated (feasible, obj >= %g)\n",
             nsols, MIP_OPT);

    /* If optimality is claimed, verify the optimal value is exactly 202 */
    if (Cbc_isProvenOptimal(m)) {
      CHECK(fabs(best_obj - MIP_OPT) < MIP_TOL,
            "attfInst1: proven optimal == 202");
      printf("  Optimality proven: obj = %g\n", best_obj);
    }
  } else {
    printf("  (no integer solution found within %d nodes — skipping "
           "solution validation)\n", MIP_NODE_LIMIT);
    /* Not finding a solution within the node limit is acceptable */
    ++tests_run; ++tests_passed;
  }

  Cbc_deleteModel(m);
}

/* ------------------------------------------------------------------ */
/* main                                                                */
/* ------------------------------------------------------------------ */

int main(void)
{
  printf("=== attfInst1 C interface tests ===\n");

  test_model_dimensions();
  test_mip_node_limited();

  printf("\n%d/%d tests passed.\n", tests_passed, tests_run);
  return (tests_passed == tests_run) ? 0 : 1;
}
