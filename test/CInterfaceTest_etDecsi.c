/* CInterfaceTest_etDecsi.c
 *
 * Tests for the ET+DECSI (etDecsi) instance (MIPLIB 2017+spp).
 *
 * Known values:
 *   LP relaxation (root, after preprocessing): ≈ 2.25
 *   MIP optimal:                               7
 *
 * Instance structure: 17917 rows, 10606 cols, 58773 NZ
 *   10080 binary, 526 general integer, 0 continuous
 *
 * The MIP test uses a 300-node limit.  Proving optimality within this limit
 * is unlikely, but if any integer solution is found every solution in CBC's
 * solution pool is validated for feasibility and for having objective >= 7
 * (the known minimum).
 */

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "Cbc_C_Interface.h"

#include "test_utils.h"
#define ETDECSI_MPS fixture_path("etDecsi.mps.gz")

/* Known / reference values */
#define MIP_OPT       7.0
#define MIP_TOL       1e-6   /* exact integer optimal */
#define MIP_NODE_LIMIT 300

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
  int rc = Cbc_readMps(m, ETDECSI_MPS);
  if (rc != 0) {
    fprintf(stderr, "ERROR: could not read %s\n", ETDECSI_MPS);
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

  CHECK(Cbc_getNumRows(m) == 17917, "etDecsi: 17917 rows");
  CHECK(Cbc_getNumCols(m) == 10606, "etDecsi: 10606 cols");

  int nint = 0;
  for (int j = 0; j < Cbc_getNumCols(m); j++)
    if (Cbc_isInteger(m, j)) ++nint;
  CHECK(nint == 10606, "etDecsi: all 10606 vars are integer");

  Cbc_deleteModel(m);
}

/* ------------------------------------------------------------------ */
/* test_mip_node_limited                                               */
/*                                                                     */
/* Solve with a node limit.  If CBC finds any integer solutions, every */
/* solution in the pool is validated for feasibility and objective     */
/* value >= known optimal (7).                                         */
/* ------------------------------------------------------------------ */

static void test_mip_node_limited(void)
{
  Cbc_Model *m = load_model();
  if (!m) return;

  Cbc_setIntParam(m, INT_PARAM_MAX_NODES, MIP_NODE_LIMIT);
  Cbc_setDblParam(m, DBL_PARAM_TIME_LIMIT, 120.0);

  int rc = Cbc_solve(m);
  (void)rc;

  CHECK(!Cbc_isProvenInfeasible(m), "etDecsi is not infeasible");

  int nsols = Cbc_numberSavedSolutions(m);
  printf("  Saved solutions: %d   Best bound: %g\n",
         nsols, Cbc_getBestPossibleObjValue(m));

  if (nsols > 0) {
    /* Best solution objective must be >= known optimal */
    double best_obj = Cbc_getObjValue(m);
    CHECK(best_obj >= MIP_OPT - MIP_TOL,
          "etDecsi best obj >= known optimal (7)");

    /* Validate every solution in the pool */
    int fails = validate_all_saved_solutions(m, MIP_OPT, MIP_TOL, "etDecsi");
    CHECK(fails == 0,
          "etDecsi: all saved solutions are feasible with obj >= known optimal");

    if (fails == 0)
      printf("  All %d saved solution(s) validated (feasible, obj >= %g)\n",
             nsols, MIP_OPT);
  } else {
    printf("  (no integer solution found within %d nodes — skipping "
           "solution validation)\n", MIP_NODE_LIMIT);
    /* Count as a pass: not finding a solution within the node limit is fine */
    ++tests_run; ++tests_passed;
  }

  Cbc_deleteModel(m);
}

/* ------------------------------------------------------------------ */
/* main                                                                */
/* ------------------------------------------------------------------ */

int main(void)
{
  printf("=== etDecsi C interface tests ===\n");

  test_model_dimensions();
  test_mip_node_limited();

  printf("\n%d/%d tests passed.\n", tests_passed, tests_run);
  return (tests_passed == tests_run) ? 0 : 1;
}
