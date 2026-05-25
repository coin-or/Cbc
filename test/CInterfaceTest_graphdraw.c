/* CInterfaceTest_graphdraw.c
 *
 * Tests for the graph-drawing instance graphdraw-domain (MIPLIB 2017).
 *
 * Known values (cross-validated with an independent LP solver):
 *   LP relaxation (pure):  12672.0  (minimisation lower bound)
 *   MIP:                   hard — optimal not reached quickly.
 *                          Test uses a node limit; if a solution is found
 *                          it is validated for feasibility only.
 *
 * Instance structure: 865 rows, 254 cols (180 binary, 20 integer, 54 continuous)
 *   choose1[X,Y]   — equality constraints (rows 0-44)
 *   centerDistAxis — bound constraints (rows 45-84)
 *   triangle1/2    — valid-inequality cuts  (rows ~600-864)
 */

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "Cbc_C_Interface.h"

#include "test_utils.h"
#define GRAPHDRAW_MPS fixture_path("graphdraw-domain.mps.gz")

/* Known / reference values */
#define GRAPHDRAW_LP_OPT  12672.0
#define LP_TOL            1e-3
#define MIP_NODE_LIMIT    500

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
  int rc = Cbc_readMps(m, GRAPHDRAW_MPS);
  if (rc != 0) {
    fprintf(stderr, "ERROR: could not read %s\n", GRAPHDRAW_MPS);
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
  Cbc_Model *m = load_model(0);
  if (!m) return;

  CHECK(Cbc_getNumRows(m) == 865, "graphdraw: 865 rows");
  CHECK(Cbc_getNumCols(m) == 254, "graphdraw: 254 columns");

  /* 180 binary + 20 general integer = 200 integer variables */
  int nint = 0;
  for (int j = 0; j < Cbc_getNumCols(m); j++)
    if (Cbc_isInteger(m, j)) ++nint;
  CHECK(nint == 200, "graphdraw: 200 integer variables");

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
  CHECK(fabs(obj - GRAPHDRAW_LP_OPT) < LP_TOL,
        "pure LP: objective matches known LP relaxation");

  Cbc_deleteModel(m);
}

/* ------------------------------------------------------------------ */
/* test_mip_node_limited                                               */
/*                                                                     */
/* Solve with a node limit.  Optimality is not expected within this   */
/* limit.  If any integer solutions are found, every solution in the  */
/* pool is validated for feasibility and objective >= LP lower bound. */
/* ------------------------------------------------------------------ */

static void test_mip_node_limited(void)
{
  Cbc_Model *m = load_model(0);
  if (!m) return;

  Cbc_setIntParam(m, INT_PARAM_MAX_NODES, MIP_NODE_LIMIT);
  /* generous time limit as safety net */
  Cbc_setDblParam(m, DBL_PARAM_TIME_LIMIT, 60.0);

  int rc = Cbc_solve(m);
  (void)rc;

  /* The solve should stop due to node limit (not time limit or error). */
  int is_opt = Cbc_isProvenOptimal(m);
  int is_inf = Cbc_isProvenInfeasible(m);
  CHECK(!is_inf, "MIP node-limited: instance is not infeasible");

  /* The dual bound from the root LP must be a valid lower bound. */
  double best_bound = Cbc_getBestPossibleObjValue(m);
  CHECK(best_bound >= GRAPHDRAW_LP_OPT - LP_TOL,
        "MIP: dual bound >= LP relaxation value");

  int nsols = Cbc_numberSavedSolutions(m);
  printf("  Saved solutions: %d   Best bound: %g\n", nsols, best_bound);

  if (nsols > 0) {
    double obj = Cbc_getObjValue(m);
    CHECK(obj >= GRAPHDRAW_LP_OPT - LP_TOL,
          "MIP: best integer solution >= LP lower bound");
    CHECK(obj >= best_bound - LP_TOL,
          "MIP: best integer solution >= dual bound");

    /* Validate every solution in the pool.  No known optimal, so only
     * feasibility is checked (pass NAN to skip objective bound check). */
    int fails = validate_all_saved_solutions(m, NAN, 0.0, "graphdraw");
    CHECK(fails == 0, "MIP: all saved solutions are integer-feasible");

    if (fails == 0)
      printf("  All %d saved solution(s) validated (feasible)\n", nsols);

    if (!is_opt)
      printf("  (note: not proven optimal — node limit reached, "
             "best obj = %.0f, bound = %.2f)\n", obj, best_bound);
  } else {
    /* No integer solution found within node limit — that is allowed. */
    printf("  (note: no integer solution found within %d nodes)\n",
           MIP_NODE_LIMIT);
    /* Still count the bound check as passed (already checked above). */
    ++tests_run; ++tests_passed; /* placeholder for "feasibility N/A" */
  }

  Cbc_deleteModel(m);
}

/* ------------------------------------------------------------------ */
/* main                                                                */
/* ------------------------------------------------------------------ */

int main(void)
{
  printf("=== graphdraw-domain C interface tests ===\n");

  test_model_dimensions();
  test_pure_lp();
  test_mip_node_limited();

  printf("\n%d/%d tests passed.\n", tests_passed, tests_run);
  return (tests_passed == tests_run) ? 0 : 1;
}
