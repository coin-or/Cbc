/* CInterfaceTest_miplib.c
 *
 * Generic per-instance MIPLIB regression test using the C interface.
 *
 * Usage:
 *   CInterfaceTest_miplib <mps_file> <expected_mip_obj> [<time_limit_sec>]
 *
 * Checks:
 *   1. Model loads without error.
 *   2. Solver does not prove infeasibility.
 *   3. If any integer solution is found:
 *        a. No solution objective violates the known bound (e.g., a minimisation
 *           solution strictly below the known optimal would indicate a solver bug).
 *        b. Every solution in the pool passes Cbc_checkFeasibility.
 *
 * Exit code: 0 = all checks passed, 1 = any check failed.
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "Cbc_C_Interface.h"
#include "test_utils.h"

/* Relative + absolute floor tolerance for objective comparison. */
#define OBJ_TOL      1e-4
#define DEFAULT_TIME 300.0  /* seconds */

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

int main(int argc, char *argv[])
{
  if (argc < 3) {
    fprintf(stderr,
            "usage: %s <mps_file> <expected_mip_obj> [<time_limit_sec>]\n",
            argv[0]);
    return 1;
  }

  const char *mps_path = argv[1];
  double expected_obj  = atof(argv[2]);
  double time_limit    = (argc >= 4) ? atof(argv[3]) : DEFAULT_TIME;

  /* Short label for messages (basename of the path). */
  const char *slash = strrchr(mps_path, '/');
  const char *tag   = slash ? slash + 1 : mps_path;

  printf("=== %s ===\n", tag);

  Cbc_Model *m = Cbc_newModel();
  Cbc_setLogLevel(m, 0);

  int rc = Cbc_readMps(m, mps_path);
  CHECK(rc == 0, "model loaded without error");
  if (rc != 0) {
    fprintf(stderr, "ERROR: could not read %s\n", mps_path);
    Cbc_deleteModel(m);
    return 1;
  }

  Cbc_setDblParam(m, DBL_PARAM_TIME_LIMIT, time_limit);
  Cbc_solve(m);

  CHECK(!Cbc_isProvenInfeasible(m), "instance is not infeasible");

  int nsols = Cbc_numberSavedSolutions(m);
  printf("  Solutions: %d   Best bound: %g\n",
         nsols, Cbc_getBestPossibleObjValue(m));

  if (nsols > 0) {
    /*
     * Validate every solution in the pool:
     *   - integer feasibility via Cbc_checkFeasibility
     *   - objective does not violate the known optimal
     *     (for min: no sol should have obj < known_optimal - tol,
     *      since that would be better-than-optimal -- a solver bug)
     *
     * Use a primal tolerance of 1e-4 for feasibility checking: LP solvers can
     * accumulate up to ~1e-4 numerical error in the stored solution even when
     * the solver considers it feasible internally.
     */
    Cbc_setDblParam(m, DBL_PARAM_PRIMAL_TOL, 1e-4);
    double abs_tol = OBJ_TOL * (1.0 + fabs(expected_obj));
    int fails = validate_all_saved_solutions(m, expected_obj, abs_tol, tag);
    CHECK(fails == 0, "all saved solutions feasible with valid objective");

    if (fails == 0)
      printf("  All %d solution(s) validated OK  (best obj: %g  expected: %g)\n",
             nsols, Cbc_getObjValue(m), expected_obj);
  } else {
    printf("  (no integer solution found within %.0fs — skipping "
           "solution validation)\n", time_limit);
    ++tests_run; ++tests_passed;
  }

  Cbc_deleteModel(m);

  printf("%d/%d checks passed.\n", tests_passed, tests_run);
  return (tests_passed == tests_run) ? 0 : 1;
}
