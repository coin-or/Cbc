/* test_utils.h — shared helpers for MIPster C interface tests.
 *
 * Provides:
 *   fixture_path()               — locate test MPS fixture files
 *   validate_all_saved_solutions() — check every solution in CBC's pool
 *
 * Search order for fixture_path(fname):
 *   1. $MIPSTER_FIXTURE_DIR/<fname>           (runtime env override)
 *   2. FIXTURE_DIR/fixtures/<fname>           (compile-time build-tree path)
 *   3. ./fixtures/<fname>                     (cwd fallback)
 */

#ifndef MIPSTER_TEST_UTILS_H
#define MIPSTER_TEST_UTILS_H

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "Cbc_C_Interface.h"

#ifndef FIXTURE_DIR
#  define FIXTURE_DIR "."
#endif

static char _fp_buf[4096];

/* Returns a path to the fixture file <fname> (e.g. "nursesched-sprint02.mps.gz").
 * Uses a static buffer — consume before calling again. */
static const char *fixture_path(const char *fname)
{
  const char *env = getenv("MIPSTER_FIXTURE_DIR");
  if (env && env[0])
    snprintf(_fp_buf, sizeof(_fp_buf), "%s/%s", env, fname);
  else
    snprintf(_fp_buf, sizeof(_fp_buf), FIXTURE_DIR "/fixtures/%s", fname);
  return _fp_buf;
}

/* validate_all_saved_solutions — check every solution in CBC's pool.
 *
 * For each solution i in [0, Cbc_numberSavedSolutions(m)):
 *   1. Verifies integer feasibility via Cbc_checkFeasibility().
 *   2. If known_optimal is not NaN, verifies the solution objective is on the
 *      correct side of known_optimal (within obj_tol):
 *        - minimisation (ObjSense > 0): obj >= known_optimal - obj_tol
 *        - maximisation (ObjSense < 0): obj <= known_optimal + obj_tol
 *
 * tag: short label printed in diagnostic messages (e.g., the instance name).
 * Returns the number of failures (0 = all solutions passed). */
static int validate_all_saved_solutions(
    Cbc_Model *m,
    double known_optimal,
    double obj_tol,
    const char *tag)
{
  int n = Cbc_numberSavedSolutions(m);
  int failures = 0;
  double sense = Cbc_getObjSense(m); /* +1 min, -1 max */
  int check_obj = !isnan(known_optimal);

  for (int i = 0; i < n; i++) {
    const double *x = Cbc_savedSolution(m, i);
    double obj = Cbc_savedSolutionObj(m, i);

    double maxViolRow, maxViolCol;
    int rowIdx, colIdx;
    char feas = Cbc_checkFeasibility(m, x,
                                     &maxViolRow, &rowIdx,
                                     &maxViolCol, &colIdx);
    if (!feas) {
      fprintf(stderr,
              "[%s] FAIL saved sol %d infeasible: "
              "maxViolRow=%g (row %d)  maxViolCol=%g (col %d)\n",
              tag, i, maxViolRow, rowIdx, maxViolCol, colIdx);
      failures++;
    }

    if (check_obj) {
      if (sense > 0 && obj < known_optimal - obj_tol) {
        fprintf(stderr,
                "[%s] FAIL saved sol %d: obj %g < known_optimal %g (tol %g)\n",
                tag, i, obj, known_optimal, obj_tol);
        failures++;
      } else if (sense < 0 && obj > known_optimal + obj_tol) {
        fprintf(stderr,
                "[%s] FAIL saved sol %d: obj %g > known_optimal %g (tol %g)\n",
                tag, i, obj, known_optimal, obj_tol);
        failures++;
      }
    }
  }

  return failures;
}

#endif /* MIPSTER_TEST_UTILS_H */
