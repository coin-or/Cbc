/* CInterfaceTest_trdta5581.c
 *
 * Tests for the trdta5581 instance (integer infeasible).
 *
 * Known properties:
 *   LP relaxation:   infeasible (LP itself is infeasible)
 *   MIP:             infeasible — bound propagation detects this in ~0.007 s
 *
 * Instance structure: 147771 rows, 107140 cols (104715 binary, 2425 integer)
 * Bound propagation proves infeasibility at row classTeacherPr(320,0,1184939)
 * in round 2 without needing an LP solve.
 */

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

#include "Cbc_C_Interface.h"

#ifndef FIXTURE_DIR
#  define FIXTURE_DIR "."
#endif
#define TRDTA_MPS FIXTURE_DIR "/fixtures/trdta5581.mps.gz"

static int tests_run    = 0;
static int tests_passed = 0;

#define CHECK(cond, msg)                                              \
  do {                                                               \
    ++tests_run;                                                     \
    if (cond) {                                                      \
      ++tests_passed;                                                \
    } else {                                                         \
      fprintf(stderr, "FAIL [%s:%d]: %s\n", __FILE__, __LINE__, msg); \
    }                                                                \
  } while (0)

static Cbc_Model *load_model(void)
{
  Cbc_Model *m = Cbc_newModel();
  Cbc_setLogLevel(m, 0);
  int rc = Cbc_readMps(m, TRDTA_MPS);
  assert(rc == 0);
  return m;
}

/* ------------------------------------------------------------------ */
static void test_model_dimensions(void)
{
  printf("test_model_dimensions\n");
  Cbc_Model *m = load_model();

  CHECK(Cbc_getNumRows(m) == 147771,   "expected 147771 rows");
  CHECK(Cbc_getNumCols(m) == 107140,   "expected 107140 cols");

  int nbin = 0, nint = 0;
  for (int j = 0; j < Cbc_getNumCols(m); ++j) {
    if (Cbc_isInteger(m, j)) {
      const double *lb = Cbc_getColLower(m);
      const double *ub = Cbc_getColUpper(m);
      if (lb[j] == 0.0 && ub[j] == 1.0) nbin++; else nint++;
    }
  }
  CHECK(nbin == 104715, "expected 104715 binary vars");
  CHECK(nint ==   2425, "expected 2425 general integer vars");

  Cbc_deleteModel(m);
}

/* ------------------------------------------------------------------ */
static void test_bp_infeasibility(void)
{
  printf("test_bp_infeasibility\n");
  Cbc_Model *m = load_model();

  /* Use default BP settings (milpbt); just cap time and disable MIP overhead */
  Cbc_setIntParam(m, INT_PARAM_CGRAPH,          0);
  Cbc_setIntParam(m, INT_PARAM_CLIQUE_MERGING,  0);
  Cbc_setDblParam(m, DBL_PARAM_TIME_LIMIT,     30.0);

  Cbc_solve(m);

  CHECK(Cbc_isProvenInfeasible(m), "MIP must be proven infeasible via BP");
  /* If infeasible was proven there should be no solution */
  CHECK(!Cbc_isProvenOptimal(m),   "must not be proven optimal");

  Cbc_deleteModel(m);
}

/* ------------------------------------------------------------------ */
int main(void)
{
  printf("=== trdta5581 C-interface tests ===\n");

  test_model_dimensions();
  test_bp_infeasibility();

  printf("=== %d / %d tests passed ===\n", tests_passed, tests_run);
  return (tests_passed == tests_run) ? 0 : 1;
}
