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
#include <time.h>

#include "Cbc_C_Interface.h"

#ifndef FIXTURE_DIR
#  define FIXTURE_DIR "."
#endif

static char _fp_buf[4096];
static char _sp_buf[4096];

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

/* Returns the path to a reference solution file for <base> (e.g. "j3041_1"),
 * or NULL if the file does not exist.  Solutions live in fixtures/solutions/.
 * Uses $MIPSTER_SOLUTIONS_DIR if set, otherwise FIXTURE_DIR/fixtures/solutions/.
 * Uses a static buffer — consume before calling again. */
static const char *solution_path(const char *base)
{
  const char *env = getenv("MIPSTER_SOLUTIONS_DIR");
  if (env && env[0])
    snprintf(_sp_buf, sizeof(_sp_buf), "%s/%s.sol", env, base);
  else
    snprintf(_sp_buf, sizeof(_sp_buf), FIXTURE_DIR "/fixtures/solutions/%s.sol", base);

  FILE *f = fopen(_sp_buf, "r");
  if (!f)
    return NULL;
  fclose(f);
  return _sp_buf;
}

/* validate_all_saved_solutions — check every solution in CBC's pool.
 *
 * For each solution i in [0, Cbc_numberSavedSolutions(m)):
 *   1. Feasibility: Cbc_checkFeasibility() — bounds, integrality, rows.
 *   2. Objective consistency: computed c·x must match the stored
 *      Cbc_savedSolutionObj() within a relative tolerance of 1e-6.
 *      Catches internal inconsistencies (solver stores a wrong obj value).
 *   3. Bound check (if known_optimal is not NaN): stored obj must not be
 *      better than the certified optimum (within obj_tol):
 *        - minimisation (ObjSense > 0): obj >= known_optimal - obj_tol
 *        - maximisation (ObjSense < 0): obj <= known_optimal + obj_tol
 *
 * tag: short label for diagnostic messages (e.g., the instance name).
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
  int check_bound = !isnan(known_optimal);
  int ncols = Cbc_getNumCols(m);
  const double *objCoef = Cbc_getObjCoefficients(m);

  for (int i = 0; i < n; i++) {
    const double *x = Cbc_savedSolution(m, i);
    double reported_obj = Cbc_savedSolutionObj(m, i);

    /* 1. Feasibility */
    double maxViolRow, maxViolCol;
    int rowIdx, colIdx;
    if (!Cbc_checkFeasibility(m, x, &maxViolRow, &rowIdx,
                                    &maxViolCol, &colIdx)) {
      fprintf(stderr,
              "[%s] FAIL saved sol %d infeasible: "
              "maxViolRow=%g (row %d)  maxViolCol=%g (col %d)\n",
              tag, i, maxViolRow, rowIdx, maxViolCol, colIdx);
      failures++;
    }

    /* 2. Objective consistency: computed c·x vs reported obj */
    double computed_obj = 0.0;
    for (int j = 0; j < ncols; j++)
      computed_obj += objCoef[j] * x[j];
    double obj_rel_tol = 1e-6 * fmax(1.0, fabs(reported_obj));
    if (fabs(computed_obj - reported_obj) > obj_rel_tol) {
      fprintf(stderr,
              "[%s] FAIL saved sol %d: computed obj %g != reported obj %g"
              " (diff=%g)\n",
              tag, i, computed_obj, reported_obj,
              computed_obj - reported_obj);
      failures++;
    }

    /* 3. Bound check: pool solution must not beat the certified optimum */
    if (check_bound) {
      if (sense > 0 && reported_obj < known_optimal - obj_tol) {
        fprintf(stderr,
                "[%s] FAIL saved sol %d: obj %g < known_optimal %g (tol %g)\n",
                tag, i, reported_obj, known_optimal, obj_tol);
        failures++;
      } else if (sense < 0 && reported_obj > known_optimal + obj_tol) {
        fprintf(stderr,
                "[%s] FAIL saved sol %d: obj %g > known_optimal %g (tol %g)\n",
                tag, i, reported_obj, known_optimal, obj_tol);
        failures++;
      }
    }
  }

  return failures;
}

/* ------------------------------------------------------------------ */
/* Performance measurement                                              */
/* ------------------------------------------------------------------ */

typedef struct {
  const char *name;
  int         pass;
  int         is_optimal;
  long        nodes;
  double      wall_time;  /* solve wall-clock seconds */
  double      obj;        /* displayed objective (after any sign flip); 1e30 = no solution */
} MipPerfRecord;

static double perf_wall_time(void)
{
  struct timespec ts;
  clock_gettime(CLOCK_MONOTONIC, &ts);
  return (double)ts.tv_sec + ts.tv_nsec * 1e-9;
}

/* Print a per-instance table + aggregate summary.
 * title: short label, e.g. "Set Covering". */
static void print_perf_summary(const MipPerfRecord *recs, int n, const char *title)
{
  printf("\n=== %s — Performance Summary ===\n", title);
  printf("  %-42s %7s  %8s  %12s  %s\n",
         "Instance", "Nodes", "Time(s)", "Obj", "Status");
  printf("  %-42s %7s  %8s  %12s  %s\n",
         "------------------------------------------", "-------", "--------",
         "------------", "----------");

  long   total_nodes = 0;
  double total_time  = 0.0;
  int    n_passed = 0, n_failed = 0, n_optimal = 0;

  for (int i = 0; i < n; i++) {
    const MipPerfRecord *r = &recs[i];
    const char *status = !r->pass    ? "FAIL"    :
                         r->is_optimal ? "optimal" : "feasible";
    if (r->obj < 1e29)
      printf("  %-42s %7ld  %8.3f  %12.2f  %s\n",
             r->name, r->nodes, r->wall_time, r->obj, status);
    else
      printf("  %-42s %7ld  %8.3f  %12s  %s\n",
             r->name, r->nodes, r->wall_time, "-", status);
    total_nodes += r->nodes;
    total_time  += r->wall_time;
    if (r->pass) n_passed++; else n_failed++;
    if (r->is_optimal) n_optimal++;
  }

  printf("  %-42s %7s  %8s  %12s  %s\n",
         "------------------------------------------", "-------", "--------",
         "------------", "----------");
  printf("  Total: %d instances  %d passed  %d failed  %d proven optimal\n",
         n, n_passed, n_failed, n_optimal);
  if (total_nodes > 0 && total_time > 1e-9)
    printf("  Nodes: %ld  Time: %.3fs  Throughput: %.0f nodes/s  [machine index]\n",
           total_nodes, total_time, (double)total_nodes / total_time);
  else
    printf("  Nodes: %ld  Time: %.3fs  [all solved at root]\n",
           total_nodes, total_time);
}

/* Compute the effective MIP objective tolerance after a solve.
 * Accounts for the solver's ratioGap and allowableGap settings, which
 * determine how far from the true optimal the solver may stop and still
 * declare isProvenOptimal() == 1.  The returned tolerance is the absolute
 * deviation that should be accepted as "correct" for certified_opt.
 *
 * base_tol: minimum absolute tolerance for floating-point noise (e.g. 1e-6).
 *
 * Formula: max(base_tol, allowableGap, ratioGap * |certified_opt|)
 *
 * Example: ratioGap=0.0001 (MIPster default), certified_opt=3200022
 *   → effective tol ≈ 320  (any integer within 320 of optimal is accepted)
 */
static double mip_obj_tol(Cbc_Model *m, double certified_opt, double base_tol)
{
  double ratio_gap = Cbc_getDblParam(m, DBL_PARAM_GAP_RATIO);
  double abs_gap   = Cbc_getDblParam(m, DBL_PARAM_ALLOWABLE_GAP);
  double gap_tol   = fmax(abs_gap, ratio_gap * fabs(certified_opt));
  return fmax(base_tol, gap_tol);
}

/* Builder for mip_diag: loads a model from an MPS file.
 * userdata must be a (const char *) path to the .mps or .mps.gz file. */
static Cbc_Model *build_mps_model(void *userdata)
{
  const char *path = (const char *)userdata;
  Cbc_Model *m = Cbc_newModel();
  if (Cbc_readMps(m, path) != 0) {
    Cbc_deleteModel(m);
    return NULL;
  }
  return m;
}

#endif /* MIPSTER_TEST_UTILS_H */
