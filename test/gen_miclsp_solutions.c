/* gen_miclsp_solutions.c — generate reference optimal solution files for all
 * MICLSP test fixtures.
 *
 * For each fixture instance this program:
 *   1. Builds the MICLSP model via the C API (same builder as the test suite).
 *   2. Solves with preprocessing disabled and a caller-supplied time limit so
 *      the solution variable names match the original column names exactly
 *      (required by -debugCuts / CbcMipStart::read).
 *   3. If the solve is proven optimal AND the objective matches the certified
 *      value, writes a .sol file in the format expected by -debugCuts:
 *          <col_idx> <col_name> <value>
 *      one line per non-zero variable.
 *
 * Usage:
 *   gen_miclsp_solutions [outdir [timelimit_seconds]]
 *
 * Defaults:
 *   outdir        = fixtures/solutions   (relative to cwd; created if absent)
 *   timelimit     = 1800 seconds
 *
 * On success the exit code is 0.  A non-zero exit indicates that at least one
 * instance that should be solvable was not.
 */

#include "miclsp_builder.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>


static int write_solution(Cbc_Model *m, const char *outpath)
{
  FILE *fp = fopen(outpath, "w");
  if (!fp) {
    fprintf(stderr, "  ERROR: cannot open %s for writing\n", outpath);
    return 0;
  }
  int ncols = Cbc_getNumCols(m);
  const double *sol = Cbc_getColSolution(m);
  char name[256];
  for (int i = 0; i < ncols; i++) {
    if (fabs(sol[i]) > 1e-9) {
      Cbc_getColName(m, i, name, sizeof(name));
      fprintf(fp, "%d %s %.15g\n", i, name, sol[i]);
    }
  }
  fclose(fp);
  return 1;
}

int main(int argc, char *argv[])
{
  const char *outdir = "fixtures/solutions";
  double timelimit = 1800.0;

  if (argc >= 2)
    outdir = argv[1];
  if (argc >= 3)
    timelimit = atof(argv[2]);

#ifdef _WIN32
  _mkdir(outdir);
#else
  mkdir(outdir, 0755);
#endif

  printf("gen_miclsp_solutions: outdir=%s  timelimit=%.0fs\n\n", outdir, timelimit);

  int n_optimal = 0, n_timelimit = 0, n_wrong = 0;

  for (int k = 0; k < N_MICLSP_INSTANCES; k++) {
    const MiclspInstance *inst = &MICLSP_INSTANCES[k];
    const char *name = MICLSP_NAMES[k];

    printf("[%2d/%d] %-40s  opt=%d  ", k + 1, N_MICLSP_INSTANCES, name, inst->opt);
    fflush(stdout);

    Cbc_Model *m = build_miclsp(inst);
    Cbc_setMaximumSeconds(m, timelimit);
    Cbc_setParameter(m, "preprocess", "off");

    Cbc_solve(m);

    double obj = Cbc_getObjValue(m);
    int proven = Cbc_isProvenOptimal(m);

    if (proven && fabs(obj - (double)inst->opt) < 0.5) {
      char outpath[512];
      snprintf(outpath, sizeof(outpath), "%s/%s.sol", outdir, name);
      if (write_solution(m, outpath)) {
        printf("OPTIMAL  obj=%.0f  -> %s\n", obj, outpath);
        n_optimal++;
      } else {
        n_wrong++;
      }
    } else if (!proven && fabs(obj - (double)inst->opt) < 0.5) {
      /* Found the certified optimal value but didn't prove it (time/node limit).
       * Save the solution anyway — it is a valid integer feasible point at the
       * certified optimal value, suitable for use as a debugCuts reference. */
      char outpath[512];
      snprintf(outpath, sizeof(outpath), "%s/%s.sol", outdir, name);
      if (write_solution(m, outpath)) {
        printf("TIME LIMIT (obj=certified)  obj=%.0f  -> %s\n", obj, outpath);
        n_optimal++;
      } else {
        n_wrong++;
      }
    } else if (!proven) {
      printf("TIME/NODE LIMIT  obj=%.0f  (certified=%d) — no solution saved\n", obj, inst->opt);
      n_timelimit++;
    } else {
      printf("WRONG OPTIMAL  obj=%.0f  (certified=%d) — skipping\n", obj, inst->opt);
      n_wrong++;
    }

    Cbc_deleteModel(m);
  }

  printf("\n=== Summary: %d optimal  %d time/node limit  %d wrong ===\n",
    n_optimal, n_timelimit, n_wrong);

  return (n_wrong > 0) ? 1 : 0;
}
