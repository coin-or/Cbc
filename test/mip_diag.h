/* mip_diag.h — shared diagnostic helper for MIP tests.
 *
 * When a MIP solve incorrectly proves optimality at the wrong value, call
 * mip_diag_wrong_optimal() to re-run with each default-enabled feature
 * disabled one at a time.  The output identifies which feature is responsible
 * for the incorrect pruning.
 *
 * Usage:
 *   static Cbc_Model *my_builder(void *ud) {
 *     MyData *d = (MyData *)ud;
 *     Cbc_Model *m = Cbc_newModel();
 *     // ... build model structure only (no time/node/cut/search params) ...
 *     return m;
 *   }
 *
 *   if (is_proven && fabs(obj - certified_opt) > 1.0) {
 *     mip_diag_wrong_optimal(my_builder, data, certified_opt, 900);
 *     // default debugCuts run (all features at default settings)
 *     mip_diag_debug_cuts(my_builder, data, certified_opt, 900, sol, NULL, NULL);
 *     // explicit nodeBoundProp=on run to isolate cut-generator soundness
 *     mip_diag_debug_cuts(my_builder, data, certified_opt, 900, sol,
 *                         "nodeBoundProp", "on");
 *   }
 *
 * Each diagnostic run applies both the time limit and the node limit;
 * whichever is reached first stops the run.
 */

#ifndef MIP_DIAG_H
#define MIP_DIAG_H

#include "Cbc_C_Interface.h"
#include <math.h>
#include <stdio.h>

typedef struct {
  const char *label;   /* human-readable run name           */
  const char *param;   /* Cbc_setParameter name             */
  const char *value;   /* value to set                      */
  const char *param2;  /* optional 2nd parameter (NULL if unused) */
  const char *value2;  /* optional 2nd value                */
  const char *param3;  /* optional 3rd parameter (NULL if unused) */
  const char *value3;  /* optional 3rd value                */
} MipDiagConfig;

/* Diagnostic configurations — only features that are ON by default.
 *
 * Infrastructure (default on):
 *   cgraph        : conflict graph (default "on")
 *   clqstr        : clique strengthening / clique merging (default "before")
 *   boundPropLevel: bound propagation before initial LP solve (default "milpbt")
 *   nodeBoundProp : knapsack bound propagation at B&B nodes (default "on")
 *   preprocess    : presolve / postsolve pipeline
 *
 * Cut families (default "ifmove"):
 *   gomory, knapsack, flow, mixed (MIR), probing, reduce2, twoMir, zeroHalf
 *
 * Tolerances (test whether tight FP tolerances change the outcome — relevant
 * for big-M formulations on platforms with FMA / different rounding):
 *   primalTolerance, integerTolerance, dualTolerance
 *
 * Also one combined run that turns off all standard cuts at once, and one
 * that turns off bound propagation at root AND nodes simultaneously.
 */
static const MipDiagConfig MIP_DIAG_CONFIGS[] = {
  /* infrastructure */
  { "no-cgraph",         "cgraph",          "off",   NULL,            NULL,  NULL,           NULL  },
  { "no-clqstr",         "clqstr",          "off",   NULL,            NULL,  NULL,           NULL  },
  { "no-boundproplevel", "boundPropL",      "off",   NULL,            NULL,  NULL,           NULL  },
  { "no-nodeboundprop",  "nodeBoundProp",   "off",   NULL,            NULL,  NULL,           NULL  },
  { "no-all-boundprop",  "boundPropL",      "off",   "nodeBoundProp", "off", NULL,           NULL  },
  { "no-preprocess",     "preprocess",      "off",   NULL,            NULL,  NULL,           NULL  },
  /* all cuts combined */
  { "no-all-cuts",       "cuts",            "off",   NULL,            NULL,  NULL,           NULL  },
  /* individual cuts */
  { "no-gomory",         "gomory",          "off",   NULL,            NULL,  NULL,           NULL  },
  { "no-knapsack",       "knapsack",        "off",   NULL,            NULL,  NULL,           NULL  },
  { "no-flow",           "flow",            "off",   NULL,            NULL,  NULL,           NULL  },
  { "no-mir",            "mixed",           "off",   NULL,            NULL,  NULL,           NULL  },
  { "no-probing",        "probing",         "off",   NULL,            NULL,  NULL,           NULL  },
  { "no-redsplit2",      "reduce2",         "off",   NULL,            NULL,  NULL,           NULL  },
  { "no-twomir",         "two",             "off",   NULL,            NULL,  NULL,           NULL  },
  { "no-zerohalf",       "zero",            "off",   NULL,            NULL,  NULL,           NULL  },
  /* tolerance tightening — relevant for big-M formulations (e.g. JSSP)
   * where FMA / FP rounding may interact with the default tolerances. */
  { "tight-prim-tol",    "primalTolerance", "1e-9",  NULL,            NULL,  NULL,           NULL  },
  { "tight-int-tol",     "integerTolerance","1e-9",  NULL,            NULL,  NULL,           NULL  },
  { "tight-all-tol",     "primalTolerance", "1e-9",  "integerTolerance","1e-9","dualTolerance","1e-9" },
  /* Integer-objective auto-detect probes:
   *
   * Cbc auto-detects integer objectives and sets cutoffIncrement ≈ 0.999
   * (CbcModel.cpp:1126). The fathoming check at CbcModel.cpp:11480 then
   * prunes any node whose LP-bound exceeds bestObj - 0.999 (i.e. anything
   * within ~0.001 of the incumbent). FP noise from FMA on Apple Silicon
   * could push the bound just above this threshold and wrongly prune the
   * optimum.
   *
   * tight-cutoff-inc: set increment=1.0 (slightly stricter than the auto-
   *   detected 0.999). Auto-detect won't override (since 0.999 < 1.0), so
   *   the user value sticks. Tests whether the 0.001 safety margin matters.
   * no-fake-obj: set more2MipOptions=lessused (bit 536870912) which, as a
   *   side-effect at CbcModel.cpp:503-505, disables createFake/randomCost
   *   for problems with few integer-objective variables (JSSP makespan
   *   formulation hits this code path with iType=7).
   */
  { "tight-cutoff-inc",  "increment",       "1.0",   NULL,            NULL,  NULL,           NULL  },
  { "no-fake-obj",       "more2!MipOptions","lessused", NULL,         NULL,  NULL,           NULL  },
};

#define MIP_N_DIAG_CONFIGS \
  ((int)(sizeof(MIP_DIAG_CONFIGS) / sizeof(MIP_DIAG_CONFIGS[0])))

/* Write the current best solution of m to path in mipstart-compatible format.
 * Returns 1 on success, 0 on failure. */
static int mip_diag_write_sol(Cbc_Model *m, const char *path)
{
  FILE *f = fopen(path, "w");
  if (!f) {
    fprintf(stderr, "[DIAG] cannot write solution to %s\n", path);
    return 0;
  }
  int n = Cbc_getNumCols(m);
  double obj = Cbc_getObjValue(m);
  const double *sol = Cbc_getColSolution(m);
  const double *objCoef = Cbc_getObjCoefficients(m);
  fprintf(f, "MIP solution - objective value %.10g\n", obj);
  char colname[256];
  for (int i = 0; i < n; i++) {
    Cbc_getColName(m, i, colname, sizeof(colname));
    fprintf(f, "%6d %-24s %24.15g %24.15g\n", i, colname, sol[i], objCoef[i]);
  }
  fclose(f);
  return 1;
}

/* Run all diagnostic configurations to locate the culprit feature.
 *
 * builder       : creates a fresh Cbc_Model (model structure only — no
 *                 time/node/cut/search parameters).  The diagnostic frees it.
 * userdata      : passed verbatim to builder.
 * certified_opt : known-correct MIP optimum.
 * time_limit_sec: per-run wall-clock limit in seconds (0 = no limit).
 */
static void mip_diag_wrong_optimal(
    Cbc_Model *(*builder)(void *userdata),
    void *userdata,
    double certified_opt,
    int time_limit_sec)
{
  printf("\n  [DIAG] Wrong optimal — running %d diagnostic configurations\n",
         MIP_N_DIAG_CONFIGS);
  printf("  [DIAG] certified_opt=%.6g  time_limit=%ds  (no node limit)\n\n",
         certified_opt, time_limit_sec);

  for (int k = 0; k < MIP_N_DIAG_CONFIGS; k++) {
    const MipDiagConfig *cfg = &MIP_DIAG_CONFIGS[k];
    if (cfg->param3)
      printf("  [DIAG %2d/%d] %-22s  (%s=%s, %s=%s, %s=%s)\n",
             k + 1, MIP_N_DIAG_CONFIGS, cfg->label,
             cfg->param,  cfg->value,
             cfg->param2, cfg->value2,
             cfg->param3, cfg->value3);
    else if (cfg->param2)
      printf("  [DIAG %2d/%d] %-22s  (%s=%s, %s=%s)\n",
             k + 1, MIP_N_DIAG_CONFIGS, cfg->label,
             cfg->param,  cfg->value,
             cfg->param2, cfg->value2);
    else
      printf("  [DIAG %2d/%d] %-22s  (%s=%s)\n",
             k + 1, MIP_N_DIAG_CONFIGS, cfg->label, cfg->param, cfg->value);
    fflush(stdout);

    Cbc_Model *m = builder(userdata);
    if (!m) {
      printf("  [DIAG %2d/%d] %-22s  → builder returned NULL, skipping\n",
             k + 1, MIP_N_DIAG_CONFIGS, cfg->label);
      continue;
    }
    Cbc_setLogLevel(m, 0);

    if (time_limit_sec > 0)
      Cbc_setDblParam(m, DBL_PARAM_TIME_LIMIT, (double)time_limit_sec);

    Cbc_setParameter(m, cfg->param,  cfg->value);
    if (cfg->param2) Cbc_setParameter(m, cfg->param2, cfg->value2);
    if (cfg->param3) Cbc_setParameter(m, cfg->param3, cfg->value3);
    Cbc_solve(m);

    int is_proven  = Cbc_isProvenOptimal(m);
    int time_hit   = Cbc_isSecondsLimitReached(m);
    int node_hit   = Cbc_isNodeLimitReached(m);
    int nsaved     = Cbc_numberSavedSolutions(m);
    double obj     = Cbc_getObjValue(m);
    double bound   = Cbc_getBestPossibleObjValue(m);

    if (is_proven) {
      /* Use the solver's own gap tolerances to judge correctness */
      double ratio_gap = Cbc_getDblParam(m, DBL_PARAM_GAP_RATIO);
      double abs_gap   = Cbc_getDblParam(m, DBL_PARAM_ALLOWABLE_GAP);
      double tol = fmax(abs_gap, ratio_gap * fabs(certified_opt));
      tol = fmax(tol, 1e-6);  /* always allow floating-point noise */
      if (fabs(obj - certified_opt) > tol)
        printf("  [DIAG %2d/%d] %-22s  → WRONG opt=%.10g (expected %.10g, tol=%.4g)"
               " — BUG ALSO PRESENT\n",
               k + 1, MIP_N_DIAG_CONFIGS, cfg->label, obj, certified_opt, tol);
      else
        printf("  [DIAG %2d/%d] %-22s  → OK  obj=%.10g certified=%.10g (tol=%.4g)\n",
               k + 1, MIP_N_DIAG_CONFIGS, cfg->label, obj, certified_opt, tol);
    } else if (nsaved > 0) {
      double gap = (obj > 1e-10) ? 100.0 * (obj - bound) / obj : 0.0;
      double rel_err = fabs(obj - certified_opt) / fmax(fabs(certified_opt), 1.0);
      int correct = (rel_err <= 1e-3);
      printf("  [DIAG %2d/%d] %-22s  → not proven  obj=%.10g  bound=%.10g"
             "  gap=%.4f%%  (%s)%s\n",
             k + 1, MIP_N_DIAG_CONFIGS, cfg->label, obj, bound, gap,
             time_hit ? "time limit" : node_hit ? "node limit" : "stopped",
             correct ? "  *** LEAD: found correct obj — this feature is SUSPECT ***" : "");
    } else {
      printf("  [DIAG %2d/%d] %-22s  → no solution found (%s)\n",
             k + 1, MIP_N_DIAG_CONFIGS, cfg->label,
             time_hit ? "time limit" : node_hit ? "node limit" : "stopped");
    }
    fflush(stdout);

    Cbc_deleteModel(m);
  }

  printf("\n  [DIAG] OK = feature not responsible for the bug.\n"
         "  [DIAG] WRONG = bug present without this feature too.\n"
         "  [DIAG] LEAD = found correct obj (not proven): this feature is suspect.\n");
}

/* Run a solve with OsiRowCutDebugger active to identify exactly which cut
 * incorrectly excludes the certified optimal solution.
 *
 * ref_sol_path  : path to a .sol file containing the certified optimal
 *                 solution (in mipster -writeSolution format or in the
 *                 "idx colname value" format used by gen_miclsp_solutions).
 *                 Every cut generated during the solve will be checked
 *                 against this solution; any cut that excludes it prints a
 *                 diagnostic message to stdout.
 *
 * extra_param / extra_value : optional extra parameter to set before solving
 *                 (e.g. "nodeBoundProp", "on").  Pass NULL for default config.
 *
 * If no "bad row" lines appear in the output the row-cut generators are NOT
 * the source of the bug; the problem likely lies in bound propagation
 * (column-bound changes / node pruning) rather than in row cuts.
 */
static void mip_diag_debug_cuts(
    Cbc_Model *(*builder)(void *userdata),
    void *userdata,
    double certified_opt,
    int time_limit_sec,
    const char *ref_sol_path,
    const char *extra_param,
    const char *extra_value)
{
  if (!ref_sol_path) {
    printf("\n  [DIAG] debugCuts: no reference solution path provided,"
           " skipping.\n");
    return;
  }

  if (extra_param)
    printf("\n  [DIAG] debugCuts pass (%s=%s) — reference solution: %s\n",
           extra_param, extra_value, ref_sol_path);
  else
    printf("\n  [DIAG] debugCuts pass (default config) — reference solution:"
           " %s\n", ref_sol_path);
  printf("  [DIAG] time_limit=%ds  no node limit"
         " — any invalid cut will be flagged below.\n", time_limit_sec);
  printf("  [DIAG] If no 'bad row' lines appear: row-cut generators are NOT"
         " the source of the bug.\n\n");
  fflush(stdout);

  Cbc_Model *m = builder(userdata);
  if (!m) {
    printf("  [DIAG] debugCuts: builder returned NULL, skipping.\n");
    return;
  }
  /* Keep log level 1 so OsiRowCutDebugger messages reach stdout. */
  Cbc_setLogLevel(m, 1);

  if (time_limit_sec > 0)
    Cbc_setDblParam(m, DBL_PARAM_TIME_LIMIT, (double)time_limit_sec);

  if (extra_param)
    Cbc_setParameter(m, extra_param, extra_value);

  /* Activate the row-cut debugger for the given solution file. */
  Cbc_setParameter(m, "debugCuts", ref_sol_path);

  Cbc_solve(m);

  int is_proven = Cbc_isProvenOptimal(m);
  double obj    = Cbc_getObjValue(m);

  printf("\n  [DIAG] debugCuts result: %s  obj=%.10g  (certified=%.10g)\n",
         is_proven ? "proven" : "not proven", obj, certified_opt);
  fflush(stdout);

  Cbc_deleteModel(m);
}

#endif /* MIP_DIAG_H */
