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
  const char *label;  /* human-readable run name */
  const char *param;  /* Cbc_setParameter name    */
  const char *value;  /* value to set             */
} MipDiagConfig;

/* Diagnostic configurations — only features that are ON by default.
 *
 * Infrastructure (default on):
 *   cgraph        : conflict graph (default "on")
 *   clqstr        : clique strengthening / clique merging (default "before")
 *   boundPropLevel: bound propagation before initial LP solve (default "milpbt")
 *   nodeBoundProp : knapsack bound propagation at B&B nodes (default "on")
 *
 * Cut families (default "ifmove"):
 *   gomory, knapsack, flow, mixed (MIR), probing, reduce2, twoMir, zeroHalf
 *
 * Also one combined run that turns off all standard cuts at once.
 */
static const MipDiagConfig MIP_DIAG_CONFIGS[] = {
  /* infrastructure */
  { "no-cgraph",        "cgraph",        "off" },
  { "no-clqstr",        "clqstr",        "off" },
  { "no-boundproplevel","boundPropL",    "off" },
  { "no-nodeboundprop", "nodeBoundProp", "off" },
  /* all cuts combined */
  { "no-all-cuts",      "cuts",          "off" },
  /* individual cuts */
  { "no-gomory",        "gomory",        "off" },
  { "no-knapsack",      "knapsack",      "off" },
  { "no-flow",          "flow",          "off" },
  { "no-mir",           "mixed",         "off" },
  { "no-probing",       "probing",       "off" },
  { "no-redsplit2",     "reduce2",       "off" },
  { "no-twomir",        "two",           "off" },
  { "no-zerohalf",      "zero",          "off" },
};

#define MIP_N_DIAG_CONFIGS \
  ((int)(sizeof(MIP_DIAG_CONFIGS) / sizeof(MIP_DIAG_CONFIGS[0])))

/* Run all diagnostic configurations to locate the culprit feature.
 *
 * builder       : creates a fresh Cbc_Model (model structure only — no
 *                 time/node/cut/search parameters).  The diagnostic frees it.
 * userdata      : passed verbatim to builder.
 * certified_opt : known-correct MIP optimum.
 * time_limit_sec: per-run wall-clock limit in seconds (0 = no limit).
 * node_limit    : per-run B&B node limit (0 = no limit).
 */
static void mip_diag_wrong_optimal(
    Cbc_Model *(*builder)(void *userdata),
    void *userdata,
    int certified_opt,
    int time_limit_sec)
{
  printf("\n  [DIAG] Wrong optimal — running %d diagnostic configurations\n",
         MIP_N_DIAG_CONFIGS);
  printf("  [DIAG] certified_opt=%d  time_limit=%ds  (no node limit)\n\n",
         certified_opt, time_limit_sec);

  for (int k = 0; k < MIP_N_DIAG_CONFIGS; k++) {
    const MipDiagConfig *cfg = &MIP_DIAG_CONFIGS[k];
    printf("  [DIAG %2d/%d] %-22s  (%s=%s)\n",
           k + 1, MIP_N_DIAG_CONFIGS, cfg->label, cfg->param, cfg->value);
    fflush(stdout);

    Cbc_Model *m = builder(userdata);
    Cbc_setLogLevel(m, 0);

    if (time_limit_sec > 0)
      Cbc_setDblParam(m, DBL_PARAM_TIME_LIMIT, (double)time_limit_sec);

    Cbc_setParameter(m, cfg->param, cfg->value);
    Cbc_solve(m);

    int is_proven  = Cbc_isProvenOptimal(m);
    int time_hit   = Cbc_isSecondsLimitReached(m);
    int node_hit   = Cbc_isNodeLimitReached(m);
    int nsaved     = Cbc_numberSavedSolutions(m);
    double obj     = Cbc_getObjValue(m);
    double bound   = Cbc_getBestPossibleObjValue(m);

    if (is_proven) {
      if (fabs(obj - (double)certified_opt) > 1.0)
        printf("  [DIAG %2d/%d] %-22s  → WRONG opt=%.0f (expected %d)"
               " — BUG ALSO PRESENT\n",
               k + 1, MIP_N_DIAG_CONFIGS, cfg->label, obj, certified_opt);
      else
        printf("  [DIAG %2d/%d] %-22s  → OK  obj=%.0f certified=%d\n",
               k + 1, MIP_N_DIAG_CONFIGS, cfg->label, obj, certified_opt);
    } else if (nsaved > 0) {
      double gap = (obj > 1e-10) ? 100.0 * (obj - bound) / obj : 0.0;
      int correct = (fabs(obj - (double)certified_opt) <= 1.0);
      printf("  [DIAG %2d/%d] %-22s  → not proven  obj=%.0f  bound=%.0f"
             "  gap=%.1f%%  (%s)%s\n",
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
    int certified_opt,
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

  printf("\n  [DIAG] debugCuts result: %s  obj=%.0f  (certified=%d)\n",
         is_proven ? "proven" : "not proven", obj, certified_opt);
  fflush(stdout);

  Cbc_deleteModel(m);
}

#endif /* MIP_DIAG_H */
