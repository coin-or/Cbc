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
 *     mip_diag_wrong_optimal(my_builder, data, certified_opt, 900, 10000);
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
 *   cgraph       : conflict graph (default "on")
 *   clqstr       : clique strengthening preprocessing (default "before")
 *   nodeBoundProp: knapsack bound propagation at B&B nodes (default "on")
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
    int time_limit_sec,
    int node_limit)
{
  printf("\n  [DIAG] Wrong optimal — running %d diagnostic configurations\n",
         MIP_N_DIAG_CONFIGS);
  printf("  [DIAG] certified_opt=%d  time_limit=%ds  node_limit=%d\n\n",
         certified_opt, time_limit_sec, node_limit);

  for (int k = 0; k < MIP_N_DIAG_CONFIGS; k++) {
    const MipDiagConfig *cfg = &MIP_DIAG_CONFIGS[k];
    printf("  [DIAG %2d/%d] %-22s  (%s=%s) ... ",
           k + 1, MIP_N_DIAG_CONFIGS, cfg->label, cfg->param, cfg->value);
    fflush(stdout);

    Cbc_Model *m = builder(userdata);
    Cbc_setLogLevel(m, 0);

    if (time_limit_sec > 0)
      Cbc_setDblParam(m, DBL_PARAM_TIME_LIMIT, (double)time_limit_sec);
    if (node_limit > 0)
      Cbc_setIntParam(m, INT_PARAM_MAX_NODES, node_limit);

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
        printf("WRONG opt=%.0f (expected %d) — BUG ALSO PRESENT\n",
               obj, certified_opt);
      else
        printf("OK  obj=%.0f matches certified=%d\n", obj, certified_opt);
    } else if (nsaved > 0) {
      double gap = (obj > 1e-10) ? 100.0 * (obj - bound) / obj : 0.0;
      printf("not proven  obj=%.0f  bound=%.0f  gap=%.1f%%  (%s)\n",
             obj, bound, gap,
             time_hit ? "time limit" : node_hit ? "node limit" : "stopped");
    } else {
      printf("no solution found (%s)\n",
             time_hit ? "time limit" : node_hit ? "node limit" : "stopped");
    }

    Cbc_deleteModel(m);
  }

  printf("\n  [DIAG] Configurations reporting OK narrow down the culprit.\n"
         "  [DIAG] If a run returns WRONG, that feature is not responsible.\n");
}

#endif /* MIP_DIAG_H */
