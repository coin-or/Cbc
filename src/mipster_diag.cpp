/*
 * mipster_diag — standalone wrong-result diagnostic tool for MIPster
 *
 * Diagnoses wrong-optimal and wrong-infeasible MIP results by:
 *   1. Confirming the bug is reproducible within the given time limit.
 *   2. Auto-obtaining a reference solution (if -sol is not given) via a
 *      progressive preprocess-off/nodeBoundProp-off/cuts-off ladder.
 *   3. Running a -debugCuts pass which prints:
 *        "bad row"                  — for any row cut that excludes the
 *                                     certified optimal solution
 *        "nodeBoundProp BAD FIXING" — for any bound fixing in
 *                                     CbcBoundPropagation that contradicts
 *                                     the certified optimal solution
 *      If neither type of message appears, the bug lies outside row cuts and
 *      bound propagation; use -scan to run a full feature-disabling sweep.
 *   4. (Optional, -scan) Running 13 feature-disabling configurations to
 *      narrow down the culprit when debugCuts alone is not conclusive.
 *
 * Usage:
 *   mipster_diag <problem.{mps,lp}[.gz]> -opt <value> [OPTIONS]
 *
 *   -opt N        Certified optimal value (required)
 *   -sol FILE     Known-good .sol file for the -debugCuts pass;
 *                 if omitted the tool attempts to auto-obtain one
 *   -time N       Per-run wall-clock time limit in seconds (default: 60)
 *   -max          Maximisation problem (default: minimisation)
 *   -scan         Also run 13 feature-disabling configurations to locate
 *                 the culprit (slow: up to 13 × time-limit seconds extra)
 *   -p PAR VAL    Extra solver parameter applied to every run.
 *                 Repeat to add multiple parameters.
 *                 Use to reproduce bugs that only trigger with specific settings,
 *                 e.g.: -p nodeBoundProp on -p nodeBoundPropMinDepth 2
 *
 * To generate a reference solution manually:
 *   mipster problem.mps -preprocess off -solve -solu /tmp/ref.sol
 *   mipster_diag problem.mps -opt <N> -sol /tmp/ref.sol
 */

#include "Cbc_C_Interface.h"
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fcntl.h>
#include <unistd.h>
#include <vector>

/* ── Diagnostic feature configurations ──────────────────────── */

struct DiagConfig {
  const char *label;
  const char *param;
  const char *value;
};

static const DiagConfig CONFIGS[] = {
  /* infrastructure */
  { "no-cgraph",        "cgraph",        "off" },
  { "no-clqstr",        "clqstr",        "off" },
  { "no-boundproplevel","boundPropL",    "off" },
  { "no-nodeboundprop", "nodeBoundProp", "off" },
  /* all cuts combined */
  { "no-all-cuts",      "cuts",          "off" },
  /* individual cut families */
  { "no-gomory",        "gomory",        "off" },
  { "no-knapsack",      "knapsack",      "off" },
  { "no-flow",          "flow",          "off" },
  { "no-mir",           "mixed",         "off" },
  { "no-probing",       "probing",       "off" },
  { "no-redsplit2",     "reduce2",       "off" },
  { "no-twomir",        "two",           "off" },
  { "no-zerohalf",      "zero",          "off" },
};
static const int N_CONFIGS = (int)(sizeof(CONFIGS) / sizeof(CONFIGS[0]));

/* ── Suppress solver stdout noise during non-diagnostic solves ── */

/* Redirect both stdout and stderr to /dev/null, run the solve, then restore.
 * BAD FIXING and bad row use plain printf (stdout) — they would be silenced.
 * Only use this for baseline / scan / auto-sol runs; use Cbc_solve()
 * directly for the debugCuts pass where we want those messages. */
static void silent_solve(Cbc_Model *m)
{
  fflush(stdout);
  fflush(stderr);
  int saved_out = dup(STDOUT_FILENO);
  int saved_err = dup(STDERR_FILENO);
  int devnull = open("/dev/null", O_WRONLY);
  dup2(devnull, STDOUT_FILENO);
  dup2(devnull, STDERR_FILENO);
  close(devnull);

  Cbc_solve(m);

  fflush(stdout);
  fflush(stderr);
  dup2(saved_out, STDOUT_FILENO);
  dup2(saved_err, STDERR_FILENO);
  close(saved_out);
  close(saved_err);
}

/* ── Result classification ───────────────────────────────────── */

enum ResultKind { RK_OK, RK_LEAD, RK_WRONG_OPT, RK_INFEASIBLE,
                  RK_ABANDONED, RK_NOSOL };

static ResultKind classify(Cbc_Model *m, double cert)
{
  if (Cbc_isAbandoned(m))
    return RK_ABANDONED;
  if (Cbc_isProvenOptimal(m))
    return (fabs(Cbc_getObjValue(m) - cert) <= 0.5) ? RK_OK : RK_WRONG_OPT;
  if (Cbc_isProvenInfeasible(m))
    return RK_INFEASIBLE;
  if (Cbc_bestSolution(m))
    return (fabs(Cbc_getObjValue(m) - cert) <= 0.5) ? RK_LEAD : RK_NOSOL;
  return RK_NOSOL;
}

static const char *rk_name(ResultKind k)
{
  switch (k) {
  case RK_OK:         return "OK";
  case RK_LEAD:       return "LEAD";
  case RK_WRONG_OPT:  return "WRONG_OPT";
  case RK_INFEASIBLE: return "INFEASIBLE";
  case RK_ABANDONED:  return "ABANDONED";
  case RK_NOSOL:      return "NOSOL";
  }
  return "?";
}

/* ── Build a fresh model ─────────────────────────────────────── */

struct ExtraParam { const char *param; const char *value; };

struct Opts {
  const char *filename;
  int         is_max;
  int         time_limit;
  int         log_level;
  std::vector<ExtraParam> extra; /* applied to every run */
};

static Cbc_Model *make_model(const Opts &o)
{
  Cbc_Model *m = Cbc_newModel();
  Cbc_setLogLevel(m, o.log_level);

  char base[1024];
  snprintf(base, sizeof(base), "%s", o.filename);
  char *gz = strstr(base, ".gz");
  if (gz) *gz = '\0';
  const char *ext = strrchr(base, '.');
  int rc = (ext && strcmp(ext, ".lp") == 0)
           ? Cbc_readLp(m, o.filename)
           : Cbc_readMps(m, o.filename);
  if (rc) {
    fprintf(stderr, "mipster_diag: failed to read '%s'\n", o.filename);
    Cbc_deleteModel(m);
    return nullptr;
  }
  if (o.is_max)
    Cbc_setObjSense(m, -1.0);
  if (o.time_limit > 0)
    Cbc_setDblParam(m, DBL_PARAM_TIME_LIMIT, (double)o.time_limit);
  for (const auto &ep : o.extra)
    Cbc_setParameter(m, ep.param, ep.value);
  return m;
}

/* ── Write incumbent to a .sol file ──────────────────────────── */

static bool write_sol(Cbc_Model *m, const char *path)
{
  const double *sol = Cbc_bestSolution(m);
  if (!sol) return false;
  FILE *f = fopen(path, "w");
  if (!f) return false;
  int n = Cbc_getNumCols(m);
  char name[512];
  for (int j = 0; j < n; j++) {
    if (fabs(sol[j]) > 1e-10) {
      Cbc_getColName(m, j, name, sizeof(name));
      fprintf(f, "%d %s %g\n", j, name, sol[j]);
    }
  }
  fclose(f);
  return true;
}

/* ── Feature scan ────────────────────────────────────────────── */

/* Runs all CONFIGS. If any run produces a correct incumbent and
 * sol_out_path != NULL, writes the first such solution to sol_out_path.
 * Returns true if a solution was written. */
static bool run_feature_scan(const Opts &base, double cert,
                             const char *sol_out_path)
{
  printf("\n  Feature scan: %d configs, %ds time limit each\n\n",
         N_CONFIGS, base.time_limit);

  int n_ok = 0, n_lead = 0;
  bool wrote_sol = false;

  for (int k = 0; k < N_CONFIGS; k++) {
    const DiagConfig &cfg = CONFIGS[k];
    Cbc_Model *m = make_model(base);
    if (!m) {
      printf("  [%2d/%d] %-22s → LOAD FAILED\n", k + 1, N_CONFIGS, cfg.label);
      continue;
    }
    Cbc_setParameter(m, cfg.param, cfg.value);
    silent_solve(m);

    ResultKind rk = classify(m, cert);
    double obj   = Cbc_getObjValue(m);
    double bound = Cbc_getBestPossibleObjValue(m);
    int t_hit    = Cbc_isSecondsLimitReached(m);
    int n_hit    = Cbc_isNodeLimitReached(m);
    const char *stop = t_hit ? "time" : n_hit ? "nodes" : "stopped";

    switch (rk) {
    case RK_OK:
      printf("  [%2d/%d] %-22s → OK          obj=%.6g\n",
             k + 1, N_CONFIGS, cfg.label, obj);
      n_ok++;
      break;
    case RK_LEAD:
      printf("  [%2d/%d] %-22s → LEAD        obj=%.6g  bound=%.6g  (%s)"
             "  *** found certified opt — SUSPECT ***\n",
             k + 1, N_CONFIGS, cfg.label, obj, bound, stop);
      n_lead++;
      break;
    case RK_WRONG_OPT:
      printf("  [%2d/%d] %-22s → WRONG_OPT   obj=%.6g  certified=%.6g"
             "  — bug also present\n",
             k + 1, N_CONFIGS, cfg.label, obj, cert);
      break;
    case RK_INFEASIBLE:
      printf("  [%2d/%d] %-22s → INFEASIBLE               — bug also present\n",
             k + 1, N_CONFIGS, cfg.label);
      break;
    case RK_ABANDONED:
      printf("  [%2d/%d] %-22s → ABANDONED\n", k + 1, N_CONFIGS, cfg.label);
      break;
    case RK_NOSOL:
      printf("  [%2d/%d] %-22s → NOSOL       bound=%.6g  (%s)\n",
             k + 1, N_CONFIGS, cfg.label, bound, stop);
      break;
    }
    fflush(stdout);

    if (!wrote_sol && sol_out_path && (rk == RK_OK || rk == RK_LEAD))
      wrote_sol = write_sol(m, sol_out_path);

    Cbc_deleteModel(m);
  }

  printf("\n  Scan summary: %d OK (culprit identified), %d LEAD (suspect)\n",
         n_ok, n_lead);
  return wrote_sol;
}

/* ── Auto-obtain reference solution ─────────────────────────── */

/* Tries progressively more conservative configs to find a correct
 * incumbent and write it to sol_out_path. Returns true on success. */
static bool auto_obtain_sol(const Opts &base, double cert,
                            const char *sol_out_path)
{
  /* Each entry accumulates the previous params */
  static const char *ladder[][2] = {
    { "preprocess",    "off" },
    { "nodeBoundProp", "off" },
    { "cuts",          "off" },
  };
  static const int N_LADDER = 3;

  printf("\n  Auto-obtaining reference solution (trying %d progressively "
         "conservative configs)...\n", N_LADDER);
  fflush(stdout);

  for (int n = 1; n <= N_LADDER; n++) {
    Cbc_Model *m = make_model(base);
    if (!m) break;
    for (int i = 0; i < n; i++)
      Cbc_setParameter(m, ladder[i][0], ladder[i][1]);
    silent_solve(m);

    ResultKind rk = classify(m, cert);
    if (rk == RK_OK || rk == RK_LEAD) {
      bool ok = write_sol(m, sol_out_path);
      Cbc_deleteModel(m);
      if (ok) {
        printf("  Auto-sol: found correct incumbent (ladder level %d), "
               "written to %s\n", n, sol_out_path);
        return true;
      }
    }
    Cbc_deleteModel(m);
  }
  printf("  Auto-sol: no reference solution found.\n"
         "  Provide one manually:  mipster problem.mps -preprocess off "
         "-solve -solu /tmp/ref.sol\n");
  return false;
}

/* ── debugCuts pass ──────────────────────────────────────────── */

static void run_debug_cuts(const Opts &base, const char *sol_path, double cert)
{
  printf("\n  debugCuts pass — reference: %s\n"
         "  Watch for:\n"
         "    'bad row N lb <= act <= ub'           invalid row cut\n"
         "    'nodeBoundProp BAD FIXING (phase,rnd)' wrong bound fixing\n"
         "  No such lines → row cuts and nodeBoundProp are NOT the source.\n"
         "  Use -scan to run feature-disabling sweep.\n\n",
         sol_path);
  fflush(stdout);

  Cbc_Model *m = make_model(base);  /* log_level stays 0; BAD FIXING uses printf */
  if (!m) return;
  Cbc_setParameter(m, "debugCuts", sol_path);

  /* Suppress stderr noise (▶ Non-default parameters) while keeping stdout
   * visible so BAD FIXING and bad row messages print through. */
  fflush(stderr);
  int saved_err = dup(STDERR_FILENO);
  int devnull = open("/dev/null", O_WRONLY);
  dup2(devnull, STDERR_FILENO);
  close(devnull);

  Cbc_solve(m);

  fflush(stderr);
  dup2(saved_err, STDERR_FILENO);
  close(saved_err);

  printf("\n  debugCuts result: %s  obj=%.6g  (certified=%.6g)\n",
         rk_name(classify(m, cert)), Cbc_getObjValue(m), cert);
  fflush(stdout);
  Cbc_deleteModel(m);
}

/* ── Usage ───────────────────────────────────────────────────── */

static void usage(const char *prog)
{
  fprintf(stderr,
    "Usage: %s <problem.{mps,lp}[.gz]> -opt <value> [OPTIONS]\n"
    "\n"
    "Diagnose wrong MIPster results (wrong-optimal or wrong-infeasible).\n"
    "By default: confirms the bug, then runs a -debugCuts pass that prints\n"
    "'bad row' (invalid row cut) and 'nodeBoundProp BAD FIXING' (wrong bound\n"
    "fixing) messages — enough to identify most soundness bugs.\n"
    "\n"
    "Required:\n"
    "  -opt N        Certified optimal value\n"
    "\n"
    "Options:\n"
    "  -sol FILE     Known-good .sol file for the debugCuts pass\n"
    "                (if omitted the tool attempts to auto-obtain one)\n"
    "  -time N       Per-run wall-clock time limit in seconds (default: 60)\n"
    "  -max          Maximisation problem (default: minimisation)\n"
    "  -scan         Also run %d feature-disabling configs to narrow down\n"
    "                the culprit (slow: up to %d × time-limit extra seconds)\n"
    "  -p PAR VAL    Extra solver parameter applied to every run (repeatable).\n"
    "                Use to reproduce bugs triggered by specific settings.\n"
    "\n"
    "Quick workflow:\n"
    "  # Default: fast debugCuts pass (a few minutes at most)\n"
    "  %s problem.mps -opt <N> [-sol /tmp/ref.sol] [-time 120]\n"
    "  %s problem.mps -opt <N> -p nodeBoundProp on -p nodeBoundPropMinDepth 2\n"
    "  # Full feature scan if debugCuts is not conclusive:\n"
    "  %s problem.mps -opt <N> -sol /tmp/ref.sol -scan\n",
    prog, N_CONFIGS, N_CONFIGS, prog, prog, prog);
}

/* ── main ────────────────────────────────────────────────────── */

int main(int argc, char **argv)
{
  const char *filename   = nullptr;
  const char *sol_file   = nullptr;
  double      cert       = 0.0;
  bool        has_opt    = false;
  int         time_limit = 60;
  bool        is_max     = false;
  bool        do_scan    = false;
  std::vector<ExtraParam> extra;

  for (int i = 1; i < argc; i++) {
    if (!strcmp(argv[i], "-opt") && i + 1 < argc)
      { cert = atof(argv[++i]); has_opt = true; }
    else if (!strcmp(argv[i], "-sol") && i + 1 < argc)
      sol_file = argv[++i];
    else if (!strcmp(argv[i], "-time") && i + 1 < argc)
      time_limit = atoi(argv[++i]);
    else if (!strcmp(argv[i], "-max"))
      is_max = true;
    else if (!strcmp(argv[i], "-scan"))
      do_scan = true;
    else if (!strcmp(argv[i], "-p") && i + 2 < argc) {
      ExtraParam ep = { argv[i + 1], argv[i + 2] };
      extra.push_back(ep);
      i += 2;
    }
    else if (!strcmp(argv[i], "-h") || !strcmp(argv[i], "--help"))
      { usage(argv[0]); return 0; }
    else if (argv[i][0] != '-')
      filename = argv[i];
    else {
      fprintf(stderr, "Unknown option: %s\n\n", argv[i]);
      usage(argv[0]);
      return 1;
    }
  }

  if (!filename || !has_opt) { usage(argv[0]); return 1; }

  printf("mipster_diag\n");
  printf("  Problem:       %s\n", filename);
  printf("  Certified opt: %.6g  (%s)\n", cert, is_max ? "max" : "min");
  printf("  Time limit:    %ds per run\n", time_limit);
  if (sol_file) printf("  Sol file:      %s\n", sol_file);
  if (do_scan)  printf("  Feature scan:  enabled (%d configs)\n", N_CONFIGS);
  for (const auto &ep : extra)
    printf("  Extra param:   %s %s\n", ep.param, ep.value);

  Opts base = { filename, is_max ? 1 : 0, time_limit, /*log_level=*/0, extra };

  /* ── Phase 1: Baseline — confirm the bug ── */
  printf("\n  Baseline run (default settings, %ds)...\n", time_limit);
  fflush(stdout);
  {
    Cbc_Model *m = make_model(base);
    if (!m) return 1;
    silent_solve(m);
    ResultKind rk = classify(m, cert);
    double obj    = Cbc_getObjValue(m);
    double bound  = Cbc_getBestPossibleObjValue(m);
    Cbc_deleteModel(m);

    if (rk == RK_OK) {
      printf("  Baseline: CORRECT (obj=%.6g) — bug not reproduced within %ds.\n"
             "  Try a larger -time N.\n", obj, time_limit);
      return 0;
    }
    switch (rk) {
    case RK_WRONG_OPT:
      printf("  Baseline: WRONG_OPT  obj=%.6g  certified=%.6g — confirmed.\n",
             obj, cert);
      break;
    case RK_INFEASIBLE:
      printf("  Baseline: WRONG_INFEASIBLE — confirmed.\n");
      break;
    case RK_LEAD:
      printf("  Baseline: NOSOL/partial  obj=%.6g  bound=%.6g — "
             "may need more time.\n", obj, bound);
      break;
    default:
      printf("  Baseline: %s  bound=%.6g — proceeding.\n", rk_name(rk), bound);
      break;
    }
  }

  char scan_sol_path[256] = "";
  snprintf(scan_sol_path, sizeof(scan_sol_path),
           "/tmp/mipster_diag_%d.sol", (int)getpid());
  bool cleanup_sol = false;

  /* ── Phase 2: Auto-obtain reference solution if not provided ── */
  if (!sol_file) {
    cleanup_sol = auto_obtain_sol(base, cert, scan_sol_path);
    if (cleanup_sol)
      sol_file = scan_sol_path;
  }

  /* ── Phase 3: debugCuts pass ── */
  if (sol_file)
    run_debug_cuts(base, sol_file, cert);
  else
    printf("\n  debugCuts: skipped — no reference solution available.\n"
           "  Run: mipster problem.mps -preprocess off -solve -solu /tmp/ref.sol\n"
           "  Then rerun with -sol /tmp/ref.sol\n");

  /* ── Phase 4 (optional): Feature scan ── */
  if (do_scan) {
    bool have_scan_sol = run_feature_scan(base, cert,
                                          sol_file ? nullptr : scan_sol_path);
    if (!sol_file && have_scan_sol) {
      sol_file = scan_sol_path;
      cleanup_sol = true;
    }
  } else {
    printf("\n  Tip: if debugCuts output above is not conclusive, rerun with -scan\n"
           "  to run %d feature-disabling configurations.\n", N_CONFIGS);
  }

  if (cleanup_sol)
    unlink(scan_sol_path);

  return 0;
}
