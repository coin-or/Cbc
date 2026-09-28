#!/usr/bin/env python3
"""
analyze_lp_params.py — Analyse LP relaxation parameter tuning experiments.

Reads lp_results.csv, averages over random seeds, then produces:
  - Per-param overall rankings (shifted geomean, mean time, solve rate)
  - Which params are fastest most often
  - Params substantially faster than a baseline for ≥N instances
  - Error/timeout summary: which instances and params had failures
  - Error classification: TRIVIALLY_OPTIMAL, CRASH (SIGSEGV/SIGABRT), WRONG_RESULT

Error categories:
  TRIVIALLY_OPTIMAL  STATUS=ERROR in CSV but check_result=yes — the LP was solved
                     without any iterations (problem trivially feasible at LP root);
                     the solver wrote a valid solution but printed no "✔ Optimal" line.
  CRASH_SIGSEGV      Solver received SIGSEGV — log contains "Signal SIGSEGV caught".
  CRASH_SIGABRT      Solver received SIGABRT — log contains "Signal SIGABRT caught".
  CRASH_OTHER        Non-zero exit code with no check_result and no recognised signal.
  WRONG_RESULT       STATUS=OPTIMAL but objective value exceeds the best known LP
                     relaxation value by more than --wrong-tol (relative).

Usage:
    python3 analyze_lp_params.py --dir <experiment_dir> [OPTIONS]

Options:
    --dir DIR            Experiment directory (must contain lp_results.csv)
    --timelimit T        LP time limit in seconds (default: 7200)
    --penalty-mult M     Penalty multiplier for failures (default: 2.0)
    --shift S            Shift for shifted geometric mean in seconds (default: 1.0)
    --speedup-thresh T   Comma-separated speedup ratios for "substantially better"
                         (default: 1.1,1.25,1.5,2.0)
    --min-instances N    Comma-separated minimum-instance thresholds
                         (default: 10,50,100)
    --baseline PARAM     Baseline param tag for comparisons (default: dual_default)
    --wrong-tol TOL      Relative tolerance for flagging wrong results (default: 1e-4)
    --out FILE           Output report file (default: lp_analysis.txt in exp dir)
    --csv-out FILE       Per-instance-param avg times CSV
                         (default: lp_avg_times.csv in exp dir)
"""

import argparse
import csv
import math
import os
import re
import sys
from collections import defaultdict


# ── Helpers ───────────────────────────────────────────────────────────────────

def shifted_geomean(values, shift):
    """exp(mean(log(v + shift))) - shift.  Returns nan for empty input."""
    if not values:
        return float("nan")
    return math.exp(sum(math.log(v + shift) for v in values) / len(values)) - shift


def pct(num, denom):
    return 100.0 * num / denom if denom else float("nan")


def fmt_nan(v, fmt=".3f"):
    return format(v, fmt) if not math.isnan(v) else "N/A"


# ── Output helpers ────────────────────────────────────────────────────────────

class Report:
    def __init__(self):
        self._lines = []

    def raw(self, s=""):
        self._lines.append(s)

    def h1(self, title):
        self._lines += ["", "=" * 74, f"  {title}", "=" * 74]

    def h2(self, title):
        self._lines += ["", f"── {title} " + "─" * max(0, 70 - len(title))]

    def table(self, headers, rows, fmt=None):
        if fmt is None:
            fmt = ["<"] * len(headers)
        widths = [len(h) for h in headers]
        for row in rows:
            for i, cell in enumerate(row):
                widths[i] = max(widths[i], len(str(cell)))
        sep = "  ".join("-" * w for w in widths)
        hdr = "  ".join(f"{h:{fmt[i]}{widths[i]}}" for i, h in enumerate(headers))
        self._lines.append(hdr)
        self._lines.append(sep)
        for row in rows:
            self._lines.append("  ".join(
                f"{str(c):{fmt[i]}{widths[i]}}" for i, c in enumerate(row)))

    def text(self):
        return "\n".join(self._lines)


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--dir", required=True)
    parser.add_argument("--timelimit", type=float, default=7200.0)
    parser.add_argument("--penalty-mult", type=float, default=2.0)
    parser.add_argument("--shift", type=float, default=1.0)
    parser.add_argument("--speedup-thresh", default="1.1,1.25,1.5,2.0")
    parser.add_argument("--min-instances", default="10,50,100")
    parser.add_argument("--baseline", default="dual_default")
    parser.add_argument("--wrong-tol", type=float, default=1e-4,
                        help="Relative tolerance for flagging wrong results")
    parser.add_argument("--out", default=None)
    parser.add_argument("--csv-out", default=None)
    parser.add_argument("--exclude-prefix", default="guess_",
                        help="Exclude param tags starting with this prefix (default: guess_)")
    args = parser.parse_args()

    exp_dir = args.dir
    csv_file = os.path.join(exp_dir, "lp_results.csv")
    if not os.path.isfile(csv_file):
        sys.exit(f"Error: {csv_file} not found")

    penalty = args.timelimit * args.penalty_mult
    speedup_thresholds = [float(x) for x in args.speedup_thresh.split(",")]
    min_inst_thresholds = [int(x) for x in args.min_instances.split(",")]
    out_file = args.out or os.path.join(exp_dir, "lp_analysis.txt")
    csv_out_file = args.csv_out or os.path.join(exp_dir, "lp_avg_times.csv")

    # ── Load CSV ──────────────────────────────────────────────────────────────
    # Columns: instance, param_tag, seed, status, obj_from_log, obj_from_check,
    #          wall_seconds, check_result, exit_code
    raw = defaultdict(list)   # (instance, param_tag) -> [(seed, time, status)]
    instances_set = set()
    params_set = set()
    all_rows = []
    exclude_prefix = args.exclude_prefix or ""

    with open(csv_file, newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            inst = row["instance"]
            param = row["param_tag"]
            if exclude_prefix and param.startswith(exclude_prefix):
                continue
            seed = row["seed"]
            status = row["status"]
            try:
                wall = float(row["wall_seconds"])
            except (ValueError, KeyError):
                wall = penalty
            raw[(inst, param)].append((seed, wall, status))
            instances_set.add(inst)
            params_set.add(param)
            all_rows.append(row)

    instances = sorted(instances_set)
    params = sorted(params_set)
    n_inst = len(instances)
    n_params = len(params)
    n_total = len(all_rows)

    # ── Classify errors by inspecting log files ───────────────────────────────
    # Categories:
    #   TRIVIALLY_OPTIMAL  — check_result=yes despite STATUS=ERROR; no iterations needed
    #   CRASH_SIGSEGV      — "Signal SIGSEGV caught" in log
    #   CRASH_SIGABRT      — "Signal SIGABRT caught" in log
    #   CRASH_OTHER        — non-zero exit, no check_result, no recognised signal
    #   WRONG_RESULT       — STATUS=OPTIMAL but obj > best-known by > wrong_tol

    print("Classifying errors (reading log files for crash cases)...",
          file=sys.stderr, flush=True)

    # First pass: build best-known objective per instance from valid OPTIMAL runs
    best_obj = {}   # instance -> lowest valid obj_from_check
    for row in all_rows:
        if row["status"] != "OPTIMAL":
            continue
        chk = row.get("check_result", "")
        if not chk.startswith("yes"):
            continue
        try:
            obj = float(row["obj_from_check"])
        except (ValueError, TypeError):
            continue
        inst = row["instance"]
        if inst not in best_obj or obj < best_obj[inst]:
            best_obj[inst] = obj

    # Second pass: classify each row
    error_class = {}   # (inst, param, seed) -> category string

    for row in all_rows:
        inst  = row["instance"]
        param = row["param_tag"]
        seed  = row["seed"]
        status = row["status"]
        chk   = row.get("check_result", "")
        key   = (inst, param, seed)
        try:
            exit_code = int(row.get("exit_code", 0))
        except (ValueError, TypeError):
            exit_code = 0

        if status == "OPTIMAL":
            # Check for wrong result
            try:
                obj = float(row["obj_from_check"])
            except (ValueError, TypeError):
                error_class[key] = "OK"
                continue
            b = best_obj.get(inst)
            if b is not None:
                tol_abs = args.wrong_tol * max(1.0, abs(b))
                if obj > b + tol_abs:
                    rel = (obj - b) / max(abs(b), 1e-10)
                    error_class[key] = f"WRONG_RESULT(rel={rel:.2e})"
                else:
                    error_class[key] = "OK"
            else:
                error_class[key] = "OK"

        elif status == "ERROR":
            if chk.startswith("yes"):
                # Solved trivially — no LP iterations, no "✔ Optimal" line
                error_class[key] = "TRIVIALLY_OPTIMAL"
            else:
                # Read log file to detect crash signal
                prefix = f"{inst}_{param}_s{seed}_fpp"
                log_path = os.path.join(exp_dir, f"{prefix}.log")
                crash_cat = "CRASH_OTHER"
                if os.path.isfile(log_path):
                    try:
                        with open(log_path, errors="replace") as lf:
                            log_text = lf.read()
                        if "Signal SIGSEGV caught" in log_text:
                            crash_cat = "CRASH_SIGSEGV"
                        elif "Signal SIGABRT caught" in log_text:
                            crash_cat = "CRASH_SIGABRT"
                    except OSError:
                        pass
                error_class[key] = crash_cat

        elif status in ("TIMEOUT", "TIMEOUT_KILLED"):
            error_class[key] = "TIMEOUT"
        else:
            error_class[key] = status

    # Aggregate error classes per (inst, param) — for multi-seed view
    # Also remap raw STATUS=ERROR to real category for avg_time computation:
    # TRIVIALLY_OPTIMAL rows should be treated as OPTIMAL for timing purposes.

    # ── Compute per-(instance, param) averages ────────────────────────────────
    # TRIVIALLY_OPTIMAL treated as OPTIMAL; all other non-OPTIMAL → penalty.
    # avg_time   = mean across all seeds (including penalised).
    # all_solved = True if every seed was OPTIMAL or TRIVIALLY_OPTIMAL.
    # any_solved = True if at least one seed was solved.

    avg_time   = {}
    all_solved = {}
    any_solved = {}

    for (inst, param), runs in raw.items():
        times, n_opt = [], 0
        for seed, wall, status in runs:
            eclass = error_class.get((inst, param, seed), status)
            # WRONG_RESULT has status=OPTIMAL but objective outside tolerance —
            # treat as failure (penalise) because the result cannot be trusted.
            actually_solved = (
                (status == "OPTIMAL" and not eclass.startswith("WRONG"))
                or eclass == "TRIVIALLY_OPTIMAL"
            )
            if actually_solved:
                times.append(wall)
                n_opt += 1
            else:
                times.append(penalty)
        avg_time[(inst, param)]   = sum(times) / len(times)
        all_solved[(inst, param)] = (n_opt == len(runs))
        any_solved[(inst, param)] = (n_opt > 0)

    # Write per-instance-param averages
    with open(csv_out_file, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["instance", "param_tag", "avg_wall_seconds",
                    "all_seeds_solved", "any_seed_solved"])
        for inst in instances:
            for param in params:
                key = (inst, param)
                if key in avg_time:
                    w.writerow([inst, param,
                                 f"{avg_time[key]:.4f}",
                                 all_solved.get(key, False),
                                 any_solved.get(key, False)])

    # ── Per-param summary stats ───────────────────────────────────────────────
    param_stats = {}
    for param in params:
        times_all = []
        n_all_s = n_any_s = 0
        for inst in instances:
            key = (inst, param)
            t = avg_time.get(key, penalty)
            times_all.append(t)
            if all_solved.get(key, False):
                n_all_s += 1
            if any_solved.get(key, False):
                n_any_s += 1
        param_stats[param] = {
            "n_solved":    n_all_s,
            "n_any_solved": n_any_s,
            "sgm_all":     shifted_geomean(times_all, args.shift),
            "mean_all":    sum(times_all) / len(times_all),
        }

    # ── Best-param per instance ───────────────────────────────────────────────
    best_param_counts = defaultdict(int)
    best_param_map = {}
    for inst in instances:
        best = min(params, key=lambda p: (avg_time.get((inst, p), penalty), p))
        best_param_counts[best] += 1
        best_param_map[inst] = best

    # ── Substantially-better-than-baseline counts ────────────────────────────
    baseline = args.baseline
    if baseline not in params_set:
        print(f"Warning: baseline '{baseline}' not found; using '{params[0]}'",
              file=sys.stderr)
        baseline = params[0]

    speedup_counts = defaultdict(lambda: defaultdict(int))
    for inst in instances:
        base_t = avg_time.get((inst, baseline), penalty)
        for param in params:
            pt = avg_time.get((inst, param), penalty)
            speedup = base_t / pt if pt > 0 else float("inf")
            for thr in speedup_thresholds:
                if speedup >= thr:
                    speedup_counts[param][thr] += 1

    # ── Pre-compute category summary stats (needed for report header) ──────────
    cat_counts_pre   = defaultdict(int)
    wrong_detail_pre = []
    n_triv_pre       = 0
    for row in all_rows:
        cat = error_class.get((row["instance"], row["param_tag"], row["seed"]),
                              row["status"])
        cat_counts_pre[cat] += 1
        if cat.startswith("WRONG"):
            wrong_detail_pre.append(row)
        if cat == "TRIVIALLY_OPTIMAL":
            n_triv_pre += 1
    n_crashes_hdr  = sum(cat_counts_pre[c] for c in cat_counts_pre if c.startswith("CRASH"))
    n_wrong_hdr    = len(wrong_detail_pre)
    n_timeouts_hdr = (cat_counts_pre.get("TIMEOUT", 0) +
                      cat_counts_pre.get("TIMEOUT_KILLED", 0))

    # ── Build report ──────────────────────────────────────────────────────────
    rpt = Report()

    rpt.h1(f"LP Relaxation Parameter Analysis  —  {os.path.basename(exp_dir)}")
    rpt.raw(f"  Instances  : {n_inst}")
    rpt.raw(f"  Parameters : {n_params}")
    rpt.raw(f"  Total runs : {n_total}  (seeds: {sorted({r['seed'] for r in all_rows})})")
    rpt.raw(f"  Time limit : {args.timelimit:.0f}s  |  Penalty (failures): {penalty:.0f}s")
    rpt.raw(f"  Baseline   : {baseline}")
    rpt.raw(f"  SGM shift  : {args.shift}s")
    rpt.raw(f"  Crashes    : {n_crashes_hdr}"
            f"  |  Wrong results: {n_wrong_hdr}"
            f"  |  Timeouts: {n_timeouts_hdr}"
            f"  |  Trivially-optimal: {n_triv_pre}")

    # ── TABLE 1: Overall ranking — SGM over all instances ────────────────────
    rpt.h2("Overall ranking — shifted geomean over ALL instances (penalty for failures)")
    rpt.raw(f"  Shift = {args.shift}s.  Lower is better.  Timeouts/failures → penalty time.")
    rpt.raw("")
    sorted_all = sorted(params, key=lambda p: param_stats[p]["sgm_all"])
    rows = []
    for rank, param in enumerate(sorted_all, 1):
        s = param_stats[param]
        rows.append((rank, param,
                     fmt_nan(s["sgm_all"]),
                     s["n_solved"],
                     f"{pct(s['n_solved'], n_inst):.1f}%",
                     best_param_counts.get(param, 0)))
    rpt.table(
        ["Rank", "Param", "SGM(all)s", "#Solved", "%Solved", "#Best"],
        rows,
        fmt=[">", "<", ">", ">", ">", ">"])

    # ── TABLE 2: How often each param is the fastest ──────────────────────────
    rpt.h2("How often each param is the FASTEST (lowest avg time, incl. penalties)")
    rpt.raw("")
    rows = [(p, best_param_counts.get(p, 0),
             f"{pct(best_param_counts.get(p, 0), n_inst):.1f}%")
            for p in sorted(params, key=lambda p: -best_param_counts.get(p, 0))]
    rpt.table(["Param", "#Best", "%Best"], rows, fmt=["<", ">", ">"])

    # ── TABLE 3: Substantially faster than baseline ───────────────────────────
    rpt.h2(f"Count of instances where param is SUBSTANTIALLY FASTER than {baseline}")
    rpt.raw(f"  Speedup = {baseline}_avg_time / param_avg_time  ≥ threshold")
    rpt.raw("")
    thr_hdrs = [f"≥{t:.2f}×" for t in speedup_thresholds]
    fmt = ["<"] + [">"] * len(speedup_thresholds)
    rows = []
    for param in sorted(params, key=lambda p: -speedup_counts[p][speedup_thresholds[0]]):
        rows.append([param] + [speedup_counts[param][thr] for thr in speedup_thresholds])
    rpt.table(["Param"] + thr_hdrs, rows, fmt=fmt)

    # ── TABLE 4: Params with wins ≥ N instances per threshold ────────────────
    for thr in speedup_thresholds:
        for min_inst in min_inst_thresholds:
            candidates = [(p, speedup_counts[p][thr])
                          for p in params if speedup_counts[p][thr] >= min_inst]
            if not candidates:
                continue
            rpt.h2(f"Params with ≥{min_inst} instances at speedup ≥{thr:.2f}× vs {baseline}")
            rpt.raw("")
            candidates.sort(key=lambda x: -x[1])
            rows = [(p, cnt, f"{pct(cnt, n_inst):.1f}%",
                     fmt_nan(param_stats[p]["sgm_all"]),
                     param_stats[p]["n_solved"])
                    for p, cnt in candidates]
            rpt.table(["Param", "#Instances", "%Total", "SGM(all)s", "#Solved"],
                      rows, fmt=["<", ">", ">", ">", ">"])

    # ── TABLE 5: Best-param speedup distribution vs baseline ─────────────────
    rpt.h2(f"Distribution of best-param speedup vs {baseline} across instances")
    rpt.raw("  Best param = fastest avg time per instance.")
    rpt.raw("")
    dist_thrs = [1.0, 1.1, 1.25, 1.5, 2.0, 3.0, 5.0, 10.0]
    dist_counts = defaultdict(int)
    for inst in instances:
        best_t = avg_time.get((inst, best_param_map[inst]), penalty)
        base_t = avg_time.get((inst, baseline), penalty)
        sp = base_t / best_t if best_t > 0 else float("inf")
        for thr in dist_thrs:
            if sp >= thr:
                dist_counts[thr] += 1
    rows = [(f"≥{t:.1f}×", dist_counts[t], f"{pct(dist_counts[t], n_inst):.1f}%")
            for t in dist_thrs]
    rpt.table(["Speedup threshold", "#Instances", "%Total"], rows, fmt=["<", ">", ">"])

    # ── TABLE 6: Head-to-head vs baseline ─────────────────────────────────────
    for ref in ([baseline] + (["dual_default"] if baseline != "dual_default" and
                               "dual_default" in params_set else [])):
        rpt.h2(f"Head-to-head: each param vs {ref}  (geometric mean time ratio)")
        rpt.raw("  ratio < 1.0 = param FASTER;  > 1.0 = param SLOWER")
        rpt.raw("  Computed only over instances where BOTH solved in all seeds.")
        rpt.raw("")
        hh_rows = []
        for param in params:
            if param == ref:
                continue
            ratios, n_both, n_win, n_lose = [], 0, 0, 0
            for inst in instances:
                if (all_solved.get((inst, param), False) and
                        all_solved.get((inst, ref), False)):
                    pt = avg_time[(inst, param)]
                    rt = avg_time[(inst, ref)]
                    n_both += 1
                    ratio = pt / rt if rt > 0 else float("inf")
                    ratios.append(ratio)
                    if pt < rt:
                        n_win += 1
                    elif rt < pt:
                        n_lose += 1
            if not ratios:
                continue
            gmr = math.exp(sum(math.log(r) for r in ratios) / len(ratios))
            hh_rows.append((param, f"{gmr:.4f}", n_both,
                             n_win, n_lose, f"{pct(n_win, n_both):.1f}%"))
        hh_rows.sort(key=lambda r: float(r[1]))
        rpt.table(
            ["Param", "GeoMean ratio", "#Both solved", "#Param wins",
             f"#{ref} wins", "%Param wins"],
            hh_rows, fmt=["<", ">", ">", ">", ">", ">"])

    # ── DOMINANCE / WORTHLESS PARAM ANALYSIS ─────────────────────────────────
    rpt.h1("Dominance Analysis")
    rpt.raw("")
    rpt.raw("  Identifies parameter settings that add little or no value to the racing")
    rpt.raw("  portfolio: every instance they would win is already covered by another param.")

    # Build numpy-style structures from avg_time dict for efficient computation
    import math as _math
    param_list = list(params)
    inst_list  = list(instances)
    ni, np_ = len(inst_list), len(param_list)
    T_mat = [[avg_time.get((inst, p), penalty) for p in param_list]
             for inst in inst_list]  # list[ni][np_]

    def _sgm_vec(col_min):
        return math.exp(sum(math.log(v + 1.0) for v in col_min) / len(col_min)) - 1.0

    # Greedy portfolio construction
    selected_idx = []
    remaining_idx = list(range(np_))
    cur_min = [penalty] * ni
    greedy_curve = []
    while remaining_idx:
        best_sgm, best_j = float("inf"), -1
        for j in remaining_idx:
            candidate = [min(cur_min[i], T_mat[i][j]) for i in range(ni)]
            s = _sgm_vec(candidate)
            if s < best_sgm:
                best_sgm, best_j = s, j
        selected_idx.append(best_j)
        remaining_idx.remove(best_j)
        cur_min = [min(cur_min[i], T_mat[i][best_j]) for i in range(ni)]
        greedy_curve.append((len(selected_idx), param_list[best_j], best_sgm))

    rpt.h2("Greedy portfolio construction — marginal value of each added param")
    rpt.raw("")
    rpt.raw("  K=1 is the single best param. Each subsequent row adds the param that")
    rpt.raw("  gives the largest SGM improvement to the racing portfolio (min over K).")
    rpt.raw("")
    gcurve_rows = []
    prev_sgm = None
    for k, name, s in greedy_curve:
        delta = s - prev_sgm if prev_sgm is not None else 0.0
        pct_gain = abs(delta) / prev_sgm * 100.0 if prev_sgm else 0.0
        gcurve_rows.append((k, name, f"{s:.3f}", f"{delta:+.3f}", f"{pct_gain:.1f}%"))
        prev_sgm = s
    rpt.table(["K", "Param added", "PortfolioSGM", "Delta", "%Gain"],
              gcurve_rows, fmt=[">", "<", ">", ">", ">"])

    # Marginal removal cost: SGM of full portfolio minus each param
    full_min = [min(T_mat[i]) for i in range(ni)]
    full_sgm_val = _sgm_vec(full_min)

    removal_cost = {}
    for j, p in enumerate(param_list):
        rest_min = [min(T_mat[i][k] for k in range(np_) if k != j) for i in range(ni)]
        removal_cost[p] = _sgm_vec(rest_min) - full_sgm_val

    # Domination fraction: % instances where some other param is ≥10% faster
    dom_frac = {}
    best_dominator = {}
    for j, p in enumerate(param_list):
        dominated_count = 0
        best_df, best_dn = 0.0, ""
        for k, q in enumerate(param_list):
            if k == j:
                continue
            wins = sum(1 for i in range(ni) if T_mat[i][k] <= T_mat[i][j] * 0.90)
            if wins > best_df * ni:
                best_df, best_dn = wins / ni, q
        dom_frac[p] = sum(
            1 for i in range(ni)
            if any(T_mat[i][k] <= T_mat[i][j] * 0.90 for k in range(np_) if k != j)
        ) / ni
        best_dominator[p] = (best_df, best_dn)

    # Unique wins: instances where param is the ONLY one that solved (below penalty)
    unique_wins_count = {}
    for j, p in enumerate(param_list):
        n_solvers = [sum(1 for k in range(np_) if T_mat[i][k] < penalty)
                     for i in range(ni)]
        unique_wins_count[p] = sum(
            1 for i in range(ni)
            if T_mat[i][j] < penalty and n_solvers[i] == 1
        )

    rpt.h2("Per-param dominance and portfolio contribution")
    rpt.raw("")
    rpt.raw("  Dom%      : % of instances where at least one other param is ≥10% faster.")
    rpt.raw("  BestDom%  : % dominated by the single strongest competitor.")
    rpt.raw("  UniqueWins: instances where this param is the ONLY one that solves.")
    rpt.raw("  RemovalSGM: SGM increase when param removed from full portfolio (higher = more valuable).")
    rpt.raw("")
    dom_rows = sorted(param_list, key=lambda p: removal_cost[p])
    dom_table_rows = []
    for p in dom_rows:
        bd_pct, bd_name = best_dominator[p]
        ind_sgm = param_stats[p]["sgm_all"]
        dom_table_rows.append((
            p,
            f"{fmt_nan(ind_sgm, '.1f')}",
            f"{dom_frac[p]*100:.1f}%",
            f"{bd_pct*100:.1f}%",
            bd_name,
            unique_wins_count[p],
            f"{removal_cost[p]:+.4f}s",
        ))
    rpt.table(
        ["Param", "IndSGM", "Dom%", "BestDom%", "BestDominator", "UniqueWins", "RemovalSGM"],
        dom_table_rows,
        fmt=["<", ">", ">", ">", "<", ">", ">"])

    # Worthless params: removal cost < 0.01s (absolute) — negligible contribution
    worthless_thresh = 0.01
    worthless_params = [p for p in param_list if removal_cost[p] < worthless_thresh]
    worthless_params.sort(key=lambda p: removal_cost[p])

    rpt.h2(f"Worthless params — removal from full portfolio costs < {worthless_thresh:.3f}s")
    rpt.raw("")
    if worthless_params:
        rpt.raw(f"  These {len(worthless_params)} params are fully covered by the rest of the portfolio.")
        rpt.raw(f"  Removing them has negligible effect on racing performance.")
        rpt.raw(f"  Consider removing from lp_params.txt to reduce experiment cost.")
        rpt.raw("")
        w_rows = [(p, f"{removal_cost[p]:+.4f}s",
                   f"{param_stats[p]['sgm_all']:.1f}s",
                   f"{dom_frac[p]*100:.1f}%",
                   best_dominator[p][1])
                  for p in worthless_params]
        rpt.table(["Param", "RemovalCost", "IndSGM", "Dom%", "PrimaryDominator"],
                  w_rows, fmt=["<", ">", ">", ">", "<"])
    else:
        rpt.raw("  (no clearly worthless params found)")

    # ── ERROR / TIMEOUT SECTION ───────────────────────────────────────────────
    rpt.h1("Errors, Crashes, and Wrong Results")

    # Build aggregated counts per category
    cat_counts = defaultdict(int)          # category -> total run count
    cat_by_param = defaultdict(lambda: defaultdict(int))   # param -> cat -> count
    cat_by_inst  = defaultdict(lambda: defaultdict(int))   # inst  -> cat -> count
    # For crash/wrong-result detail: list of (inst, param, seed, category)
    crash_detail  = []
    wrong_detail  = []
    triv_by_param = defaultdict(int)       # param -> count of TRIVIALLY_OPTIMAL

    for row in all_rows:
        inst  = row["instance"]
        param = row["param_tag"]
        seed  = row["seed"]
        cat   = error_class.get((inst, param, seed), row["status"])
        cat_counts[cat] += 1
        cat_by_param[param][cat] += 1
        cat_by_inst[inst][cat] += 1
        if cat.startswith("CRASH"):
            crash_detail.append((inst, param, seed, cat))
        elif cat.startswith("WRONG"):
            wrong_detail.append((inst, param, seed, cat))
        elif cat == "TRIVIALLY_OPTIMAL":
            triv_by_param[param] += 1

    n_crashes = sum(1 for c in cat_counts if c.startswith("CRASH")
                    for _ in range(cat_counts[c]))
    n_wrong   = len(wrong_detail)
    n_triv    = cat_counts["TRIVIALLY_OPTIMAL"]
    n_timeouts = cat_counts.get("TIMEOUT", 0) + cat_counts.get("TIMEOUT_KILLED", 0)

    # ── Header stats ──────────────────────────────────────────────────────────
    rpt.h2("Summary counts by category")
    rpt.raw("")
    rpt.raw("  Categories:")
    rpt.raw("    OK                — solved OPTIMAL, objective matches best known")
    rpt.raw("    TRIVIALLY_OPTIMAL — STATUS=ERROR in CSV but solution is valid; LP solved")
    rpt.raw("                        without iterations (no '✔ Optimal' line emitted)")
    rpt.raw("    CRASH_SIGSEGV     — solver received SIGSEGV (segmentation fault)")
    rpt.raw("    CRASH_SIGABRT     — solver received SIGABRT (assertion failure / abort)")
    rpt.raw("    CRASH_OTHER       — non-zero exit, no recognised signal")
    rpt.raw("    WRONG_RESULT      — OPTIMAL reported but obj > best-known by "
            f">{args.wrong_tol:.0e} rel")
    rpt.raw("    TIMEOUT(_KILLED)  — hit time limit")
    rpt.raw("")
    rpt.table(["Category", "Count", "%Total"],
              [(cat, cnt, f"{pct(cnt, n_total):.2f}%")
               for cat, cnt in sorted(cat_counts.items(), key=lambda x: -x[1])],
              fmt=["<", ">", ">"])

    # ── CRASHES ───────────────────────────────────────────────────────────────
    rpt.h2("Crashes — per param")
    rpt.raw("")
    crash_param_rows = []
    for param in sorted(params):
        n_segv  = cat_by_param[param].get("CRASH_SIGSEGV", 0)
        n_abort = cat_by_param[param].get("CRASH_SIGABRT", 0)
        n_other = cat_by_param[param].get("CRASH_OTHER", 0)
        total   = n_segv + n_abort + n_other
        if total == 0:
            continue
        crash_param_rows.append((param, n_segv, n_abort, n_other, total))
    if crash_param_rows:
        crash_param_rows.sort(key=lambda r: -r[4])
        rpt.table(["Param", "#SIGSEGV", "#SIGABRT", "#Other", "#Total"],
                  crash_param_rows, fmt=["<", ">", ">", ">", ">"])
    else:
        rpt.raw("  (no crashes)")

    rpt.h2("Crashes — per instance (instances with ≥1 crash)")
    rpt.raw("")
    crash_inst_rows = []
    for inst in sorted(instances):
        n_segv  = cat_by_inst[inst].get("CRASH_SIGSEGV", 0)
        n_abort = cat_by_inst[inst].get("CRASH_SIGABRT", 0)
        n_other = cat_by_inst[inst].get("CRASH_OTHER", 0)
        total   = n_segv + n_abort + n_other
        if total == 0:
            continue
        # Which params crash on this instance?
        crashing_params = sorted({p for (i, p, s, c) in crash_detail
                                   if i == inst and c.startswith("CRASH")})
        crash_inst_rows.append((inst, n_segv, n_abort, n_other, total,
                                 ", ".join(crashing_params)))
    if crash_inst_rows:
        crash_inst_rows.sort(key=lambda r: -r[4])
        rpt.table(["Instance", "#SIGSEGV", "#SIGABRT", "#Other", "#Total", "Crashing params"],
                  crash_inst_rows, fmt=["<", ">", ">", ">", ">", "<"])
    else:
        rpt.raw("  (no crashes)")

    rpt.h2("Crash detail: instance × param × seed × signal")
    rpt.raw("")
    if crash_detail:
        crash_detail.sort(key=lambda r: (r[3], r[0], r[1]))
        rpt.table(["Instance", "Param", "Seed", "Signal"],
                  crash_detail, fmt=["<", "<", "<", "<"])
    else:
        rpt.raw("  (none)")

    # ── WRONG RESULTS ─────────────────────────────────────────────────────────
    rpt.h2(f"Wrong results — OPTIMAL reported but obj > best-known by >{args.wrong_tol:.0e}")
    rpt.raw("")
    if wrong_detail:
        # Summarise by instance: how many runs wrong, what is the deviation?
        wrong_inst_summary = defaultdict(list)
        for inst, param, seed, cat in wrong_detail:
            wrong_inst_summary[inst].append((param, seed, cat))

        wrong_inst_rows = []
        for inst in sorted(wrong_inst_summary.keys()):
            entries = wrong_inst_summary[inst]
            # extract max relative error from category string "WRONG_RESULT(rel=X)"
            rels = []
            for _, _, cat in entries:
                m = re.search(r"rel=([\d.e+\-]+)", cat)
                if m:
                    try:
                        rels.append(float(m.group(1)))
                    except ValueError:
                        pass
            max_rel = max(rels) if rels else float("nan")
            affected_params = sorted({p for p, _, _ in entries})
            wrong_inst_rows.append((inst, len(entries), f"{max_rel:.2e}",
                                    f"{best_obj.get(inst, float('nan')):.6g}",
                                    ", ".join(affected_params)))
        rpt.table(["Instance", "#Wrong runs", "Max rel err", "Best known obj",
                   "Affected params (sample)"],
                  wrong_inst_rows, fmt=["<", ">", ">", ">", "<"])

        rpt.raw("")
        rpt.h2("Wrong result detail: instance × param × seed × deviation")
        rpt.raw("")
        wrong_detail_sorted = sorted(wrong_detail, key=lambda r: (r[0], r[1], r[2]))
        rpt.table(["Instance", "Param", "Seed", "Deviation"],
                  wrong_detail_sorted, fmt=["<", "<", "<", "<"])
    else:
        rpt.raw("  (none — all OPTIMAL results match best-known objective)")

    # ── TRIVIALLY OPTIMAL ─────────────────────────────────────────────────────
    rpt.h2("Trivially optimal runs (STATUS=ERROR in CSV but solution is valid)")
    rpt.raw("  These are solved at LP root without iterations — not genuine errors.")
    rpt.raw("")
    if n_triv:
        triv_inst = sorted({inst for (inst, p, s, c) in
                             [(r["instance"], r["param_tag"], r["seed"],
                               error_class.get((r["instance"], r["param_tag"], r["seed"]), ""))
                              for r in all_rows]
                             if c == "TRIVIALLY_OPTIMAL"})
        rpt.raw(f"  Total trivially-optimal runs : {n_triv}")
        rpt.raw(f"  Distinct instances affected  : {len(triv_inst)}")
        rpt.raw(f"  Instances: {', '.join(triv_inst)}")

        rpt.raw("")
        triv_rows = [(p, cnt, f"{pct(cnt, 3):.0f}% of seeds")
                     for p, cnt in sorted(triv_by_param.items(), key=lambda x: -x[1])]
        rpt.table(["Param", "#Trivially-optimal runs", ""],
                  triv_rows, fmt=["<", ">", "<"])
    else:
        rpt.raw("  (none)")

    # ── TIMEOUTS ──────────────────────────────────────────────────────────────
    rpt.h2("Timeouts — per param")
    rpt.raw("")
    timeout_param_rows = []
    for param in sorted(params):
        n_to = (cat_by_param[param].get("TIMEOUT", 0) +
                cat_by_param[param].get("TIMEOUT_KILLED", 0))
        if n_to == 0:
            continue
        timeout_param_rows.append((param, n_to, f"{pct(n_to, n_inst * 3):.1f}%"))
    if timeout_param_rows:
        timeout_param_rows.sort(key=lambda r: -r[1])
        rpt.table(["Param", "#Timeouts", "%Of runs"], timeout_param_rows,
                  fmt=["<", ">", ">"])
    else:
        rpt.raw("  (no timeouts)")

    rpt.h2("Timeouts — per instance")
    rpt.raw("")
    timeout_inst_rows = []
    for inst in sorted(instances):
        n_to = (cat_by_inst[inst].get("TIMEOUT", 0) +
                cat_by_inst[inst].get("TIMEOUT_KILLED", 0))
        if n_to == 0:
            continue
        timed_params = sorted({p for r in all_rows
                                if r["instance"] == inst and
                                r["status"] in ("TIMEOUT", "TIMEOUT_KILLED")
                                for p in [r["param_tag"]]})
        timeout_inst_rows.append((inst, n_to, ", ".join(timed_params)))
    if timeout_inst_rows:
        timeout_inst_rows.sort(key=lambda r: -r[1])
        rpt.table(["Instance", "#Timeouts", "Params with timeouts"],
                  timeout_inst_rows, fmt=["<", ">", "<"])
    else:
        rpt.raw("  (no timeouts)")

    # ── Write & print ─────────────────────────────────────────────────────────
    output = rpt.text()
    with open(out_file, "w") as f:
        f.write(output)
        f.write("\n")

    print(output)
    print(f"\n→ Report saved to  : {out_file}")
    print(f"→ Avg times CSV    : {csv_out_file}")


if __name__ == "__main__":
    main()
