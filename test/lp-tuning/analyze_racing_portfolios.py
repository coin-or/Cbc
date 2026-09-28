#!/usr/bin/env python3
"""
analyze_racing_portfolios.py — Racing LP portfolio analysis.

For the "racingLP" feature in CBC, multiple CLP solves run in parallel with
different parameter settings and the first to finish wins.  This script
simulates that behaviour: given a portfolio of K parameter configurations,
the effective time per instance is min(t_1, ..., t_K).

Reads lp_avg_times.csv (produced by analyze_lp_params.py) and finds the
best combination of K parameters to race together.

Strategy
--------
  • Greedy construction: start empty, at each step add the param that gives
    the lowest portfolio SGM over all instances.  This is O(K × P) evaluations,
    fast for any practical K.
  • Exhaustive search: for K ≤ --exhaustive-k (default 4), enumerate all
    C(P, K) combinations and report the true optimum.  Compare with greedy
    to check quality.  C(29, 4) ≈ 23 k combinations — runs in < 1 s with
    numpy.
  • Performance curve: table of SGM / mean / solve-rate for each portfolio
    size from 1 to --max-k.
  • Instance coverage map: for the greedy portfolio at each K, which param
    is actually fastest per instance.

Metrics
-------
  racing_time(inst, portfolio) = min(avg_time[inst][p]  for p in portfolio)
                                 (penalty if no param has data for this inst)
  SGM = shifted geometric mean over racing_time values across all instances
  solve_rate = fraction of instances where ≥1 param in portfolio solved
               (all_seeds_solved = True in avg_times CSV)

Usage
-----
  python3 analyze_racing_portfolios.py \\
      --dir /home/haroldo/experiments/cbc/lp_relax_2026_05_12_noblas/ \\
      [--timelimit 14400] [--penalty-mult 2.0] [--shift 1.0] \\
      [--max-k 8] [--exhaustive-k 4] [--baseline dual_default] \\
      [--out report.txt] [--csv-out portfolio_curve.csv]
"""

import argparse
import csv
import itertools
import math
import os
import sys
from collections import defaultdict

try:
    import numpy as np
    HAS_NUMPY = True
except ImportError:
    HAS_NUMPY = False


# ── Helpers ───────────────────────────────────────────────────────────────────

def shifted_geomean(values, shift):
    """exp(mean(log(v + shift))) - shift.  Returns nan for empty."""
    if not values:
        return float("nan")
    return math.exp(sum(math.log(v + shift) for v in values) / len(values)) - shift


def pct(num, denom):
    return 100.0 * num / denom if denom else float("nan")


def fmt_nan(v, fmt=".3f"):
    return format(v, fmt) if not math.isnan(v) else "N/A"


# ── Report builder ────────────────────────────────────────────────────────────

class Report:
    def __init__(self):
        self._lines = []

    def raw(self, s=""):
        self._lines.append(s)

    def h1(self, title):
        self._lines += ["", "=" * 76, f"  {title}", "=" * 76]

    def h2(self, title):
        self._lines += ["", f"── {title} " + "─" * max(0, 72 - len(title))]

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


# ── Portfolio evaluation ──────────────────────────────────────────────────────

def portfolio_metrics(inst_list, time_mat, solved_mat, portfolio_indices, shift, penalty):
    """
    Compute SGM, mean, and solve-rate for a portfolio (list of column indices).

    Parameters
    ----------
    time_mat    : list of lists, shape (n_inst, n_params) — avg times (penalised)
    solved_mat  : list of lists, shape (n_inst, n_params) — bool, all seeds solved
    portfolio_indices : list of int (column indices)
    """
    racing_times = []
    n_solved = 0
    for i in range(len(inst_list)):
        t_min = min(time_mat[i][p] for p in portfolio_indices)
        racing_times.append(t_min)
        if any(solved_mat[i][p] for p in portfolio_indices):
            n_solved += 1

    sgm = shifted_geomean(racing_times, shift)
    mean = sum(racing_times) / len(racing_times)
    solve_rate = pct(n_solved, len(inst_list))
    return sgm, mean, solve_rate, n_solved


def portfolio_metrics_np(time_np, solved_np, portfolio_indices, shift, penalty):
    """Numpy-accelerated variant for exhaustive search."""
    sub_t = time_np[:, portfolio_indices]
    sub_s = solved_np[:, portfolio_indices]
    racing_times = sub_t.min(axis=1)
    n_solved = int(sub_s.any(axis=1).sum())
    n_inst = time_np.shape[0]
    sgm = math.exp(float(np.log(racing_times + shift).mean())) - shift
    mean = float(racing_times.mean())
    solve_rate = pct(n_solved, n_inst)
    return sgm, mean, solve_rate, n_solved


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--dir", required=True,
                        help="Experiment directory (must contain lp_avg_times.csv)")
    parser.add_argument("--avg-csv", default=None,
                        help="Path to avg-times CSV (default: lp_avg_times.csv in --dir). "
                             "Use lp_avg_times_full.csv for the 70-param extended set.")
    parser.add_argument("--timelimit", type=float, default=7200.0,
                        help="LP time limit used in experiment (default: 7200)")
    parser.add_argument("--penalty-mult", type=float, default=2.0,
                        help="Penalty multiplier for failures (default: 2.0)")
    parser.add_argument("--shift", type=float, default=1.0,
                        help="Shift for shifted geometric mean in seconds (default: 1.0)")
    parser.add_argument("--max-k", type=int, default=8,
                        help="Maximum portfolio size for greedy construction (default: 8)")
    parser.add_argument("--exhaustive-k", type=int, default=4,
                        help="Maximum K for exhaustive optimal search (default: 4)")
    parser.add_argument("--baseline", default="dual_default",
                        help="Baseline single param for comparison (default: dual_default)")
    parser.add_argument("--out", default=None,
                        help="Output report file (default: lp_racing_analysis.txt in exp dir)")
    parser.add_argument("--csv-out", default=None,
                        help="Portfolio curve CSV (default: lp_racing_curve.csv in exp dir)")
    parser.add_argument("--exclude-prefix", default="guess_",
                        help="Exclude param tags starting with this prefix (default: guess_)")
    args = parser.parse_args()

    exp_dir = args.exp_dir = args.dir
    avg_csv = args.avg_csv if args.avg_csv else os.path.join(exp_dir, "lp_avg_times.csv")
    if not os.path.isfile(avg_csv):
        sys.exit(f"Error: {avg_csv} not found.  Run analyze_lp_params.py first.")

    penalty = args.timelimit * args.penalty_mult
    exclude_prefix = args.exclude_prefix or ""
    out_file = args.out or os.path.join(exp_dir, "lp_racing_analysis.txt")
    csv_out  = args.csv_out or os.path.join(exp_dir, "lp_racing_curve.csv")

    # ── Load avg times ────────────────────────────────────────────────────────
    # Columns: instance, param_tag, avg_wall_seconds, all_seeds_solved, any_seed_solved
    raw = {}        # (instance, param) -> (avg_time, all_solved, any_solved)
    instances_set = set()
    params_set    = set()

    with open(avg_csv, newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            inst  = row["instance"]
            param = row["param_tag"]
            if exclude_prefix and param.startswith(exclude_prefix):
                continue
            try:
                t = float(row["avg_wall_seconds"])
            except (ValueError, KeyError):
                t = penalty
            all_s = row.get("all_seeds_solved", "False").strip() == "True"
            any_s = row.get("any_seed_solved",  "False").strip() == "True"
            raw[(inst, param)] = (t, all_s, any_s)
            instances_set.add(inst)
            params_set.add(param)

    instances = sorted(instances_set)
    params    = sorted(params_set)
    n_inst    = len(instances)
    n_params  = len(params)

    param_idx = {p: i for i, p in enumerate(params)}
    inst_idx  = {inst: i for i, inst in enumerate(instances)}

    # Build matrices (n_inst × n_params), fill missing with penalty
    time_mat   = [[penalty] * n_params for _ in range(n_inst)]
    solved_mat = [[False]   * n_params for _ in range(n_inst)]

    for (inst, param), (t, all_s, any_s) in raw.items():
        r = inst_idx[inst]
        c = param_idx[param]
        time_mat[r][c]   = t
        solved_mat[r][c] = all_s   # use all_seeds_solved for conservatism

    # Numpy arrays for fast exhaustive search
    if HAS_NUMPY:
        time_np   = np.array(time_mat, dtype=np.float64)
        solved_np = np.array(solved_mat, dtype=bool)
        eval_fn   = lambda idxs: portfolio_metrics_np(time_np, solved_np, idxs,
                                                       args.shift, penalty)
    else:
        eval_fn   = lambda idxs: portfolio_metrics(instances, time_mat, solved_mat,
                                                    idxs, args.shift, penalty)
        print("Warning: numpy not available — exhaustive search will be slower.",
              file=sys.stderr)

    # Baseline single-param performance
    if args.baseline in param_idx:
        b_idx = param_idx[args.baseline]
        baseline_sgm, baseline_mean, baseline_sr, baseline_ns = eval_fn([b_idx])
    else:
        print(f"Warning: baseline '{args.baseline}' not found.", file=sys.stderr)
        baseline_sgm, baseline_mean, baseline_sr, baseline_ns = float("nan"), float("nan"), 0.0, 0

    # ── Greedy portfolio construction ─────────────────────────────────────────
    print("Building greedy portfolio...", file=sys.stderr, flush=True)

    greedy_portfolio = []   # list of param indices, in order added
    greedy_curve    = []    # (k, sgm, mean, solve_rate, n_solved, param_added, sgm_gain_pct)

    remaining = list(range(n_params))
    prev_sgm  = float("inf")

    for k in range(1, min(args.max_k, n_params) + 1):
        best_sgm, best_param_idx = float("inf"), None
        for p in remaining:
            candidate = greedy_portfolio + [p]
            sgm, _, _, _ = eval_fn(candidate)
            if sgm < best_sgm:
                best_sgm = sgm
                best_param_idx = p

        greedy_portfolio.append(best_param_idx)
        remaining.remove(best_param_idx)
        sgm, mean, sr, ns = eval_fn(greedy_portfolio)
        gain = (prev_sgm - sgm) / prev_sgm * 100 if prev_sgm < float("inf") else float("nan")
        prev_sgm = sgm
        greedy_curve.append((k, sgm, mean, sr, ns, params[best_param_idx], gain))

    # ── Exhaustive optimal search for small K ─────────────────────────────────
    max_ex = min(args.exhaustive_k, n_params)
    exhaustive_results = {}  # k -> (best_sgm, best_combo_indices)

    for k in range(1, max_ex + 1):
        n_combos = math.comb(n_params, k)
        print(f"Exhaustive K={k}: {n_combos} combinations...", file=sys.stderr, flush=True)
        best_sgm, best_combo = float("inf"), None
        for combo in itertools.combinations(range(n_params), k):
            sgm, _, _, _ = eval_fn(list(combo))
            if sgm < best_sgm:
                best_sgm = sgm
                best_combo = combo
        exhaustive_results[k] = (best_sgm, best_combo)

    # ── Per-instance coverage at each greedy portfolio size ───────────────────
    def coverage_map(portfolio_indices):
        """For each instance, which param is fastest in the portfolio?"""
        coverage = {}
        for i, inst in enumerate(instances):
            times = [(time_mat[i][p], p) for p in portfolio_indices]
            best_t, best_p = min(times)
            coverage[inst] = (params[best_p], best_t, solved_mat[i][best_p])
        return coverage

    # Build report
    rpt = Report()
    rpt.h1(f"Racing LP Portfolio Analysis  —  {os.path.basename(os.path.normpath(exp_dir))}")
    rpt.raw(f"  Instances      : {n_inst}")
    rpt.raw(f"  Params pool    : {n_params}  ({', '.join(params)})")
    rpt.raw(f"  Time limit     : {args.timelimit:.0f}s")
    rpt.raw(f"  Penalty        : {penalty:.0f}s (= {args.penalty_mult}× timelimit)")
    rpt.raw(f"  SGM shift      : {args.shift}s")
    rpt.raw(f"  Exhaustive up  : K ≤ {max_ex}  (optimal guarantees)")
    rpt.raw(f"  Greedy up to   : K ≤ {len(greedy_curve)}")
    rpt.raw(f"  Baseline param : {args.baseline}  "
            f"(SGM={fmt_nan(baseline_sgm, '.3f')}s, solve={baseline_sr:.1f}%)")
    rpt.raw(f"  Numpy          : {'yes (fast)' if HAS_NUMPY else 'no (slow)'}")

    # ── TABLE 1: Greedy portfolio construction curve ──────────────────────────
    rpt.h2("Greedy portfolio construction — performance curve")
    rpt.raw("  At each step the param that maximally reduces SGM is added.")
    rpt.raw("  SGM = shifted geometric mean of racing_time = min(t_p, p∈portfolio) per instance.")
    rpt.raw(f"  Baseline (K=1, {args.baseline}): SGM={fmt_nan(baseline_sgm, '.3f')}s  "
            f"solve={baseline_sr:.1f}%")
    rpt.raw("")

    rows = []
    k1_sgm = greedy_curve[0][1] if greedy_curve else float("nan")
    for k, sgm, mean, sr, ns, param_added, gain in greedy_curve:
        vs_k1  = f"{(k1_sgm - sgm) / k1_sgm * 100:.1f}%" if k > 1 else "—"
        vs_base = (f"{(baseline_sgm - sgm) / baseline_sgm * 100:.1f}%"
                   if not math.isnan(baseline_sgm) else "N/A")
        step_gain = f"{gain:.1f}%" if not math.isnan(gain) else "—"
        rows.append((k, f"{sgm:.3f}", vs_k1, vs_base, f"{sr:.1f}%", ns, step_gain, param_added))

    rpt.table(
        ["K", "SGM(s)", "vs_K1", "vs_base", "Solve%", "#Solved", "ΔStep%", "Added param"],
        rows,
        fmt=[">", ">", ">", ">", ">", ">", ">", "<"])

    rpt.raw("")
    rpt.raw("  Greedy portfolio compositions:")
    running = []
    for k, sgm, mean, sr, ns, param_added, gain in greedy_curve:
        running.append(param_added)
        rpt.raw(f"    K={k}: [{', '.join(running)}]")

    # ── TABLE 2: Exhaustive optimal portfolios ────────────────────────────────
    rpt.h2(f"Exhaustive optimal portfolios (K ≤ {max_ex})")
    rpt.raw("  These are globally optimal — no combination of K params beats them.")
    rpt.raw("")

    ex_rows = []
    for k, (ex_sgm, ex_combo) in sorted(exhaustive_results.items()):
        ex_params = [params[i] for i in ex_combo]
        # Compare with greedy
        greedy_sgm = greedy_curve[k - 1][1] if k <= len(greedy_curve) else float("nan")
        greedy_params = [params[i] for i in greedy_portfolio[:k]]
        if abs(ex_sgm - greedy_sgm) < 1e-9:
            greedy_note = "✓ greedy optimal"
        else:
            diff = (greedy_sgm - ex_sgm) / ex_sgm * 100
            greedy_note = f"greedy +{diff:.2f}% worse"
        ex_rows.append((k, f"{ex_sgm:.4f}",
                        f"{greedy_sgm:.4f}", greedy_note,
                        ", ".join(ex_params)))
    rpt.table(
        ["K", "Opt SGM(s)", "Greedy SGM(s)", "Greedy quality", "Optimal params"],
        ex_rows, fmt=[">", ">", ">", "<", "<"])

    rpt.raw("")
    rpt.raw("  Optimal compositions:")
    for k, (ex_sgm, ex_combo) in sorted(exhaustive_results.items()):
        ex_params = [params[i] for i in ex_combo]
        rpt.raw(f"    K={k}: [{', '.join(ex_params)}]  SGM={ex_sgm:.4f}s")

    # ── TABLE 3: Marginal value of each param in the greedy portfolio ─────────
    rpt.h2("Marginal contribution of each param added (greedy order)")
    rpt.raw("  Shows how much each addition improves racing performance.")
    rpt.raw("")

    marginal_rows = []
    for k, sgm, mean, sr, ns, param_added, gain in greedy_curve:
        # Count instances where this param is now the portfolio winner
        cov = coverage_map(greedy_portfolio[:k])
        n_this_fastest = sum(1 for inst, (p, t, s) in cov.items() if p == param_added)
        pct_fastest = pct(n_this_fastest, n_inst)
        marginal_rows.append((
            k, param_added,
            f"{gain:.1f}%" if not math.isnan(gain) else "—",
            f"{sgm:.3f}",
            f"{sr:.1f}%", ns,
            n_this_fastest, f"{pct_fastest:.1f}%"))

    rpt.table(
        ["K", "Param added", "ΔSGM%", "New SGM(s)", "Solve%", "#Solved",
         "#Fastest inst", "%Fastest"],
        marginal_rows,
        fmt=[">", "<", ">", ">", ">", ">", ">", ">"])

    # ── TABLE 4: Per-param "uniqueness" — instances only this param solves ────
    rpt.h2("Per-param uniqueness — instances where this is the ONLY solver")
    rpt.raw("  'Unique' = all other params timed out / failed on this instance.")
    rpt.raw("  High uniqueness → param is hard to drop from any portfolio.")
    rpt.raw("")

    uniqueness_rows = []
    for param in params:
        p_col = param_idx[param]
        unique = 0
        fastest_solo = 0
        for i in range(n_inst):
            others = [j for j in range(n_params) if j != p_col]
            # Unique solver: all others failed (time >= penalty)
            if solved_mat[i][p_col] and all(not solved_mat[i][j] for j in others):
                unique += 1
            # Fastest on instance considering only successful solvers
            if solved_mat[i][p_col]:
                my_t = time_mat[i][p_col]
                others_solved = [time_mat[i][j] for j in others if solved_mat[i][j]]
                if not others_solved or my_t <= min(others_solved):
                    fastest_solo += 1
        uniqueness_rows.append((param, unique, f"{pct(unique, n_inst):.1f}%",
                                 fastest_solo, f"{pct(fastest_solo, n_inst):.1f}%"))
    uniqueness_rows.sort(key=lambda r: (-r[1], -r[3]))
    rpt.table(
        ["Param", "#Unique wins", "%Unique", "#Fastest overall", "%Fastest overall"],
        uniqueness_rows, fmt=["<", ">", ">", ">", ">"])

    # ── TABLE 5: Coverage map for greedy K=1..max shown ──────────────────────
    show_ks = [k for k in [1, 2, 3, 4, 5] if k <= len(greedy_curve)]
    rpt.h2(f"Per-param racing contributions for greedy portfolios (K={show_ks})")
    rpt.raw("  Shows how many instances each param 'wins' (is fastest) in the portfolio.")
    rpt.raw("")

    cov_header = ["Param"] + [f"K={k}" for k in show_ks]
    cov_fmt    = ["<"] + [">"] * len(show_ks)
    cov_rows   = []
    for param in params:
        row = [param]
        for k in show_ks:
            cov = coverage_map(greedy_portfolio[:k])
            n_win = sum(1 for inst, (p, t, s) in cov.items() if p == param)
            row.append(n_win if n_win > 0 else ".")
        cov_rows.append(row)
    # Sort by K=1 or K=2 wins descending
    cov_rows.sort(key=lambda r: -int(r[1]) if r[1] != "." else 0)
    rpt.table(cov_header, cov_rows, fmt=cov_fmt)

    # ── TABLE 6: Instances most improved by racing (vs best single param) ─────
    rpt.h2("Instances most improved by racing (K=2 portfolio vs best single param)")
    rpt.raw("  speedup = best_single_param_time / K=2_racing_time")
    rpt.raw("")

    if len(greedy_curve) >= 2:
        port2 = greedy_portfolio[:2]
        inst_speedups = []
        for i, inst in enumerate(instances):
            # Best individual (any param)
            best_solo = min(time_mat[i])
            racing2   = min(time_mat[i][p] for p in port2)
            speedup   = best_solo / racing2 if racing2 > 0 else float("inf")
            inst_speedups.append((inst, speedup, best_solo, racing2,
                                  solved_mat[i][port2[0]], solved_mat[i][port2[1]]))
        inst_speedups.sort(key=lambda r: -r[1])
        top_rows = []
        for inst, sp, best_t, race_t, s0, s1 in inst_speedups[:30]:
            note = "only_p2" if (not s0 and s1) else ("only_p1" if (s0 and not s1) else "both")
            top_rows.append((inst, f"{sp:.2f}×",
                              f"{best_t:.2f}s", f"{race_t:.2f}s", note))
        rpt.table(
            ["Instance", "Speedup", "BestSolo(s)", "Racing(s)", "Coverage"],
            top_rows, fmt=["<", ">", ">", ">", "<"])

    # ── TABLE 7: Instances where ALL params fail (hard instances) ────────────
    rpt.h2("Hard instances — no param solves within timelimit")
    rpt.raw("")
    hard = []
    for i, inst in enumerate(instances):
        if not any(solved_mat[i]):
            n_timeout = sum(1 for p in range(n_params) if not solved_mat[i][p])
            best_t = min(time_mat[i])
            hard.append((inst, n_timeout, f"{best_t:.0f}s"))
    if hard:
        hard.sort(key=lambda r: -r[1])
        rpt.table(["Instance", "#Params failed", "Best time (penalised)"],
                  hard, fmt=["<", ">", ">"])
        rpt.raw(f"\n  Total hard instances (not solved by any param): {len(hard)}")
    else:
        rpt.raw("  (all instances solved by at least one param)")

    # ── Summary / recommendation ──────────────────────────────────────────────
    rpt.h1("Summary / Recommendation")
    rpt.raw("")
    rpt.raw("  Racing portfolio recommendations (greedy construction):")
    rpt.raw("")
    running = []
    for k, sgm, mean, sr, ns, param_added, gain in greedy_curve[:args.max_k]:
        running.append(param_added)
        vs_base = (f"{(baseline_sgm - sgm) / baseline_sgm * 100:.1f}%"
                   if not math.isnan(baseline_sgm) else "N/A")
        rpt.raw(f"    K={k}:  SGM={sgm:.3f}s  solve={sr:.1f}%  vs_baseline={vs_base}")
        rpt.raw(f"         params: {running}")
        rpt.raw("")

    if exhaustive_results:
        rpt.raw("  Exhaustive-verified optimal (guaranteed globally optimal):")
        rpt.raw("")
        for k, (ex_sgm, ex_combo) in sorted(exhaustive_results.items()):
            ex_params = [params[i] for i in ex_combo]
            vs_base = (f"{(baseline_sgm - ex_sgm) / baseline_sgm * 100:.1f}%"
                       if not math.isnan(baseline_sgm) else "N/A")
            rpt.raw(f"    K={k}:  SGM={ex_sgm:.4f}s  vs_baseline={vs_base}")
            rpt.raw(f"         params: {ex_params}")
            rpt.raw("")

    # ── Write portfolio curve CSV ─────────────────────────────────────────────
    with open(csv_out, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["k", "construction", "sgm_s", "mean_s", "solve_pct",
                    "n_solved", "param_added", "portfolio"])
        running = []
        for k, sgm, mean, sr, ns, param_added, gain in greedy_curve:
            running.append(param_added)
            w.writerow([k, "greedy", f"{sgm:.4f}", f"{mean:.4f}",
                        f"{sr:.2f}", ns, param_added, "|".join(running)])
        for k, (ex_sgm, ex_combo) in sorted(exhaustive_results.items()):
            ex_params = [params[i] for i in ex_combo]
            ex_sgm_v, ex_mean_v, ex_sr, ex_ns = eval_fn(list(ex_combo))
            w.writerow([k, "exhaustive_opt", f"{ex_sgm:.4f}", f"{ex_mean_v:.4f}",
                        f"{ex_sr:.2f}", ex_ns, "",
                        "|".join(ex_params)])

    # ── Output ────────────────────────────────────────────────────────────────
    output = rpt.text()
    with open(out_file, "w") as f:
        f.write(output)
        f.write("\n")

    print(output)
    print(f"\n  Report saved to : {out_file}", file=sys.stderr)
    print(f"  CSV saved to    : {csv_out}", file=sys.stderr)


if __name__ == "__main__":
    main()
