#!/usr/bin/env python3
"""
kfold_racing_portfolio.py — K-fold validation for racing LP portfolios.

For a given portfolio size K_race (e.g. 2 or 3), evaluates how much
improvement we get by finding the optimal K_race-param racing portfolio on
the TRAINING folds and measuring the SGM speedup on the HELD-OUT test fold.

This answers: "If we hard-code the best K=2 (or K=3) racing config, what
out-of-sample speedup can we expect vs dual_default?"

Method
------
  For each fold:
    1. Find the exhaustive-optimal K_race portfolio on training instances.
    2. Simulate racing on test instances: racing_time = min(t_p for p in portfolio).
    3. Measure SGM(racing_time) / SGM(dual_default_time) on test instances.
  Report mean ± std speedup across folds.

Usage
-----
  python3 lp_tune_scripts/kfold_racing_portfolio.py \\
      --avg-csv ~/experiments/cbc/lp_relax_2026_05_15_noblas/lp_avg_times_full.csv \\
      --timelimit 10800 \\
      --k-race 2 3 \\
      --n-folds 10 \\
      --baseline dual_default
"""

import argparse
import csv
import itertools
import math
import random
import sys

import numpy as np


def shifted_geomean(values, shift=1.0):
    if len(values) == 0:
        return float("nan")
    return math.exp(sum(math.log(max(v, 1e-9) + shift) for v in values) / len(values)) - shift


def racing_time(inst_times, portfolio):
    """Effective time = min over portfolio params; penalty if all missing."""
    times = [inst_times[p] for p in portfolio if p in inst_times]
    return min(times) if times else float("inf")


def best_portfolio_exhaustive(instances, times, params, k, penalty):
    """Find optimal k-param portfolio (minimise SGM) by exhaustive search."""
    best_sgm = float("inf")
    best_combo = None
    n = len(instances)
    shift = 1.0

    # Pre-build numpy matrix for speed: shape (n_instances, n_params)
    idx = {p: i for i, p in enumerate(params)}
    T = np.full((n, len(params)), penalty)
    for j, inst in enumerate(instances):
        for p, i in idx.items():
            if p in times[inst]:
                T[j, i] = times[inst][p]

    log_T = np.log(T + shift)

    for combo in itertools.combinations(range(len(params)), k):
        combo = list(combo)
        racing = T[:, combo].min(axis=1)
        sgm = math.exp(np.log(racing + shift).mean()) - shift
        if sgm < best_sgm:
            best_sgm = sgm
            best_combo = combo

    return [params[i] for i in best_combo], best_sgm


def eval_portfolio(instances, times, portfolio, baseline, penalty, shift=1.0):
    racing_times  = [racing_time(times[inst], portfolio) for inst in instances]
    baseline_times = [times[inst].get(baseline, penalty) for inst in instances]
    sgm_racing   = shifted_geomean(racing_times,   shift)
    sgm_baseline = shifted_geomean(baseline_times, shift)
    speedup = sgm_baseline / sgm_racing if sgm_racing > 0 else float("nan")
    return sgm_racing, sgm_baseline, speedup


def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--avg-csv", required=True,
                        help="Path to lp_avg_times_full.csv")
    parser.add_argument("--timelimit", type=float, default=10800.0)
    parser.add_argument("--penalty-mult", type=float, default=2.0)
    parser.add_argument("--shift", type=float, default=1.0)
    parser.add_argument("--k-race", type=int, nargs="+", default=[2, 3],
                        help="Portfolio sizes to evaluate (default: 2 3)")
    parser.add_argument("--n-folds", type=int, default=10)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--baseline", default="dual_default")
    parser.add_argument("--exclude-prefix", default="guess_")
    args = parser.parse_args()

    penalty = args.timelimit * args.penalty_mult
    shift   = args.shift

    # ── Load data ─────────────────────────────────────────────────────────────
    times = {}   # inst -> {param -> avg_wall_seconds}
    params_set = set()
    with open(args.avg_csv, newline="") as f:
        for row in csv.DictReader(f):
            inst  = row["instance"]
            param = row["param_tag"]
            if args.exclude_prefix and param.startswith(args.exclude_prefix):
                continue
            t = float(row["avg_wall_seconds"])
            times.setdefault(inst, {})[param] = t
            params_set.add(param)

    instances = sorted(times.keys())
    params    = sorted(params_set)
    n         = len(instances)
    print(f"Instances: {n},  Params: {len(params)},  Folds: {args.n_folds}")
    print(f"Baseline: {args.baseline}  Timelimit: {args.timelimit}s  Penalty: {penalty}s")
    print()

    # Fill missing entries with penalty
    for inst in instances:
        for p in params:
            times[inst].setdefault(p, penalty)

    # ── K-fold split ──────────────────────────────────────────────────────────
    rng = random.Random(args.seed)
    shuffled = instances[:]
    rng.shuffle(shuffled)
    fold_size = n // args.n_folds
    folds = []
    for i in range(args.n_folds):
        start = i * fold_size
        end   = start + fold_size if i < args.n_folds - 1 else n
        folds.append(shuffled[start:end])

    # ── Evaluate each K_race ──────────────────────────────────────────────────
    for k_race in args.k_race:
        print(f"{'='*60}")
        print(f"  Racing portfolio size K={k_race}")
        print(f"{'='*60}")

        fold_speedups   = []
        fold_sgm_racing = []
        fold_sgm_base   = []
        fold_portfolios = []

        for fold_idx in range(args.n_folds):
            test_insts  = folds[fold_idx]
            train_insts = [inst for i, fold in enumerate(folds)
                           if i != fold_idx for inst in fold]

            # Find best portfolio on training set
            portfolio, train_sgm = best_portfolio_exhaustive(
                train_insts, times, params, k_race, penalty)

            # Evaluate on test set
            sgm_r, sgm_b, speedup = eval_portfolio(
                test_insts, times, portfolio, args.baseline, penalty, shift)

            fold_speedups.append(speedup)
            fold_sgm_racing.append(sgm_r)
            fold_sgm_base.append(sgm_b)
            fold_portfolios.append(portfolio)

            print(f"  Fold {fold_idx+1:2d}: train_SGM={train_sgm:.3f}s  "
                  f"test_SGM={sgm_r:.3f}s  baseline={sgm_b:.3f}s  "
                  f"speedup={speedup:.3f}x")
            print(f"           portfolio: {portfolio}")

        mean_sp  = float(np.mean(fold_speedups))
        std_sp   = float(np.std(fold_speedups))
        min_sp   = float(np.min(fold_speedups))
        max_sp   = float(np.max(fold_speedups))
        mean_sgm = float(np.mean(fold_sgm_racing))
        mean_base= float(np.mean(fold_sgm_base))

        print()
        print(f"  K={k_race} summary:")
        print(f"    Mean speedup : {mean_sp:.3f}x  (std={std_sp:.3f}, "
              f"min={min_sp:.3f}, max={max_sp:.3f})")
        print(f"    Mean test SGM: {mean_sgm:.3f}s  (baseline={mean_base:.3f}s)")
        print()

        # Portfolio frequency — which params appear most often
        from collections import Counter
        param_counts = Counter(p for port in fold_portfolios for p in port)
        print(f"  Portfolio stability (how often each param appears across folds):")
        for param, count in param_counts.most_common():
            bar = "█" * count + "░" * (args.n_folds - count)
            print(f"    {count:2d}/{args.n_folds}  {bar}  {param}")
        print()

    # ── Full-data optimal for reference ──────────────────────────────────────
    print(f"{'='*60}")
    print("  Full-data optimal portfolios (in-sample reference)")
    print(f"{'='*60}")
    for k_race in args.k_race:
        portfolio, sgm = best_portfolio_exhaustive(
            instances, times, params, k_race, penalty)
        sgm_r, sgm_b, speedup = eval_portfolio(
            instances, times, portfolio, args.baseline, penalty, shift)
        print(f"  K={k_race}: SGM={sgm:.3f}s  speedup={speedup:.3f}x  {portfolio}")
    print()


if __name__ == "__main__":
    main()
