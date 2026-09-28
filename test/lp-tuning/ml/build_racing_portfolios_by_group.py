#!/usr/bin/env python3
"""
build_racing_portfolios_by_group.py — Group-based racing portfolio analysis.

Pipeline
--------
  1. Train the current 12-rep RF recommender on all instances.
  2. Predict the recommended config (group) for each instance.
  3. Group instances by their predicted recommendation.
  4. For each group, exhaustively find the best K=2 and K=3 racing portfolio
     over the 12 cluster representatives, minimising SGM(min-time-in-portfolio).
  5. Print a summary table: group → best pair + triple + per-group SGM gains.
  6. Compare vs global K=2 / K=3 portfolios and single-param recommender.
  7. Write a JSON with the group→portfolio mapping (input for the full k-fold).

The group→portfolio mapping is the core artifact:
  - At inference time: predict group (same RF as current recommender)
    → look up hardcoded portfolio for that group → use it for racing.
  - No new model is needed; only the lookup table changes.

Usage
-----
  python3 lp_tune_scripts/build_racing_portfolios_by_group.py \\
      [--avg-csv PATH]   [--features PATH]  [--out-json PATH]
      [--timelimit 10800] [--search-all]    [--baseline dual_default]

Options
-------
  --avg-csv     Path to lp_avg_times_full.csv  (default: auto-detect from EXP_DIR)
  --features    Path to features.csv           (default: ~/inst/miplib/2017+spp/features.csv)
  --out-json    Output JSON file for group→portfolio map  (default: racing_portfolios_by_group.json)
  --timelimit   LP timelimit in seconds  (default: 10800)
  --search-all  Search all 70 params instead of just the 12 cluster reps
  --baseline    Baseline param for speedup comparison  (default: dual_default)
  --penalty-mult   Penalty multiplier over timelimit  (default: 2.0)
  --shift       SGM shift  (default: 1.0)
"""

import argparse
import csv
import itertools
import json
import math
import os
import sys
from collections import defaultdict

import numpy as np
from sklearn.ensemble import RandomForestClassifier

# ---------------------------------------------------------------------------
# Defaults — same as build_lp_param_scorer.py
# ---------------------------------------------------------------------------
EXP_DIR      = os.path.expanduser("~/experiments/cbc/lp_relax_2026_05_15_noblas")
FEATURES_CSV = os.path.expanduser("~/inst/miplib/2017+spp/features.csv")

REPS = [
    'dual_pertv72', 'dual_pesteep_psi1_pertv61', 'dual_pesteep_scaling_off',
    'primal_idiot30', 'dual_pesteep_psineg1_pertv61', 'primal_sprint',
    'dual_pertv58', 'primal_idiot10', 'dual_pesteep_pertv58',
    'primal_idiot50', 'primal_idiot500', 'primal_idiot60',
]

RF_PARAMS = dict(n_estimators=150, max_depth=8, random_state=42, n_jobs=-1)

DROP_PARAMS = {
    "primal_PEsteep_idiot100", "primal_idiot200_pertv61",
    "primal_idiot200_pertvm1483", "primal_idiot200_scaling_equi",
    "primal_idiot300_pertvm1483", "primal_idiot40_pertvm1483",
}


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def shifted_geomean(values, shift=1.0):
    if not values:
        return float("nan")
    return math.exp(sum(math.log(max(v, 1e-9) + shift) for v in values) / len(values)) - shift


def load_avg_times(avg_csv, penalty):
    """Returns times[inst][param] = avg_wall_seconds, filling missing with penalty."""
    times = defaultdict(dict)
    with open(avg_csv, newline="") as f:
        for row in csv.DictReader(f):
            inst, param = row["instance"], row["param_tag"]
            if param not in DROP_PARAMS:
                times[inst][param] = float(row["avg_wall_seconds"])
    return times


def load_features(features_path):
    """Returns feature_names (list), inst_features dict[inst] = np.array."""
    feature_names = None
    inst_features = {}
    with open(features_path, newline="") as f:
        for i, row in enumerate(csv.reader(f)):
            if i == 0:
                feature_names = row[1:]
                continue
            vals = []
            for v in row[1:]:
                try:
                    vals.append(float(v))
                except (ValueError, TypeError):
                    vals.append(float("nan"))
            inst_features[row[0]] = np.array(vals, dtype=np.float64)
    return feature_names, inst_features


def build_Xy(instances, times, params, penalty, inst_features, medians=None):
    """Build feature matrix X and best-param labels y for a set of instances."""
    X_rows, y_rows, valid = [], [], []
    for name in instances:
        if name not in inst_features:
            continue
        row_times = [times[name].get(p, penalty) for p in params]
        best_p = params[int(np.argmin(row_times))]
        X_rows.append(inst_features[name])
        y_rows.append(best_p)
        valid.append(name)
    if not X_rows:
        return None, None, [], None
    X = np.vstack(X_rows)
    if medians is None:
        medians = np.nanmedian(X, axis=0)
    nan_mask = np.isnan(X)
    X[nan_mask] = np.take(medians, np.where(nan_mask)[1])
    return X, np.array(y_rows), valid, medians


def best_portfolio_exhaustive(instances, times, params, k, penalty, shift=1.0):
    """Exhaustive search: best k-param portfolio minimising SGM(min_time)."""
    n = len(instances)
    idx = {p: i for i, p in enumerate(params)}
    T = np.full((n, len(params)), penalty)
    for j, inst in enumerate(instances):
        for p, i in idx.items():
            T[j, i] = times[inst].get(p, penalty)

    best_sgm, best_combo = float("inf"), None
    for combo in itertools.combinations(range(len(params)), k):
        combo_l = list(combo)
        racing = T[:, combo_l].min(axis=1)
        sgm = math.exp(np.log(racing + shift).mean()) - shift
        if sgm < best_sgm:
            best_sgm = sgm
            best_combo = combo_l
    return [params[i] for i in best_combo], best_sgm


def simulate_racing(instances, times, portfolio, penalty):
    """Simulated racing time = min time over portfolio per instance."""
    return [min(times[inst].get(p, penalty) for p in portfolio)
            for inst in instances]


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--avg-csv",
                        default=os.path.join(EXP_DIR, "lp_avg_times_full.csv"))
    parser.add_argument("--features",    default=FEATURES_CSV)
    parser.add_argument("--out-json",    default="racing_portfolios_by_group.json")
    parser.add_argument("--timelimit",   type=float, default=10800.0)
    parser.add_argument("--penalty-mult",type=float, default=2.0)
    parser.add_argument("--shift",       type=float, default=1.0)
    parser.add_argument("--baseline",    default="dual_default")
    parser.add_argument("--search-all",  action="store_true",
                        help="Search all 70 params for portfolios (not just 12 reps)")
    args = parser.parse_args()

    penalty = args.timelimit * args.penalty_mult
    shift   = args.shift

    # ── Load ─────────────────────────────────────────────────────────────────
    print("Loading features ...", flush=True)
    feat_names, inst_features = load_features(args.features)
    print(f"  {len(inst_features)} instances, {len(feat_names)} features")

    print("Loading avg times ...", flush=True)
    times = load_avg_times(args.avg_csv, penalty)
    all_instances = sorted(times.keys())
    all_params    = sorted({p for d in times.values() for p in d})
    # Fill missing with penalty
    for inst in all_instances:
        for p in all_params:
            times[inst].setdefault(p, penalty)
    print(f"  {len(all_instances)} instances, {len(all_params)} params")

    search_params = all_params if args.search_all else REPS
    print(f"  Portfolio search over: {len(search_params)} configs"
          f"  (C({len(search_params)},2)={len(search_params)*(len(search_params)-1)//2} pairs,"
          f"  C({len(search_params)},3)={len(list(itertools.combinations(search_params,3)))} triples)")

    common = sorted(set(all_instances) & set(inst_features.keys()))
    print(f"  Instances with features: {len(common)}")
    print()

    # ── Step 1: train RF recommender on all instances ────────────────────────
    print("Training single-param RF recommender on all instances ...", flush=True)
    X_all, y_all, valid_all, medians = build_Xy(
        common, times, REPS, penalty, inst_features)
    clf = RandomForestClassifier(**RF_PARAMS)
    clf.fit(X_all, y_all)

    # Predict recommendation (group) for each instance
    preds = clf.predict(X_all)
    inst_group = {inst: pred for inst, pred in zip(valid_all, preds)}
    print(f"  Trained on {len(valid_all)} instances")

    # ── Step 2: group instances by recommendation ────────────────────────────
    groups = defaultdict(list)
    for inst, grp in inst_group.items():
        groups[grp].append(inst)
    print(f"\nGroups ({len(groups)} unique recommendations):")
    for grp in sorted(groups, key=lambda g: -len(groups[g])):
        print(f"  {grp:<40}  {len(groups[grp]):3d} instances")
    print()

    # ── Step 3: exhaustive best portfolio per group ───────────────────────────
    print("Finding best K=2 and K=3 portfolios per group (exhaustive) ...", flush=True)
    group_portfolios = {}
    for grp in sorted(groups):
        grp_insts = groups[grp]
        pair,  sgm_pair  = best_portfolio_exhaustive(
            grp_insts, times, search_params, 2, penalty, shift)
        triple, sgm_triple = best_portfolio_exhaustive(
            grp_insts, times, search_params, 3, penalty, shift)
        # Baseline SGM for this group
        base_times = [times[inst].get(args.baseline, penalty) for inst in grp_insts]
        sgm_base   = shifted_geomean(base_times, shift)
        # Single-param recommender SGM for this group
        single_times = [times[inst].get(grp, penalty) for inst in grp_insts]
        sgm_single   = shifted_geomean(single_times, shift)

        group_portfolios[grp] = {
            "n_instances": len(grp_insts),
            "pair":   pair,
            "triple": triple,
            "sgm_baseline": sgm_base,
            "sgm_single":   sgm_single,
            "sgm_pair":     sgm_pair,
            "sgm_triple":   sgm_triple,
        }
        sp_single = sgm_base / sgm_single if sgm_single > 0 else float("nan")
        sp_pair   = sgm_base / sgm_pair   if sgm_pair   > 0 else float("nan")
        sp_triple = sgm_base / sgm_triple if sgm_triple > 0 else float("nan")
        print(f"\n  Group: {grp}  ({len(grp_insts)} instances)")
        print(f"    Baseline  SGM: {sgm_base:.3f}s")
        print(f"    Single    SGM: {sgm_single:.3f}s  ({sp_single:.3f}x)")
        print(f"    Best K=2  SGM: {sgm_pair:.3f}s  ({sp_pair:.3f}x)  → {pair}")
        print(f"    Best K=3  SGM: {sgm_triple:.3f}s  ({sp_triple:.3f}x)  → {triple}")

    # ── Step 4: overall personalized portfolio evaluation ────────────────────
    print("\n\n" + "="*70)
    print("Overall performance (all instances, personalised racing vs global)")
    print("="*70)

    personalized_pair_times   = []
    personalized_triple_times = []
    single_times_all          = []
    baseline_times_all        = []

    for inst in valid_all:
        grp = inst_group[inst]
        pair   = group_portfolios[grp]["pair"]
        triple = group_portfolios[grp]["triple"]
        personalized_pair_times.append(min(times[inst].get(p, penalty) for p in pair))
        personalized_triple_times.append(min(times[inst].get(p, penalty) for p in triple))
        single_times_all.append(times[inst].get(grp, penalty))
        baseline_times_all.append(times[inst].get(args.baseline, penalty))

    # Global best K=2 and K=3 portfolios (in-sample, for comparison)
    print("\nFinding global K=2 and K=3 (in-sample) ...", flush=True)
    global_pair,   sgm_gp = best_portfolio_exhaustive(
        valid_all, times, search_params, 2, penalty, shift)
    global_triple, sgm_gt = best_portfolio_exhaustive(
        valid_all, times, search_params, 3, penalty, shift)
    global_pair_times   = simulate_racing(valid_all, times, global_pair,   penalty)
    global_triple_times = simulate_racing(valid_all, times, global_triple, penalty)

    sgm_base    = shifted_geomean(baseline_times_all,        shift)
    sgm_single  = shifted_geomean(single_times_all,          shift)
    sgm_gp_eval = shifted_geomean(global_pair_times,         shift)
    sgm_gt_eval = shifted_geomean(global_triple_times,       shift)
    sgm_pp      = shifted_geomean(personalized_pair_times,   shift)
    sgm_pt      = shifted_geomean(personalized_triple_times, shift)

    def speedup(base, val):
        return base / val if val > 0 else float("nan")

    rows = [
        ("dual_default (baseline)",     sgm_base,    1.000),
        ("single recommender",          sgm_single,  speedup(sgm_base, sgm_single)),
        ("global K=2 racing",           sgm_gp_eval, speedup(sgm_base, sgm_gp_eval)),
        ("global K=3 racing",           sgm_gt_eval, speedup(sgm_base, sgm_gt_eval)),
        ("personalized K=2 (group)",    sgm_pp,      speedup(sgm_base, sgm_pp)),
        ("personalized K=3 (group)",    sgm_pt,      speedup(sgm_base, sgm_pt)),
    ]
    header = f"  {'Method':<35}  {'SGM (s)':>8}  {'Speedup':>8}"
    print("\n" + header)
    print("  " + "-" * 55)
    for name, sgm, sp in rows:
        print(f"  {name:<35}  {sgm:8.3f}  {sp:8.3f}x")

    print(f"\n  Note: global K=2={global_pair}  K=3={global_triple}")
    print(f"  Note: personalized uses group-specific best portfolio (in-sample).")
    print(f"  K-fold validation in kfold_racing_recommender.py gives out-of-sample estimate.")

    # ── Step 5: write JSON output ─────────────────────────────────────────────
    out_data = {
        "reps":             REPS,
        "search_params":    search_params,
        "baseline":         args.baseline,
        "global_pair":      global_pair,
        "global_triple":    global_triple,
        "groups":           group_portfolios,
    }
    with open(args.out_json, "w") as f:
        json.dump(out_data, f, indent=2)
    print(f"\nGroup→portfolio mapping written to: {args.out_json}")


if __name__ == "__main__":
    main()
