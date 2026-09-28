#!/usr/bin/env python3
"""
kfold_racing_recommender.py — K-fold validation of feature-based racing portfolio
recommendation.

Strategy
--------
  The current single-param recommender predicts the best LP config per instance.
  We use that as a "group" signal to pick personalized racing portfolios:

    Training fold:
      1. Train single-param RF (12 cluster reps) → predict group for each
         training instance.
      2. Group training instances by predicted group.
      3. For each group G, exhaustively find best K=2 portfolio P2(G) and
         K=3 portfolio P3(G) over the search_params set.
      4. Label each training instance with P2(G) and P3(G) (its group's portfolio).
      5. Train pair_clf RF:   features → P2 label (string).
         Train triple_clf RF: features → P3 label (string).

    Test fold:
      1. Apply pair_clf   → predicted pair   per test instance.
      2. Apply triple_clf → predicted triple per test instance.
      3. Simulate racing: racing_time = min(time over predicted portfolio).
      4. Compute SGM speedup vs baseline.

  Additionally evaluates:
    - dual_default (baseline)
    - global K=2 (best pair on training set, same for all)
    - global K=3 (best triple on training set, same for all)
    - single-param recommender (current approach)
    - personalized K=2 (feature-based, this script)
    - personalized K=3 (feature-based, this script)

Usage
-----
  python3 lp_tune_scripts/kfold_racing_recommender.py \\
      [--avg-csv PATH]  [--features PATH]  [--partitions DIR] \\
      [--timelimit 10800]  [--search-all]  [--baseline dual_default]

Options
-------
  --avg-csv       Path to lp_avg_times_full.csv
  --features      Path to features.csv
  --partitions    Directory with fold_00.txt … fold_09.txt
  --timelimit     LP timelimit (default: 10800)
  --search-all    Search all 70 params for portfolios (not just 12 reps)
  --baseline      Baseline param (default: dual_default)
  --penalty-mult  Penalty multiplier (default: 2.0)
  --shift         SGM shift in seconds (default: 1.0)
  --out           Output report file (default: kfold_racing_recommender.txt)
"""

import argparse
import csv
import itertools
import math
import os
import sys
from collections import defaultdict

import numpy as np
from sklearn.ensemble import RandomForestClassifier

# ---------------------------------------------------------------------------
# Defaults
# ---------------------------------------------------------------------------
EXP_DIR      = os.path.expanduser("~/experiments/cbc/lp_relax_2026_05_15_noblas")
FEATURES_CSV = os.path.expanduser("~/inst/miplib/2017+spp/features.csv")
PARTITIONS   = os.path.expanduser("~/inst/miplib/2017+spp/partitions")

REPS = [
    'dual_pertv72', 'dual_pesteep_psi1_pertv61', 'dual_pesteep_scaling_off',
    'primal_idiot30', 'dual_pesteep_psineg1_pertv61', 'primal_sprint',
    'dual_pertv58', 'primal_idiot10', 'dual_pesteep_pertv58',
    'primal_idiot50', 'primal_idiot500', 'primal_idiot60',
]

DROP_PARAMS = {
    "primal_PEsteep_idiot100", "primal_idiot200_pertv61",
    "primal_idiot200_pertvm1483", "primal_idiot200_scaling_equi",
    "primal_idiot300_pertvm1483", "primal_idiot40_pertvm1483",
}

RF_SINGLE = dict(n_estimators=150, max_depth=8,  random_state=42, n_jobs=-1)
RF_RACING = dict(n_estimators=150, max_depth=8,  random_state=42, n_jobs=-1)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def shifted_geomean(values, shift=1.0):
    if not values:
        return float("nan")
    return math.exp(sum(math.log(max(v, 1e-9) + shift) for v in values) / len(values)) - shift


def speedup_str(base, val):
    if val <= 0 or math.isnan(val):
        return "N/A"
    r = base / val
    return f"{r:.3f}x {'↑' if r > 1.005 else ('↓' if r < 0.995 else '─')}"


def load_avg_times(avg_csv, penalty):
    times = defaultdict(dict)
    with open(avg_csv, newline="") as f:
        for row in csv.DictReader(f):
            inst, param = row["instance"], row["param_tag"]
            if param not in DROP_PARAMS:
                times[inst][param] = float(row["avg_wall_seconds"])
    return times


def load_features(features_path):
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


def load_folds(partitions_dir):
    folds = {}
    for fname in sorted(os.listdir(partitions_dir)):
        if not fname.startswith("fold_") or not fname.endswith(".txt"):
            continue
        idx = int(fname[5:7])
        with open(os.path.join(partitions_dir, fname)) as f:
            folds[idx] = [l.strip() for l in f if l.strip()]
    return folds


def impute_X(X, medians):
    X = X.copy()
    nan_mask = np.isnan(X)
    if nan_mask.any():
        X[nan_mask] = np.take(medians, np.where(nan_mask)[1])
    return X


def build_single_Xy(instances, times, params, penalty, inst_features, medians=None):
    """Features + best-single-param label."""
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
    X = impute_X(X, medians)
    return X, np.array(y_rows), valid, medians


def portfolio_tag(portfolio):
    """Canonical string label for a portfolio (sorted params joined with '|')."""
    return "|".join(sorted(portfolio))


def best_portfolio_exhaustive(instances, times, params, k, penalty, shift=1.0):
    """Find best k-param portfolio (min SGM) by exhaustive search. Returns (portfolio, sgm)."""
    n = len(instances)
    idx = {p: i for i, p in enumerate(params)}
    T = np.full((n, len(params)), penalty)
    for j, inst in enumerate(instances):
        for p, i in idx.items():
            T[j, i] = times[inst].get(p, penalty)

    best_sgm, best_combo = float("inf"), None
    for combo in itertools.combinations(range(len(params)), k):
        racing = T[:, list(combo)].min(axis=1)
        sgm = math.exp(np.log(racing + shift).mean()) - shift
        if sgm < best_sgm:
            best_sgm, best_combo = sgm, list(combo)
    return [params[i] for i in best_combo], best_sgm


def build_racing_Xy(instances, inst_group, group_portfolios_k, inst_features, medians):
    """Build features + portfolio-label for the racing classifier."""
    X_rows, y_rows, valid = [], [], []
    for inst in instances:
        if inst not in inst_features or inst not in inst_group:
            continue
        grp = inst_group[inst]
        if grp not in group_portfolios_k:
            continue
        X_rows.append(inst_features[inst])
        y_rows.append(portfolio_tag(group_portfolios_k[grp]))
        valid.append(inst)
    if not X_rows:
        return None, None, []
    X = np.vstack(X_rows)
    X = impute_X(X, medians)
    return X, np.array(y_rows), valid


def simulate_racing(instances, times, portfolio_per_inst, penalty):
    """Racing time per instance using per-instance portfolio (list of param names)."""
    return [
        min(times[inst].get(p, penalty) for p in portfolio_per_inst[inst])
        for inst in instances
        if inst in portfolio_per_inst
    ]


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--avg-csv",
                        default=os.path.join(EXP_DIR, "lp_avg_times_full.csv"))
    parser.add_argument("--features",    default=FEATURES_CSV)
    parser.add_argument("--partitions",  default=PARTITIONS)
    parser.add_argument("--timelimit",   type=float, default=10800.0)
    parser.add_argument("--penalty-mult",type=float, default=2.0)
    parser.add_argument("--shift",       type=float, default=1.0)
    parser.add_argument("--baseline",    default="dual_default")
    parser.add_argument("--search-all",  action="store_true",
                        help="Search all params for portfolios, not just 12 reps")
    parser.add_argument("--out",         default=None)
    args = parser.parse_args()

    penalty  = args.timelimit * args.penalty_mult
    shift    = args.shift

    # ── Load ─────────────────────────────────────────────────────────────────
    print("Loading features ...", flush=True)
    feat_names, inst_features = load_features(args.features)
    print(f"  {len(inst_features)} instances, {len(feat_names)} features")

    print("Loading avg times ...", flush=True)
    times = load_avg_times(args.avg_csv, penalty)
    all_instances = sorted(times.keys())
    all_params    = sorted({p for d in times.values() for p in d})
    for inst in all_instances:
        for p in all_params:
            times[inst].setdefault(p, penalty)
    print(f"  {len(all_instances)} instances, {len(all_params)} params")

    folds = load_folds(args.partitions)
    k_folds = len(folds)
    print(f"  {k_folds} folds from {args.partitions}")

    search_params = all_params if args.search_all else REPS
    n_pairs   = len(search_params) * (len(search_params) - 1) // 2
    n_triples = sum(1 for _ in itertools.combinations(search_params, 3))
    print(f"  Portfolio search over {len(search_params)} params "
          f"({n_pairs} pairs, {n_triples} triples)")
    print()

    out_file = args.out or os.path.join(EXP_DIR, "kfold_racing_recommender.txt")

    common = sorted(set(all_instances) & set(inst_features))

    # ── Per-fold evaluation ──────────────────────────────────────────────────
    # Accumulators for pooled SGM (each entry = one test-fold instance's time)
    acc = {
        "baseline":        [],
        "single":          [],
        "global_k2":       [],
        "global_k3":       [],
        "personalized_k2": [],
        "personalized_k3": [],
    }
    fold_results   = []   # per-fold summary rows
    inst_records   = []   # per-instance detail rows for personalization analysis

    for fold_idx in sorted(folds):
        test_names  = [n for n in folds[fold_idx] if n in inst_features and n in times]
        train_names = [n for fi, insts in folds.items()
                       if fi != fold_idx for n in insts
                       if n in inst_features and n in times]
        if not test_names or not train_names:
            continue

        print(f"Fold {fold_idx:2d}: {len(train_names)} train, "
              f"{len(test_names)} test ...", flush=True)

        # ── 1. Train single-param RF on training fold ─────────────────────
        X_tr, y_tr, tr_valid, tr_medians = build_single_Xy(
            train_names, times, REPS, penalty, inst_features)
        single_clf = RandomForestClassifier(**RF_SINGLE)
        single_clf.fit(X_tr, y_tr)

        # ── 2. Predict group for each TRAINING instance ───────────────────
        single_train_preds = single_clf.predict(X_tr)
        inst_group_train = {inst: pred
                            for inst, pred in zip(tr_valid, single_train_preds)}
        from collections import defaultdict as _dd
        groups_tr = _dd(list)
        for inst, grp in inst_group_train.items():
            groups_tr[grp].append(inst)

        print(f"         {len(groups_tr)} groups in training fold: "
              + ", ".join(f"{g}({len(v)})" for g, v in sorted(groups_tr.items())),
              flush=True)

        # ── 3. Exhaustive best K=2/K=3 per training group ────────────────
        group_best_pair   = {}  # grp → list of param names
        group_best_triple = {}
        for grp, grp_insts in groups_tr.items():
            pair, _   = best_portfolio_exhaustive(
                grp_insts, times, search_params, 2, penalty, shift)
            triple, _ = best_portfolio_exhaustive(
                grp_insts, times, search_params, 3, penalty, shift)
            group_best_pair[grp]   = pair
            group_best_triple[grp] = triple

        # ── 4. Build racing training labels ───────────────────────────────
        X_pair,   y_pair,   valid_pair   = build_racing_Xy(
            tr_valid, inst_group_train, group_best_pair,   inst_features, tr_medians)
        X_triple, y_triple, valid_triple = build_racing_Xy(
            tr_valid, inst_group_train, group_best_triple, inst_features, tr_medians)

        # ── 5. Train pair/triple classifiers ─────────────────────────────
        pair_clf   = RandomForestClassifier(**RF_RACING)
        triple_clf = RandomForestClassifier(**RF_RACING)
        pair_clf.fit(X_pair, y_pair)
        triple_clf.fit(X_triple, y_triple)

        # ── 6. Global K=2 / K=3 from training set ─────────────────────────
        global_pair,   _ = best_portfolio_exhaustive(
            tr_valid, times, search_params, 2, penalty, shift)
        global_triple, _ = best_portfolio_exhaustive(
            tr_valid, times, search_params, 3, penalty, shift)

        print(f"         global K=2: {global_pair}", flush=True)
        print(f"         global K=3: {global_triple}", flush=True)

        # ── 7. Evaluate on test fold ──────────────────────────────────────
        # Build test feature matrix
        X_te_rows, te_valid = [], []
        for name in test_names:
            X_te_rows.append(inst_features[name])
            te_valid.append(name)
        X_te = impute_X(np.vstack(X_te_rows), tr_medians)

        # Predict single param (group)
        single_preds_te = single_clf.predict(X_te)

        # Predict pair/triple portfolios
        pair_preds_te   = pair_clf.predict(X_te)
        triple_preds_te = triple_clf.predict(X_te)

        # Decode pair/triple labels back to list of params
        def decode_portfolio(label):
            return label.split("|")

        portfolio_pair_per_inst   = {inst: decode_portfolio(label)
                                     for inst, label in zip(te_valid, pair_preds_te)}
        portfolio_triple_per_inst = {inst: decode_portfolio(label)
                                     for inst, label in zip(te_valid, triple_preds_te)}

        # Collect times
        fold_base_times    = [times[inst].get(args.baseline, penalty) for inst in te_valid]
        fold_single_times  = [times[inst].get(pred, penalty)
                               for inst, pred in zip(te_valid, single_preds_te)]
        fold_gk2_times     = [min(times[inst].get(p, penalty) for p in global_pair)
                               for inst in te_valid]
        fold_gk3_times     = [min(times[inst].get(p, penalty) for p in global_triple)
                               for inst in te_valid]
        fold_pk2_times     = simulate_racing(te_valid, times, portfolio_pair_per_inst,   penalty)
        fold_pk3_times     = simulate_racing(te_valid, times, portfolio_triple_per_inst, penalty)

        sgm_base   = shifted_geomean(fold_base_times,   shift)
        sgm_single = shifted_geomean(fold_single_times, shift)
        sgm_gk2    = shifted_geomean(fold_gk2_times,    shift)
        sgm_gk3    = shifted_geomean(fold_gk3_times,    shift)
        sgm_pk2    = shifted_geomean(fold_pk2_times,    shift)
        sgm_pk3    = shifted_geomean(fold_pk3_times,    shift)

        fold_results.append({
            "fold":       fold_idx,
            "n_test":     len(te_valid),
            "sgm_base":   sgm_base,
            "sgm_single": sgm_single,
            "sgm_gk2":    sgm_gk2,
            "sgm_gk3":    sgm_gk3,
            "sgm_pk2":    sgm_pk2,
            "sgm_pk3":    sgm_pk3,
        })

        # Accumulate for pooled SGM
        acc["baseline"].extend(fold_base_times)
        acc["single"].extend(fold_single_times)
        acc["global_k2"].extend(fold_gk2_times)
        acc["global_k3"].extend(fold_gk3_times)
        acc["personalized_k2"].extend(fold_pk2_times)
        acc["personalized_k3"].extend(fold_pk3_times)

        # Per-instance records for personalization analysis
        for i, inst in enumerate(te_valid):
            inst_records.append({
                "instance":   inst,
                "base":       fold_base_times[i],
                "single":     fold_single_times[i],
                "gk2":        fold_gk2_times[i],
                "gk3":        fold_gk3_times[i],
                "pk2":        fold_pk2_times[i],
                "pk3":        fold_pk3_times[i],
                "pk2_label":  pair_preds_te[i],
                "pk3_label":  triple_preds_te[i],
                "gk2_label":  portfolio_tag(global_pair),
                "gk3_label":  portfolio_tag(global_triple),
            })

        print(f"         SGM: base={sgm_base:.3f} single={sgm_single:.3f} "
              f"gK2={sgm_gk2:.3f} gK3={sgm_gk3:.3f} "
              f"pK2={sgm_pk2:.3f} pK3={sgm_pk3:.3f}", flush=True)
        print(f"         Speedup vs base: single={speedup_str(sgm_base,sgm_single)} "
              f"gK2={speedup_str(sgm_base,sgm_gk2)} "
              f"gK3={speedup_str(sgm_base,sgm_gk3)} "
              f"pK2={speedup_str(sgm_base,sgm_pk2)} "
              f"pK3={speedup_str(sgm_base,sgm_pk3)}", flush=True)

    # ── Pooled SGM summary ────────────────────────────────────────────────────
    lines = []
    lines.append("")
    lines.append("=" * 72)
    lines.append("  K-FOLD RACING RECOMMENDER — RESULTS")
    lines.append("=" * 72)
    lines.append(f"  Experiment : {EXP_DIR}")
    lines.append(f"  Folds      : {k_folds}  |  Instances: {len(common)}")
    lines.append(f"  Search set : {len(search_params)} params"
                 f" ({'all' if args.search_all else '12 reps'})")
    lines.append(f"  Baseline   : {args.baseline}  |  Timelimit: {args.timelimit}s")
    lines.append("")

    # Per-fold table
    lines.append("  Per-fold SGM (seconds):")
    hdr = (f"  {'Fold':>4}  {'N':>4}  {'base':>7}  {'single':>7}  "
           f"{'gK2':>7}  {'gK3':>7}  {'pK2':>7}  {'pK3':>7}")
    lines.append(hdr)
    lines.append("  " + "-" * 64)
    for r in fold_results:
        lines.append(
            f"  {r['fold']:>4}  {r['n_test']:>4}  "
            f"{r['sgm_base']:>7.3f}  {r['sgm_single']:>7.3f}  "
            f"{r['sgm_gk2']:>7.3f}  {r['sgm_gk3']:>7.3f}  "
            f"{r['sgm_pk2']:>7.3f}  {r['sgm_pk3']:>7.3f}"
        )

    # Pooled summary
    lines.append("")
    lines.append("  Pooled SGM across all test folds:")
    lines.append("")
    sgm_base_pool = shifted_geomean(acc["baseline"],        shift)
    summary_methods = [
        ("dual_default (baseline)",    acc["baseline"]),
        ("single recommender",         acc["single"]),
        ("global K=2 racing",          acc["global_k2"]),
        ("global K=3 racing",          acc["global_k3"]),
        ("personalized K=2 (group)",   acc["personalized_k2"]),
        ("personalized K=3 (group)",   acc["personalized_k3"]),
    ]
    col_w = max(len(m) for m, _ in summary_methods) + 2
    hdr2 = f"  {'Method':<{col_w}}  {'Pooled SGM':>11}  {'Speedup':>10}"
    lines.append(hdr2)
    lines.append("  " + "-" * (col_w + 26))
    for name, times_list in summary_methods:
        sgm = shifted_geomean(times_list, shift)
        sp  = speedup_str(sgm_base_pool, sgm)
        lines.append(f"  {name:<{col_w}}  {sgm:>11.3f}  {sp:>10}")

    lines.append("")
    lines.append("  Interpretation:")
    lines.append("    pK2 vs gK2: does personalized portfolio beat global for K=2?")
    lines.append("    pK3 vs gK3: does personalized portfolio beat global for K=3?")
    lines.append("    pK2 vs single: does K=2 racing beat single-thread recommender?")
    lines.append("    pK3 vs single: does K=3 racing beat single-thread recommender?")
    lines.append("")
    lines.append("  If pK2 >> single: racing with personalized configs beats single method.")
    lines.append("  If pK2 ≈ gK2: grouping adds little; global portfolio is sufficient.")

    # ── Personalization gain analysis ─────────────────────────────────────────
    lines.append("")
    lines.append("=" * 72)
    lines.append("  PERSONALIZATION GAIN ANALYSIS (pK vs gK, per instance)")
    lines.append("=" * 72)

    for k_label, gk_key, pk_key, lbl_gk, lbl_pk in [
        ("K=2", "gk2", "pk2", "gk2_label", "pk2_label"),
        ("K=3", "gk3", "pk3", "gk3_label", "pk3_label"),
    ]:
        helped, hurt, same = [], [], []
        diff_helped, diff_hurt = [], []
        for r in inst_records:
            gk, pk = r[gk_key], r[pk_key]
            diff = gk - pk  # positive = personalized was faster
            if r[lbl_gk] == r[lbl_pk]:
                same.append(r["instance"])
            elif diff > 1e-6:
                helped.append((r["instance"], gk, pk, diff))
                diff_helped.append(diff)
            else:
                hurt.append((r["instance"], gk, pk, diff))
                diff_hurt.append(abs(diff))

        n_tot = len(inst_records)
        lines.append("")
        lines.append(f"  {k_label}: personalization helped {len(helped)}/{n_tot}, "
                     f"hurt {len(hurt)}/{n_tot}, "
                     f"same portfolio {len(same)}/{n_tot}")

        if diff_helped:
            avg_gain = sum(diff_helped) / len(diff_helped)
            lines.append(f"       avg time saved when helped : {avg_gain:.3f}s")
        if diff_hurt:
            avg_loss = sum(diff_hurt) / len(diff_hurt)
            lines.append(f"       avg time lost  when hurt   : {avg_loss:.3f}s")

        # Top-10 instances where personalization helped most
        helped_sorted = sorted(helped, key=lambda x: -x[3])[:10]
        if helped_sorted:
            lines.append(f"       Top instances where pK{k_label[-1]} helped:")
            lines.append(f"         {'Instance':<35}  {'gK time':>8}  {'pK time':>8}  {'saved':>7}")
            lines.append(f"         " + "-" * 63)
            for inst, gk, pk, d in helped_sorted:
                lines.append(f"         {inst:<35}  {gk:>8.3f}  {pk:>8.3f}  {d:>7.3f}s")

        # Top-10 instances where personalization hurt most
        hurt_sorted = sorted(hurt, key=lambda x: x[3])[:10]
        if hurt_sorted:
            lines.append(f"       Top instances where pK{k_label[-1]} hurt:")
            lines.append(f"         {'Instance':<35}  {'gK time':>8}  {'pK time':>8}  {'extra':>7}")
            lines.append(f"         " + "-" * 63)
            for inst, gk, pk, d in hurt_sorted:
                lines.append(f"         {inst:<35}  {gk:>8.3f}  {pk:>8.3f}  {abs(d):>7.3f}s")

    lines.append("")

    report = "\n".join(lines)
    print(report)

    with open(out_file, "w") as f:
        f.write(report + "\n")
    print(f"\nReport written to: {out_file}", file=sys.stderr)


if __name__ == "__main__":
    main()
