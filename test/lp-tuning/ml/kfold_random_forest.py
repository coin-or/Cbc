#!/usr/bin/env python3
"""
kfold_random_forest.py — K-fold evaluation using scikit-learn RandomForestClassifier.

Labels instances with the best LP param (lowest penalised avg_wall_seconds) and trains
a Random Forest to predict that label from instance features. Out-of-sample test times
are collected across all 10 folds for the pooled SGM comparison.

Key differences vs single DTC:
  - Ensemble of trees — less variance, better generalisation
  - Hard labels (best param per instance) as classification target, same as sklearn DTC
  - Uses max_features='sqrt' by default (standard for RF classification)

Usage:
    python3 kfold_random_forest.py --dir EXPERIMENT_DIR [options]

Options:
    --dir DIR              Experiment directory (must contain lp_avg_times.csv)
    --features PATH        Features CSV (default: ~/inst/miplib/2017+spp/features.csv)
    --partitions DIR       Fold files directory
                           (default: ~/inst/miplib/2017+spp/partitions)
    --n-estimators LIST    Comma-separated n_estimators values (default: 50,100,200,500)
    --max-depth LIST       Comma-separated max_depth values (default: None,4,8)
                           Use 0 for unlimited depth
    --timelimit T          LP time limit in seconds (default: 10800)
    --penalty-mult M       Penalty multiplier (default: 2.0)
    --shift S              SGM shift in seconds (default: 1.0)
    --baseline PARAM       Baseline param for comparison (default: cbc_default)
    --min-samples-leaf N   Minimum samples per leaf (default: 5)
    --out FILE             Output report file (default: kfold_random_forest.txt in exp dir)
"""

import argparse
import csv
import math
import os
import sys
from collections import defaultdict, Counter

import numpy as np
from sklearn.ensemble import RandomForestClassifier


# ── Helpers ────────────────────────────────────────────────────────────────────

def shifted_geomean(values, shift=1.0):
    if not values:
        return float("nan")
    return math.exp(sum(math.log(v + shift) for v in values) / len(values)) - shift


def fmt(v, decimals=3):
    return f"{v:.{decimals}f}" if not math.isnan(v) else "N/A"


def speedup_str(baseline_sgm, param_sgm):
    if math.isnan(baseline_sgm) or math.isnan(param_sgm) or param_sgm == 0:
        return "N/A"
    ratio = baseline_sgm / param_sgm
    return f"{ratio:.3f}x {'faster' if ratio >= 1.0 else '(slower)'}"


class Report:
    def __init__(self):
        self._lines = []

    def raw(self, s=""):
        self._lines.append(s)

    def h1(self, title):
        self._lines += ["", "=" * 74, f"  {title}", "=" * 74]

    def h2(self, title):
        self._lines += ["", f"── {title} " + "─" * max(0, 70 - len(title))]

    def table(self, headers, rows, align=None):
        if align is None:
            align = ["<"] * len(headers)
        widths = [len(h) for h in headers]
        for row in rows:
            for i, cell in enumerate(row):
                widths[i] = max(widths[i], len(str(cell)))
        sep = "  ".join("-" * w for w in widths)
        hdr = "  ".join(f"{h:{align[i]}{widths[i]}}" for i, h in enumerate(headers))
        self._lines.append(hdr)
        self._lines.append(sep)
        for row in rows:
            self._lines.append("  ".join(
                f"{str(c):{align[i]}{widths[i]}}" for i, c in enumerate(row)))

    def text(self):
        return "\n".join(self._lines)


# ── Data loading ───────────────────────────────────────────────────────────────

def load_features(features_path):
    feature_names = None
    inst_features = {}
    with open(features_path, newline="") as f:
        for i, row in enumerate(csv.reader(f)):
            if i == 0:
                feature_names = row[1:]
                continue
            name = row[0]
            vals = []
            for v in row[1:]:
                try:
                    vals.append(float(v))
                except (ValueError, TypeError):
                    vals.append(float("nan"))
            inst_features[name] = np.array(vals, dtype=np.float64)
    return feature_names, inst_features


def load_avg_times(avg_csv, penalty):
    times = defaultdict(dict)
    with open(avg_csv, newline="") as f:
        for row in csv.DictReader(f):
            times[row["instance"]][row["param_tag"]] = float(row["avg_wall_seconds"])
    instances = sorted(times)
    params    = sorted({p for d in times.values() for p in d})
    return times, instances, params


def load_folds(partitions_dir):
    folds = {}
    for fname in sorted(os.listdir(partitions_dir)):
        if not fname.startswith("fold_") or not fname.endswith(".txt"):
            continue
        idx = int(fname[5:7])
        with open(os.path.join(partitions_dir, fname)) as f:
            folds[idx] = [l.strip() for l in f if l.strip()]
    return folds


def best_param_for(inst_name, times, params, penalty):
    best_p, best_t = None, float("inf")
    for p in params:
        t = times[inst_name].get(p, penalty)
        if t < best_t:
            best_t, best_p = t, p
    return best_p


def build_Xy(inst_names, times, params, penalty, inst_features, train_medians=None):
    X_rows, y_rows, valid = [], [], []
    for name in inst_names:
        if name not in inst_features:
            continue
        X_rows.append(inst_features[name])
        y_rows.append(best_param_for(name, times, params, penalty))
        valid.append(name)

    if not X_rows:
        return None, None, [], None

    X = np.vstack(X_rows)
    if train_medians is None:
        col_medians = np.nanmedian(X, axis=0)
    else:
        col_medians = train_medians
    nan_mask = np.isnan(X)
    X[nan_mask] = np.take(col_medians, np.where(nan_mask)[1])

    return X, np.array(y_rows), valid, col_medians


# ── Main ───────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--dir",              required=True)
    parser.add_argument("--features",
                        default=os.path.expanduser("~/inst/miplib/2017+spp/features.csv"))
    parser.add_argument("--partitions",
                        default=os.path.expanduser("~/inst/miplib/2017+spp/partitions"))
    parser.add_argument("--n-estimators",     default="50,100,200,500")
    parser.add_argument("--max-depth",        default="0,4,8",
                        help="Comma-separated max_depth values; use 0 for unlimited")
    parser.add_argument("--timelimit",        type=float, default=10800.0)
    parser.add_argument("--penalty-mult",     type=float, default=2.0)
    parser.add_argument("--shift",            type=float, default=1.0)
    parser.add_argument("--baseline",         default="cbc_default")
    parser.add_argument("--min-samples-leaf", type=int,   default=5)
    parser.add_argument("--out",              default=None)
    args = parser.parse_args()

    exp_dir   = args.dir
    penalty   = args.timelimit * args.penalty_mult
    shift     = args.shift
    baseline  = args.baseline
    n_est_list = [int(v) for v in args.n_estimators.split(",")]
    max_depth_list = [None if v == "0" else int(v) for v in args.max_depth.split(",")]
    out_file  = args.out or os.path.join(exp_dir, "kfold_random_forest.txt")
    avg_csv   = os.path.join(exp_dir, "lp_avg_times.csv")

    print("Loading features ...", flush=True)
    feature_names, inst_features = load_features(args.features)
    print(f"  {len(inst_features)} instances, {len(feature_names)} features")

    print("Loading avg times ...", flush=True)
    times, all_instances, params = load_avg_times(avg_csv, penalty)
    print(f"  {len(all_instances)} instances, {len(params)} params")

    folds = load_folds(args.partitions)
    k = len(folds)
    print(f"  {k} folds")

    rep = Report()
    rep.h1("K-Fold Evaluation: scikit-learn RandomForestClassifier")
    rep.raw(f"  Experiment      : {exp_dir}")
    rep.raw(f"  Folds           : {k}")
    rep.raw(f"  n_estimators    : {n_est_list}")
    rep.raw(f"  max_depth       : {['unlimited' if d is None else d for d in max_depth_list]}")
    rep.raw(f"  Min-samples-leaf: {args.min_samples_leaf}  |  max_features : sqrt")
    rep.raw(f"  Baseline        : {baseline}  |  Penalty : {penalty:.0f}s  |  SGM shift : {shift}s")
    rep.raw("")
    rep.raw("  Label per instance = param with lowest avg_wall_seconds (penalised).")
    rep.raw("  RF learns to predict that label from instance features (Gini criterion).")
    rep.raw("  NaN features imputed with per-fold training-set column median (no leakage).")

    # Run all (n_est, max_depth) combinations
    configs = [(n, d) for d in max_depth_list for n in n_est_list]
    summary_rows = []

    for max_depth in max_depth_list:
        depth_label = "unlimited" if max_depth is None else str(max_depth)
        rep.h2(f"max_depth={depth_label} — varying n_estimators")
        rep.raw("")

        n_est_summary = []

        for n_est in n_est_list:
            config_label = f"n={n_est} depth={depth_label}"
            print(f"  {config_label} ...", flush=True)

            all_tree_times, all_base_times = [], []

            for fold_idx in sorted(folds):
                test_names  = folds[fold_idx]
                train_names = [inst for fi, insts in folds.items()
                               if fi != fold_idx for inst in insts]

                X_tr, y_tr, tr_valid, tr_medians = build_Xy(
                    train_names, times, params, penalty, inst_features)
                if X_tr is None:
                    continue

                clf = RandomForestClassifier(
                    n_estimators=n_est,
                    max_depth=max_depth,
                    min_samples_leaf=args.min_samples_leaf,
                    max_features="sqrt",
                    random_state=42,
                    n_jobs=-1,
                )
                clf.fit(X_tr, y_tr)

                X_te_rows, te_valid = [], []
                for name in test_names:
                    if name not in inst_features:
                        continue
                    X_te_rows.append(inst_features[name])
                    te_valid.append(name)

                if not X_te_rows:
                    continue

                X_te = np.vstack(X_te_rows)
                nan_mask = np.isnan(X_te)
                X_te[nan_mask] = np.take(tr_medians, np.where(nan_mask)[1])

                predictions = clf.predict(X_te)

                for name, rec_param in zip(te_valid, predictions):
                    all_tree_times.append(times[name].get(rec_param, penalty))
                    all_base_times.append(times[name].get(baseline,  penalty))

            pooled_tree = shifted_geomean(all_tree_times, shift)
            pooled_base = shifted_geomean(all_base_times, shift)
            n_est_summary.append((n_est, pooled_tree, pooled_base))
            summary_rows.append((depth_label, n_est, pooled_tree, pooled_base))
            print(f"    pooled SGM: {fmt(pooled_tree)} vs baseline {fmt(pooled_base)}"
                  f"  →  {speedup_str(pooled_base, pooled_tree)}")

        rep.table(
            ["n_estimators", "Pooled SGM (RF)", "Pooled SGM (baseline)", "Speedup"],
            [[n, fmt(st), fmt(sb), speedup_str(sb, st)] for n, st, sb in n_est_summary],
            align=["<", ">", ">", "<"],
        )

    # ── Grand summary ──────────────────────────────────────────────────────────
    rep.h2("Grand summary: all (max_depth, n_estimators) combinations")
    rep.raw("")
    rep.table(
        ["max_depth", "n_estimators", "Pooled SGM (RF)", "Pooled SGM (baseline)", "Speedup"],
        [[d, n, fmt(st), fmt(sb), speedup_str(sb, st)] for d, n, st, sb in summary_rows],
        align=["<", ">", ">", ">", "<"],
    )

    # ── Comparison with best DTC and fbps ──────────────────────────────────────
    rep.h2("Comparison with best single methods (same experiment, same folds)")
    rep.raw("")
    best_rf = min(summary_rows, key=lambda r: r[2])
    rep.raw(f"  Best RF config  : max_depth={best_rf[0]}, n_estimators={best_rf[1]}")
    rep.raw(f"  Best RF SGM     : {fmt(best_rf[2])}"
            f"  →  {speedup_str(best_rf[3], best_rf[2])}")
    rep.raw("")
    cmp_data = [
        ("cbc_default (baseline)",    13.234, "─"),
        ("Best single param",         12.923, "─"),
        ("fbps DTree depth=4",        12.316, "fbps/kfold_dtree.py"),
        ("sklearn DTC depth=4",       12.057, "kfold_sklearn_dtree.py"),
        (f"Best RF ({best_rf[0]},n={best_rf[1]})", best_rf[2], "kfold_random_forest.py"),
    ]
    rep.table(
        ["Method", "Pooled SGM", "Source"],
        [[m, fmt(s), src] for m, s, src in cmp_data],
        align=["<", ">", "<"],
    )

    # ── Feature importances from best config (full-data fit) ───────────────────
    rep.h2("Top-20 feature importances — best RF config (full-data fit)")
    rep.raw("")
    X_all, y_all, all_valid, _ = build_Xy(
        all_instances, times, params, penalty, inst_features)
    best_clf = RandomForestClassifier(
        n_estimators=best_rf[1],
        max_depth=None if best_rf[0] == "unlimited" else int(best_rf[0]),
        min_samples_leaf=args.min_samples_leaf,
        max_features="sqrt",
        random_state=42,
        n_jobs=-1,
    )
    best_clf.fit(X_all, y_all)
    importances = best_clf.feature_importances_
    top20 = sorted(zip(feature_names, importances), key=lambda x: -x[1])[:20]
    rep.table(
        ["Feature", "Importance"],
        [[name, f"{imp:.4f}"] for name, imp in top20],
        align=["<", ">"],
    )

    # ── Write & print ──────────────────────────────────────────────────────────
    report_text = rep.text()
    with open(out_file, "w") as f:
        f.write(report_text + "\n")
    print(report_text)
    print(f"\nReport written to: {out_file}", file=sys.stderr)


if __name__ == "__main__":
    main()
