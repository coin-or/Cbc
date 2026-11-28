#!/usr/bin/env python3
"""
kfold_xgboost.py — K-fold evaluation using XGBoost (XGBClassifier).

Labels each instance with the best LP param (lowest penalised avg_wall_seconds)
and trains XGBoost to predict that label from instance features.
Out-of-sample test times are pooled across all 10 folds for the SGM comparison.

Scans a grid of n_estimators × max_depth × learning_rate to find the best config.

Usage:
    python3 kfold_xgboost.py --dir EXPERIMENT_DIR [options]

Options:
    --dir DIR              Experiment directory (must contain lp_avg_times.csv)
    --features PATH        Features CSV (default: ~/inst/miplib/2017+spp/features.csv)
    --partitions DIR       Fold files directory
                           (default: ~/inst/miplib/2017+spp/partitions)
    --n-estimators LIST    Comma-separated n_estimators (default: 50,100,200,500)
    --max-depth LIST       Comma-separated max_depth (default: 3,4,6,8)
    --learning-rate LIST   Comma-separated eta values (default: 0.1,0.3)
    --subsample F          Row subsampling ratio (default: 0.8)
    --colsample F          Column subsampling per tree (default: 0.8)
    --timelimit T          LP time limit in seconds (default: 10800)
    --penalty-mult M       Penalty multiplier (default: 2.0)
    --shift S              SGM shift in seconds (default: 1.0)
    --baseline PARAM       Baseline param for comparison (default: cbc_default)
    --min-child-weight N   Min child weight (default: 5)
    --out FILE             Output report file (default: kfold_xgboost.txt in exp dir)
"""

import argparse
import csv
import math
import os
import sys
from collections import defaultdict

import numpy as np
from xgboost import XGBClassifier
from sklearn.preprocessing import LabelEncoder


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
    col_medians = np.nanmedian(X, axis=0) if train_medians is None else train_medians
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
    parser.add_argument("--max-depth",        default="3,4,6,8")
    parser.add_argument("--learning-rate",    default="0.1,0.3")
    parser.add_argument("--subsample",        type=float, default=0.8)
    parser.add_argument("--colsample",        type=float, default=0.8)
    parser.add_argument("--min-child-weight", type=int,   default=5)
    parser.add_argument("--timelimit",        type=float, default=10800.0)
    parser.add_argument("--penalty-mult",     type=float, default=2.0)
    parser.add_argument("--shift",            type=float, default=1.0)
    parser.add_argument("--baseline",         default="cbc_default")
    parser.add_argument("--out",              default=None)
    args = parser.parse_args()

    exp_dir   = args.dir
    penalty   = args.timelimit * args.penalty_mult
    shift     = args.shift
    baseline  = args.baseline
    n_est_list = [int(v) for v in args.n_estimators.split(",")]
    depth_list = [int(v) for v in args.max_depth.split(",")]
    lr_list    = [float(v) for v in args.learning_rate.split(",")]
    out_file   = args.out or os.path.join(exp_dir, "kfold_xgboost.txt")
    avg_csv    = os.path.join(exp_dir, "lp_avg_times.csv")

    print("Loading features ...", flush=True)
    feature_names, inst_features = load_features(args.features)
    print(f"  {len(inst_features)} instances, {len(feature_names)} features")

    print("Loading avg times ...", flush=True)
    times, all_instances, params = load_avg_times(avg_csv, penalty)
    print(f"  {len(all_instances)} instances, {len(params)} params")

    folds = load_folds(args.partitions)
    k = len(folds)
    print(f"  {k} folds")

    # Fit a global LabelEncoder on all possible class names so all folds use
    # consistent integer encoding (required by XGBoost).
    all_labels = sorted({best_param_for(n, times, params, penalty)
                         for n in all_instances if n in inst_features})
    le = LabelEncoder().fit(all_labels)
    n_classes = len(le.classes_)
    print(f"  {n_classes} label classes")

    rep = Report()
    rep.h1("K-Fold Evaluation: XGBoost (XGBClassifier)")
    rep.raw(f"  Experiment      : {exp_dir}")
    rep.raw(f"  Folds           : {k}")
    rep.raw(f"  n_estimators    : {n_est_list}")
    rep.raw(f"  max_depth       : {depth_list}")
    rep.raw(f"  learning_rate   : {lr_list}")
    rep.raw(f"  subsample       : {args.subsample}  |  colsample_bytree : {args.colsample}")
    rep.raw(f"  min_child_weight: {args.min_child_weight}")
    rep.raw(f"  Baseline        : {baseline}  |  Penalty : {penalty:.0f}s  |  SGM shift : {shift}s")
    rep.raw("")
    rep.raw("  Label per instance = param with lowest avg_wall_seconds (penalised).")
    rep.raw("  XGBoost predicts that label from instance features (softmax objective).")
    rep.raw("  NaN features passed as-is (XGBoost handles them natively).")

    summary_rows = []   # (lr, depth, n_est, pooled_sgm, pooled_base)

    configs = [(lr, d, n) for lr in lr_list for d in depth_list for n in n_est_list]
    total = len(configs)
    print(f"\n  Running {total} configs × {k} folds ...\n", flush=True)

    for ci, (lr, depth, n_est) in enumerate(configs, 1):
        label = f"lr={lr} depth={depth} n={n_est}"
        print(f"  [{ci:>3}/{total}] {label} ...", end=" ", flush=True)

        all_tree_times, all_base_times = [], []

        for fold_idx in sorted(folds):
            test_names  = folds[fold_idx]
            train_names = [inst for fi, insts in folds.items()
                           if fi != fold_idx for inst in insts]

            X_tr, y_tr_str, tr_valid, tr_medians = build_Xy(
                train_names, times, params, penalty, inst_features)
            if X_tr is None:
                continue
            y_tr = le.transform(y_tr_str)

            clf = XGBClassifier(
                n_estimators=n_est,
                max_depth=depth,
                learning_rate=lr,
                subsample=args.subsample,
                colsample_bytree=args.colsample,
                min_child_weight=args.min_child_weight,
                objective="multi:softmax",
                num_class=n_classes,
                eval_metric="mlogloss",
                random_state=42,
                n_jobs=-1,
                verbosity=0,
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
            # Impute NaN with training medians (no leakage)
            nan_mask = np.isnan(X_te)
            X_te[nan_mask] = np.take(tr_medians, np.where(nan_mask)[1])

            pred_ints   = clf.predict(X_te)
            predictions = le.inverse_transform(pred_ints.astype(int))

            for name, rec_param in zip(te_valid, predictions):
                all_tree_times.append(times[name].get(rec_param, penalty))
                all_base_times.append(times[name].get(baseline,  penalty))

        pooled_tree = shifted_geomean(all_tree_times, shift)
        pooled_base = shifted_geomean(all_base_times, shift)
        summary_rows.append((lr, depth, n_est, pooled_tree, pooled_base))
        print(f"SGM {fmt(pooled_tree)}  {speedup_str(pooled_base, pooled_tree)}")

    # ── Summary table sorted by SGM ────────────────────────────────────────────
    rep.h2("All configs — sorted by pooled out-of-sample SGM (best first)")
    rep.raw("")
    sorted_rows = sorted(summary_rows, key=lambda r: r[3])
    rep.table(
        ["lr", "max_depth", "n_estimators", "Pooled SGM (XGB)", "Baseline SGM", "Speedup"],
        [[lr, d, n, fmt(st), fmt(sb), speedup_str(sb, st)]
         for lr, d, n, st, sb in sorted_rows],
        align=["<", ">", ">", ">", ">", "<"],
    )

    # ── Best config breakdown by learning rate ─────────────────────────────────
    rep.h2("Best SGM per learning_rate")
    rep.raw("")
    for lr in lr_list:
        best = min((r for r in summary_rows if r[0] == lr), key=lambda r: r[3])
        rep.raw(f"  lr={lr}: best SGM {fmt(best[3])} "
                f"(depth={best[1]}, n={best[2]})  →  {speedup_str(best[4], best[3])}")

    # ── Comparison with other methods ──────────────────────────────────────────
    best_xgb = min(summary_rows, key=lambda r: r[3])
    rep.h2("Comparison with previously evaluated methods (same folds)")
    rep.raw("")
    cmp_data = [
        ("cbc_default (baseline)",          13.234),
        ("Best single param",               12.923),
        ("fbps DTree depth=4",              12.316),
        ("sklearn DTC depth=4",             12.057),
        ("Random Forest (depth=8, n=200)",  10.346),
        (f"XGBoost (lr={best_xgb[0]}, depth={best_xgb[1]}, n={best_xgb[2]})",
         best_xgb[3]),
    ]
    rep.table(
        ["Method", "Pooled SGM", "Speedup vs baseline"],
        [[m, fmt(s), speedup_str(13.234, s)] for m, s in cmp_data],
        align=["<", ">", "<"],
    )

    # ── Feature importances — best config, full-data fit ──────────────────────
    rep.h2(f"Top-20 feature importances — best XGB config (full-data fit)")
    rep.raw("")
    X_all, y_all_str, all_valid, _ = build_Xy(
        all_instances, times, params, penalty, inst_features)
    y_all = le.transform(y_all_str)
    best_clf = XGBClassifier(
        n_estimators=best_xgb[2],
        max_depth=best_xgb[1],
        learning_rate=best_xgb[0],
        subsample=args.subsample,
        colsample_bytree=args.colsample,
        min_child_weight=args.min_child_weight,
        objective="multi:softmax",
        num_class=n_classes,
        eval_metric="mlogloss",
        random_state=42,
        n_jobs=-1,
        verbosity=0,
    )
    best_clf.fit(X_all, y_all)
    importances = best_clf.feature_importances_
    top20 = sorted(zip(feature_names, importances), key=lambda x: -x[1])[:20]
    rep.table(
        ["Feature", "Importance (gain)"],
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
