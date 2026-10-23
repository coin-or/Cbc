#!/usr/bin/env python3
"""
kfold_lightgbm.py — K-fold evaluation using LightGBM (LGBMClassifier).

Labels each instance with the best LP param (lowest penalised avg_wall_seconds)
and trains LightGBM to predict that label from instance features.
Out-of-sample test times are pooled across all 10 folds for the SGM comparison.

LightGBM differences vs XGBoost:
  - Leaf-wise (best-first) tree growth instead of level-wise
  - Histogram-based splits: faster, native NaN handling
  - num_leaves controls model complexity (more direct than max_depth)

Usage:
    python3 kfold_lightgbm.py --dir EXPERIMENT_DIR [options]

Options:
    --dir DIR              Experiment directory (must contain lp_avg_times.csv)
    --features PATH        Features CSV (default: ~/inst/miplib/2017+spp/features.csv)
    --partitions DIR       Fold files directory
                           (default: ~/inst/miplib/2017+spp/partitions)
    --n-estimators LIST    Comma-separated n_estimators (default: 50,100,200)
    --num-leaves LIST      Comma-separated num_leaves (default: 15,31,63,127)
    --learning-rate LIST   Comma-separated learning rates (default: 0.05,0.1,0.3)
    --subsample F          Row subsampling ratio (default: 0.8)
    --colsample F          Column subsampling per tree (default: 0.8)
    --min-child-samples N  Min samples in a leaf (default: 5)
    --timelimit T          LP time limit in seconds (default: 10800)
    --penalty-mult M       Penalty multiplier (default: 2.0)
    --shift S              SGM shift in seconds (default: 1.0)
    --baseline PARAM       Baseline param for comparison (default: cbc_default)
    --out FILE             Output report file (default: kfold_lightgbm.txt in exp dir)
"""

import argparse
import csv
import math
import os
import sys
import warnings
from collections import defaultdict

import numpy as np
from lightgbm import LGBMClassifier

warnings.filterwarnings("ignore", category=UserWarning)


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


def build_Xy(inst_names, times, params, penalty, inst_features):
    """NaN values are kept as-is — LightGBM handles them natively."""
    X_rows, y_rows, valid = [], [], []
    for name in inst_names:
        if name not in inst_features:
            continue
        X_rows.append(inst_features[name])
        y_rows.append(best_param_for(name, times, params, penalty))
        valid.append(name)
    if not X_rows:
        return None, None, []
    return np.vstack(X_rows), np.array(y_rows), valid


# ── Main ───────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--dir",                required=True)
    parser.add_argument("--features",
                        default=os.path.expanduser("~/inst/miplib/2017+spp/features.csv"))
    parser.add_argument("--partitions",
                        default=os.path.expanduser("~/inst/miplib/2017+spp/partitions"))
    parser.add_argument("--n-estimators",       default="50,100,200")
    parser.add_argument("--num-leaves",         default="15,31,63,127")
    parser.add_argument("--learning-rate",      default="0.05,0.1,0.3")
    parser.add_argument("--subsample",          type=float, default=0.8)
    parser.add_argument("--colsample",          type=float, default=0.8)
    parser.add_argument("--min-child-samples",  type=int,   default=5)
    parser.add_argument("--timelimit",          type=float, default=10800.0)
    parser.add_argument("--penalty-mult",       type=float, default=2.0)
    parser.add_argument("--shift",              type=float, default=1.0)
    parser.add_argument("--baseline",           default="cbc_default")
    parser.add_argument("--out",                default=None)
    args = parser.parse_args()

    exp_dir        = args.dir
    penalty        = args.timelimit * args.penalty_mult
    shift          = args.shift
    baseline       = args.baseline
    n_est_list     = [int(v)   for v in args.n_estimators.split(",")]
    leaves_list    = [int(v)   for v in args.num_leaves.split(",")]
    lr_list        = [float(v) for v in args.learning_rate.split(",")]
    out_file       = args.out or os.path.join(exp_dir, "kfold_lightgbm.txt")
    avg_csv        = os.path.join(exp_dir, "lp_avg_times.csv")

    print("Loading features ...", flush=True)
    feature_names, inst_features = load_features(args.features)
    print(f"  {len(inst_features)} instances, {len(feature_names)} features")

    print("Loading avg times ...", flush=True)
    times, all_instances, params = load_avg_times(avg_csv, penalty)
    print(f"  {len(all_instances)} instances, {len(params)} params")

    folds = load_folds(args.partitions)
    k = len(folds)
    print(f"  {k} folds")

    all_labels = sorted({best_param_for(n, times, params, penalty)
                         for n in all_instances if n in inst_features})
    label_to_int = {l: i for i, l in enumerate(all_labels)}
    int_to_label = {i: l for l, i in label_to_int.items()}
    n_classes = len(all_labels)
    print(f"  {n_classes} label classes")

    configs = [(lr, nl, n) for lr in lr_list for nl in leaves_list for n in n_est_list]
    total = len(configs)

    rep = Report()
    rep.h1("K-Fold Evaluation: LightGBM (LGBMClassifier)")
    rep.raw(f"  Experiment       : {exp_dir}")
    rep.raw(f"  Folds            : {k}")
    rep.raw(f"  n_estimators     : {n_est_list}")
    rep.raw(f"  num_leaves       : {leaves_list}")
    rep.raw(f"  learning_rate    : {lr_list}")
    rep.raw(f"  subsample        : {args.subsample}  |  colsample_bytree : {args.colsample}")
    rep.raw(f"  min_child_samples: {args.min_child_samples}")
    rep.raw(f"  Baseline         : {baseline}  |  Penalty : {penalty:.0f}s  |  SGM shift : {shift}s")
    rep.raw("")
    rep.raw("  Label per instance = param with lowest avg_wall_seconds (penalised).")
    rep.raw("  LightGBM predicts that label from instance features (softmax objective).")
    rep.raw("  NaN features handled natively by LightGBM (no imputation needed).")

    summary_rows = []
    print(f"\n  Running {total} configs × {k} folds ...\n", flush=True)

    for ci, (lr, num_leaves, n_est) in enumerate(configs, 1):
        label = f"lr={lr} leaves={num_leaves} n={n_est}"
        print(f"  [{ci:>3}/{total}] {label} ...", end=" ", flush=True)

        all_tree_times, all_base_times = [], []

        for fold_idx in sorted(folds):
            test_names  = folds[fold_idx]
            train_names = [inst for fi, insts in folds.items()
                           if fi != fold_idx for inst in insts]

            X_tr, y_tr_str, tr_valid = build_Xy(
                train_names, times, params, penalty, inst_features)
            if X_tr is None:
                continue
            y_tr = np.array([label_to_int[l] for l in y_tr_str])

            clf = LGBMClassifier(
                n_estimators=n_est,
                num_leaves=num_leaves,
                learning_rate=lr,
                subsample=args.subsample,
                colsample_bytree=args.colsample,
                min_child_samples=args.min_child_samples,
                objective="multiclass",
                num_class=n_classes,
                random_state=42,
                n_jobs=-1,
                verbose=-1,
            )
            clf.fit(X_tr, y_tr)

            X_te, _, te_valid = build_Xy(
                test_names, times, params, penalty, inst_features)
            if X_te is None:
                continue

            pred_ints   = clf.predict(X_te)
            predictions = [int_to_label[int(p)] for p in pred_ints]

            for name, rec_param in zip(te_valid, predictions):
                all_tree_times.append(times[name].get(rec_param, penalty))
                all_base_times.append(times[name].get(baseline,  penalty))

        pooled_tree = shifted_geomean(all_tree_times, shift)
        pooled_base = shifted_geomean(all_base_times, shift)
        summary_rows.append((lr, num_leaves, n_est, pooled_tree, pooled_base))
        print(f"SGM {fmt(pooled_tree)}  {speedup_str(pooled_base, pooled_tree)}")

    # ── Summary sorted by SGM ──────────────────────────────────────────────────
    rep.h2("All configs — sorted by pooled out-of-sample SGM (best first)")
    rep.raw("")
    sorted_rows = sorted(summary_rows, key=lambda r: r[3])
    rep.table(
        ["lr", "num_leaves", "n_estimators", "Pooled SGM (LGBM)", "Baseline SGM", "Speedup"],
        [[lr, nl, n, fmt(st), fmt(sb), speedup_str(sb, st)]
         for lr, nl, n, st, sb in sorted_rows],
        align=["<", ">", ">", ">", ">", "<"],
    )

    # ── Best per learning rate ─────────────────────────────────────────────────
    rep.h2("Best SGM per learning_rate")
    rep.raw("")
    for lr in lr_list:
        best = min((r for r in summary_rows if r[0] == lr), key=lambda r: r[3])
        rep.raw(f"  lr={lr}: best SGM {fmt(best[3])} "
                f"(leaves={best[1]}, n={best[2]})  →  {speedup_str(best[4], best[3])}")

    # ── Comparison with all methods ────────────────────────────────────────────
    best_lgbm = min(summary_rows, key=lambda r: r[3])
    rep.h2("Comparison with all evaluated methods (same folds)")
    rep.raw("")
    cmp_data = [
        ("cbc_default (baseline)",          13.234),
        ("Best single param",               12.923),
        ("fbps DTree depth=4",              12.316),
        ("sklearn DTC depth=4",             12.057),
        ("Random Forest (depth=8, n=200)",  10.346),
        ("XGBoost (lr=0.1, depth=6, n=200)", 10.009),
        (f"LightGBM (lr={best_lgbm[0]}, leaves={best_lgbm[1]}, n={best_lgbm[2]})",
         best_lgbm[3]),
    ]
    rep.table(
        ["Method", "Pooled SGM", "Speedup vs baseline"],
        [[m, fmt(s), speedup_str(13.234, s)] for m, s in cmp_data],
        align=["<", ">", "<"],
    )

    # ── Feature importances — best config, full-data fit ──────────────────────
    rep.h2(f"Top-20 feature importances — best LGBM config (full-data fit)")
    rep.raw("")
    X_all, y_all_str, all_valid = build_Xy(
        all_instances, times, params, penalty, inst_features)
    y_all = np.array([label_to_int[l] for l in y_all_str])
    best_clf = LGBMClassifier(
        n_estimators=best_lgbm[2],
        num_leaves=best_lgbm[1],
        learning_rate=best_lgbm[0],
        subsample=args.subsample,
        colsample_bytree=args.colsample,
        min_child_samples=args.min_child_samples,
        objective="multiclass",
        num_class=n_classes,
        random_state=42,
        n_jobs=-1,
        verbose=-1,
    )
    best_clf.fit(X_all, y_all)
    importances = best_clf.feature_importances_
    top20 = sorted(zip(feature_names, importances), key=lambda x: -x[1])[:20]
    rep.table(
        ["Feature", "Importance (split count)"],
        [[name, str(imp)] for name, imp in top20],
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
