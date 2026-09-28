#!/usr/bin/env python3
"""
kfold_multioutput_regression.py — K-fold evaluation using multi-output regression.

Instead of predicting the single best param label (classification), this script
trains models to predict the solve time for ALL params simultaneously.  At inference
the recommended param is argmin of the predicted time vector.

This formulation:
  - Uses all time information (no hard label; near-ties contribute less gradient)
  - Penalises catastrophic wrong choices more than near-tie wrong choices
  - Target = log1p(avg_wall_seconds) per param, transformed back for SGM

Methods evaluated:
  - RandomForestRegressor     — natively multi-output (one forest, all outputs)
  - MultiOutputRegressor(XGBRegressor)    — 24 independent XGB boosters
  - MultiOutputRegressor(LGBMRegressor)   — 24 independent LGBM boosters

Usage:
    python3 kfold_multioutput_regression.py --dir EXPERIMENT_DIR [options]

Options:
    --dir DIR           Experiment directory (must contain lp_avg_times.csv)
    --features PATH     Features CSV (default: ~/inst/miplib/2017+spp/features.csv)
    --partitions DIR    Fold files directory
                        (default: ~/inst/miplib/2017+spp/partitions)
    --methods LIST      Comma-separated: rf,xgb,lgbm  (default: rf,xgb,lgbm)
    --timelimit T       LP time limit in seconds (default: 10800)
    --penalty-mult M    Penalty multiplier (default: 2.0)
    --shift S           SGM shift in seconds (default: 1.0)
    --baseline PARAM    Baseline param for comparison (default: cbc_default)
    --out FILE          Output report file (default: kfold_multioutput_regression.txt)
"""

import argparse
import csv
import math
import os
import sys
import warnings
from collections import defaultdict

import numpy as np
from sklearn.ensemble import RandomForestRegressor
from sklearn.multioutput import MultiOutputRegressor

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


def build_XY(inst_names, times, params, penalty, inst_features, train_medians=None):
    """
    Build feature matrix X and target matrix Y.
    Y[i, j] = log1p(avg_wall_seconds for instance i, param j).
    Returns X (n×d), Y (n×p), valid_names, col_medians_used.
    """
    X_rows, Y_rows, valid = [], [], []
    for name in inst_names:
        if name not in inst_features:
            continue
        X_rows.append(inst_features[name])
        y_row = [math.log1p(times[name].get(p, penalty)) for p in params]
        Y_rows.append(y_row)
        valid.append(name)

    if not X_rows:
        return None, None, [], None

    X = np.vstack(X_rows)
    col_medians = np.nanmedian(X, axis=0) if train_medians is None else train_medians
    nan_mask = np.isnan(X)
    X[nan_mask] = np.take(col_medians, np.where(nan_mask)[1])

    return X, np.array(Y_rows, dtype=np.float64), valid, col_medians


def evaluate_predictions(Y_pred, te_valid, params, times, penalty, baseline):
    """Given predicted log1p-time matrix, return actual tree times and baseline times."""
    tree_times, base_times = [], []
    for i, name in enumerate(te_valid):
        best_param_idx = int(np.argmin(Y_pred[i]))
        rec_param = params[best_param_idx]
        tree_times.append(times[name].get(rec_param, penalty))
        base_times.append(times[name].get(baseline, penalty))
    return tree_times, base_times


# ── Grid configs ───────────────────────────────────────────────────────────────

def rf_configs():
    """(label, kwargs) pairs for RandomForestRegressor."""
    return [
        ("n=50  depth=8",   dict(n_estimators=50,  max_depth=8,    min_samples_leaf=5)),
        ("n=100 depth=8",   dict(n_estimators=100, max_depth=8,    min_samples_leaf=5)),
        ("n=200 depth=8",   dict(n_estimators=200, max_depth=8,    min_samples_leaf=5)),
        ("n=50  depth=∞",   dict(n_estimators=50,  max_depth=None, min_samples_leaf=5)),
        ("n=100 depth=∞",   dict(n_estimators=100, max_depth=None, min_samples_leaf=5)),
        ("n=200 depth=∞",   dict(n_estimators=200, max_depth=None, min_samples_leaf=5)),
    ]


def xgb_configs():
    """(label, kwargs) pairs for XGBRegressor."""
    from xgboost import XGBRegressor
    cfgs = []
    for lr in [0.1, 0.3]:
        for depth in [4, 6]:
            for n in [100, 200]:
                label = f"lr={lr} d={depth} n={n}"
                cfgs.append((label, dict(
                    n_estimators=n, max_depth=depth, learning_rate=lr,
                    subsample=0.8, colsample_bytree=0.8, min_child_weight=5,
                    objective="reg:squarederror", random_state=42,
                    n_jobs=-1, verbosity=0,
                )))
    return cfgs


def lgbm_configs():
    """(label, kwargs) pairs for LGBMRegressor."""
    from lightgbm import LGBMRegressor
    cfgs = []
    for lr in [0.05, 0.1]:
        for leaves in [31, 63]:
            for n in [100, 200]:
                label = f"lr={lr} leaves={leaves} n={n}"
                cfgs.append((label, dict(
                    n_estimators=n, num_leaves=leaves, learning_rate=lr,
                    subsample=0.8, colsample_bytree=0.8, min_child_samples=5,
                    objective="regression", random_state=42,
                    n_jobs=-1, verbose=-1,
                )))
    return cfgs


# ── Run one method grid ────────────────────────────────────────────────────────

def run_method(method_name, configs, make_estimator_fn,
               folds, times, params, penalty, inst_features, baseline, shift):
    """
    Run all configs for one method via k-fold CV.
    Returns list of (label, pooled_sgm, pooled_base).
    """
    results = []
    total = len(configs)
    print(f"\n  {'─'*60}", flush=True)
    print(f"  {method_name}  ({total} configs × {len(folds)} folds)", flush=True)

    for ci, (label, kwargs) in enumerate(configs, 1):
        print(f"  [{ci:>2}/{total}] {label} ...", end=" ", flush=True)

        all_tree_times, all_base_times = [], []

        for fold_idx in sorted(folds):
            test_names  = folds[fold_idx]
            train_names = [inst for fi, insts in folds.items()
                           if fi != fold_idx for inst in insts]

            X_tr, Y_tr, tr_valid, tr_medians = build_XY(
                train_names, times, params, penalty, inst_features)
            if X_tr is None:
                continue

            clf = make_estimator_fn(kwargs)
            clf.fit(X_tr, Y_tr)

            X_te, _, te_valid, _ = build_XY(
                test_names, times, params, penalty, inst_features, tr_medians)
            if X_te is None:
                continue

            Y_pred = clf.predict(X_te)
            tt, bt = evaluate_predictions(
                Y_pred, te_valid, params, times, penalty, baseline)
            all_tree_times.extend(tt)
            all_base_times.extend(bt)

        pooled_tree = shifted_geomean(all_tree_times, shift)
        pooled_base = shifted_geomean(all_base_times, shift)
        results.append((label, pooled_tree, pooled_base))
        print(f"SGM {fmt(pooled_tree)}  {speedup_str(pooled_base, pooled_tree)}")

    return results


# ── Main ───────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--dir",          required=True)
    parser.add_argument("--features",
                        default=os.path.expanduser("~/inst/miplib/2017+spp/features.csv"))
    parser.add_argument("--partitions",
                        default=os.path.expanduser("~/inst/miplib/2017+spp/partitions"))
    parser.add_argument("--methods",      default="rf,xgb,lgbm")
    parser.add_argument("--timelimit",    type=float, default=10800.0)
    parser.add_argument("--penalty-mult", type=float, default=2.0)
    parser.add_argument("--shift",        type=float, default=1.0)
    parser.add_argument("--baseline",     default="cbc_default")
    parser.add_argument("--out",          default=None)
    args = parser.parse_args()

    exp_dir  = args.dir
    penalty  = args.timelimit * args.penalty_mult
    shift    = args.shift
    baseline = args.baseline
    methods  = [m.strip() for m in args.methods.split(",")]
    out_file = args.out or os.path.join(exp_dir, "kfold_multioutput_regression.txt")
    avg_csv  = os.path.join(exp_dir, "lp_avg_times.csv")

    print("Loading features ...", flush=True)
    feature_names, inst_features = load_features(args.features)
    print(f"  {len(inst_features)} instances, {len(feature_names)} features")

    print("Loading avg times ...", flush=True)
    times, all_instances, params = load_avg_times(avg_csv, penalty)
    print(f"  {len(all_instances)} instances, {len(params)} params")

    folds = load_folds(args.partitions)
    k = len(folds)
    print(f"  {k} folds, {len(params)} regression targets (one per param)")

    rep = Report()
    rep.h1("K-Fold Evaluation: Multi-Output Regression")
    rep.raw(f"  Experiment  : {exp_dir}")
    rep.raw(f"  Folds       : {k}  |  Targets : {len(params)} params")
    rep.raw(f"  Target      : log1p(avg_wall_seconds) per param — argmin at inference")
    rep.raw(f"  Baseline    : {baseline}  |  Penalty : {penalty:.0f}s  |  SGM shift : {shift}s")
    rep.raw("")
    rep.raw("  For each test instance: predict time vector → pick param with min predicted time.")
    rep.raw("  Evaluated on actual (not predicted) times for SGM comparison.")

    all_results = {}   # method_name -> list of (label, sgm, base_sgm)

    # ── Random Forest ──────────────────────────────────────────────────────────
    if "rf" in methods:
        def make_rf(kwargs):
            return RandomForestRegressor(random_state=42, n_jobs=-1, **kwargs)

        all_results["RF"] = run_method(
            "Random Forest (multi-output, native)",
            rf_configs(), make_rf,
            folds, times, params, penalty, inst_features, baseline, shift)

    # ── XGBoost ────────────────────────────────────────────────────────────────
    if "xgb" in methods:
        from xgboost import XGBRegressor

        def make_xgb(kwargs):
            return MultiOutputRegressor(XGBRegressor(**kwargs), n_jobs=1)

        all_results["XGBoost"] = run_method(
            "XGBoost (MultiOutputRegressor, 1 booster per param)",
            xgb_configs(), make_xgb,
            folds, times, params, penalty, inst_features, baseline, shift)

    # ── LightGBM ───────────────────────────────────────────────────────────────
    if "lgbm" in methods:
        from lightgbm import LGBMRegressor

        def make_lgbm(kwargs):
            return MultiOutputRegressor(LGBMRegressor(**kwargs), n_jobs=1)

        all_results["LightGBM"] = run_method(
            "LightGBM (MultiOutputRegressor, 1 booster per param)",
            lgbm_configs(), make_lgbm,
            folds, times, params, penalty, inst_features, baseline, shift)

    # ── Per-method summary tables ──────────────────────────────────────────────
    for method_name, results in all_results.items():
        rep.h2(f"{method_name} — all configs sorted by SGM")
        rep.raw("")
        sorted_r = sorted(results, key=lambda r: r[1])
        rep.table(
            ["Config", "Pooled SGM", "Baseline SGM", "Speedup"],
            [[lbl, fmt(st), fmt(sb), speedup_str(sb, st)]
             for lbl, st, sb in sorted_r],
            align=["<", ">", ">", "<"],
        )

    # ── Grand comparison ───────────────────────────────────────────────────────
    rep.h2("Grand comparison — best of each method")
    rep.raw("")

    # Classification best results (from prior runs)
    prior = [
        ("cbc_default (baseline)",                 13.234, "—"),
        ("Best single param",                      12.923, "classification"),
        ("fbps DTree depth=4",                     12.316, "classification"),
        ("sklearn DTC depth=4",                    12.057, "classification"),
        ("RF classify (depth=8, n=200)",           10.346, "classification"),
        ("XGBoost classify (lr=0.1,d=6,n=200)",   10.009, "classification"),
        ("LightGBM classify (lr=0.05,l=31,n=200)", 10.010, "classification"),
    ]
    rows = [[m, fmt(s), src, "—"] for m, s, src in prior]

    for method_name, results in all_results.items():
        best = min(results, key=lambda r: r[1])
        rows.append([
            f"{method_name} regress ({best[0]})",
            fmt(best[1]),
            "regression",
            speedup_str(best[2], best[1]),
        ])

    rep.table(
        ["Method", "Pooled SGM", "Framing", "Speedup vs baseline"],
        rows,
        align=["<", ">", "<", "<"],
    )

    # ── Write & print ──────────────────────────────────────────────────────────
    report_text = rep.text()
    with open(out_file, "w") as f:
        f.write(report_text + "\n")
    print(report_text)
    print(f"\nReport written to: {out_file}", file=sys.stderr)


if __name__ == "__main__":
    main()
