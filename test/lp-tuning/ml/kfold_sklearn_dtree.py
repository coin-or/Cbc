#!/usr/bin/env python3
"""
kfold_sklearn_dtree.py — K-fold evaluation using scikit-learn DecisionTreeClassifier.

Compares with the fbps-based kfold_dtree.py. Key differences:
  - fbps tree   : minimises sum of avg_wall_seconds in each leaf (regression on times)
  - sklearn DTC : predicts the *best* param label per instance (Gini classification)

Label assignment: for each training instance, the label is the param with the lowest
avg_wall_seconds (already penalised for failures/timeouts). The DTC learns to predict
that label from instance features.

At test time, the predicted param is looked up in the actual avg_wall_seconds table
to compute the real solve time for SGM comparison.

Note: sklearn's standard DecisionTreeClassifier uses axis-aligned (univariate) splits,
identical to fbps in that regard. The difference is purely in the split criterion
(Gini impurity vs sum-of-times).

Usage:
    python3 kfold_sklearn_dtree.py --dir EXPERIMENT_DIR [options]

Options:
    --dir DIR           Experiment directory (must contain lp_avg_times.csv)
    --features PATH     Features CSV (default: ~/inst/miplib/2017+spp/features.csv)
    --partitions DIR    Fold files directory
                        (default: ~/inst/miplib/2017+spp/partitions)
    --depths LIST       Comma-separated max_depth values (default: 1,2,3,4)
    --timelimit T       LP time limit in seconds (default: 10800)
    --penalty-mult M    Penalty multiplier (default: 2.0)
    --shift S           SGM shift in seconds (default: 1.0)
    --baseline PARAM    Baseline param for comparison (default: cbc_default)
    --criterion C       DTC split criterion: gini or entropy (default: gini)
    --min-samples-leaf N  Minimum samples per leaf (default: 5)
    --out FILE          Output report file (default: kfold_sklearn_dtree.txt in exp dir)
"""

import argparse
import csv
import io
import math
import os
import sys
from collections import defaultdict, Counter

import numpy as np
from sklearn.tree import DecisionTreeClassifier, export_text


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
    """
    Returns:
        feature_names : list of str
        inst_features : dict  name -> np.ndarray (float64, NaN for missing)
    """
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
    """Returns times[instance][param] = penalised avg wall seconds, plus sorted lists."""
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
    """Return the param with lowest avg time for this instance."""
    best_p, best_t = None, float("inf")
    for p in params:
        t = times[inst_name].get(p, penalty)
        if t < best_t:
            best_t, best_p = t, p
    return best_p


# ── Build X / y matrices for a set of instances ───────────────────────────────

def build_Xy(inst_names, times, params, penalty, inst_features, feature_names):
    """
    Build feature matrix X and label vector y for sklearn.
    Instances without feature data are skipped.
    Returns X (n×d), y (n,), valid_names (list of n instance names).
    """
    X_rows, y_rows, valid = [], [], []
    for name in inst_names:
        if name not in inst_features:
            continue
        feats = inst_features[name]
        # Replace NaN with column median — computed per-call (small datasets)
        # We'll impute after collecting all rows
        X_rows.append(feats)
        y_rows.append(best_param_for(name, times, params, penalty))
        valid.append(name)

    if not X_rows:
        return None, None, []

    X = np.vstack(X_rows)
    # Simple NaN imputation: replace with column median of this subset
    col_medians = np.nanmedian(X, axis=0)
    nan_mask = np.isnan(X)
    X[nan_mask] = np.take(col_medians, np.where(nan_mask)[1])

    return X, np.array(y_rows), valid


# ── Main ───────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--dir",          required=True)
    parser.add_argument("--features",
                        default=os.path.expanduser("~/inst/miplib/2017+spp/features.csv"))
    parser.add_argument("--partitions",
                        default=os.path.expanduser("~/inst/miplib/2017+spp/partitions"))
    parser.add_argument("--depths",       default="1,2,3,4")
    parser.add_argument("--timelimit",    type=float, default=10800.0)
    parser.add_argument("--penalty-mult", type=float, default=2.0)
    parser.add_argument("--shift",        type=float, default=1.0)
    parser.add_argument("--baseline",     default="cbc_default")
    parser.add_argument("--criterion",    default="gini", choices=["gini", "entropy"])
    parser.add_argument("--min-samples-leaf", type=int, default=5)
    parser.add_argument("--out",          default=None)
    args = parser.parse_args()

    exp_dir   = args.dir
    penalty   = args.timelimit * args.penalty_mult
    shift     = args.shift
    baseline  = args.baseline
    depths    = [int(d) for d in args.depths.split(",")]
    out_file  = args.out or os.path.join(exp_dir, "kfold_sklearn_dtree.txt")
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
    rep.h1("K-Fold Evaluation: scikit-learn DecisionTreeClassifier")
    rep.raw(f"  Experiment      : {exp_dir}")
    rep.raw(f"  Folds           : {k}  |  Depths : {depths}")
    rep.raw(f"  Criterion       : {args.criterion}  |  Min-samples-leaf : {args.min_samples_leaf}")
    rep.raw(f"  Baseline        : {baseline}  |  Penalty : {penalty:.0f}s  |  SGM shift : {shift}s")
    rep.raw("")
    rep.raw("  Label per instance = param with lowest avg_wall_seconds (penalised).")
    rep.raw("  DTC learns to predict that label from instance features (Gini criterion).")
    rep.raw("  NaN features imputed with per-subset column median.")

    baseline_idx_in_params = params.index(baseline) if baseline in params else None

    summary_rows = []

    for depth in depths:
        rep.h2(f"Depth {depth} — per-fold cross-validation")
        rep.raw("")

        fold_headers = ["Fold", "N", "Train N", "#Leaves",
                        "Test SGM (tree)", "Test SGM (baseline)", "Speedup"]
        fold_align   = ["<", ">", ">", ">", ">", ">", "<"]
        fold_rows_   = []

        all_tree_times, all_base_times = [], []

        for fold_idx in sorted(folds):
            test_names  = folds[fold_idx]
            train_names = [inst for fi, insts in folds.items()
                           if fi != fold_idx for inst in insts]

            # Build train matrix
            X_tr, y_tr, tr_valid = build_Xy(
                train_names, times, params, penalty, inst_features, feature_names)
            if X_tr is None:
                print(f"  depth={depth} fold={fold_idx:02d}: no training data, skipping")
                continue

            # Train
            clf = DecisionTreeClassifier(
                max_depth=depth,
                criterion=args.criterion,
                min_samples_leaf=args.min_samples_leaf,
                random_state=42,
            )
            clf.fit(X_tr, y_tr)
            n_leaves = clf.get_n_leaves()

            # Build test matrix — use train median for NaN imputation to avoid leakage
            X_te_rows, te_valid = [], []
            for name in test_names:
                if name not in inst_features:
                    continue
                X_te_rows.append(inst_features[name])
                te_valid.append(name)

            if not X_te_rows:
                continue

            X_te = np.vstack(X_te_rows)
            # Impute NaN using training-set column medians (no leakage)
            tr_medians = np.nanmedian(X_tr, axis=0)
            nan_mask = np.isnan(X_te)
            X_te[nan_mask] = np.take(tr_medians, np.where(nan_mask)[1])

            predictions = clf.predict(X_te)

            tree_times, base_times = [], []
            for name, rec_param in zip(te_valid, predictions):
                t_tree = times[name].get(rec_param, penalty)
                t_base = times[name].get(baseline,  penalty)
                tree_times.append(t_tree)
                base_times.append(t_base)

            all_tree_times.extend(tree_times)
            all_base_times.extend(base_times)

            sgm_t = shifted_geomean(tree_times, shift)
            sgm_b = shifted_geomean(base_times, shift)

            print(f"  depth={depth} fold={fold_idx:02d}: "
                  f"train={len(tr_valid)} test={len(te_valid)} "
                  f"leaves={n_leaves}  SGM {fmt(sgm_t)} vs {fmt(sgm_b)}")

            fold_rows_.append([
                f"{fold_idx:02d}", len(te_valid), len(tr_valid), n_leaves,
                fmt(sgm_t), fmt(sgm_b), speedup_str(sgm_b, sgm_t),
            ])

        rep.table(fold_headers, fold_rows_, align=fold_align)

        pooled_tree = shifted_geomean(all_tree_times, shift)
        pooled_base = shifted_geomean(all_base_times, shift)
        rep.raw("")
        rep.raw(f"  Pooled test-set SGM  (sklearn DTC depth={depth}) : {fmt(pooled_tree)}")
        rep.raw(f"  Pooled test-set SGM  (baseline)                  : {fmt(pooled_base)}")
        rep.raw(f"  Pooled speedup                                   : "
                f"{speedup_str(pooled_base, pooled_tree)}")

        summary_rows.append((depth, pooled_tree, pooled_base))

    # ── Summary ────────────────────────────────────────────────────────────────
    rep.h2("Summary: pooled out-of-sample speedup by depth")
    rep.raw("")
    rep.table(
        ["Depth", "Pooled SGM (sklearn DTC)", "Pooled SGM (baseline)", "Speedup"],
        [[d, fmt(st), fmt(sb), speedup_str(sb, st)] for d, st, sb in summary_rows],
        align=["<", ">", ">", "<"],
    )

    # ── Comparison with fbps ───────────────────────────────────────────────────
    fbps_results = {1: 13.327, 2: 12.685, 3: 12.399, 4: 12.316}
    rep.h2("Comparison with fbps decision tree (same experiment, same folds)")
    rep.raw("")
    rep.raw("  fbps optimises sum-of-times in each leaf (regression objective).")
    rep.raw("  sklearn DTC predicts best-param label per instance (Gini classification).\n")
    cmp_rows = []
    for d, st, sb in summary_rows:
        fbps_sgm = fbps_results.get(d)
        fbps_str = fmt(fbps_sgm) if fbps_sgm else "N/A"
        winner = "fbps" if (fbps_sgm and fbps_sgm <= st) else "sklearn"
        cmp_rows.append([d, fmt(st), fbps_str, fmt(sb), winner])
    rep.table(
        ["Depth", "sklearn SGM", "fbps SGM", "Baseline SGM", "Winner"],
        cmp_rows, align=["<", ">", ">", ">", "<"],
    )

    # ── In-sample full-data DTC + rule printout ────────────────────────────────
    rep.h2("In-sample full-data sklearn DTC (optimistic / upper-bound reference)")
    rep.raw("")
    insample_rows = []
    best_depth_tree = None
    best_depth_feat_names = None

    for depth in depths:
        X_all, y_all, all_valid = build_Xy(
            all_instances, times, params, penalty, inst_features, feature_names)
        clf = DecisionTreeClassifier(
            max_depth=depth, criterion=args.criterion,
            min_samples_leaf=args.min_samples_leaf, random_state=42)
        clf.fit(X_all, y_all)
        preds = clf.predict(X_all)
        tree_t = [times[n].get(p, penalty) for n, p in zip(all_valid, preds)]
        base_t = [times[n].get(baseline, penalty) for n in all_valid]
        sgm_t = shifted_geomean(tree_t, shift)
        sgm_b = shifted_geomean(base_t, shift)
        insample_rows.append([depth, clf.get_n_leaves(),
                               fmt(sgm_t), fmt(sgm_b), speedup_str(sgm_b, sgm_t)])
        if depth == max(depths):
            best_depth_tree = clf
            best_depth_feat_names = feature_names

    rep.table(
        ["Depth", "#Leaves", "In-sample SGM", "Baseline SGM", "Speedup"],
        insample_rows, align=["<", ">", ">", ">", "<"],
    )

    # Print the deepest in-sample tree rules
    if best_depth_tree is not None:
        rep.h2(f"Decision rules — depth {max(depths)} in-sample tree")
        rep.raw("")
        tree_text = export_text(best_depth_tree,
                                feature_names=best_depth_feat_names,
                                max_depth=max(depths))
        for line in tree_text.splitlines():
            rep.raw("  " + line)

        rep.h2(f"Label distribution at leaves — depth {max(depths)} in-sample tree")
        rep.raw("")
        leaf_ids = best_depth_tree.apply(X_all)
        from collections import defaultdict as dd
        leaf_preds = dd(list)
        for lid, pred in zip(leaf_ids, preds):
            leaf_preds[lid].append(pred)
        for lid, pred_list in sorted(leaf_preds.items()):
            freq = Counter(pred_list)
            top = freq.most_common(1)[0]
            rep.raw(f"  Leaf {lid:>3d} ({len(pred_list):>3d} instances): "
                    f"{top[0]}  ({top[1]}/{len(pred_list)} agree)")

    # ── Write & print ──────────────────────────────────────────────────────────
    report_text = rep.text()
    with open(out_file, "w") as f:
        f.write(report_text + "\n")
    print(report_text)
    print(f"\nReport written to: {out_file}", file=sys.stderr)


if __name__ == "__main__":
    main()
