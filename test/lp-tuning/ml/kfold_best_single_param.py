#!/usr/bin/env python3
"""
kfold_best_single_param.py — K-fold evaluation of the "best single parameter" strategy.

For each fold k (0..K-1):
  - Training set  : instances from all folds except k
  - Test set      : instances in fold k
  - Selection     : pick the param with the lowest shifted-geomean (SGM) on the
                    training set, among params that return correct results
  - Evaluation    : measure SGM of the selected param on the test set, compared
                    to the baseline (default: cbc_default)

This gives an unbiased (out-of-sample) estimate of how much we gain by switching
from the current default to the best single alternative LP parameter.

Additionally reports:
  - Oracle upper bound: what SGM would be if we knew the best param per instance
    (virtual per-instance oracle, in-sample on full set)
  - Full-data selection: best param using all 380 instances (in-sample, optimistic)

Usage:
    python3 kfold_best_single_param.py --dir EXPERIMENT_DIR [options]

Options:
    --dir DIR           Experiment directory (must contain lp_avg_times.csv)
    --partitions DIR    Directory with fold_NN.txt files
                        (default: ~/inst/miplib/2017+spp/partitions)
    --timelimit T       LP time limit in seconds (default: 10800)
    --penalty-mult M    Penalty multiplier (default: 2.0)
    --shift S           SGM shift in seconds (default: 1.0)
    --baseline PARAM    Baseline param for comparison (default: cbc_default)
    --out FILE          Output report file (default: kfold_best_single.txt in exp dir)
"""

import argparse
import csv
import math
import os
import sys
from collections import defaultdict


# ── Helpers ────────────────────────────────────────────────────────────────────

def shifted_geomean(values, shift=1.0):
    if not values:
        return float("nan")
    return math.exp(sum(math.log(v + shift) for v in values) / len(values)) - shift


def fmt(v, decimals=3):
    if math.isnan(v):
        return "N/A"
    return f"{v:.{decimals}f}"


def speedup_str(baseline_sgm, param_sgm):
    """Return '1.23x faster' or '0.98x (slower)' string."""
    if math.isnan(baseline_sgm) or math.isnan(param_sgm) or param_sgm == 0:
        return "N/A"
    ratio = baseline_sgm / param_sgm
    if ratio >= 1.0:
        return f"{ratio:.3f}x faster"
    else:
        return f"{ratio:.3f}x  (slower)"


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


# ── Core logic ─────────────────────────────────────────────────────────────────

def load_avg_times(csv_path, penalty):
    """
    Returns:
        times[instance][param] = penalised avg wall seconds
        params                 = sorted list of param tags
        instances              = sorted list of instance names
    """
    times = defaultdict(dict)
    with open(csv_path, newline="") as f:
        for row in csv.DictReader(f):
            inst  = row["instance"]
            param = row["param_tag"]
            try:
                t = float(row["avg_wall_seconds"])
            except (ValueError, KeyError):
                t = penalty
            times[inst][param] = t

    instances = sorted(times)
    params    = sorted({p for d in times.values() for p in d})
    return times, instances, params


def load_folds(partitions_dir):
    """Returns dict: fold_index -> [instance_name, ...]"""
    folds = {}
    for fname in sorted(os.listdir(partitions_dir)):
        if not fname.startswith("fold_") or not fname.endswith(".txt"):
            continue
        idx = int(fname[5:7])
        path = os.path.join(partitions_dir, fname)
        with open(path) as f:
            instances = [line.strip() for line in f if line.strip()]
        folds[idx] = instances
    return folds


def best_param_on(instance_subset, times, params, penalty, shift):
    """Select the param with the lowest SGM over instance_subset."""
    best_p, best_sgm = None, float("inf")
    for p in params:
        vals = [times[inst].get(p, penalty) for inst in instance_subset]
        sgm  = shifted_geomean(vals, shift)
        if sgm < best_sgm:
            best_sgm, best_p = sgm, p
    return best_p, best_sgm


def sgm_on(instance_subset, times, param, penalty, shift):
    vals = [times[inst].get(param, penalty) for inst in instance_subset]
    return shifted_geomean(vals, shift)


def solve_rate(instance_subset, times, param, penalty):
    solved = sum(1 for inst in instance_subset
                 if times[inst].get(param, penalty) < penalty)
    return solved / len(instance_subset) if instance_subset else float("nan")


# ── Main ───────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--dir",         required=True,
                        help="Experiment directory containing lp_avg_times.csv")
    parser.add_argument("--partitions",
                        default=os.path.expanduser("~/inst/miplib/2017+spp/partitions"),
                        help="Directory with fold_NN.txt partition files")
    parser.add_argument("--timelimit",   type=float, default=10800.0)
    parser.add_argument("--penalty-mult", type=float, default=2.0)
    parser.add_argument("--shift",       type=float, default=1.0)
    parser.add_argument("--baseline",    default="cbc_default")
    parser.add_argument("--out",         default=None)
    args = parser.parse_args()

    exp_dir     = args.dir
    penalty     = args.timelimit * args.penalty_mult
    shift       = args.shift
    baseline    = args.baseline
    out_file    = args.out or os.path.join(exp_dir, "kfold_best_single.txt")

    avg_csv = os.path.join(exp_dir, "lp_avg_times.csv")
    if not os.path.isfile(avg_csv):
        sys.exit(f"Error: {avg_csv} not found.\n"
                 "Run analyze_lp_params.py first to generate lp_avg_times.csv.")

    if not os.path.isdir(args.partitions):
        sys.exit(f"Error: partitions directory not found: {args.partitions}")

    times, all_instances, params = load_avg_times(avg_csv, penalty)
    folds = load_folds(args.partitions)
    k = len(folds)

    # Filter instances to those present in both the CSV and the partition files
    fold_instances = {inst for fold in folds.values() for inst in fold}
    missing = fold_instances - set(all_instances)
    if missing:
        print(f"Warning: {len(missing)} instances in folds not found in CSV "
              f"(they will use penalty): {sorted(missing)[:5]}...", file=sys.stderr)

    rep = Report()
    rep.h1("K-Fold Evaluation: Best Single LP Parameter")
    rep.raw(f"  Experiment : {exp_dir}")
    rep.raw(f"  Avg-times  : {avg_csv}")
    rep.raw(f"  Folds      : {k}  ({len(all_instances)} total instances, "
            f"{len(fold_instances)} in partition files)")
    rep.raw(f"  Parameters : {len(params)}")
    rep.raw(f"  Baseline   : {baseline}")
    rep.raw(f"  Time limit : {args.timelimit:.0f}s  "
            f"Penalty : {penalty:.0f}s  "
            f"SGM shift : {shift}s")

    # ── Per-fold cross-validation ──────────────────────────────────────────────
    rep.h2("Per-fold cross-validation results")
    rep.raw("")
    rep.raw("  For each fold: best param is selected on the 9-fold training set,")
    rep.raw("  then evaluated (SGM) on the held-out test fold.\n")

    fold_headers = ["Fold", "Test N", "Selected param (train)",
                    "Train SGM", "Test SGM (selected)", "Test SGM (baseline)",
                    "Speedup vs baseline"]
    fold_align   = ["<", ">", "<", ">", ">", ">", "<"]
    fold_rows    = []

    test_sgms_selected = []
    test_sgms_baseline = []
    selected_params    = []

    for fold_idx in sorted(folds):
        test_set  = [inst for inst in folds[fold_idx] if inst in times]
        train_set = [inst for fold_j, fold_insts in folds.items()
                     if fold_j != fold_idx
                     for inst in fold_insts
                     if inst in times]

        best_p, train_sgm = best_param_on(train_set, times, params, penalty, shift)

        test_sgm_sel  = sgm_on(test_set, times, best_p,  penalty, shift)
        test_sgm_base = sgm_on(test_set, times, baseline, penalty, shift)

        test_sgms_selected.append(test_sgm_sel)
        test_sgms_baseline.append(test_sgm_base)
        selected_params.append(best_p)

        fold_rows.append([
            f"{fold_idx:02d}",
            len(test_set),
            best_p,
            fmt(train_sgm),
            fmt(test_sgm_sel),
            fmt(test_sgm_base),
            speedup_str(test_sgm_base, test_sgm_sel),
        ])

    rep.table(fold_headers, fold_rows, align=fold_align)

    # ── Aggregate cross-validation summary ────────────────────────────────────
    rep.h2("Aggregate cross-validation summary")

    # Concatenate all test-fold instances to compute a single pooled SGM
    all_test_times_sel  = []
    all_test_times_base = []
    for fold_idx in sorted(folds):
        test_set = [inst for inst in folds[fold_idx] if inst in times]
        # Use the param selected for THIS fold
        best_p = selected_params[fold_idx]
        for inst in test_set:
            all_test_times_sel.append(times[inst].get(best_p,  penalty))
            all_test_times_base.append(times[inst].get(baseline, penalty))

    pooled_sgm_sel  = shifted_geomean(all_test_times_sel,  shift)
    pooled_sgm_base = shifted_geomean(all_test_times_base, shift)
    mean_fold_sgm_sel  = sum(test_sgms_selected) / k
    mean_fold_sgm_base = sum(test_sgms_baseline) / k

    from collections import Counter
    param_freq = Counter(selected_params)
    most_selected = param_freq.most_common()

    rep.raw("")
    rep.raw(f"  Pooled test-set SGM  (selected param) : {fmt(pooled_sgm_sel)}")
    rep.raw(f"  Pooled test-set SGM  (baseline)       : {fmt(pooled_sgm_base)}")
    rep.raw(f"  Pooled speedup                        : "
            f"{speedup_str(pooled_sgm_base, pooled_sgm_sel)}")
    rep.raw("")
    rep.raw(f"  Mean per-fold SGM (selected)  : {fmt(mean_fold_sgm_sel)}")
    rep.raw(f"  Mean per-fold SGM (baseline)  : {fmt(mean_fold_sgm_base)}")
    rep.raw("")
    rep.raw("  Param selected per fold (frequency):")
    for p, cnt in most_selected:
        rep.raw(f"    {p:<35s} selected in {cnt}/{k} folds")

    # ── In-sample full-data selection (optimistic reference) ──────────────────
    rep.h2("In-sample full-data selection (optimistic / upper-bound reference)")
    rep.raw("")
    rep.raw("  Best param selected using ALL instances (no hold-out).")
    rep.raw("  This is optimistic — the same data was used for selection and evaluation.\n")

    best_p_full, best_sgm_full = best_param_on(
        [i for i in all_instances if i in times], times, params, penalty, shift)
    baseline_sgm_full = sgm_on(
        [i for i in all_instances if i in times], times, baseline, penalty, shift)

    rep.raw(f"  Best param (full data) : {best_p_full}")
    rep.raw(f"  SGM (best param, all)  : {fmt(best_sgm_full)}")
    rep.raw(f"  SGM (baseline, all)    : {fmt(baseline_sgm_full)}")
    rep.raw(f"  In-sample speedup      : {speedup_str(baseline_sgm_full, best_sgm_full)}")

    # ── Full ranking on training-all (for reference) ──────────────────────────
    rep.h2("Full ranking of all params on complete instance set (in-sample)")
    rep.raw("")

    ranked = []
    inst_list = [i for i in all_instances if i in times]
    for p in params:
        sgm = sgm_on(inst_list, times, p, penalty, shift)
        sr  = solve_rate(inst_list, times, p, penalty)
        ranked.append((sgm, p, sr))
    ranked.sort()

    rank_rows = []
    for rank, (sgm, p, sr) in enumerate(ranked, 1):
        tag = ""
        if p == baseline:
            tag = " ← baseline"
        elif p == best_p_full:
            tag = " ← in-sample best"
        rank_rows.append([rank, p, fmt(sgm), f"{sr*100:.1f}%", tag])

    rep.table(["Rank", "Param", "SGM(all)s", "%Solved", ""],
              rank_rows, align=["<", "<", ">", ">", "<"])

    # ── Per-fold solve rates ──────────────────────────────────────────────────
    rep.h2("Test-fold solve rates (selected param vs baseline)")
    rep.raw("")
    sr_rows = []
    for fold_idx in sorted(folds):
        test_set = [inst for inst in folds[fold_idx] if inst in times]
        best_p   = selected_params[fold_idx]
        sr_sel   = solve_rate(test_set, times, best_p,   penalty)
        sr_base  = solve_rate(test_set, times, baseline, penalty)
        sr_rows.append([
            f"{fold_idx:02d}", len(test_set), best_p,
            f"{sr_sel*100:.1f}%", f"{sr_base*100:.1f}%",
        ])
    rep.table(["Fold", "N", "Selected param", "%Solved (selected)", "%Solved (baseline)"],
              sr_rows, align=["<", ">", "<", ">", ">"])

    # ── Write & print ──────────────────────────────────────────────────────────
    report_text = rep.text()
    with open(out_file, "w") as f:
        f.write(report_text + "\n")

    print(report_text)
    print(f"\nReport written to: {out_file}", file=sys.stderr)


if __name__ == "__main__":
    main()
