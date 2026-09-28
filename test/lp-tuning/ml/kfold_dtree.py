#!/usr/bin/env python3
"""
kfold_dtree.py — K-fold evaluation of feature-based decision tree LP parameter selection.

For tree depths 1, 2 and 3:
  For each fold k (0..K-1):
    - Training set : instances from all other K-1 folds
    - Test set     : instances in fold k
    - Selection    : build a DTree on the training set (uses fbps greedy algorithm)
    - Evaluation   : traverse the tree for each test instance using its features
                     to get the recommended param; compare SGM vs baseline

Key interpretation notes:
  - The DTree optimises total sum-of-times (not SGM) on the training set.
    Each leaf recommends the param with the lowest sum of avg_wall_seconds
    over its training instances.
  - When evaluating on test instances, we traverse the tree using feature
    values to reach a leaf and use that leaf's recommended param.
  - avg_wall_seconds already contains penalty times for failed/timed-out runs
    (produced by analyze_lp_params.py).
  - We report SGM (shifted geometric mean) for comparisons, consistent with
    the single-param analysis.

Usage:
    python3 kfold_dtree.py --dir EXPERIMENT_DIR [options]

Options:
    --dir DIR           Experiment directory (must contain lp_avg_times.csv)
    --features PATH     Features CSV (default: ~/inst/miplib/2017+spp/features.csv)
    --partitions DIR    Directory with fold_NN.txt files
                        (default: ~/inst/miplib/2017+spp/partitions)
    --depths LIST       Comma-separated depths to evaluate (default: 1,2,3)
    --min-leaf N        Minimum instances per leaf node (default: 20)
    --timelimit T       LP time limit in seconds (default: 10800)
    --penalty-mult M    Penalty multiplier (default: 2.0)
    --shift S           SGM shift in seconds (default: 1.0)
    --baseline PARAM    Baseline param for comparison (default: cbc_default)
    --out FILE          Output report file (default: kfold_dtree.txt in exp dir)
"""

import argparse
import csv
import math
import os
import sys
from collections import defaultdict

# ── fbps library ───────────────────────────────────────────────────────────────
_SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(_SCRIPT_DIR, '..', 'fbps', 'fbps'))

from instance import InstanceSet, Instance
from psetting import PSetting
from dtree import DTree, Node


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

def build_results(iset, avg_times_csv, exclude_prefix, penalty):
    """
    Build a minimal Results object from lp_avg_times.csv.
    avg_wall_seconds already encodes penalty for failed runs.
    """
    psetting_map = {}
    psettings = []

    with open(avg_times_csv, newline="") as f:
        for row in csv.DictReader(f):
            p = row["param_tag"]
            if exclude_prefix and p.startswith(exclude_prefix):
                continue
            if p not in psetting_map:
                ps = PSetting(p)
                ps.idx = len(psettings)
                psettings.append(ps)
                psetting_map[p] = ps

    # Initialise every instance result to penalty
    for inst in iset.instances:
        inst.results = [penalty] * len(psettings)

    with open(avg_times_csv, newline="") as f:
        for row in csv.DictReader(f):
            p = row["param_tag"]
            if exclude_prefix and p.startswith(exclude_prefix):
                continue
            name = row["instance"]
            if name not in iset.instByName:
                continue
            inst = iset.instByName[name]
            inst.results[psetting_map[p].idx] = float(row["avg_wall_seconds"])

    class _Results:
        pass
    r = _Results()
    r.psettings = psettings
    r.psettingByName = psetting_map
    return r


def load_folds(partitions_dir):
    folds = {}
    for fname in sorted(os.listdir(partitions_dir)):
        if not fname.startswith("fold_") or not fname.endswith(".txt"):
            continue
        idx = int(fname[5:7])
        with open(os.path.join(partitions_dir, fname)) as f:
            folds[idx] = [line.strip() for line in f if line.strip()]
    return folds


# ── Subset InstanceSet ─────────────────────────────────────────────────────────

def subset_iset(full_iset, names):
    """
    Return a new InstanceSet restricted to `names` (a set of instance name strings).
    Shares the same features list and feature metadata as full_iset.
    get_branching_values_feature() iterates self.instances dynamically, so
    the subset will compute valid branching values for the subset automatically.
    """
    sub = InstanceSet()                         # empty — no CSV loading
    sub.features = full_iset.features           # shared (read-only in tree building)
    sub.featureValues = full_iset.featureValues # shared
    sub.branchingOptions = full_iset.branchingOptions  # shared (not used in greedy_branch)
    sub.instances = [inst for inst in full_iset.instances if inst.name in names]
    sub.instByName = {inst.name: inst for inst in sub.instances}
    return sub


# ── Tree traversal ─────────────────────────────────────────────────────────────

def traverse(node, inst):
    """
    Walk the decision tree for a single instance using its feature values.
    Returns the recommended PSetting name (string).
    """
    if not node.children_nodes:
        return node.bestPS[0].setting
    feat_val = inst.features[node.branch_feat_idx]
    branch_val = node.branch_value
    # mirror the same comparison used in InstanceSet.branch()
    try:
        go_left = feat_val <= branch_val
    except TypeError:
        go_left = str(feat_val) <= str(branch_val)
    return traverse(node.children_nodes[0 if go_left else 1], inst)


# ── Per-fold tree evaluation ───────────────────────────────────────────────────

def eval_fold(fold_idx, folds, full_iset, results, depth, min_leaf, penalty,
              shift, baseline, param_to_idx):
    """
    Train DTree on all-but-fold-k, evaluate on fold-k.
    Returns (test_sgm_tree, test_sgm_baseline, leaf_assignments dict)
    """
    test_names  = set(folds[fold_idx])
    train_names = {inst for fi, insts in folds.items()
                   if fi != fold_idx for inst in insts}

    train_iset = subset_iset(full_iset, train_names)

    tree = DTree(train_iset, results, depth, min_leaf, default_setting="")
    tree.build()

    # Evaluate on test instances
    test_insts = [full_iset.instByName[n] for n in test_names
                  if n in full_iset.instByName]

    baseline_idx = param_to_idx.get(baseline)

    tree_times = []
    base_times = []
    assignments = {}   # instance_name -> recommended param

    for inst in test_insts:
        rec_param = traverse(tree.root, inst)
        rec_idx   = param_to_idx.get(rec_param)

        t_tree = inst.results[rec_idx]  if rec_idx  is not None else penalty
        t_base = inst.results[baseline_idx] if baseline_idx is not None else penalty

        tree_times.append(t_tree)
        base_times.append(t_base)
        assignments[inst.name] = rec_param

    sgm_tree = shifted_geomean(tree_times, shift)
    sgm_base = shifted_geomean(base_times, shift)
    return sgm_tree, sgm_base, assignments


# ── In-sample full-data tree (optimistic reference) ───────────────────────────

def eval_full_insample(full_iset, results, depth, min_leaf, penalty,
                       shift, baseline, param_to_idx):
    """Build tree on all instances, evaluate on all instances (optimistic)."""
    tree = DTree(full_iset, results, depth, min_leaf, default_setting="")
    tree.build()

    baseline_idx = param_to_idx.get(baseline)
    tree_times, base_times = [], []
    assignments = {}

    for inst in full_iset.instances:
        rec_param = traverse(tree.root, inst)
        rec_idx   = param_to_idx.get(rec_param)
        tree_times.append(inst.results[rec_idx]  if rec_idx  is not None else penalty)
        base_times.append(inst.results[baseline_idx] if baseline_idx is not None else penalty)
        assignments[inst.name] = rec_param

    return (shifted_geomean(tree_times, shift),
            shifted_geomean(base_times, shift),
            tree, assignments)


# ── Main ───────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--dir",         required=True)
    parser.add_argument("--features",
                        default=os.path.expanduser("~/inst/miplib/2017+spp/features.csv"))
    parser.add_argument("--partitions",
                        default=os.path.expanduser("~/inst/miplib/2017+spp/partitions"))
    parser.add_argument("--depths",      default="1,2,3")
    parser.add_argument("--min-leaf",    type=int, default=20)
    parser.add_argument("--timelimit",   type=float, default=10800.0)
    parser.add_argument("--penalty-mult", type=float, default=2.0)
    parser.add_argument("--shift",       type=float, default=1.0)
    parser.add_argument("--baseline",    default="cbc_default")
    parser.add_argument("--out",         default=None)
    args = parser.parse_args()

    exp_dir   = args.dir
    penalty   = args.timelimit * args.penalty_mult
    shift     = args.shift
    baseline  = args.baseline
    min_leaf  = args.min_leaf
    depths    = [int(d) for d in args.depths.split(",")]
    out_file  = args.out or os.path.join(exp_dir, "kfold_dtree.txt")
    avg_csv   = os.path.join(exp_dir, "lp_avg_times.csv")

    for path, label in [(avg_csv, "lp_avg_times.csv"),
                        (args.features, "features CSV"),
                        (args.partitions, "partitions dir")]:
        if not os.path.exists(path):
            sys.exit(f"Error: {label} not found: {path}")

    # ── Load data ──────────────────────────────────────────────────────────────
    print("Loading features ...", flush=True)
    full_iset = InstanceSet(args.features)
    print(f"  {len(full_iset.instances)} instances, {len(full_iset.features)} features")

    print("Loading results ...", flush=True)
    results = build_results(full_iset, avg_csv, "guess_", penalty)
    param_to_idx = {ps.setting: ps.idx for ps in results.psettings}
    print(f"  {len(results.psettings)} parameter settings")

    folds = load_folds(args.partitions)
    k = len(folds)
    print(f"  {k} folds loaded")

    rep = Report()
    rep.h1("K-Fold Evaluation: Feature-Based Decision Tree (LP Parameter Selection)")
    rep.raw(f"  Experiment : {exp_dir}")
    rep.raw(f"  Folds      : {k}  |  Min-leaf : {min_leaf}  |  Depths : {depths}")
    rep.raw(f"  Baseline   : {baseline}  |  Penalty : {penalty:.0f}s  |  SGM shift : {shift}s")

    # ── Cross-validation per depth ─────────────────────────────────────────────
    # summary_rows[depth] = (pooled_sgm_tree, pooled_sgm_base, speedup)
    summary_rows = []

    for depth in depths:
        rep.h2(f"Depth {depth} — per-fold cross-validation")
        rep.raw("")

        fold_headers = ["Fold", "N", "Train N", "#Leaves",
                        "Test SGM (tree)", "Test SGM (baseline)", "Speedup"]
        fold_align   = ["<", ">", ">", ">", ">", ">", "<"]
        fold_rows_   = []

        all_tree_times = []
        all_base_times = []

        for fold_idx in sorted(folds):
            test_names  = set(folds[fold_idx])
            train_names = {inst for fi, insts in folds.items()
                           if fi != fold_idx for inst in insts}

            train_iset = subset_iset(full_iset, train_names)
            train_insts_found = len(train_iset.instances)

            print(f"  depth={depth} fold={fold_idx:02d}: "
                  f"train={train_insts_found} test={len(test_names)} ...",
                  end=" ", flush=True)

            tree = DTree(train_iset, results, depth, min_leaf, default_setting="")
            tree.build()
            n_leaves = len(tree.leafs)
            print(f"{n_leaves} leaves")

            # Evaluate on test
            baseline_idx = param_to_idx.get(baseline)
            tree_times, base_times = [], []

            for name in sorted(test_names):
                if name not in full_iset.instByName:
                    continue
                inst = full_iset.instByName[name]
                rec_param = traverse(tree.root, inst)
                rec_idx   = param_to_idx.get(rec_param)

                tree_times.append(inst.results[rec_idx]  if rec_idx  is not None else penalty)
                base_times.append(inst.results[baseline_idx] if baseline_idx is not None else penalty)

            all_tree_times.extend(tree_times)
            all_base_times.extend(base_times)

            sgm_t = shifted_geomean(tree_times, shift)
            sgm_b = shifted_geomean(base_times, shift)

            fold_rows_.append([
                f"{fold_idx:02d}", len(test_names), train_insts_found, n_leaves,
                fmt(sgm_t), fmt(sgm_b), speedup_str(sgm_b, sgm_t),
            ])

        rep.table(fold_headers, fold_rows_, align=fold_align)

        pooled_tree = shifted_geomean(all_tree_times, shift)
        pooled_base = shifted_geomean(all_base_times, shift)
        rep.raw("")
        rep.raw(f"  Pooled test-set SGM  (tree depth={depth}) : {fmt(pooled_tree)}")
        rep.raw(f"  Pooled test-set SGM  (baseline)           : {fmt(pooled_base)}")
        rep.raw(f"  Pooled speedup                            : "
                f"{speedup_str(pooled_base, pooled_tree)}")

        summary_rows.append((depth, pooled_tree, pooled_base))

    # ── Summary across depths ──────────────────────────────────────────────────
    rep.h2("Summary: pooled out-of-sample speedup by depth")
    rep.raw("")

    sum_headers = ["Depth", "Pooled SGM (tree)", "Pooled SGM (baseline)", "Speedup vs baseline"]
    sum_align   = ["<", ">", ">", "<"]
    sum_table   = []
    for depth, st, sb in summary_rows:
        sum_table.append([depth, fmt(st), fmt(sb), speedup_str(sb, st)])
    rep.table(sum_headers, sum_table, align=sum_align)

    # ── In-sample full-data trees (optimistic reference) ──────────────────────
    rep.h2("In-sample full-data trees (optimistic / upper-bound reference)")
    rep.raw("")
    rep.raw("  Trees trained and evaluated on ALL 380 instances — no hold-out.")
    rep.raw("  Results are optimistic; use k-fold pooled SGM for unbiased estimates.\n")

    insample_headers = ["Depth", "#Leaves", "In-sample SGM (tree)",
                        "In-sample SGM (baseline)", "In-sample speedup"]
    insample_align   = ["<", ">", ">", ">", "<"]
    insample_rows_   = []

    for depth in depths:
        print(f"  In-sample tree depth={depth} ...", end=" ", flush=True)
        sgm_t, sgm_b, tree, _ = eval_full_insample(
            full_iset, results, depth, min_leaf, penalty, shift, baseline, param_to_idx)
        print(f"{len(tree.leafs)} leaves")
        insample_rows_.append([depth, len(tree.leafs),
                                fmt(sgm_t), fmt(sgm_b), speedup_str(sgm_b, sgm_t)])

    rep.table(insample_headers, insample_rows_, align=insample_align)

    # ── In-sample leaf detail for depth 3 ─────────────────────────────────────
    rep.h2("In-sample leaf assignments — depth 3 (for reference)")
    rep.raw("")
    print("  Building depth-3 in-sample tree for leaf detail ...", flush=True)
    _, _, tree3, assigns3 = eval_full_insample(
        full_iset, results, 3, min_leaf, penalty, shift, baseline, param_to_idx)

    import io
    buf = io.StringIO()
    _print_tree(tree3.root, full_iset.features, out=buf)
    for line in buf.getvalue().splitlines():
        rep.raw("  " + line)

    rep.h2("IF-THEN rules — depth 3 (in-sample tree)")
    rep.raw("")
    buf2 = io.StringIO()
    _write_rules(tree3.root, full_iset.features, out=buf2)
    for line in buf2.getvalue().splitlines():
        rep.raw("  " + line)

    # ── Write & print ──────────────────────────────────────────────────────────
    report_text = rep.text()
    with open(out_file, "w") as f:
        f.write(report_text + "\n")

    print(report_text)
    print(f"\nReport written to: {out_file}", file=sys.stderr)


# ── Tree printing helpers ──────────────────────────────────────────────────────

def _fmt_val(v):
    if isinstance(v, float):
        return f"{v:.6g}"
    return str(v)


def _print_tree(node, features, indent=0, out=None):
    if out is None:
        out = sys.stdout
    prefix = "  " * indent
    if not node.children_nodes:
        ps = node.bestPS[0]
        n = len(node.instances)
        total = sum(inst.results[ps.idx] for inst in node.instances)
        avg   = total / n if n else 0
        out.write(f"{prefix}→ RECOMMEND: {ps.setting}\n")
        out.write(f"{prefix}  Instances: {n}  Avg time: {avg:.3f}s\n")
    else:
        feat  = features[node.branch_feat_idx]
        bval  = _fmt_val(node.branch_value)
        out.write(f"{prefix}Split: {feat} ≤ {bval}\n")
        for i, child in enumerate(node.children_nodes):
            label = f"≤{bval}" if i == 0 else f">{bval}"
            out.write(f"{prefix}  [{label}] ({len(child.instances)} instances)\n")
            _print_tree(child, features, indent + 2, out)


def _write_rules(node, features, path=None, out=None):
    if path is None:
        path = []
    if out is None:
        out = sys.stdout
    if not node.children_nodes:
        cond = " AND ".join(path) if path else "(all instances)"
        out.write(f"IF {cond}\n   THEN {node.bestPS[0].setting}"
                  f"  ({len(node.instances)} training instances)\n\n")
        return
    feat = features[node.branch_feat_idx]
    bval = _fmt_val(node.branch_value)
    for i, child in enumerate(node.children_nodes):
        cond = f"{feat} {'≤' if i == 0 else '>'} {bval}"
        _write_rules(child, features, path + [cond], out)


if __name__ == "__main__":
    main()
