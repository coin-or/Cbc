#!/usr/bin/env python3
"""
run_lp_dtree.py — Feature-based decision tree for LP parameter selection.

Reads lp_avg_times.csv (produced by analyze_lp_params.py) and instance features,
then builds an optimal decision tree (via the fbps greedy algorithm) that maps
instance characteristics to the best LP parameter setting.

The tree minimises the sum of solve times across all instances in each leaf,
using the penalty time (2×timelimit) for timed-out runs.

Usage:
    python3 run_lp_dtree.py --dir EXPERIMENT_DIR [options]

Key options:
    --features PATH     Path to features CSV (default: ~/inst/miplib/2017+spp/features.csv)
    --max-depth N       Maximum tree depth (default: 3)
    --min-leaf N        Minimum instances per leaf (default: 20)
    --baseline PARAM    Baseline param tag for speedup reporting (default: dual_default)
    --out FILE          Output text file (default: lp_dtree.txt in experiment dir)
    --svg FILE          Output SVG file (default: lp_dtree.svg in experiment dir)
    --exclude-prefix S  Exclude params starting with S (default: guess_)
"""

import sys
import os
import argparse
import csv
import math
from collections import defaultdict
from time import process_time

# Add fbps library path
_SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(_SCRIPT_DIR, 'fbps', 'fbps'))

from instance import InstanceSet, Instance
from psetting import PSetting
from results import Results, FillStrategy
from dtree import DTree, Node


def shifted_geomean(vals, shift=1.0):
    if not vals:
        return float("nan")
    return math.exp(sum(math.log(v + shift) for v in vals) / len(vals)) - shift


def format_val(v):
    if v > 10000:
        return f"{v:.0f}"
    if v > 100:
        return f"{v:.1f}"
    return f"{v:.3f}"


def print_tree(node, features, indent=0, out=None):
    """Recursively print the decision tree in human-readable format."""
    prefix = "  " * indent
    if out is None:
        out = sys.stdout

    if not node.children_nodes:
        # Leaf node
        ps = node.bestPS[0]
        n = len(node.instances)
        total_time = sum(inst.results[ps.idx] for inst in node.instances)
        avg_time = total_time / n if n else 0
        out.write(f"{prefix}→ RECOMMEND: {ps.setting}\n")
        out.write(f"{prefix}  Instances: {n}  |  Avg time: {format_val(avg_time)}s  |  Total: {format_val(total_time)}s\n")
    else:
        feat = features[node.branch_feat_idx]
        bval = node.branch_value
        bval_str = f"{bval:.6g}" if isinstance(bval, float) else str(bval)
        out.write(f"{prefix}Split on: {feat} ≤ {bval_str}\n")
        for i, child in enumerate(node.children_nodes):
            branch_label = f"≤ {bval_str}" if i == 0 else f"> {bval_str}"
            out.write(f"{prefix}  [{branch_label}]  ({len(child.instances)} instances)\n")
            print_tree(child, features, indent + 2, out)


def compute_speedup(tree, baseline_param, iset, results):
    """Compute total time with tree vs baseline param."""
    if baseline_param not in results.psettingByName:
        return None, None, None
    bp_idx = results.psettingByName[baseline_param].idx
    baseline_total = sum(inst.results[bp_idx] for inst in iset.instances)

    tree_total = 0.0
    for leaf in tree.leafs:
        ps_idx = leaf.bestPS[0].idx
        for inst in leaf.instances:
            tree_total += inst.results[ps_idx]

    speedup = baseline_total / tree_total if tree_total > 0 else float("inf")
    return baseline_total, tree_total, speedup


def write_leaf_summary(tree, features, results, penalty, out):
    """Write a summary table of all leaves."""
    out.write("\n── Leaf Summary ─────────────────────────────────────────────────────\n")
    out.write(f"{'Leaf':<6} {'Param':<40} {'#Inst':>6} {'SolvedPct':>10} {'AvgTime':>9} {'SGM':>9}\n")
    out.write("─" * 85 + "\n")
    for i, leaf in enumerate(tree.leafs):
        ps = leaf.bestPS[0]
        ps_idx = ps.idx
        times = [inst.results[ps_idx] for inst in leaf.instances]
        n_solved = sum(1 for t in times if t < penalty * 0.99)
        avg_t = sum(times) / len(times) if times else 0
        sgm = shifted_geomean(times)
        pct = 100 * n_solved / len(leaf.instances)
        out.write(f"{i+1:<6} {ps.setting:<40} {len(leaf.instances):>6} {pct:>9.1f}% {avg_t:>9.1f} {sgm:>9.1f}\n")


def write_split_conditions(node, features, path=None, out=None):
    """Print the conditions (feature splits) leading to each leaf."""
    if path is None:
        path = []
    if out is None:
        out = sys.stdout
    if not node.children_nodes:
        ps = node.bestPS[0]
        conditions = " AND ".join(path) if path else "(root)"
        out.write(f"  IF {conditions}\n  THEN {ps.setting}  ({len(node.instances)} instances)\n\n")
        return
    feat = features[node.branch_feat_idx]
    bval = node.branch_value
    bval_str = f"{bval:.6g}" if isinstance(bval, float) else str(bval)
    for i, child in enumerate(node.children_nodes):
        if i == 0:
            cond = f"{feat} ≤ {bval_str}"
        else:
            cond = f"{feat} > {bval_str}"
        write_split_conditions(child, features, path + [cond], out)


def build_results_from_csv(iset, avg_times_csv, exclude_prefix, penalty):
    """
    Build a Results-compatible object directly from lp_avg_times.csv,
    bypassing the fbps CSV reader (which expects only 3 columns and skips header).
    We load avg_wall_seconds directly — already penalized for timeouts.
    """
    # Collect psettings
    psetting_map = {}
    psettings = []

    # First pass: collect all param tags
    with open(avg_times_csv, newline="") as f:
        for i, row in enumerate(csv.DictReader(f)):
            p = row["param_tag"]
            if exclude_prefix and p.startswith(exclude_prefix):
                continue
            if p not in psetting_map:
                ps = PSetting(p)
                ps.idx = len(psettings)
                psettings.append(ps)
                psetting_map[p] = ps

    # Initialize instance results to penalty (handle missing combos)
    for inst in iset.instances:
        inst.results = [penalty] * len(psettings)

    # Second pass: fill in actual times
    with open(avg_times_csv, newline="") as f:
        for row in csv.DictReader(f):
            inst_name = row["instance"]
            p = row["param_tag"]
            if exclude_prefix and p.startswith(exclude_prefix):
                continue
            if inst_name not in iset.instByName:
                continue
            inst = iset.instByName[inst_name]
            ps = psetting_map[p]
            t = float(row["avg_wall_seconds"])
            inst.results[ps.idx] = t

    # Attach to a minimal Results-like object
    class _Results:
        pass
    r = _Results()
    r.psettings = psettings
    r.psettingByName = psetting_map
    return r


def generate_svg(tree, features, svg_path):
    """Generate graphviz SVG of the decision tree."""
    try:
        from graphviz import Digraph
    except ImportError:
        print("  (graphviz not installed — skipping SVG)")
        return

    def add_nodes(dot, node, parent_id=None, edge_label=""):
        node_id = node.node_id
        if not node.children_nodes:
            ps = node.bestPS[0]
            n = len(node.instances)
            label = f"{ps.setting}\n{n} instances"
            dot.node(node_id, label, shape="folder", style="filled", fillcolor="lightyellow")
        else:
            feat = features[node.branch_feat_idx]
            bval = node.branch_value
            bval_str = f"{bval:.5g}" if isinstance(bval, float) else str(bval)
            label = f"{feat}\n≤ {bval_str}?"
            dot.node(node_id, label, shape="box3d", style="filled", fillcolor="lightblue2")
        if parent_id:
            dot.edge(parent_id, node_id, label=edge_label)
        for i, child in enumerate(node.children_nodes):
            bval = node.branch_value
            bval_str = f"{bval:.5g}" if isinstance(bval, float) else str(bval)
            el = f"≤{bval_str}" if i == 0 else f">{bval_str}"
            add_nodes(dot, child, node_id, el)

    dot = Digraph(format="svg")
    dot.attr("node", shape="box", fontsize="11")
    dot.attr("graph", rankdir="TB")
    add_nodes(dot, tree.root)

    # graphviz.render writes to file; use source to write SVG directly
    base = svg_path.replace(".svg", "")
    dot.render(base, cleanup=True)
    # dot.render produces base.svg when format="svg"
    rendered = base + ".svg"
    if os.path.exists(rendered) and rendered != svg_path:
        os.rename(rendered, svg_path)
    print(f"  SVG saved to: {svg_path}")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dir", required=True, help="Experiment directory (contains lp_avg_times.csv)")
    ap.add_argument("--features", default=os.path.expanduser("~/inst/miplib/2017+spp/features.csv"))
    ap.add_argument("--max-depth", type=int, default=3)
    ap.add_argument("--min-leaf", type=int, default=20,
                    help="Minimum instances per leaf node (default: 20)")
    ap.add_argument("--baseline", default="dual_default")
    ap.add_argument("--timelimit", type=float, default=7200.0)
    ap.add_argument("--penalty-mult", type=float, default=2.0)
    ap.add_argument("--out", default=None)
    ap.add_argument("--svg", default=None)
    ap.add_argument("--exclude-prefix", default="guess_")
    args = ap.parse_args()

    exp_dir = args.dir
    avg_times_csv = os.path.join(exp_dir, "lp_avg_times.csv")
    if not os.path.isfile(avg_times_csv):
        sys.exit(f"Error: {avg_times_csv} not found — run analyze_lp_params.py first")
    if not os.path.isfile(args.features):
        sys.exit(f"Error: features file not found: {args.features}")

    penalty = args.timelimit * args.penalty_mult
    out_path = args.out or os.path.join(exp_dir, "lp_dtree.txt")
    svg_path = args.svg or os.path.join(exp_dir, "lp_dtree.svg")
    exp_name = os.path.basename(os.path.normpath(exp_dir))

    print(f"Loading features from: {args.features}")
    t0 = process_time()
    iset = InstanceSet(args.features)
    print(f"  {len(iset.instances)} instances, {len(iset.features)} features  ({process_time()-t0:.1f}s)")

    print(f"Loading results from: {avg_times_csv}")
    t0 = process_time()
    results = build_results_from_csv(iset, avg_times_csv, args.exclude_prefix, penalty)
    n_params = len(results.psettings)
    print(f"  {n_params} parameter settings  ({process_time()-t0:.1f}s)")

    print(f"Building decision tree (max_depth={args.max_depth}, min_leaf={args.min_leaf}) ...")
    t0 = process_time()
    tree = DTree(iset, results, args.max_depth, args.min_leaf, args.baseline
                 if args.baseline in results.psettingByName else "")
    tree.build()
    elapsed = process_time() - t0
    print(f"  Done in {elapsed:.1f}s  ({len(tree.leafs)} leaves)")

    baseline_total, tree_total, speedup = compute_speedup(tree, args.baseline, iset, results)

    # ── Write text report ───────────────────────────────────────────────────
    lines = []
    w = lambda s="": lines.append(s)

    w(f"LP Decision Tree Analysis  —  {exp_name}")
    w("=" * 70)
    w(f"  Instances  : {len(iset.instances)}")
    w(f"  Features   : {len(iset.features)}")
    w(f"  Parameters : {n_params}")
    w(f"  Max depth  : {args.max_depth}")
    w(f"  Min leaf   : {args.min_leaf}")
    w(f"  Penalty    : {penalty:.0f}s  (timelimit {args.timelimit:.0f}s × {args.penalty_mult})")
    w(f"  Baseline   : {args.baseline}")
    if baseline_total is not None:
        w(f"  Total time baseline : {baseline_total:.1f}s")
        w(f"  Total time with tree: {tree_total:.1f}s")
        w(f"  Speedup vs baseline : {speedup:.3f}×")
    w()

    w("── Decision Tree Structure ───────────────────────────────────────────")
    import io
    buf = io.StringIO()
    print_tree(tree.root, iset.features, out=buf)
    w(buf.getvalue())

    w("── Rule Summary (IF-THEN conditions) ────────────────────────────────")
    w()
    buf2 = io.StringIO()
    write_split_conditions(tree.root, iset.features, out=buf2)
    w(buf2.getvalue())

    buf3 = io.StringIO()
    write_leaf_summary(tree, iset.features, results, penalty, buf3)
    w(buf3.getvalue())

    w()
    w("── Per-leaf instance lists ───────────────────────────────────────────")
    for i, leaf in enumerate(tree.leafs):
        ps = leaf.bestPS[0]
        w(f"\nLeaf {i+1}: {ps.setting}  ({len(leaf.instances)} instances)")
        times = [inst.results[ps.idx] for inst in leaf.instances]
        for inst, t in sorted(zip(leaf.instances, times), key=lambda x: x[1]):
            solved = "✓" if t < penalty * 0.99 else "✗"
            w(f"  {solved} {inst.name:<45} {format_val(t)}s")

    report = "\n".join(lines)
    with open(out_path, "w") as f:
        f.write(report)
    print(f"\n→ Text report: {out_path}")
    print(report[:2000])  # preview

    print(f"\nGenerating SVG ...")
    generate_svg(tree, iset.features, svg_path)


if __name__ == "__main__":
    main()
