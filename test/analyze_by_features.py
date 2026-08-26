#!/usr/bin/env python3
"""
analyze_by_features.py — Break a compare_benchmarks-style comparison down by
matrix/instance features (set partitioning/packing/covering structure, pure
binary vs. mixed-integer, instance size, density, presence of general
integer variables, ...) instead of only reporting one aggregate number.

Motivation: an aggregate avg_cost/avg_pgap/avg_dgap across all instances can
hide the fact that a build is much stronger on some instance structures than
others (e.g. dominant on binary/SPP-structured or very large instances, but
actually weaker -- fewer solved, worse gap, *slower* -- on instances with
many general-integer variables). This script reuses compare_benchmarks.py's
own loading/classification/gap/cost logic (so results are always consistent
with `compare-benchmarks`), joins each instance against a matrix-features
table, and reports the same metrics compare-benchmarks reports (n, n_solved,
avg_cost, avg primal/dual gap -- gap_cap-penalized for outright failures, an
average runtime formatted as HhMMm) per feature-defined group, for each
experiment directory given.

Feature file: any CSV/TSV with a "Name" column (instance name) and the
~208-column feature vector Cbc's own `-writeFeatures -csvFeatures <file>`
produces (percBin, percInteger, rPercPartitioning/rPercPacking/rPercCovering,
rows, cols, density, genInt, ...). Two are auto-detected if --features-file
is not given (first match wins):
    Cbc/test/mip-sanity-data/features.tsv   (this repo's own generated set)
    ~/inst/miplib/2017+spp/features.csv     (broader MIPLIB feature set, if present)
If an instance used in the experiment is missing from the feature file, it
is simply excluded from every feature-based group (but still counted in the
"ALL" row) -- regenerate/extend the feature file with:
    cbc <instance>.mps.gz -writeFeatures -csvFeatures out.csv -exit

Usage:
    analyze_by_features.py <dir1> <dir2> [... <dirN>] [options]
    analyze_by_features.py label1=<dir1> label2=<dir2> [options]
    analyze_by_features.py --list experiments.txt [options]

Options:
    --weights FILE         JSON weights file (default: same as compare_benchmarks.py)
    --bks-tsv FILE         bks.tsv for results.tsv-style dirs lacking an embedded BKS
    --features-file FILE   instance-features CSV/TSV (see above for auto-detection)
    --groups FILE           JSON file overriding the built-in group definitions
                            (see BUILTIN_GROUPS below for the expected shape)
    --no-color
    -o/--output FILE       write the report to FILE instead of stdout
"""

import sys
import os
import json
import argparse

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import compare_benchmarks as cb

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

DEFAULT_FEATURES_CANDIDATES = [
    os.path.expanduser("~/inst/miplib/2017+spp/features.csv"),
    os.path.join(SCRIPT_DIR, "mip-sanity-data", "features.tsv"),
]

# ─────────────────────────────────────────────────────────────────────────────
# Group definitions
# ─────────────────────────────────────────────────────────────────────────────
# Each group is (name, pandas query expression string evaluated against the
# joined/feature dataframe via DataFrame.eval-style df.query()). Kept as
# plain query strings (not lambdas) so --groups can override them from JSON
# without embedding Python code.

BUILTIN_GROUPS = [
    ("ALL", None),
    ("percBin>=95% (near-pure binary)", "percBin >= 95"),
    ("percBin<95%", "percBin < 95"),
    ("SPP/pack/cov rows>=90%", "_sppc >= 90"),
    ("SPP/pack/cov rows<90%", "_sppc < 90"),
    ("Pure SPP/Cov + pure binary (0/1 matrix)", "_sppc >= 99.9 and percBin >= 99.9"),
    ("Everything else (not pure SPP/Cov+bin)", "not (_sppc >= 99.9 and percBin >= 99.9)"),
    ("percInteger>=95% (integer-heavy)", "percInteger >= 95"),
    ("percInteger<95%", "percInteger < 95"),
    ("Pure integer (100%)", "percInteger >= 100"),
    ("Mixed (has continuous vars)", "percInteger < 100"),
    ("Knapsack-dominant rows (>=30%)", "_kp >= 30"),
    ("Knapsack rows<30%", "_kp < 30"),
    ("density>=5%", "density >= 5"),
    ("density<5%", "density < 5"),
    ("rows<500 (small)", "rows < 500"),
    ("500<=rows<5000 (medium)", "rows >= 500 and rows < 5000"),
    ("5000<=rows<50000 (large)", "rows >= 5000 and rows < 50000"),
    ("rows>=50000 (huge)", "rows >= 50000"),
    ("cols<500", "cols < 500"),
    ("500<=cols<5000", "cols >= 500 and cols < 5000"),
    ("5000<=cols<50000", "cols >= 5000 and cols < 50000"),
    ("cols>=50000", "cols >= 50000"),
    ("Has genInt vars (general integer)", "genInt > 0"),
    ("No genInt vars (bin+cont only)", "genInt == 0"),
]


def fmt_hms(seconds):
    if pd.isna(seconds):
        return "  n/a"
    total = int(round(seconds))
    h, rem = divmod(total, 3600)
    m, _ = divmod(rem, 60)
    return f"{h:d}h{m:02d}m"


def load_features(path):
    if path.endswith(".tsv"):
        feat = pd.read_csv(path, sep="\t")
    else:
        feat = pd.read_csv(path)
    feat["Name"] = feat["Name"].astype(str).str.strip()
    feat = feat.drop_duplicates(subset="Name", keep="first").set_index("Name")
    feat["_sppc"] = (
        feat.get("rPercPartitioning", 0).fillna(0)
        + feat.get("rPercPacking", 0).fillna(0)
        + feat.get("rPercCovering", 0).fillna(0)
    )
    feat["_kp"] = feat.get("rPercKnapsack", 0).fillna(0) + feat.get(
        "rPercIntegerKnapsack", 0
    ).fillna(0)
    return feat


def load_experiments(dirs, labels, weights, bks_tsv):
    bks_lookup = cb.load_bks_tsv(bks_tsv) if bks_tsv and os.path.exists(bks_tsv) else {}
    plo, phi = weights["primal_gap_clip"]
    dlo, dhi = weights["dual_gap_clip"]
    dfs = {}
    for label, d in zip(labels, dirs):
        df, _kind = cb.detect_and_load(d, bks_lookup)
        df["instance"] = df["instance"].astype(str).str.strip()
        df = cb.enrich(df, weights)
        df["pgap_clip"] = df["primal_gap_bks"].clip(plo, phi)
        df["dgap_clip"] = df["dual_gap_bks"].clip(dlo, dhi)
        dfs[label] = df.set_index("instance")
    return dfs


def build_joined(dfs, labels, feat):
    """Inner-join all experiments on instance, then left-join features."""
    joined = None
    keep_cols = [
        "category", "cost", "pgap_clip", "dgap_clip", "elapsed_s",
    ]
    for label in labels:
        sub = dfs[label][keep_cols].add_suffix(f"__{label}")
        joined = sub if joined is None else joined.join(sub, how="inner")
    joined = joined.join(feat, how="left")
    return joined


def summarize_group(joined, labels, name, query):
    sub = joined if query is None else joined.query(query)
    n = len(sub)
    row = {"name": name, "n": n}
    if n == 0:
        for label in labels:
            row[label] = None
        return row
    for label in labels:
        cat = sub[f"category__{label}"]
        n_solved = int(((cat == "OPTIMAL") | (cat == "INFEASIBLE_CONFIRMED")).sum())
        row[label] = {
            "n_solved": n_solved,
            "avg_cost": sub[f"cost__{label}"].mean(),
            "avg_pgap": sub[f"pgap_clip__{label}"].mean(skipna=True),
            "avg_dgap": sub[f"dgap_clip__{label}"].mean(skipna=True),
            "avg_t": sub[f"elapsed_s__{label}"].mean(skipna=True),
        }
    return row


def print_report(rows, labels, feature_coverage, out):
    out.write(f"=== analyze_by_features: {len(labels)} experiment(s), "
              f"feature coverage {feature_coverage} ===\n\n")
    header_cols = []
    for label in labels:
        header_cols.append(f"{label}: solved / cost / pgap / dgap / avg_t")
    out.write(f"{'group':42s} {'n':>4s}  " + "   |   ".join(header_cols) + "\n")
    out.write("-" * (42 + 6 + len(labels) * 46) + "\n")
    for row in rows:
        name, n = row["name"], row["n"]
        parts = []
        for label in labels:
            v = row[label]
            if v is None:
                parts.append(f"{'n/a':>44s}")
                continue
            parts.append(
                f"{v['n_solved']:3d} / {v['avg_cost']:6.1f} / "
                f"{v['avg_pgap']:5.1f}% / {v['avg_dgap']:5.1f}% / {fmt_hms(v['avg_t'])}"
            )
        out.write(f"{name:42s} {n:4d}  " + "   |   ".join(parts) + "\n")


def parse_args():
    p = argparse.ArgumentParser(
        description="Break down a Cbc benchmark comparison by instance/matrix features.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("dirs", nargs="*", help="experiment dirs, optionally label=dir")
    p.add_argument("--list", help="file with one dir (or label=dir) per line")
    p.add_argument("--weights", default=cb.DEFAULT_WEIGHTS_PATH, help="weights JSON file")
    p.add_argument("--bks-tsv", default=cb.DEFAULT_BKS_TSV, help="bks.tsv fallback")
    p.add_argument("--features-file", default=None,
                    help="instance features CSV/TSV (auto-detected if omitted)")
    p.add_argument("--groups", default=None,
                    help="JSON file overriding built-in group definitions "
                         "(list of [name, query] pairs; query may be null for 'all')")
    p.add_argument("--no-color", action="store_true", help="(accepted, currently no-op)")
    p.add_argument("-o", "--output", help="write report to file instead of stdout")
    return p.parse_args()


def main():
    args = parse_args()
    dirs, labels = cb.collect_dirs_and_labels(args)
    if len(dirs) < 1:
        print("error: need at least one experiment directory", file=sys.stderr)
        sys.exit(2)

    weights = cb.load_weights(args.weights)

    features_file = args.features_file
    if features_file is None:
        for cand in DEFAULT_FEATURES_CANDIDATES:
            if os.path.exists(cand):
                features_file = cand
                break
    if features_file is None:
        print("error: no features file found/given; pass --features-file", file=sys.stderr)
        sys.exit(2)
    feat = load_features(features_file)

    dfs = load_experiments(dirs, labels, weights, args.bks_tsv)
    joined = build_joined(dfs, labels, feat)

    total_n = len(joined)
    with_feat = int(joined["cols"].notna().sum()) if "cols" in joined.columns else 0
    feature_coverage = f"{with_feat}/{total_n} instances"

    if args.groups:
        with open(args.groups) as fh:
            groups = [tuple(g) for g in json.load(fh)]
    else:
        groups = BUILTIN_GROUPS

    rows = [summarize_group(joined, labels, name, query) for name, query in groups]

    out = open(args.output, "w") if args.output else sys.stdout
    try:
        print_report(rows, labels, feature_coverage, out)
    finally:
        if args.output:
            out.close()


if __name__ == "__main__":
    main()
