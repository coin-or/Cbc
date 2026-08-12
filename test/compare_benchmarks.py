#!/usr/bin/env python3
"""
compare_benchmarks.py — Compare N>=2 Cbc benchmark experiment directories.

A from-scratch port of mipster's compare_multi_experiments.py data model
(load/classify/gap/cost/rank) to cbc-workspace, with three changes driven by
the specific failure modes that make this kind of comparison easy to get
subtly wrong:

  1. Costs/weights are externalized to a small JSON file (see
     compare_benchmarks_weights.json next to this script) instead of being
     hardcoded, so you can ask "how sensitive is this ranking to how harshly
     we penalize 'no solution found'?" without editing Python.
  2. Status classification is deliberately conservative about what counts as
     "found an integer-feasible solution". Cbc still writes an "Objective
     value:" line (and a .sol file) for a *fractional* LP relaxation when no
     integer-feasible incumbent exists ("Stopped ... (no integer solution -
     continuous used)"); a naive `objective != NA` check would silently
     treat that as a real solution. This script trusts the harness's own
     solution_found/proven_infeasible/timed_out flags when present (which
     are derived from Cbc's actual status text, not just "was a number
     printed"), and cross-checks them against status text/objective
     presence, warning loudly on any inconsistency (see --verify).
  3. Primal gap (objective vs BKS) and dual gap (dual bound vs BKS) are
     tracked as two separate numbers, not folded into a single "gap" -- a
     run that stopped with a great primal solution but a weak bound (or vice
     versa) looks very different under the two metrics, and conflating them
     hides that.

Supports two input layouts, auto-detected per directory from the header of
whichever file is present:

  * summary.tsv  (mipster-style "full benchmark run" output; self-contained,
    carries its own BKS in an 'expected' column) — this is the format used by
    /home/haroldo/experiments/cbc/<run>/summary.tsv.
  * results.tsv  (cbc-workspace's own ./test / run-mip-sanity-tests output;
    no embedded BKS -- looked up from Cbc/test/mip-sanity-data/bks.tsv by
    instance name instead).

Usage:
    compare_benchmarks.py <dir1> <dir2> [... <dirN>] [options]
    compare_benchmarks.py trust5=<dir1> trust10=<dir2> [options]
    compare_benchmarks.py --list experiments.txt [options]

experiments.txt: one directory per line; optionally "label=<dir>".

Options:
    --weights FILE       JSON weights/cost file (default: compare_benchmarks_weights.json
                          next to this script)
    --bks-tsv FILE       bks.tsv to use for results.tsv-style dirs lacking an
                          embedded BKS (default: Cbc/test/mip-sanity-data/bks.tsv)
    --top-n N            Instances shown in the biggest-divergence spotlight [default: 20]
    --gap-tol PCT        Minimum |gap delta| (pp) to flag an instance as changed [default: 0.01]
    --verify             Cross-check solution_found/proven_infeasible/timed_out
                          against status text + objective presence; print any
                          inconsistencies found (data-quality check on the
                          input files themselves, not on Cbc)
    --no-color           Disable ANSI colors
    -o/--output FILE     Write the text report to FILE instead of stdout
"""

import sys
import os
import re
import json
import argparse
import math

import pandas as pd
import numpy as np

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
DEFAULT_WEIGHTS_PATH = os.path.join(SCRIPT_DIR, "compare_benchmarks_weights.json")
DEFAULT_BKS_TSV = os.path.join(SCRIPT_DIR, "mip-sanity-data", "bks.tsv")

NA = float("nan")


# ─────────────────────────────────────────────────────────────────────────────
# Weights
# ─────────────────────────────────────────────────────────────────────────────

DEFAULT_WEIGHTS = {
    "cost_solved": 0.0,
    "gap_cap": 100.0,
    "cost_no_solution": 200.0,
    "cost_overtime": 300.0,
    "cost_error": 250.0,
    "cost_wrong": 1000.0,
    "cost_unknown": 200.0,
    "primal_gap_clip": [-50.0, 200.0],
    "dual_gap_clip": [-50.0, 200.0],
    "tie_break": "n_solved_desc",
}


def load_weights(path):
    weights = dict(DEFAULT_WEIGHTS)
    if path and os.path.exists(path):
        with open(path) as fh:
            raw = json.load(fh)
        for k, v in raw.items():
            if k.startswith("_"):
                continue
            weights[k] = v
    return weights


# ─────────────────────────────────────────────────────────────────────────────
# Status classification
# ─────────────────────────────────────────────────────────────────────────────
#
# Every status string we've observed in the wild (across both the mipster
# full-benchmark harness and cbc-workspace's own run-mip-sanity-tests)
# collapses into one of these buckets. "concluded" means the search produced
# a *certificate*: a proven optimum or a proven infeasibility -- not merely
# "stopped with something to show for it".

CATEGORIES = [
    "OPTIMAL",             # proven optimal
    "INFEASIBLE_CONFIRMED",  # proven infeasible, matches expectation
    "TIMEOUT_WITH_SOL",    # hit a limit, real integer-feasible incumbent exists
    "TIMEOUT_NO_SOL",      # hit a limit, no integer-feasible incumbent
    "OVERTIME",            # ignored its own limit, hard-killed by the harness
    "ERROR",               # crashed / ASAN error / no parsable result at all
    "WRONG",               # validator/harness caught an incorrect claim
    "UNKNOWN",             # unrecognized status text (fallback, should be rare)
]

_WRONG_MARKERS = ("WRONG", "INVALID_SOL")
_ERROR_MARKERS = ("ERROR", "CRASH")


def classify_status(status, solution_found, proven_infeasible, timed_out):
    """
    Classify a single run's raw status text into one of CATEGORIES, using the
    harness's own boolean flags (when available) rather than re-deriving
    "was a solution found" purely from whether an objective number is
    present -- Cbc happily prints an objective for a fractional relaxation
    too, so that check alone is not sufficient (see module docstring).
    """
    s = str(status).strip()

    if any(m in s for m in _WRONG_MARKERS):
        return "WRONG"
    if any(m in s for m in _ERROR_MARKERS):
        return "ERROR"

    if s.startswith("OVERTIME"):
        return "OVERTIME"

    if s.startswith("SOLVED"):
        if proven_infeasible == 1 or "(inf)" in s:
            return "INFEASIBLE_CONFIRMED"
        return "OPTIMAL"

    if s.startswith("TIMEOUT") or s.startswith("MEMLIMIT"):
        if solution_found == 1:
            return "TIMEOUT_WITH_SOL"
        return "TIMEOUT_NO_SOL"

    if s.startswith("NO_SOLUTION"):
        return "TIMEOUT_NO_SOL"

    return "UNKNOWN"


# ─────────────────────────────────────────────────────────────────────────────
# Loading: summary.tsv (mipster-style, self-contained BKS)
# ─────────────────────────────────────────────────────────────────────────────

def _to_num(series):
    return pd.to_numeric(series.replace("-", NA).infer_objects(copy=False), errors="coerce")


def _parse_gap_field(text):
    text = str(text).strip()
    if not text or text == "-" or text == "nan":
        return NA
    if text.startswith(">"):
        return 100.0
    if text.endswith("%"):
        text = text[:-1]
    try:
        return min(float(text), 100.0)
    except ValueError:
        return NA


def load_summary_tsv(path):
    df = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    df.columns = [c.strip() for c in df.columns]
    df["instance"] = df["instance"].str.strip()

    for col in ("objective", "dual_bound", "expected"):
        if col not in df.columns:
            df[col] = "-"
    for col in ("solution_found", "proven_infeasible", "timed_out"):
        if col not in df.columns:
            df[col] = NA

    out = pd.DataFrame({
        "instance":          df["instance"],
        "status":            df["status"].str.strip(),
        "objective":         _to_num(df["objective"]),
        "dual_bound":        _to_num(df["dual_bound"]),
        "bks":               _to_num(df["expected"]),
        "elapsed_s":         pd.to_numeric(df.get("elapsed_s", NA), errors="coerce"),
        "gap_field":         df.get("gap_field", "").apply(_parse_gap_field),
        "solution_found":    pd.to_numeric(df["solution_found"], errors="coerce"),
        "proven_infeasible": pd.to_numeric(df["proven_infeasible"], errors="coerce"),
        "timed_out":         pd.to_numeric(df["timed_out"], errors="coerce"),
    })
    return out


# ─────────────────────────────────────────────────────────────────────────────
# Loading: results.tsv (cbc-workspace ./test style, no embedded BKS)
# ─────────────────────────────────────────────────────────────────────────────

def load_results_tsv(path, bks_lookup):
    df = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    df.columns = [c.strip() for c in df.columns]
    df["instance"] = df["instance"].str.strip()

    obj = _to_num(df["obj"])
    bound = _to_num(df["bound"])
    gap_pct = _to_num(df["gap_pct"])
    is_optimal = pd.to_numeric(df["is_optimal"], errors="coerce").fillna(0)

    # results.tsv's "status" column is the *validator* outcome
    # (PASS/FAIL/OVERTIME/ERROR), not a Cbc status string -- remap it onto
    # the same vocabulary used by classify_status() so both loaders share
    # one classification path. Whether a genuine integer-feasible solution
    # was found is judged from obj being non-empty (results.tsv only ever
    # populates "obj" from Cbc's "Objective value:" line, which Cbc omits
    # entirely -- printing "No feasible solution found" instead -- whenever
    # there is no incumbent; see run-mip-sanity-tests' parse_cbc_log).
    solution_found = obj.notna().astype(int)

    remapped_status = []
    proven_infeasible = []
    for i in range(len(df)):
        v = df["status"].iloc[i].strip()
        opt = is_optimal.iloc[i] == 1
        has_sol = solution_found.iloc[i] == 1
        if v == "PASS":
            if opt:
                remapped_status.append("SOLVED")
                proven_infeasible.append(0)
            elif not has_sol and bound.isna().iloc[i] is False and pd.isna(obj.iloc[i]):
                # PASS with no solution and no claimed optimum only happens
                # for a confirmed proven-infeasible run.
                remapped_status.append("SOLVED(inf)")
                proven_infeasible.append(1)
            else:
                remapped_status.append("TIMEOUT(gap)" if has_sol else "TIMEOUT(no_sol)")
                proven_infeasible.append(0)
        elif v == "OVERTIME":
            remapped_status.append("OVERTIME")
            proven_infeasible.append(0)
        elif v == "ERROR":
            remapped_status.append("ERROR")
            proven_infeasible.append(0)
        elif v == "FAIL":
            remapped_status.append("WRONG")
            proven_infeasible.append(0)
        else:
            remapped_status.append(v)
            proven_infeasible.append(0)

    bks = df["instance"].map(bks_lookup) if bks_lookup else pd.Series([NA] * len(df))

    out = pd.DataFrame({
        "instance":          df["instance"],
        "status":            pd.Series(remapped_status),
        "objective":         obj,
        "dual_bound":        bound,
        "bks":               bks.reset_index(drop=True),
        "elapsed_s":         pd.to_numeric(df["elapsed_s"], errors="coerce"),
        "gap_field":         gap_pct,
        "solution_found":    solution_found,
        "proven_infeasible": pd.Series(proven_infeasible),
        "timed_out":         NA,
    })
    return out


def load_bks_tsv(path):
    """Returns dict instance -> signed BKS objective (NaN for infeasible/NA)."""
    lookup = {}
    if not os.path.exists(path):
        return lookup
    bdf = pd.read_csv(path, sep="\t", dtype=str)
    for _, row in bdf.iterrows():
        inst = row["instance"].strip()
        status = row.get("status", "").strip()
        obj = row.get("objective", "NA")
        if status == "infeasible" or obj in ("NA", "", None):
            lookup[inst] = NA
        else:
            try:
                lookup[inst] = float(obj)
            except ValueError:
                lookup[inst] = NA
    return lookup


def detect_and_load(exp_dir, bks_lookup):
    summary_path = os.path.join(exp_dir, "summary.tsv")
    results_path = os.path.join(exp_dir, "results.tsv")
    if os.path.exists(summary_path):
        return load_summary_tsv(summary_path), "summary.tsv"
    if os.path.exists(results_path):
        return load_results_tsv(results_path, bks_lookup), "results.tsv"
    raise FileNotFoundError(
        f"Neither summary.tsv nor results.tsv found in {exp_dir!r}"
    )


# ─────────────────────────────────────────────────────────────────────────────
# Gaps + cost
# ─────────────────────────────────────────────────────────────────────────────

def enrich(df, weights):
    df = df.copy()
    df["category"] = df.apply(
        lambda r: classify_status(
            r["status"], r["solution_found"], r["proven_infeasible"], r["timed_out"]
        ),
        axis=1,
    )
    df["concluded"] = df["category"].isin(["OPTIMAL", "INFEASIBLE_CONFIRMED"])
    df["has_solution"] = df["category"].isin(["OPTIMAL", "TIMEOUT_WITH_SOL"])

    def primal_gap(row):
        if row["category"] == "INFEASIBLE_CONFIRMED":
            return NA
        bks, obj = row["bks"], row["objective"]
        if pd.isna(bks) or pd.isna(obj) or not row["has_solution"]:
            return NA
        if abs(bks) <= 1e-10:
            return 0.0 if abs(obj) <= 1e-10 else NA
        return (obj - bks) / abs(bks) * 100.0

    def dual_gap(row):
        if row["category"] == "INFEASIBLE_CONFIRMED":
            return NA
        bks, dual = row["bks"], row["dual_bound"]
        if pd.isna(bks) or pd.isna(dual):
            return NA
        if abs(bks) <= 1e-10:
            return 0.0 if abs(dual) <= 1e-10 else NA
        return (bks - dual) / abs(bks) * 100.0

    df["primal_gap_bks"] = df.apply(primal_gap, axis=1)
    df["dual_gap_bks"] = df.apply(dual_gap, axis=1)

    gap_cap = weights["gap_cap"]

    def cost(row):
        cat = row["category"]
        if cat in ("OPTIMAL", "INFEASIBLE_CONFIRMED"):
            return weights["cost_solved"]
        if cat == "WRONG":
            return weights["cost_wrong"]
        if cat == "OVERTIME":
            return weights["cost_overtime"]
        if cat == "ERROR":
            return weights["cost_error"]
        if cat == "TIMEOUT_NO_SOL":
            return weights["cost_no_solution"]
        if cat == "TIMEOUT_WITH_SOL":
            g = row["gap_field"]
            if pd.isna(g):
                # solution found but gap unreported/unknown -- try to derive
                # from the primal/dual bounds directly before falling back.
                bks_dual = row["dual_gap_bks"]
                if not pd.isna(bks_dual):
                    g = max(0.0, bks_dual)
                else:
                    g = gap_cap
            return min(abs(g), gap_cap)
        return weights["cost_unknown"]

    df["cost"] = df.apply(cost, axis=1)
    return df


# ─────────────────────────────────────────────────────────────────────────────
# Cross-check (--verify): does the harness's own bookkeeping add up?
# ─────────────────────────────────────────────────────────────────────────────

def verify_consistency(df, label):
    """
    Sanity-checks the loaded/normalized rows against each other, independent
    of any weighting choice:
      - objective present  <=>  solution_found == 1 (when the flag exists)
      - proven_infeasible == 1  =>  objective absent
      - category == INFEASIBLE_CONFIRMED  =>  bks says infeasible or unknown,
        never "optimal" (that combination is a validator bug, not a solver
        outcome -- cbc_validate_sol already hard-errors on it, so it should
        never reach here as INFEASIBLE_CONFIRMED)
    Returns a list of human-readable issue strings (empty if all clean).
    """
    issues = []
    for _, row in df.iterrows():
        inst = row["instance"]
        has_sf_flag = not pd.isna(row["solution_found"])
        obj_present = not pd.isna(row["objective"])
        if has_sf_flag:
            sf = row["solution_found"] == 1
            if sf and not obj_present:
                issues.append(f"[{label}] {inst}: solution_found=1 but objective is empty")
            if obj_present and not sf and row["category"] not in ("WRONG",):
                issues.append(f"[{label}] {inst}: objective={row['objective']} present but solution_found=0")
        if row["proven_infeasible"] == 1 and obj_present:
            issues.append(f"[{label}] {inst}: proven_infeasible=1 but objective={row['objective']} present")
        if row["category"] == "INFEASIBLE_CONFIRMED" and not pd.isna(row["bks"]):
            issues.append(
                f"[{label}] {inst}: classified INFEASIBLE_CONFIRMED but bks={row['bks']} "
                f"looks like a finite objective, not an infeasible instance"
            )
    return issues


# ─────────────────────────────────────────────────────────────────────────────
# Merge + ranking
# ─────────────────────────────────────────────────────────────────────────────

def merge_all(dfs, labels):
    merged = None
    for i, df in enumerate(dfs):
        cols = {
            "status": f"status_{i}", "category": f"category_{i}",
            "cost": f"cost_{i}", "elapsed_s": f"elapsed_s_{i}",
            "objective": f"objective_{i}", "dual_bound": f"dual_bound_{i}",
            "primal_gap_bks": f"primal_gap_bks_{i}", "dual_gap_bks": f"dual_gap_bks_{i}",
            "concluded": f"concluded_{i}",
        }
        sub = df[["instance", "bks"] + list(cols.keys())].rename(columns=cols)
        if merged is None:
            merged = sub
        else:
            merged = merged.merge(sub.drop(columns=["bks"]), on="instance", how="outer")
    return merged


def compute_rankings(dfs, labels, merged, weights):
    n = len(dfs)
    cost_cols = [f"cost_{i}" for i in range(n)]
    n_total = len(merged)

    # An experiment that never ran an instance (missing row after the outer
    # join) is penalized as if it had produced no solution at all -- it must
    # never look *better* than a real (bad) result just because it's absent.
    pen = merged[cost_cols].fillna(weights["cost_no_solution"])
    pen_best = pen.min(axis=1)

    lo, hi = weights["primal_gap_clip"]
    dlo, dhi = weights["dual_gap_clip"]

    records = []
    for i, df in enumerate(dfs):
        cost_pen = pen[f"cost_{i}"]
        # NOTE: n_solved counts *proven* outcomes (category, not cost) --
        # a TIMEOUT_WITH_SOL row can also incidentally have cost==0.0 (a
        # reported 0.00% gap without an optimality certificate), and that
        # must not be conflated with a confirmed-optimal/-infeasible result.
        n_solved = int(df["concluded"].sum())
        n_wrong = int((df["category"] == "WRONG").sum())
        n_error = int((df["category"] == "ERROR").sum())
        n_overtime = int((df["category"] == "OVERTIME").sum())
        n_no_sol = int((df["category"] == "TIMEOUT_NO_SOL").sum())
        n_timeout_sol = int((df["category"] == "TIMEOUT_WITH_SOL").sum())

        pgb = df["primal_gap_bks"].clip(lo, hi)
        dgb = df["dual_gap_bks"].clip(dlo, dhi)

        n_wins = int((cost_pen == pen_best).sum())
        tie_count = pen.eq(pen_best, axis=0).sum(axis=1)
        n_sole_wins = int(((cost_pen == pen_best) & (tie_count == 1)).sum())

        records.append({
            "label": labels[i],
            "n_instances": len(df),
            "n_total_union": n_total,
            "avg_cost": cost_pen.mean(),
            "total_cost": cost_pen.sum(),
            "n_solved": n_solved,
            "n_timeout_with_sol": n_timeout_sol,
            "n_no_solution": n_no_sol,
            "n_overtime": n_overtime,
            "n_error": n_error,
            "n_wrong": n_wrong,
            "avg_primal_gap_bks": pgb.mean(skipna=True),
            "avg_dual_gap_bks": dgb.mean(skipna=True),
            "n_wins": n_wins,
            "n_sole_wins": n_sole_wins,
        })

    rdf = pd.DataFrame(records)
    rdf = rdf.sort_values(["avg_cost", "n_solved"], ascending=[True, False]).reset_index(drop=True)
    rdf["rank"] = range(1, len(rdf) + 1)
    return rdf


# ─────────────────────────────────────────────────────────────────────────────
# CLI plumbing
# ─────────────────────────────────────────────────────────────────────────────

def collect_dirs_and_labels(args):
    entries = []
    if args.list:
        with open(args.list) as fh:
            for line in fh:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                entries.append(line)
    else:
        entries = args.dirs

    dirs, labels = [], []
    for e in entries:
        if "=" in e and not os.path.exists(e):
            label, d = e.split("=", 1)
        else:
            d = e
            label = os.path.basename(os.path.normpath(d))
        dirs.append(d)
        labels.append(label)
    return dirs, labels


def parse_args():
    p = argparse.ArgumentParser(
        description="Compare N>=2 Cbc benchmark experiment directories.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("dirs", nargs="*", help="experiment dirs, optionally label=dir")
    p.add_argument("--list", help="file with one dir (or label=dir) per line")
    p.add_argument("--weights", default=DEFAULT_WEIGHTS_PATH, help="weights JSON file")
    p.add_argument("--bks-tsv", default=DEFAULT_BKS_TSV, help="bks.tsv for results.tsv-style dirs")
    p.add_argument("--top-n", type=int, default=20)
    p.add_argument("--gap-tol", type=float, default=0.01)
    p.add_argument("--verify", action="store_true", help="cross-check input-file bookkeeping")
    p.add_argument("--no-color", action="store_true")
    p.add_argument("-o", "--output", help="write report to file instead of stdout")
    return p.parse_args()


class Colors:
    def __init__(self, enabled):
        self.BOLD = "\033[1m" if enabled else ""
        self.RESET = "\033[0m" if enabled else ""
        self.DIM = "\033[2m" if enabled else ""
        self.RED = "\033[31m" if enabled else ""
        self.GREEN = "\033[32m" if enabled else ""
        self.YELLOW = "\033[33m" if enabled else ""
        self.CYAN = "\033[36m" if enabled else ""


def fmt(v, nd=2, pct=False):
    if v is None or (isinstance(v, float) and math.isnan(v)):
        return "-"
    s = f"{v:.{nd}f}"
    return s + "%" if pct else s


def main():
    args = parse_args()
    out = []

    def emit(s=""):
        out.append(s)

    dirs, labels = collect_dirs_and_labels(args)
    if len(dirs) < 2:
        print("Error: need at least 2 experiment directories to compare.", file=sys.stderr)
        return 1

    C = Colors(enabled=(not args.no_color) and sys.stdout.isatty() and not args.output)
    weights = load_weights(args.weights)
    bks_lookup = load_bks_tsv(args.bks_tsv)

    dfs, formats = [], []
    for d, lbl in zip(dirs, labels):
        raw, fmt_name = detect_and_load(d, bks_lookup)
        df = enrich(raw, weights)
        dfs.append(df)
        formats.append(fmt_name)

    if args.verify:
        all_issues = []
        for df, lbl in zip(dfs, labels):
            all_issues.extend(verify_consistency(df, lbl))
        emit(f"{C.BOLD}=== --verify: input bookkeeping cross-check ==={C.RESET}")
        if all_issues:
            emit(f"{C.YELLOW}{len(all_issues)} inconsistenc(y/ies) found:{C.RESET}")
            for issue in all_issues:
                emit(f"  {C.YELLOW}⚠{C.RESET}  {issue}")
        else:
            emit(f"{C.GREEN}No inconsistencies found.{C.RESET}")
        emit()

    merged = merge_all(dfs, labels)
    rankings = compute_rankings(dfs, labels, merged, weights)

    emit(f"{C.BOLD}=== compare_benchmarks: {len(dirs)} experiment(s) ==={C.RESET}")
    for d, lbl, fmt_name, df in zip(dirs, labels, formats, dfs):
        emit(f"  {C.CYAN}{lbl}{C.RESET}: {d}  ({fmt_name}, {len(df)} instances)")
    emit(f"  Weights file: {args.weights}")
    emit()

    emit(f"{C.BOLD}=== Ranking (lower avg cost = better; weights: solved={weights['cost_solved']}, "
         f"no_sol={weights['cost_no_solution']}, overtime={weights['cost_overtime']}, "
         f"error={weights['cost_error']}, wrong={weights['cost_wrong']}) ==={C.RESET}")
    header = (f"  {'#':>2} {'label':<20} {'avg_cost':>9} {'solved':>7} {'to_sol':>7} "
              f"{'no_sol':>7} {'over':>5} {'err':>4} {'wrong':>6} {'avg_pgap':>9} {'avg_dgap':>9} {'wins':>5}")
    emit(header)
    emit("  " + "─" * (len(header) - 2))
    for _, r in rankings.iterrows():
        emit(
            f"  {r['rank']:>2} {r['label']:<20} {fmt(r['avg_cost']):>9} {r['n_solved']:>7} "
            f"{r['n_timeout_with_sol']:>7} {r['n_no_solution']:>7} {r['n_overtime']:>5} "
            f"{r['n_error']:>4} {r['n_wrong']:>6} {fmt(r['avg_primal_gap_bks'], pct=True):>9} "
            f"{fmt(r['avg_dual_gap_bks'], pct=True):>9} {r['n_wins']:>5}"
        )
    emit()
    best = rankings.iloc[0]
    emit(f"  {C.GREEN}{C.BOLD}Winner: {best['label']}{C.RESET} "
         f"(avg cost {fmt(best['avg_cost'])} vs {' / '.join(fmt(v) for v in rankings['avg_cost'][1:])})")
    emit()

    # ── Per-instance divergence spotlight (biggest cost gap, any pair) ──────
    n = len(dfs)
    cost_cols = [f"cost_{i}" for i in range(n)]
    spread = merged[cost_cols].fillna(weights["cost_no_solution"]).max(axis=1) - \
        merged[cost_cols].fillna(weights["cost_no_solution"]).min(axis=1)
    merged = merged.assign(_spread=spread)
    top = merged.sort_values("_spread", ascending=False).head(args.top_n)

    emit(f"{C.BOLD}=== Top {min(args.top_n, len(top))} instances by cost spread across experiments ==={C.RESET}")
    hdr2 = "  " + f"{'instance':<40}" + "".join(f"{lbl[:14]:>16}" for lbl in labels)
    emit(hdr2)
    for _, row in top.iterrows():
        line = f"  {row['instance'][:40]:<40}"
        for i in range(n):
            c = row.get(f"cost_{i}")
            cat = row.get(f"category_{i}", "-")
            cell = f"{fmt(c)}({cat[:3] if isinstance(cat, str) else '-'})"
            line += f"{cell:>16}"
        emit(line)
    emit()

    # ── Instance-set bookkeeping ─────────────────────────────────────────────
    only_in = {}
    for i, lbl in enumerate(labels):
        mask = merged[f"cost_{i}"].isna()
        if mask.any():
            only_in[lbl] = int((~mask).sum())
    if len(set(len(df) for df in dfs)) > 1:
        emit(f"{C.YELLOW}Note: instance sets differ across experiments "
             f"({', '.join(f'{lbl}: {len(df)}' for lbl, df in zip(labels, dfs))}); "
             f"missing instances penalized at cost_no_solution.{C.RESET}")
        emit()

    text = "\n".join(out)
    if args.output:
        with open(args.output, "w") as fh:
            fh.write(text + "\n")
    else:
        print(text)
    return 0


if __name__ == "__main__":
    sys.exit(main())
