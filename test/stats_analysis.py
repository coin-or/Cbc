#!/usr/bin/env python3
"""
stats_analysis.py — Detailed time/efficiency breakdown for ONE Cbc benchmark
experiment, from the per-instance stats CSV files produced by Cbc's own
`-csvStatistics <file> -writeStatistics` (see Cbc/src/CbcSolverStatistics.cpp).

This is a from-scratch port of mipster's stats_analysis.py to cbc-workspace,
adapted for upstream Cbc's own CSV schema, which differs from mipster's
custom fork in two important ways (see "Porting notes" below):

  1. Column names.  mipster's fork writes prefixed columns like
     `cut_<generator>_time/_cuts/_calls/_avgnz` and `heur_<name>_time/_sols/
     _execs`. Upstream Cbc's own -csvStatistics (CbcSolverStatistics::writeCsv)
     instead writes one fixed, canonical set of columns:
       Name,result,time,sys,elapsed,objective,continuous,lp_seconds,tightened,
       cut_time,nodes,iterations,rows,columns,processed_rows,processed_columns,
       cgraph_time,cgraph_density,clqstr_extended,clqstr_dominated,clqstr_time,
       coefstr_changed,coefstr_rows,coefstr_time,
       rowred_fixed,rowred_duplicate,rowred_parallel,rowred_time,
       <Gen>_cuts,<Gen>_time  (one pair per canonical cut generator name,
                               e.g. "Clique_cuts,Clique_time" -- no "_calls"
                               or "_avgnz" column for any generator), and
       runtime_options.
     This script's CUT_GENERATORS table and column lookups are rewritten to
     match that scheme directly (no prefix, no calls/avgnz).

  2. No per-heuristic data at all.  Upstream Cbc's -csvStatistics does not
     instrument individual heuristics (RINS/FeasibilityPump/diving/...) the
     way mipster's fork does -- there is no heur_*_time/_sols/_execs column
     to read. The "Heuristics" and "Wasted heuristic time" sections are kept
     (for forward-compatibility / in case a future Cbc build adds this
     instrumentation) but degrade gracefully to an explicit "not available in
     this Cbc build" message instead of mipster's generic "no data found".

  3. New sections not present in mipster's version, covering phases this
     Cbc fork tracks that mipster's stats.csv schema doesn't have:
     conflict-graph construction (cgraph_time/density), clique strengthening
     (clqstr_*), coefficient tightening (coefstr_*), and row reduction
     (rowred_*) -- all pre-root-LP preprocessing steps. These get their own
     "Presolve / tightening phases" section instead of being silently
     absorbed into "Other".

  4. No external "rich" dependency. mipster's version renders tables with
     the `rich` package; this port uses the same plain/ANSI-color table
     style already established by compare_benchmarks.py in this repo, so it
     needs nothing beyond pandas (already a dependency here).

Usage:
    stats_analysis.py --outdir /path/to/exp/
    stats_analysis.py --outdir /path/to/exp/ --top 20 --no-color
    stats_analysis.py --statsfile /path/to/stats.csv --top 15
    stats_analysis.py --outdir /path/to/exp/ --cut-breakdown
"""

import argparse
import math
import sys
from pathlib import Path

import pandas as pd

# ─────────────────────────────────────────────────────────────────────────────
# Canonical cut generators (matches CbcSolverStatistics.cpp's
# canonicalGeneratorNames(), which is the fixed column set Cbc itself writes)
# ─────────────────────────────────────────────────────────────────────────────

CUT_GENERATORS = [
    ("Probing",               "Probing"),
    ("Gomory",                "Gomory (root)"),
    ("GomoryL1",              "Gomory L1"),
    ("GomoryL2",              "Gomory L2"),
    ("Gomory(2)",             "Gomory (tree)"),
    ("Knapsack",              "Knapsack"),
    ("Clique",                "Clique (BK)"),
    ("OddWheel",               "OddWheel"),
    ("MixedIntegerRounding2", "MIR2"),
    ("FlowCover",             "FlowCover"),
    ("TwoMirCuts",            "TwoMIR"),
    ("TwoMirCutsL1",          "TwoMIR L1"),
    ("TwoMirCutsL2",          "TwoMIR L2"),
    ("LiftAndProject",        "LiftAndProject"),
    ("ZeroHalf",              "ZeroHalf"),
    ("Reduce-and-split",      "ReduceSplit"),
    ("Reduce-and-split(2)",   "ReduceSplit2"),
]

# Kept for forward-compatibility: no upstream Cbc build currently populates
# these columns (see module docstring point 2). If a future -csvStatistics
# version adds per-heuristic instrumentation matching this naming scheme,
# the Heuristics/Wasted-heuristics sections below will pick it up
# automatically instead of needing a rewrite.
HEURISTICS = [
    ("feasibility_pump", "FeasibilityPump"),
    ("feasibilityjump",  "FeasibilityJump"),
    ("rins",             "RINS"),
    ("rens",             "RENS"),
    ("rounding",         "Rounding"),
    ("greedy_cover",     "GreedyCover"),
    ("greedy_equality",  "GreedyEquality"),
    ("divecoefficient",  "DiveCoefficient"),
    ("diveguided",       "DiveGuided"),
    ("divevectorlength", "DiveVectorLength"),
    ("vnd",              "VND"),
    ("naive",            "Naive"),
]

NA = float("nan")


# ─────────────────────────────────────────────────────────────────────────────
# Small helpers (numeric extraction, formatting -- no external deps)
# ─────────────────────────────────────────────────────────────────────────────

def _safe_sum(df, col, default=0.0):
    if col in df.columns:
        return pd.to_numeric(df[col], errors="coerce").fillna(0.0).sum()
    return default


def _safe_col(df, col):
    if col in df.columns:
        return pd.to_numeric(df[col], errors="coerce").fillna(0.0)
    return pd.Series([0.0] * len(df), index=df.index)


def load_stats(path):
    df = pd.read_csv(path)
    df.columns = [c.strip() for c in df.columns]
    if df.columns[0] != "Name":
        df = df.rename(columns={df.columns[0]: "Name"})
    df = df.rename(columns={"Name": "instance"})
    return df


def fmt_time(seconds):
    if seconds != seconds:  # NaN
        return "-"
    if seconds < 0.001:
        return "0"
    if seconds >= 3600:
        return f"{seconds/3600:.2f}h"
    if seconds >= 60:
        return f"{seconds/60:.1f}min"
    return f"{seconds:.2f}s"


def fmt_pct(frac):
    return f"{frac * 100:.1f}%"


def bar(frac, width=20):
    frac = max(0.0, min(1.0, frac)) if frac == frac else 0.0
    filled = round(frac * width)
    return "█" * filled + "░" * (width - filled)


class Colors:
    def __init__(self, enabled):
        self.BOLD = "\033[1m" if enabled else ""
        self.RESET = "\033[0m" if enabled else ""
        self.DIM = "\033[2m" if enabled else ""
        self.RED = "\033[31m" if enabled else ""
        self.YELLOW = "\033[33m" if enabled else ""
        self.GREEN = "\033[32m" if enabled else ""
        self.CYAN = "\033[36m" if enabled else ""
        self.MAGENTA = "\033[35m" if enabled else ""
        self.BLUE = "\033[34m" if enabled else ""


class Report:
    """Accumulates plain-text lines (with optional ANSI color) for one run."""

    def __init__(self, colors):
        self.C = colors
        self.lines = []

    def emit(self, s=""):
        self.lines.append(s)

    def title(self, s):
        self.emit(f"{self.C.BOLD}{self.C.CYAN}{s}{self.C.RESET}")

    def table(self, title, headers, widths, rows, footer=None):
        """rows: list of tuples of pre-formatted strings, one per column."""
        self.title(title)
        header_line = "  " + "".join(f"{h:>{w}}" if i else f"{h:<{w}}"
                                      for i, (h, w) in enumerate(zip(headers, widths)))
        self.emit(header_line)
        self.emit("  " + "─" * (sum(widths)))
        for row in rows:
            line = "  " + "".join(f"{c:>{w}}" if i else f"{c:<{w}}"
                                   for i, (c, w) in enumerate(zip(row, widths)))
            self.emit(line)
        if footer:
            self.emit("  " + "─" * (sum(widths)))
            self.emit("  " + "".join(f"{c:>{w}}" if i else f"{c:<{w}}"
                                      for i, (c, w) in enumerate(zip(footer, widths))))
        self.emit()

    def text(self, s):
        self.emit(s)


# ─────────────────────────────────────────────────────────────────────────────
# Section 1: Time budget breakdown
# ─────────────────────────────────────────────────────────────────────────────

def section_time_budget(df, rpt):
    total_time = _safe_sum(df, "time")
    lp_time = _safe_sum(df, "lp_seconds")
    cut_time = _safe_sum(df, "cut_time")
    cgraph_time = _safe_sum(df, "cgraph_time")
    presolve_time = (
        _safe_sum(df, "clqstr_time") + _safe_sum(df, "coefstr_time") + _safe_sum(df, "rowred_time")
    )

    accounted = lp_time + cut_time + cgraph_time + presolve_time
    other_time = max(total_time - accounted, 0.0)

    if total_time <= 0:
        rpt.text(f"{rpt.C.DIM}No timing data available.{rpt.C.RESET}")
        return

    phases = [
        ("LP solve", lp_time),
        ("Cut generation", cut_time),
        ("Conflict graph", cgraph_time),
        ("Presolve/tightening", presolve_time),
        ("Other (heuristics, branching, ...)", other_time),
    ]
    rows = []
    for name, t in phases:
        frac = t / total_time if total_time > 0 else 0.0
        rows.append((name, fmt_time(t), fmt_pct(frac), bar(frac)))
    rpt.table(
        "⏱  Time budget (summed across all instances)",
        ["Phase", "Total time", "Share", "Distribution"],
        [38, 12, 8, 22],
        rows,
        footer=("Total", fmt_time(total_time), "100%", ""),
    )


# ─────────────────────────────────────────────────────────────────────────────
# Section 2: Presolve / tightening phases (new vs mipster's schema)
# ─────────────────────────────────────────────────────────────────────────────

def section_presolve_tightening(df, rpt):
    have_any = any(
        c in df.columns
        for c in ("clqstr_time", "coefstr_time", "rowred_time", "cgraph_time")
    )
    if not have_any:
        return

    rows = []
    if "cgraph_time" in df.columns:
        t = _safe_sum(df, "cgraph_time")
        dens = pd.to_numeric(df.get("cgraph_density", pd.Series(dtype=float)), errors="coerce").mean()
        rows.append(("Conflict graph build", fmt_time(t), "-", "-",
                      f"avg density {dens:.4f}" if dens == dens else "-"))
    if "clqstr_time" in df.columns:
        t = _safe_sum(df, "clqstr_time")
        ext = int(_safe_sum(df, "clqstr_extended"))
        dom = int(_safe_sum(df, "clqstr_dominated"))
        rows.append(("Clique strengthening", fmt_time(t), f"{ext:,}", f"{dom:,}",
                      "constraints extended / dominated"))
    if "coefstr_time" in df.columns:
        t = _safe_sum(df, "coefstr_time")
        chg = int(_safe_sum(df, "coefstr_changed"))
        rws = int(_safe_sum(df, "coefstr_rows"))
        rows.append(("Coefficient tightening", fmt_time(t), f"{chg:,}", f"{rws:,}",
                      "coefficients reduced / rows touched"))
    if "rowred_time" in df.columns:
        t = _safe_sum(df, "rowred_time")
        fixed = int(_safe_sum(df, "rowred_fixed"))
        dup = int(_safe_sum(df, "rowred_duplicate"))
        par = int(_safe_sum(df, "rowred_parallel"))
        rows.append(("Row reduction", fmt_time(t), f"{fixed:,}", f"{dup:,}",
                      f"fixed / duplicate (+ {par:,} parallel)"))

    if not rows:
        return
    rpt.table(
        "🔧 Presolve / tightening phases (pre-root-LP)",
        ["Phase", "Total time", "Metric A", "Metric B", "Notes"],
        [26, 12, 10, 10, 40],
        rows,
    )


# ─────────────────────────────────────────────────────────────────────────────
# Section 3: Cut generator ranking
# ─────────────────────────────────────────────────────────────────────────────

def section_cuts(df, rpt, top_n):
    rows = []
    for gen, label in CUT_GENERATORS:
        time_col = f"{gen}_time"
        cuts_col = f"{gen}_cuts"
        t = _safe_sum(df, time_col)
        cuts = _safe_sum(df, cuts_col)
        if t > 0 or cuts > 0:
            rows.append({
                "label": label, "time": t, "cuts": int(cuts),
                "cuts_per_sec": cuts / t if t > 0 else 0.0,
            })

    if not rows:
        rpt.text(f"{rpt.C.DIM}No cut generator data found.{rpt.C.RESET}")
        return

    rows.sort(key=lambda r: r["time"], reverse=True)
    max_time = rows[0]["time"] if rows else 1.0

    out_rows = []
    for r in rows[:top_n]:
        frac = r["time"] / max_time if max_time > 0 else 0.0
        out_rows.append((
            r["label"], fmt_time(r["time"]), bar(frac),
            f"{r['cuts']:,}", f"{r['cuts_per_sec']:.2f}",
        ))
    rpt.table(
        f"✂  Cut generators — top {min(top_n, len(rows))} by time",
        ["Generator", "Total time", "Bar", "Total cuts", "Cuts/sec"],
        [22, 12, 22, 11, 10],
        out_rows,
    )


# ─────────────────────────────────────────────────────────────────────────────
# Section 4: Heuristic ranking (degrades gracefully -- see module docstring)
# ─────────────────────────────────────────────────────────────────────────────

def section_heuristics(df, rpt, top_n):
    rows = []
    for pfx, label in HEURISTICS:
        t = _safe_sum(df, f"heur_{pfx}_time")
        sols = _safe_sum(df, f"heur_{pfx}_sols")
        execs = _safe_sum(df, f"heur_{pfx}_execs")
        if t > 0 or sols > 0:
            rows.append({"label": label, "time": t, "sols": int(sols), "execs": int(execs)})

    if not rows:
        rpt.title("💡 Heuristics")
        rpt.text(
            f"  {rpt.C.DIM}Not available: this Cbc build's -csvStatistics does not "
            f"instrument individual heuristics (no heur_*_time/_sols/_execs columns). "
            f"Heuristic time is folded into the 'Other' bucket in the time budget above."
            f"{rpt.C.RESET}"
        )
        rpt.emit()
        return

    rows.sort(key=lambda r: r["time"], reverse=True)
    max_time = rows[0]["time"] if rows else 1.0
    out_rows = []
    for r in rows[:top_n]:
        frac = r["time"] / max_time if max_time > 0 else 0.0
        out_rows.append((r["label"], fmt_time(r["time"]), bar(frac),
                         f"{r['sols']:,}", f"{r['execs']:,}"))
    rpt.table(
        f"💡 Heuristics — top {min(top_n, len(rows))} by time",
        ["Heuristic", "Total time", "Bar", "Solutions", "Executions"],
        [22, 12, 22, 10, 11],
        out_rows,
    )


# ─────────────────────────────────────────────────────────────────────────────
# Section 5: Worst instances
# ─────────────────────────────────────────────────────────────────────────────

def section_worst_instances(df, rpt, top_n):
    if "time" not in df.columns:
        return
    df = df.copy()
    df["_time"] = _safe_col(df, "time")
    df["_lp_time"] = _safe_col(df, "lp_seconds")
    df["_cut_time"] = _safe_col(df, "cut_time")
    df["_cgraph_time"] = _safe_col(df, "cgraph_time")
    df["_presolve_time"] = _safe_col(df, "clqstr_time") + _safe_col(df, "coefstr_time") + _safe_col(df, "rowred_time")

    worst = df.nlargest(top_n, "_time")
    rows = []
    for _, r in worst.iterrows():
        nodes_val = r.get("nodes", "")
        try:
            nodes_str = f"{int(float(nodes_val)):,}"
        except (ValueError, TypeError):
            nodes_str = str(nodes_val)
        rows.append((
            str(r.get("instance", "?"))[:30], fmt_time(r["_time"]), fmt_time(r["_lp_time"]),
            fmt_time(r["_cut_time"]), fmt_time(r["_cgraph_time"]), fmt_time(r["_presolve_time"]), nodes_str,
        ))
    rpt.table(
        f"🐢 Top-{min(top_n, len(worst))} slowest instances",
        ["Instance", "Total", "LP", "Cuts", "CGraph", "Presolve", "Nodes"],
        [32, 9, 9, 9, 9, 10, 12],
        rows,
    )


# ─────────────────────────────────────────────────────────────────────────────
# Section 6: Wasted cut time (time spent, no cuts produced)
# ─────────────────────────────────────────────────────────────────────────────

def section_wasted_cuts(df, rpt, top_n):
    gen_rows = []
    for gen, label in CUT_GENERATORS:
        t = _safe_sum(df, f"{gen}_time")
        cuts = _safe_sum(df, f"{gen}_cuts")
        if t <= 0:
            continue
        # No "_calls" column exists upstream (see module docstring point 1),
        # so "ran but produced nothing" is judged from time>0 & cuts==0
        # directly, per-instance, instead of via a separate calls counter.
        zero_mask = (_safe_col(df, f"{gen}_time") > 0) & (_safe_col(df, f"{gen}_cuts") == 0)
        zero_cut_time = _safe_col(df, f"{gen}_time")[zero_mask].sum()
        zero_cut_insts = int(zero_mask.sum())
        time_per_cut = t / cuts if cuts > 0 else float("inf")
        gen_rows.append({
            "label": label, "gen": gen, "time": t, "cuts": int(cuts),
            "time_per_cut": time_per_cut, "zero_cut_insts": zero_cut_insts,
            "zero_cut_time": zero_cut_time,
        })

    if not gen_rows:
        return

    gen_rows.sort(key=lambda r: r["zero_cut_time"], reverse=True)
    max_wasted = gen_rows[0]["zero_cut_time"] if gen_rows else 1.0

    out_rows = []
    for r in gen_rows[:top_n]:
        wasted_frac = r["zero_cut_time"] / max_wasted if max_wasted > 0 else 0.0
        tpc = r["time_per_cut"]
        tpc_str = f"{tpc:.4f}" if tpc < float("inf") else "inf"
        out_rows.append((
            r["label"], fmt_time(r["zero_cut_time"]), bar(wasted_frac),
            f"{r['zero_cut_insts']:,}", fmt_time(r["time"]), tpc_str,
        ))
    rpt.table(
        "⚠  Cut generators: time spent where no cuts were produced",
        ["Generator", "Wasted time", "Bar", "# Insts(0 cuts)", "Total time", "Time/cut (s)"],
        [22, 13, 22, 16, 12, 13],
        out_rows,
    )

    spotlight_gens = [r for r in gen_rows if r["zero_cut_time"] > 0][:3]
    for gen in spotlight_gens:
        pfx, label = gen["gen"], gen["label"]
        t_col, c_col = f"{pfx}_time", f"{pfx}_cuts"
        mask = (_safe_col(df, t_col) > 0) & (_safe_col(df, c_col) == 0)
        sub = df[mask].copy()
        if sub.empty:
            continue
        sub["_wt"] = _safe_col(sub, t_col)
        worst = sub.nlargest(min(top_n, 8), "_wt")
        rows = []
        for _, r in worst.iterrows():
            total_t = pd.to_numeric(r.get("time", 0), errors="coerce")
            total_t = 0.0 if total_t != total_t else float(total_t)
            rows.append((str(r.get("instance", "?"))[:30], fmt_time(r["_wt"]), fmt_time(total_t)))
        rpt.table(
            f"  ↳ {label}: top instances with 0 cuts, most time spent",
            ["Instance", "Time (wasted)", "Total time"],
            [32, 15, 12],
            rows,
        )


# ─────────────────────────────────────────────────────────────────────────────
# Section 7: Wasted heuristic time (degrades gracefully; see docstring)
# ─────────────────────────────────────────────────────────────────────────────

def section_wasted_heuristics(df, rpt, top_n):
    heur_rows = []
    for pfx, label in HEURISTICS:
        t = _safe_sum(df, f"heur_{pfx}_time")
        execs = _safe_sum(df, f"heur_{pfx}_execs")
        if t <= 0 or execs <= 0:
            continue
        zero_sol_mask = (_safe_col(df, f"heur_{pfx}_execs") > 0) & (_safe_col(df, f"heur_{pfx}_sols") == 0)
        zero_sol_time = _safe_col(df, f"heur_{pfx}_time")[zero_sol_mask].sum()
        heur_rows.append({"label": label, "pfx": pfx, "time": t, "zero_sol_time": zero_sol_time,
                           "zero_sol_insts": int(zero_sol_mask.sum())})

    if not heur_rows:
        return  # section_heuristics() already explained why, no need to repeat

    heur_rows.sort(key=lambda r: r["zero_sol_time"], reverse=True)
    max_wasted = heur_rows[0]["zero_sol_time"] if heur_rows else 1.0
    out_rows = []
    for r in heur_rows[:top_n]:
        wasted_frac = r["zero_sol_time"] / max_wasted if max_wasted > 0 else 0.0
        out_rows.append((r["label"], fmt_time(r["zero_sol_time"]), bar(wasted_frac),
                         f"{r['zero_sol_insts']:,}", fmt_time(r["time"])))
    rpt.table(
        "⚠  Heuristics: time spent where no solution was found",
        ["Heuristic", "Wasted time", "Bar", "# Insts(0 sols)", "Total time"],
        [22, 13, 22, 16, 12],
        out_rows,
    )


# ─────────────────────────────────────────────────────────────────────────────
# Section 8: Per-instance cut breakdown (optional, verbose: --cut-breakdown)
# ─────────────────────────────────────────────────────────────────────────────

def section_cut_breakdown_by_instance(df, rpt, top_n):
    if "cut_time" not in df.columns:
        return
    df = df.copy()
    df["_cut_time"] = _safe_col(df, "cut_time")
    worst = df.nlargest(top_n, "_cut_time")

    active_gens = [(gen, label) for gen, label in CUT_GENERATORS if _safe_sum(df, f"{gen}_time") > 0]
    if not active_gens:
        return

    headers = ["Instance", "Total"] + [label[:14] for _, label in active_gens]
    widths = [24, 9] + [14] * len(active_gens)
    rows = []
    for _, r in worst.iterrows():
        total_ct = r["_cut_time"]
        if total_ct <= 0:
            break
        row = [str(r.get("instance", "?"))[:22], fmt_time(total_ct)]
        for gen, _label in active_gens:
            t = pd.to_numeric(r.get(f"{gen}_time", 0), errors="coerce")
            t = 0.0 if t != t else float(t)
            pct = t / total_ct if total_ct > 0 else 0.0
            row.append(f"{fmt_time(t)}({pct*100:.0f}%)" if t > 0 else "-")
        rows.append(tuple(row))

    rpt.table(f"✂  Cut time breakdown — top-{min(top_n, len(worst))} instances by cut time",
              headers, widths, rows)


# ─────────────────────────────────────────────────────────────────────────────
# Main
# ─────────────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Detailed time/efficiency analysis from a single Cbc benchmark "
                    "experiment's stats CSV files (Cbc's own -csvStatistics/-writeStatistics).",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    grp = parser.add_mutually_exclusive_group(required=True)
    grp.add_argument("--outdir", help="Experiment output directory (looks for stats.csv, or combines *.stats.csv)")
    grp.add_argument("--statsfile", help="Direct path to a stats.csv file")
    parser.add_argument("--top", type=int, default=15, metavar="N",
                        help="Number of top entries to show per table [default: 15]")
    parser.add_argument("--no-color", action="store_true", help="Suppress ANSI colors")
    parser.add_argument("--cut-breakdown", action="store_true",
                        help="Also show per-instance cut-type breakdown table")
    parser.add_argument("-o", "--output", help="Write the plain-text report to FILE "
                        "(default: <outdir>/stats_analysis.txt when --outdir is given, else stdout)")
    args = parser.parse_args()

    if args.statsfile:
        stats_path = Path(args.statsfile)
    else:
        outdir = Path(args.outdir)
        stats_path = outdir / "stats.csv"
        if not stats_path.exists():
            parts = sorted(outdir.glob("*.stats.csv"))
            if not parts:
                print(f"Error: no stats.csv or *.stats.csv files found in {outdir}", file=sys.stderr)
                return 1
            combined = pd.concat([pd.read_csv(p) for p in parts], ignore_index=True)
            stats_path = outdir / "stats.csv"
            combined.to_csv(stats_path, index=False)
            print(f"Combined {len(parts)} stats files -> {stats_path}")

    if not stats_path.exists():
        print(f"Error: {stats_path} not found", file=sys.stderr)
        return 1

    df = load_stats(stats_path)
    n = len(df)

    out_to_file = bool(args.output) or bool(args.outdir)
    C = Colors(enabled=not args.no_color and sys.stdout.isatty() and not (args.output))
    rpt = Report(C)

    label = Path(args.outdir).name if args.outdir else stats_path.name
    rpt.emit(f"{C.BOLD}=== Stats analysis: {label} ==={C.RESET}")
    rpt.emit(f"  {n} instance(s)   {stats_path}")
    rpt.emit()

    section_time_budget(df, rpt)
    section_presolve_tightening(df, rpt)
    section_cuts(df, rpt, args.top)
    section_heuristics(df, rpt, args.top)
    section_wasted_cuts(df, rpt, args.top)
    section_wasted_heuristics(df, rpt, args.top)
    section_worst_instances(df, rpt, args.top)
    if args.cut_breakdown:
        section_cut_breakdown_by_instance(df, rpt, args.top)

    text = "\n".join(rpt.lines)
    out_path = args.output
    if not out_path and args.outdir:
        out_path = str(Path(args.outdir) / "stats_analysis.txt")

    if out_path:
        with open(out_path, "w") as fh:
            fh.write(text + "\n")
        print(text)
        print(f"\nAnalysis saved to: {out_path}")
    else:
        print(text)
    return 0


if __name__ == "__main__":
    sys.exit(main())
