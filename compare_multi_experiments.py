#!/usr/bin/env python3
"""
compare_multi_experiments.py — Compare N≥2 CBC experiment directories.

Usage:
    python3 compare_multi_experiments.py <dir1> <dir2> ... [options]
    python3 compare_multi_experiments.py --list exp_list.txt [options]

exp_list.txt: one directory per line; optionally "label=<dir>" to set a short name.

Options:
    -o / --output FILE      Output PDF  [default: compare_multi_TIMESTAMP.pdf]
    -l / --labels L1 L2…    Short labels for each experiment (must match count)
    --list FILE             Read experiment dirs from a file
    --top-n N               Number of best/worst instances per spotlight [default: 20]

Cost system (lower is better):
    Solved to optimality / proven infeasible  →    0  pts
    Timeout with feasible solution            →  gap%  (0–100)
    No feasible solution                      →  200  pts  (sentinel)
    Overtime (killed past wallclock limit)    →  300  pts  (sentinel)

Ranking uses total cost as the primary metric; ties broken by # solved (descending).
"""

import sys
import os
import re
import argparse
import textwrap
import datetime

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.gridspec as gridspec
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.ticker import MaxNLocator
from matplotlib.lines import Line2D
from matplotlib.colors import to_rgba

# ─────────────────────────────────────────────────────────────────────────────
# Palette & style
# ─────────────────────────────────────────────────────────────────────────────

PALETTE = {
    "solved":   "#2ecc71",
    "timeout":  "#e67e22",
    "overtime": "#e74c3c",
    "infeas":   "#9b59b6",
    "other":    "#95a5a6",
    "blue":     "#2980b9",
    "dark":     "#2c3e50",
    "light":    "#ecf0f1",
    "accent":   "#f39c12",
    "teal":     "#1abc9c",
}

# Up to 10 distinct experiment colours from matplotlib's tab10 colourmap
_EXP_CMAP = plt.get_cmap("tab10")

def exp_color(i):
    """Return a distinct colour for experiment index i."""
    return _EXP_CMAP(i % 10)


plt.rcParams.update({
    "font.family":        "DejaVu Sans",
    "font.size":          9,
    "axes.spines.top":    False,
    "axes.spines.right":  False,
    "axes.grid":          True,
    "grid.alpha":         0.3,
    "axes.labelsize":     9,
    "axes.titlesize":     10,
    "axes.titleweight":   "bold",
    "figure.facecolor":   "white",
    "axes.facecolor":     "white",
})


# ─────────────────────────────────────────────────────────────────────────────
# Data loading  (copied from compare_experiments.py with minor adaptations)
# ─────────────────────────────────────────────────────────────────────────────

def load_summary(exp_dir):
    path = os.path.join(exp_dir, "summary.tsv")
    rows = []
    with open(path) as fh:
        header_line = next(fh, None)
        headers = [h.strip() for h in (header_line or "").rstrip("\n").split("\t")]
        has_dual = "dual_bound" in headers
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 6:
                continue
            if has_dual:
                rows.append({
                    "instance":   parts[0].strip(),
                    "status":     parts[1].strip(),
                    "objective":  parts[2].strip(),
                    "dual_bound": parts[3].strip() if len(parts) > 3 else "-",
                    "expected":   parts[4].strip() if len(parts) > 4 else "-",
                    "elapsed_s":  parts[5].strip() if len(parts) > 5 else "",
                    "threads":    parts[6].strip() if len(parts) > 6 else "",
                    "gap_field":  parts[7].strip() if len(parts) > 7 else "",
                })
            else:
                rows.append({
                    "instance":   parts[0].strip(),
                    "status":     parts[1].strip(),
                    "objective":  parts[2].strip(),
                    "dual_bound": "-",
                    "expected":   parts[3].strip() if len(parts) > 3 else "-",
                    "elapsed_s":  parts[4].strip() if len(parts) > 4 else "",
                    "threads":    parts[5].strip() if len(parts) > 5 else "",
                    "gap_field":  parts[6].strip() if len(parts) > 6 else "",
                })

    df = pd.DataFrame(rows, columns=[
        "instance", "status", "objective", "dual_bound", "expected",
        "elapsed_s", "threads", "gap_field",
    ])

    # Keep only instances with file-level evidence of execution:
    #   a .log file, a .result file, or a row in stats.csv
    files_in_dir = set(os.listdir(exp_dir))
    executed = (
        {f[:-4] for f in files_in_dir if f.endswith(".log")}
        | {f[:-7] for f in files_in_dir if f.endswith(".result")}
    )
    stats_path = os.path.join(exp_dir, "stats.csv")
    if os.path.exists(stats_path):
        try:
            sc = pd.read_csv(stats_path)
            if len(sc.columns) > 0:
                executed |= set(sc.iloc[:, 0].astype(str).str.strip())
        except Exception:
            pass
    df = df[df["instance"].isin(executed)].reset_index(drop=True)

    def parse_gap_field(text):
        text = str(text).strip()
        if not text or text == "-":
            return None
        if text.startswith(">"):
            return 100.0
        if text.endswith("%"):
            text = text[:-1]
        try:
            return min(float(text), 100.0)
        except ValueError:
            return None

    df["elapsed_s"]  = pd.to_numeric(df["elapsed_s"],  errors="coerce")
    df["objective"]  = pd.to_numeric(df["objective"],  errors="coerce")
    df["dual_bound"] = pd.to_numeric(df["dual_bound"], errors="coerce")
    df["expected"]   = pd.to_numeric(df["expected"],   errors="coerce")

    def categorise(s):
        s = s.strip()
        if s.startswith("SOLVED"):
            return "Solved"
        if "INFEASIBLE" in s and "WRONG" not in s:
            return "Solved"
        if s.startswith("TIMEOUT") and "no_sol" in s:
            return "Timeout (no sol)"
        if s.startswith("TIMEOUT"):
            return "Timeout (gap)"
        if s == "OVERTIME":
            return "Overtime"
        return "Other"

    df["category"] = df["status"].apply(categorise)

    def compute_gap(row):
        s = row["status"]
        if s.startswith("SOLVED"):
            return 0.0
        if "INFEASIBLE" in s and "WRONG" not in s:
            return 0.0
        if s == "OVERTIME":
            return 300.0
        if "WRONG" in s:
            return 200.0
        gf = parse_gap_field(row.get("gap_field", ""))
        if gf is not None:
            return gf
        m = re.search(r"gap=([0-9.]+)%", s)
        if m:
            return min(float(m.group(1)), 100.0)
        if "no_sol" in s or (row["objective"] != row["objective"]):
            return 200.0
        # Compute from primal/dual bounds: min(100, (obj-dual)/|dual|*100)
        dual = row["dual_bound"]
        obj  = row["objective"]
        if dual == dual and obj == obj:
            if abs(dual) > 1e-10:
                return min(abs(obj - dual) / abs(dual) * 100, 100.0)
            else:
                return 0.0 if abs(obj) <= 1e-10 else 100.0
        # Solution found but no gap info: conservative cap at 100
        if obj == obj:
            return 100.0
        return 200.0

    df["gap_pct"] = df.apply(compute_gap, axis=1)

    # ── Primal and dual bound gaps relative to BKS ───────────────────────────
    # primal_gap_bks = (obj − bks) / |bks| × 100  (min; <0 = CBC beat BKS)
    # dual_gap_bks   = (bks − dual) / |bks| × 100  (min; 0 = proven optimal)
    # Special case: BKS = 0 and obj = 0  →  gap = 0 (feasibility / cost-zero problems)
    def _is_infeasible(s):
        return ("INFEASIBLE" in s or "(inf)" in s) and "WRONG" not in s

    def primal_gap_bks(row):
        bks, obj = row["expected"], row["objective"]
        if bks != bks:                              # NaN: no BKS data
            return float("nan")
        if _is_infeasible(row["status"]):
            return float("nan")
        if obj != obj:                              # no solution found
            return float("nan")
        if abs(bks) <= 1e-10:                       # BKS = 0
            return 0.0 if abs(obj) <= 1e-10 else float("nan")
        return (obj - bks) / abs(bks) * 100

    def dual_gap_bks(row):
        bks, dual = row["expected"], row["dual_bound"]
        if bks != bks:
            return float("nan")
        if _is_infeasible(row["status"]):
            return float("nan")
        if dual != dual:
            return float("nan")
        if abs(bks) <= 1e-10:
            return 0.0 if abs(dual) <= 1e-10 else float("nan")
        return (bks - dual) / abs(bks) * 100

    df["primal_gap_bks"] = df.apply(primal_gap_bks, axis=1)
    df["dual_gap_bks"]   = df.apply(dual_gap_bks,   axis=1)

    # ── Unified score (used for all rankings and analysis) ───────────────────
    # cost = general CBC gap = min(100, (obj-dual)/|dual|*100)  if solution found
    # cost = 200                                                 if no solution / WRONG
    # cost = 300                                                 if OVERTIME (killed past wallclock)
    # cost = 0                                                   if proven optimal / infeasible
    def bks_cost(row):
        s = row["status"]
        if s == "OVERTIME":
            return 300.0
        if "WRONG" in s:
            return 200.0
        if _is_infeasible(s):
            return 0.0
        return row["gap_pct"]

    df["cost"] = df.apply(bks_cost, axis=1)
    return df


def load_meta(exp_dir):
    meta = {
        "binary": "n/a", "instances": "?", "timelimit": "?",
        "parallel": "?", "started": "?", "finished": "?",
        "cbc_opts": "", "git_commits": {},
    }
    path = os.path.join(exp_dir, "report.txt")
    if not os.path.exists(path):
        return meta
    with open(path) as fh:
        text = fh.read()
    for key, pat in [
        ("binary",    r"^Binary:\s+(.+)"),
        ("instances", r"^Instances:\s+(\d+)"),
        ("timelimit", r"^Timelimit:\s+(\S+)"),
        ("parallel",  r"^Parallel:\s+(\d+)"),
        ("started",   r"^Started:\s+(.+)"),
        ("finished",  r"^Finished:\s+(.+)"),
        ("cbc_opts",  r"^CBC opts:\s+(.+)"),
    ]:
        m = re.search(pat, text, re.MULTILINE)
        if m:
            meta[key] = m.group(1).strip()
    for proj in ("CoinUtils", "Osi", "Clp", "Cgl", "Cbc"):
        m = re.search(rf"{proj}\s+git:\s+([0-9a-f]{{7,}})", text, re.IGNORECASE)
        if m:
            meta["git_commits"][proj] = m.group(1)
    return meta


# ─────────────────────────────────────────────────────────────────────────────
# Merge helper
# ─────────────────────────────────────────────────────────────────────────────

def merge_all(dfs, labels):
    """
    Outer-join all per-experiment DataFrames on 'instance'.
    Each experiment contributes columns suffixed by its label index:
      cost_0, cost_1, ..., status_0, ..., elapsed_s_0, ...
      primal_gap_bks_0, ..., dual_gap_bks_0, ...

    Also adds:
      best_cost     — minimum cost across all experiments (virtual best)
      median_cost   — median cost across all experiments
    """
    merged = None
    for i, (df, lbl) in enumerate(zip(dfs, labels)):
        cols = {
            "cost":         f"cost_{i}",
            "status":        f"status_{i}",
            "category":      f"category_{i}",
            "elapsed_s":     f"elapsed_s_{i}",
            "objective":     f"objective_{i}",
            "dual_bound":    f"dual_bound_{i}",
            "gap_pct":       f"gap_pct_{i}",
            "primal_gap_bks": f"primal_gap_bks_{i}",
            "dual_gap_bks":   f"dual_gap_bks_{i}",
        }
        sub = df[["instance", "expected"] + list(cols.keys())].rename(columns=cols)
        if merged is None:
            merged = sub
        else:
            # expected should be the same across experiments for the same instance
            sub_noexp = sub.drop(columns=["expected"])
            merged = merged.merge(sub_noexp, on="instance", how="outer")

    cost_cols = [f"cost_{i}" for i in range(len(labels))]
    # best_cost ignores NaN (only counts experiments that ran the instance)
    merged["best_cost"]   = merged[cost_cols].min(axis=1)
    merged["median_cost"] = merged[cost_cols].median(axis=1)
    # penalized view: experiments that didn't run an instance get cost=200
    pen = merged[cost_cols].fillna(200.0)
    merged["best_cost_pen"] = pen.min(axis=1)
    return merged


# ─────────────────────────────────────────────────────────────────────────────
# Rankings computation
# ─────────────────────────────────────────────────────────────────────────────

def compute_rankings(dfs, labels, merged):
    """
    Build a DataFrame with one row per experiment, columns:
      label, n_inst_own, n_inst_total, avg_cost, total_cost,
      n_solved, n_with_sol, pct_solved,
      avg_gap, avg_primal_gap_bks, avg_dual_gap_bks,
      median_dual_gap_bks, n_wins, n_sole_wins, rank

    Instances an experiment did not run (crashed, missing) receive the maximum
    penalty cost of 200.  All metrics are computed over the full union of
    instances so that experiments with different instance sets are comparable.
    """
    n = len(dfs)
    cost_cols = [f"cost_{i}" for i in range(n)]
    n_total = len(merged)

    # Penalized cost matrix: NaN → 200 (crash / no run = no solution)
    pen = merged[cost_cols].fillna(200.0)
    pen_best = pen.min(axis=1)

    records = []
    for i, (df, lbl) in enumerate(zip(dfs, labels)):
        cost_pen = pen[f"cost_{i}"]

        total    = cost_pen.sum()
        avg_cost = cost_pen.mean()
        n_solved = int((cost_pen == 0).sum())
        n_sol    = int((cost_pen < 200).sum())

        # BKS-relative gaps: computed only over instances the experiment ran
        pgb = df["primal_gap_bks"].clip(-50, 200)
        dgb = df["dual_gap_bks"].clip(-50, 200)
        avg_pgb = pgb.mean(skipna=True)
        avg_dgb = dgb.mean(skipna=True)
        med_dgb = dgb.median(skipna=True)

        # Wins using penalized costs (missing run = cost 200, never beats a real result)
        n_wins = int((cost_pen == pen_best).sum())
        tie_ok = pen.eq(pen_best, axis=0).sum(axis=1)
        n_sole = int(((cost_pen == pen_best) & (tie_ok == 1)).sum())

        records.append({
            "label":               lbl,
            "n_inst_own":          len(df),
            "n_inst_total":        n_total,
            "avg_cost":            avg_cost,
            "total_cost":          total,
            "n_solved":            n_solved,
            "n_with_sol":          n_sol,
            "pct_solved":          n_solved / n_total * 100,
            "avg_gap":             avg_cost,
            "avg_gap_timeout":     avg_cost,   # kept for compat
            "median_gap_timeout":  float("nan"),
            "avg_primal_gap_bks":  avg_pgb,
            "avg_dual_gap_bks":    avg_dgb,
            "median_dual_gap_bks": med_dgb,
            "n_wins":              n_wins,
            "n_sole_wins":         n_sole,
        })

    rdf = pd.DataFrame(records)
    # Rank: lowest avg_cost = rank 1; break ties by n_solved descending
    rdf = rdf.sort_values(
        ["avg_cost", "n_solved"],
        ascending=[True, False]
    ).reset_index(drop=True)
    rdf["rank"] = range(1, len(rdf) + 1)
    return rdf


# ─────────────────────────────────────────────────────────────────────────────
# Layout helpers
# ─────────────────────────────────────────────────────────────────────────────

_INTRO_INSTANCE_SET = """\
INSTANCE SET
The benchmark contains 358 MIP instances from the combined set
~/inst/miplib/2017+spp, drawn from two sources:

  • MIPLIB 2017 (240 instances)
    The standard MIP benchmark library (miplib.zib.de).  Covers a broad
    range of real-world and academic problems across many application domains.

  • Additional set-structure instances (118 instances)
    Instances where a high proportion of constraints are of set-packing,
    set-partitioning, or set-covering type.  These problems are structurally
    amenable to conflict-graph techniques (clique cuts, odd-wheel cuts, BK
    clique separation) and are the primary target of recent CBC improvements.
    Drawn from the MIPLIB 2017 extended collection and additional public
    libraries (CORAL, MIPLIB 2010, MIPLIB 2003).

All instances are minimisation problems solved with a fixed time limit
(typically 300 s).  The best known solution (BKS) for each instance is
sourced from ~/inst/miplib/2017+spp/bks.tsv.
"""

_INTRO_METRICS = """\
METRICS & COST DEFINITION

  Cost  (primary ranking metric)
    ┌─────────────────────────────────────────────────────────────────┐
    │  gap   = min(100,  (obj − dual) / |dual| × 100 %)              │
    │  cost = gap                               if solution found    │
    │  cost = 200                               if no solution found │
    │  cost = 300                               if OVERTIME (killed) │
    └─────────────────────────────────────────────────────────────────┘
    0 means the run closed the primal–dual gap to zero (proven optimal).
    Values in (0, 100] represent the residual optimality gap at
    termination, capped at 100 %.  200 is a penalty for failing to
    find any feasible solution within the time limit.  300 is a higher
    penalty for OVERTIME runs that ignored the time limit and were
    hard-killed by the wallclock.

  Primal Gap BKS  =  (obj − BKS) / |BKS| × 100 %
    Distance between CBC's best incumbent and the best known solution.
    Negative values mean CBC found a better solution than the stored BKS
    (stale reference or a new record).

  Dual Gap BKS  =  (BKS − dual) / |BKS| × 100 %
    Distance between the dual lower bound and the BKS.
    0 % means the run proved the BKS optimal.
    Large values indicate the solver is still far from proving optimality.

  Gap (CBC-internal)  =  (obj − dual) / |dual| × 100 %
    Primal–dual gap as reported by CBC at termination.
    Measures proof-of-optimality progress within a single run but is
    not directly comparable across instances without a common reference.
"""


# ─────────────────────────────────────────────────────────────────────────────
# Page 0 — Introduction
# ─────────────────────────────────────────────────────────────────────────────

def page_intro(pdf, dirs, labels, dfs, metas):
    """
    Introductory page: experiment list, instance set description,
    and metric / cost definitions.
    """
    fig = plt.figure(figsize=(11, 8.5))
    fig.patch.set_facecolor("#fafafa")

    # ── Title block ──────────────────────────────────────────────────────────
    fig.text(0.5, 0.96,
             f"CBC Multi-Experiment Comparison  ({len(dirs)} experiments)",
             ha="center", va="top",
             fontsize=16, fontweight="bold", color=PALETTE["dark"])
    fig.text(0.5, 0.918,
             f"Generated by compare_multi_experiments.py",
             ha="center", va="top", fontsize=8, color="#7f8c8d", style="italic")

    # ── Thin rule ─────────────────────────────────────────────────────────────
    rule = fig.add_axes([0.05, 0.898, 0.90, 0.002])
    rule.set_facecolor(PALETTE["dark"])
    rule.axis("off")

    # ── Experiment table ──────────────────────────────────────────────────────
    ax_tbl = fig.add_axes([0.05, 0.715, 0.90, 0.175])
    ax_tbl.axis("off")
    col_hdrs = ["#", "Label", "Directory", "Instances", "Timelimit", "Binary"]
    rows = []
    cell_colors = []
    for i, (d, lbl, df, meta) in enumerate(zip(dirs, labels, dfs, metas)):
        n_inst = len(df)
        rows.append([
            str(i + 1), lbl, d,
            str(n_inst),
            meta.get("timelimit", "?"),
            os.path.basename(meta.get("binary", "n/a")),
        ])
        cell_colors.append(
            [f"#{int(255 - i*12):02x}{int(249 - i*8):02x}{int(196 + i*10):02x}"] +
            ["white"] * (len(col_hdrs) - 1)
        )

    tbl = ax_tbl.table(
        cellText=rows,
        colLabels=col_hdrs,
        cellColours=cell_colors,
        cellLoc="left",
        loc="upper center",
        bbox=[0.0, 0.0, 1.0, 1.0],
    )
    tbl.auto_set_font_size(False)
    tbl.set_fontsize(8)
    for (r, c), cell in tbl.get_celld().items():
        cell.set_edgecolor("#dddddd")
        if r == 0:
            cell.set_facecolor(PALETTE["dark"])
            cell.set_text_props(color="white", fontweight="bold")
        cell.set_linewidth(0.5)
        if c == 2:  # directory column — smaller font
            cell.set_text_props(fontsize=6.5)

    # ── Thin rule ─────────────────────────────────────────────────────────────
    rule2 = fig.add_axes([0.05, 0.708, 0.90, 0.001])
    rule2.set_facecolor("#cccccc")
    rule2.axis("off")

    # ── Two-column body ───────────────────────────────────────────────────────
    ax_left = fig.add_axes([0.05, 0.07, 0.43, 0.625])
    ax_left.axis("off")
    ax_left.add_patch(plt.Rectangle((0, 0), 1, 1,
                                     facecolor="#f0f4f8", edgecolor="#c5cfe0",
                                     linewidth=1, transform=ax_left.transAxes))
    ax_left.text(0.04, 0.97, _INTRO_INSTANCE_SET.strip(),
                 transform=ax_left.transAxes,
                 ha="left", va="top", fontsize=7.8, color=PALETTE["dark"],
                 fontfamily="monospace", linespacing=1.55)

    ax_right = fig.add_axes([0.52, 0.07, 0.44, 0.625])
    ax_right.axis("off")
    ax_right.add_patch(plt.Rectangle((0, 0), 1, 1,
                                      facecolor="#f5f5f0", edgecolor="#c5c5b0",
                                      linewidth=1, transform=ax_right.transAxes))
    ax_right.text(0.04, 0.97, _INTRO_METRICS.strip(),
                  transform=ax_right.transAxes,
                  ha="left", va="top", fontsize=7.8, color=PALETTE["dark"],
                  fontfamily="monospace", linespacing=1.55)

    # ── Footer ────────────────────────────────────────────────────────────────
    fig.text(0.5, 0.02,
             "Generated by compare_multi_experiments.py  |  CBC — COIN-OR Branch & Cut",
             ha="center", va="bottom", fontsize=7, color="#aaa")

    pdf.savefig(fig, bbox_inches="tight")
    plt.close(fig)



def page_title(fig, title, subtitle=""):
    fig.text(0.5, 0.97, title, ha="center", va="top",
             fontsize=14, fontweight="bold", color=PALETTE["dark"])
    if subtitle:
        fig.text(0.5, 0.942, subtitle, ha="center", va="top",
                 fontsize=8, color="#7f8c8d")


def exp_linestyle(i):
    """Return a distinct linestyle for experiment index i (cycles through 5 styles)."""
    styles = ["-", "--", ":", "-.", (0, (3, 1, 1, 1))]
    return styles[i % len(styles)]


def exp_legend_handles(labels):
    """Return Line2D list for a legend (shows both color and linestyle)."""
    return [Line2D([0], [0], color=exp_color(i), linewidth=1.8,
                   linestyle=exp_linestyle(i), label=lbl)
            for i, lbl in enumerate(labels)]


# ─────────────────────────────────────────────────────────────────────────────
# Page 1 — Executive Summary
# ─────────────────────────────────────────────────────────────────────────────

def page_executive_summary(pdf, dfs, labels, metas, rankings, merged):
    n_exp = len(labels)
    best_label = rankings.iloc[0]["label"]

    fig = plt.figure(figsize=(11, 8.5))
    page_title(fig, "Multi-Experiment Comparison — Executive Summary",
               subtitle=(
                   f"{n_exp} experiments  |  "
                   f"Overall winner (lowest avg cost per instance): {best_label}  |  "
                   f"Cost: 0=optimal  0–100=gap  200=no solution  300=overtime  (lower is better)"
               ))

    # ── Ranked metrics table ─────────────────────────────────────────────────
    ax_tbl = fig.add_axes([0.02, 0.47, 0.96, 0.44])
    ax_tbl.axis("off")

    col_labels = [
        "Rank", "Label", "Avg\nCost", "Own\nInst", "#Solved", "#With Sol",
        "%Solved", "Avg\nGap", "Avg Primal\nGap BKS", "Avg Dual\nGap BKS", "Wins", "Sole Wins",
    ]
    rows_data   = []
    cell_colors = []
    for _, row in rankings.iterrows():
        rows_data.append([
            f"#{row['rank']}",
            row["label"],
            f"{row['avg_cost']:.2f}",
            f"{row['n_inst_own']}/{row['n_inst_total']}",
            f"{row['n_solved']}",
            f"{row['n_with_sol']}",
            f"{row['pct_solved']:.1f}%",
            f"{row['avg_gap']:.1f}%" if not np.isnan(row["avg_gap"]) else "—",
            f"{row['avg_primal_gap_bks']:.1f}%" if not np.isnan(row["avg_primal_gap_bks"]) else "—",
            f"{row['avg_dual_gap_bks']:.1f}%" if not np.isnan(row["avg_dual_gap_bks"]) else "—",
            f"{row['n_wins']}",
            f"{row['n_sole_wins']}",
        ])
        base = "#fff9c4" if row["rank"] == 1 else "white"
        cell_colors.append([base] * len(col_labels))

    tbl = ax_tbl.table(
        cellText=rows_data,
        colLabels=col_labels,
        cellColours=cell_colors,
        cellLoc="center",
        loc="upper center",
        bbox=[0.0, 0.0, 1.0, 1.0],
    )
    tbl.auto_set_font_size(False)
    tbl.set_fontsize(8)
    for (r, c), cell in tbl.get_celld().items():
        cell.set_linewidth(0.4)
        if r == 0:
            cell.set_text_props(fontweight="bold")
            cell.set_facecolor("#d5d8dc")
        # Colour the label cell with experiment colour
        if r > 0 and c == 1:
            exp_idx = next(
                (j for j, lbl in enumerate(labels)
                 if lbl == rows_data[r - 1][1]), None)
            if exp_idx is not None:
                cell.set_facecolor(to_rgba(exp_color(exp_idx), alpha=0.25))

    ax_tbl.set_title("Rankings — all experiments", loc="left",
                     fontsize=9, fontweight="bold", pad=4)

    # ── Mini score-band stacked bar (one bar per experiment, ordered by rank) ─
    ax_bar = fig.add_axes([0.05, 0.07, 0.60, 0.30])
    bands = [
        ("Optimal",      lambda s: s == 0,               PALETTE["solved"]),
        ("Tight (0–10]", lambda s: (s > 0) & (s <= 10),  PALETTE["teal"]),
        ("Mod (10–50]",  lambda s: (s > 10) & (s <= 50), PALETTE["timeout"]),
        ("Large (50–100]", lambda s: (s > 50) & (s < 200),          PALETTE["accent"]),
        ("No sol",         lambda s: (s >= 200) & (s < 300),        PALETTE["overtime"]),
        ("Overtime",       lambda s: s >= 300,                       "#8e1a0e"),
    ]
    x = np.arange(n_exp)
    rank_order = list(rankings["label"])
    bottoms = np.zeros(n_exp)
    for _, mask_fn, col in bands:
        vals = []
        for lbl in rank_order:
            idx = labels.index(lbl)
            vals.append(mask_fn(dfs[idx]["cost"]).sum())
        vals = np.array(vals, dtype=float)
        ax_bar.bar(x, vals, bottom=bottoms, color=col, edgecolor="white",
                   linewidth=0.4)
        bottoms += vals

    ax_bar.set_xticks(x)
    ax_bar.set_xticklabels(rank_order, rotation=20, ha="right", fontsize=8)
    ax_bar.set_ylabel("# instances")
    ax_bar.set_title("Cost Band Breakdown  (ranked left→right)")
    ax_bar.legend(handles=[mpatches.Patch(color=col, label=nm)
                            for nm, _, col in bands],
                  fontsize=7, loc="upper right")

    # ── Meta summary text (rightmost panel) ──────────────────────────────────
    ax_meta = fig.add_axes([0.68, 0.07, 0.30, 0.30])
    ax_meta.axis("off")
    lines = []
    for lbl, meta in zip(labels, metas):
        cbc_git = meta["git_commits"].get("Cbc", "?")
        lines.append(f"[{lbl}]")
        lines.append(f"  Timelimit: {meta['timelimit']}s")
        lines.append(f"  Instances: {meta['instances']}")
        lines.append(f"  Cbc git:   {cbc_git[:8] if cbc_git != '?' else '?'}")
        lines.append("")
    ax_meta.text(0.0, 1.0, "\n".join(lines), va="top", ha="left",
                 fontsize=7, family="monospace", transform=ax_meta.transAxes,
                 color=PALETTE["dark"])
    ax_meta.set_title("Experiment Info", loc="left", fontsize=8,
                      fontweight="bold", pad=4)

    pdf.savefig(fig, bbox_inches="tight")
    plt.close(fig)


# ─────────────────────────────────────────────────────────────────────────────
# Page 2 — Score Profiles & Analysis
# ─────────────────────────────────────────────────────────────────────────────

def page_cost_profiles(pdf, dfs, labels, merged, rankings):
    n_exp = len(labels)
    cost_cols = [f"cost_{i}" for i in range(n_exp)]

    fig = plt.figure(figsize=(11, 8.5))
    page_title(fig, "Cost Profiles & Analysis",
               subtitle="Cost: 0=optimal  0–100=gap  200=no solution  300=overtime  (lower is better)")
    gs = gridspec.GridSpec(2, 2, figure=fig,
                           top=0.90, bottom=0.08, left=0.08, right=0.97,
                           hspace=0.46, wspace=0.36)

    # ── A: Cumulative score CDF ───────────────────────────────────────────────
    ax_a = fig.add_subplot(gs[0, 0])
    for i, (df, lbl) in enumerate(zip(dfs, labels)):
        srt = np.sort(df["cost"].values)
        ax_a.step(srt, np.arange(1, len(srt) + 1),
                  color=exp_color(i), linewidth=1.8,
                  linestyle=exp_linestyle(i),
                  where="post", label=lbl)
    ax_a.axvline(100, color="#ddd", linewidth=0.8, linestyle=":")
    ax_a.axvline(200, color="#ddd", linewidth=0.8, linestyle=":")
    ax_a.axvline(300, color="#ddd", linewidth=0.8, linestyle=":")
    ax_a.set_xlabel("Cost")
    ax_a.set_ylabel("# instances with cost ≤ x")
    ax_a.set_title("Cumulative Cost Profile\n(higher & left = better)")
    ax_a.legend(fontsize=7, framealpha=0.85)
    ax_a.set_xlim(-5, 315)

    # ── B: Score band breakdown grouped bar ───────────────────────────────────
    ax_b = fig.add_subplot(gs[0, 1])
    band_specs = [
        ("Optimal",         lambda s: s == 0),
        ("Tight\n(0–10]",   lambda s: (s > 0) & (s <= 10)),
        ("Mod\n(10–50]",    lambda s: (s > 10) & (s <= 50)),
        ("Large\n(50–100]", lambda s: (s > 50) & (s < 200)),
        ("No sol\n(200)",   lambda s: (s >= 200) & (s < 300)),
        ("Overtime\n(300)", lambda s: s >= 300),
    ]
    band_colors = [PALETTE["solved"], PALETTE["teal"], PALETTE["timeout"],
                   PALETTE["accent"], PALETTE["overtime"], "#8e1a0e"]
    x = np.arange(len(band_specs))
    w = 0.8 / n_exp
    for i, (df, lbl) in enumerate(zip(dfs, labels)):
        vals = [mask(df["cost"]).sum() for _, mask in band_specs]
        bars = ax_b.bar(x + (i - n_exp / 2 + 0.5) * w, vals, w,
                        color=exp_color(i), edgecolor="white", linewidth=0.4,
                        label=lbl)
    ax_b.set_xticks(x)
    ax_b.set_xticklabels([nm for nm, _ in band_specs], fontsize=7.5)
    ax_b.set_ylabel("# instances")
    ax_b.set_title("Cost Band Breakdown")
    ax_b.legend(fontsize=7, framealpha=0.85)

    # ── C: Pairwise score differences vs best-ranked experiment ──────────────
    ax_c = fig.add_subplot(gs[1, 0])
    # Use the top-ranked experiment as baseline (find its index in original labels list)
    best_label_c = rankings.iloc[0]["label"]
    best_idx = labels.index(best_label_c) if best_label_c in labels else 0
    # Sorted bar: for each other experiment, sort instances by Δcost = exp - best
    for i, (df, lbl) in enumerate(zip(dfs, labels)):
        if i == best_idx:
            continue
        base_col = f"cost_{best_idx}"
        exp_col  = f"cost_{i}"
        valid = merged[[base_col, exp_col]].dropna()
        delta = (valid[exp_col] - valid[base_col]).values
        delta_sorted = np.sort(delta)
        pos = np.arange(len(delta_sorted))
        ax_c.step(pos, delta_sorted, color=exp_color(i),
                  linewidth=1.4, where="post",
                  linestyle=exp_linestyle(i),
                  label=f"{lbl} − {labels[best_idx]}")
    ax_c.axhline(0, color=PALETTE["dark"], linewidth=0.9, linestyle="--")
    ax_c.set_xlabel("Instances (sorted by Δ cost)")
    ax_c.set_ylabel("Δ cost vs best-ranked experiment")
    ax_c.set_title(f"Cost Differences vs {labels[best_idx]}\n"
                   f"(below 0 = better than {labels[best_idx]})")
    ax_c.legend(fontsize=7, framealpha=0.85)
    ax_c.set_ylim(-210, 210)

    # ── D: Head-to-head win matrix (heatmap) ─────────────────────────────────
    ax_d = fig.add_subplot(gs[1, 1])
    # wins[i, j] = # instances where exp i has STRICTLY better score than exp j
    win_matrix = np.zeros((n_exp, n_exp), dtype=int)
    for i in range(n_exp):
        for j in range(n_exp):
            if i == j:
                win_matrix[i, j] = -1  # diagonal (unused)
                continue
            sc_i = merged[f"cost_{i}"]
            sc_j = merged[f"cost_{j}"]
            valid = sc_i.notna() & sc_j.notna()
            win_matrix[i, j] = int((sc_i[valid] < sc_j[valid]).sum())

    # Show as text in a table-style heatmap
    im = ax_d.imshow(win_matrix.astype(float), cmap="RdYlGn", aspect="auto",
                     vmin=0, vmax=merged.shape[0])
    for ii in range(n_exp):
        for jj in range(n_exp):
            if ii == jj:
                ax_d.text(jj, ii, "—", ha="center", va="center", fontsize=8,
                          color="#888")
            else:
                ax_d.text(jj, ii, str(win_matrix[ii, jj]),
                          ha="center", va="center", fontsize=7.5,
                          color="white" if win_matrix[ii, jj] > merged.shape[0] * 0.5 else "black")
    ax_d.set_xticks(range(n_exp))
    ax_d.set_yticks(range(n_exp))
    ax_d.set_xticklabels(labels, rotation=25, ha="right", fontsize=7.5)
    ax_d.set_yticklabels(labels, fontsize=7.5)
    ax_d.set_title("Head-to-Head Wins\n(row wins against column)")
    plt.colorbar(im, ax=ax_d, shrink=0.7, label="# wins")

    pdf.savefig(fig, bbox_inches="tight")
    plt.close(fig)


# ─────────────────────────────────────────────────────────────────────────────
# Page 3 — BKS Gap Analysis
# ─────────────────────────────────────────────────────────────────────────────

def page_bks_gaps(pdf, dfs, labels):
    CLIP_LO, CLIP_HI = -30.0, 200.0

    fig = plt.figure(figsize=(11, 8.5))
    page_title(fig, "Primal & Dual Bound Gaps vs Best-Known Solution (BKS)",
               subtitle=(
                   "primal gap = (obj − BKS)/|BKS|×100  |  "
                   "dual gap = (BKS − dual)/|BKS|×100  |  "
                   f"CDFs clipped to [{CLIP_LO:.0f}%, {CLIP_HI:.0f}%]"
               ))
    gs = gridspec.GridSpec(2, 2, figure=fig,
                           top=0.89, bottom=0.08, left=0.08, right=0.97,
                           hspace=0.44, wspace=0.34)

    # ── A: Primal gap BKS CDF (timeout instances with solution) ──────────────
    ax_a = fig.add_subplot(gs[0, 0])
    for i, (df, lbl) in enumerate(zip(dfs, labels)):
        raw = df.loc[
            df["primal_gap_bks"].notna() & (df["gap_pct"] > 0),
            "primal_gap_bks"
        ].values
        if len(raw) > 0:
            srt = np.sort(np.clip(raw, CLIP_LO, CLIP_HI))
            ax_a.step(srt, np.arange(1, len(srt) + 1),
                      color=exp_color(i), linewidth=1.8,
                      linestyle=exp_linestyle(i),
                      where="post", label=f"{lbl} (n={len(raw)})")
    ax_a.axvline(0, color="#bbb", linewidth=0.8, linestyle=":")
    ax_a.set_xlabel("Primal gap vs BKS (%)\n(negative = CBC beat BKS)")
    ax_a.set_ylabel("# instances")
    ax_a.set_title("Primal Gap CDF\n(timeout instances with solution)")
    ax_a.legend(fontsize=7, framealpha=0.85)

    # ── B: Dual gap BKS CDF (all non-infeasible) ─────────────────────────────
    ax_b = fig.add_subplot(gs[0, 1])
    for i, (df, lbl) in enumerate(zip(dfs, labels)):
        raw = df.loc[df["dual_gap_bks"].notna(), "dual_gap_bks"].values
        if len(raw) > 0:
            srt = np.sort(np.clip(raw, CLIP_LO, CLIP_HI))
            proven = (raw <= 1e-6).sum()
            ax_b.step(srt, np.arange(1, len(srt) + 1),
                      color=exp_color(i), linewidth=1.8,
                      linestyle=exp_linestyle(i),
                      where="post", label=f"{lbl} (n={len(raw)}, {proven} opt.)")
    ax_b.axvline(0, color="#bbb", linewidth=0.8, linestyle=":")
    ax_b.set_xlabel("Dual gap vs BKS (%)\n(0 = proven optimal)")
    ax_b.set_ylabel("# instances")
    ax_b.set_title("Dual Gap CDF\n(all non-infeasible instances)")
    ax_b.legend(fontsize=7, framealpha=0.85)

    # ── C: Summary table — avg/median BKS gaps ───────────────────────────────
    ax_c = fig.add_subplot(gs[1, 0])
    ax_c.axis("off")
    col_labels = ["Label", "n (primal)", "Med primal\ngap BKS",
                  "n (dual)", "Med dual\ngap BKS", "Proven\nopt."]
    rows_data = []
    for i, (df, lbl) in enumerate(zip(dfs, labels)):
        pgb = df.loc[df["primal_gap_bks"].notna() & (df["gap_pct"] > 0),
                     "primal_gap_bks"]
        dgb = df.loc[df["dual_gap_bks"].notna(), "dual_gap_bks"]
        med_pgb = np.nanmedian(np.clip(pgb.values, CLIP_LO, CLIP_HI)) if len(pgb) > 0 else float("nan")
        med_dgb = np.nanmedian(np.clip(dgb.values, CLIP_LO, CLIP_HI)) if len(dgb) > 0 else float("nan")
        proven = (df["dual_gap_bks"] <= 1e-6).sum()
        rows_data.append([
            lbl,
            str(len(pgb)),
            f"{med_pgb:.1f}%" if not np.isnan(med_pgb) else "—",
            str(len(dgb)),
            f"{med_dgb:.1f}%" if not np.isnan(med_dgb) else "—",
            str(proven),
        ])
    tbl = ax_c.table(cellText=rows_data, colLabels=col_labels,
                     cellLoc="center", loc="center",
                     bbox=[0.0, 0.0, 1.0, 1.0])
    tbl.auto_set_font_size(False)
    tbl.set_fontsize(8)
    for (r, cc), cell in tbl.get_celld().items():
        cell.set_linewidth(0.4)
        if r == 0:
            cell.set_text_props(fontweight="bold")
            cell.set_facecolor("#d5d8dc")
        elif cc == 0:
            exp_idx = next((j for j, lbl in enumerate(labels)
                            if lbl == rows_data[r - 1][0]), None)
            if exp_idx is not None:
                cell.set_facecolor(to_rgba(exp_color(exp_idx), alpha=0.25))
    ax_c.set_title("BKS Gap Summary", loc="left", fontsize=9,
                   fontweight="bold", pad=4)

    # ── D: % instances within gap thresholds (BKS proximity bar chart) ───────
    ax_d = fig.add_subplot(gs[1, 1])
    thresholds = [0.0, 1.0, 5.0, 10.0, 20.0]
    threshold_labels = ["=0%", "≤1%", "≤5%", "≤10%", "≤20%"]
    x = np.arange(len(thresholds))
    w = 0.8 / len(labels)
    for i, (df, lbl) in enumerate(zip(dfs, labels)):
        dgb = df["dual_gap_bks"].clip(CLIP_LO, CLIP_HI)
        total_valid = dgb.notna().sum()
        if total_valid == 0:
            continue
        pcts = []
        for thr in thresholds:
            if thr == 0.0:
                pcts.append((dgb <= 1e-6).sum() / total_valid * 100)
            else:
                pcts.append((dgb.notna() & (dgb <= thr)).sum() / total_valid * 100)
        ax_d.bar(x + (i - len(labels) / 2 + 0.5) * w, pcts, w,
                 color=exp_color(i), edgecolor="white", linewidth=0.4,
                 label=lbl)
    ax_d.set_xticks(x)
    ax_d.set_xticklabels(threshold_labels, fontsize=8)
    ax_d.set_ylabel("% instances")
    ax_d.set_title("% Instances Within Dual Gap Threshold\n(vs BKS, higher = better)")
    ax_d.legend(fontsize=7, framealpha=0.85)

    pdf.savefig(fig, bbox_inches="tight")
    plt.close(fig)


# ─────────────────────────────────────────────────────────────────────────────
# Page 4 — Category Heatmap
# ─────────────────────────────────────────────────────────────────────────────

def page_category_heatmap(pdf, dfs, labels, merged):
    """
    Show average cost per experiment × instance category.
    Also show % solved per experiment × category.
    """
    n_exp = len(labels)

    # Instance categories from 'category_0' (use first experiment's classification)
    cat_col = "category_0"
    if cat_col not in merged.columns:
        return
    categories = [c for c in merged[cat_col].dropna().unique()
                  if c not in ("", "Other")]
    categories = sorted(categories)
    if not categories:
        return

    fig = plt.figure(figsize=(11, 8.5))
    page_title(fig, "Performance by Instance Category",
               subtitle="Average cost per experiment × category  |  lower is better")
    gs = gridspec.GridSpec(1, 2, figure=fig,
                           top=0.89, bottom=0.10, left=0.08, right=0.97,
                           wspace=0.44)

    # ── Left: Avg score heatmap ───────────────────────────────────────────────
    ax_a = fig.add_subplot(gs[0, 0])
    cost_mat = np.full((n_exp, len(categories)), float("nan"))
    for i in range(n_exp):
        sc_col = f"cost_{i}"
        cat_col_i = f"category_{i}"
        col_to_use = cat_col_i if cat_col_i in merged.columns else cat_col
        for j, cat in enumerate(categories):
            mask = merged[col_to_use] == cat
            vals = merged.loc[mask, sc_col].dropna()
            if len(vals) > 0:
                cost_mat[i, j] = vals.mean()

    im_a = ax_a.imshow(cost_mat, cmap="RdYlGn_r", aspect="auto",
                       vmin=0, vmax=100)
    for ii in range(n_exp):
        for jj in range(len(categories)):
            v = cost_mat[ii, jj]
            if not np.isnan(v):
                ax_a.text(jj, ii, f"{v:.0f}",
                          ha="center", va="center", fontsize=7,
                          color="white" if v > 60 else "black")
    ax_a.set_xticks(range(len(categories)))
    ax_a.set_xticklabels(categories, rotation=30, ha="right", fontsize=7.5)
    ax_a.set_yticks(range(n_exp))
    ax_a.set_yticklabels(labels, fontsize=8)
    ax_a.set_title("Avg Cost by Category\n(red=high/bad, green=low/good)")
    plt.colorbar(im_a, ax=ax_a, shrink=0.7, label="avg cost")

    # ── Right: % solved heatmap ───────────────────────────────────────────────
    ax_b = fig.add_subplot(gs[0, 1])
    solved_mat = np.full((n_exp, len(categories)), float("nan"))
    for i in range(n_exp):
        sc_col = f"cost_{i}"
        cat_col_i = f"category_{i}"
        col_to_use = cat_col_i if cat_col_i in merged.columns else cat_col
        for j, cat in enumerate(categories):
            mask = merged[col_to_use] == cat
            vals = merged.loc[mask, sc_col].dropna()
            if len(vals) > 0:
                solved_mat[i, j] = (vals == 0).sum() / len(vals) * 100

    im_b = ax_b.imshow(solved_mat, cmap="RdYlGn", aspect="auto",
                       vmin=0, vmax=100)
    for ii in range(n_exp):
        for jj in range(len(categories)):
            v = solved_mat[ii, jj]
            if not np.isnan(v):
                ax_b.text(jj, ii, f"{v:.0f}%",
                          ha="center", va="center", fontsize=7,
                          color="white" if v < 30 else "black")
    ax_b.set_xticks(range(len(categories)))
    ax_b.set_xticklabels(categories, rotation=30, ha="right", fontsize=7.5)
    ax_b.set_yticks(range(n_exp))
    ax_b.set_yticklabels(labels, fontsize=8)
    ax_b.set_title("% Solved by Category\n(green=high/good, red=low/bad)")
    plt.colorbar(im_b, ax=ax_b, shrink=0.7, label="% solved")

    pdf.savefig(fig, bbox_inches="tight")
    plt.close(fig)


# ─────────────────────────────────────────────────────────────────────────────
# Page 5+ — Per-Experiment Spotlight  (2 experiments per page)
# ─────────────────────────────────────────────────────────────────────────────

def page_experiment_spotlight(pdf, dfs, labels, rankings, merged, top_n=20):
    """
    For each experiment: top-N instances where it outperforms the median most,
    and top-N where it falls behind most.

    Two experiments are shown per page side-by-side.
    """
    n_exp = len(labels)
    cost_cols = [f"cost_{i}" for i in range(n_exp)]

    # Pre-compute per-instance median score and per-experiment relative performance
    merged = merged.copy()
    merged["median_cost"] = merged[cost_cols].median(axis=1)

    # Group into pages of 2
    rank_order = list(rankings["label"])
    pages_of_two = [rank_order[k:k + 2] for k in range(0, n_exp, 2)]

    for page_labels in pages_of_two:
        n_on_page = len(page_labels)
        fig = plt.figure(figsize=(11, 8.5))
        page_title(fig, "Per-Experiment Spotlight: Strengths & Weaknesses",
                   subtitle=(
                       "Δ cost = this experiment − median  "
                       "| negative (green bars) = better than median  "
                       "| positive (red bars) = worse than median"
                   ))

        gs = gridspec.GridSpec(1, n_on_page, figure=fig,
                               top=0.90, bottom=0.07, left=0.06, right=0.97,
                               wspace=0.38)

        for col_pos, lbl in enumerate(page_labels):
            exp_idx = labels.index(lbl)
            sc_col  = f"cost_{exp_idx}"

            valid = merged[["instance", sc_col, "median_cost"]].dropna().copy()
            valid["delta"] = valid[sc_col] - valid["median_cost"]

            # Filter: only instances where there is SOME variance across experiments
            # (skip instances where all experiments have the same score)
            all_same = merged[cost_cols].nunique(axis=1) == 1
            valid = valid[~all_same.reindex(valid.index, fill_value=False)]

            improvements  = valid.nsmallest(top_n, "delta")   # delta << 0 = best
            degradations  = valid.nlargest(top_n, "delta")    # delta >> 0 = worst

            ax = fig.add_subplot(gs[0, col_pos])

            # Combine: improvements on top (most negative first), then degradations
            impr_sorted = improvements.sort_values("delta")
            degr_sorted = degradations.sort_values("delta", ascending=False)
            combined = pd.concat([impr_sorted, degr_sorted], ignore_index=True)

            # Insert a small separator gap between improvements and degradations
            names  = list(combined["instance"])
            deltas = list(combined["delta"])
            n_impr = len(impr_sorted)
            # Insert blank row as separator
            names.insert(n_impr, "─" * 14)
            deltas.insert(n_impr, 0.0)

            y  = np.arange(len(names))
            colors = []
            for k, d in enumerate(deltas):
                if names[k].startswith("─"):
                    colors.append("#cccccc")
                elif d < 0:
                    colors.append(PALETTE["solved"])    # green = improvement
                else:
                    colors.append(PALETTE["overtime"])  # red = degradation

            ax.barh(y, deltas, color=colors, edgecolor="none", height=0.75)
            ax.axvline(0, color=PALETTE["dark"], linewidth=0.8, linestyle="--")

            short_names = [nm[:28] for nm in names]
            ax.set_yticks(y)
            ax.set_yticklabels(short_names, fontsize=5.5)
            ax.invert_yaxis()
            ax.set_xlabel("Δ cost (this − median)")

            n_improve_total = (valid["delta"] < -1).sum()
            n_degrade_total = (valid["delta"] > 1).sum()

            exp_col = exp_color(exp_idx)
            ax.set_title(
                f"{lbl}\n"
                f"Better on {n_improve_total} inst. / worse on {n_degrade_total} inst.\n"
                f"(showing top {top_n} each)",
                color=exp_col, fontsize=9,
            )
            ax.grid(axis="x", alpha=0.3)
            ax.set_axisbelow(True)

            # Add total cost and rank to the panel
            rank_row = rankings[rankings["label"] == lbl]
            if len(rank_row) > 0:
                r = rank_row.iloc[0]
                summary_txt = (
                    f"Rank #{r['rank']}  |  Avg cost: {r['avg_cost']:.2f}  "
                    f"({r['n_inst_own']}/{r['n_inst_total']} inst. run)\n"
                    f"Solved: {r['n_solved']} ({r['pct_solved']:.0f}% of total)  |  "
                    f"Wins: {r['n_wins']}  Sole wins: {r['n_sole_wins']}"
                )
                ax.text(0.98, 0.01, summary_txt, ha="right", va="bottom",
                        transform=ax.transAxes, fontsize=6.5,
                        color=PALETTE["dark"],
                        bbox=dict(boxstyle="round,pad=0.3", facecolor="#f0f0f0",
                                  edgecolor="#ccc", linewidth=0.6))

            # Legend only once
            if col_pos == 0:
                legend_items = [
                    mpatches.Patch(color=PALETTE["solved"],  label="Better than median"),
                    mpatches.Patch(color=PALETTE["overtime"], label="Worse than median"),
                ]
                ax.legend(handles=legend_items, fontsize=7, loc="lower right")

        pdf.savefig(fig, bbox_inches="tight")
        plt.close(fig)


# ─────────────────────────────────────────────────────────────────────────────
# Page N — Instance Detail Table
# ─────────────────────────────────────────────────────────────────────────────

def page_detail_table(pdf, labels, merged, rankings):
    """
    Per-instance table showing costs from all experiments.
    Sorted by: best_cost ascending (easiest/best first), then median_cost.
    Rows are colour-coded by which experiment wins on that instance.
    """
    n_exp = len(labels)
    cost_cols = [f"cost_{i}" for i in range(n_exp)]
    rank_order = list(rankings["label"])

    # Build sorted display frame
    disp = merged[["instance", "expected"] + cost_cols + ["best_cost"]].copy()
    disp = disp.sort_values(["best_cost", "instance"], na_position="last")

    ROWS_PER_PAGE = 45
    pages = [disp.iloc[k:k + ROWS_PER_PAGE]
             for k in range(0, len(disp), ROWS_PER_PAGE)]

    def fmt_cost(v):
        if v != v or np.isnan(v):
            return "—"
        if v >= 300:
            return "OT"
        if v >= 200:
            return "NS"
        return f"{v:.1f}"

    # Column order: instance, expected, then scores in rank order
    rank_indices = [labels.index(lbl) for lbl in rank_order]
    col_headers  = ["Instance", "BKS"] + [f"#{j+1}\n{rank_order[j]}" for j in range(n_exp)]

    for pi, chunk in enumerate(pages):
        fig = plt.figure(figsize=(11, 8.5))
        page_title(fig,
                   f"Instance Detail  ({pi+1}/{len(pages)})",
                   subtitle=f"Columns: instance | BKS | costs in rank order  "
                             f"| green cell = best on that instance")
        ax = fig.add_axes([0.01, 0.03, 0.98, 0.89])
        ax.axis("off")

        cell_text   = []
        cell_colors = []

        for _, row in chunk.iterrows():
            scores_raw = [row[f"cost_{ri}"] for ri in rank_indices]
            best_s = min((s for s in scores_raw if not np.isnan(s)), default=float("nan"))

            bks_val = row["expected"]
            bks_str = f"{bks_val:.4g}" if not np.isnan(bks_val) else "—"

            text_row = [str(row["instance"])[:30], bks_str]
            col_row  = ["#f4f4f4", "#f4f4f4"]
            for s in scores_raw:
                text_row.append(fmt_cost(s))
                if np.isnan(s):
                    col_row.append("#f4f4f4")
                elif abs(s - best_s) < 0.05 and best_s < 200:
                    col_row.append("#d5f5e3")   # green — winner/tied
                elif s >= 300:
                    col_row.append("#c0392b22") # dark red — overtime
                elif s >= 200:
                    col_row.append("#fdecea")   # light red — no solution
                elif s > best_s + 20:
                    col_row.append("#fef9e7")   # light amber — significantly worse
                else:
                    col_row.append("#f9f9f9")

            cell_text.append(text_row)
            cell_colors.append(col_row)

        tbl = ax.table(
            cellText=cell_text,
            colLabels=col_headers,
            cellColours=cell_colors,
            cellLoc="center",
            loc="upper center",
            bbox=[0.0, 0.0, 1.0, 1.0],
        )
        tbl.auto_set_font_size(False)
        tbl.set_fontsize(6.5)
        for (r, c), cell in tbl.get_celld().items():
            cell.set_linewidth(0.3)
            if r == 0:
                cell.set_text_props(fontweight="bold")
                cell.set_facecolor("#d5d8dc")

        pdf.savefig(fig, bbox_inches="tight")
        plt.close(fig)


# ─────────────────────────────────────────────────────────────────────────────
# Main
# ─────────────────────────────────────────────────────────────────────────────

def parse_args():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("dirs", nargs="*", metavar="DIR",
                        help="Experiment directories to compare")
    parser.add_argument("--list", dest="list_file", metavar="FILE",
                        help="File listing experiment dirs (one per line; 'label=dir')")
    parser.add_argument("-o", "--output", metavar="FILE",
                        help="Output PDF path")
    parser.add_argument("-l", "--labels", nargs="+", metavar="LABEL",
                        help="Short labels for each experiment")
    parser.add_argument("--top-n", type=int, default=20, metavar="N",
                        help="Top-N instances per spotlight panel [default: 20]")
    return parser.parse_args()


def collect_dirs_and_labels(args):
    """Return (dirs, labels) from CLI args or --list file."""
    dirs, labels = [], []

    if args.list_file:
        with open(args.list_file) as fh:
            for line in fh:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                if "=" in line:
                    lbl, d = line.split("=", 1)
                    labels.append(lbl.strip())
                    dirs.append(d.strip())
                else:
                    dirs.append(line)
                    labels.append(os.path.basename(line.rstrip("/")))

    for d in args.dirs:
        dirs.append(d)
        labels.append(os.path.basename(d.rstrip("/")))

    if args.labels:
        if len(args.labels) != len(dirs):
            print(f"ERROR: --labels has {len(args.labels)} entries but "
                  f"{len(dirs)} dirs were specified", file=sys.stderr)
            sys.exit(1)
        labels = args.labels

    return dirs, labels


def main():
    args = parse_args()
    dirs, labels = collect_dirs_and_labels(args)

    if len(dirs) < 2:
        print("ERROR: need at least 2 experiment directories", file=sys.stderr)
        sys.exit(1)

    ts = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    if args.output:
        out_pdf = args.output
    else:
        out_pdf = f"compare_multi_{ts}.pdf"

    print(f"Comparing {len(dirs)} experiments:")

    dfs   = []
    metas = []
    for i, (d, lbl) in enumerate(zip(dirs, labels)):
        print(f"  [{i+1}/{len(dirs)}] {lbl}  ({d})")
        df   = load_summary(d)
        meta = load_meta(d)
        print(f"         {len(df)} instances")
        dfs.append(df)
        metas.append(meta)

    n_common = len(set.intersection(*[set(df["instance"]) for df in dfs]))
    print(f"  {n_common} instances common to all experiments\n")

    print("  Building merged table …")
    merged = merge_all(dfs, labels)

    print("  Computing rankings …")
    rankings = compute_rankings(dfs, labels, merged)

    print("\n  Rankings (by avg cost per instance, missing/crashed=200, overtime=300):")
    for _, row in rankings.iterrows():
        print(f"    #{row['rank']:2d}  {row['label']:<30s}  "
              f"avg_cost={row['avg_cost']:.2f}  "
              f"solved={row['n_solved']}/{row['n_inst_total']} ({row['pct_solved']:.1f}%)  "
              f"ran={row['n_inst_own']}  wins={row['n_wins']}")

    print(f"\n  Writing PDF → {out_pdf} …")
    with PdfPages(out_pdf) as pdf:
        d = pdf.infodict()
        d["Title"]   = f"CBC Multi-Experiment Comparison ({len(dirs)} experiments)"
        d["Creator"] = "compare_multi_experiments.py"

        print("  Page 0: Introduction")
        page_intro(pdf, dirs, labels, dfs, metas)

        print("  Page 1: Executive summary")
        page_executive_summary(pdf, dfs, labels, metas, rankings, merged)

        print("  Page 2: Cost profiles & head-to-head")
        page_cost_profiles(pdf, dfs, labels, merged, rankings)

        print("  Page 3: BKS gap analysis")
        page_bks_gaps(pdf, dfs, labels)

        print("  Page 4: Category heatmap")
        page_category_heatmap(pdf, dfs, labels, merged)

        n_spotlight_pages = (len(labels) + 1) // 2
        print(f"  Pages 5–{4 + n_spotlight_pages}: Per-experiment spotlights")
        page_experiment_spotlight(pdf, dfs, labels, rankings, merged,
                                  top_n=args.top_n)

        print(f"  Pages {5 + n_spotlight_pages}+: Instance detail tables")
        page_detail_table(pdf, labels, merged, rankings)

    print(f"\nDone.  → {out_pdf}")


if __name__ == "__main__":
    main()
