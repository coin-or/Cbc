#!/usr/bin/env python3
"""
make_lp_report.py — Generate a PDF report from LP relaxation parameter tuning results.

Reads the outputs of analyze_lp_params.py and analyze_racing_portfolios.py and
produces a multi-page PDF with charts and tables summarising the discoveries.

Usage:
    python3 make_lp_report.py <experiment_dir> [output.pdf] [OPTIONS]

    # Run analyses first if CSVs are missing:
    python3 analyze_lp_params.py --dir <exp_dir> --timelimit 7200
    python3 analyze_racing_portfolios.py --dir <exp_dir> --timelimit 7200

Options:
    --timelimit T       LP time limit used in experiment (default: 7200)
    --penalty-mult M    Penalty multiplier (default: 2.0)
    --shift S           SGM shift in seconds (default: 1.0)
    --baseline PARAM    Baseline param for comparison (default: dual_default)
    --top N             Number of top params to highlight (default: 8)
    --exhaustive-k K    Max K for exhaustive portfolio search (default: 4)
    --max-k K           Max portfolio size to show (default: 8)
    --dtree-depth N     Decision tree max depth (default: 3)
    --no-dtree          Skip decision tree page
    --features PATH     Instance features CSV (default: ~/inst/miplib/2017+spp/features.csv)
    --exclude-prefix S  Exclude params starting with S from dtree (default: guess_)
"""

import argparse
import csv
import itertools
import math
import os
import sys
import tempfile
from collections import defaultdict

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.ticker import LogLocator, NullFormatter, FuncFormatter
from matplotlib.lines import Line2D

# fbps decision tree library (sibling directory)
_SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(_SCRIPT_DIR, "fbps", "fbps"))
try:
    from instance import InstanceSet
    from psetting import PSetting
    from dtree import DTree
    _FBPS_AVAILABLE = True
except ImportError:
    _FBPS_AVAILABLE = False

# ── Style ─────────────────────────────────────────────────────────────────────

PALETTE = {
    "blue":   "#2C7BB6",
    "red":    "#D7191C",
    "green":  "#1A9641",
    "orange": "#FDAE61",
    "purple": "#762A83",
    "grey":   "#888888",
    "teal":   "#00929F",
    "yellow": "#F4D03F",
}
CMAP_QUAL = [PALETTE["blue"], PALETTE["red"], PALETTE["green"], PALETTE["orange"],
             PALETTE["purple"], PALETTE["teal"], PALETTE["yellow"], PALETTE["grey"]]

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 9,
    "axes.titlesize": 10,
    "axes.labelsize": 9,
    "xtick.labelsize": 8,
    "ytick.labelsize": 8,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "figure.dpi": 150,
})


# ── Helpers ───────────────────────────────────────────────────────────────────

def shifted_geomean(values, shift):
    if not values:
        return float("nan")
    return math.exp(sum(math.log(v + shift) for v in values) / len(values)) - shift


def short_name(param, maxlen=22):
    return param if len(param) <= maxlen else param[:maxlen - 1] + "…"


def page_title(fig, text, subtitle=""):
    fig.text(0.5, 0.97, text, ha="center", va="top", fontsize=13, fontweight="bold")
    if subtitle:
        fig.text(0.5, 0.94, subtitle, ha="center", va="top", fontsize=9,
                 color=PALETTE["grey"])


def add_page_footer(fig, exp_name, page_num):
    fig.text(0.5, 0.01, f"{exp_name}  •  page {page_num}",
             ha="center", fontsize=7, color=PALETTE["grey"])


def hbar_chart(ax, labels, values, color=PALETTE["blue"], baseline_val=None,
               baseline_label=None, xlabel="", title="", log=False):
    """Horizontal bar chart, sorted descending."""
    y = range(len(labels))
    bars = ax.barh(list(y), values, color=color, alpha=0.82, height=0.65)
    ax.set_yticks(list(y))
    ax.set_yticklabels(labels, fontsize=8)
    ax.invert_yaxis()
    ax.set_xlabel(xlabel)
    ax.set_title(title, pad=4)
    if log:
        ax.set_xscale("log")
    if baseline_val is not None:
        ax.axvline(baseline_val, color=PALETTE["red"], lw=1.5, ls="--",
                   label=baseline_label or "baseline")
        ax.legend(fontsize=7, loc="lower right")
    # Value labels on bars
    for bar, val in zip(bars, values):
        if math.isnan(val):
            continue
        ax.text(bar.get_width() * (1.02 if not log else 1.05),
                bar.get_y() + bar.get_height() / 2,
                f"{val:.2f}" if val < 1000 else f"{val:.0f}",
                va="center", ha="left", fontsize=7)
    ax.margins(x=0.15)


# ── Data loading ──────────────────────────────────────────────────────────────

def load_avg_times(exp_dir, penalty, exclude_prefix="guess_"):
    """Load lp_avg_times.csv → DataFrames and dicts."""
    csv_path = os.path.join(exp_dir, "lp_avg_times.csv")
    if not os.path.isfile(csv_path):
        return None, None, None, None
    df = pd.read_csv(csv_path)
    if exclude_prefix:
        df = df[~df["param_tag"].str.startswith(exclude_prefix)]
    df["avg_wall_seconds"] = pd.to_numeric(df["avg_wall_seconds"], errors="coerce").fillna(penalty)
    df["all_seeds_solved"] = df["all_seeds_solved"].astype(str).str.strip() == "True"
    df["any_seed_solved"]  = df["any_seed_solved"].astype(str).str.strip() == "True"
    instances = sorted(df["instance"].unique())
    params    = sorted(df["param_tag"].unique())
    # Pivot: rows=instances, cols=params
    time_piv   = df.pivot_table("avg_wall_seconds", "instance", "param_tag").reindex(
        index=instances, columns=params).fillna(penalty)
    solved_piv = df.pivot_table("all_seeds_solved", "instance", "param_tag",
                                aggfunc="max").reindex(
        index=instances, columns=params).fillna(False).astype(bool)
    return df, time_piv, solved_piv, instances, params


def load_racing_curve(exp_dir):
    """Load lp_racing_curve.csv if available."""
    csv_path = os.path.join(exp_dir, "lp_racing_curve.csv")
    if not os.path.isfile(csv_path):
        return None
    return pd.read_csv(csv_path)


def load_raw_results(exp_dir):
    """Load lp_results.csv for status breakdown."""
    csv_path = os.path.join(exp_dir, "lp_results.csv")
    if not os.path.isfile(csv_path):
        return None
    return pd.read_csv(csv_path)


# ── Portfolio helpers ─────────────────────────────────────────────────────────

def portfolio_sgm(time_np, solved_np, port_cols, shift):
    sub_t = time_np[:, port_cols]
    racing = sub_t.min(axis=1)
    return math.exp(float(np.log(racing + shift).mean())) - shift


def portfolio_solve_rate(solved_np, port_cols):
    return float(solved_np[:, port_cols].any(axis=1).mean()) * 100


def greedy_portfolio(time_np, solved_np, params, shift, max_k):
    remaining = list(range(len(params)))
    portfolio = []
    curve = []  # (k, sgm, solve_rate, param_added)
    for k in range(1, max_k + 1):
        best_sgm, best_p = float("inf"), None
        for p in remaining:
            sgm = portfolio_sgm(time_np, solved_np, portfolio + [p], shift)
            if sgm < best_sgm:
                best_sgm, best_p = sgm, p
        portfolio.append(best_p)
        remaining.remove(best_p)
        sr = portfolio_solve_rate(solved_np, portfolio)
        curve.append((k, best_sgm, sr, params[best_p]))
    return curve, portfolio


def exhaustive_optimal(time_np, solved_np, params, shift, max_k):
    results = {}
    for k in range(1, max_k + 1):
        best_sgm, best_combo = float("inf"), None
        for combo in itertools.combinations(range(len(params)), k):
            sgm = portfolio_sgm(time_np, solved_np, list(combo), shift)
            if sgm < best_sgm:
                best_sgm, best_combo = sgm, combo
        results[k] = (best_sgm, [params[i] for i in best_combo])
    return results


# ── Pages ─────────────────────────────────────────────────────────────────────

def page_overview(pdf, exp_name, raw_df, time_piv, solved_piv, penalty, timelimit, args):
    """Page 1 — Overview: instance count, solve rates, status breakdown."""
    fig = plt.figure(figsize=(11, 8.5))
    page_title(fig, f"LP Relaxation Experiment — {exp_name}",
               f"Instances: {time_piv.shape[0]}  |  Params: {time_piv.shape[1]}"
               f"  |  Timelimit: {timelimit}s  |  Penalty: {penalty:.0f}s")

    gs = fig.add_gridspec(2, 3, left=0.07, right=0.97, top=0.90, bottom=0.08,
                          hspace=0.45, wspace=0.38)

    params    = list(time_piv.columns)
    instances = list(time_piv.index)
    n_inst    = len(instances)
    shift     = args.shift

    # ── (0,0)+(0,1): Solve rate per param (horizontal bar) ───────────────────
    ax_sr = fig.add_subplot(gs[0, :2])
    solve_rates = [(p, solved_piv[p].mean() * 100) for p in params]
    solve_rates.sort(key=lambda x: -x[1])
    snames = [short_name(p) for p, _ in solve_rates]
    svals  = [v for _, v in solve_rates]
    baseline_sr = next((v for p, v in solve_rates if p == args.baseline), None)
    hbar_chart(ax_sr, snames, svals, color=PALETTE["blue"],
               baseline_val=baseline_sr,
               baseline_label=f"{args.baseline} ({baseline_sr:.1f}%)" if baseline_sr else None,
               xlabel="Solve rate (%, all seeds solved)", title="Solve rate per parameter")

    # ── (0,2): Status breakdown pie ──────────────────────────────────────────
    ax_pie = fig.add_subplot(gs[0, 2])
    if raw_df is not None:
        counts = raw_df["status"].value_counts()
        labels = list(counts.index)
        sizes  = list(counts.values)
        colors_pie = [PALETTE["green"] if l == "OPTIMAL"
                      else PALETTE["orange"] if l in ("TIMEOUT", "TIMEOUT_KILLED")
                      else PALETTE["red"] if l == "ERROR"
                      else PALETTE["grey"] for l in labels]
        wedges, texts, autotexts = ax_pie.pie(
            sizes, labels=None, autopct="%1.1f%%", colors=colors_pie,
            startangle=90, pctdistance=0.75)
        for at in autotexts:
            at.set_fontsize(7)
        ax_pie.legend(wedges, [f"{l} ({v})" for l, v in zip(labels, sizes)],
                      loc="lower center", bbox_to_anchor=(0.5, -0.28),
                      fontsize=7, ncol=1)
        ax_pie.set_title("Run status breakdown", pad=4)
    else:
        ax_pie.text(0.5, 0.5, "lp_results.csv\nnot found", ha="center", va="center")
        ax_pie.axis("off")

    # ── (1,0)+(1,1): SGM(all) bar chart ──────────────────────────────────────
    ax_sgm = fig.add_subplot(gs[1, :2])
    sgm_vals = []
    for p in params:
        times = list(time_piv[p])
        sgm_vals.append((p, shifted_geomean(times, shift)))
    sgm_vals.sort(key=lambda x: x[1])
    snames2 = [short_name(p) for p, _ in sgm_vals]
    svals2  = [v for _, v in sgm_vals]
    baseline_sgm = next((v for p, v in sgm_vals if p == args.baseline), None)
    hbar_chart(ax_sgm, snames2, svals2, color=PALETTE["teal"],
               baseline_val=baseline_sgm,
               baseline_label=f"{args.baseline}" if baseline_sgm else None,
               xlabel="SGM solve time (s, penalty for failures)",
               title=f"SGM over ALL instances (shift={shift}s, penalty={penalty:.0f}s)")

    # ── (1,2): Instance difficulty histogram ─────────────────────────────────
    ax_hist = fig.add_subplot(gs[1, 2])
    # Best time per instance across all params
    best_times = time_piv.min(axis=1)
    solved_mask = solved_piv.any(axis=1)
    solved_times = best_times[solved_mask]
    bins = np.logspace(np.log10(0.001), np.log10(timelimit), 30)
    ax_hist.hist(solved_times, bins=bins, color=PALETTE["blue"], alpha=0.8, edgecolor="white")
    ax_hist.set_xscale("log")
    ax_hist.set_xlabel("Best solve time (s) — solved instances")
    ax_hist.set_ylabel("# Instances")
    ax_hist.set_title(f"Instance difficulty distribution\n"
                      f"({solved_mask.sum()} solved / {n_inst} total)")
    ax_hist.xaxis.set_major_formatter(FuncFormatter(lambda x, _: f"{x:.3g}"))

    add_page_footer(fig, exp_name, 1)
    pdf.savefig(fig, bbox_inches="tight")
    plt.close(fig)


def page_param_rankings(pdf, exp_name, time_piv, solved_piv, penalty, args):
    """Page 2 — Param rankings: SGM(all), solve rates, speedup distribution, freq winner."""
    fig = plt.figure(figsize=(11, 8.5))
    page_title(fig, "Parameter Rankings", exp_name)

    params    = list(time_piv.columns)
    instances = list(time_piv.index)
    shift     = args.shift
    n_inst    = len(instances)

    gs = fig.add_gridspec(2, 2, left=0.18, right=0.97, top=0.90, bottom=0.08,
                          hspace=0.45, wspace=0.42)

    # ── SGM(all) ranking ─────────────────────────────────────────────────────
    ax1 = fig.add_subplot(gs[0, 0])
    sgm_all = {p: shifted_geomean(list(time_piv[p]), shift) for p in params}
    sorted_params = sorted(params, key=lambda p: sgm_all[p])
    names1 = [short_name(p) for p in sorted_params]
    vals1  = [sgm_all[p] for p in sorted_params]
    baseline_val = sgm_all.get(args.baseline)
    colors1 = [PALETTE["orange"] if p == args.baseline else PALETTE["teal"]
                for p in sorted_params]
    ax1.barh(range(len(names1)), vals1, color=colors1, alpha=0.85, height=0.65)
    ax1.set_yticks(range(len(names1)))
    ax1.set_yticklabels(names1, fontsize=7)
    ax1.invert_yaxis()
    ax1.set_xlabel("SGM (s)")
    ax1.set_title("SGM — all instances\n(lower = better; failures → penalty)", pad=4)
    if baseline_val:
        ax1.axvline(baseline_val, color=PALETTE["red"], lw=1.2, ls="--")
    patch = mpatches.Patch(color=PALETTE["orange"], label="baseline")
    ax1.legend(handles=[patch], fontsize=7, loc="lower right")

    # ── Solve rate ranking ────────────────────────────────────────────────────
    ax2 = fig.add_subplot(gs[0, 1])
    solve_rates = {p: solved_piv[p].mean() * 100 for p in params}
    sorted_sr = sorted(params, key=lambda p: solve_rates[p], reverse=True)
    names2 = [short_name(p) for p in sorted_sr]
    vals2  = [solve_rates[p] for p in sorted_sr]
    colors2 = [PALETTE["orange"] if p == args.baseline else PALETTE["blue"]
                for p in sorted_sr]
    ax2.barh(range(len(names2)), vals2, color=colors2, alpha=0.85, height=0.65)
    ax2.set_yticks(range(len(names2)))
    ax2.set_yticklabels(names2, fontsize=7)
    ax2.invert_yaxis()
    ax2.set_xlabel("Solve rate (%)")
    ax2.set_title("Solve rate (all seeds solved)\n(higher = better; failures counted against)", pad=4)
    ax2.axvline(solve_rates.get(args.baseline, 0), color=PALETTE["red"], lw=1.2, ls="--")
    patch2 = mpatches.Patch(color=PALETTE["orange"], label="baseline")
    ax2.legend(handles=[patch2], fontsize=7, loc="lower right")

    # ── Speedup vs baseline scatter (each param vs baseline, per instance) ────
    ax3 = fig.add_subplot(gs[1, 0])
    if args.baseline in params:
        base_times = time_piv[args.baseline].values
        top_params = sorted_params[:args.top]  # top params by SGM(all)
        for idx, p in enumerate(top_params[:5]):
            if p == args.baseline:
                continue
            p_times = time_piv[p].values
            # Only include instances where both solved
            mask = solved_piv[args.baseline].values & solved_piv[p].values
            if mask.sum() < 3:
                continue
            speedups = base_times[mask] / np.maximum(p_times[mask], 0.001)
            label = short_name(p, 18)
            ax3.hist(np.log10(speedups), bins=30, alpha=0.5,
                     label=f"{label}", color=CMAP_QUAL[idx % len(CMAP_QUAL)])
        ax3.axvline(0, color="black", lw=1, ls="--")
        ax3.set_xlabel("log₁₀(speedup vs baseline)")
        ax3.set_ylabel("# Instances")
        ax3.set_title(f"Speedup distribution vs {short_name(args.baseline)}\n"
                      "(>0 = faster than baseline)", pad=4)
        ax3.legend(fontsize=6, loc="upper left")
    else:
        ax3.text(0.5, 0.5, f"Baseline '{args.baseline}'\nnot found", ha="center", va="center")
        ax3.axis("off")

    # ── % Fastest per instance bar chart (top params) ─────────────────────────
    ax4 = fig.add_subplot(gs[1, 1])
    best_counts = defaultdict(int)
    for inst in instances:
        best_p = time_piv.loc[inst].idxmin()
        best_counts[best_p] += 1
    top_best = sorted(best_counts.items(), key=lambda x: -x[1])[:15]
    names4 = [short_name(p) for p, _ in top_best]
    vals4  = [c / n_inst * 100 for _, c in top_best]
    colors4 = [PALETTE["orange"] if p == args.baseline else PALETTE["green"]
                for p, _ in top_best]
    ax4.barh(range(len(names4)), vals4, color=colors4, alpha=0.85, height=0.65)
    ax4.set_yticks(range(len(names4)))
    ax4.set_yticklabels(names4, fontsize=7)
    ax4.invert_yaxis()
    ax4.set_xlabel("% instances where this param is fastest")
    ax4.set_title("Most frequent winner\n(% of instances where param had lowest avg time)", pad=4)
    for i, v in enumerate(vals4):
        ax4.text(v + 0.3, i, f"{v:.1f}%", va="center", fontsize=7)

    add_page_footer(fig, exp_name, 2)
    pdf.savefig(fig, bbox_inches="tight")
    plt.close(fig)


def page_performance_profiles(pdf, exp_name, time_piv, solved_piv, penalty, args):
    """Page 3 — Performance profiles (time ratio CDF) for top params."""
    fig = plt.figure(figsize=(11, 8.5))
    page_title(fig, "Performance Profiles", exp_name)

    params    = list(time_piv.columns)
    instances = list(time_piv.index)
    shift     = args.shift

    gs = fig.add_gridspec(2, 2, left=0.10, right=0.97, top=0.90, bottom=0.08,
                          hspace=0.45, wspace=0.38)

    # Rank params by SGM(all)
    sgm_all = {p: shifted_geomean(list(time_piv[p]), shift) for p in params}
    ranked  = sorted(params, key=lambda p: sgm_all[p])
    top     = ranked[:args.top]

    # ── (0,0)+(0,1): Performance profile (Dolan-Moré style, normalized to ALL instances) ──
    ax_pp = fig.add_subplot(gs[0, :])
    n_total = len(instances)
    for idx, p in enumerate(top):
        times_p = []
        for inst in instances:
            if not solved_piv.loc[inst, p]:
                continue
            t_p = time_piv.loc[inst, p]
            best_t = time_piv.loc[inst, ranked].min()
            ratio  = t_p / max(best_t, 0.001)
            times_p.append(ratio)
        times_p.sort()
        # Normalize y to ALL instances so curves plateau at the solve rate
        ys = [(i + 1) / n_total for i in range(len(times_p))]
        label = short_name(p, 20)
        ax_pp.step(times_p, ys, where="post", label=label,
                   color=CMAP_QUAL[idx % len(CMAP_QUAL)], lw=1.5, alpha=0.85)

    ax_pp.set_xscale("log")
    ax_pp.set_xlabel("Time ratio τ  (param_time / best_time_for_instance)")
    ax_pp.set_ylabel("Fraction of ALL instances with ratio ≤ τ")
    ax_pp.set_title(f"Performance profile — top {args.top} params\n"
                    "(higher = better; curves plateau at solve rate — unsolved instances never appear)", pad=4)
    ax_pp.set_xlim(left=1.0)
    ax_pp.set_ylim(0, 1.05)
    ax_pp.axvline(1.0, color="black", lw=0.8, ls="--", alpha=0.5)
    ax_pp.legend(fontsize=7, loc="lower right", ncol=2)
    ax_pp.grid(True, which="both", alpha=0.2)

    # ── (1,0): Head-to-head scatter: #1 vs #2 param ──────────────────────────
    ax_scat = fig.add_subplot(gs[1, 0])
    if len(ranked) >= 2:
        p1, p2 = ranked[0], ranked[1]
        xs = time_piv[p1].values
        ys2 = time_piv[p2].values
        both = solved_piv[p1].values & solved_piv[p2].values
        c_both   = PALETTE["blue"]
        c_p1only = PALETTE["green"]
        c_p2only = PALETTE["red"]
        # Scatter: solved by both
        ax_scat.scatter(xs[both], ys2[both], s=10, alpha=0.5,
                        color=c_both, label="both solved")
        # p1 solved, p2 timeout → show at penalty y
        p1_only = solved_piv[p1].values & ~solved_piv[p2].values
        if p1_only.any():
            ax_scat.scatter(xs[p1_only], ys2[p1_only], s=10, alpha=0.4,
                            color=c_p1only, marker="^",
                            label=f"only {short_name(p1,12)} solved")
        p2_only = ~solved_piv[p1].values & solved_piv[p2].values
        if p2_only.any():
            ax_scat.scatter(xs[p2_only], ys2[p2_only], s=10, alpha=0.4,
                            color=c_p2only, marker="v",
                            label=f"only {short_name(p2,12)} solved")
        lo = min(xs[both].min(), ys2[both].min()) * 0.8 if both.any() else 0.001
        hi = max(xs[both].max(), ys2[both].max()) * 1.2 if both.any() else 1000
        ax_scat.plot([lo, hi], [lo, hi], "k--", lw=0.8, alpha=0.5)
        ax_scat.set_xscale("log"); ax_scat.set_yscale("log")
        ax_scat.set_xlabel(f"{short_name(p1)} (s)")
        ax_scat.set_ylabel(f"{short_name(p2)} (s)")
        ax_scat.set_title(f"Head-to-head: rank 1 vs rank 2\n"
                          f"(below diagonal = {short_name(p2,14)} faster)", pad=4)
        ax_scat.legend(fontsize=7)
        ax_scat.grid(True, which="both", alpha=0.15)

    # ── (1,1): CDF of solve times for top params ─────────────────────────────
    ax_cdf = fig.add_subplot(gs[1, 1])
    max_t = max(time_piv[p][solved_piv[p]].max()
                for p in top if solved_piv[p].any()) * 1.05
    for idx, p in enumerate(top):
        t_vals = sorted(time_piv[p][solved_piv[p]].values)
        if not t_vals:
            continue
        ys3 = [(i + 1) / len(instances) for i in range(len(t_vals))]
        ax_cdf.step(t_vals, ys3, where="post",
                    label=short_name(p, 18),
                    color=CMAP_QUAL[idx % len(CMAP_QUAL)], lw=1.4, alpha=0.85)
    ax_cdf.set_xscale("log")
    ax_cdf.set_xlabel("Solve time (s)")
    ax_cdf.set_ylabel("Fraction of ALL instances solved")
    ax_cdf.set_title("Solve time CDF (fraction of all instances solved by time t)", pad=4)
    ax_cdf.set_ylim(0, 1.05)
    ax_cdf.legend(fontsize=6, loc="upper left", ncol=2)
    ax_cdf.grid(True, which="both", alpha=0.2)

    add_page_footer(fig, exp_name, 3)
    pdf.savefig(fig, bbox_inches="tight")
    plt.close(fig)


def page_racing_portfolios(pdf, exp_name, time_piv, solved_piv, args,
                           greedy_curve, exhaustive_results):
    """Page 4 — Racing portfolio analysis: SGM curve, coverage, contributions."""
    fig = plt.figure(figsize=(11, 8.5))
    page_title(fig, "Racing LP Portfolio Analysis",
               f"racingLP: K solvers run in parallel; effective time = min(t₁,…,tₖ)")

    gs = fig.add_gridspec(2, 2, left=0.10, right=0.97, top=0.88, bottom=0.08,
                          hspace=0.50, wspace=0.38)

    shift     = args.shift
    params    = list(time_piv.columns)
    instances = list(time_piv.index)
    n_inst    = len(instances)

    # ── (0,0): SGM vs K curve ─────────────────────────────────────────────────
    ax_curve = fig.add_subplot(gs[0, 0])
    ks_g    = [row[0] for row in greedy_curve]
    sgms_g  = [row[1] for row in greedy_curve]
    srs_g   = [row[2] for row in greedy_curve]

    ax_curve.plot(ks_g, sgms_g, "o-", color=PALETTE["blue"], lw=2,
                  ms=6, label="Greedy portfolio")

    # Exhaustive overlay
    if exhaustive_results:
        ks_e  = sorted(exhaustive_results.keys())
        sgms_e = [exhaustive_results[k][0] for k in ks_e]
        ax_curve.plot(ks_e, sgms_e, "s--", color=PALETTE["red"], lw=1.5,
                      ms=5, label="Exhaustive optimal", zorder=5)

    # Baseline K=1 reference
    base_sgm = shifted_geomean(list(time_piv[args.baseline]), shift) if args.baseline in params else None
    if base_sgm:
        ax_curve.axhline(base_sgm, color=PALETTE["grey"], lw=1, ls=":",
                         label=f"baseline ({short_name(args.baseline)})")

    # Annotate each point with param added
    for k, sgm, sr, p_added in greedy_curve:
        ax_curve.annotate(short_name(p_added, 14),
                          xy=(k, sgm), xytext=(4, 4),
                          textcoords="offset points", fontsize=6, color=PALETTE["blue"])

    ax_curve.set_yscale("log")
    ax_curve.set_xlabel("Portfolio size K (# parallel runners)")
    ax_curve.set_ylabel("SGM solve time (s, log scale)")
    ax_curve.set_title("Racing portfolio: SGM vs portfolio size\n"
                       "(lower = better; annotated with param added at each step)", pad=4)
    ax_curve.set_xticks(ks_g)
    ax_curve.legend(fontsize=7)
    ax_curve.grid(True, which="both", alpha=0.2)

    # ── (0,1): Solve rate vs K ────────────────────────────────────────────────
    ax_sr = fig.add_subplot(gs[0, 1])
    ax_sr.plot(ks_g, srs_g, "o-", color=PALETTE["green"], lw=2, ms=6)
    if base_sgm and args.baseline in params:
        base_sr = solved_piv[args.baseline].mean() * 100
        ax_sr.axhline(base_sr, color=PALETTE["grey"], lw=1, ls=":",
                      label=f"baseline ({base_sr:.1f}%)")
    for k, _, sr, p_added in greedy_curve:
        ax_sr.annotate(f"{sr:.1f}%", xy=(k, sr), xytext=(2, 4),
                       textcoords="offset points", fontsize=7)
    ax_sr.set_xlabel("Portfolio size K")
    ax_sr.set_ylabel("Solve rate (%)")
    ax_sr.set_title("Racing portfolio: solve rate vs K\n"
                    "(fraction of instances solved by ≥1 param)", pad=4)
    ax_sr.set_xticks(ks_g)
    ax_sr.set_ylim(0, 105)
    ax_sr.legend(fontsize=7)
    ax_sr.grid(True, alpha=0.2)

    # ── (1,0): Marginal SGM improvement at each step ──────────────────────────
    ax_marg = fig.add_subplot(gs[1, 0])
    gains = []
    prev = None
    for k, sgm, sr, p_added in greedy_curve:
        if prev is not None:
            gains.append((k, (prev - sgm) / prev * 100, p_added))
        prev = sgm
    if gains:
        ks_m    = [g[0] for g in gains]
        g_vals  = [g[1] for g in gains]
        p_names = [short_name(g[2], 16) for g in gains]
        bars = ax_marg.bar(ks_m, g_vals, color=PALETTE["orange"], alpha=0.85)
        ax_marg.set_xticks(ks_m)
        ax_marg.set_xticklabels([f"K={k}\n{n}" for k, n in zip(ks_m, p_names)],
                                 fontsize=7)
        ax_marg.set_xlabel("Portfolio size K (param added)")
        ax_marg.set_ylabel("ΔSGM% (marginal improvement)")
        ax_marg.set_title("Marginal SGM gain at each step\n"
                           "(where does adding a runner stop helping?)", pad=4)
        for bar, val in zip(bars, g_vals):
            ax_marg.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.3,
                         f"{val:.1f}%", ha="center", fontsize=7)
        ax_marg.grid(True, axis="y", alpha=0.2)

    # ── (1,1): Per-param contribution table at key K values ──────────────────
    ax_tbl = fig.add_subplot(gs[1, 1])
    ax_tbl.axis("off")

    # Build contribution table: how many instances each param wins at K=1..5
    show_ks = [k for k in [1, 2, 3, 4, 5] if k <= len(greedy_curve)]
    port_by_k = {}
    running = []
    param_to_idx = {p: i for i, p in enumerate(params)}
    time_np   = time_piv.values
    solved_np = solved_piv.values

    for k, sgm, sr, p_added in greedy_curve:
        running.append(p_added)
        port_by_k[k] = list(running)

    table_data = []
    col_labels = ["Param"] + [f"K={k}" for k in show_ks]
    for p in (port_by_k.get(max(show_ks), []) if show_ks else []):
        row = [short_name(p, 20)]
        for k in show_ks:
            port = port_by_k.get(k, [])
            port_idx = [param_to_idx[pp] for pp in port if pp in param_to_idx]
            if not port_idx:
                row.append(".")
                continue
            sub_t = time_np[:, port_idx]
            winner_col = sub_t.argmin(axis=1)
            p_col = param_to_idx.get(p)
            n_wins = sum(1 for i, wc in enumerate(winner_col)
                         if port_idx[wc] == p_col and solved_np[i, p_col])
            row.append(str(n_wins) if n_wins > 0 else ".")
        table_data.append(row)

    if table_data:
        tbl = ax_tbl.table(
            cellText=table_data,
            colLabels=col_labels,
            cellLoc="center", loc="center",
            bbox=[0, 0, 1, 1])
        tbl.auto_set_font_size(False)
        tbl.set_fontsize(8)
        # Header styling
        for j in range(len(col_labels)):
            tbl[0, j].set_facecolor(PALETTE["blue"])
            tbl[0, j].set_text_props(color="white", fontweight="bold")
        ax_tbl.set_title("Instances 'won' per param in greedy portfolio\n"
                         "(# instances where this param was fastest)", pad=4)

    add_page_footer(fig, exp_name, 4)
    pdf.savefig(fig, bbox_inches="tight")
    plt.close(fig)


# ── Decision tree helpers ──────────────────────────────────────────────────────

def _build_dtree_results(iset, avg_times_csv, exclude_prefix, penalty):
    """Build a minimal results object from lp_avg_times.csv for the fbps DTree."""
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

    for inst in iset.instances:
        inst.results = [penalty] * len(psettings)

    with open(avg_times_csv, newline="") as f:
        for row in csv.DictReader(f):
            p = row["param_tag"]
            if exclude_prefix and p.startswith(exclude_prefix):
                continue
            name = row["instance"]
            if name not in iset.instByName or p not in psetting_map:
                continue
            iset.instByName[name].results[psetting_map[p].idx] = float(row["avg_wall_seconds"])

    class _R:
        pass
    r = _R()
    r.psettings = psettings
    r.psettingByName = psetting_map
    return r


def _render_dtree_png(tree, features, png_path):
    """Render decision tree to PNG via graphviz."""
    from graphviz import Digraph

    def _add(dot, node, parent_id=None, edge_label=""):
        nid = node.node_id
        if not node.children_nodes:
            ps = node.bestPS[0]
            n = len(node.instances)
            label = f"{ps.setting}\n({n} inst)"
            dot.node(nid, label, shape="box", style="filled,rounded",
                     fillcolor="#FFFACD", fontsize="10")
        else:
            feat = features[node.branch_feat_idx]
            bval = node.branch_value
            bval_s = f"{bval:.4g}" if isinstance(bval, float) else str(bval)
            label = f"{feat}\n≤ {bval_s}?"
            dot.node(nid, label, shape="diamond", style="filled",
                     fillcolor="#B0D4F1", fontsize="10")
        if parent_id:
            dot.edge(parent_id, nid, label=edge_label, fontsize="9")
        for i, child in enumerate(node.children_nodes):
            bval = node.branch_value
            bval_s = f"{bval:.4g}" if isinstance(bval, float) else str(bval)
            el = f" ≤{bval_s}" if i == 0 else f" >{bval_s}"
            _add(dot, child, nid, el)

    dot = Digraph(format="png")
    dot.attr("graph", rankdir="TB", dpi="150", bgcolor="white")
    dot.attr("node", fontname="DejaVu Sans")
    dot.attr("edge", fontname="DejaVu Sans")
    _add(dot, tree.root)

    base = png_path.replace(".png", "")
    dot.render(base, cleanup=True)
    rendered = base + ".png"
    if rendered != png_path and os.path.exists(rendered):
        os.rename(rendered, png_path)


def page_decision_tree(pdf, exp_name, exp_dir, penalty, args, page_num):
    """Decision tree page — build fbps tree, render PNG, embed + leaf table."""
    fig = plt.figure(figsize=(11, 8.5))
    page_title(fig, "Feature-Based LP Parameter Decision Tree", exp_name)

    if not _FBPS_AVAILABLE:
        fig.text(0.5, 0.5, "fbps library not found — skipping decision tree page.",
                 ha="center", va="center", fontsize=12, color=PALETTE["red"])
        add_page_footer(fig, exp_name, page_num)
        pdf.savefig(fig, bbox_inches="tight")
        plt.close(fig)
        return

    features_csv = args.features
    avg_times_csv = os.path.join(exp_dir, "lp_avg_times.csv")
    exclude_prefix = args.exclude_prefix

    if not os.path.isfile(features_csv):
        fig.text(0.5, 0.5, f"Features file not found:\n{features_csv}",
                 ha="center", va="center", fontsize=11, color=PALETTE["red"])
        add_page_footer(fig, exp_name, page_num)
        pdf.savefig(fig, bbox_inches="tight")
        plt.close(fig)
        return

    # Build tree (suppress stdout chatter from fbps)
    import io as _io
    import contextlib
    buf = _io.StringIO()
    with contextlib.redirect_stdout(buf):
        iset = InstanceSet(features_csv)
        results = _build_dtree_results(iset, avg_times_csv, exclude_prefix, penalty)
        baseline_str = args.baseline if args.baseline in results.psettingByName else ""
        tree = DTree(iset, results, args.dtree_depth, args.dtree_min_leaf, baseline_str)
        tree.build()

    n_leaves = len(tree.leafs)
    baseline_total = tree.default_time if baseline_str else None
    tree_total = sum(
        inst.results[leaf.bestPS[0].idx]
        for leaf in tree.leafs for inst in leaf.instances)
    speedup = (baseline_total / tree_total) if baseline_total else None

    # Render PNG to temp file
    with tempfile.NamedTemporaryFile(suffix=".png", delete=False) as tf:
        png_path = tf.name
    try:
        _render_dtree_png(tree, iset.features, png_path)
        img_ok = os.path.exists(png_path)
    except Exception as e:
        print(f"  Warning: graphviz render failed: {e}", file=sys.stderr)
        img_ok = False

    # ── Layout: tree image (top 65%) + leaf summary table (bottom 35%) ────────
    gs = fig.add_gridspec(2, 1, top=0.91, bottom=0.04, hspace=0.08,
                          height_ratios=[3, 1.4])

    # Tree image
    ax_img = fig.add_subplot(gs[0])
    ax_img.axis("off")
    if img_ok:
        import matplotlib.image as mpimg
        img = mpimg.imread(png_path)
        ax_img.imshow(img, aspect="equal", interpolation="lanczos")
        os.unlink(png_path)
    else:
        ax_img.text(0.5, 0.5, "Tree rendering unavailable (graphviz error)",
                    ha="center", va="center", fontsize=11)

    # Subtitle with speedup info
    subtitle_parts = [f"depth={args.dtree_depth}", f"min_leaf={args.dtree_min_leaf}",
                      f"{n_leaves} leaves", f"{len(iset.instances)} instances",
                      f"{len(results.psettings)} params"]
    if speedup:
        subtitle_parts.append(f"speedup vs {args.baseline}: {speedup:.2f}×")
    ax_img.set_title("  |  ".join(subtitle_parts), fontsize=8, color=PALETTE["grey"], pad=3)

    # Leaf summary table
    ax_tbl = fig.add_subplot(gs[1])
    ax_tbl.axis("off")

    col_labels = ["Leaf", "Recommended param", "#Instances", "Solve%",
                  "Avg time (s)", "vs baseline"]
    bp_idx = results.psettingByName[args.baseline].idx if args.baseline in results.psettingByName else None
    table_rows = []
    for i, leaf in enumerate(tree.leafs):
        ps = leaf.bestPS[0]
        times = [inst.results[ps.idx] for inst in leaf.instances]
        n = len(leaf.instances)
        n_solved = sum(1 for t in times if t < penalty * 0.99)
        avg_t = sum(times) / n
        if bp_idx is not None:
            bl_times = [inst.results[bp_idx] for inst in leaf.instances]
            bl_avg = sum(bl_times) / n
            sp_str = f"{bl_avg/avg_t:.2f}×" if avg_t > 0 else "—"
        else:
            sp_str = "—"
        table_rows.append([
            str(i + 1),
            ps.setting,
            str(n),
            f"{100*n_solved/n:.1f}%",
            f"{avg_t:.1f}",
            sp_str,
        ])

    tbl = ax_tbl.table(
        cellText=table_rows,
        colLabels=col_labels,
        cellLoc="center", loc="center",
        bbox=[0, 0, 1, 1])
    tbl.auto_set_font_size(False)
    tbl.set_fontsize(7.5)
    tbl.auto_set_column_width(list(range(len(col_labels))))
    for j in range(len(col_labels)):
        tbl[0, j].set_facecolor(PALETTE["blue"])
        tbl[0, j].set_text_props(color="white", fontweight="bold")
    # Left-align the param column
    for row_i in range(1, len(table_rows) + 1):
        tbl[row_i, 1].set_text_props(ha="left")

    add_page_footer(fig, exp_name, page_num)
    pdf.savefig(fig, bbox_inches="tight")
    plt.close(fig)


def page_discoveries(pdf, exp_name, time_piv, solved_piv, penalty, args,
                     greedy_curve, exhaustive_results, page_num=6):
    """Page 5 — Summary: key discoveries in plain text + key numbers."""
    fig = plt.figure(figsize=(11, 8.5))
    page_title(fig, "Key Discoveries & Recommendations", exp_name)

    shift     = args.shift
    params    = list(time_piv.columns)
    instances = list(time_piv.index)
    n_inst    = len(instances)

    sgm_all = {p: shifted_geomean(list(time_piv[p]), shift) for p in params}
    ranked  = sorted(params, key=lambda p: sgm_all[p])
    base_sgm = sgm_all.get(args.baseline, float("nan"))
    base_sr  = solved_piv[args.baseline].mean() * 100 if args.baseline in params else float("nan")

    # Build key facts
    best_single        = ranked[0]
    best_single_sgm    = sgm_all[best_single]
    best_single_sr     = solved_piv[best_single].mean() * 100
    best_single_vs_base = (base_sgm - best_single_sgm) / base_sgm * 100

    facts = []

    facts.append(("Baseline",
                  f"{args.baseline}  →  SGM={base_sgm:.1f}s,  solve rate={base_sr:.1f}%"))

    facts.append(("Best single param",
                  f"{best_single}  →  SGM={best_single_sgm:.1f}s,  "
                  f"solve rate={best_single_sr:.1f}%  "
                  f"({best_single_vs_base:+.1f}% vs baseline)"))

    if len(greedy_curve) >= 2:
        k2_sgm, k2_sr = greedy_curve[1][1], greedy_curve[1][2]
        p2 = greedy_curve[1][3]
        gain_k2 = (base_sgm - k2_sgm) / base_sgm * 100
        facts.append(("Racing K=2 (greedy)",
                      f"[{best_single}, {p2}]  →  SGM={k2_sgm:.1f}s,  "
                      f"solve rate={k2_sr:.1f}%  ({gain_k2:+.1f}% vs baseline)"))
    if len(greedy_curve) >= 3:
        k3_sgm, k3_sr = greedy_curve[2][1], greedy_curve[2][2]
        p3 = greedy_curve[2][3]
        gain_k3 = (base_sgm - k3_sgm) / base_sgm * 100
        port3 = [greedy_curve[i][3] for i in range(3)]
        facts.append(("Racing K=3 (greedy)",
                      f"{port3}  →  SGM={k3_sgm:.1f}s,  "
                      f"solve rate={k3_sr:.1f}%  ({gain_k3:+.1f}% vs baseline)"))

    if exhaustive_results:
        for k, (ex_sgm, ex_params) in sorted(exhaustive_results.items()):
            gain = (base_sgm - ex_sgm) / base_sgm * 100
            greedy_sgm_k = greedy_curve[k - 1][1] if k <= len(greedy_curve) else float("nan")
            greedy_diff = (greedy_sgm_k - ex_sgm) / ex_sgm * 100 if not math.isnan(greedy_sgm_k) else float("nan")
            greedy_note = (f"greedy ≈ optimal" if abs(greedy_diff) < 0.01
                           else f"greedy +{greedy_diff:.2f}% suboptimal")
            facts.append((f"Optimal K={k} (exhaustive)",
                          f"{ex_params}  →  SGM={ex_sgm:.2f}s  "
                          f"({gain:+.1f}% vs baseline,  {greedy_note})"))

    # Diminishing returns
    if len(greedy_curve) >= 3:
        gains = []
        prev = greedy_curve[0][1]
        for k, sgm, sr, p_added in greedy_curve[1:]:
            gains.append((k, (prev - sgm) / prev * 100))
            prev = sgm
        dr_k = next((k for k, g in gains if g < 5.0), None)
        if dr_k:
            facts.append(("Diminishing returns",
                          f"ΔStep% drops below 5% at K={dr_k}. "
                          f"K={dr_k - 1} is likely the sweet spot."))

    # Hard instances
    any_solved = solved_piv.any(axis=1)
    n_hard = (~any_solved).sum()
    facts.append(("Hard instances",
                  f"{n_hard} instance(s) not solved by any param within timelimit."))

    # Draw as styled text blocks
    ax = fig.add_axes([0.05, 0.05, 0.90, 0.85])
    ax.axis("off")

    y = 0.94
    line_h = 0.088
    for label, text in facts:
        # Label (bold, coloured)
        ax.text(0.0, y, f"▸ {label}:", fontsize=10, fontweight="bold",
                color=PALETTE["blue"], transform=ax.transAxes, va="top")
        y -= 0.032
        # Wrapped text
        import textwrap
        wrapped = textwrap.fill(text, width=100)
        ax.text(0.03, y, wrapped, fontsize=9, color="black",
                transform=ax.transAxes, va="top", family="monospace")
        y -= line_h + 0.004
        if y < 0.04:
            break

    add_page_footer(fig, exp_name, page_num)
    pdf.savefig(fig, bbox_inches="tight")
    plt.close(fig)


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("exp_dir", help="Experiment directory")
    parser.add_argument("output", nargs="?", default=None,
                        help="Output PDF path (default: <exp_dir>/lp_report.pdf)")
    parser.add_argument("--timelimit",    type=float, default=7200.0,
                        help="LP time limit used in experiment (default: 7200)")
    parser.add_argument("--penalty-mult", type=float, default=2.0)
    parser.add_argument("--shift",        type=float, default=1.0)
    parser.add_argument("--baseline",     default="dual_default")
    parser.add_argument("--top",          type=int,   default=8)
    parser.add_argument("--exhaustive-k", type=int,   default=4)
    parser.add_argument("--max-k",        type=int,   default=8)
    parser.add_argument("--dtree-depth",  type=int,   default=3,
                        help="Decision tree max depth (default: 3)")
    parser.add_argument("--dtree-min-leaf", type=int, default=20,
                        help="Decision tree minimum instances per leaf (default: 20)")
    parser.add_argument("--no-dtree",     action="store_true",
                        help="Skip the decision tree page")
    parser.add_argument("--features",
                        default=os.path.expanduser("~/inst/miplib/2017+spp/features.csv"),
                        help="Instance features CSV for decision tree")
    parser.add_argument("--exclude-prefix", default="guess_",
                        help="Exclude params starting with this prefix (default: guess_)")
    args = parser.parse_args()

    exp_dir = args.exp_dir
    penalty = args.timelimit * args.penalty_mult
    out_pdf = args.output or os.path.join(exp_dir, "lp_report.pdf")
    exp_name = os.path.basename(os.path.normpath(exp_dir))

    # ── Load data ─────────────────────────────────────────────────────────────
    result = load_avg_times(exp_dir, penalty, args.exclude_prefix)
    if result is None or result[0] is None:
        sys.exit("Error: lp_avg_times.csv not found. Run analyze_lp_params.py first.")
    df, time_piv, solved_piv, instances, params = result

    raw_df = load_raw_results(exp_dir)

    print(f"Loaded: {len(instances)} instances, {len(params)} params", file=sys.stderr)

    # ── Compute portfolio analysis ────────────────────────────────────────────
    time_np   = time_piv.values
    solved_np = solved_piv.values

    print("Building greedy portfolio...", file=sys.stderr, flush=True)
    greedy_curve, greedy_port_idx = greedy_portfolio(
        time_np, solved_np, params, args.shift, args.max_k)

    print(f"Running exhaustive search K≤{args.exhaustive_k}...", file=sys.stderr, flush=True)
    exhaustive_results = exhaustive_optimal(
        time_np, solved_np, params, args.shift, args.exhaustive_k)

    # ── Generate PDF ──────────────────────────────────────────────────────────
    print(f"Generating {out_pdf} ...", file=sys.stderr, flush=True)
    with PdfPages(out_pdf) as pdf:
        page_overview(pdf, exp_name, raw_df, time_piv, solved_piv,
                      penalty, args.timelimit, args)
        page_param_rankings(pdf, exp_name, time_piv, solved_piv, penalty, args)
        page_performance_profiles(pdf, exp_name, time_piv, solved_piv, penalty, args)
        page_racing_portfolios(pdf, exp_name, time_piv, solved_piv, args,
                               greedy_curve, exhaustive_results)
        page_num = 5
        if not args.no_dtree:
            if not _FBPS_AVAILABLE:
                print("  Warning: fbps library not available, skipping decision tree page.",
                      file=sys.stderr)
            else:
                print("Building decision tree...", file=sys.stderr, flush=True)
                page_decision_tree(pdf, exp_name, exp_dir, penalty, args, page_num)
                page_num += 1
        page_discoveries(pdf, exp_name, time_piv, solved_piv, penalty, args,
                         greedy_curve, exhaustive_results, page_num)

    print(f"✓ Report saved to {out_pdf}", file=sys.stderr)
    print(out_pdf)


if __name__ == "__main__":
    main()
