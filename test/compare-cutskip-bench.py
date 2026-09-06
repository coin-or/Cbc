#!/usr/bin/env python3
"""Compare two mip-root-replay benchmark runs (e.g. base vs. adaptive
root cut-generator skip) produced by bench-adaptive-cutskip, using a
size-normalized "gap closed per second" efficiency metric instead of raw
bound deltas (which aren't comparable across instances of very different
tightness/scale).

Dual-bound metric (per instance, sense-generic -- works for both min and
max because both numerator and denominator flip sign together):

    dualGapClosed = clamp( (rootBound - lpBound) / (bks - lpBound), 0, 1 )

  lpBound = the *pre-cut* root LP relaxation bound, taken from the
            fixture's own <name>.root.meta ("objValue" line) -- this is
            not something the replay run at hand controls, so it's the
            correct common baseline for comparing two cut-generator
            configurations against each other. bks==lpBound (LP already
            tight) is defined as dualGapClosed=1 rather than a division
            by zero.

Primal-bound metric: the standard MIPLIB 2017 "primal gap" (no natural
"worst possible" anchor like the LP relaxation exists on the primal
side, so this uses the usual normalization by objective magnitude
instead of by an LP bound):

    primalGap = 1                                   if no solution or
                                                       sign(best) != sign(bks)
              = |best - bks| / max(|best|, |bks|, eps)   otherwise

  primalGapClosed = 1 - primalGap   (so higher is better, like the dual
  metric)

Efficiency (the metric requested): computed per instance as
gapClosed / bbTime, then aggregated with a *geometric mean* across
instances (shifted by +1 to tolerate zeros), not an arithmetic mean --
an arithmetic mean of ratios lets one very fast or very slow instance
dominate the average; geometric mean is the standard choice here (same
reasoning as MIPLIB/SCIP aggregate benchmark scores).

Usage:
    ./compare-cutskip-bench.py BASE.tsv OTHER.tsv \
        [--fixture-dir=PATH] [--bks=PATH] [--label-base=NAME] [--label-other=NAME]
"""
import argparse
import csv
import math
import os
import sys
from collections import namedtuple

EPS = 1e-10


def read_tsv(path):
    rows = {}
    with open(path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            rows[row["instance"]] = row
    return rows


def read_bks(path):
    bks = {}
    if not os.path.exists(path):
        return bks
    with open(path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            try:
                bks[row["instance"]] = (float(row["objective"]), row.get("sense", "min"))
            except (KeyError, ValueError):
                continue
    return bks


def read_lp_bound(fixture_dir, instance, tag="root"):
    meta_path = os.path.join(fixture_dir, f"{instance}.{tag}.meta")
    if not os.path.exists(meta_path):
        return None
    with open(meta_path) as f:
        for line in f:
            parts = line.split()
            if len(parts) == 2 and parts[0] == "objValue":
                return float(parts[1])
    return None


def to_float(s):
    if s is None:
        return None
    s = s.strip()
    if s in ("none", "n/a", "", "ERROR"):
        return None
    try:
        return float(s)
    except ValueError:
        return None


Metrics = namedtuple("Metrics", "instance dualGapClosed primalGapClosed bbTime "
                                 "dualEff primalEff bound best bks lpBound")


def compute_metrics(row, bks_map, fixture_dir, tag):
    inst = row["instance"]
    bbTime = to_float(row.get("bbTime"))
    bound = to_float(row.get("bound"))
    best = to_float(row.get("best"))
    if bbTime is None or bound is None:
        return None  # ERROR row (timeout/crash) -- excluded, reported separately

    bks_entry = bks_map.get(inst)
    bks_val = bks_entry[0] if bks_entry else None
    lp_bound = read_lp_bound(fixture_dir, inst, tag)

    dual_gap_closed = None
    if bks_val is not None and lp_bound is not None:
        denom = bks_val - lp_bound
        if abs(denom) < EPS:
            dual_gap_closed = 1.0
        else:
            dual_gap_closed = (bound - lp_bound) / denom
            dual_gap_closed = max(0.0, min(1.0, dual_gap_closed))

    primal_gap_closed = None
    if bks_val is not None:
        if best is None:
            primal_gap = 1.0
        elif (best >= 0) != (bks_val >= 0) and abs(bks_val) > EPS:
            primal_gap = 1.0
        else:
            primal_gap = abs(best - bks_val) / max(abs(best), abs(bks_val), EPS)
            primal_gap = min(1.0, primal_gap)
        primal_gap_closed = 1.0 - primal_gap

    dual_eff = (dual_gap_closed / bbTime) if (dual_gap_closed is not None and bbTime > 0) else None
    primal_eff = (primal_gap_closed / bbTime) if (primal_gap_closed is not None and bbTime > 0) else None

    return Metrics(inst, dual_gap_closed, primal_gap_closed, bbTime, dual_eff, primal_eff,
                    bound, best, bks_val, lp_bound)


def geomean(values):
    values = [v for v in values if v is not None]
    if not values:
        return None
    # Shift by +1 to tolerate exact zeros (a generator that closes 0% gap
    # in some time still contributes a finite, very small efficiency
    # rather than collapsing the whole geometric mean to zero).
    logs = [math.log(v + 1.0) for v in values]
    return math.exp(sum(logs) / len(logs)) - 1.0


def mean(values):
    values = [v for v in values if v is not None]
    if not values:
        return None
    return sum(values) / len(values)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("base_tsv")
    ap.add_argument("other_tsv")
    ap.add_argument("--fixture-dir", default=os.path.expanduser(
        "~/instances/mip-sanity-data/rootFixtures"))
    ap.add_argument("--bks", default=os.path.join(
        os.path.dirname(os.path.abspath(__file__)), "mip-sanity-data", "bks.tsv"))
    ap.add_argument("--tag", default="root")
    ap.add_argument("--label-base", default="base")
    ap.add_argument("--label-other", default="other")
    ap.add_argument("--top-n", type=int, default=15)
    args = ap.parse_args()

    base_rows = read_tsv(args.base_tsv)
    other_rows = read_tsv(args.other_tsv)
    bks_map = read_bks(args.bks)

    common = sorted(set(base_rows) & set(other_rows))
    only_base = sorted(set(base_rows) - set(other_rows))
    only_other = sorted(set(other_rows) - set(base_rows))

    base_metrics = {}
    other_metrics = {}
    base_errors = []
    other_errors = []
    for inst in common:
        bm = compute_metrics(base_rows[inst], bks_map, args.fixture_dir, args.tag)
        om = compute_metrics(other_rows[inst], bks_map, args.fixture_dir, args.tag)
        if bm is None:
            base_errors.append(inst)
        else:
            base_metrics[inst] = bm
        if om is None:
            other_errors.append(inst)
        else:
            other_metrics[inst] = om

    both_ok = sorted(set(base_metrics) & set(other_metrics))

    def col(d, name, field):
        return [getattr(d[i], field) for i in name if getattr(d[i], field) is not None]

    print(f"Instances compared: {len(common)}  (base-only: {len(only_base)}, "
          f"other-only: {len(only_other)})")
    print(f"Errors (excluded from aggregates): base={len(base_errors)} "
          f"other={len(other_errors)}  both-ok={len(both_ok)}")
    print()

    header = f"{'':22s} {args.label_base:>14s} {args.label_other:>14s} {'delta':>10s}"
    print(header)
    print("-" * len(header))

    def report(label, base_vals, other_vals, pct=True, agg=mean):
        b = agg(base_vals)
        o = agg(other_vals)
        if b is None or o is None:
            print(f"{label:22s} {'n/a':>14s} {'n/a':>14s} {'n/a':>10s}")
            return
        d = o - b
        fmt = "{:.4f}" if not pct else "{:.2%}"
        print(f"{label:22s} {fmt.format(b):>14s} {fmt.format(o):>14s} {fmt.format(d):>10s}")

    report("Mean dual gap closed", col(base_metrics, both_ok, "dualGapClosed"),
           col(other_metrics, both_ok, "dualGapClosed"))
    report("Mean primal gap closed", col(base_metrics, both_ok, "primalGapClosed"),
           col(other_metrics, both_ok, "primalGapClosed"))
    report("Mean bbTime (s)", col(base_metrics, both_ok, "bbTime"),
           col(other_metrics, both_ok, "bbTime"), pct=False)
    report("Geomean dual eff (%/s)", col(base_metrics, both_ok, "dualEff"),
           col(other_metrics, both_ok, "dualEff"), agg=geomean)
    report("Geomean primal eff (%/s)", col(base_metrics, both_ok, "primalEff"),
           col(other_metrics, both_ok, "primalEff"), agg=geomean)

    print()
    print(f"Top {args.top_n} instances by dual-efficiency delta (other - base):")
    deltas = []
    for inst in both_ok:
        be = base_metrics[inst].dualEff
        oe = other_metrics[inst].dualEff
        if be is None or oe is None:
            continue
        deltas.append((oe - be, inst, be, oe))
    deltas.sort(reverse=True)
    print(f"{'instance':40s} {'base eff':>12s} {'other eff':>12s} {'delta':>12s}")
    for d, inst, be, oe in deltas[:args.top_n]:
        print(f"{inst:40s} {be:12.4f} {oe:12.4f} {d:12.4f}")
    print("... worst:")
    for d, inst, be, oe in deltas[-args.top_n:]:
        print(f"{inst:40s} {be:12.4f} {oe:12.4f} {d:12.4f}")


if __name__ == "__main__":
    main()
