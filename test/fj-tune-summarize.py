#!/usr/bin/env python3
"""fj-tune-summarize.py — aggregate per-config TSVs from fj-tune-sweep.sh into
a single summary table: how often a root-only (or N-node) replay finds *any*
incumbent, and how good it is relative to bks.tsv.

Usage: fj-tune-summarize.py <outdir> <configs.tsv>

Reads <outdir>/<config>.tsv for every config named in configs.tsv (columns:
instance, found, best, bound, bks, bbtime, optimal) and prints one summary
row per config to stdout (tab-separated).
"""
import sys
import os
import math


def read_configs(path):
    names = []
    with open(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            name = line.split("\t")[0]
            names.append(name)
    return names


def fnum(s):
    try:
        v = float(s)
        if math.isnan(v) or math.isinf(v):
            return None
        return v
    except (ValueError, TypeError):
        return None


def summarize(path):
    n = 0
    n_found = 0
    n_optimal = 0
    gaps = []
    times = []
    with open(path) as f:
        header = f.readline()
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 7:
                continue
            _inst, found, best, bound, bks, bbtime, optimal = parts[:7]
            n += 1
            found = found == "1"
            if found:
                n_found += 1
            if optimal == "1":
                n_optimal += 1
            t = fnum(bbtime)
            if t is not None:
                times.append(t)
            bv = fnum(best)
            kv = fnum(bks)
            if found and bv is not None and kv is not None:
                # Normalized primal gap |obj-bks| / max(|obj|,|bks|,eps), the
                # same convention compare_benchmarks.py uses elsewhere in this
                # repo -- bounded to [0,100]%, unlike a pure relative-to-bks
                # gap which blows up on near-zero-BKS instances (e.g. the
                # markshare family, bks in {1,14}; mushroom-best, bks~0.055).
                denom = max(abs(bv), abs(kv), 1e-9)
                gaps.append(100.0 * abs(bv - kv) / denom)
    return {
        "n": n,
        "n_found": n_found,
        "n_optimal": n_optimal,
        "avg_gap_bks_pct": (sum(gaps) / len(gaps)) if gaps else float("nan"),
        "n_gap_samples": len(gaps),
        "avg_bbtime_s": (sum(times) / len(times)) if times else float("nan"),
        "total_bbtime_s": sum(times) if times else 0.0,
    }


def main():
    if len(sys.argv) != 3:
        print(__doc__, file=sys.stderr)
        return 2
    outdir, configs_path = sys.argv[1], sys.argv[2]
    names = read_configs(configs_path)
    cols = ["config", "n", "found", "found_pct", "optimal",
            "avg_gap_bks_pct", "n_gap_samples", "avg_bbtime_s", "total_bbtime_s"]
    print("\t".join(cols))
    for name in names:
        path = os.path.join(outdir, name + ".tsv")
        if not os.path.exists(path):
            continue
        s = summarize(path)
        found_pct = 100.0 * s["n_found"] / s["n"] if s["n"] else float("nan")
        print("\t".join(str(x) for x in [
            name, s["n"], s["n_found"], f"{found_pct:.2f}", s["n_optimal"],
            f"{s['avg_gap_bks_pct']:.4f}" if not math.isnan(s["avg_gap_bks_pct"]) else "nan",
            s["n_gap_samples"],
            f"{s['avg_bbtime_s']:.4f}" if not math.isnan(s["avg_bbtime_s"]) else "nan",
            f"{s['total_bbtime_s']:.2f}",
        ]))
    return 0


if __name__ == "__main__":
    sys.exit(main())
