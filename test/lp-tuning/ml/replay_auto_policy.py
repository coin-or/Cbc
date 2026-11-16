#!/usr/bin/env python3
"""Replay the compiled sequential auto scorer against its historical dataset.

Requires a C++ compiler and installed Osi headers (pkg-config osi, or
--include-dir). No Python packages or LP solves are required. Timings are frozen
historical measurements, not measurements of the current executable.
"""

import argparse
import csv
import ctypes
import json
import math
import os
from pathlib import Path
import shlex
import subprocess
import tempfile


FEATURES = {
    "rNzPercRowsCovering": "OFnzPercRowsCovering",
    "colsLeast4Nz": "OFcolsLeast4Nz",
    "percRowsLeast4096Nz": "OFpercRowsLeast4096Nz",
    "rNzPercRowsHubImplication": "OFnzPercRowsHubImplication",
    "percBin": "OFpercBin",
    "cols": "OFcols",
    "rNzpercRowsPartPacking": "OFnzpercRowsPartPacking",
    "rowsLeast8Nz": "OFrowsLeast8Nz",
    "objStdDev": "OFobjStdDev",
    "rNzrowsPartPacking": "OFnzrowsPartPacking",
}


def main():
    repo = Path(__file__).resolve().parents[3]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", type=Path,
                        default=Path(__file__).with_name("auto_policy_results.csv"))
    parser.add_argument("--include-dir", type=Path)
    args = parser.parse_args()
    include = (["-I" + str(args.include_dir)] if args.include_dir else
               shlex.split(subprocess.check_output(
                   ["pkg-config", "--cflags", "osi"], text=True)))
    rows = list(csv.DictReader(args.data.open()))
    if not rows or len({r["instance"] for r in rows}) != len(rows):
        raise ValueError("Expected nonempty, unique instance rows")
    eligible = [r for r in rows if r["rule"] != "bypass"]
    baseline, policy, selected = [], [], []
    failures_before = failures_after = new_failures = 0
    counts = [0] * 5
    with tempfile.TemporaryDirectory(prefix="cbc-auto-replay-") as directory:
        directory = Path(directory)
        wrapper = directory / "wrapper.cpp"
        wrapper.write_text(
            '#include "CbcLpParamScorer.hpp"\n#include "OsiFeatures.hpp"\n'
            'extern "C" int featureCount() { return OFCount; }\n'
            'extern "C" int featureIndex(int i) { const int indices[] = {'
            + ",".join(FEATURES.values()) + '}; return indices[i]; }\n'
            'extern "C" const char *recommend(const double *x) {'
            'return cbcRecommendLpParam(x); }\n')
        library = directory / "scorer.so"
        subprocess.run(
            shlex.split(os.environ.get("CXX", "c++"))
            + ["-std=c++11", "-shared", "-fPIC", "-I" + str(repo / "src")]
            + include + [str(wrapper), str(repo / "src/CbcLpParamScorer.cpp"),
                         "-o", str(library)], check=True)
        scorer = ctypes.CDLL(str(library))
        scorer.featureCount.restype = ctypes.c_int
        scorer.featureIndex.argtypes = [ctypes.c_int]
        scorer.featureIndex.restype = ctypes.c_int
        scorer.recommend.argtypes = [ctypes.POINTER(ctypes.c_double)]
        scorer.recommend.restype = ctypes.c_char_p
        vector_type = ctypes.c_double * scorer.featureCount()
        for row in eligible:
            vector = vector_type()
            for i, name in enumerate(FEATURES):
                vector[scorer.featureIndex(i)] = float(row[name])
            actual = scorer.recommend(vector).decode()
            if actual != row["policy_tag"]:
                raise ValueError(f"{row['instance']}: {actual} != {row['policy_tag']}")
            counts[int(row["rule"])] += 1
            before, after = [], []
            for seed in (1, 2, 3):
                before.append(float(row[f"baseline_s{seed}"]))
                after.append(float(row[f"policy_s{seed}"]))
                b = int(row[f"baseline_valid_s{seed}"])
                a = int(row[f"policy_valid_s{seed}"])
                if b not in (0, 1) or a not in (0, 1):
                    raise ValueError("Invalid validation flag")
                failures_before += 1 - b
                failures_after += 1 - a
                new_failures += b and not a
            if any(not math.isfinite(t) or t < 0 for t in before + after):
                raise ValueError("Invalid penalized runtime")
            baseline.append(sum(before) / 3)
            policy.append(sum(after) / 3)
            selected.append(row["rule"] != "0")
    wins = sum(s and (b + 1) / (p + 1) >= 2 and b - p >= 1
               for b, p, s in zip(baseline, policy, selected))
    regressions = sum(s and (p + 1) / (b + 1) > 1.1 and p - b >= 1
                      for b, p, s in zip(baseline, policy, selected))
    sgm = lambda times: math.expm1(sum(map(math.log1p, times)) / len(times))
    print(json.dumps({
        "instances": len(eligible), "preselection_bypasses": len(rows) - len(eligible),
        "fallback_and_rule_counts": counts,
        "big_wins": wins, "significant_regressions": regressions,
        "new_failed_seeds": new_failures,
        "failed_seeds_baseline": failures_before, "failed_seeds_policy": failures_after,
        "baseline_sgm_s": sgm(baseline), "policy_sgm_s": sgm(policy),
        "baseline_penalized_s": sum(baseline), "policy_penalized_s": sum(policy),
        "saved_penalized_s": sum(baseline) - sum(policy),
    }, indent=2))


if __name__ == "__main__":
    main()
