#!/usr/bin/env python3
"""summarize_lp_results.py — Cross-validate LP relaxation experiment results.

Reads lp_results.csv, determines the best feasible LP objective per instance,
and flags runs that reported OPTIMAL but whose objective is worse than the
best known feasible solution (fake optimals).

Usage:
    python3 summarize_lp_results.py <experiment_dir> [--tol 1e-6] [--csv] [--verbose]

Output:
    - Per-instance summary: best obj, number of runs, fake optimals
    - Per-run detail with FAKE_OPTIMAL / OK / INFEASIBLE / TIMEOUT tags
    - Writes <experiment_dir>/lp_summary.csv with annotated results
"""

import argparse
import csv
import os
import sys
from collections import defaultdict


def parse_check_result(check_str):
    """Parse 'yes;primal=1.07e-14;dual=3.13e-13' into (feasible, primal_err, dual_err)."""
    if not check_str or check_str == "NA":
        return None, None, None
    parts = check_str.split(";")
    feasible = parts[0].strip() if parts else None
    primal_err = None
    dual_err = None
    for p in parts[1:]:
        p = p.strip()
        if p.startswith("primal="):
            try:
                primal_err = float(p.split("=", 1)[1])
            except ValueError:
                pass
        elif p.startswith("dual="):
            try:
                dual_err = float(p.split("=", 1)[1])
            except ValueError:
                pass
    return feasible, primal_err, dual_err


def check_optimal(check_str):
    """Value of the 'optimal=' field ('yes'/'no'), None if absent/NA (older CSVs)."""
    for p in (check_str or "").split(";"):
        p = p.strip()
        if p.startswith("optimal="):
            v = p.split("=", 1)[1]
            return None if v == "NA" else v
    return None


def check_valid(check_str):
    """Solution passed -checkSolution: primal feasible and not flagged non-optimal."""
    return (parse_check_result(check_str)[0] == "yes"
            and check_optimal(check_str) != "no")


def parse_obj(val):
    """Parse objective value, return float or None."""
    if not val or val == "NA":
        return None
    try:
        return float(val)
    except ValueError:
        return None


def main():
    parser = argparse.ArgumentParser(description="Summarize LP relaxation experiment results")
    parser.add_argument("exp_dir", help="Experiment directory containing lp_results.csv")
    parser.add_argument("--tol", type=float, default=1e-6,
                        help="Relative tolerance for comparing objectives (default: 1e-6)")
    parser.add_argument("--csv", action="store_true",
                        help="Only write CSV output, no console report")
    parser.add_argument("--verbose", "-v", action="store_true",
                        help="Show per-run details")
    args = parser.parse_args()

    csv_path = os.path.join(args.exp_dir, "lp_results.csv")
    if not os.path.isfile(csv_path):
        print(f"Error: {csv_path} not found", file=sys.stderr)
        sys.exit(1)

    # Read all rows
    rows = []
    with open(csv_path, newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            rows.append(row)

    if not rows:
        print("No results found in CSV", file=sys.stderr)
        sys.exit(1)

    # Group by instance
    by_instance = defaultdict(list)
    for row in rows:
        by_instance[row["instance"]].append(row)

    # For each instance, find the best feasible objective
    # "Valid" = check_result starts with "yes", optimal!=no, AND status is OPTIMAL
    # Best = lowest obj (CBC minimizes internally)
    instance_best = {}  # instance -> best_obj (float)
    for inst, inst_rows in by_instance.items():
        best = None
        for r in inst_rows:
            valid = check_valid(r.get("check_result", ""))
            obj = parse_obj(r.get("obj_from_check"))
            if obj is None:
                obj = parse_obj(r.get("obj_from_log"))
            if valid and obj is not None:
                if best is None or obj < best:
                    best = obj
        instance_best[inst] = best

    # Annotate each row
    annotated = []
    fake_count = 0
    total_optimal = 0
    for row in rows:
        inst = row["instance"]
        status = row.get("status", "")
        feasible, primal_err, dual_err = parse_check_result(row.get("check_result", ""))
        obj = parse_obj(row.get("obj_from_check"))
        if obj is None:
            obj = parse_obj(row.get("obj_from_log"))
        best = instance_best.get(inst)

        if status == "OPTIMAL":
            total_optimal += 1
            if feasible != "yes":
                tag = "FAKE_OPTIMAL(infeasible)"
                fake_count += 1
            elif check_optimal(row.get("check_result", "")) == "no":
                tag = "FAKE_OPTIMAL(dual_infeasible)"
                fake_count += 1
            elif best is not None and obj is not None:
                # Check if this obj is worse than best known
                tol = args.tol * max(1.0, abs(best))
                if obj > best + tol:
                    tag = "FAKE_OPTIMAL"
                    fake_count += 1
                else:
                    tag = "OK"
            else:
                tag = "OK"
        elif status in ("TIMEOUT", "TIMEOUT_KILLED"):
            tag = "TIMEOUT"
        elif status == "INFEASIBLE":
            tag = "INFEASIBLE"
        elif status == "ERROR":
            tag = "ERROR"
        else:
            tag = status

        gap_to_best = ""
        if best is not None and obj is not None and abs(best) > 1e-10:
            gap_to_best = f"{(obj - best) / abs(best) * 100:.6f}%"
        elif best is not None and obj is not None:
            gap_to_best = f"{obj - best:.2e}"

        annotated.append({
            **row,
            "validation": tag,
            "best_obj": best if best is not None else "NA",
            "gap_to_best": gap_to_best if gap_to_best else "NA",
        })

    # Write annotated CSV
    out_csv = os.path.join(args.exp_dir, "lp_summary.csv")
    fieldnames = list(rows[0].keys()) + ["validation", "best_obj", "gap_to_best"]
    with open(out_csv, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(annotated)

    if args.csv:
        return

    # Console report
    n_instances = len(by_instance)
    n_runs = len(rows)
    n_fake = fake_count

    print(f"\n{'='*70}")
    print(f"  LP Relaxation Results Summary")
    print(f"{'─'*70}")
    print(f"  Experiment:     {args.exp_dir}")
    print(f"  Total runs:     {n_runs}")
    print(f"  Instances:      {n_instances}")
    print(f"  Optimal runs:   {total_optimal}")
    print(f"  Fake optimals:  {n_fake}")
    print(f"  Tolerance:      {args.tol}")
    print(f"{'='*70}\n")

    # Per-instance summary
    print(f"  {'Instance':<35} {'Runs':>5} {'Opt':>4} {'Fake':>5} {'Best Obj':>18} {'Feasible':>8}")
    print(f"  {'─'*35} {'─'*5} {'─'*4} {'─'*5} {'─'*18} {'─'*8}")

    for inst in sorted(by_instance.keys()):
        inst_rows = by_instance[inst]
        n = len(inst_rows)
        n_opt = sum(1 for r in inst_rows if r.get("status") == "OPTIMAL")
        n_feas = sum(1 for r in inst_rows
                     if check_valid(r.get("check_result", "")))
        inst_fakes = sum(1 for a in annotated
                         if a["instance"] == inst and "FAKE" in a["validation"])
        best = instance_best.get(inst)
        best_str = f"{best:.8g}" if best is not None else "NA"
        print(f"  {inst:<35} {n:>5} {n_opt:>4} {inst_fakes:>5} {best_str:>18} {n_feas:>8}")

    # Show fake optimals detail
    fakes = [a for a in annotated if "FAKE" in a["validation"]]
    if fakes:
        print(f"\n{'─'*70}")
        print(f"  FAKE OPTIMALS ({len(fakes)}):")
        print(f"{'─'*70}")
        for a in fakes:
            print(f"  {a['instance']:<30} {a['param_tag']:<25} s{a['seed']}")
            print(f"    obj={a.get('obj_from_check', 'NA'):>18}  "
                  f"best={a['best_obj']:>18}  "
                  f"gap={a['gap_to_best']:>12}  "
                  f"{a['validation']}")

    if args.verbose:
        print(f"\n{'─'*70}")
        print(f"  ALL RUNS:")
        print(f"{'─'*70}")
        for a in annotated:
            print(f"  {a['instance']:<30} {a['param_tag']:<20} s{a['seed']}  "
                  f"obj={a.get('obj_from_check', 'NA'):>16}  "
                  f"t={a['wall_seconds']:>8}s  "
                  f"{a['validation']}")

    print(f"\n  Summary CSV: {out_csv}")
    print()


if __name__ == "__main__":
    main()
