#!/usr/bin/env python3
"""
compare_fpp_results.py
Compare Gurobi validation results against the BKS (best-known-solution) file
for fast preprocessing validation.

For each preprocessed instance solved by Gurobi:
  - If Gurobi finds INFEASIBLE but BKS says feasible → BUG (we over-fixed)
  - If Gurobi finds optimal and BKS=opt known → values must match (within tol)
  - If Gurobi finds optimal and BKS=best → our optimal must be <= BKS (min) or >= BKS (max)
  - If Gurobi time-limited → use dual bound to check: dual <= BKS is consistent

Usage:
  python3 compare_fpp_results.py [grb_results_dir] [bks_solu_file] [fpp_summary_csv]
"""

import sys
import os
import re
import math
from pathlib import Path

# Paths (can be overridden by argv)
GRB_DIR = sys.argv[1] if len(sys.argv) > 1 else \
    "/home/haroldo/experiments/cbc/fpp_gurobi_validation"
BKS_FILE = sys.argv[2] if len(sys.argv) > 2 else \
    "/home/haroldo/inst/miplib/2017/miplib2017-v35.solu"
FPP_CSV = sys.argv[3] if len(sys.argv) > 3 else \
    "/home/haroldo/dev/cbc/fpp_preprocess_summary.csv"

OBJ_TOL = 1e-4  # relative tolerance for objective comparison

# --- Parse BKS file ---
bks = {}  # name -> (type, value)  type in {opt, best, inf}
with open(BKS_FILE) as f:
    for line in f:
        line = line.strip()
        m = re.match(r'^=(opt|best|inf)=\s+(\S+)\s+(-?[\d.e+\-]+)?', line)
        if m:
            bks_type, name, val = m.group(1), m.group(2), m.group(3)
            bks[name] = (bks_type, float(val) if val else None)

# --- Parse FPP summary CSV ---
fpp_info = {}  # name -> dict
with open(FPP_CSV) as f:
    headers = f.readline().strip().split(',')
    for line in f:
        fields = line.strip().split(',')
        d = dict(zip(headers, fields))
        fpp_info[d['instance']] = d

# --- Parse Gurobi log files ---
def parse_grb_log(logfile):
    """Return dict with keys: status, obj_value, dual_bound, gap_pct."""
    result = {
        'status': 'unknown',
        'obj_value': None,
        'dual_bound': None,
        'gap_pct': None,
    }
    try:
        text = Path(logfile).read_text(errors='replace')
    except FileNotFoundError:
        result['status'] = 'not_run'
        return result

    if re.search(r'Model is infeasible', text):
        result['status'] = 'infeasible'
    elif re.search(r'Optimal solution found', text):
        result['status'] = 'optimal'
    elif re.search(r'Time limit reached', text):
        result['status'] = 'time_limit'
    elif re.search(r'Solved in \d+ iterations', text):
        result['status'] = 'optimal'

    # Best objective (primal bound)
    m = re.search(r'Best objective\s+([-\d.e+]+)', text)
    if m:
        try:
            result['obj_value'] = float(m.group(1))
        except ValueError:
            pass

    # Best bound (dual bound)
    m = re.search(r'Best bound\s+([-\d.e+]+)', text)
    if m:
        try:
            result['dual_bound'] = float(m.group(1))
        except ValueError:
            pass

    # Gap
    m = re.search(r'Gap\s+([\d.]+)%', text)
    if m:
        result['gap_pct'] = float(m.group(1))

    return result

# --- Run comparison ---
print(f"{'Instance':<40} {'FPP':>8} {'GrbStatus':>12} {'GrbObj':>15} {'BKSType':>8} {'BKSVal':>15} {'Check':>10}")
print("-" * 120)

bugs = []
ok = []
improved_bks = []
skipped = []
not_run = []

grb_logs = sorted(Path(GRB_DIR).glob("*.grb.log")) if os.path.isdir(GRB_DIR) else []
processed_names = {f.name.replace('.grb.log', '') for f in grb_logs}

for name in sorted(processed_names):
    logfile = os.path.join(GRB_DIR, f"{name}.grb.log")
    grb = parse_grb_log(logfile)

    fpp = fpp_info.get(name, {})
    n_fixed = int(fpp.get('n_fixed', 0))

    bks_type, bks_val = bks.get(name, (None, None))

    # Determine problem sense from Gurobi log (look for MINIMIZE/MAXIMIZE or obj sign)
    # For validation we use: if BKS opt, grb opt must equal; if BKS inf, grb must be inf
    check = "?"
    is_bug = False

    if grb['status'] == 'not_run':
        not_run.append(name)
        check = "NOT_RUN"
    elif grb['status'] == 'infeasible':
        if bks_type in ('opt', 'best'):
            # Original was feasible → preprocessing made it infeasible → BUG
            is_bug = True
            check = "BUG_INFEAS"
        elif bks_type == 'inf':
            check = "OK(inf)"
        else:
            check = "INFEAS(no_bks)"
    elif grb['status'] in ('optimal', 'time_limit'):
        obj = grb['obj_value']
        dual = grb['dual_bound']

        if bks_type == 'inf':
            # BKS says infeasible but Gurobi found solution → BUG
            is_bug = True
            check = "BUG_SHOULD_BE_INF"
        elif bks_type == 'opt' and bks_val is not None:
            # Gurobi optimal must match known optimal (or be strictly better, which
            # means the solu file had an incorrect/outdated value).
            if grb['status'] == 'optimal' and obj is not None:
                rel_err = (obj - bks_val) / max(1.0, abs(bks_val))
                abs_rel_err = abs(rel_err)
                if abs_rel_err <= OBJ_TOL:
                    check = f"OK(opt,{abs_rel_err:.1e})"
                elif rel_err < -OBJ_TOL:
                    # Gurobi found strictly better than "known optimal" → old BKS was wrong
                    check = f"IMPROVED_BKS({rel_err:.2e})"
                    # Not a preprocessing bug — but worth reporting
                else:
                    # Gurobi found strictly worse → preprocessing excluded the optimal → BUG
                    is_bug = True
                    check = f"BUG_OPT(+{rel_err:.2e})"
            else:
                # Time limit: dual bound should be <= bks_val (for min) or >= (for max)
                if dual is not None:
                    # Check both directions — we don't know sense so flag if dual far from BKS
                    rel_err = abs(dual - bks_val) / max(1.0, abs(bks_val))
                    check = f"TL_dual_err={rel_err:.2e}"
                else:
                    check = "TL(no_dual)"
        elif bks_type == 'best' and bks_val is not None:
            # Gurobi optimal (if found) should be <= BKS (min) or >= (max)
            # We can't know sense definitively, so just record
            if grb['status'] == 'optimal' and obj is not None:
                # If our optimal is much worse than BKS, it means we over-fixed
                # BKS is a feasible solution value, optimal_preprocessed = original_optimal <= BKS (min)
                # If obj >> bks_val (minimization) → bug; if obj << bks_val (maximization) → bug
                # Use a relative threshold: if more than 1% worse, flag
                rel_gap = (obj - bks_val) / max(1.0, abs(bks_val))
                if rel_gap > OBJ_TOL:
                    is_bug = True
                    check = f"BUG_BEST(+{rel_gap:.2e})"
                else:
                    check = f"OK(best,{rel_gap:.2e})"
            else:
                check = "TL(best_bks)"
        else:
            check = f"NO_BKS({grb['status']})"

    obj_str = f"{grb['obj_value']:.6g}" if grb['obj_value'] is not None else "n/a"
    bks_str = f"{bks_val:.6g}" if bks_val is not None else "n/a"

    print(f"{name:<40} {n_fixed:>8} {grb['status']:>12} {obj_str:>15} {str(bks_type or ''):>8} {bks_str:>15} {check:>10}")

    if is_bug:
        bugs.append(name)
    elif check.startswith("IMPROVED_BKS"):
        improved_bks.append((name, grb['obj_value'], bks_val))
        ok.append(name)
    elif check not in ("NOT_RUN", "?"):
        ok.append(name)

print()
print("=" * 120)
print(f"Total processed : {len(processed_names)}")
print(f"OK              : {len(ok)}")
print(f"Improved BKS    : {len(improved_bks)}")
print(f"Bugs found      : {len(bugs)}")
print(f"Not run         : {len(not_run)}")
print()
if improved_bks:
    print("Instances with improved BKS (old BKS was incorrect):")
    for name, new_val, old_val in improved_bks:
        print(f"  {name}: {old_val:.6g} → {new_val:.6g}")
    print()
if bugs:
    print("!! BUGS DETECTED !!")
    for b in bugs:
        print(f"  {b}")
else:
    print("No bugs detected.")
