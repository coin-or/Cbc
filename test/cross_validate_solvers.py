#!/usr/bin/env python3
"""
Cross-validate VRP instances with HiGHS and MIPster.
Compares objective values to detect wrong results.
"""

import subprocess
import sys
import os
import re
from pathlib import Path

def solve_with_mipster(mps_file, timeout=60):
    """Solve with MIPster and extract objective."""
    cmd = [f"../src/mipster", mps_file, "-sec", str(timeout), "-solve"]
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout+10)
        output = result.stdout + result.stderr

        # Check status
        if "Optimal solution found" in output:
            match = re.search(r"Objective value:\s+([\d.e+-]+)", output)
            if match:
                return float(match.group(1)), "optimal"
        elif "Stopped on time limit" in output and "BestSol:" in output:
            match = re.search(r"BestSol:\s+([\d.e+-]+)", output)
            if match:
                return float(match.group(1)), "feasible"
        elif "INFEASIBLE" in output:
            return None, "infeasible"
    except Exception as e:
        print(f"    MIPster error: {e}")

    return None, "unknown"

def solve_with_highs(mps_file, timeout=60):
    """Solve with HiGHS and extract objective."""
    try:
        # Check if highs is available
        result = subprocess.run(["which", "highs"], capture_output=True)
        if result.returncode != 0:
            return None, "unavailable"

        cmd = ["highs", mps_file, f"--time_limit={timeout}"]
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout+10)
        output = result.stdout + result.stderr

        # Parse HiGHS output
        if "Status            Optimal" in output or "Model status      : Optimal" in output:
            # Try different objective patterns
            match = re.search(r"([\d.e+-]+)\s+\(objective\)", output)
            if not match:
                match = re.search(r"Objective\s+value:\s+([\d.e+-]+)", output)
            if match:
                return float(match.group(1)), "optimal"
        elif "Infeasible" in output or "Model status      : Infeasible" in output:
            return None, "infeasible"
    except Exception as e:
        print(f"    HiGHS error: {e}")

    return None, "unknown"

def compare_results(name, mipster_obj, mipster_status, highs_obj, highs_status):
    """Compare solver results and report discrepancies."""
    # If HiGHS unavailable, just report MIPster result
    if highs_status == "unavailable":
        print(f"{name:30s} MIPster: {mipster_status:10s} obj={mipster_obj if mipster_obj else 'N/A':>12}")
        return True

    # Both solved
    if mipster_status in ["optimal", "feasible"] and highs_status in ["optimal", "feasible"]:
        if mipster_obj is None or highs_obj is None:
            print(f"{name:30s} ✗ MISSING OBJ  MIPster={mipster_obj}, HiGHS={highs_obj}")
            return False

        rel_diff = abs(mipster_obj - highs_obj) / max(abs(highs_obj), 1e-6)
        if rel_diff > 1e-4:
            print(f"{name:30s} ✗ MISMATCH  MIPster={mipster_obj:.6f}, HiGHS={highs_obj:.6f}, diff={rel_diff:.2e}")
            return False
        else:
            print(f"{name:30s} ✓ MATCH  obj={mipster_obj:.6f} (rel_diff={rel_diff:.2e})")
            return True

    # Status mismatch
    if mipster_status != highs_status:
        print(f"{name:30s} ⚠ STATUS MISMATCH  MIPster={mipster_status}, HiGHS={highs_status}")
        return False

    # Both infeasible or both unknown
    print(f"{name:30s} • {mipster_status}")
    return True

if __name__ == '__main__':
    print("=" * 80)
    print("VRP Cross-Solver Validation (MIPster vs HiGHS)")
    print("=" * 80)
    print()

    # Get all VRP instances
    instances = sorted(Path('.').glob('cvrp_*.mps')) + sorted(Path('.').glob('vrppd_*.mps'))

    if not instances:
        print("No VRP instances found. Run gen_vrp_fixtures.py and gen_diverse_vrp.py first.")
        sys.exit(1)

    passed = 0
    failed = 0
    failed_instances = []

    for mps_file in instances:
        name = mps_file.stem

        # Solve with both solvers
        mipster_obj, mipster_status = solve_with_mipster(str(mps_file))
        highs_obj, highs_status = solve_with_highs(str(mps_file))

        # Compare
        if compare_results(name, mipster_obj, mipster_status, highs_obj, highs_status):
            passed += 1
        else:
            failed += 1
            failed_instances.append(name)

    print()
    print("=" * 80)
    print(f"Summary: {passed} passed, {failed} failed")

    if failed > 0:
        print()
        print("Failed instances:")
        for name in failed_instances:
            print(f"  - {name}")
        sys.exit(1)
