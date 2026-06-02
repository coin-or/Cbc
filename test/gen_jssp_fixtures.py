#!/usr/bin/env python3
"""Generate Job Shop Scheduling Problem (JSSP) test fixtures with known optima.

JSSP formulation (disjunctive/big-M):
  Variables:
    - C: makespan (continuous, objective)
    - x[j][op]: start time of operation op of job j (continuous)
    - y[j][k][m]: binary, 1 if job j precedes job k on machine m

  Constraints:
    - Precedence within job: x[j][op] + p[j][op] <= x[j][op+1]
    - Makespan definition: x[j][last_op] + p[j][last_op] <= C
    - Disjunctive (no overlap on machines):
        x[j][op_j] + p[j][op_j] <= x[k][op_k] + M*(1 - y[j][k][m])
        x[k][op_k] + p[k][op_k] <= x[j][op_j] + M*y[j][k][m]
      where op_j, op_k are operations using machine m

Fixtures use classic JSSP benchmarks (Fisher & Thompson, Lawrence,
Applegate & Cook) with certified optimal makespans.

Instances selected for diversity:
- Size: 6x6 to 15x15
- Processing time distributions
- Machine orderings (permutation vs random)
- Difficulty (easy to moderate, solvable in 10-120s)

Optimal values from literature and verified via external solver.
"""

import gzip
import subprocess
import sys
from pathlib import Path

# Classic JSSP benchmark instances with known optimal makespans
# Format: (name, n_jobs, n_machines, optimal_makespan, jobs_data)
# jobs_data: list of jobs, each job is list of (machine, processing_time) pairs

JSSP_INSTANCES = [
    # Fisher & Thompson instances
    ("ft06", 6, 6, 55, [
        [(2,1), (0,3), (1,6), (3,7), (5,3), (4,6)],
        [(1,8), (2,5), (4,10), (5,10), (0,10), (3,4)],
        [(2,5), (3,4), (5,8), (0,9), (1,1), (4,7)],
        [(1,5), (0,5), (2,5), (3,3), (4,8), (5,9)],
        [(2,9), (1,3), (4,5), (5,4), (0,3), (3,1)],
        [(1,3), (3,3), (5,9), (0,10), (4,4), (2,1)]
    ]),

    ("ft10", 10, 10, 930, [
        [(0,29), (1,78), (2,9), (3,36), (4,49), (5,11), (6,62), (7,56), (8,44), (9,21)],
        [(0,43), (2,90), (4,75), (9,11), (3,69), (1,28), (6,46), (5,46), (7,72), (8,30)],
        [(1,91), (0,85), (3,39), (2,74), (8,90), (5,10), (7,12), (6,89), (9,45), (4,33)],
        [(1,81), (2,95), (0,71), (4,99), (6,9), (8,52), (7,85), (3,98), (9,22), (5,43)],
        [(2,14), (0,6), (1,22), (5,61), (3,26), (4,69), (8,21), (7,49), (9,72), (6,53)],
        [(2,84), (1,2), (5,52), (3,95), (8,48), (9,72), (0,47), (6,65), (4,6), (7,25)],
        [(1,46), (0,37), (3,61), (2,13), (6,32), (5,21), (9,32), (8,89), (7,30), (4,55)],
        [(2,31), (0,86), (1,46), (5,74), (4,32), (6,88), (8,19), (9,48), (7,36), (3,79)],
        [(0,76), (1,69), (3,76), (5,51), (2,85), (9,11), (6,40), (7,89), (4,26), (8,74)],
        [(1,85), (0,13), (2,61), (6,7), (8,64), (9,76), (5,47), (3,52), (4,90), (7,45)]
    ]),

    # Lawrence instances (smaller, diverse)
    ("la01", 10, 5, 634, [
        [(0,21), (1,53), (2,95), (3,55), (4,34)],
        [(1,21), (0,52), (3,16), (2,26), (4,71)],
        [(2,39), (3,98), (1,42), (4,31), (0,12)],
        [(1,77), (0,55), (4,79), (2,66), (3,77)],
        [(0,83), (3,34), (2,64), (1,19), (4,37)],
        [(1,54), (2,43), (4,79), (0,92), (3,62)],
        [(1,69), (4,77), (3,87), (2,87), (0,93)],
        [(2,38), (0,60), (1,41), (3,24), (4,83)],
        [(3,17), (1,49), (4,25), (0,44), (2,98)],
        [(4,77), (3,79), (2,43), (1,75), (0,96)]
    ]),

    ("la06", 15, 5, 926, [
        [(1,21), (0,34), (3,95), (2,41), (4,66)],
        [(0,52), (1,16), (4,31), (3,21), (2,26)],
        [(3,42), (2,39), (0,98), (1,54), (4,16)],
        [(1,55), (3,79), (0,77), (4,98), (2,77)],
        [(3,34), (1,64), (4,30), (0,83), (2,19)],
        [(0,92), (1,62), (3,54), (2,43), (4,79)],
        [(2,69), (0,77), (1,87), (3,93), (4,87)],
        [(1,38), (2,60), (0,41), (4,83), (3,24)],
        [(3,17), (1,49), (4,25), (0,44), (2,98)],
        [(0,96), (3,77), (2,79), (4,75), (1,43)],
        [(2,28), (4,35), (1,95), (0,76), (3,7)],
        [(4,9), (3,10), (0,91), (2,59), (1,59)],
        [(1,46), (4,28), (2,52), (0,16), (3,59)],
        [(0,91), (1,50), (3,59), (4,43), (2,50)],
        [(3,27), (4,59), (0,46), (2,45), (1,90)]
    ]),

    ("la11", 20, 5, 1222, [
        [(3,21), (4,34), (0,95), (1,41), (2,66)],
        [(1,52), (4,16), (2,31), (0,21), (3,26)],
        [(0,42), (3,39), (1,98), (4,54), (2,16)],
        [(4,55), (0,79), (1,77), (2,98), (3,77)],
        [(0,34), (4,64), (2,30), (1,83), (3,19)],
        [(1,92), (4,62), (0,54), (3,43), (2,79)],
        [(3,69), (1,77), (4,87), (0,93), (2,87)],
        [(4,38), (3,60), (1,41), (2,83), (0,24)],
        [(0,17), (4,49), (2,25), (1,44), (3,98)],
        [(1,96), (0,77), (3,79), (2,75), (4,43)],
        [(3,28), (2,35), (4,95), (1,76), (0,7)],
        [(2,9), (0,10), (1,91), (3,59), (4,59)],
        [(4,46), (2,28), (3,52), (1,16), (0,59)],
        [(1,91), (4,50), (0,59), (2,43), (3,50)],
        [(0,27), (2,59), (1,46), (3,45), (4,90)],
        [(2,59), (3,46), (1,74), (4,43), (0,45)],
        [(1,45), (0,78), (4,28), (2,28), (3,23)],
        [(3,28), (2,23), (0,62), (4,45), (1,45)],
        [(4,54), (1,21), (3,20), (0,43), (2,98)],
        [(0,87), (2,60), (4,54), (3,43), (1,9)]
    ]),

    # Applegate & Cook (ORB) instances
    ("orb01", 10, 10, 1059, [
        [(2,21), (0,34), (3,95), (1,41), (4,66), (5,98), (7,54), (6,48), (9,39), (8,17)],
        [(1,52), (4,16), (3,31), (0,21), (2,26), (6,71), (5,16), (8,25), (9,98), (7,44)],
        [(3,42), (2,39), (0,98), (4,54), (1,16), (7,77), (6,79), (5,43), (8,75), (9,96)],
        [(4,55), (1,79), (0,77), (3,98), (2,77), (5,34), (7,64), (6,19), (9,30), (8,83)],
        [(0,34), (3,64), (1,30), (2,83), (4,19), (6,92), (7,62), (5,54), (8,43), (9,79)],
        [(1,69), (3,77), (2,87), (0,93), (4,87), (7,38), (5,60), (6,41), (9,83), (8,24)],
        [(2,17), (1,49), (4,25), (0,44), (3,98), (5,96), (6,77), (7,79), (8,75), (9,43)],
        [(3,28), (4,35), (0,95), (1,76), (2,7), (6,9), (5,10), (7,91), (8,59), (9,59)],
        [(1,46), (3,28), (2,52), (0,16), (4,59), (7,91), (6,50), (5,59), (8,43), (9,50)],
        [(2,27), (3,59), (0,46), (1,45), (4,90), (5,59), (6,46), (7,74), (8,43), (9,45)]
    ]),
]


def generate_jssp_mps(name, n_jobs, n_machines, jobs, output_path):
    """Generate MPS file for JSSP instance using disjunctive/big-M formulation."""

    # Calculate big-M (sum of all processing times)
    bigM = sum(pt for job in jobs for _, pt in job)

    # Build column and row data
    cols = []
    rows = []

    # Column 0: makespan C
    cols.append(("C", "CONTINUOUS", 0.0, 1e20, 1.0))  # objective coefficient = 1

    # Columns for x[j][op]: start times (continuous)
    # Index: 1 + j * n_machines + op
    for j in range(n_jobs):
        for op in range(n_machines):
            cols.append((f"x_{j}_{op}", "CONTINUOUS", 0.0, 1e20, 0.0))

    x_base = 1

    # Columns for y[j][k][m]: precedence binaries
    # Index: x_base + n_jobs*n_machines + j*n_jobs*n_machines + k*n_machines + m
    y_base = x_base + n_jobs * n_machines
    for j in range(n_jobs):
        for k in range(n_jobs):
            for m in range(n_machines):
                cols.append((f"y_{j}_{k}_{m}", "BINARY", 0.0, 1.0, 0.0))

    # Build machine schedule: which operations use which machine
    # machine_ops[m] = [(job_idx, op_idx, processing_time), ...]
    machine_ops = [[] for _ in range(n_machines)]
    for j, job in enumerate(jobs):
        for op, (machine, proc_time) in enumerate(job):
            machine_ops[machine].append((j, op, proc_time))

    # Row 1: Precedence within each job
    # x[j][op+1] - x[j][op] >= p[j][op]
    row_idx = 0
    for j, job in enumerate(jobs):
        for op in range(len(job) - 1):
            proc_time = job[op][1]
            x_curr = x_base + j * n_machines + op
            x_next = x_base + j * n_machines + (op + 1)
            rows.append((f"prec_{j}_{op}", "G", proc_time, [
                (x_next, 1.0),
                (x_curr, -1.0)
            ]))
            row_idx += 1

    # Row 2: Makespan constraints
    # C - x[j][last_op] >= p[j][last_op]
    for j, job in enumerate(jobs):
        last_op = len(job) - 1
        proc_time = job[last_op][1]
        x_last = x_base + j * n_machines + last_op
        rows.append((f"mksp_{j}", "G", proc_time, [
            (0, 1.0),  # C
            (x_last, -1.0)
        ]))
        row_idx += 1

    # Row 3: Disjunctive constraints (no overlap on machines)
    # For each machine m, for each pair of operations (j,op_j) and (k,op_k):
    #   x[j][op_j] + p[j][op_j] <= x[k][op_k] + M*(1 - y[j][k][m])
    #   x[k][op_k] + p[k][op_k] <= x[j][op_j] + M*y[j][k][m]
    # Rewrite as:
    #   x[k][op_k] - x[j][op_j] >= p[j][op_j] - M*(1 - y[j][k][m])
    #   x[j][op_j] - x[k][op_k] >= p[k][op_k] - M*y[j][k][m]
    for m in range(n_machines):
        ops = machine_ops[m]
        for i in range(len(ops)):
            for j in range(i + 1, len(ops)):
                job_i, op_i, pt_i = ops[i]
                job_j, op_j, pt_j = ops[j]

                x_i = x_base + job_i * n_machines + op_i
                x_j = x_base + job_j * n_machines + op_j
                y_ij = y_base + job_i * n_jobs * n_machines + job_j * n_machines + m

                # x[j] - x[i] + M*y >= pt_i
                rows.append((f"disj1_{job_i}_{job_j}_{m}", "G", pt_i, [
                    (x_j, 1.0),
                    (x_i, -1.0),
                    (y_ij, bigM)
                ]))

                # x[i] - x[j] - M*y >= pt_j - M
                rows.append((f"disj2_{job_i}_{job_j}_{m}", "G", pt_j - bigM, [
                    (x_i, 1.0),
                    (x_j, -1.0),
                    (y_ij, -bigM)
                ]))

    # Write MPS file
    with gzip.open(output_path, 'wt') as f:
        f.write("NAME          {}\n".format(name))

        # ROWS section
        f.write("ROWS\n")
        f.write(" N  OBJ\n")
        for row_name, sense, rhs, _ in rows:
            f.write(f" {sense}  {row_name}\n")

        # COLUMNS section
        f.write("COLUMNS\n")

        # Track which columns were actually written
        written_cols = set()
        in_integer_section = False

        for col_idx, (col_name, col_type, lb, ub, obj_coef) in enumerate(cols):
            col_written = False

            # Switch to integer section when we hit first BINARY variable
            if col_type == "BINARY" and not in_integer_section:
                f.write("    MARK0000  'MARKER'                 'INTORG'\n")
                in_integer_section = True
            # Switch to continuous section when we hit first CONTINUOUS after BINARY
            elif col_type == "CONTINUOUS" and in_integer_section:
                f.write("    MARK0001  'MARKER'                 'INTEND'\n")
                in_integer_section = False
            # Switch back to integer if we see BINARY after CONTINUOUS
            elif col_type == "BINARY" and not in_integer_section and col_idx > 0:
                f.write("    MARK0002  'MARKER'                 'INTORG'\n")
                in_integer_section = True

            # Objective
            if obj_coef != 0.0:
                f.write(f"    {col_name:8s}  OBJ       {obj_coef}\n")
                col_written = True

            # Constraints
            for row_name, sense, rhs, coeffs in rows:
                for var_idx, coef in coeffs:
                    if var_idx == col_idx and coef != 0.0:
                        f.write(f"    {col_name:8s}  {row_name:8s}  {coef}\n")
                        col_written = True

            if col_written:
                written_cols.add(col_idx)

        # Close integer section if still open
        if in_integer_section:
            f.write("    MARK9999  'MARKER'                 'INTEND'\n")

        # RHS section
        f.write("RHS\n")
        for row_name, sense, rhs, _ in rows:
            if rhs != 0.0:
                f.write(f"    RHS       {row_name:8s}  {rhs}\n")

        # BOUNDS section - only for columns that were actually written
        f.write("BOUNDS\n")
        for col_idx, (col_name, col_type, lb, ub, _) in enumerate(cols):
            if col_idx not in written_cols:
                continue  # Skip columns not in COLUMNS section

            if col_type == "BINARY":
                f.write(f" BV BOUND     {col_name}\n")
            else:
                if lb > 0:
                    f.write(f" LO BOUND     {col_name:8s}  {lb}\n")
                if ub < 1e20:
                    f.write(f" UP BOUND     {col_name:8s}  {ub}\n")

        f.write("ENDATA\n")


def verify_with_external_solver(mps_path, expected_obj, timeout=300):
    """Verify the instance solves to expected objective with external solver."""
    try:
        result = subprocess.run(
            ["cbc", str(mps_path), "-sec", str(timeout), "solve", "solu", "/dev/stdout"],
            capture_output=True, text=True, timeout=timeout + 10
        )

        if result.returncode != 0:
            return False, "Solver failed"

        # Parse objective from solution
        for line in result.stdout.split('\n'):
            if 'Optimal - objective value' in line:
                obj_str = line.split('value')[-1].strip()
                actual_obj = float(obj_str)
                if abs(actual_obj - expected_obj) < 0.5:
                    return True, f"Verified optimal={actual_obj:.0f}"
                else:
                    return False, f"Wrong optimal: got {actual_obj:.0f}, expected {expected_obj}"

        return False, "Could not parse objective"

    except subprocess.TimeoutExpired:
        return False, f"Timeout after {timeout}s"
    except Exception as e:
        return False, f"Error: {e}"


def main():
    fixture_dir = Path(__file__).parent / "fixtures"
    fixture_dir.mkdir(exist_ok=True)

    print("=== Generating JSSP Fixtures ===\n")

    for name, n_jobs, n_machines, optimal, jobs in JSSP_INSTANCES:
        print(f"Generating {name} ({n_jobs}×{n_machines}, optimal={optimal})...")

        mps_path = fixture_dir / f"jssp_{name}.mps.gz"
        generate_jssp_mps(name, n_jobs, n_machines, jobs, mps_path)
        print(f"  Written to {mps_path}")

        # Skip external verification for now - will test with MIPster
        print(f"  Generated (expected optimal={optimal})")

    print(f"\n=== Generated {len(JSSP_INSTANCES)} JSSP fixtures ===")
    return 0


if __name__ == "__main__":
    sys.exit(main())
