#!/usr/bin/env python3
"""
Generate Quadratic Assignment Problem (QAP) fixtures.

Problem:
  minimize   sum_i sum_j f[i][j] * d[x[i]][x[j]]
  subject to sum_j x[i][j] = 1   for all i (each facility assigned once)
             sum_i x[i][j] = 1   for all j (each location used once)
             x[i][j] in {0,1}

Linearization (Kaufman-Broeckx):
  minimize   sum_i sum_j sum_k sum_l f[i][j] * d[k][l] * y[i][j][k][l]
  subject to sum_j x[i][j] = 1   for all i
             sum_i x[i][j] = 1   for all j
             y[i][j][k][l] >= x[i][k] + x[j][l] - 1   for all i,j,k,l with i<j
             x[i][j], y[i][j][k][l] in {0,1}

Instances are solved with external MIP solver to obtain certified optimal values.
Outputs MPS files for MIPster integration tests.
"""

import sys
import math
import random
from pathlib import Path


def generate_qap_instance(n, seed, flow_type="uniform", dist_type="euclidean"):
    """Generate a QAP instance with n facilities/locations.

    Args:
        n: problem size (number of facilities and locations)
        seed: random seed
        flow_type: "uniform", "sparse", "hub" (flow matrix pattern)
        dist_type: "euclidean", "manhattan", "grid" (distance matrix pattern)
    """
    rng = random.Random(seed)

    # Flow matrix f[i][j] (facility i to facility j)
    f = [[0] * n for _ in range(n)]
    if flow_type == "uniform":
        for i in range(n):
            for j in range(n):
                if i != j:
                    f[i][j] = rng.randint(1, 50)
    elif flow_type == "sparse":
        # Only 30% of pairs have flow
        for i in range(n):
            for j in range(n):
                if i != j and rng.random() < 0.3:
                    f[i][j] = rng.randint(5, 100)
    elif flow_type == "hub":
        # Hub-and-spoke: facility 0 is hub with high flow
        for i in range(1, n):
            f[0][i] = f[i][0] = rng.randint(50, 100)
        for i in range(1, n):
            for j in range(i+1, n):
                f[i][j] = f[j][i] = rng.randint(1, 10)

    # Distance matrix d[k][l] (location k to location l)
    d = [[0] * n for _ in range(n)]
    if dist_type == "euclidean":
        # Random 2D coordinates
        coords = [(rng.randint(0, 100), rng.randint(0, 100)) for _ in range(n)]
        for k in range(n):
            for l in range(n):
                if k != l:
                    dx = coords[k][0] - coords[l][0]
                    dy = coords[k][1] - coords[l][1]
                    d[k][l] = int(math.sqrt(dx*dx + dy*dy))
    elif dist_type == "manhattan":
        coords = [(rng.randint(0, 100), rng.randint(0, 100)) for _ in range(n)]
        for k in range(n):
            for l in range(n):
                if k != l:
                    d[k][l] = abs(coords[k][0] - coords[l][0]) + abs(coords[k][1] - coords[l][1])
    elif dist_type == "grid":
        # Grid layout: n locations arranged in sqrt(n) x sqrt(n) grid
        side = int(math.sqrt(n))
        for k in range(n):
            for l in range(n):
                if k != l:
                    row_k, col_k = k // side, k % side
                    row_l, col_l = l // side, l % side
                    d[k][l] = abs(row_k - row_l) + abs(col_k - col_l)

    return f, d


def write_qap_mps(filename, n, f, d):
    """Write QAP instance in MPS format (linearized formulation).

    We use a compact linearization with auxiliary variables y[i][k][j][l] only for i < j.
    Variables:
      x[i][k]: facility i assigned to location k
      y[i][k][j][l]: x[i][k] AND x[j][l] for i < j
    """
    with open(filename, 'w') as fp:
        fp.write("NAME          QAP\n")

        # ROWS section
        fp.write("ROWS\n")
        fp.write(" N  OBJ\n")

        # Assignment constraints: each facility to one location
        for i in range(n):
            fp.write(f" E  FACIL{i}\n")

        # Assignment constraints: each location used once
        for k in range(n):
            fp.write(f" E  LOC{k}\n")

        # McCormick linearization: y[i][k][j][l] >= x[i][k] + x[j][l] - 1
        for i in range(n):
            for k in range(n):
                for j in range(i+1, n):
                    for l in range(n):
                        if f[i][j] * d[k][l] > 0:  # Only if non-zero cost
                            fp.write(f" G  MC_{i}_{k}_{j}_{l}\n")

        # COLUMNS section
        fp.write("COLUMNS\n")

        # x[i][k] variables
        for i in range(n):
            for k in range(n):
                var_name = f"x_{i}_{k}"
                # Assignment constraints
                fp.write(f"    {var_name:<10s}  FACIL{i}  1.0\n")
                fp.write(f"    {var_name:<10s}  LOC{k}  1.0\n")
                # McCormick constraints (appear in y terms)

        # y[i][k][j][l] variables (auxiliary for linearization)
        for i in range(n):
            for k in range(n):
                for j in range(i+1, n):
                    for l in range(n):
                        cost = f[i][j] * d[k][l]
                        if cost > 0:
                            var_name = f"y_{i}_{k}_{j}_{l}"
                            # Objective coefficient
                            fp.write(f"    {var_name:<10s}  OBJ       {cost}\n")
                            # McCormick constraint: y >= x[i][k] + x[j][l] - 1
                            fp.write(f"    {var_name:<10s}  MC_{i}_{k}_{j}_{l}  1.0\n")

        # Add x terms to McCormick constraints
        for i in range(n):
            for k in range(n):
                for j in range(i+1, n):
                    for l in range(n):
                        if f[i][j] * d[k][l] > 0:
                            x_ik = f"x_{i}_{k}"
                            x_jl = f"x_{j}_{l}"
                            fp.write(f"    {x_ik:<10s}  MC_{i}_{k}_{j}_{l}  -1.0\n")
                            fp.write(f"    {x_jl:<10s}  MC_{i}_{k}_{j}_{l}  -1.0\n")

        # RHS section
        fp.write("RHS\n")
        for i in range(n):
            fp.write(f"    RHS1      FACIL{i}  1.0\n")
        for k in range(n):
            fp.write(f"    RHS1      LOC{k}  1.0\n")
        for i in range(n):
            for k in range(n):
                for j in range(i+1, n):
                    for l in range(n):
                        if f[i][j] * d[k][l] > 0:
                            fp.write(f"    RHS1      MC_{i}_{k}_{j}_{l}  -1.0\n")

        # BOUNDS section
        fp.write("BOUNDS\n")
        # All variables binary by default
        for i in range(n):
            for k in range(n):
                fp.write(f" BV BND1      x_{i}_{k}\n")
        for i in range(n):
            for k in range(n):
                for j in range(i+1, n):
                    for l in range(n):
                        if f[i][j] * d[k][l] > 0:
                            fp.write(f" BV BND1      y_{i}_{k}_{j}_{l}\n")

        fp.write("ENDATA\n")


def solve_qap_external(n, f, d, timeout=60):
    """Solve QAP with external solver to get certified optimal value."""
    try:
        import gurobipy as gp
        from gurobipy import GRB
    except ImportError:
        print("WARNING: External solver not available, skipping solve", file=sys.stderr)
        return None, None

    env = gp.Env(empty=True)
    env.setParam("OutputFlag", 0)
    env.setParam("TimeLimit", timeout)
    env.start()

    m = gp.Model(env=env)

    # Variables: x[i][k] = 1 if facility i assigned to location k
    x = {}
    for i in range(n):
        for k in range(n):
            x[i,k] = m.addVar(vtype=GRB.BINARY, name=f"x_{i}_{k}")

    # Auxiliary variables for linearization
    y = {}
    for i in range(n):
        for k in range(n):
            for j in range(i+1, n):
                for l in range(n):
                    if f[i][j] * d[k][l] > 0:
                        y[i,k,j,l] = m.addVar(vtype=GRB.BINARY, name=f"y_{i}_{k}_{j}_{l}")

    # Objective: minimize flow * distance
    obj = gp.QuadExpr()
    for i in range(n):
        for k in range(n):
            for j in range(i+1, n):
                for l in range(n):
                    if f[i][j] * d[k][l] > 0:
                        obj += f[i][j] * d[k][l] * y[i,k,j,l]
    m.setObjective(obj, GRB.MINIMIZE)

    # Constraints: each facility assigned once
    for i in range(n):
        m.addConstr(gp.quicksum(x[i,k] for k in range(n)) == 1, name=f"facil_{i}")

    # Constraints: each location used once
    for k in range(n):
        m.addConstr(gp.quicksum(x[i,k] for i in range(n)) == 1, name=f"loc_{k}")

    # McCormick linearization: y[i,k,j,l] >= x[i,k] + x[j,l] - 1
    for i in range(n):
        for k in range(n):
            for j in range(i+1, n):
                for l in range(n):
                    if f[i][j] * d[k][l] > 0:
                        m.addConstr(y[i,k,j,l] >= x[i,k] + x[j,l] - 1)

    m.optimize()

    if m.Status == GRB.OPTIMAL:
        obj_val = m.ObjVal
        env.dispose()
        return obj_val, "optimal"
    elif m.Status == GRB.TIME_LIMIT:
        obj_val = m.ObjVal if m.SolCount > 0 else None
        env.dispose()
        return obj_val, "timelimit"
    else:
        env.dispose()
        return None, f"status_{m.Status}"


# Test specifications: (n, seed, flow_type, dist_type, timeout)
SPECS = [
    # Small instances for quick CI tests
    (4,  42,  "uniform", "euclidean", 10),
    (5,  42,  "uniform", "euclidean", 30),
    (5,  137, "sparse",  "manhattan", 30),
    (6,  42,  "uniform", "grid",      60),
    (6,  137, "hub",     "euclidean", 60),

    # Medium instances
    (7,  42,  "uniform", "euclidean", 120),
    (7,  137, "sparse",  "grid",      120),
]


def main():
    import gzip

    fixture_dir = Path(__file__).parent / "fixtures"
    fixture_dir.mkdir(exist_ok=True)

    results = []

    for n, seed, flow_type, dist_type, timeout in SPECS:
        name = f"qap_n{n}_s{seed}_{flow_type}_{dist_type}"
        print(f"Generating {name}...", flush=True)

        f, d = generate_qap_instance(n, seed, flow_type, dist_type)

        # Write MPS
        mps_path = f"/tmp/{name}.mps"
        write_qap_mps(mps_path, n, f, d)

        # Solve with external solver
        obj, status = solve_qap_external(n, f, d, timeout)

        if status == "optimal":
            print(f"  {name}: obj={obj:.2f} (optimal)")
            results.append((name, n, obj, timeout))

            # Compress to fixtures/
            gz_path = fixture_dir / f"{name}.mps.gz"
            with open(mps_path, 'rb') as f_in:
                with gzip.open(gz_path, 'wb') as f_out:
                    f_out.writelines(f_in)
            print(f"  Written to {gz_path}")
        else:
            print(f"  {name}: {status} (skipping)")

    # Generate C test file header
    print("\n" + "="*70)
    print("Add to test/CInterfaceTest_qap.c:")
    print("="*70)
    print("typedef struct {")
    print("  const char *name;")
    print("  double expected_obj;")
    print("  int timeout_sec;")
    print("} QapTestCase;\n")
    print("static const QapTestCase qap_test_cases[] = {")
    for name, n, obj, timeout in results:
        print(f'  {{"{name}", {obj:.2f}, {timeout}}},')
    print("};")


if __name__ == "__main__":
    main()
