#!/usr/bin/env python3
"""
Generate Set Covering Problem (SCP) fixtures.

Problem:
  minimize   sum_j c[j] * x[j]
  subject to sum_{j: i in S[j]} x[j] >= 1   for all i (element i must be covered)
             x[j] in {0,1}

where S[j] is the set of elements covered by set j.

This is a classic pure binary problem where:
  - Clique cuts should be very effective (conflicts between sets)
  - Odd-wheel cuts can strengthen formulation
  - Preprocessing can identify dominated sets

Instances are generated with various patterns:
  - uniform: random coverage
  - sparse: low overlap between sets
  - dense: high overlap
  - hierarchical: some sets subsume others

Solved with external solver to obtain certified optimal values.
"""

import sys
import math
import random
import gzip
from pathlib import Path


def generate_scp_instance(n_elements, n_sets, seed, pattern="uniform", density=0.3, cost_type="integer"):
    """Generate a set covering instance.

    Args:
        n_elements: number of elements to cover
        n_sets: number of available sets
        seed: random seed
        pattern: coverage pattern ("uniform", "sparse", "dense", "hierarchical")
        density: coverage density (prob that set j covers element i)
        cost_type: "integer", "continuous", "large", "fractional"
    """
    rng = random.Random(seed)

    # Set costs with diversity
    if cost_type == "integer":
        costs = [rng.randint(1, 100) for _ in range(n_sets)]
    elif cost_type == "continuous":
        costs = [round(rng.uniform(1.0, 100.0), 2) for _ in range(n_sets)]
    elif cost_type == "large":
        costs = [rng.randint(100, 1000) for _ in range(n_sets)]
    elif cost_type == "fractional":
        costs = [round(rng.uniform(0.1, 10.0), 2) for _ in range(n_sets)]
    else:
        costs = [rng.randint(1, 100) for _ in range(n_sets)]

    # Coverage matrix: covers[j][i] = 1 if set j covers element i
    covers = [[0] * n_elements for _ in range(n_sets)]

    if pattern == "uniform":
        for j in range(n_sets):
            for i in range(n_elements):
                if rng.random() < density:
                    covers[j][i] = 1

    elif pattern == "sparse":
        # Each set covers a small random subset (20% of elements)
        for j in range(n_sets):
            k = max(1, int(0.2 * n_elements))
            elements = rng.sample(range(n_elements), k)
            for i in elements:
                covers[j][i] = 1

    elif pattern == "dense":
        # Each set covers a large random subset (70% of elements)
        for j in range(n_sets):
            k = max(1, int(0.7 * n_elements))
            elements = rng.sample(range(n_elements), k)
            for i in elements:
                covers[j][i] = 1

    elif pattern == "hierarchical":
        # Some sets subsume others; test dominated column detection
        for j in range(n_sets):
            # Base coverage
            for i in range(n_elements):
                if rng.random() < density:
                    covers[j][i] = 1
        # Make every 5th set a superset of previous one
        for j in range(5, n_sets, 5):
            if j > 0:
                for i in range(n_elements):
                    if covers[j-1][i]:
                        covers[j][i] = 1

    # Ensure every element is covered by at least one set
    for i in range(n_elements):
        if sum(covers[j][i] for j in range(n_sets)) == 0:
            j = rng.randint(0, n_sets - 1)
            covers[j][i] = 1

    return costs, covers


def write_scp_mps(filename, n_elements, n_sets, costs, covers):
    """Write SCP instance in MPS format."""
    with open(filename, 'w') as fp:
        fp.write("NAME          SETCOVER\n")

        # ROWS section
        fp.write("ROWS\n")
        fp.write(" N  OBJ\n")
        for i in range(n_elements):
            fp.write(f" G  COVER_{i}\n")

        # COLUMNS section
        fp.write("COLUMNS\n")
        for j in range(n_sets):
            var_name = f"x_{j}"
            # Objective coefficient
            fp.write(f"    {var_name:<10s}  OBJ       {costs[j]}\n")
            # Coverage constraints
            for i in range(n_elements):
                if covers[j][i]:
                    fp.write(f"    {var_name:<10s}  COVER_{i}  1.0\n")

        # RHS section
        fp.write("RHS\n")
        for i in range(n_elements):
            fp.write(f"    RHS1      COVER_{i}  1.0\n")

        # BOUNDS section
        fp.write("BOUNDS\n")
        for j in range(n_sets):
            fp.write(f" BV BND1      x_{j}\n")

        fp.write("ENDATA\n")


def solve_scp_external(n_elements, n_sets, costs, covers, timeout=60):
    """Solve SCP with external solver to get certified optimal value."""
    try:
        import gurobipy as gp
        from gurobipy import GRB
    except ImportError:
        print("WARNING: external solver not available, skipping solve", file=sys.stderr)
        return None, None

    env = gp.Env(empty=True)
    env.setParam("OutputFlag", 0)
    env.setParam("TimeLimit", timeout)
    env.start()

    m = gp.Model(env=env)

    # Variables: x[j] = 1 if set j is selected
    x = [m.addVar(vtype=GRB.BINARY, obj=costs[j], name=f"x_{j}")
         for j in range(n_sets)]

    # Constraints: each element covered at least once
    for i in range(n_elements):
        m.addConstr(gp.quicksum(covers[j][i] * x[j] for j in range(n_sets)) >= 1,
                   name=f"cover_{i}")

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


# Test specifications: (n_elements, n_sets, seed, pattern, density, cost_type, timeout)
SPECS = [
    # Small instances for quick CI - diverse cost types
    (10,  15,  42,  "uniform",      0.3, "integer",     10),
    (15,  20,  42,  "sparse",       0.3, "continuous",  15),
    (20,  30,  42,  "dense",        0.3, "integer",     20),
    (20,  25,  137, "hierarchical", 0.3, "fractional",  20),

    # Medium instances - diverse cost types
    (30,  50,  42,  "uniform",      0.3, "large",       60),
    (40,  60,  137, "sparse",       0.3, "integer",     60),
    (50,  80,  42,  "dense",        0.3, "continuous",  90),
]


def main():
    fixture_dir = Path(__file__).parent / "fixtures"
    fixture_dir.mkdir(exist_ok=True)

    results = []

    for n_elem, n_sets, seed, pattern, density, cost_type, timeout in SPECS:
        name = f"scp_e{n_elem}_s{n_sets}_sd{seed}_{pattern}"
        print(f"Generating {name} (costs={cost_type})...", flush=True)

        costs, covers = generate_scp_instance(n_elem, n_sets, seed, pattern, density, cost_type)

        # Write MPS
        mps_path = f"/tmp/{name}.mps"
        write_scp_mps(mps_path, n_elem, n_sets, costs, covers)

        # Solve with external solver
        obj, status = solve_scp_external(n_elem, n_sets, costs, covers, timeout)

        if status == "optimal":
            print(f"  {name}: obj={obj:.2f} (optimal)")
            results.append((name, n_elem, n_sets, obj, timeout))

            # Compress to fixtures/
            gz_path = fixture_dir / f"{name}.mps.gz"
            with open(mps_path, 'rb') as f_in:
                with gzip.open(gz_path, 'wb') as f_out:
                    f_out.writelines(f_in)
            print(f"  Written to {gz_path}")
        else:
            print(f"  {name}: {status} (skipping)")

    # Generate C test file snippet
    print("\n" + "="*70)
    print("Add to test/CInterfaceTest_setcover.c:")
    print("="*70)
    print("typedef struct {")
    print("  const char *name;")
    print("  double expected_obj;")
    print("  int timeout_sec;")
    print("} SetCoverTestCase;\n")
    print("static const SetCoverTestCase setcover_test_cases[] = {")
    for name, n_elem, n_sets, obj, timeout in results:
        print(f'  {{"{name}", {obj:.2f}, {timeout}}},')
    print("};")


if __name__ == "__main__":
    main()
