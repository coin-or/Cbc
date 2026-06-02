#!/usr/bin/env python3
"""
Generate Set Packing Problem (SPP) fixtures.

Problem:
  maximize   sum_j c[j] * x[j]
  subject to sum_{j: i in S[j]} x[j] <= 1   for all i (element i in at most one set)
             x[j] in {0,1}

where S[j] is the set of elements in set j.

This is a classic pure binary problem where:
  - Clique cuts are ESSENTIAL (conflicts form maximal cliques)
  - Odd-wheel cuts strengthen LP relaxation
  - Conflict graph is the core structure

Set packing is equivalent to:
  - Maximum weighted independent set in the conflict graph
  - Maximum weighted clique in the complement graph

Instances are generated with various conflict patterns:
  - random: random element overlaps
  - k_clique: k-partite structure
  - geometric: geometric intersection patterns
  - structured: special conflict patterns (e.g., scheduling)

Solved with external solver to obtain certified optimal values.
"""

import sys
import math
import random
import gzip
from pathlib import Path


def generate_spp_instance(n_elements, n_sets, seed, pattern="random", k=3, value_type="integer"):
    """Generate a set packing instance.

    Args:
        n_elements: number of elements (universe size)
        n_sets: number of candidate sets
        seed: random seed
        pattern: conflict pattern ("random", "k_partite", "geometric", "scheduling")
        k: parameter for pattern (e.g., clique size for k_partite)
        value_type: "integer", "continuous", "large", "fractional"
    """
    rng = random.Random(seed)

    # Set values (weights) with diversity
    if value_type == "integer":
        values = [rng.randint(1, 100) for _ in range(n_sets)]
    elif value_type == "continuous":
        values = [round(rng.uniform(1.0, 100.0), 2) for _ in range(n_sets)]
    elif value_type == "large":
        values = [rng.randint(100, 1000) for _ in range(n_sets)]
    elif value_type == "fractional":
        values = [round(rng.uniform(0.1, 10.0), 2) for _ in range(n_sets)]
    else:
        values = [rng.randint(1, 100) for _ in range(n_sets)]

    # Set memberships: members[j] = list of elements in set j
    members = [[] for _ in range(n_sets)]

    if pattern == "random":
        # Each set contains k random elements
        for j in range(n_sets):
            size = rng.randint(2, min(k, n_elements))
            members[j] = rng.sample(range(n_elements), size)

    elif pattern == "k_partite":
        # k-partite structure: sets form k cliques
        # Sets within same clique share elements (conflict)
        sets_per_clique = n_sets // k
        for clique_id in range(k):
            # Common elements for this clique
            common = rng.sample(range(n_elements), min(3, n_elements))
            for i in range(sets_per_clique):
                j = clique_id * sets_per_clique + i
                if j < n_sets:
                    members[j] = common + rng.sample(
                        [e for e in range(n_elements) if e not in common],
                        min(2, n_elements - len(common))
                    )

    elif pattern == "geometric":
        # Geometric pattern: intervals on a line
        # Set j = [start[j], end[j])
        # Conflicts when intervals overlap
        for j in range(n_sets):
            duration = rng.randint(2, min(5, n_elements))
            start = rng.randint(0, n_elements - 1)
            end = min(start + duration, n_elements)
            members[j] = list(range(start, end))
            if not members[j]:  # Empty set
                members[j] = [rng.randint(0, n_elements-1)]

    elif pattern == "scheduling":
        # Scheduling pattern: jobs with time windows
        # Each set represents a job covering time slots
        # Jobs conflict if time windows overlap
        for j in range(n_sets):
            duration = rng.randint(2, min(5, n_elements))
            start = rng.randint(0, max(1, n_elements - duration))
            members[j] = list(range(start, min(start + duration, n_elements)))

    # Ensure non-empty sets
    for j in range(n_sets):
        if not members[j]:
            members[j] = [rng.randint(0, n_elements - 1)]

    return values, members


def write_spp_mps(filename, n_elements, n_sets, values, members):
    """Write SPP instance in MPS format."""
    with open(filename, 'w') as fp:
        fp.write("NAME          SETPACK\n")

        # ROWS section
        fp.write("ROWS\n")
        fp.write(" N  OBJ\n")
        for i in range(n_elements):
            fp.write(f" L  PACK_{i}\n")

        # COLUMNS section
        fp.write("COLUMNS\n")
        for j in range(n_sets):
            var_name = f"x_{j}"
            # Objective coefficient (negate for maximization)
            fp.write(f"    {var_name:<10s}  OBJ       {-values[j]}\n")
            # Packing constraints
            for i in members[j]:
                fp.write(f"    {var_name:<10s}  PACK_{i}  1.0\n")

        # RHS section
        fp.write("RHS\n")
        for i in range(n_elements):
            fp.write(f"    RHS1      PACK_{i}  1.0\n")

        # BOUNDS section
        fp.write("BOUNDS\n")
        for j in range(n_sets):
            fp.write(f" BV BND1      x_{j}\n")

        fp.write("ENDATA\n")


def solve_spp_external(n_elements, n_sets, values, members, timeout=60):
    """Solve SPP with external solver to get certified optimal value."""
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
    x = [m.addVar(vtype=GRB.BINARY, obj=values[j], name=f"x_{j}")
         for j in range(n_sets)]

    m.setObjective(gp.quicksum(values[j] * x[j] for j in range(n_sets)), GRB.MAXIMIZE)

    # Constraints: each element in at most one set
    for i in range(n_elements):
        # Find sets containing element i
        sets_with_i = [j for j in range(n_sets) if i in members[j]]
        if sets_with_i:
            m.addConstr(gp.quicksum(x[j] for j in sets_with_i) <= 1,
                       name=f"pack_{i}")

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


# Test specifications: (n_elements, n_sets, seed, pattern, k, value_type, timeout)
SPECS = [
    # Small instances for quick CI - diverse value types
    (10,  15,  42,  "random",     3, "integer",     10),
    (12,  18,  42,  "k_partite",  3, "continuous",  15),
    (15,  20,  137, "geometric",  3, "integer",     20),
    (20,  25,  42,  "scheduling", 4, "fractional",  20),

    # Medium instances - diverse value types
    (30,  50,  42,  "random",     4, "large",       60),
    (35,  60,  137, "k_partite",  4, "integer",     60),
    (40,  70,  42,  "geometric",  5, "continuous",  90),
    (50,  80,  137, "scheduling", 5, "integer",     90),
]


def main():
    fixture_dir = Path(__file__).parent / "fixtures"
    fixture_dir.mkdir(exist_ok=True)

    results = []

    for n_elem, n_sets, seed, pattern, k, value_type, timeout in SPECS:
        name = f"spp_e{n_elem}_s{n_sets}_sd{seed}_{pattern}"
        print(f"Generating {name} (values={value_type})...", flush=True)

        values, members = generate_spp_instance(n_elem, n_sets, seed, pattern, k, value_type)

        # Write MPS
        mps_path = f"/tmp/{name}.mps"
        write_spp_mps(mps_path, n_elem, n_sets, values, members)

        # Solve with external solver
        obj, status = solve_spp_external(n_elem, n_sets, values, members, timeout)

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
    print("Add to test/CInterfaceTest_setpacking.c:")
    print("="*70)
    print("typedef struct {")
    print("  const char *name;")
    print("  double expected_obj;")
    print("  int timeout_sec;")
    print("} SetPackingTestCase;\n")
    print("static const SetPackingTestCase setpacking_test_cases[] = {")
    for name, n_elem, n_sets, obj, timeout in results:
        print(f'  {{"{name}", {obj:.2f}, {timeout}}},')
    print("};")


if __name__ == "__main__":
    main()
