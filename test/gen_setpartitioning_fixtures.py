#!/usr/bin/env python3
"""
Generate Set Partitioning Problem (SPP) fixtures.

Problem:
  minimize   sum_j c[j] * x[j]
  subject to sum_{j: i in S[j]} x[j] = 1   for all i (element i in exactly one set)
             x[j] in {0,1}

where S[j] is the set of elements in set j.

Set partitioning is stricter than set covering (= vs >=) and set packing (<= vs =).
This structure appears in:
  - Crew scheduling (each flight covered by exactly one crew)
  - Vehicle routing (each customer visited exactly once)
  - Graph partitioning

Key properties:
  - Pure binary problem
  - Clique cuts effective on conflict structure
  - Often highly degenerate (many optimal solutions)
  - Preprocessing can identify implied equalities

Instances are generated with various patterns:
  - uniform: random partitions
  - overlapping: sets share elements (needs auxiliary sets)
  - structured: special patterns (e.g., graph partitioning)

Solved with external solver to obtain certified optimal values.
"""

import sys
import math
import random
import gzip
from pathlib import Path


def generate_spp_instance(n_elements, n_sets, seed, pattern="uniform"):
    """Generate a set partitioning instance.

    Args:
        n_elements: number of elements (universe size)
        n_sets: number of candidate sets
        seed: random seed
        pattern: generation pattern ("uniform", "overlapping", "structured")
    """
    rng = random.Random(seed)

    # Set costs
    costs = [rng.randint(1, 100) for _ in range(n_sets)]

    # Set memberships: members[j] = list of elements in set j
    members = [[] for _ in range(n_sets)]

    if pattern == "uniform":
        # Generate random partitions; ensure each element appears in multiple sets
        # so the problem is non-trivial
        for i in range(n_elements):
            # Each element appears in 2-4 sets
            n_appearances = rng.randint(2, min(4, n_sets))
            sets_for_i = rng.sample(range(n_sets), n_appearances)
            for j in sets_for_i:
                members[j].append(i)

    elif pattern == "overlapping":
        # Start with a valid partition, then add alternative covers
        partition = list(range(n_elements))
        rng.shuffle(partition)

        # First n_sets/2 sets form a partition
        sets_per_part = max(1, n_sets // 2)
        elem_per_set = max(1, n_elements // sets_per_part)

        for j in range(sets_per_part):
            start = j * elem_per_set
            end = min((j+1) * elem_per_set, n_elements)
            members[j] = partition[start:end]

        # Remaining sets are random covers
        for j in range(sets_per_part, n_sets):
            size = rng.randint(2, max(3, n_elements // 3))
            members[j] = rng.sample(range(n_elements), size)

    elif pattern == "structured":
        # Graph partitioning pattern: partition nodes of a grid graph
        # Each set is a contiguous region
        side = int(math.sqrt(n_elements)) + 1
        grid = list(range(n_elements))

        for j in range(n_sets):
            # Random region: start point + grow in 4 directions
            if j < n_sets // 2:
                # Small contiguous regions
                size = rng.randint(2, min(5, n_elements))
                start = rng.randint(0, n_elements - 1)
                members[j] = [start]
                candidates = [start]
                while len(members[j]) < size and candidates:
                    current = candidates.pop(0)
                    # Add neighbors (4-connected grid)
                    row, col = current // side, current % side
                    for dr, dc in [(-1,0), (1,0), (0,-1), (0,1)]:
                        nr, nc = row + dr, col + dc
                        if 0 <= nr < side and 0 <= nc < side:
                            neighbor = nr * side + nc
                            if neighbor < n_elements and neighbor not in members[j]:
                                members[j].append(neighbor)
                                candidates.append(neighbor)
            else:
                # Random covers
                size = rng.randint(2, max(3, n_elements // 4))
                members[j] = rng.sample(range(n_elements), size)

    # Ensure every element appears in at least one set
    for i in range(n_elements):
        appears = any(i in members[j] for j in range(n_sets))
        if not appears:
            j = rng.randint(0, n_sets - 1)
            members[j].append(i)

    return costs, members


def write_spp_mps(filename, n_elements, n_sets, costs, members):
    """Write set partitioning instance in MPS format."""
    with open(filename, 'w') as fp:
        fp.write("NAME          SETPART\n")

        # ROWS section
        fp.write("ROWS\n")
        fp.write(" N  OBJ\n")
        for i in range(n_elements):
            fp.write(f" E  PART_{i}\n")  # Equality constraint

        # COLUMNS section
        fp.write("COLUMNS\n")
        for j in range(n_sets):
            var_name = f"x_{j}"
            # Objective coefficient
            fp.write(f"    {var_name:<10s}  OBJ       {costs[j]}\n")
            # Partitioning constraints
            for i in members[j]:
                fp.write(f"    {var_name:<10s}  PART_{i}  1.0\n")

        # RHS section
        fp.write("RHS\n")
        for i in range(n_elements):
            fp.write(f"    RHS1      PART_{i}  1.0\n")

        # BOUNDS section
        fp.write("BOUNDS\n")
        for j in range(n_sets):
            fp.write(f" BV BND1      x_{j}\n")

        fp.write("ENDATA\n")


def solve_spp_external(n_elements, n_sets, costs, members, timeout=60):
    """Solve set partitioning with external solver to get certified optimal value."""
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

    # Constraints: each element in exactly one set
    for i in range(n_elements):
        # Find sets containing element i
        sets_with_i = [j for j in range(n_sets) if i in members[j]]
        if sets_with_i:
            m.addConstr(gp.quicksum(x[j] for j in sets_with_i) == 1,
                       name=f"part_{i}")
        else:
            # Infeasible if element not covered
            print(f"WARNING: element {i} not covered by any set", file=sys.stderr)
            env.dispose()
            return None, "infeasible"

    m.optimize()

    if m.Status == GRB.OPTIMAL:
        obj_val = m.ObjVal
        env.dispose()
        return obj_val, "optimal"
    elif m.Status == GRB.TIME_LIMIT:
        obj_val = m.ObjVal if m.SolCount > 0 else None
        env.dispose()
        return obj_val, "timelimit"
    elif m.Status == GRB.INFEASIBLE:
        env.dispose()
        return None, "infeasible"
    else:
        env.dispose()
        return None, f"status_{m.Status}"


# Test specifications: (n_elements, n_sets, seed, pattern, timeout)
SPECS = [
    # Small instances for quick CI
    (10,  20,  42,  "uniform",     10),
    (15,  30,  42,  "overlapping", 20),
    (20,  40,  137, "structured",  30),

    # Medium instances
    (30,  60,  42,  "uniform",     60),
    (40,  80,  137, "overlapping", 90),
    (50, 100,  42,  "structured",  120),
]


def main():
    fixture_dir = Path(__file__).parent / "fixtures"
    fixture_dir.mkdir(exist_ok=True)

    results = []

    for n_elem, n_sets, seed, pattern, timeout in SPECS:
        name = f"sprt_e{n_elem}_s{n_sets}_sd{seed}_{pattern}"
        print(f"Generating {name}...", flush=True)

        costs, members = generate_spp_instance(n_elem, n_sets, seed, pattern)

        # Write MPS
        mps_path = f"/tmp/{name}.mps"
        write_spp_mps(mps_path, n_elem, n_sets, costs, members)

        # Solve with external solver
        obj, status = solve_spp_external(n_elem, n_sets, costs, members, timeout)

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
    print("Add to test/CInterfaceTest_setpartitioning.c:")
    print("="*70)
    print("typedef struct {")
    print("  const char *name;")
    print("  double expected_obj;")
    print("  int timeout_sec;")
    print("} SetPartitioningTestCase;\n")
    print("static const SetPartitioningTestCase setpartitioning_test_cases[] = {")
    for name, n_elem, n_sets, obj, timeout in results:
        print(f'  {{"{name}", {obj:.2f}, {timeout}}},')
    print("};")


if __name__ == "__main__":
    main()
