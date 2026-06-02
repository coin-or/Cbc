#!/usr/bin/env python3
"""
Generate Bin Packing with Conflicts (BPC) fixtures.

Problem:
  minimize   sum_b y[b]
  subject to sum_b x[i][b] = 1              for all i (item i assigned to one bin)
             sum_i w[i] * x[i][b] <= C * y[b]  for all b (capacity constraint)
             x[i][b] + x[j][b] <= y[b]      for all conflicts (i,j) in E
             x[i][b], y[b] in {0,1}

where:
  - w[i] is the weight of item i
  - C is the bin capacity
  - E is the set of conflicting item pairs (cannot be in same bin)
  - y[b] = 1 if bin b is used

This combines classical bin packing with conflict constraints, appearing in:
  - Task scheduling with resource conflicts
  - Storage assignment with incompatibility
  - Warehouse packing with hazardous materials
  - Graph coloring with weighted vertices

Key properties:
  - Pure binary problem
  - Clique cuts ESSENTIAL for conflict structure
  - Knapsack cuts effective for capacity constraints
  - Combination of packing + covering constraints

Conflict patterns:
  - random: random item pairs conflict
  - clique: items form k-cliques (mutual exclusion sets)
  - bipartite: two types of items, conflicts across types
  - geometric: items conflict if "close" in some space

Solved with external solver to obtain certified optimal values.
"""

import sys
import math
import random
import gzip
from pathlib import Path


def generate_bpc_instance(n_items, capacity, seed, conflict_pattern="random",
                          conflict_density=0.2, weight_pattern="uniform"):
    """Generate a bin packing with conflicts instance.

    Args:
        n_items: number of items to pack
        capacity: bin capacity
        seed: random seed
        conflict_pattern: conflict graph structure
        conflict_density: fraction of possible conflicts that exist
        weight_pattern: item weight distribution
    """
    rng = random.Random(seed)

    # Item weights
    if weight_pattern == "uniform":
        weights = [rng.randint(capacity // 10, capacity // 3) for _ in range(n_items)]
    elif weight_pattern == "diverse":
        # Mix of small, medium, large items
        weights = []
        for i in range(n_items):
            if i % 3 == 0:
                weights.append(rng.randint(capacity // 20, capacity // 10))  # small
            elif i % 3 == 1:
                weights.append(rng.randint(capacity // 5, capacity // 3))    # medium
            else:
                weights.append(rng.randint(capacity // 2, 2*capacity // 3))  # large
    else:  # "tight"
        # Most items close to capacity/3 for tighter packing
        weights = [rng.randint(capacity // 4, capacity // 3) for _ in range(n_items)]

    # Conflict edges: set of (i,j) pairs with i < j
    conflicts = set()

    if conflict_pattern == "random":
        # Random conflicts with given density
        for i in range(n_items):
            for j in range(i+1, n_items):
                if rng.random() < conflict_density:
                    conflicts.add((i, j))

    elif conflict_pattern == "clique":
        # Items partitioned into k cliques (items in same clique conflict)
        k = max(2, int(math.sqrt(n_items)))
        items_per_clique = n_items // k
        for clique_id in range(k):
            start = clique_id * items_per_clique
            end = min((clique_id + 1) * items_per_clique, n_items)
            items = list(range(start, end))
            # All pairs within clique conflict
            for idx1 in range(len(items)):
                for idx2 in range(idx1+1, len(items)):
                    conflicts.add((items[idx1], items[idx2]))

    elif conflict_pattern == "bipartite":
        # Items split into two types; conflicts only across types
        split = n_items // 2
        type_a = list(range(split))
        type_b = list(range(split, n_items))
        for i in type_a:
            for j in type_b:
                if rng.random() < conflict_density * 2:  # compensate for reduced pairs
                    conflicts.add((i, j))

    elif conflict_pattern == "geometric":
        # Items have positions; conflicts if distance < threshold
        positions = [(rng.random() * 100, rng.random() * 100) for _ in range(n_items)]
        threshold = 100 * math.sqrt(conflict_density / math.pi)  # area-based density
        for i in range(n_items):
            for j in range(i+1, n_items):
                dx = positions[i][0] - positions[j][0]
                dy = positions[i][1] - positions[j][1]
                dist = math.sqrt(dx*dx + dy*dy)
                if dist < threshold:
                    conflicts.add((i, j))

    return weights, list(conflicts)


def write_bpc_mps(filename, n_items, capacity, weights, conflicts):
    """Write bin packing with conflicts in MPS format.

    We use a 2-index formulation with a fixed number of bins.
    Upper bound on bins needed: n_items (trivial worst case).
    Better bound: sum(weights) / capacity + k where k accounts for conflicts.
    """
    # Upper bound on bins: sum of weights / capacity (rounded up) + conflict slack
    n_bins = max(n_items, (sum(weights) + capacity - 1) // capacity + len(conflicts) // 10)
    n_bins = min(n_bins, n_items)  # Cap at n_items

    with open(filename, 'w') as fp:
        fp.write("NAME          BPC\n")

        # ROWS section
        fp.write("ROWS\n")
        fp.write(" N  OBJ\n")

        # Assignment: each item in exactly one bin
        for i in range(n_items):
            fp.write(f" E  ASSIGN_{i}\n")

        # Capacity: total weight in bin b <= capacity * y[b]
        for b in range(n_bins):
            fp.write(f" L  CAP_{b}\n")

        # Conflict: x[i][b] + x[j][b] <= y[b] for each conflict (i,j)
        for i, j in conflicts:
            for b in range(n_bins):
                fp.write(f" L  CONF_{i}_{j}_{b}\n")

        # COLUMNS section
        fp.write("COLUMNS\n")

        # Item assignment variables x[i][b]
        for i in range(n_items):
            for b in range(n_bins):
                var_name = f"x_{i}_{b}"
                # Assignment constraint
                fp.write(f"    {var_name:<10s}  ASSIGN_{i}  1.0\n")
                # Capacity constraint
                fp.write(f"    {var_name:<10s}  CAP_{b}  {weights[i]}\n")
                # Conflict constraints
                for j in range(n_items):
                    if i < j and (i, j) in conflicts:
                        fp.write(f"    {var_name:<10s}  CONF_{i}_{j}_{b}  1.0\n")
                    elif j < i and (j, i) in conflicts:
                        fp.write(f"    {var_name:<10s}  CONF_{j}_{i}_{b}  1.0\n")

        # Bin usage variables y[b] - write all coefficients in one block per variable
        for b in range(n_bins):
            var_name = f"y_{b}"
            # Objective
            fp.write(f"    {var_name:<10s}  OBJ       1.0\n")
            # Capacity constraint
            fp.write(f"    {var_name:<10s}  CAP_{b}  {-capacity}\n")
            # Conflict constraints
            for i, j in conflicts:
                fp.write(f"    {var_name:<10s}  CONF_{i}_{j}_{b}  -1.0\n")

        # RHS section
        fp.write("RHS\n")
        for i in range(n_items):
            fp.write(f"    RHS1      ASSIGN_{i}  1.0\n")
        # Capacity and conflict RHS are 0 (already in canonical form)

        # BOUNDS section
        fp.write("BOUNDS\n")
        for b in range(n_bins):
            fp.write(f" BV BND1      y_{b}\n")
        for i in range(n_items):
            for b in range(n_bins):
                fp.write(f" BV BND1      x_{i}_{b}\n")

        fp.write("ENDATA\n")


def solve_bpc_external(n_items, capacity, weights, conflicts, timeout=120):
    """Solve bin packing with conflicts using external solver."""
    try:
        import gurobipy as gp
        from gurobipy import GRB
    except ImportError:
        print("WARNING: external solver not available, skipping solve", file=sys.stderr)
        return None, None, None

    n_bins = max(n_items, (sum(weights) + capacity - 1) // capacity + len(conflicts) // 10)
    n_bins = min(n_bins, n_items)

    env = gp.Env(empty=True)
    env.setParam("OutputFlag", 0)
    env.setParam("TimeLimit", timeout)
    env.start()

    m = gp.Model(env=env)

    # Variables
    y = [m.addVar(vtype=GRB.BINARY, obj=1.0, name=f"y_{b}") for b in range(n_bins)]
    x = [[m.addVar(vtype=GRB.BINARY, name=f"x_{i}_{b}")
          for b in range(n_bins)] for i in range(n_items)]

    # Constraints: each item in exactly one bin
    for i in range(n_items):
        m.addConstr(gp.quicksum(x[i][b] for b in range(n_bins)) == 1, name=f"assign_{i}")

    # Capacity constraints
    for b in range(n_bins):
        m.addConstr(gp.quicksum(weights[i] * x[i][b] for i in range(n_items)) <= capacity * y[b],
                   name=f"cap_{b}")

    # Conflict constraints
    for i, j in conflicts:
        for b in range(n_bins):
            m.addConstr(x[i][b] + x[j][b] <= y[b], name=f"conf_{i}_{j}_{b}")

    m.optimize()

    if m.Status == GRB.OPTIMAL:
        obj_val = m.ObjVal
        env.dispose()
        return obj_val, "optimal", n_bins
    elif m.Status == GRB.TIME_LIMIT:
        obj_val = m.ObjVal if m.SolCount > 0 else None
        env.dispose()
        return obj_val, "timelimit", n_bins
    else:
        env.dispose()
        return None, f"status_{m.Status}", n_bins


# Test specifications: (n_items, capacity, seed, conflict_pattern, conflict_density, weight_pattern, timeout)
SPECS = [
    # Small instances for quick CI
    (10,  100,  42,  "random",    0.2, "uniform", 20),
    (12,  150,  42,  "clique",    0.3, "diverse", 30),
    (15,  100,  137, "bipartite", 0.3, "uniform", 30),
    (15,  120,  42,  "geometric", 0.2, "tight",   40),

    # Medium instances
    (20,  150,  42,  "random",    0.2, "diverse", 90),
    (25,  200,  137, "clique",    0.3, "uniform", 120),
    (30,  180,  42,  "geometric", 0.2, "tight",   120),
]


def main():
    fixture_dir = Path(__file__).parent / "fixtures"
    fixture_dir.mkdir(exist_ok=True)

    results = []

    for n_items, capacity, seed, conf_pat, conf_dens, weight_pat, timeout in SPECS:
        name = f"bpc_n{n_items}_c{capacity}_sd{seed}_{conf_pat}_{weight_pat}"
        print(f"Generating {name}...", flush=True)

        weights, conflicts = generate_bpc_instance(
            n_items, capacity, seed, conf_pat, conf_dens, weight_pat
        )

        # Write MPS
        mps_path = f"/tmp/{name}.mps"
        write_bpc_mps(mps_path, n_items, capacity, weights, conflicts)

        # Solve with external solver
        obj, status, n_bins = solve_bpc_external(n_items, capacity, weights, conflicts, timeout)

        if status == "optimal":
            print(f"  {name}: obj={obj:.0f} bins (optimal, n_bins={n_bins}, conflicts={len(conflicts)})")
            results.append((name, n_items, capacity, obj, timeout))

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
    print("Add to test/CInterfaceTest_bpc.c:")
    print("="*70)
    print("typedef struct {")
    print("  const char *name;")
    print("  double expected_obj;")
    print("  int timeout_sec;")
    print("} BpcTestCase;\n")
    print("static const BpcTestCase bpc_test_cases[] = {")
    for name, n_items, capacity, obj, timeout in results:
        print(f'  {{"{name}", {obj:.0f}, {timeout}}},')
    print("};")


if __name__ == "__main__":
    main()
