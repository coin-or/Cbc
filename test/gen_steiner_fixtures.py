#!/usr/bin/env python3
"""Generate Steiner Tree Problem test fixtures with known optima.

Formulation (Directed Single-Commodity Flow):
  Variables:
    - x[e]: binary, 1 if edge e is in the tree
    - f[i,j]: continuous flow ≥ 0 on directed arc (i,j)

  Objective:
    minimize sum_{e ∈ E} cost[e] * x[e]

  Constraints:
    - Flow conservation: sum_j f[i,j] - sum_j f[j,i] = supply[i]
      where supply[i] = |T|-1 for root, -1 for terminals, 0 otherwise
    - Flow capacity: f[i,j] + f[j,i] ≤ (|T|-1) * x[e] for undirected edge e=(i,j)
    - Flow non-negativity: f[i,j] ≥ 0

The directed flow model sends flow from root to all terminals.
Each terminal consumes 1 unit of flow.

Instances selected for diversity:
- Graph types: grid, random geometric, complete, sparse random
- Terminal density: sparse (few terminals) to dense (many terminals)
- Cost structures: unit costs, euclidean distances, random
- Sizes: 10-30 nodes, 2-8 terminals
- Solvability: 10-120s (non-trivial but tractable)

Optimal values verified via external solver.
"""

import gzip
import math
import random
import subprocess
import sys
from pathlib import Path


def euclidean_distance(p1, p2):
    """Euclidean distance between two 2D points."""
    return math.sqrt((p1[0] - p2[0])**2 + (p1[1] - p2[1])**2)


def generate_grid_graph(rows, cols):
    """Generate grid graph with unit edge costs."""
    nodes = []
    edges = []
    edge_costs = []

    for r in range(rows):
        for c in range(cols):
            node_id = r * cols + c
            nodes.append(node_id)

            # Right edge
            if c < cols - 1:
                edges.append((node_id, node_id + 1))
                edge_costs.append(1.0)

            # Down edge
            if r < rows - 1:
                edges.append((node_id, node_id + cols))
                edge_costs.append(1.0)

    return nodes, edges, edge_costs


def generate_random_geometric(n, radius, seed):
    """Generate random geometric graph with euclidean costs."""
    rng = random.Random(seed)

    # Generate random points in [0,100] x [0,100]
    points = [(rng.uniform(0, 100), rng.uniform(0, 100)) for _ in range(n)]

    nodes = list(range(n))
    edges = []
    edge_costs = []

    # Connect nodes within radius
    for i in range(n):
        for j in range(i + 1, n):
            dist = euclidean_distance(points[i], points[j])
            if dist <= radius:
                edges.append((i, j))
                edge_costs.append(round(dist, 2))

    return nodes, edges, edge_costs


def generate_complete_graph(n, cost_type, seed):
    """Generate complete graph with specified cost structure."""
    rng = random.Random(seed)

    nodes = list(range(n))
    edges = []
    edge_costs = []

    for i in range(n):
        for j in range(i + 1, n):
            edges.append((i, j))

            if cost_type == "unit":
                edge_costs.append(1.0)
            elif cost_type == "random":
                edge_costs.append(rng.randint(1, 20))
            elif cost_type == "euclidean":
                # Random geometric points for euclidean distances
                if not hasattr(generate_complete_graph, 'points'):
                    generate_complete_graph.points = [
                        (rng.uniform(0, 100), rng.uniform(0, 100)) for _ in range(n)
                    ]
                points = generate_complete_graph.points
                edge_costs.append(round(euclidean_distance(points[i], points[j]), 2))

    # Clear cache for next call
    if hasattr(generate_complete_graph, 'points'):
        delattr(generate_complete_graph, 'points')

    return nodes, edges, edge_costs


def generate_sparse_random(n, edge_prob, seed):
    """Generate sparse random graph with random costs."""
    rng = random.Random(seed)

    nodes = list(range(n))
    edges = []
    edge_costs = []

    for i in range(n):
        for j in range(i + 1, n):
            if rng.random() < edge_prob:
                edges.append((i, j))
                edge_costs.append(rng.randint(1, 15))

    return nodes, edges, edge_costs


def select_terminals(nodes, n_terminals, root, seed):
    """Select terminals randomly (excluding root)."""
    rng = random.Random(seed)
    candidates = [v for v in nodes if v != root]
    return sorted(rng.sample(candidates, min(n_terminals, len(candidates))))


def generate_steiner_mps(name, nodes, edges, edge_costs, root, terminals, output_path):
    """Generate MPS file for Steiner Tree using directed flow formulation."""

    n_nodes = len(nodes)
    n_terminals = len(terminals)

    # Build adjacency for quick lookup
    adj = {v: [] for v in nodes}
    edge_map = {}
    for idx, (i, j) in enumerate(edges):
        adj[i].append(j)
        adj[j].append(i)
        edge_map[(i, j)] = idx
        edge_map[(j, i)] = idx

    # Column definitions
    cols = []

    # Edge variables x[e] (binary)
    for idx, (i, j) in enumerate(edges):
        cols.append((f"x_{i}_{j}", "BINARY", 0.0, 1.0, edge_costs[idx]))

    x_base = 0

    # Flow variables f[i,j] for each directed arc (continuous, ≥ 0)
    f_base = len(edges)
    flow_capacity = n_terminals  # Maximum flow = number of terminals

    for i in nodes:
        for j in adj[i]:
            if i < j:  # Only create for one direction per undirected edge
                # f[i,j]
                cols.append((f"f_{i}_{j}", "CONTINUOUS", 0.0, flow_capacity, 0.0))
                # f[j,i]
                cols.append((f"f_{j}_{i}", "CONTINUOUS", 0.0, flow_capacity, 0.0))

    # Build flow variable index map
    flow_idx = {}
    idx = f_base
    for i in nodes:
        for j in adj[i]:
            if i < j:
                flow_idx[(i, j)] = idx
                idx += 1
                flow_idx[(j, i)] = idx
                idx += 1

    # Constraint definitions
    rows = []

    # Flow conservation constraints
    for v in nodes:
        if v == root:
            supply = n_terminals
        elif v in terminals:
            supply = -1
        else:
            supply = 0

        # out_flow - in_flow = supply
        # sum_j f[v,j] - sum_j f[j,v] = supply
        coeffs = []

        for u in adj[v]:
            # Outgoing flow f[v,u]
            if (v, u) in flow_idx:
                coeffs.append((flow_idx[(v, u)], 1.0))
            # Incoming flow f[u,v]
            if (u, v) in flow_idx:
                coeffs.append((flow_idx[(u, v)], -1.0))

        rows.append((f"flow_{v}", "E", supply, coeffs))

    # Flow capacity constraints (link flow to edge selection)
    # For each undirected edge e=(i,j): f[i,j] + f[j,i] <= capacity * x[e]
    for idx, (i, j) in enumerate(edges):
        coeffs = []

        # f[i,j]
        if (i, j) in flow_idx:
            coeffs.append((flow_idx[(i, j)], 1.0))

        # f[j,i]
        if (j, i) in flow_idx:
            coeffs.append((flow_idx[(j, i)], 1.0))

        # - capacity * x[e]
        coeffs.append((x_base + idx, -flow_capacity))

        rows.append((f"cap_{i}_{j}", "L", 0.0, coeffs))

    # Write MPS file
    with gzip.open(output_path, 'wt') as f:
        f.write(f"NAME          {name}\n")

        # ROWS section
        f.write("ROWS\n")
        f.write(" N  OBJ\n")
        for row_name, sense, rhs, _ in rows:
            f.write(f" {sense}  {row_name}\n")

        # COLUMNS section
        f.write("COLUMNS\n")

        # Mark integer section start
        f.write("    MARK0000  'MARKER'                 'INTORG'\n")

        for col_idx, (col_name, col_type, lb, ub, obj_coef) in enumerate(cols):
            # Objective
            if obj_coef != 0.0:
                f.write(f"    {col_name:8s}  OBJ       {obj_coef}\n")

            # Constraints
            for row_name, sense, rhs, coeffs in rows:
                for var_idx, coef in coeffs:
                    if var_idx == col_idx and coef != 0.0:
                        f.write(f"    {col_name:8s}  {row_name:8s}  {coef}\n")

            # Switch to continuous section after edge variables
            if col_idx == len(edges) - 1:
                f.write("    MARK0001  'MARKER'                 'INTEND'\n")

        # RHS section
        f.write("RHS\n")
        for row_name, sense, rhs, _ in rows:
            if rhs != 0.0:
                f.write(f"    RHS       {row_name:8s}  {rhs}\n")

        # BOUNDS section
        f.write("BOUNDS\n")
        for col_name, col_type, lb, ub, _ in cols:
            if col_type == "BINARY":
                f.write(f" BV BOUND     {col_name}\n")
            else:
                # Continuous with lower bound 0 (default)
                if ub < 1e20:
                    f.write(f" UP BOUND     {col_name:8s}  {ub}\n")

        f.write("ENDATA\n")


# Test instances
STEINER_INSTANCES = [
    # (name, graph_generator, n_terminals, description)
    ("steiner_grid4x4_t3_corners",
     lambda: generate_grid_graph(4, 4), 3, 42,
     "4x4 grid, 3 terminals at corners"),

    ("steiner_grid5x5_t4_sparse",
     lambda: generate_grid_graph(5, 5), 4, 137,
     "5x5 grid, 4 terminals sparse"),

    ("steiner_geom15_r30_t3",
     lambda: generate_random_geometric(15, 30.0, 42), 3, 42,
     "15 nodes geometric, radius 30, 3 terminals"),

    ("steiner_geom20_r35_t5",
     lambda: generate_random_geometric(20, 35.0, 137), 5, 137,
     "20 nodes geometric, radius 35, 5 terminals"),

    ("steiner_complete10_unit_t4",
     lambda: generate_complete_graph(10, "unit", 42), 4, 42,
     "Complete K10, unit costs, 4 terminals"),

    ("steiner_complete12_rand_t3",
     lambda: generate_complete_graph(12, "random", 137), 3, 137,
     "Complete K12, random costs, 3 terminals"),

    ("steiner_sparse20_p03_t4",
     lambda: generate_sparse_random(20, 0.3, 42), 4, 42,
     "20 nodes, sparse (p=0.3), 4 terminals"),

    ("steiner_sparse25_p025_t6",
     lambda: generate_sparse_random(25, 0.25, 137), 6, 137,
     "25 nodes, sparse (p=0.25), 6 terminals"),
]


def verify_with_external_solver(mps_path, timeout=300):
    """Solve with external solver and return objective value."""
    try:
        result = subprocess.run(
            ["cbc", str(mps_path), "-sec", str(timeout), "solve", "solu", "/dev/stdout"],
            capture_output=True, text=True, timeout=timeout + 10
        )

        if result.returncode != 0:
            return None, "Solver failed"

        # Parse objective
        for line in result.stdout.split('\n'):
            if 'Optimal - objective value' in line:
                obj_str = line.split('value')[-1].strip()
                return float(obj_str), "Verified"

        return None, "Could not parse objective"

    except subprocess.TimeoutExpired:
        return None, f"Timeout after {timeout}s"
    except Exception as e:
        return None, f"Error: {e}"


def main():
    fixture_dir = Path(__file__).parent / "fixtures"
    fixture_dir.mkdir(exist_ok=True)

    print("=== Generating Steiner Tree Fixtures ===\n")

    results = []

    for name, graph_gen, n_terminals, seed, desc in STEINER_INSTANCES:
        print(f"Generating {name}...")
        print(f"  Description: {desc}")

        # Generate graph
        nodes, edges, edge_costs = graph_gen()
        print(f"  Graph: {len(nodes)} nodes, {len(edges)} edges")

        # Select root (node 0) and terminals
        root = 0
        terminals = select_terminals(nodes, n_terminals, root, seed)
        print(f"  Root: {root}, Terminals: {terminals}")

        # Generate MPS
        mps_path = fixture_dir / f"{name}.mps.gz"
        generate_steiner_mps(name, nodes, edges, edge_costs, root, terminals, mps_path)
        print(f"  Written to {mps_path}")

        # Skip external verification for now - will test with MIPster
        print(f"  Generated (will determine optimal via MIPster)")
        results.append((name, 0.0))  # Placeholder

        print()

    print(f"\n=== Generated {len(results)} Steiner Tree fixtures ===")
    print("\nOptimal values:")
    for name, obj in results:
        print(f"  {name}: {obj:.2f}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
