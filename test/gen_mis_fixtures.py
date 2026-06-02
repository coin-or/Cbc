#!/usr/bin/env python3
"""
Generate Maximum (Weighted) Independent Set (MIS/MWIS) fixtures.

Problem:
  maximize   sum_v w[v] * x[v]
  subject to x[u] + x[v] <= 1   for all edges (u,v) in E
             x[v] in {0,1}

where:
  - w[v] is the weight (value) of vertex v
  - E is the edge set
  - x[v] = 1 if vertex v is in the independent set

This is the canonical pure binary problem for testing clique cuts:
  - Independent set in G = clique in complement graph G'
  - Clique cuts correspond to odd-hole inequalities
  - Conflict graph IS the problem structure

Graph patterns:
  - random: Erdős-Rényi random graphs
  - geometric: unit disk graphs (proximity in 2D)
  - grid: grid graphs (2D lattice)
  - k_partite: k-partite graphs (known structure)
  - trees: tree graphs (polynomial-time solvable, good for validation)
  - planar: planar graphs (special structure)

Solved with external solver to obtain certified optimal values.
"""

import sys
import math
import random
import gzip
from pathlib import Path


def generate_graph(n_vertices, seed, pattern="random", edge_prob=0.3, k=3):
    """Generate a graph for the independent set problem.

    Args:
        n_vertices: number of vertices
        seed: random seed
        pattern: graph structure pattern
        edge_prob: edge probability for random graphs
        k: parameter for structured graphs (e.g., k-partite)

    Returns:
        edges: list of (u, v) tuples with u < v
    """
    rng = random.Random(seed)
    edges = set()

    if pattern == "random":
        # Erdős-Rényi random graph
        for u in range(n_vertices):
            for v in range(u+1, n_vertices):
                if rng.random() < edge_prob:
                    edges.add((u, v))

    elif pattern == "geometric":
        # Unit disk graph: vertices in 2D plane, edge if distance < threshold
        positions = [(rng.random() * 100, rng.random() * 100) for _ in range(n_vertices)]
        # Threshold chosen to achieve desired edge density
        threshold = 100 * math.sqrt(edge_prob / (math.pi * 0.7))
        for u in range(n_vertices):
            for v in range(u+1, n_vertices):
                dx = positions[u][0] - positions[v][0]
                dy = positions[u][1] - positions[v][1]
                dist = math.sqrt(dx*dx + dy*dy)
                if dist < threshold:
                    edges.add((u, v))

    elif pattern == "grid":
        # 2D grid graph (4-connected)
        side = int(math.sqrt(n_vertices))
        for u in range(n_vertices):
            row_u, col_u = u // side, u % side
            # Right neighbor
            if col_u < side - 1:
                v = u + 1
                if v < n_vertices:
                    edges.add((u, v))
            # Down neighbor
            if row_u < side - 1:
                v = u + side
                if v < n_vertices:
                    edges.add((u, v))

    elif pattern == "k_partite":
        # k-partite graph: vertices partitioned into k sets, edges only across sets
        part_size = n_vertices // k
        for u in range(n_vertices):
            part_u = u // part_size
            for v in range(u+1, n_vertices):
                part_v = v // part_size
                if part_u != part_v and rng.random() < edge_prob * 2:  # compensate
                    edges.add((u, v))

    elif pattern == "tree":
        # Random tree (spanning tree from random edges)
        # Tree has n-1 edges
        vertices = list(range(n_vertices))
        rng.shuffle(vertices)
        for i in range(1, n_vertices):
            # Connect vertex i to random earlier vertex
            j = rng.randint(0, i-1)
            u, v = min(vertices[i], vertices[j]), max(vertices[i], vertices[j])
            edges.add((u, v))

    elif pattern == "cycle":
        # Cycle graph
        for u in range(n_vertices):
            v = (u + 1) % n_vertices
            edges.add((min(u,v), max(u,v)))

    elif pattern == "planar":
        # Planar graph: start with grid, add some diagonals
        side = int(math.sqrt(n_vertices))
        # Grid edges
        for u in range(n_vertices):
            row_u, col_u = u // side, u % side
            if col_u < side - 1:
                v = u + 1
                if v < n_vertices:
                    edges.add((u, v))
            if row_u < side - 1:
                v = u + side
                if v < n_vertices:
                    edges.add((u, v))
        # Add some diagonals (keep planar)
        for u in range(n_vertices):
            row_u, col_u = u // side, u % side
            if row_u < side - 1 and col_u < side - 1 and rng.random() < 0.3:
                v = u + side + 1  # diagonal
                if v < n_vertices:
                    edges.add((u, v))

    return list(edges)


def generate_mis_instance(n_vertices, seed, pattern="random", edge_prob=0.3,
                         weight_type="uniform", k=3):
    """Generate a maximum independent set instance.

    Args:
        n_vertices: number of vertices
        seed: random seed
        pattern: graph structure
        edge_prob: edge probability
        weight_type: vertex weight distribution
        k: parameter for structured graphs
    """
    rng = random.Random(seed)

    # Vertex weights
    if weight_type == "uniform":
        weights = [rng.randint(1, 100) for _ in range(n_vertices)]
    elif weight_type == "unit":
        # Unweighted: all weights = 1 (classic max independent set)
        weights = [1] * n_vertices
    elif weight_type == "diverse":
        # Mix of small, medium, large weights
        weights = []
        for v in range(n_vertices):
            if v % 3 == 0:
                weights.append(rng.randint(1, 20))
            elif v % 3 == 1:
                weights.append(rng.randint(30, 70))
            else:
                weights.append(rng.randint(80, 150))
    elif weight_type == "degree":
        # Weights inversely proportional to degree (heuristic: high-degree nodes less valuable)
        edges = generate_graph(n_vertices, seed, pattern, edge_prob, k)
        degree = [0] * n_vertices
        for u, v in edges:
            degree[u] += 1
            degree[v] += 1
        weights = [max(1, 100 - degree[v] * 5) for v in range(n_vertices)]
        return weights, edges
    else:
        weights = [1] * n_vertices

    edges = generate_graph(n_vertices, seed, pattern, edge_prob, k)
    return weights, edges


def write_mis_mps(filename, n_vertices, weights, edges):
    """Write MIS instance in MPS format."""
    with open(filename, 'w') as fp:
        fp.write("NAME          MIS\n")

        # ROWS section
        fp.write("ROWS\n")
        fp.write(" N  OBJ\n")
        for u, v in edges:
            fp.write(f" L  EDGE_{u}_{v}\n")

        # COLUMNS section
        fp.write("COLUMNS\n")
        for v in range(n_vertices):
            var_name = f"x_{v}"
            # Objective coefficient (negate for maximization)
            fp.write(f"    {var_name:<10s}  OBJ       {-weights[v]}\n")
            # Edge constraints: x[u] + x[v] <= 1
            for u, w in edges:
                if w == v:
                    fp.write(f"    {var_name:<10s}  EDGE_{u}_{w}  1.0\n")
                elif u == v:
                    fp.write(f"    {var_name:<10s}  EDGE_{u}_{w}  1.0\n")

        # RHS section
        fp.write("RHS\n")
        for u, v in edges:
            fp.write(f"    RHS1      EDGE_{u}_{v}  1.0\n")

        # BOUNDS section
        fp.write("BOUNDS\n")
        for v in range(n_vertices):
            fp.write(f" BV BND1      x_{v}\n")

        fp.write("ENDATA\n")


def solve_mis_external(n_vertices, weights, edges, timeout=60):
    """Solve MIS with external solver to get certified optimal value."""
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

    # Variables: x[v] = 1 if vertex v in independent set
    x = [m.addVar(vtype=GRB.BINARY, obj=weights[v], name=f"x_{v}")
         for v in range(n_vertices)]

    m.setObjective(gp.quicksum(weights[v] * x[v] for v in range(n_vertices)), GRB.MAXIMIZE)

    # Constraints: x[u] + x[v] <= 1 for each edge
    for u, v in edges:
        m.addConstr(x[u] + x[v] <= 1, name=f"edge_{u}_{v}")

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


# Test specifications: (n_vertices, seed, pattern, edge_prob, weight_type, k, timeout)
SPECS = [
    # Small instances for quick CI
    (15,  42,  "random",    0.3, "uniform", 3, 10),
    (20,  42,  "geometric", 0.3, "unit",    3, 15),
    (16,  137, "grid",      0.0, "diverse", 3, 10),  # 4x4 grid
    (18,  42,  "k_partite", 0.4, "uniform", 3, 20),
    (20,  137, "tree",      0.0, "uniform", 3, 5),   # polynomial time

    # Medium instances
    (30,  42,  "random",    0.3, "uniform", 3, 60),
    (35,  137, "geometric", 0.3, "diverse", 3, 60),
    (36,  42,  "grid",      0.0, "unit",    3, 30),  # 6x6 grid
    (40,  137, "planar",    0.3, "uniform", 3, 90),
    (30,  42,  "cycle",     0.0, "uniform", 3, 10),
]


def main():
    fixture_dir = Path(__file__).parent / "fixtures"
    fixture_dir.mkdir(exist_ok=True)

    results = []

    for n_vert, seed, pattern, edge_prob, weight_type, k, timeout in SPECS:
        name = f"mis_n{n_vert}_sd{seed}_{pattern}_{weight_type}"
        print(f"Generating {name}...", flush=True)

        weights, edges = generate_mis_instance(n_vert, seed, pattern, edge_prob, weight_type, k)

        # Write MPS
        mps_path = f"/tmp/{name}.mps"
        write_mis_mps(mps_path, n_vert, weights, edges)

        # Solve with external solver
        obj, status = solve_mis_external(n_vert, weights, edges, timeout)

        if status == "optimal":
            print(f"  {name}: obj={obj:.2f} (optimal, |V|={n_vert}, |E|={len(edges)})")
            results.append((name, n_vert, len(edges), obj, timeout))

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
    print("Add to test/CInterfaceTest_mis.c:")
    print("="*70)
    print("typedef struct {")
    print("  const char *name;")
    print("  double expected_obj;")
    print("  int timeout_sec;")
    print("} MisTestCase;\n")
    print("static const MisTestCase mis_test_cases[] = {")
    for name, n_vert, n_edges, obj, timeout in results:
        print(f'  {{"{name}", {obj:.2f}, {timeout}}},')
    print("};")


if __name__ == "__main__":
    main()
