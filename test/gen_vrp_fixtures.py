#!/usr/bin/env python3
"""
Generate diverse CVRP and VRPPD test fixtures with various characteristics.

Characteristics explored:
  - Graph density (sparse vs dense)
  - Demand patterns (uniform, varied, clustered)
  - Capacity tightness (loose, medium, tight)
  - Geographic structure (clustered, uniform, grid)
  - Problem size (small, medium, large)
  - Multi-vehicle coordination needs

All instances use the correct capacity formulation:
  u[i] - u[j] + Q*x[i,j] <= Q - d[j]
"""

import random
import math


def euclidean_distance(p1, p2):
    """Euclidean distance between two points."""
    return math.sqrt((p1[0] - p2[0])**2 + (p1[1] - p2[1])**2)


def generate_positions(n, pattern="uniform", seed=None):
    """
    Generate node positions.

    Args:
      n: number of positions
      pattern: "uniform", "clustered", "grid", "ring"
      seed: random seed

    Returns:
      dict {node_id: (x, y)}
    """
    if seed:
        random.seed(seed)

    positions = {}

    if pattern == "uniform":
        # Uniformly distributed in 100x100 square
        for i in range(n):
            positions[i] = (random.uniform(0, 100), random.uniform(0, 100))

    elif pattern == "clustered":
        # 2-4 clusters of nodes
        n_clusters = random.randint(2, 4)
        cluster_centers = [(random.uniform(20, 80), random.uniform(20, 80))
                          for _ in range(n_clusters)]
        for i in range(n):
            center = cluster_centers[i % n_clusters]
            positions[i] = (center[0] + random.gauss(0, 10),
                           center[1] + random.gauss(0, 10))

    elif pattern == "grid":
        # Arrange in grid
        side = int(math.ceil(math.sqrt(n)))
        spacing = 100 / (side + 1)
        for i in range(n):
            row, col = divmod(i, side)
            positions[i] = ((col + 1) * spacing, (row + 1) * spacing)

    elif pattern == "ring":
        # Arrange in ring around depot
        for i in range(n):
            angle = 2 * math.pi * i / n
            radius = 50 + random.uniform(-5, 5)
            positions[i] = (50 + radius * math.cos(angle),
                           50 + radius * math.sin(angle))

    return positions


def compute_cost_matrix(positions, cost_type="euclidean"):
    """
    Compute cost matrix from positions.

    Args:
      positions: dict {node_id: (x, y)}
      cost_type: "euclidean", "manhattan", "rounded"
    """
    costs = {}
    nodes = sorted(positions.keys())

    for i in nodes:
        for j in nodes:
            if i == j:
                continue

            if cost_type == "euclidean":
                costs[(i,j)] = euclidean_distance(positions[i], positions[j])
            elif cost_type == "manhattan":
                costs[(i,j)] = abs(positions[i][0] - positions[j][0]) + \
                              abs(positions[i][1] - positions[j][1])
            elif cost_type == "rounded":
                dist = euclidean_distance(positions[i], positions[j])
                costs[(i,j)] = round(dist)

    return costs


def generate_demands(n, pattern="uniform", total=None, seed=None):
    """
    Generate customer demands.

    Args:
      n: number of customers
      pattern: "uniform", "varied", "skewed", "few_large"
      total: target total demand (if None, generated)
      seed: random seed

    Returns:
      dict {customer_id: demand} (customer IDs start at 1)
    """
    if seed:
        random.seed(seed)

    demands = {}

    if pattern == "uniform":
        # All demands similar
        base = 5 if total is None else total // n
        for i in range(1, n+1):
            demands[i] = base + random.randint(-1, 1)

    elif pattern == "varied":
        # Wide range of demands
        for i in range(1, n+1):
            demands[i] = random.randint(2, 10)

    elif pattern == "skewed":
        # Most small, few large
        for i in range(1, n+1):
            if random.random() < 0.2:  # 20% large
                demands[i] = random.randint(8, 12)
            else:
                demands[i] = random.randint(2, 4)

    elif pattern == "few_large":
        # Few very large demands
        for i in range(1, n+1):
            if random.random() < 0.1:  # 10% very large
                demands[i] = random.randint(15, 20)
            else:
                demands[i] = random.randint(2, 5)

    # Adjust to match target total if specified
    if total is not None:
        current = sum(demands.values())
        if current > 0:
            scale = total / current
            for i in demands:
                demands[i] = max(1, round(demands[i] * scale))

    return demands


def write_cvrp_mps(name, K, customers, Q, demands, costs, desc):
    """Write CVRP instance to MPS file."""

    N = len(customers)

    print(f"Generating {name}.mps")
    print(f"  {desc}")
    print(f"  Vehicles: {K}, Customers: {N}, Capacity: {Q}")
    print(f"  Demand: total={sum(demands[i] for i in customers)}, "
          f"min={min(demands.values())}, max={max(demands.values())}, "
          f"avg={sum(demands.values())/len(demands):.1f}")

    with open(f"{name}.mps", "w") as f:
        f.write(f"NAME          {name}\n")
        f.write("ROWS\n")
        f.write(" N  OBJ\n")

        for i in customers:
            f.write(f" E  OUT{i}\n")
            f.write(f" E  IN{i}\n")
            f.write(f" E  FLOW{i}\n")

        f.write(" L  FLEET\n")

        cap_idx = 0
        cap_map = {}
        for i in [0] + customers:
            for j in customers:
                if i == j:
                    continue
                cap_map[(i,j)] = cap_idx
                f.write(f" L  CAP{cap_idx}\n")
                cap_idx += 1

        for i in customers:
            f.write(f" G  DLBO{i}\n")
            f.write(f" L  DUBO{i}\n")

        f.write("COLUMNS\n")
        f.write("    MARK0000  'MARKER'                 'INTORG'\n")

        for i in [0] + customers:
            for j in [0] + customers:
                if i == j:
                    continue
                v = f"x{i}_{j}"

                f.write(f"    {v:12s}  OBJ       {costs.get((i,j), 9999):.2f}\n")

                if i > 0:
                    f.write(f"    {v:12s}  OUT{i}     1\n")
                if j > 0:
                    f.write(f"    {v:12s}  IN{j}      1\n")
                if i > 0:
                    f.write(f"    {v:12s}  FLOW{i}    1\n")
                if j > 0:
                    f.write(f"    {v:12s}  FLOW{j}   -1\n")
                if i == 0:
                    f.write(f"    {v:12s}  FLEET      1\n")

                if (i,j) in cap_map:
                    idx = cap_map[(i,j)]
                    f.write(f"    {v:12s}  CAP{idx}    {Q}\n")

        f.write("    MARK0000  'MARKER'                 'INTEND'\n")

        for i in customers:
            v = f"u{i}"
            for (arc_i, arc_j) in cap_map:
                idx = cap_map[(arc_i, arc_j)]
                if i == arc_i and arc_i > 0:
                    f.write(f"    {v:12s}  CAP{idx}     1\n")
                elif i == arc_j:
                    f.write(f"    {v:12s}  CAP{idx}    -1\n")
            f.write(f"    {v:12s}  DLBO{i}    1\n")
            f.write(f"    {v:12s}  DUBO{i}    1\n")

        f.write("RHS\n")
        for i in customers:
            f.write(f"    RHS1      OUT{i}      1\n")
            f.write(f"    RHS1      IN{i}       1\n")
        f.write(f"    RHS1      FLEET      {K}\n")

        for (i, j), idx in cap_map.items():
            rhs = Q - demands[j]
            f.write(f"    RHS1      CAP{idx}     {rhs}\n")

        for i in customers:
            f.write(f"    RHS1      DLBO{i}    {demands[i]}\n")
            f.write(f"    RHS1      DUBO{i}    {Q}\n")

        f.write("BOUNDS\n")
        f.write("ENDATA\n")

    print(f"✓ Generated {name}.mps\n")


def write_vrppd_mps(name, K, requests, Q, pickup_demand, costs, desc):
    """Write VRPPD instance to MPS file."""

    R = requests
    customers = list(range(1, 2*R+1))
    demands = {}
    for r in range(1, R+1):
        demands[r] = pickup_demand
        demands[r+R] = -pickup_demand

    print(f"Generating {name}.mps (VRPPD)")
    print(f"  {desc}")
    print(f"  Vehicles: {K}, Requests: {R}, Capacity: {Q}")
    print(f"  Pickup demand: {pickup_demand}")

    with open(f"{name}.mps", "w") as f:
        f.write(f"NAME          {name}\n")
        f.write("ROWS\n")
        f.write(" N  OBJ\n")

        for i in customers:
            f.write(f" E  OUT{i}\n")
            f.write(f" E  IN{i}\n")
            f.write(f" E  FLOW{i}\n")

        f.write(" L  FLEET\n")

        cap_idx = 0
        cap_map = {}
        for i in [0] + customers:
            for j in customers:
                if i == j:
                    continue
                cap_map[(i,j)] = cap_idx
                f.write(f" L  CAP{cap_idx}\n")
                cap_idx += 1

        for i in customers:
            f.write(f" G  DLBO{i}\n")
            f.write(f" L  DUBO{i}\n")

        for r in range(1, R+1):
            f.write(f" L  PREC{r}\n")

        f.write("COLUMNS\n")
        f.write("    MARK0000  'MARKER'                 'INTORG'\n")

        for i in [0] + customers:
            for j in [0] + customers:
                if i == j:
                    continue
                v = f"x{i}_{j}"

                f.write(f"    {v:12s}  OBJ       {costs.get((i,j), 9999):.2f}\n")

                if i > 0:
                    f.write(f"    {v:12s}  OUT{i}     1\n")
                if j > 0:
                    f.write(f"    {v:12s}  IN{j}      1\n")
                if i > 0:
                    f.write(f"    {v:12s}  FLOW{i}    1\n")
                if j > 0:
                    f.write(f"    {v:12s}  FLOW{j}   -1\n")
                if i == 0:
                    f.write(f"    {v:12s}  FLEET      1\n")

                if (i,j) in cap_map:
                    idx = cap_map[(i,j)]
                    f.write(f"    {v:12s}  CAP{idx}    {Q}\n")

                for r in range(1, R+1):
                    if i == r+R and j == r:
                        f.write(f"    {v:12s}  PREC{r}     1\n")

        f.write("    MARK0000  'MARKER'                 'INTEND'\n")

        for i in customers:
            v = f"u{i}"
            for (arc_i, arc_j) in cap_map:
                idx = cap_map[(arc_i, arc_j)]
                if i == arc_i and arc_i > 0:
                    f.write(f"    {v:12s}  CAP{idx}     1\n")
                elif i == arc_j:
                    f.write(f"    {v:12s}  CAP{idx}    -1\n")
            f.write(f"    {v:12s}  DLBO{i}    1\n")
            f.write(f"    {v:12s}  DUBO{i}    1\n")

        f.write("RHS\n")
        for i in customers:
            f.write(f"    RHS1      OUT{i}      1\n")
            f.write(f"    RHS1      IN{i}       1\n")
        f.write(f"    RHS1      FLEET      {K}\n")

        for (i, j), idx in cap_map.items():
            rhs = Q - demands[j]
            f.write(f"    RHS1      CAP{idx}     {rhs}\n")

        for i in customers:
            f.write(f"    RHS1      DLBO{i}    {demands[i]}\n")
            f.write(f"    RHS1      DUBO{i}    {Q}\n")

        for r in range(1, R+1):
            f.write(f"    RHS1      PREC{r}     0\n")

        f.write("BOUNDS\n")
        f.write("ENDATA\n")

    print(f"✓ Generated {name}.mps\n")


if __name__ == '__main__':
    print("=" * 70)
    print("CVRP TEST FIXTURES - Diverse Characteristics")
    print("=" * 70)
    print()

    # ========== CVRP: Graph Density ==========
    print("--- Graph Density ---")

    # Sparse graph (grid layout)
    pos = generate_positions(9, pattern="grid", seed=100)
    costs = compute_cost_matrix(pos, cost_type="rounded")
    demands = generate_demands(9, pattern="uniform", total=30, seed=100)
    write_cvrp_mps("cvrp_sparse_grid", K=2, customers=list(range(1, 10)),
                   Q=20, demands=demands, costs=costs,
                   desc="Sparse: grid layout, clear structure")

    # Dense graph (uniform random)
    pos = generate_positions(8, pattern="uniform", seed=101)
    costs = compute_cost_matrix(pos, cost_type="euclidean")
    demands = generate_demands(8, pattern="uniform", total=28, seed=101)
    write_cvrp_mps("cvrp_dense_uniform", K=2, customers=list(range(1, 9)),
                   Q=20, demands=demands, costs=costs,
                   desc="Dense: uniform distribution")

    # ========== CVRP: Demand Patterns ==========
    print("--- Demand Patterns ---")

    # Uniform demands
    pos = generate_positions(10, pattern="uniform", seed=200)
    costs = compute_cost_matrix(pos, cost_type="rounded")
    demands = {i: 5 for i in range(1, 11)}
    write_cvrp_mps("cvrp_demand_uniform", K=2, customers=list(range(1, 11)),
                   Q=30, demands=demands, costs=costs,
                   desc="Uniform demands: all customers similar")

    # Skewed demands (most small, few large)
    demands = generate_demands(10, pattern="skewed", seed=201)
    write_cvrp_mps("cvrp_demand_skewed", K=3, customers=list(range(1, 11)),
                   Q=25, demands=demands, costs=costs,
                   desc="Skewed demands: 80% small, 20% large")

    # Few very large demands
    demands = generate_demands(10, pattern="few_large", seed=202)
    write_cvrp_mps("cvrp_demand_outliers", K=3, customers=list(range(1, 11)),
                   Q=30, demands=demands, costs=costs,
                   desc="Outlier demands: 10% very large customers")

    # ========== CVRP: Capacity Tightness ==========
    print("--- Capacity Tightness ---")

    pos = generate_positions(12, pattern="clustered", seed=300)
    costs = compute_cost_matrix(pos, cost_type="rounded")
    demands = generate_demands(12, pattern="varied", seed=300)
    total_demand = sum(demands.values())

    # Very tight (forces many vehicles)
    Q_tight = max(demands.values()) + 2
    K_tight = (total_demand + Q_tight - 1) // Q_tight
    write_cvrp_mps("cvrp_capacity_tight", K=K_tight, customers=list(range(1, 13)),
                   Q=Q_tight, demands=demands, costs=costs,
                   desc=f"Very tight capacity: Q={Q_tight}, forces {K_tight} vehicles")

    # Medium (requires planning)
    Q_medium = total_demand // 3 + 5
    write_cvrp_mps("cvrp_capacity_medium", K=3, customers=list(range(1, 13)),
                   Q=Q_medium, demands=demands, costs=costs,
                   desc=f"Medium capacity: Q={Q_medium}, needs smart routing")

    # Loose (mostly unconstrained)
    Q_loose = total_demand
    write_cvrp_mps("cvrp_capacity_loose", K=1, customers=list(range(1, 13)),
                   Q=Q_loose, demands=demands, costs=costs,
                   desc=f"Loose capacity: Q={Q_loose}, single vehicle possible")

    # ========== CVRP: Geographic Structure ==========
    print("--- Geographic Structure ---")

    # Clustered customers
    pos = generate_positions(15, pattern="clustered", seed=400)
    costs = compute_cost_matrix(pos, cost_type="euclidean")
    demands = generate_demands(15, pattern="uniform", total=45, seed=400)
    write_cvrp_mps("cvrp_geo_clustered", K=3, customers=list(range(1, 16)),
                   Q=20, demands=demands, costs=costs,
                   desc="Clustered geography: 2-4 customer clusters")

    # Ring layout
    pos = generate_positions(12, pattern="ring", seed=401)
    costs = compute_cost_matrix(pos, cost_type="rounded")
    demands = generate_demands(12, pattern="uniform", total=36, seed=401)
    write_cvrp_mps("cvrp_geo_ring", K=2, customers=list(range(1, 13)),
                   Q=25, demands=demands, costs=costs,
                   desc="Ring geography: customers around depot")

    # ========== CVRP: Problem Size ==========
    print("--- Problem Size ---")

    # Small
    pos = generate_positions(6, pattern="uniform", seed=500)
    costs = compute_cost_matrix(pos, cost_type="rounded")
    demands = generate_demands(6, pattern="uniform", total=18, seed=500)
    write_cvrp_mps("cvrp_size_small", K=1, customers=list(range(1, 7)),
                   Q=20, demands=demands, costs=costs,
                   desc="Small: 6 customers, quick to solve")

    # Large
    pos = generate_positions(20, pattern="uniform", seed=501)
    costs = compute_cost_matrix(pos, cost_type="rounded")
    demands = generate_demands(20, pattern="varied", seed=501)
    write_cvrp_mps("cvrp_size_large", K=4, customers=list(range(1, 21)),
                   Q=30, demands=demands, costs=costs,
                   desc="Large: 20 customers, challenging")

    # ========== VRPPD ==========
    print("=" * 70)
    print("VRPPD TEST FIXTURES")
    print("=" * 70)
    print()

    # Small VRPPD with uniform layout
    pos = generate_positions(5, pattern="uniform", seed=600)
    costs = compute_cost_matrix(pos, cost_type="rounded")
    write_vrppd_mps("vrppd_small_uniform", K=1, requests=2, Q=15,
                    pickup_demand=5, costs=costs,
                    desc="Small uniform: 2 requests, easy")

    # Medium VRPPD with clustered layout
    pos = generate_positions(7, pattern="clustered", seed=601)
    costs = compute_cost_matrix(pos, cost_type="euclidean")
    write_vrppd_mps("vrppd_medium_clustered", K=2, requests=3, Q=12,
                    pickup_demand=4, costs=costs,
                    desc="Medium clustered: 3 requests, 2 vehicles")

    # Tight capacity VRPPD
    pos = generate_positions(9, pattern="ring", seed=602)
    costs = compute_cost_matrix(pos, cost_type="rounded")
    write_vrppd_mps("vrppd_tight_ring", K=3, requests=4, Q=10,
                    pickup_demand=4, costs=costs,
                    desc="Tight ring: 4 requests, tight capacity")

    # Large VRPPD
    pos = generate_positions(13, pattern="uniform", seed=603)
    costs = compute_cost_matrix(pos, cost_type="euclidean")
    write_vrppd_mps("vrppd_large_uniform", K=3, requests=6, Q=18,
                    pickup_demand=6, costs=costs,
                    desc="Large uniform: 6 requests, challenging")

    print("=" * 70)
    print("SUMMARY")
    print("=" * 70)
    print("\nGenerated fixtures covering:")
    print("  ✓ Graph density (sparse/dense)")
    print("  ✓ Demand patterns (uniform/skewed/outliers)")
    print("  ✓ Capacity tightness (tight/medium/loose)")
    print("  ✓ Geographic structure (clustered/ring/uniform)")
    print("  ✓ Problem size (small/medium/large)")
    print("  ✓ CVRP and VRPPD variants")
    print("\nAll instances use correct capacity formulation.")
    print("Ready for comprehensive testing!")
