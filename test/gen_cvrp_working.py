#!/usr/bin/env python3
"""
Generate CVRP (Capacitated Vehicle Routing Problem) and VRPPD (with Pickup-Delivery)
test fixtures with various characteristics.

Correct capacity constraint formulation:
  u[i] - u[j] + Q*x[i,j] <= Q - d[j]

When x[i,j]=1: u[j] >= u[i] + d[j] (load changes by demand)
When x[i,j]=0: u[i] - u[j] <= Q - d[j] (non-binding)
"""

import random

def generate_cvrp(name, K, customers, Q, demands, costs, desc):
    """
    Generate CVRP instance.

    Args:
      name: output filename (without .mps)
      K: number of vehicles
      customers: list of customer IDs (depot is 0)
      Q: vehicle capacity
      demands: dict {customer_id: demand}
      costs: dict {(i,j): cost}
      desc: description printed to console
    """
    N = len(customers)

    print(f"Generating {name}.mps")
    print(f"  {desc}")
    print(f"  Vehicles: {K}, Customers: {N}, Capacity: {Q}")
    print(f"  Total demand: {sum(demands[i] for i in customers)}")

    with open(f"{name}.mps", "w") as f:
        f.write(f"NAME          {name}\n")
        f.write("ROWS\n")
        f.write(" N  OBJ\n")

        # Each customer: out-degree = 1
        for i in customers:
            f.write(f" E  OUT{i}\n")

        # Each customer: in-degree = 1
        for i in customers:
            f.write(f" E  IN{i}\n")

        # Flow conservation
        for i in customers:
            f.write(f" E  FLOW{i}\n")

        # Fleet limit
        f.write(" L  FLEET\n")

        # Capacity constraints: u[i] - u[j] + Q*x[i,j] <= Q - d[j]
        cap_idx = 0
        cap_map = {}
        for i in [0] + customers:
            for j in customers:
                if i == j:
                    continue
                cap_map[(i,j)] = cap_idx
                f.write(f" L  CAP{cap_idx}\n")
                cap_idx += 1

        # Demand bounds
        for i in customers:
            f.write(f" G  DLBO{i}\n")
            f.write(f" L  DUBO{i}\n")

        f.write("COLUMNS\n")
        f.write("    MARK0000  'MARKER'                 'INTORG'\n")

        # x[i,j] variables
        for i in [0] + customers:
            for j in [0] + customers:
                if i == j:
                    continue
                v = f"x{i}_{j}"

                f.write(f"    {v:12s}  OBJ       {costs.get((i,j), 9999)}\n")

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

        # u[i] variables
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


def generate_vrppd(name, K, requests, Q, pickup_demand, costs, desc):
    """
    Generate VRPPD instance.

    Args:
      name: output filename
      K: number of vehicles
      requests: number of pickup-delivery pairs
      Q: capacity
      pickup_demand: demand size for pickups (delivery is negative)
      costs: dict {(i,j): cost}
      desc: description

    Nodes: 0=depot, 1..R=pickups, R+1..2R=deliveries
    """
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

        # Capacity
        cap_idx = 0
        cap_map = {}
        for i in [0] + customers:
            for j in customers:
                if i == j:
                    continue
                cap_map[(i,j)] = cap_idx
                f.write(f" L  CAP{cap_idx}\n")
                cap_idx += 1

        # Demand bounds
        for i in customers:
            f.write(f" G  DLBO{i}\n")
            f.write(f" L  DUBO{i}\n")

        # Precedence: forbid delivery→pickup for each request
        for r in range(1, R+1):
            f.write(f" L  PREC{r}\n")

        f.write("COLUMNS\n")
        f.write("    MARK0000  'MARKER'                 'INTORG'\n")

        for i in [0] + customers:
            for j in [0] + customers:
                if i == j:
                    continue
                v = f"x{i}_{j}"

                f.write(f"    {v:12s}  OBJ       {costs.get((i,j), 9999)}\n")

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

                # Precedence: forbid delivery→pickup arcs
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
    # Instance 1: Small tight capacity (needs both vehicles)
    customers = [1, 2, 3, 4]
    demands = {1: 3, 2: 4, 3: 5, 4: 6}
    costs = {(i,j): 10 + abs(i-j)*5 for i in [0]+customers for j in [0]+customers if i!=j}
    generate_cvrp("cvrp_tight", K=2, customers=customers, Q=10, demands=demands, costs=costs,
                  desc="Tight capacity - total demand 18, Q=10, needs 2 vehicles")

    # Instance 2: Loose capacity (single vehicle sufficient)
    generate_cvrp("cvrp_loose", K=1, customers=customers, Q=20, demands=demands, costs=costs,
                  desc="Loose capacity - single vehicle can serve all")

    # Instance 3: Larger instance with mixed demands
    random.seed(42)
    customers = list(range(1, 11))
    demands = {i: random.randint(2, 8) for i in customers}
    costs = {(i,j): 10 + 5*((i-j)**2)**0.5 for i in [0]+customers for j in [0]+customers if i!=j}
    generate_cvrp("cvrp_medium", K=3, customers=customers, Q=15, demands=demands, costs=costs,
                  desc="Medium instance - 10 customers, 3 vehicles, mixed demands")

    # Instance 4: VRPPD small
    costs_pd = {(i,j): 10 + abs(i-j)*5 for i in range(5) for j in range(5) if i!=j}
    generate_vrppd("vrppd_small", K=1, requests=2, Q=20, pickup_demand=5, costs=costs_pd,
                   desc="Small VRPPD - 2 requests, 1 vehicle")

    # Instance 5: VRPPD tight
    generate_vrppd("vrppd_tight", K=2, requests=3, Q=10, pickup_demand=4, costs=costs_pd,
                   desc="Tight VRPPD - 3 requests, 2 vehicles, tight capacity")

    print("All instances generated successfully!")
