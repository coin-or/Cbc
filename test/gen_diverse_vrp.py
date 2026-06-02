#!/usr/bin/env python3
"""
Generate VRP instances with truly diverse characteristics AND objective scales.

Key improvements over gen_vrp_fixtures.py:
- Different cost scales (small, medium, large values)
- More varied problem structures
- Clearer differentiation between instances
"""

import sys
sys.path.append('.')
from gen_vrp_fixtures import (generate_positions, compute_cost_matrix,
                               generate_demands, write_cvrp_mps, write_vrppd_mps)

def scale_costs(costs, scale):
    """Scale all costs by a factor."""
    return {(i,j): costs[(i,j)] * scale for (i,j) in costs}

if __name__ == '__main__':
    print("=" * 70)
    print("DIVERSE VRP INSTANCES - Different Scales and Structures")
    print("=" * 70)
    print()

    # ========== Small Scale (costs 1-20) ==========
    print("--- Small Scale Instances ---")

    pos = generate_positions(5, pattern="uniform", seed=1000)
    costs = compute_cost_matrix(pos, cost_type="rounded")
    costs = scale_costs(costs, 0.1)  # Scale down to 1-20 range
    demands = {1: 3, 2: 4, 3: 5, 4: 2, 5: 3}
    write_cvrp_mps("cvrp_small_scale", K=1, customers=[1,2,3,4,5],
                   Q=20, demands=demands, costs=costs,
                   desc="Small scale: costs 1-20, easy problem")

    # ========== Medium Scale (costs 10-100) ==========
    print("--- Medium Scale Instances ---")

    pos = generate_positions(8, pattern="grid", seed=1001)
    costs = compute_cost_matrix(pos, cost_type="euclidean")
    demands = generate_demands(8, pattern="uniform", total=32, seed=1001)
    write_cvrp_mps("cvrp_medium_scale", K=2, customers=list(range(1,9)),
                   Q=20, demands=demands, costs=costs,
                   desc="Medium scale: costs 10-100, grid structure")

    # ========== Large Scale (costs 100-500) ==========
    print("--- Large Scale Instances ---")

    pos = generate_positions(10, pattern="clustered", seed=1002)
    costs = compute_cost_matrix(pos, cost_type="euclidean")
    costs = scale_costs(costs, 5.0)  # Scale up
    demands = generate_demands(10, pattern="skewed", seed=1002)
    write_cvrp_mps("cvrp_large_scale", K=3, customers=list(range(1,11)),
                   Q=25, demands=demands, costs=costs,
                   desc="Large scale: costs 100-500, clustered")

    # ========== Asymmetric Costs ==========
    print("--- Asymmetric Costs ---")

    pos = generate_positions(7, pattern="ring", seed=1003)
    costs = compute_cost_matrix(pos, cost_type="euclidean")
    # Make asymmetric: one direction is 1.5x more expensive
    asym_costs = {}
    for (i,j) in costs:
        if i < j:
            asym_costs[(i,j)] = costs[(i,j)] * 1.5
            asym_costs[(j,i)] = costs[(i,j)]
        elif i > j and (i,j) not in asym_costs:
            asym_costs[(i,j)] = costs[(i,j)]
    demands = {i: 4 for i in range(1,8)}
    write_cvrp_mps("cvrp_asymmetric", K=2, customers=list(range(1,8)),
                   Q=18, demands=demands, costs=asym_costs,
                   desc="Asymmetric costs: c(i,j) != c(j,i)")

    # ========== High Capacity Utilization ==========
    print("--- High Capacity Utilization ---")

    customers = list(range(1, 13))
    demands = {i: 6 + (i % 3) for i in customers}  # Demands: 6,7,8,6,7,8,...
    total_demand = sum(demands.values())
    Q = (total_demand // 3) + 1  # Very tight: ~3 vehicles needed
    pos = generate_positions(12, pattern="uniform", seed=1004)
    costs = compute_cost_matrix(pos, cost_type="rounded")
    write_cvrp_mps("cvrp_high_utilization", K=3, customers=customers,
                   Q=Q, demands=demands, costs=costs,
                   desc=f"High utilization: Q={Q}, total demand={total_demand}")

    # ========== Many Small Deliveries ==========
    print("--- Many Small Deliveries ---")

    customers = list(range(1, 16))
    demands = {i: 2 for i in customers}  # All small
    pos = generate_positions(15, pattern="clustered", seed=1005)
    costs = compute_cost_matrix(pos, cost_type="euclidean")
    write_cvrp_mps("cvrp_many_small", K=2, customers=customers,
                   Q=20, demands=demands, costs=costs,
                   desc="Many small: 15 customers with demand=2")

    # ========== Few Large Deliveries ==========
    print("--- Few Large Deliveries ---")

    customers = list(range(1, 6))
    demands = {i: 15 + i for i in customers}  # Large: 16,17,18,19,20
    pos = generate_positions(5, pattern="uniform", seed=1006)
    costs = compute_cost_matrix(pos, cost_type="rounded")
    write_cvrp_mps("cvrp_few_large", K=5, customers=customers,
                   Q=25, demands=demands, costs=costs,
                   desc="Few large: 5 customers with large demands (16-20)")

    # ========== VRPPD with Different Scales ==========
    print("\n--- VRPPD Diverse Scale ---")

    # Small scale VRPPD
    pos = generate_positions(5, pattern="grid", seed=2000)
    costs = compute_cost_matrix(pos, cost_type="rounded")
    costs = scale_costs(costs, 0.2)
    write_vrppd_mps("vrppd_small_scale", K=1, requests=2, Q=12,
                    pickup_demand=4, costs=costs,
                    desc="Small scale VRPPD: costs 2-30")

    # Large scale VRPPD with high demand
    pos = generate_positions(9, pattern="uniform", seed=2001)
    costs = compute_cost_matrix(pos, cost_type="euclidean")
    costs = scale_costs(costs, 3.0)
    write_vrppd_mps("vrppd_large_scale", K=2, requests=4, Q=20,
                    pickup_demand=8, costs=costs,
                    desc="Large scale VRPPD: costs 30-450, high demand")

    # Tight VRPPD with many requests
    pos = generate_positions(11, pattern="clustered", seed=2002)
    costs = compute_cost_matrix(pos, cost_type="euclidean")
    write_vrppd_mps("vrppd_many_requests", K=3, requests=5, Q=15,
                    pickup_demand=5, costs=costs,
                    desc="Many requests VRPPD: 5 requests, tight coordination")

    print("\n" + "=" * 70)
    print("Generated 10 diverse instances with varied:")
    print("  ✓ Cost scales (small: 1-20, medium: 10-100, large: 100-500)")
    print("  ✓ Problem structures (symmetric/asymmetric, utilization)")
    print("  ✓ Delivery patterns (many small vs few large)")
    print("  ✓ VRPPD variations (small/large scale, many requests)")
    print("=" * 70)
