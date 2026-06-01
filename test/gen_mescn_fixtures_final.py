#!/usr/bin/env python3
"""Final MESCN fixture generators for CI."""

import random
from collections import namedtuple

MESCNInstance = namedtuple('MESCNInstance', [
    'name', 'F', 'W', 'R', 'P', 'T',
    'setup_cost_f', 'setup_cost_w', 'prod_cost',
    'transport_fw', 'transport_wr',
    'holding_f', 'holding_w',
    'cap_factory', 'cap_warehouse',
    'batch_size', 'truck_size', 'demand',
    'fully_connected', 'fw_arcs', 'wr_arcs'
])

def generate_moderate_sparse(seed, F=3, W=4, R=5, P=2, T=3):
    """Sparse network instance - requires branching."""
    random.seed(seed)
    name = f"moderate_sparse_s{seed}"

    setup_cost_f = [random.randint(50, 100) for _ in range(F * P)]
    setup_cost_w = [random.randint(35, 70) for _ in range(W * P)]
    prod_cost = [random.randint(12, 28) for _ in range(F * P)]
    transport_fw = [random.randint(15, 40) for _ in range(F * W * P)]
    transport_wr = [random.randint(10, 30) for _ in range(W * R * P)]
    holding_f = [random.randint(4, 10) for _ in range(F * P)]
    holding_w = [random.randint(4, 10) for _ in range(W * P)]

    batch_size = [random.choice([7, 11, 13]) for _ in range(P)]
    truck_size = 19

    demand = []
    for r in range(R):
        for p in range(P):
            for t in range(T):
                base = random.randint(3, 7)
                # MUST be multiple of truck_size for integer feasibility
                d = base * truck_size
                demand.append(d)

    max_demand_per_period = max(
        sum(demand[r * P * T + p * T + t]
            for r in range(R) for p in range(P))
        for t in range(T)
    )

    total_cap_needed = int(max_demand_per_period * 1.12)

    cap_factory = []
    cap_factory.append(int(total_cap_needed * 0.35))
    cap_factory.append(int(total_cap_needed * 0.30))
    remaining = total_cap_needed - sum(cap_factory)
    for f in range(2, F):
        if f < F - 1:
            cap = int(remaining * 0.40)
        else:
            cap = remaining - sum(cap_factory[2:])
        cap_factory.append(max(max(batch_size) * 2, cap))

    cap_warehouse = [int(cap_factory[w % F] * random.uniform(0.45, 0.75))
                     for w in range(W)]

    fw_arcs = []
    for f in range(F):
        n_connections = max(2, int(W * 0.60))
        connected_warehouses = random.sample(range(W), n_connections)
        for w in connected_warehouses:
            fw_arcs.append((f, w))

    wr_arcs = []
    for w in range(W):
        n_connections = max(2, int(R * 0.60))
        connected_retailers = random.sample(range(R), n_connections)
        for r in connected_retailers:
            wr_arcs.append((w, r))

    # Ensure connectivity
    reachable_retailers = set()
    for w, r in wr_arcs:
        for f, ww in fw_arcs:
            if ww == w:
                reachable_retailers.add(r)
                break

    for r in range(R):
        if r not in reachable_retailers:
            for w in range(W):
                for f, ww in fw_arcs:
                    if ww == w:
                        wr_arcs.append((w, r))
                        reachable_retailers.add(r)
                        break
                if r in reachable_retailers:
                    break

    return MESCNInstance(
        name=name, F=F, W=W, R=R, P=P, T=T,
        setup_cost_f=setup_cost_f, setup_cost_w=setup_cost_w, prod_cost=prod_cost,
        transport_fw=transport_fw, transport_wr=transport_wr,
        holding_f=holding_f, holding_w=holding_w,
        cap_factory=cap_factory, cap_warehouse=cap_warehouse,
        batch_size=batch_size, truck_size=truck_size, demand=demand,
        fully_connected=0, fw_arcs=fw_arcs, wr_arcs=wr_arcs
    )
