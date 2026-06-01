#!/usr/bin/env python3
"""Generate Multi-Echelon Supply Chain Network (MESCN) test fixtures.

Simplified version that doesn't require Gurobi/HiGHS for cross-validation.
Instead, generates small instances with known solutions that can be manually verified.
"""

import argparse
import sys
from dataclasses import dataclass
from typing import List, Optional, Tuple

import numpy as np


@dataclass
class MESCNInstance:
    """Multi-Echelon Supply Chain Network instance."""
    name: str
    F: int  # num factories
    W: int  # num warehouses
    R: int  # num retail stores
    P: int  # num products
    T: int  # num time periods

    # Cost parameters
    setup_cost_f: List[int]
    setup_cost_w: List[int]
    prod_cost: List[int]
    transport_fw: List[int]
    transport_wr: List[int]
    holding_f: List[int]
    holding_w: List[int]

    # Capacity and batch parameters
    cap_factory: List[int]
    cap_warehouse: List[int]
    batch_size: List[int]
    truck_size: int

    # Demand
    demand: List[int]

    # Network structure
    fw_arcs: Optional[List[Tuple[int, int]]]
    wr_arcs: Optional[List[Tuple[int, int]]]

    # Certified solution (will be computed or manually specified)
    optimal_obj: Optional[int]


def generate_tiny_feasible():
    """Generate a tiny feasible instance with obvious optimal solution.

    Structure:
    - 1 factory, 1 warehouse, 1 retailer, 1 product, 2 periods
    - Simple demand pattern: [10, 20]
    - Batch size 10, truck size 10
    - Optimal: produce 1 batch period 1, 2 batches period 2
    """
    F, W, R, P, T = 1, 1, 1, 1, 2

    # Costs (favor just-in-time production)
    setup_cost_f = [100]  # F*P=1
    setup_cost_w = [50]   # W*P=1
    prod_cost = [5]       # F*P=1
    transport_fw = [2]    # F*W*P=1
    transport_wr = [1]    # W*R*P=1
    holding_f = [10]      # F*P=1 (expensive to hold at factory)
    holding_w = [10]      # W*P=1 (expensive to hold at warehouse)

    # Capacities
    cap_factory = [500]   # F=1 (plenty of capacity)
    cap_warehouse = [500] # W=1
    batch_size = [10]     # P=1
    truck_size = 10

    # Demand
    demand = [10, 20]     # R*P*T = 1*1*2

    # Network: fully connected
    fw_arcs = None
    wr_arcs = None

    # Optimal solution (manually computed):
    # Period 0: produce 10 (1 batch), ship 10, deliver 10 -> obj = 100 + 5*10 + 2*10 + 1*10 = 180
    # Period 1: produce 20 (2 batches), ship 20, deliver 20 -> obj = 100 + 5*20 + 2*20 + 1*20 = 260
    # Total: 440
    optimal_obj = 440

    return MESCNInstance(
        name="tiny_feasible",
        F=F, W=W, R=R, P=P, T=T,
        setup_cost_f=setup_cost_f,
        setup_cost_w=setup_cost_w,
        prod_cost=prod_cost,
        transport_fw=transport_fw,
        transport_wr=transport_wr,
        holding_f=holding_f,
        holding_w=holding_w,
        cap_factory=cap_factory,
        cap_warehouse=cap_warehouse,
        batch_size=batch_size,
        truck_size=truck_size,
        demand=demand,
        fw_arcs=fw_arcs,
        wr_arcs=wr_arcs,
        optimal_obj=optimal_obj
    )


def generate_small_instance(seed: int):
    """Generate a small random feasible instance."""
    rng = np.random.RandomState(seed)

    F, W, R, P, T = 2, 2, 3, 2, 3

    # Costs
    setup_cost_f = [rng.randint(80, 150) for _ in range(F * P)]
    setup_cost_w = [rng.randint(40, 80) for _ in range(W * P)]
    prod_cost = [rng.randint(3, 10) for _ in range(F * P)]
    transport_fw = [rng.randint(1, 5) for _ in range(F * W * P)]
    transport_wr = [rng.randint(1, 3) for _ in range(W * R * P)]
    holding_f = [rng.randint(1, 3) for _ in range(F * P)]
    holding_w = [rng.randint(1, 3) for _ in range(W * P)]

    # Batch and truck
    batch_size = [rng.randint(20, 40) for _ in range(P)]
    truck_size = 100

    # Demand (designed to be satisfiable)
    base_demand = [rng.randint(30, 70) for _ in range(R * P * T)]

    # Capacities (generous to ensure feasibility)
    total_demand = sum(base_demand)
    avg_per_period = total_demand / T
    cap_factory = [int(avg_per_period / F * 1.5) for _ in range(F)]
    cap_warehouse = [int(avg_per_period / W * 2.0) for _ in range(W)]

    return MESCNInstance(
        name=f"small_{seed}",
        F=F, W=W, R=R, P=P, T=T,
        setup_cost_f=setup_cost_f,
        setup_cost_w=setup_cost_w,
        prod_cost=prod_cost,
        transport_fw=transport_fw,
        transport_wr=transport_wr,
        holding_f=holding_f,
        holding_w=holding_w,
        cap_factory=cap_factory,
        cap_warehouse=cap_warehouse,
        batch_size=batch_size,
        truck_size=truck_size,
        demand=base_demand,
        fw_arcs=None,
        wr_arcs=None,
        optimal_obj=None  # Will be solved
    )


def generate_infeasible_demand_spike():
    """Generate instance with impossible demand spike."""
    F, W, R, P, T = 2, 2, 2, 1, 3

    setup_cost_f = [100, 100]
    setup_cost_w = [50, 50]
    prod_cost = [5, 5]
    transport_fw = [2, 2, 2, 2]
    transport_wr = [1, 1, 1, 1]
    holding_f = [2, 2]
    holding_w = [2, 2]

    cap_factory = [100, 100]  # Total capacity per period: 200
    cap_warehouse = [200, 200]
    batch_size = [10]
    truck_size = 10

    # Demand spike in middle period exceeds total capacity
    # Period 0: 20, Period 1: 300 (IMPOSSIBLE), Period 2: 20
    demand = [
        10, 150, 10,  # retailer 0, product 0, periods 0-2
        10, 150, 10   # retailer 1, product 0, periods 0-2
    ]  # Period 1 total: 300 > 200 (max capacity)

    return MESCNInstance(
        name="infeas_demand_spike",
        F=F, W=W, R=R, P=P, T=T,
        setup_cost_f=setup_cost_f,
        setup_cost_w=setup_cost_w,
        prod_cost=prod_cost,
        transport_fw=transport_fw,
        transport_wr=transport_wr,
        holding_f=holding_f,
        holding_w=holding_w,
        cap_factory=cap_factory,
        cap_warehouse=cap_warehouse,
        batch_size=batch_size,
        truck_size=truck_size,
        demand=demand,
        fw_arcs=None,
        wr_arcs=None,
        optimal_obj=-1  # Infeasible
    )


def generate_infeasible_disconnected():
    """Generate instance with disconnected network."""
    F, W, R, P, T = 2, 2, 2, 1, 2

    setup_cost_f = [100, 100]
    setup_cost_w = [50, 50]
    prod_cost = [5, 5]
    transport_fw = [2, 2, 2, 2]
    transport_wr = [1, 1, 1, 1]
    holding_f = [2, 2]
    holding_w = [2, 2]

    cap_factory = [200, 200]
    cap_warehouse = [200, 200]
    batch_size = [10]
    truck_size = 10

    # Normal demand
    demand = [30, 30, 30, 30]  # 2 retailers * 1 product * 2 periods

    # Disconnected network: factory 0 can only reach warehouse 0,
    # warehouse 0 can only reach retailer 0
    # But retailer 1 has demand with no supply path
    fw_arcs = [(0, 0), (1, 1)]  # Factory 0->Warehouse 0, Factory 1->Warehouse 1
    wr_arcs = [(0, 0), (1, 1)]  # Warehouse 0->Retailer 0, Warehouse 1->Retailer 1
    # This is actually feasible! Let me make it truly disconnected:
    # Remove connection to retailer 1
    wr_arcs = [(0, 0), (1, 0)]  # Both warehouses can only reach retailer 0
    # Now retailer 1 cannot be supplied

    return MESCNInstance(
        name="infeas_disconnected",
        F=F, W=W, R=R, P=P, T=T,
        setup_cost_f=setup_cost_f,
        setup_cost_w=setup_cost_w,
        prod_cost=prod_cost,
        transport_fw=transport_fw,
        transport_wr=transport_wr,
        holding_f=holding_f,
        holding_w=holding_w,
        cap_factory=cap_factory,
        cap_warehouse=cap_warehouse,
        batch_size=batch_size,
        truck_size=truck_size,
        demand=demand,
        fw_arcs=fw_arcs,
        wr_arcs=wr_arcs,
        optimal_obj=-1  # Infeasible
    )


def write_c_header(instances: List[MESCNInstance], output_path: str):
    """Write instances to C header file."""
    lines = []
    lines.append("/* Multi-Echelon Supply Chain Network (MESCN) test fixtures.")
    lines.append(" * Generated by gen_mescn_fixtures_simple.py — do not edit manually.")
    lines.append(" */")
    lines.append("")
    lines.append("#ifndef MESCN_FIXTURES_H")
    lines.append("#define MESCN_FIXTURES_H")
    lines.append("")
    lines.append("#include <stddef.h>")
    lines.append("")

    lines.append("typedef struct {")
    lines.append("  const char *name;")
    lines.append("  int certified_optimal;  /* -1 for infeasible, -999 for unknown, >= 0 for optimal */")
    lines.append("  int F, W, R, P, T;  /* dimensions */")
    lines.append("  const int *setup_cost_f;      /* [F*P] */")
    lines.append("  const int *setup_cost_w;      /* [W*P] */")
    lines.append("  const int *prod_cost;         /* [F*P] */")
    lines.append("  const int *transport_fw;      /* [F*W*P] */")
    lines.append("  const int *transport_wr;      /* [W*R*P] */")
    lines.append("  const int *holding_f;         /* [F*P] */")
    lines.append("  const int *holding_w;         /* [W*P] */")
    lines.append("  const int *cap_factory;       /* [F] */")
    lines.append("  const int *cap_warehouse;     /* [W] */")
    lines.append("  const int *batch_size;        /* [P] */")
    lines.append("  int truck_size;")
    lines.append("  const int *demand;            /* [R*P*T] */")
    lines.append("  int fully_connected;          /* 1 if all arcs present, 0 if custom network */")
    lines.append("  int n_fw_arcs;")
    lines.append("  const int *fw_arcs;           /* [(f,w) pairs] flattened */")
    lines.append("  int n_wr_arcs;")
    lines.append("  const int *wr_arcs;           /* [(w,r) pairs] flattened */")
    lines.append("} MESCNInstance;")
    lines.append("")

    # Data arrays
    for inst in instances:
        prefix = inst.name.replace('-', '_')

        def write_array(name, arr):
            lines.append(f"static const int {prefix}_{name}[] = {{")
            lines.append("  " + ", ".join(map(str, arr)))
            lines.append("};")
            lines.append("")

        write_array("setup_cost_f", inst.setup_cost_f)
        write_array("setup_cost_w", inst.setup_cost_w)
        write_array("prod_cost", inst.prod_cost)
        write_array("transport_fw", inst.transport_fw)
        write_array("transport_wr", inst.transport_wr)
        write_array("holding_f", inst.holding_f)
        write_array("holding_w", inst.holding_w)
        write_array("cap_factory", inst.cap_factory)
        write_array("cap_warehouse", inst.cap_warehouse)
        write_array("batch_size", inst.batch_size)
        write_array("demand", inst.demand)

        # Network arcs
        if inst.fw_arcs is not None:
            lines.append(f"static const int {prefix}_fw_arcs[] = {{")
            flat = []
            for f, w in inst.fw_arcs:
                flat.extend([f, w])
            lines.append("  " + ", ".join(map(str, flat)))
            lines.append("};")
            lines.append("")

        if inst.wr_arcs is not None:
            lines.append(f"static const int {prefix}_wr_arcs[] = {{")
            flat = []
            for w, r in inst.wr_arcs:
                flat.extend([w, r])
            lines.append("  " + ", ".join(map(str, flat)))
            lines.append("};")
            lines.append("")

    # Instance array
    lines.append("static const MESCNInstance MESCN_FIXTURES[] = {")
    for inst in instances:
        prefix = inst.name.replace('-', '_')
        fully_connected = 1 if inst.fw_arcs is None else 0
        n_fw = 0 if inst.fw_arcs is None else len(inst.fw_arcs)
        n_wr = 0 if inst.wr_arcs is None else len(inst.wr_arcs)
        fw_ptr = "NULL" if inst.fw_arcs is None else f"{prefix}_fw_arcs"
        wr_ptr = "NULL" if inst.wr_arcs is None else f"{prefix}_wr_arcs"
        opt_val = inst.optimal_obj if inst.optimal_obj is not None else -999

        lines.append("  {")
        lines.append(f'    .name = "{inst.name}",')
        lines.append(f"    .certified_optimal = {opt_val},")
        lines.append(f"    .F = {inst.F}, .W = {inst.W}, .R = {inst.R}, .P = {inst.P}, .T = {inst.T},")
        lines.append(f"    .setup_cost_f = {prefix}_setup_cost_f,")
        lines.append(f"    .setup_cost_w = {prefix}_setup_cost_w,")
        lines.append(f"    .prod_cost = {prefix}_prod_cost,")
        lines.append(f"    .transport_fw = {prefix}_transport_fw,")
        lines.append(f"    .transport_wr = {prefix}_transport_wr,")
        lines.append(f"    .holding_f = {prefix}_holding_f,")
        lines.append(f"    .holding_w = {prefix}_holding_w,")
        lines.append(f"    .cap_factory = {prefix}_cap_factory,")
        lines.append(f"    .cap_warehouse = {prefix}_cap_warehouse,")
        lines.append(f"    .batch_size = {prefix}_batch_size,")
        lines.append(f"    .truck_size = {inst.truck_size},")
        lines.append(f"    .demand = {prefix}_demand,")
        lines.append(f"    .fully_connected = {fully_connected},")
        lines.append(f"    .n_fw_arcs = {n_fw},")
        lines.append(f"    .fw_arcs = {fw_ptr},")
        lines.append(f"    .n_wr_arcs = {n_wr},")
        lines.append(f"    .wr_arcs = {wr_ptr}")
        lines.append("  },")
    lines.append("};")
    lines.append("")
    lines.append(f"#define N_MESCN_FIXTURES {len(instances)}")
    lines.append("")
    lines.append("#endif /* MESCN_FIXTURES_H */")

    with open(output_path, 'w') as f:
        f.write('\n'.join(lines))

    print(f"Wrote {output_path}")


def main():
    parser = argparse.ArgumentParser(description="Generate MESCN test fixtures (simple version)")
    parser.add_argument('--output', type=str, default='mescn_fixtures.h',
                       help='Output C header file')
    args = parser.parse_args()

    instances = []

    # Tiny instance with manually computed optimal
    print("Generating tiny_feasible...")
    instances.append(generate_tiny_feasible())

    # Small random instances (optimal will be computed by test)
    for seed in [1234, 5678, 9012]:
        print(f"Generating small_{seed}...")
        instances.append(generate_small_instance(seed))

    # Infeasible instances
    print("Generating infeas_demand_spike...")
    instances.append(generate_infeasible_demand_spike())

    print("Generating infeas_disconnected...")
    instances.append(generate_infeasible_disconnected())

    # Write header
    write_c_header(instances, args.output)

    print(f"\nGenerated {len(instances)} fixtures:")
    for inst in instances:
        if inst.optimal_obj == -1:
            status = "INFEASIBLE"
        elif inst.optimal_obj is None or inst.optimal_obj == -999:
            status = "optimal=UNKNOWN (will be computed)"
        else:
            status = f"optimal={inst.optimal_obj}"
        print(f"  {inst.name:30s}  {status}")


if __name__ == '__main__':
    main()
