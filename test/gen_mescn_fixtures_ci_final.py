#!/usr/bin/env python3
"""Final MESCN CI test suite with verified integer-feasible instances."""

from gen_mescn_fixtures_simple import (
    generate_tiny_feasible,
    generate_infeasible_demand_spike,
    generate_infeasible_disconnected
)
from gen_mescn_fixtures_final import generate_moderate_sparse

# Final CI test suite: 7 instances
CI_FIXTURES = [
    # Easy: correctness test with known optimal
    generate_tiny_feasible(),           # 1F×1W×1R×1P×2T, optimal=540, 0 nodes, <1s

    # Infeasible tests
    generate_infeasible_demand_spike(), # 2F×2W×2R×1P×3T, infeasible
    generate_infeasible_disconnected(), # 2F×2W×2R×1P×2T, infeasible

    # Medium: requires branching, solves optimally
    generate_moderate_sparse(7000),     # 3F×4W×5R×2P×3T, ~9k nodes, ~10s
    generate_moderate_sparse(7001),     # 3F×4W×5R×2P×3T, ~30k nodes, ~30s

    # Hard: requires extensive branching, times out but finds solutions
    generate_moderate_sparse(7002),     # 3F×4W×5R×2P×3T, ~167k nodes, 60s timeout
    generate_moderate_sparse(7007),     # 3F×4W×5R×2P×3T, ~147k nodes, 60s timeout
]

if __name__ == '__main__':
    print("=" * 80)
    print("MESCN CI Test Suite - Final")
    print("=" * 80)
    print()

    for i, inst in enumerate(CI_FIXTURES, 1):
        status_hint = "optimal=540" if inst.name == "tiny_feasible" else (
            "infeasible" if "infeas" in inst.name else
            "optimal (~9k nodes)" if inst.name == "moderate_sparse_s7000" else
            "optimal (~30k nodes)" if inst.name == "moderate_sparse_s7001" else
            "timeout (~167k nodes)" if inst.name == "moderate_sparse_s7002" else
            "timeout (~147k nodes)")

        conn = "full" if getattr(inst, 'fully_connected', 1) else f"sparse ({len(inst.fw_arcs)}fw/{len(inst.wr_arcs)}wr)"

        print(f"{i}. {inst.name:27s} {inst.F}F×{inst.W}W×{inst.R}R×{inst.P}P×{inst.T}T  {conn:20s} → {status_hint}")

    print()
    print("Summary:")
    print("  • 1 easy instance (known optimal, 0 nodes)")
    print("  • 2 infeasible instances")
    print("  • 2 medium instances (optimal with branching, 9k-30k nodes)")
    print("  • 2 hard instances (timeout, find feasible solutions, 147k-167k nodes)")
    print()
    print("Total CI time estimate: ~150 seconds")
    print()
    print("All instances exercise:")
    print("  ✓ Binary, general integer, and continuous variables")
    print("  ✓ Flow balance constraints (multi-echelon network)")
    print("  ✓ Capacity constraints (factory + warehouse)")
    print("  ✓ Big-M setup links (binary-integer coupling)")
    print("  ✓ Batch production & truck shipment constraints")
    print("  ✓ Cut generation & primal heuristics")
