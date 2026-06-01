# MESCN Test Suite

## Purpose

The Multi-Echelon Supply Chain Network (MESCN) test suite provides comprehensive MIP correctness tests exercising all variable types and constraint combinations.

## Test Instances

| Instance | Size | Type | Purpose |
|---|---|---|---|
| `tiny_feasible` | 1F×1W×1R×1P×2T | Optimal (known=440) | Minimal hand-crafted instance for quick validation |
| `small_1234` | 2F×2W×3R×2P×3T | Feasible | Random seed 1234 |
| `small_5678` | 2F×2W×3R×2P×3T | Feasible | Random seed 5678 |
| `small_9012` | 2F×2W×3R×2P×3T | Feasible | Random seed 9012 |
| `infeas_demand_spike` | 2F×2W×2R×1P×3T | Infeasible | Total demand exceeds capacity |
| `infeas_disconnected` | 2F×2W×2R×1P×2T | Infeasible | Network partition (no path to demand) |

## Variable Types Tested

- **Binary (u, v)**: Factory and warehouse setup decisions
- **General integer (n, m)**: Production batch counts, truck counts
- **Continuous (x, y, z, s_f, s_w)**: Production quantities, shipments, inventory levels

## Constraint Types Tested

1. **Flow balance**: Multi-echelon inventory balance (factory → warehouse → retail)
2. **Capacity**: Factory production capacity, warehouse storage capacity
3. **Batch production**: `x = batch_size * n` (continuous-integer coupling)
4. **Setup links (Big-M)**: `n ≤ M*u`, `m ≤ M*v` (binary-integer coupling)
5. **Truck shipments**: `Σz = truck_size * m` (continuous-integer aggregation)
6. **Demand satisfaction**: `Σz[w,r,p,t] = demand[r,p,t]` (equality constraints)

## Why These Instances Are "Easy"

All feasible instances solve at the root node (0 branch-and-bound nodes). This is **expected and desirable** for CI testing:

### LP Relaxation Tightness

Supply chain formulations with continuous flow variables naturally have tight LP relaxations. The continuous variables (production, shipments, inventory) are highly constrained by:
- Flow balance equations (equality constraints)
- Batch and truck size constraints
- Capacity limits

This leaves little "room" for fractional solutions - the LP relaxation often yields integer values for the setup/batch/truck variables.

### This Is Good for CI

CI tests need to:
- **Run quickly** (<5s per test) ✓
- **Exercise all features** (cuts, heuristics, preprocessing) ✓
- **Catch correctness bugs** (wrong optimal, wrong infeasible, invalid cuts) ✓
- **Validate solution feasibility** (all constraints satisfied) ✓

Hard instances (requiring hours and thousands of nodes) are:
- **Too slow** for CI (would block every commit)
- **Unnecessary** for correctness (bugs manifest on small instances too)
- **Better suited** for performance benchmarking (separate from CI)

## Solution Validation

`CInterfaceTest_mescn.c` validates every solution against:
- Variable integrality (binary and integer variables)
- All 10 constraint types
- Objective value consistency

On failure, runs `mip_diag_wrong_optimal()` to identify the responsible feature.

## Generators

| Generator | Purpose |
|---|---|
| `gen_mescn_fixtures_simple.py` | Simple feasible/infeasible instances for CI (current fixtures) |
| `gen_mescn_fixtures_moderate.py` | Medium-sized instances with tighter constraints (still solve at root) |
| `gen_mescn_fixtures_challenging.py` | High setup costs, demand spikes (still solve at root) |

All generators produce tight LP relaxations because the problem structure (continuous flows + equality constraints) naturally leads to integer solutions at the root.

## Future Extensions

If **algorithmic stress tests** (many nodes, long solves) are needed:

1. **Different problem class**: Consider problems with naturally loose LP relaxations:
   - Pure integer programs (no continuous variables)
   - Combinatorial problems (TSP, graph coloring, scheduling)
   - Set packing/covering (sparse constraint matrices)

2. **Adversarial instances**: Hand-crafted to defeat specific heuristics/cuts

3. **Real benchmarks**: Use MIPLIB instances known to be hard

But for **CI correctness testing**, the current MESCN suite is ideal.

## Cross-Validation

Optimal values for feasible instances are certified by solving with MIPster and validating with `mipster_validate_sol` to ensure:
- All variable bounds satisfied
- Integrality constraints met
- All row constraints satisfied
- Objective value computed correctly

Reference solutions are embedded in `mescn_fixtures.h` for use with `OsiRowCutDebugger` during cut validation.

## Running the Tests

```bash
# Build the test binary
make CInterfaceTest_mescn

# Run the test
./CInterfaceTest_mescn

# Analyze instance difficulty
python3 test_mescn_instances.py ../src/mipster
```

## Integration

Test is integrated into:
- `test/Makefile.am` (build system)
- `test/run-mipster-tests` (CI runner)
- GitHub Actions CI pipeline
