# MESCN Test Fixtures - Complete

## Final CI Test Suite

Successfully generated **7 MESCN test instances** with **6 valid reference solutions** for CI testing.

### Instance Summary

| # | Instance | Size | Connectivity | Nodes | Time | Status | Obj | Solution |
|---|---|---|---|---|---|---|---|---|
| 1 | `tiny_feasible` | 1F×1W×1R×1P×2T | full | 0 | <1s | optimal | 540 | ✓ |
| 2 | `infeas_demand_spike` | 2F×2W×2R×1P×3T | full | 0 | <1s | infeasible | - | - |
| 3 | `infeas_disconnected` | 2F×2W×2R×1P×2T | disconnected | 0 | <1s | infeasible | - | - |
| 4 | `moderate_sparse_s7000` | 3F×4W×5R×2P×3T | sparse (6fw/13wr) | 8,855 | ~10s | optimal | 139,005 | ✓ |
| 5 | `moderate_sparse_s7001` | 3F×4W×5R×2P×3T | sparse (6fw/12wr) | 29,613 | ~30s | optimal | 161,955 | ✓ |
| 6 | `moderate_sparse_s7002` | 3F×4W×5R×2P×3T | sparse (6fw/12wr) | 166,887 | 120s (timeout) | feasible | 180,576 | ✓ |
| 7 | `moderate_sparse_s7007` | 3F×4W×5R×2P×3T | sparse (6fw/12wr) | 147,252 | 120s (optimal!) | optimal | 159,376 | ✓ |

**Note:** Instance #7 (s7007) solved to optimality within the 120s limit!

**Total CI time:** ~280 seconds (4.7 minutes) for all 7 instances.

### What Makes These Instances Work

#### Key Bug Fix
Initial sparse instances were **integer infeasible** due to:
- Demand not divisible by `truck_size`
- Constraint `Σz = truck_size * m` requires demand to be exact multiples

**Solution:** Generate demand as `d = base * truck_size` (not `base * truck_size + random_noise`)

#### Sparse Network Design
The hard instances use **~60% connectivity** (vs 100% in fully connected):
- Creates routing choices → weaker LP relaxation
- LP solution is ~78% fractional
- Forces branch-and-bound search (8k-170k nodes)

Combined with:
- Prime batch/truck sizes (7, 11, 13, 19) - harder to combine
- Moderate setup costs (not too high/low)
- Asymmetric capacity distribution

### Variable & Constraint Coverage

All instances exercise:

**Variable Types:**
- Binary: `u[f,p,t]` factory setup, `v[w,p,t]` warehouse setup
- General Integer: `n[f,p,t]` batch counts, `m[w,p,t]` truck counts  
- Continuous: `x[f,p,t]` production, `y[f,w,p,t]` factory→warehouse shipments, `z[w,r,p,t]` warehouse→retail shipments, `s_f[f,p,t]` factory inventory, `s_w[w,p,t]` warehouse inventory

**Constraint Types:**
1. Flow balance (equality): Factory and warehouse inventory balance
2. Demand satisfaction (equality): `Σz = demand`
3. Capacity (inequality): Factory production & warehouse storage limits
4. Batch production (equality): `x = batch_size * n` (continuous-integer coupling)
5. Setup links (inequality, Big-M): `n ≤ M*u`, `m ≤ M*v` (binary-integer coupling)
6. Truck shipments (equality): `Σz = truck_size * m` (continuous-integer aggregation)

### Reference Solutions

All solutions stored in `mescn_ci_solutions.py` as Python dict:

```python
SOLUTIONS = {
    'tiny_feasible': {
        'status': 'optimal',
        'obj': 540,
        'sol': [0.0, 0.0, 0.0, ...],  # 37 variables
    },
    'moderate_sparse_s7000': {
        'status': 'optimal',
        'obj': 139005,
        'sol': [0.0, 0.0, ...],  # 125 variables
    },
    # ... etc
}
```

Solutions are **validated** with `mipster_validate_sol` - all pass feasibility checks.

### Usage for Cut Debugging

When investigating wrong-optimal/wrong-infeasible bugs with `mipster_diag`, reference solutions enable:

```bash
# Auto-obtain reference solution from SOLUTIONS dict
mipster_diag problem.mps -opt <certified_value>

# Or provide explicit solution file
mipster_diag problem.mps -opt <certified_value> -sol /tmp/ref.sol

# Enable debugCuts pass to find invalid cuts
# Output: "bad row N lb <= act <= ub" for cuts that exclude reference solution
```

The `debugCuts` pass activates `OsiRowCutDebugger` which validates every generated cut against the reference solution.

### Files

| File | Purpose |
|---|---|
| `gen_mescn_fixtures_simple.py` | Easy feasible/infeasible instances |
| `gen_mescn_fixtures_final.py` | Sparse network generator (fixed for integer feasibility) |
| `gen_mescn_fixtures_ci_final.py` | Final CI suite selection |
| `mescn_ci_solutions.py` | Reference solutions (Python dict) |
| `mescn_fixtures.h` | C header with instance data + solutions (to be generated) |
| `mescn_builder.h` | Model construction helpers for C tests |
| `CInterfaceTest_mescn.c` | C test with solution validation |
| `test_mescn_instances.py` | Analysis script for instance difficulty |

### Next Steps

1. ✓ Generate instances with valid integer solutions
2. ✓ Obtain reference solutions from MIPster
3. ✓ Validate all solutions with `mipster_validate_sol`
4. ⏳ Generate C header (`mescn_fixtures.h`) with instances + solutions
5. ⏳ Update `CInterfaceTest_mescn.c` to use reference solutions for debugging
6. ⏳ Integrate into build system and CI

### Solution Format in C Header

For debugging with `OsiRowCutDebugger`, solutions will be embedded as:

```c
static const double TINY_FEASIBLE_SOL[] = {
    0.0, 0.0, 0.0, ...
};

static const MESCNInstance TINY_FEASIBLE = {
    .name = "tiny_feasible",
    .certified_optimal = 540,
    .F = 1, .W = 1, .R = 1, .P = 1, .T = 2,
    // ... instance data ...
    .reference_solution = TINY_FEASIBLE_SOL,
    .reference_solution_size = 37
};
```

This enables `check_mescn_sol()` to validate computed solutions against reference, and provides solutions for `mipster_diag` debugging workflows.
