# VRP Test Suite

Comprehensive test suite for CVRP (Capacitated Vehicle Routing Problem) and VRPPD (Vehicle Routing Problem with Pickup and Delivery).

## Quick Start

```bash
# Generate all test instances
python3 gen_vrp_fixtures.py

# Test a simple instance
mipster cvrp_loose.mps -solve

# Test with time limit
mipster cvrp_size_large.mps -sec 300 -solve
```

## What's Included

- **21 test instances** (16 CVRP + 5 VRPPD)
- **Diverse characteristics**: size, density, demands, capacity, geography
- **Range of difficulty**: from instant (< 0.01s) to challenging (> 60s)
- **All use correct formulation** with validation tools available

## Instance Categories

### CVRP (Capacitated VRP)
- **Original**: Simple baseline instances (loose, tight, medium)
- **Density**: Sparse grid vs dense uniform
- **Demands**: Uniform, skewed, outliers
- **Capacity**: Tight (forces many vehicles) to loose (single vehicle)
- **Geography**: Clustered, ring, uniform
- **Size**: Small (6) to large (20 customers)

### VRPPD (VRP with Pickup-Delivery)
- **Basic**: Simple 2-3 request instances
- **Diverse**: Various layouts (uniform, clustered, ring)
- **Challenging**: Up to 6 requests, tight capacity

## Key Technical Details

### Correct Capacity Formulation

After extensive debugging, discovered the correct formulation:

```
✓ Correct:   u[i] - u[j] + Q*x[i,j] <= Q - d[j]
✗ Failed:    u[j] - u[i] + M*x[i,j] >= d[j] + M
```

**Key differences:**
1. Inequality direction: `<=` not `>=`
2. Variable order: `u[i] - u[j]` not `u[j] - u[i]`
3. Use capacity `Q` directly
4. RHS: `Q - d[j]` not `d[j] + M`

### VRPPD Specifics

Additional requirements:
1. **Precedence constraints**: Forbid delivery→pickup arcs
   ```
   x[delivery_r, pickup_r] <= 0
   ```
2. **Non-monotonic demands**: Pickup positive, delivery negative
3. **Load bounds**: `d[i] <= u[i] <= Q`

## Files

| File | Purpose |
|------|---------|
| `gen_vrp_fixtures.py` | Comprehensive fixture generator |
| `gen_cvrp_working.py` | Original simple generator |
| `VRP_INSTANCE_CATALOG.md` | Complete instance catalog with details |
| `VRP_TEST_FIXTURES_SUMMARY.md` | Overview of characteristics covered |
| `CVRP_VRPPD_BREAKTHROUGH.md` | Technical details of formulation discovery |
| `VRP_README.md` | This file |

## Sample Results

| Instance | Type | Size | Time | Nodes | Objective |
|----------|------|------|------|-------|-----------|
| cvrp_loose | CVRP | 4 | 0.03s | 0 | 90 |
| cvrp_tight | CVRP | 4 | 0.01s | 0 | 130 |
| cvrp_medium | CVRP | 10 | 20s | 259K | 320 |
| cvrp_size_small | CVRP | 6 | 0.02s | 24 | 20158 |
| cvrp_size_large | CVRP | 20 | >45s | >200K | ? |
| vrppd_small | VRPPD | 2 req | 0.01s | 0 | 70 |
| vrppd_tight_ring | VRPPD | 4 req | 0.03s | 0 | 397 |
| vrppd_large_uniform | VRPPD | 6 req | 0.13s | 12 | 370 |

## Testing Strategy

### Quick Regression (< 1s)
```bash
mipster cvrp_loose.mps -solve
mipster cvrp_size_small.mps -solve
mipster vrppd_small.mps -solve
```

### Medium Tests (1-30s)
```bash
mipster cvrp_demand_skewed.mps -sec 30 -solve
mipster cvrp_capacity_medium.mps -sec 30 -solve
mipster vrppd_large_uniform.mps -sec 30 -solve
```

### Challenging Benchmarks (> 30s)
```bash
mipster cvrp_medium.mps -sec 300 -solve
mipster cvrp_size_large.mps -sec 300 -solve
```

## Extending the Test Suite

The generator supports:

```python
# Custom CVRP
generate_cvrp(
    name="my_instance",
    K=3,                    # vehicles
    customers=[1,2,3,4,5],  # customer IDs
    Q=25,                   # capacity
    demands={1:5, 2:3, ...},
    costs={(i,j): dist, ...},
    desc="Description"
)

# Custom VRPPD
generate_vrppd(
    name="my_vrppd",
    K=2,                    # vehicles
    requests=4,             # pickup-delivery pairs
    Q=20,                   # capacity
    pickup_demand=6,
    costs={(i,j): dist, ...},
    desc="Description"
)
```

Utility functions:
- `generate_positions(n, pattern, seed)` - Create node layouts
- `compute_cost_matrix(positions, cost_type)` - Distance matrix
- `generate_demands(n, pattern, total, seed)` - Demand patterns

## History

- Initial attempts used incorrect capacity formulation (all infeasible)
- Breakthrough: discovered correct formulation from online reference
- Generated diverse test suite covering wide range of characteristics
- All instances can be validated with `validate_vrp_solutions.sh` (feasibility) and `cross_validate_solvers.py` (cross-solver comparison)

## References

See `VRP_TESTS.md` for complete technical details including:
- Problem formulations (CVRP and VRPPD)
- Complete instance catalog
- Why the capacity constraint formulation is critical
- Validation procedures
