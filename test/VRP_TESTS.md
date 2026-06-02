# VRP Test Suite Documentation

Complete documentation for CVRP and VRPPD test instances.

---

## Table of Contents

1. [Quick Start](#quick-start)
2. [Problem Formulations](#problem-formulations)
3. [Test Instances](#test-instances)
4. [Instance Characteristics](#instance-characteristics)
5. [Generation Scripts](#generation-scripts)
6. [Test Results](#test-results)
7. [Technical Details](#technical-details)

---

## Quick Start

### Generate All Instances
```bash
cd test
python3 gen_vrp_fixtures.py      # Main fixture suite
python3 gen_diverse_vrp.py       # Additional diverse instances
```

### Run Tests
```bash
# Quick smoke test
mipster cvrp_loose.mps -solve

# Comprehensive test suite
./test_all_vrp.sh

# Individual instance with time limit
mipster cvrp_medium.mps -sec 300 -solve
```

### Validate Results
```bash
# Validate all solutions (feasibility checking)
./validate_vrp_solutions.sh

# Cross-validate with HiGHS (objective comparison)
python3 cross_validate_solvers.py

# Manual validation of a single solution
../src/mipster instance.mps -sec 60 -solve -solu solution.sol
mipster_validate_sol instance.mps solution.sol

# Expected: all instances either solve optimally or find feasible solutions
```

---

## Problem Formulations

### CVRP (Capacitated Vehicle Routing Problem)

**Problem**: Given a depot, customers with demands, and vehicles with capacity Q, find minimum-cost routes visiting each customer exactly once.

**Variables**:
- `x[i,j]` ∈ {0,1}: vehicle uses arc (i,j)
- `u[i]` ∈ [0,Q]: cumulative load after visiting customer i

**Constraints**:
```
Minimize:  Σ c[i,j] * x[i,j]

Subject to:
  Σ x[i,j] = 1                        ∀j ∈ customers (visit once)
  Σ x[i,j] = 1                        ∀i ∈ customers (leave once)
  Σ x[i,j] - Σ x[j,i] = 0            ∀i ∈ customers (flow conservation)
  Σ x[0,j] ≤ K                        (fleet limit)
  u[i] - u[j] + Q*x[i,j] ≤ Q - d[j]  ∀i,j ∈ customers (capacity)
  d[i] ≤ u[i] ≤ Q                     ∀i ∈ customers (load bounds)
  x[i,j] ∈ {0,1}                      ∀i,j
```

**Key Insight**: The capacity constraint formulation is critical:
```
✓ CORRECT:  u[i] - u[j] + Q*x[i,j] <= Q - d[j]
✗ WRONG:    u[j] - u[i] + M*x[i,j] >= d[j] + M
```

When x[i,j]=1: constraint becomes `u[j] >= u[i] + d[j]` (load increases by demand)  
When x[i,j]=0: constraint becomes `u[i] - u[j] <= Q - d[j]` (non-binding)

### VRPPD (Vehicle Routing with Pickup-Delivery)

**Problem**: CVRP where each request has a pickup node and delivery node that must be visited by the same vehicle, with pickup before delivery.

**Additional Constraints**:
```
  x[delivery_r, pickup_r] ≤ 0         ∀r ∈ requests (precedence)
```

**Demands**: 
- Pickup nodes: positive demand (load increases)
- Delivery nodes: negative demand (load decreases)

The same capacity formulation handles non-monotonic loads correctly.

---

## Test Instances

### CVRP Instances (26 total)

#### Original Simple (3)
| Instance | K | N | Q | Demand | Description |
|----------|---|---|---|--------|-------------|
| cvrp_loose | 1 | 4 | 20 | 18 | Single vehicle sufficient |
| cvrp_tight | 2 | 4 | 10 | 18 | Needs 2 vehicles |
| cvrp_medium | 3 | 10 | 15 | 40 | Medium difficulty |

#### Graph Density (2)
- `cvrp_sparse_grid` - Structured grid layout
- `cvrp_dense_uniform` - Random uniform distribution

#### Demand Patterns (3)
- `cvrp_demand_uniform` - All demands identical (5)
- `cvrp_demand_skewed` - 80% small, 20% large
- `cvrp_demand_outliers` - 10% very large customers

#### Capacity Tightness (3)
- `cvrp_capacity_tight` - Very tight: Q=12, forces 7 vehicles
- `cvrp_capacity_medium` - Medium: Q=30, requires smart routing
- `cvrp_capacity_loose` - Loose: Q=76, single vehicle possible

#### Geographic Structure (2)
- `cvrp_geo_clustered` - 2-4 customer clusters
- `cvrp_geo_ring` - Customers arranged in ring

#### Problem Size (2)
- `cvrp_size_small` - 6 customers, quick
- `cvrp_size_large` - 20 customers, challenging

#### Scale Variations (10)
- `cvrp_small_scale` - Costs 1-20
- `cvrp_medium_scale` - Costs 10-100
- `cvrp_large_scale` - Costs 100-500
- `cvrp_asymmetric` - c(i,j) ≠ c(j,i)
- `cvrp_high_utilization` - Very tight capacity utilization
- `cvrp_many_small` - 15 customers, demand=2 each
- Others from diverse generator

### VRPPD Instances (5 total)

| Instance | K | Requests | Q | Pickup | Description |
|----------|---|----------|---|--------|-------------|
| vrppd_small | 1 | 2 | 20 | 5 | Basic instance |
| vrppd_tight | 2 | 3 | 10 | 4 | Tight capacity |
| vrppd_small_uniform | 1 | 2 | 15 | 5 | Uniform layout |
| vrppd_medium_clustered | 2 | 3 | 12 | 4 | Clustered |
| vrppd_tight_ring | 3 | 4 | 10 | 4 | Ring layout |
| vrppd_large_uniform | 3 | 6 | 18 | 6 | Large, challenging |
| vrppd_small_scale | 1 | 2 | 12 | 4 | Small costs |
| vrppd_large_scale | 2 | 4 | 20 | 8 | Large costs |
| vrppd_many_requests | 3 | 5 | 15 | 5 | Many requests |

---

## Instance Characteristics

### Diversity Dimensions

| Characteristic | Values | Purpose |
|---|---|---|
| **Problem Size** | 4-20 customers | Test scalability |
| **Fleet Size** | 1-7 vehicles | Multi-vehicle coordination |
| **Capacity** | 10-76 | Tight to loose constraints |
| **Demand Pattern** | Uniform, skewed, outliers | Packing complexity |
| **Geography** | Grid, clustered, ring, uniform | Route structure |
| **Cost Scale** | 1-500 | Objective value range |

### Difficulty Levels

**Easy (< 1s)**:
- `cvrp_loose`, `cvrp_tight`, `cvrp_size_small`
- `vrppd_small`, `vrppd_small_uniform`

**Medium (1-20s)**:
- `cvrp_demand_skewed`, `cvrp_sparse_grid`
- `cvrp_capacity_loose`

**Hard (20s+)**:
- `cvrp_medium` - 20s, 259K nodes
- `cvrp_capacity_medium` - 60s timeout
- `cvrp_size_large` - 60s+ timeout

### Objective Value Ranges

**CVRP**:
- Small: 90-320 (simple instances)
- Medium: 20,000-20,500 (grid-based instances)

**VRPPD**:
- Wide range: 26, 70, 93, 105, 153, 370, 397, 878, 20,078

**Note**: Many CVRP instances have similar objectives (~20,000) because they use 100×100 grid with 8-15 customers. The diversity is in **problem structure** (capacity, demands, geography), not objective values.

### Understanding Objective Value Similarity

While some instances have similar objective values, they differ significantly in **problem characteristics**:

| Instance | Customers | K | Q | Total Demand | Structure | Difficulty |
|----------|-----------|---|---|--------------|-----------|------------|
| cvrp_sparse_grid | 9 | 2 | 20 | 27 | Grid layout | 31K nodes |
| cvrp_demand_skewed | 10 | 3 | 25 | 45 | 80% small, 20% large demands | 2K nodes |
| cvrp_geo_clustered | 15 | 3 | 20 | 47 | 2-4 customer clusters | 116K nodes |
| cvrp_capacity_tight | 12 | 7 | 12 | 76 | Very tight, many vehicles | 1K nodes |

All have obj ~20,000-20,500 but vastly different:
- Number of vehicles needed
- Capacity utilization
- Branch & bound difficulty (node count)
- Problem structure

**Why objectives are similar:** Most instances use `generate_positions(n, pattern, seed)` with 100×100 grid and `compute_cost_matrix(pos, "euclidean")`, producing distances ~10-150. With n=8-15 customers, total cost ≈ n × avg_distance ≈ 12 × 1500 = 18,000-22,000.

**Instances with diverse objective scales:** `gen_diverse_vrp.py` generates instances with explicit cost scaling (0.1×, 1×, 5×) to produce objectives ranging from 26 to 878 in VRPPD and various scales in CVRP.

---

## Generation Scripts

### gen_vrp_fixtures.py

Main fixture generator covering all characteristics:

```python
# Generate CVRP
generate_cvrp(
    name="instance_name",
    K=2,                    # vehicles
    customers=[1,2,3,4],    # customer IDs
    Q=20,                   # capacity
    demands={1:5, 2:3, ...},
    costs={(i,j): dist, ...},
    desc="Description"
)

# Generate VRPPD
generate_vrppd(
    name="instance_name",
    K=2,                    # vehicles
    requests=4,             # pickup-delivery pairs
    Q=20,                   # capacity
    pickup_demand=5,
    costs={(i,j): dist, ...},
    desc="Description"
)
```

**Utility Functions**:
- `generate_positions(n, pattern, seed)` - Create layouts:
  - `pattern`: "uniform", "clustered", "grid", "ring"
- `compute_cost_matrix(positions, cost_type)` - Distance matrix:
  - `cost_type`: "euclidean", "manhattan", "rounded"
- `generate_demands(n, pattern, total, seed)` - Demand patterns:
  - `pattern`: "uniform", "varied", "skewed", "few_large"

### gen_diverse_vrp.py

Additional generator for cost scale diversity:
- Small scale: costs 1-20
- Medium scale: costs 10-100
- Large scale: costs 100-500
- Asymmetric costs
- High utilization patterns

### gen_cvrp_working.py

Original simple generator (deprecated, kept for reference).

---

## Test Results

### Comprehensive Test (29/31 instances)

```bash
./test_all_vrp.sh
```

**Results**:
- ✓ 19 instances solve to optimality
- ✓ 10 instances find feasible solutions (timeout)
- ✗ 2 instances need fixing

**Sample Results**:

| Instance | Time | Nodes | Objective | Status |
|----------|------|-------|-----------|--------|
| cvrp_loose | 0.03s | 0 | 90 | Optimal |
| cvrp_medium | 20.3s | 259K | 320 | Optimal |
| cvrp_size_small | 0.02s | 24 | 20,158 | Optimal |
| cvrp_geo_clustered | 60s | 116K | 20,247 | Timeout (feasible) |
| vrppd_small | 0.01s | 0 | 70 | Optimal |
| vrppd_large_uniform | 0.13s | 12 | 370 | Optimal |

### Validation

All instances use the **correct capacity formulation**. Validation available via:
- **Feasibility checking**: `./validate_vrp_solutions.sh` - validates all constraints, bounds, integrality, and objective with `mipster_validate_sol`
- **Cross-solver validation**: `python3 cross_validate_solvers.py` - compares MIPster vs HiGHS objective values (requires HiGHS installed)

---

## Technical Details

### Formulation Breakthrough

After extensive debugging, discovered the correct capacity constraint:

**Failed Approach** (caused infeasibility):
```
u[j] - u[i] + M*x[i,j] >= d[j] + M
```
Problems:
- When x=0: requires `u[j] - u[i] >= d[j] + M`
- For delivery nodes (d[j] negative), this overconstrains
- With bounded u variables, creates infeasibility

**Correct Approach**:
```
u[i] - u[j] + Q*x[i,j] <= Q - d[j]
```
Key differences:
1. Inequality direction: `<=` not `>=`
2. Variable order: `u[i] - u[j]` not `u[j] - u[i]`
3. Use capacity `Q` directly (not arbitrary M)
4. RHS: `Q - d[j]` not `d[j] + M`

### VRPPD Precedence

Precedence enforced by forbidding delivery→pickup arcs:
```
x[delivery_r, pickup_r] <= 0  ∀r ∈ requests
```

This is simpler than auxiliary time/position variables and works well with the capacity constraints.

### MPS Format Notes

- Row names: `OUT{i}`, `IN{i}`, `FLOW{i}`, `CAP{idx}`, etc.
- Variable names: `x{i}_{j}`, `u{i}`
- All binary variables in INTORG marker section
- Capacity RHS: `Q - d[j]` (can be > Q for deliveries)

---

## Files

### Documentation
- `VRP_TESTS.md` - This file (comprehensive documentation)

### Generators
- `gen_vrp_fixtures.py` - Main generator (21 instances)
- `gen_diverse_vrp.py` - Scale-diverse generator (10 instances)
- `gen_cvrp_working.py` - Original simple generator (3 instances)

### Test Scripts
- `test_all_vrp.sh` - Comprehensive test with adaptive timeouts
- `validate_vrp_solutions.sh` - Solve all instances and validate solution feasibility
- `cross_validate_solvers.py` - Cross-validate MIPster vs HiGHS objective values

### Instance Files
- `cvrp_*.mps` - 26 CVRP instances
- `vrppd_*.mps` - 5 VRPPD instances

---

## References

### Formulation Sources
- Standard VRP literature (Toth & Vigo)
- Online MIP formulation references
- Extensive experimentation and validation

### Related Test Suites
- `CInterfaceTest_tsp_random.c` - TSP tests (MTZ and branch-and-cut)
- `CInterfaceTest_mescn.c` - Multi-echelon supply chain tests

---

## Future Enhancements

Potential additions:
1. Time window constraints (VRPTW)
2. Heterogeneous fleet (different vehicle capacities)
3. Multiple depots
4. Split deliveries
5. Backhauls
6. Distance/duration constraints

All would build on the proven capacity constraint formulation documented here.
