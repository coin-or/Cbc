# Pure Binary Problem Test Fixtures

This document describes the test fixture generators for pure binary (0/1) optimization problems, designed to test MIPster's clique cut generators (`CglBKClique`, `CglOddWheel`, `CglCliqueStrengthening`) and other binary-specific features.

## Overview

We have **7 problem types** with diverse instance characteristics:

| Problem | Generator | Key Features Tested |
|---------|-----------|---------------------|
| **QAP** | `gen_qap_fixtures.py` | Quadratic assignment, McCormick linearization, clique structure |
| **Set Covering** | `gen_setcover_fixtures.py` | Hitting sets, dominated columns, clique cuts |
| **Set Packing** | `gen_setpacking_fixtures.py` | Maximum weighted independent set, conflict cliques |
| **Set Partitioning** | `gen_setpartitioning_fixtures.py` | Exact covers, highly degenerate, equality constraints |
| **Bin Packing with Conflicts** | `gen_bpp_conflicts_fixtures.py` | Knapsack + conflicts, clique + capacity cuts |
| **Max Independent Set** | `gen_mis_fixtures.py` | Pure conflict graph, clique/odd-hole cuts |
| **MDKP** | `gen_mdkp_fixtures.py` | Multi-dimensional knapsack (already exists) |

---

## Diversity Dimensions

Each generator creates instances with controlled variation across these dimensions:

### 1. Size (rows/columns)
- **Small** (10-20 variables): CI quick tests (< 30s)
- **Medium** (30-60 variables): CI moderate tests (30-120s)
- **Large** (80-150 variables): Experiment-only (skipped in CI)

### 2. Density
- **Sparse** (10-30% non-zeros): Test cut sparsification
- **Medium** (30-60% non-zeros): Typical real-world density
- **Dense** (60-90% non-zeros): Stress-test clique finding

### 3. Structure
- **Random**: Baseline, no exploitable structure
- **Clustered/Clique**: k-partite, clique subgraphs
- **Geometric**: Spatial patterns (unit disk, grid)
- **Hierarchical**: Dominated sets, nested structure

### 4. Tightness
- **Loose**: LP relaxation close to optimal (< 5% gap)
- **Medium**: 5-20% integrality gap
- **Tight**: > 20% gap, requires strong cuts

### 5. Objective Distribution
- **Uniform**: random weights in [1, 100]
- **Diverse**: mix of small/medium/large weights
- **Unit**: all weights = 1 (unweighted problem)
- **Structured**: weights correlated with problem structure

---

## Problem-Specific Details

### 1. QAP (Quadratic Assignment Problem)

**File:** `gen_qap_fixtures.py`

**Formulation:** Linearized Kaufman-Broeckx (2-index with McCormick auxiliary variables)

**Diversity:**
| Dimension | Values |
|-----------|--------|
| Size | n ∈ {4, 5, 6, 7} facilities/locations |
| Flow pattern | uniform, sparse (30%), hub-and-spoke |
| Distance pattern | euclidean, manhattan, grid |
| Seed | 42, 137 (deterministic variation) |

**Why diverse:**
- Flow patterns test different conflict densities
- Distance patterns vary constraint tightness
- Small sizes ensure solvability (QAP is NP-hard, grows very hard)

**Expected solve time:** 10-120s

---

### 2. Set Covering

**File:** `gen_setcover_fixtures.py`

**Formulation:** Minimize cost to cover all elements

**Diversity:**
| Dimension | Values |
|-----------|--------|
| Size | elements ∈ {10, 15, 20, 30, 40, 50}, sets ∈ {15, 20, 30, 50, 60, 80} |
| Pattern | uniform, sparse, dense, hierarchical |
| Density | 30% (uniform), 20% (sparse), 70% (dense) |
| Seed | 42, 137 |

**Pattern details:**
- **uniform**: random 30% coverage
- **sparse**: each set covers 20% of elements (low overlap)
- **dense**: each set covers 70% of elements (high overlap → many conflicts)
- **hierarchical**: some sets subsume others (test dominated column detection)

**Expected solve time:** 10-90s

---

### 3. Set Packing

**File:** `gen_setpacking_fixtures.py`

**Formulation:** Maximize value of non-conflicting sets

**Diversity:**
| Dimension | Values |
|-----------|--------|
| Size | elements ∈ {10, 12, 15, 20, 30, 35, 40, 50}, sets ∈ {15, 18, 20, 25, 50, 60, 70, 80} |
| Pattern | random, k_partite, geometric, scheduling |
| Conflict density | 20-40% (pattern-dependent) |
| Value dist | uniform, diverse, unit, structured |
| Seed | 42, 137 |

**Pattern details:**
- **random**: random k-element sets, conflicts when elements overlap
- **k_partite**: sets partitioned into k cliques (all pairs within clique conflict)
- **geometric**: intervals on a line, conflicts when intervals overlap
- **scheduling**: time windows, conflicts on overlapping slots

**Expected solve time:** 10-90s

---

### 4. Set Partitioning

**File:** `gen_setpartitioning_fixtures.py`

**Formulation:** Minimize cost to partition elements (each in exactly one set)

**Diversity:**
| Dimension | Values |
|-----------|--------|
| Size | elements ∈ {10, 15, 20, 30, 40, 50}, sets ∈ {20, 30, 40, 60, 80, 100} |
| Pattern | uniform, overlapping, structured |
| Seed | 42, 137 |

**Pattern details:**
- **uniform**: each element appears in 2-4 sets
- **overlapping**: base partition + alternative covers
- **structured**: graph partitioning on grid (contiguous regions)

**Expected solve time:** 10-120s

**Note:** Set partitioning often highly degenerate (many optimal solutions).

---

### 5. Bin Packing with Conflicts

**File:** `gen_bpp_conflicts_fixtures.py`

**Formulation:** Minimize bins, items have weights and conflicts

**Diversity:**
| Dimension | Values |
|-----------|--------|
| Size | items ∈ {10, 12, 15, 20, 25, 30}, capacity ∈ {100, 120, 150, 180, 200} |
| Conflict pattern | random, clique, bipartite, geometric |
| Conflict density | 20-30% |
| Weight dist | uniform, diverse (small/med/large), tight |
| Seed | 42, 137 |

**Pattern details:**
- **random**: random item pairs conflict
- **clique**: items partitioned into k cliques (mutual exclusion sets)
- **bipartite**: two types, conflicts only across types
- **geometric**: items have positions, conflicts if close

**Weight distributions:**
- **uniform**: random in [capacity/10, capacity/3]
- **diverse**: mix of small (5%), medium (30%), large (65%) items
- **tight**: most items around capacity/3 (tighter packing)

**Expected solve time:** 20-120s

---

### 6. Maximum Independent Set (MIS)

**File:** `gen_mis_fixtures.py`

**Formulation:** Maximize weighted vertices with no adjacent pairs selected

**Diversity:**
| Dimension | Values |
|-----------|--------|
| Size | vertices ∈ {15, 16, 18, 20, 30, 35, 36, 40} |
| Graph pattern | random, geometric, grid, k_partite, tree, cycle, planar |
| Edge density | 30% (random/geometric), structural (others) |
| Weight dist | uniform, unit, diverse, degree-based |
| Seed | 42, 137 |

**Pattern details:**
- **random**: Erdős-Rényi G(n, p)
- **geometric**: unit disk graph (proximity in 2D)
- **grid**: 2D lattice (4-connected)
- **k_partite**: k-partite graph
- **tree**: random spanning tree (polynomial-time solvable, good validation)
- **cycle**: cycle graph
- **planar**: grid + some diagonals

**Expected solve time:** 5-90s (tree instances very fast, good for validation)

---

### 7. Multi-Dimensional Knapsack (MDKP)

**File:** `gen_mdkp_fixtures.py` (already exists)

**Formulation:** Maximize profit subject to multiple capacity constraints

**Diversity:**
| Dimension | Values |
|-----------|--------|
| Size | m ∈ {2, 3, 5} dimensions, n ∈ {20, 25, 30, 40} items |
| Tightness | α = 0.5 (capacity = 50% of total weight) |
| Seed | 42, 137 |

**Expected solve time:** < 60s

---

## Usage

### Generate All Fixtures

```bash
cd test

# Requires external MIP solver for certified optimal values
python3 gen_qap_fixtures.py
python3 gen_setcover_fixtures.py
python3 gen_setpacking_fixtures.py
python3 gen_setpartitioning_fixtures.py
python3 gen_bpp_conflicts_fixtures.py
python3 gen_mis_fixtures.py

# Existing generator
python3 gen_mdkp_fixtures.py
```

Each script:
1. Generates instances with diverse characteristics
2. Solves with external solver to obtain certified optimal values
3. Writes MPS files to `fixtures/*.mps.gz`
4. Prints C test case snippet for integration

### Add to Test Suite

Each generator prints a C test case array at the end:

```c
typedef struct {
  const char *name;
  double expected_obj;
  int timeout_sec;
} ProblemTestCase;

static const ProblemTestCase test_cases[] = {
  {"instance_name", 12345.67, 30},
  // ...
};
```

Copy this into `test/CInterfaceTest_<problem>.c` (see existing tests for full template).

---

## Design Principles

### 1. Non-Trivial but Solvable
- **Avoid**: instances that solve in < 1s (trivial, don't test anything)
- **Avoid**: instances that timeout (>300s) or don't converge
- **Target**: 10-120s solve time per instance in CI

### 2. Diverse Coverage
Each problem type covers:
- Multiple sizes (small/medium)
- Multiple structures (random/structured/geometric)
- Multiple densities (sparse/medium/dense)
- Multiple seeds (deterministic variation)

### 3. Clique-Cut Friendly
All problems chosen because they benefit from:
- Clique cuts (`CglBKClique`)
- Odd-wheel cuts (`CglOddWheel`)
- Clique strengthening (`CglCliqueStrengthening`)

### 4. Cross-Validation
All instances solved with external MIP solver to obtain **certified optimal values**.
Tests fail if MIPster:
- Claims optimal at wrong objective
- Produces infeasible solution
- Times out on known-solvable instance

---

## Testing Guidelines

### Instance Selection for CI

| Time Budget | Size Range | Patterns | Seeds |
|-------------|-----------|----------|-------|
| 10-20s | Small (10-20 vars) | All | 1 |
| 30-60s | Medium (20-40 vars) | Representative | 1-2 |
| 60-120s | Medium-Large (40-60 vars) | Hard patterns | 1 |

**CI suite target:** 8-12 instances per problem type, total CI time < 15 min.

### Instance Selection for Experiments

Use all generated fixtures (including large sizes) for:
- Performance regression testing
- Cut generator tuning
- Algorithm comparison

---

## Adding a New Problem Type

1. **Create `gen_<problem>_fixtures.py`** following existing pattern:
   ```python
   def generate_instance(n, seed, pattern, ...):
       # Generate problem data
       return data

   def write_mps(filename, data):
       # Write MPS format

   def solve_external(data, timeout):
       # Solve with external solver, return (obj, status)

   SPECS = [
       # (size, seed, pattern, ..., timeout)
   ]

   def main():
       # Generate, solve, compress, print C snippet
   ```

2. **Diversity checklist:**
   - [ ] 3+ size values (small/medium/large)
   - [ ] 3+ structural patterns
   - [ ] 2+ density levels
   - [ ] 2+ seeds for each config
   - [ ] Expected solve time 10-120s

3. **Write C test** `test/CInterfaceTest_<problem>.c`:
   - Load MPS
   - Solve with timeout
   - Validate objective + feasibility
   - Call `mip_diag_wrong_optimal()` on failure

4. **Add to build** (`test/Makefile.am`, `test/run-mipster-tests`)

---

## Maintenance

### Updating Certified Optimal Values

If algorithmic improvements in external solver or MIPster change known optimal values:

1. Re-run generator script
2. Check new optimal values are **better or equal** (never worse)
3. Update test case arrays in `CInterfaceTest_*.c`

### Removing Flaky Instances

If an instance consistently times out or is too easy:

1. Identify in CI logs (`gh run view --log`)
2. Remove from SPECS in generator script
3. Regenerate fixtures
4. Update test arrays

---

## References

- **Testing guide:** `doc/testing-guide.md`
- **Cut debugging:** `doc/cut-debugging-howto.md`
- **Existing tests:** `test/CInterfaceTest_mdkp_random.c`, `test/CInterfaceTest_tsp_random.c`
