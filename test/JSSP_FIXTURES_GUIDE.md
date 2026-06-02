# Job Shop Scheduling Problem (JSSP) Test Fixtures

## Overview

This test suite validates MIPster's handling of **disjunctive scheduling constraints** using classic JSSP benchmark instances. Unlike the pure binary problems (QAP, set packing, etc.), JSSP tests **big-M formulations** with mixed binary and continuous variables, making it an important complement to the clique-cut test suite.

## Problem Structure

**JSSP formulation (disjunctive/big-M):**

Variables:
- `C`: makespan (continuous, objective to minimize)
- `x[j][op]`: start time of operation `op` of job `j` (continuous)
- `y[j][k][m]`: binary precedence variable (1 if job `j` precedes job `k` on machine `m`)

Constraints:
1. **Precedence within job**: `x[j][op] + p[j][op] <= x[j][op+1]`
2. **Makespan definition**: `x[j][last_op] + p[j][last_op] <= C`
3. **Disjunctive (no overlap on machines)**:
   - `x[j][op_j] + p[j][op_j] <= x[k][op_k] + M*(1 - y[j][k][m])`
   - `x[k][op_k] + p[k][op_k] <= x[j][op_j] + M*y[j][k][m]`

where `M` = sum of all processing times (big-M constant).

## What Makes JSSP Interesting for Testing

1. **Different from clique-cut problems**: JSSP uses **big-M disjunctive constraints** rather than conflict graph structures, so it doesn't trigger CglBKClique, CglOddWheel, or CglCliqueStrengthening.

2. **Tests sequencing logic**: Binary variables represent ordering decisions between pairs of operations on the same machine.

3. **Mixed variable types**: Combines continuous (start times) and binary (precedence) variables.

4. **Classic benchmark suite**: Fisher & Thompson, Lawrence, and Applegate & Cook instances are well-studied with known optimal makespans.

## Test Instances

| Instance | Size (jobs×machines) | Optimal Makespan | Timeout | Source |
|----------|---------------------|------------------|---------|--------|
| jssp_ft06 | 6×6 | 55 | 30s | Fisher & Thompson |
| jssp_ft10 | 10×10 | 930 | 300s | Fisher & Thompson |
| jssp_la01 | 10×5 | 666 | 120s | Lawrence |
| jssp_la06 | 15×5 | 926 | 180s | Lawrence |
| jssp_la11 | 20×5 | 1222 | 300s | Lawrence |
| jssp_orb01 | 10×10 | 1059 | 300s | Applegate & Cook |

### Diversity

**Size diversity:**
- Small: 6×6 (ft06)
- Medium: 10×5, 10×10 (la01, ft10, orb01)
- Large: 15×5, 20×5 (la06, la11)

**Structure diversity:**
- **Fisher & Thompson (ft06, ft10)**: Balanced square instances (n×n)
- **Lawrence (la01, la06, la11)**: Rectangular instances (more jobs than machines)
- **Applegate & Cook (orb01)**: Square with diverse processing time distributions

**Processing time diversity:**
- Fisher & Thompson: Integer processing times, wide range
- Lawrence: Integer processing times, structured distributions
- Applegate & Cook: Large processing times, high variance

## Validation

Each test:
1. **Checks optimality**: `Cbc_isProvenOptimal()`
2. **Validates feasibility**: `Cbc_checkFeasibility()` (checks all precedence, disjunctive, and makespan constraints)
3. **Verifies objective**: Makespan must match known optimal value within ±0.5

## Generator

Run `gen_jssp_fixtures.py` to regenerate fixtures:

```bash
cd test
python3 gen_jssp_fixtures.py
```

The generator embeds classic benchmark data (machine orderings, processing times) and produces MPS files using the disjunctive/big-M formulation.

##References

- Fisher, H., & Thompson, G. L. (1963). "Probabilistic learning combinations of local job-shop scheduling rules"
- Lawrence, S. (1984). "Supplement to Resource constrained project scheduling"
- Applegate, D., & Cook, W. (1991). "A computational study of the job-shop scheduling problem"

Classic JSSP data available at: http://jobshop.jjvh.nl/
