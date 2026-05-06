# Diving Heuristics Analysis

*May 2026 — CBC devel*

## Overview

Diving heuristics iteratively fix integer variables and re-solve the LP relaxation,
attempting to "dive" toward a feasible integer solution. Each iteration:

1. Fixes up to `percentageToFix_` (default 20%) of integer variables that are already at integer bounds
2. Selects ONE fractional variable via the heuristic's criterion and tightens its bound
3. Re-solves the LP
4. Repeats until all integers are integral (success) or a stopping condition is hit (failure)

## Default Configuration

- **Only DivingCoefficient is ON by default**
- `percentageToFix_ = 0.2` (20% of integers fixed at bounds per iteration)
- `maxIterations_ = 100` (max dive iterations)
- `maxSimplexIterations_ = 10000` (total simplex iterations budget)
- `diveopt = 2` (d²/2^d probability schedule; always runs at root since depth=0 bypasses the check)

## Empirical Observations

### Test: air05 (7195 integers, 224 fractional at root)

| Heuristic | Iters | Frac Fixed | At-Bound | Simplex Its | Outcome |
|---|---|---|---|---|---|
| VectorLength | 10 | 10 | 2726 | 631 | Infeasible (too many vars forced to 1) |
| Guided | 101 | 101 | 2726 | 240 | Hit iter limit (fixing vars at ~0.0000) |
| Fractional | 101 | 101 | 2726 | 128 | Hit iter limit (fixing vars at ~0.0000) |
| Coefficient | 15 | 15 | 2726 | 979 | Infeasible (mixed directions, moderate fractionality) |
| LineSearch | 1 | 0 | 0 | 0 | Cannot select variable (value == rootValue) |
| PseudoCost | 9 | 9 | 2726 | 467 | Infeasible (high fractionality vars forced to 1) |

### Test: brazil3 (6606 integers after preprocessing, 898 fractional)

| Heuristic | Iters | Frac Fixed | At-Bound | Simplex Its | Outcome |
|---|---|---|---|---|---|
| Fractional | 13 | 13 | 5717 | 710 | Infeasible |
| Coefficient | 12 | 12 | 5725 | 3731 | Infeasible |
| LineSearch | 1 | 0 | 0 | 0 | Cannot select variable |

## Issues Identified

### 1. Guided and Fractional waste iterations at root

Both heuristics fix variables with fractionality < 0.001 (essentially already integer).
They run 100+ iterations fixing variables at values like 0.0000 or 0.0001, accomplishing
nothing meaningful. The LP barely changes because the perturbation is negligible.

**Root cause:** These heuristics select the variable with the smallest fractionality.
At root, many variables have tiny fractionality (numerical noise), and these get selected
first. The heuristic spends all its budget on irrelevant fixings.

### 2. LineSearch is useless at root

LineSearch computes `relDistance = fraction / (rootValue - value)`. At root, `value == rootValue`
for all variables, giving `relDistance = COIN_DBL_MAX`. No variable can be selected, so the
dive bails immediately after 1 iteration.

**Root cause:** The heuristic's criterion is fundamentally undefined at root. It only makes
sense at B&B nodes where variables have moved from their root LP values.

### 3. At-bound fixing is front-loaded

All 2726 (air05) or 5717 (brazil3) at-bound fixings happen in the first iteration.
After that, no more variables are eligible (they're already fixed). The `percentageToFix_`
limit (20%) is not the bottleneck — the number of available candidates is.

### 4. VectorLength and PseudoCost hit infeasibility quickly

These heuristics select variables with high fractionality and round them aggressively
(often all in the same direction — up to 1). In set-partitioning/covering problems,
forcing many variables to 1 quickly violates equality constraints.

### 5. No minimum fractionality threshold

There is no check to skip variables with negligible fractionality. A variable at 0.0001
is essentially integer — fixing it to 0 doesn't help the dive make progress toward
feasibility but still counts as an iteration.

## Proposed Improvements

### A. Add minimum fractionality threshold (low effort, high impact)

Skip variables with fractionality below a threshold (e.g., 0.01) in `selectVariableToBranch()`.
This would prevent Guided and Fractional from wasting iterations on near-integer variables.

```cpp
// In each selectVariableToBranch:
if (fraction < 0.01 && (1.0 - fraction) < 0.01)
    continue;  // skip near-integer variables
```

**Impact:** Guided and Fractional would focus on truly fractional variables, making their
dives more meaningful and potentially finding solutions faster.

### B. Disable LineSearch at root (low effort)

Since LineSearch cannot select any variable at root (all relative distances are infinity),
it should skip execution at root entirely rather than wasting a function call.

```cpp
// In CbcHeuristicDiveLineSearch::selectVariableToBranch or canHeuristicRun:
if (model_->getNodeCount() == 0)
    return true;  // signal all trivially roundable → dive exits gracefully
```

### C. Diversify rounding direction (medium effort, high impact)

VectorLength and PseudoCost tend to round all variables in the same direction (up to 1),
causing quick infeasibility in partitioning problems. A diversification strategy could:

- Limit consecutive same-direction roundings (e.g., after 3 consecutive round-ups, force a round-down)
- Alternate directions based on constraint structure
- Use a "balance" score that penalizes the direction that has been used more

### D. Adaptive percentageToFix_ for all modes (medium effort)

Currently, the adaptive logic (halving `percentageToFix_` after infeasibility) only
activates with `diveopt 8`. Extending it to the default mode would help:

- After infeasibility, reduce at-bound fixing to give the LP more freedom
- Reset on success

### E. New diving heuristic: DiveConflict (high effort, potentially high impact)

A new diving strategy that uses the **conflict graph** to guide variable selection:

**Criterion:** Select the fractional variable whose fixing would propagate the most
implications through the conflict graph. Variables with many conflicts (adjacent nodes
in the conflict graph) are preferred because fixing them eliminates more fractional
variables indirectly.

**Rounding direction:** Round in the direction that creates fewer conflicts with
already-fixed variables.

```
score(x_i) = numConflicts(x_i, direction) × fractionality(x_i)
```

Pick the variable with the highest score — it's both fractional AND well-connected,
so fixing it propagates maximum information.

**Rationale:** CBC already has a conflict graph infrastructure (CglBKClique, CglOddWheel).
Leveraging it for diving would combine the structural information from the constraint
matrix with the LP relaxation state.

### F. New diving heuristic: DiveObjective (medium effort)

A simpler new strategy: select the fractional variable whose fixing would improve the
objective the most per unit of fractionality:

```
score(x_i) = |obj_i| × fractionality(x_i) / columnLength(x_i)
```

Round in the objective-improving direction. This is similar to VectorLength but weights
by fractionality, avoiding the "fix nearly-integer variables" trap.

## Recommended Priority

1. **A (min fractionality threshold)** — immediate improvement, trivial to implement
2. **B (disable LineSearch at root)** — removes wasted computation
3. **D (adaptive percentageToFix_)** — helps all diving heuristics avoid infeasibility
4. **C (diversify direction)** — addresses the main failure mode
5. **E (DiveConflict)** — leverages CBC's unique conflict graph infrastructure
6. **F (DiveObjective)** — simple alternative to VectorLength
