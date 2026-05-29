# Cut Debugging in MIPster

This document explains how to detect invalid row cuts (cuts that incorrectly
exclude the optimal integer solution), what we learned about cut feasibility
and tolerances while fixing a TwoMIR bug (commit `1e871201`), and how the
`-debugCuts` parameter works internally.

---

## The Problem: Invalid Cuts

A row cut `sum_i a_i x_i >= lb` (or `<= ub`) is **invalid** if it excludes
the optimal integer solution. This can happen when:

1. **Floating-point accumulation** in the cut generator produces an RHS that
   is slightly too tight (e.g., `lb = 1.800001144` when the true value is
   `1.8`).
2. **Preprocessing remaps columns**, introducing additional rounding errors
   that compound with step 1.
3. **Scaling routines fail** to clean up the resulting near-rational RHS.

Invalid cuts cause the solver to falsely prune branches containing the true
optimal, yielding a wrong (suboptimal) proven-optimal result.

---

## The `-debugCuts` Parameter

### Usage

```sh
# First obtain a known-good solution (e.g., from a run with preprocessing off):
~/prog/cbc/bin/mipster problem.mps -preprocess off -solve -solu /tmp/known_opt.sol

# Then solve normally with debug mode active:
~/prog/cbc/bin/mipster problem.mps -debugCuts /tmp/known_opt.sol -solve
```

The `.sol` file is the standard MIPster solution file (output of `-solu`).
If any cut generated during the B&B is violated by the known solution, a
diagnostic is printed immediately:

```
bad row 1751 1.800001144 <= 1.8 <= 1.797693135e+308
```

Format: `bad row <row_idx> <rowLower> <= <solution_activity> <= <rowUpper>`

A "bad row" means the cut's bound cannot be satisfied by the known optimal,
i.e., the cut is invalid.

### What it does internally

1. Reads the `.sol` file via `CbcMipStart::read()` — lines starting with a
   digit are `idx colname value` entries.
2. Stores the solution in `debugValues_[nOrigCols]`.
3. After preprocessing, `CbcSolver.cpp` (around line 10122) maps the original
   column values through `process.originalColumns()`, fixes integer variables
   to their known values, resolves the LP, then calls
   `babModel_->solver()->activateRowCutDebugger(newValues)`.
4. `CbcModel` checks `specialOptions_ & 1` (set when debugger is active) and
   in `CHECK_KNOWN_SOLUTION` blocks (always compiled in via
   `OsiRowCutDebugger.hpp:7`) evaluates every cut against the known solution.
5. Violations are printed and — if `CoinAssert` is enabled — can abort the
   solve immediately.

### Note on `continuousSolver_`

The debugger is activated on `babModel_->solver()` (the B&B LP solver) but
NOT on `continuousSolver_`. This means the `CoinAssert(!debugger->invalidCut())`
check at `CbcModel.cpp:8532` is **not** triggered (it gates on
`continuousSolver_->getRowCutDebugger()`). The "bad row" printout from line
11288–11335 still fires and is sufficient for diagnosis.

---

## Cut Feasibility and Tolerances

### Key tolerance values

| Tolerance | Value | Used for |
|---|---|---|
| Primal feasibility | `1e-7` | LP solution is considered feasible |
| Integer feasibility | `1e-5` | Variable considered integer |
| `OsiRowCutDebugger::invalidCut()` epsilon | `1e-6` | Cut is invalid if optimal violates by > this |
| `CHECK_KNOWN_SOLUTION` check | `1e-6` | Prints "bad row" if activity < lb − 1e-6 |

A cut that is wrong by less than `1e-6` will **not** be caught by the debugger.

### How `scaleCutIntegral` works (in `OsiCuts.cpp`)

After a cut is generated, `insertIfNotDuplicateAndClean` calls
`scaleCutIntegral(elements, indices, nz, rhs, maxdelta=1e-12)` which:

1. Uses `CoinRational(x, maxdelta, maxdnom=60)` to find rational `p/q` with
   `|x - p/q| ≤ maxdelta` and `q ≤ 60`, for the RHS and all coefficients.
2. Finds the LCM of all denominators → this is the scale factor.
3. If the scaled values all round to integers (within `1e-9`), the cut is
   replaced by its integer-scaled version (e.g., `5x >= 9` instead of
   `x >= 1.8`).
4. Returns `false` if any value can't be approximated or the scale exceeds
   1000.

**Why this fails for our bug:** `maxdelta = 1e-12` is far too tight to
recognize `1.800001144 ≈ 9/5` (error = `1.144e-6 >> 1e-12`). The scaling
returns false and the raw erroneous RHS is kept.

### The TwoMIR floating-point error (fixed in commit `1e871201`)

Instance **nu25-pr12**: MIPster reported optimal=**53910** (wrong),
true optimal=**53905**.

Root cause chain:
1. Preprocessing remaps the problem to a smaller formulation.
2. `CglTwomir` generates a valid cut (mathematically) but with accumulated
   FP error: `lb = 1.800001144` instead of `1.8`.
3. `scaleCutIntegral` fails (error >> maxdelta), so the raw `lb` is kept.
4. The cut incorrectly excludes the optimal solution (violation = `1.144e-6`
   barely exceeds the `1e-6` threshold).
5. The branch containing the optimal is pruned; solver reports 53910 as optimal.

The bug is triggered only when TwoMIR AND preprocessing are both active.
Disabling either one avoids the problem because:
- Without preprocessing: the original formulation's cut is cleaner.
- Without TwoMIR: no TwoMIR cuts are generated.

**Fix (`OsiCuts.cpp`):** When `scaleCutIntegral` fails for TwoMIR cut types
(61–63), relax the cut's RHS bound by `2e-6`:
```cpp
if (!goodScale && typeCut >= 61 && typeCut <= 63) {
    if (newLb < -1.0e30)
        newUb += 2.0e-6;
    else
        newLb -= 2.0e-6;
}
```

The `2e-6` is just above the observed violation (`1.144e-6`) and within
the solver's feasibility tolerances.

---

## Systematic Diagnosis Workflow

When MIPster reports a suspiciously good (or wrong) result:

### 1. Establish a ground truth

```sh
# Disable preprocessing
~/prog/cbc/bin/mipster problem.mps -preprocess off -solve

# Disable individual cut types
~/prog/cbc/bin/mipster problem.mps -twomir off -solve
~/prog/cbc/bin/mipster problem.mps -gomory off -solve
~/prog/cbc/bin/mipster problem.mps -mir off -solve
```

If results differ, a cut type is the likely culprit. If only the
preprocessing+cut combination is wrong (not either alone), the bug is an
interaction between preprocessing's column remapping and cut generation.

### 2. Save the known good solution

```sh
~/prog/cbc/bin/mipster problem.mps -preprocess off -solve -solu /tmp/known.sol
```

### 3. Run with `-debugCuts` and look for "bad row" output

```sh
~/prog/cbc/bin/mipster problem.mps -debugCuts /tmp/known.sol -solve 2>&1 | grep "bad row"
```

If you see "bad row" lines, the identified cut type is invalid.

### 4. Narrow to the cut generator

```sh
# Enable only TwoMIR, use debugCuts
~/prog/cbc/bin/mipster problem.mps \
  -gomory off -mir off -clique off \
  -debugCuts /tmp/known.sol -solve 2>&1 | grep "bad row"
```

### 5. Add `assert`-on-bad-cut for crash diagnostics

To make the solver abort at the first invalid cut (instead of just printing),
activate the debugger also on `continuousSolver_`. Currently (as of 2026),
`babModel_->solver()` is activated but `continuousSolver_` is not — the
`CoinAssert` at `CbcModel.cpp:8532` is therefore silent.

---

## Cut Type Codes in `insertIfNotDuplicateAndClean`

| Type | Cut generator |
|---|---|
| 11–17 | Gomory |
| 21–22 | GMI |
| 31 | MIR2 |
| 41–43 | FlowCover |
| 51–52 | RedSplit2 |
| 61–63 | TwoMIR |
| 71–72 | ResidualCapacity |

Types 61–63 are the three TwoMIR insertion paths:
- **61**: MIR cut path (rejects badly-scaled cuts)
- **62**: complemented variable variant
- **63**: main TwoMIR path (the one with the bug, now fixed)
