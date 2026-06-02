# gesa3 Preprocessing Bug: Infeasible Solution Generated

**Status:** Unfixed, critical preprocessing bug  
**Affects:** Integer preprocessing (presolve)  
**Instance:** gesa3 from MIPLIB3 (1152 variables, 1368 rows)  
**Discovered:** 2026-06-02

## Symptoms

When solving gesa3 with default settings (preprocessing enabled):
- Solver produces a solution that is **structurally infeasible**
- Row constraint violation: maxViolRow=26.2141 at row 281
- Column bound violation: maxViolCol=19.8421 at col 77
- Objective value 2.79013e+07 is less than known optimal 2.7991e+07 (difference: ~900,000)

**The same instance with preprocessing disabled (`-preprocess off`) solves correctly:**
- Produces feasible solution
- Finds correct optimal objective: 2.7991e+07
- Test passes: `gesa3_o` variant (preprocessing off) PASSES

## Reproduction

```bash
# Fails with preprocessing (default)
$MIPSTER_PREFIX/bin/mipster test/fixtures/gesa3.mps.gz -solve -solu /tmp/gesa3_bad.sol
# Check feasibility - will show violations
$MIPSTER_PREFIX/bin/mipster_validate_sol test/fixtures/gesa3.mps.gz /tmp/gesa3_bad.sol

# Succeeds with preprocessing off
$MIPSTER_PREFIX/bin/mipster test/fixtures/gesa3.mps.gz -preprocess off -solve -solu /tmp/gesa3_good.sol
# Check feasibility - will pass
$MIPSTER_PREFIX/bin/mipster_validate_sol test/fixtures/gesa3.mps.gz /tmp/gesa3_good.sol
```

## Test Output

From `CInterfaceTest_miplib`:

```
[13/40] gesa3                FAIL  (6.9s)
    | [gesa3.gz] FAIL saved sol 0 infeasible: maxViolRow=26.2141 (row 281)  maxViolCol=19.8421 (col 77)
    | [gesa3.gz] FAIL saved sol 0: obj 2.79013e+07 < known_optimal 2.7991e+07 (tol 2799.1)
    | FAIL [CInterfaceTest_miplib.c:98] all saved solutions feasible with valid objective

[14/40] gesa3_o              PASS  (6.9s)    All 1 solution(s) validated OK  (best obj: 2.7991e+07  expected: 2.7991e+07)
```

The `_o` suffix variant disables preprocessing via `-preprocess off` and solves correctly.

## Instance Details

**gesa3.mps.gz:**
- Source: MIPLIB3 benchmark set
- Dimensions: 1152 columns, 1368 rows
- Compressed size: ~31KB
- Known optimal objective: 2.7991e+07
- Available at: `test/fixtures/gesa3.mps.gz`

**Related instance:** gesa2 has similar preprocessing bug but is slightly larger (1224 cols, 1392 rows).

## Root Cause

The bug is in the **integer preprocessing** (presolve) subsystem. Either:

1. **Transformation bug** - preprocessing incorrectly transforms the original problem, creating an equivalent problem that's actually different
2. **Solution mapping bug** - preprocessing correctly transforms the problem, but the solution from the preprocessed model is incorrectly mapped back to the original model's variable space

Given that:
- The violations are large (26.21, 19.84) - not numerical tolerance issues
- The objective is wrong by ~3% (900k out of 28M)
- Multiple constraints and bounds are violated simultaneously
- The bug affects both gesa2 and gesa3 (related instances from same source)

This suggests a **systematic transformation error** rather than a solution mapping bug. Likely candidates:

- Incorrect handling of special constraint types (SOS, GUB, knapsack covers)
- Bug in variable substitution or aggregation during presolve
- Incorrect coefficient updates when eliminating variables or constraints

## Comparison with Node Limit Bug

This is distinct from the node-limit-bug.md issue:
- Node limit bug: parameter passing issue, nodes processed exceed limit
- gesa3 bug: logical correctness issue, solutions are structurally wrong

## Investigation Plan

1. **Capture preprocessed model:**
   ```bash
   # Run with preprocessing, extract the preprocessed MPS
   # (May need to add debug output to save preprocessed model)
   ```

2. **Solve preprocessed model directly:**
   ```bash
   # If preprocessed model is correct, solving it with -preprocess off should give correct result
   # If it fails, the bug is in the transformation
   # If it succeeds, the bug is in solution mapping
   ```

3. **Compare original vs preprocessed:**
   - Variable mappings
   - Constraint transformations
   - Bound changes
   - Coefficient modifications

4. **Use mipster_diag for systematic feature elimination:**
   ```bash
   # This may help narrow down which preprocessing feature causes the bug
   # Though preprocessing is typically all-or-nothing
   ```

5. **Check CbcModel::integerPresolve() and CglPreProcess:**
   - `src/CbcModel.cpp` lines 19271-19709: integer preprocessing flow
   - `src/CglPreProcess.cpp`: constraint preprocessing
   - Look for variable substitution, constraint aggregation, bound propagation

## Workaround

Users encountering infeasible solutions on gesa3 or similar instances should:

```bash
# Disable preprocessing
$MIPSTER_PREFIX/bin/mipster problem.mps.gz -preprocess off -solve
```

Performance may be worse without preprocessing, but results will be correct.

## Related Files

- **Test:** `test/CInterfaceTest_miplib.c` - MIPLIB3 test suite
- **Instance:** `test/fixtures/gesa3.mps.gz` (1152 vars, 1368 rows)
- **Preprocessing code:**
  - `src/CbcModel.cpp` lines 19271-19709: `integerPresolve()`, `integerPresolveThisModel()`, `originalModel()`
  - `src/CglPreProcess.cpp`: constraint-level preprocessing
  - `src/OsiPresolve.cpp`: LP-level presolve

## References

- MIPLIB3: http://miplib.zib.de/miplib3/
- Related instance: gesa2 (similar bug, slightly larger)
- Commit: (to be added when fix is committed)
