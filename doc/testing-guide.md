# MIPster Testing Guide

Comprehensive guide for writing, running, and debugging integration tests for MIP solvers.

---

## Table of Contents

1. [Quick Start](#quick-start)
2. [Test Structure](#test-structure)
3. [Choosing Test Instances](#choosing-test-instances)
4. [Validation Methods](#validation-methods)
5. [Automated Debugging](#automated-debugging)
6. [Writing a New Test](#writing-a-new-test)
7. [CI Integration](#ci-integration)
8. [Best Practices](#best-practices)

---

## Quick Start

### Running the Full Test Suite

```sh
# Build and run all tests
cd test && make -j$(nproc) test

# Or via the test runner script
./test/run-mipster-tests
```

### Running Individual Tests

```sh
cd test

# Run a specific test binary
./CInterfaceTest_vrp fixtures/

# Run with debug library (ASan on Linux)
./run-mipster-tests --debug
```

---

## Test Structure

MIPster uses **integration tests** that verify end-to-end solver correctness on real MIP instances.

### Test Binary Anatomy

Each test (`CInterfaceTest_<name>`) follows this pattern:

```c
#include "Cbc_C_Interface.h"
#include "mip_diag.h"
#include "test_utils.h"

typedef struct {
  const char *name;
  double expected_obj;    // Certified optimal value
  int timeout_sec;
  // ... problem-specific fields
} TestCase;

int main(int argc, char **argv) {
  for each test case:
    1. Load MPS file
    2. Solve with timeout
    3. Check status (optimal/feasible/infeasible)
    4. Validate objective value
    5. Validate solution feasibility
    6. On failure: run mip_diag_wrong_optimal()
}
```

### Key Components

| Component | Purpose |
|---|---|
| `Cbc_C_Interface.h` | C API for loading/solving MIPs |
| `mip_diag.h` | Automated diagnostic passes for wrong-optimal failures |
| `test_utils.h` | Utility functions (timing, tolerances, assertions) |
| `mipster_validate_sol` | External utility for solution feasibility checking |

---

## Choosing Test Instances

### Diversity Dimensions

Good test suites cover multiple characteristics:

| Dimension | Examples | Why Important |
|---|---|---|
| **Problem Size** | 10-1000 variables | Scalability, memory |
| **Variable Types** | Binary, integer, continuous, mixed | Type-specific logic |
| **Constraint Types** | Equalities, inequalities, special structure | Presolve, cuts |
| **Tightness** | LP relaxation gap 0.1%-50% | Integrality handling |
| **Feasibility** | Optimal, feasible-only, infeasible | Status detection |
| **Difficulty** | < 1s to 60s+ | Performance baseline |
| **Special Structure** | Routing, packing, scheduling, networks | Domain-specific cuts |

### Instance Sources

**Internal fixtures** (`test/fixtures/*.mps.gz`):
- Small, fast instances for CI
- Known optimal values
- Cover specific features

**External benchmarks**:
- MIPLIB 2017: `$MIPSTER_INSTANCES/miplib/2017+spp/`
- Diverse set: `$MIPSTER_INSTANCES/super/`

### Generating Synthetic Instances

For domain-specific tests, generate instances programmatically:

```python
# Example: VRP instances with varying characteristics
def generate_vrp_instance(n_customers, capacity, demand_pattern):
    # Generate positions, demands, costs
    # Write MPS file
    return mps_path, expected_obj

# Create diverse suite
instances = [
    generate_vrp_instance(5, 20, "uniform"),
    generate_vrp_instance(10, 15, "skewed"),
    generate_vrp_instance(20, 30, "clustered"),
]
```

See `test/gen_vrp_fixtures.py` for a complete example.

---

## Validation Methods

### 1. Objective Value Validation

Compare solver objective against a certified optimal value:

```c
double expected_obj = 12345.67;
double actual_obj = Cbc_getObjValue(model);
double rel_err = fabs(actual_obj - expected_obj) / fmax(fabs(expected_obj), 1e-6);

if (rel_err > 1e-3) {
  printf("FAIL: obj=%.6f expected=%.6f (rel_err=%.2e)\n",
         actual_obj, expected_obj, rel_err);
  return 1;
}
```

**Obtaining certified values:**
- Solve with multiple solvers and verify agreement
- Use known-optimal instances from literature
- For synthetic instances, verify with exhaustive enumeration (small cases)

### 2. Solution Feasibility Validation

MIPster maintains a **solution pool** during the search. Best practice is to validate **all solutions**, not just the final incumbent.

#### Using the C API

The simplest approach is `Cbc_checkFeasibility()`, which validates:
- **Bounds**: `lb ≤ x[i] ≤ ub` for all variables
- **Integrality**: `|x[i] − round(x[i])| ≤ intTol` for integer variables
- **Constraints**: `rowLB ≤ Ax ≤ rowUB` for all rows (with `primalTol`)

```c
double maxViolRow, maxViolCol;
int rowIdx, colIdx;

const double *solution = Cbc_getColSolution(model);

char is_feasible = Cbc_checkFeasibility(model, solution, 
                                         &maxViolRow, &rowIdx,
                                         &maxViolCol, &colIdx);

if (!is_feasible) {
  if (maxViolRow > 0) {
    printf("FAIL: row %d violated by %.6f\n", rowIdx, maxViolRow);
  }
  if (maxViolCol > 0) {
    printf("FAIL: column %d violated by %.6f\n", colIdx, maxViolCol);
  }
  return 1;
}
```

#### Using mipster_validate_sol (external utility)

For comprehensive validation including objective value checking:

- **Bounds**: `lb ≤ x[i] ≤ ub` for all variables
- **Integrality**: `|x[i] − round(x[i])| ≤ 1e-5` for integer variables
- **Constraints**: `rowLB ≤ Ax ≤ rowUB` for all rows
- **Objective**: computed `c·x` matches claimed value

```c
// Get number of solutions in pool
int n_solutions = Cbc_numberSavedSolutions(model);

// Validate each solution using C API
for (int i = 0; i < n_solutions; i++) {
  const double *solution = Cbc_savedSolution(model, i);
  double sol_obj = Cbc_savedSolutionObj(model, i);
  
  double maxViolRow, maxViolCol;
  int rowIdx, colIdx;
  
  char is_feasible = Cbc_checkFeasibility(model, solution,
                                           &maxViolRow, &rowIdx,
                                           &maxViolCol, &colIdx);
  
  if (!is_feasible) {
    printf("FAIL: solution %d/%d (obj=%.2f) infeasible\n", 
           i+1, n_solutions, sol_obj);
    if (maxViolRow > 0) {
      printf("      row %d violated by %.6f\n", rowIdx, maxViolRow);
    }
    if (maxViolCol > 0) {
      printf("      column %d violated by %.6f\n", colIdx, maxViolCol);
    }
    return 1;
  }
}

printf("PASS: all %d solutions validated\n", n_solutions);
```

**Why validate all solutions:**
- Catches bugs in primal heuristics that produce infeasible solutions
- Verifies solution pool management is correct
- Ensures every reported improvement is actually feasible

**Quick validation (incumbent only):**

For simple tests, validating the final incumbent is often sufficient:

```c
// Using C API (fastest)
const double *solution = Cbc_getColSolution(model);
double maxViolRow, maxViolCol;
int rowIdx, colIdx;

if (!Cbc_checkFeasibility(model, solution, &maxViolRow, &rowIdx,
                          &maxViolCol, &colIdx)) {
  printf("FAIL: incumbent solution infeasible\n");
  return 1;
}

// Or using external validator (more comprehensive, includes objective check)
Cbc_writeSolution(model, "solution.sol");
if (!validate_solution(mps_file, "solution.sol")) {
  printf("FAIL: incumbent solution violates constraints\n");
  return 1;
}
```

**C API Functions:**
- `Cbc_checkFeasibility()` - validates bounds, integrality, and constraints
- `Cbc_getColSolution()` - retrieves incumbent solution
- `Cbc_savedSolution(model, i)` - retrieves solution `i` from pool
- `Cbc_numberSavedSolutions()` - returns pool size

See `test/CInterfaceTest_vrp.c:validate_solution()` for a fork/exec implementation using `mipster_validate_sol`.

### 3. Cross-Solver Validation

Compare objectives across multiple solvers to detect wrong results:

```python
# test/cross_validate_solvers.py
mipster_obj, mipster_status = solve_with_mipster(mps_file)
other_obj, other_status = solve_with_other_solver(mps_file)

if both optimal:
    rel_diff = abs(mipster_obj - other_obj) / max(abs(other_obj), 1e-6)
    if rel_diff > 1e-4:
        print(f"MISMATCH: MIPster={mipster_obj}, Other={other_obj}")
```

This is useful for:
- Verifying correctness of new features
- Regression testing on large instance sets
- Smoke-testing on platforms where certified values aren't available

---

## Automated Debugging

When a test reports **wrong-optimal** (claims optimal at an incorrect objective), use `mip_diag` to identify the culprit.

### How It Works

`mip_diag_wrong_optimal(mps_file, certified_obj, timeout)` runs two diagnostic passes:

1. **Feature sweep** (13 configs): Disables cut families and features one-by-one to isolate the cause
2. **Debug cuts pass**: Activates `OsiRowCutDebugger` against a reference solution

### Integration in Tests

```c
if (test_failed && expected_obj > 0.0) {
  printf("    Running diagnostic...\n");
  mip_diag_wrong_optimal(mps_file, expected_obj, timeout_sec);
}
```

### Reading Diagnostic Output

From CI logs, look for:

```
[DIAG  8/13] no-flow      → not proven  obj=8181  bound=8050  gap=1.6%  *** LEAD ***
[DIAG 12/13] no-twomir    → WRONG       obj=8204  (same wrong result)
```

- `LEAD` = this feature, when disabled, finds the correct objective (strong signal)
- `WRONG` = same wrong result persists (not the cause)
- `OK` = disabling this feature fixed the bug (it was the cause)

Then run `debugCuts` pass to find the exact invalid cut:

```sh
mipster problem.mps -debugCuts known.sol -solve 2>&1 | grep "bad row"
```

Output:
```
bad row 1234 5.0 <= 5.1 <= 5.0
```

This identifies row 1234 as an invalid cut. Cross-reference with `doc/cut-debugging-howto.md` to map the cut type code.

### Manual Diagnostic Flow

For cases not covered by automated diagnostics:

```sh
# 1. Obtain reference solution (preprocessing off)
mipster problem.mps -preprocess off -solve -solu ref.sol

# 2. Run with cut debugging
mipster problem.mps -debugCuts ref.sol -solve 2>&1 | tee debug.log

# 3. Check for violations
grep "bad row\|BAD FIXING" debug.log
```

**Full manual workflow**: `doc/cut-debugging-howto.md`

---

## Writing a New Test

### Step 1: Create Test Fixtures

Generate or collect MPS files with known solutions:

```sh
cd test

# Generate instances
python3 gen_my_problem_fixtures.py

# Compress for CI
gzip -c my_instance.mps > fixtures/my_instance.mps.gz
```

### Step 2: Write Test Binary

Create `test/CInterfaceTest_my_problem.c`:

```c
#include <stdio.h>
#include "Cbc_C_Interface.h"
#include "mip_diag.h"

typedef struct {
  const char *name;
  double expected_obj;
  int timeout_sec;
} TestCase;

static const TestCase test_cases[] = {
  {"my_easy",   100.0,  10},
  {"my_medium", 250.5,  30},
  {"my_hard",     0.0,  60},  // 0.0 = just check solvability
};

int main(int argc, char **argv) {
  if (argc < 2) {
    fprintf(stderr, "Usage: %s <fixture-dir>\n", argv[0]);
    return 1;
  }

  const char *fixture_dir = argv[1];
  int passed = 0, failed = 0;

  for (int i = 0; i < 3; i++) {
    const TestCase *tc = &test_cases[i];
    char path[512];
    snprintf(path, sizeof(path), "%s/%s.mps", fixture_dir, tc->name);

    Cbc_Model *model = Cbc_newModel();
    if (Cbc_readMps(model, path) != 0) {
      printf("FAIL: cannot read %s\n", tc->name);
      failed++;
      continue;
    }

    Cbc_setMaximumSeconds(model, tc->timeout_sec);
    Cbc_solve(model);

    int is_optimal = Cbc_isProvenOptimal(model);
    double obj = Cbc_getObjValue(model);

    if (tc->expected_obj > 0.0) {
      double rel_err = fabs(obj - tc->expected_obj) / tc->expected_obj;
      if (is_optimal && rel_err < 1e-3) {
        printf("PASS: %s obj=%.2f\n", tc->name, obj);
        passed++;
      } else {
        printf("FAIL: %s obj=%.2f expected=%.2f\n", tc->name, obj, tc->expected_obj);
        failed++;
        mip_diag_wrong_optimal(path, tc->expected_obj, tc->timeout_sec);
      }
    } else {
      // Just check solvability
      if (is_optimal || Cbc_numberSolutions(model) > 0) {
        printf("PASS: %s (feasible)\n", tc->name);
        passed++;
      } else {
        printf("FAIL: %s (no solution)\n", tc->name);
        failed++;
      }
    }

    Cbc_deleteModel(model);
  }

  printf("\n%d passed, %d failed\n", passed, failed);
  return (failed == 0) ? 0 : 1;
}
```

### Step 3: Add to Build System

Edit `test/Makefile.am`:

```makefile
bin_PROGRAMS = ... CInterfaceTest_my_problem

CInterfaceTest_my_problem_SOURCES = CInterfaceTest_my_problem.c
nodist_EXTRA_CInterfaceTest_my_problem_SOURCES = dummy.cpp
CInterfaceTest_my_problem_LDADD = ../src/libmipster.la

dist_pkgtestfixture_DATA = \
  ... \
  fixtures/my_easy.mps.gz \
  fixtures/my_medium.mps.gz \
  fixtures/my_hard.mps.gz
```

### Step 4: Add to Test Runner

Edit `test/run-mipster-tests`:

```bash
TESTS=(
  CInterfaceTest
  ...
  CInterfaceTest_my_problem
)
```

### Step 5: Rebuild and Test

```sh
# Regenerate Makefile
autoreconf --install -I BuildTools

# Reconfigure and build
./configster --opt --install

# Run your test
cd test
./CInterfaceTest_my_problem fixtures/

# Run via test runner
./run-mipster-tests
```

---

## CI Integration

Tests run automatically on:
- Linux x86_64 (manylinux_2_34)
- Linux aarch64 (manylinux_2_34)
- Windows x86_64 (MinGW/MSYS2)
- macOS x86_64 (Rosetta 2)
- macOS arm64 (Apple Silicon)

### Test Execution in CI

CI builds run the full test suite via `test/run-mipster-tests`:

1. **Debug build** (with ASan on Linux): runs first, catches memory errors
2. **Release build**: runs on all platforms for performance baseline
3. **Smoke test** (AVX2/neon variants): quick sanity check

### Accessing CI Logs

```sh
gh run list --limit 5
gh run view <run-id> --log | grep "CInterfaceTest_my_problem"
```

### Timeouts in CI

Set per-instance timeouts based on difficulty:
- Easy (< 1s): 10s timeout
- Medium (1-10s): 30s timeout
- Hard (10s+): 60s timeout

CI has a global 2-hour job timeout. Keep test suites under 15 minutes total.

---

## Best Practices

### Instance Selection

✅ **Do:**
- Include instances that **have failed in the past** (regression tests)
- Cover multiple problem types (routing, packing, scheduling, networks)
- Test edge cases (empty constraints, single-variable, infeasible)
- Use small instances (< 1000 variables) for fast CI

❌ **Don't:**
- Add instances that timeout without finding a solution (not useful for validation)
- Include 1000+ large instances (use experiment scripts instead)
- Duplicate coverage (many instances with identical characteristics)

### Objective Value Tolerances

| Value Range | Relative Tolerance |
|---|---|
| Obj > 1e6 | 1e-4 |
| Obj in [1, 1e6] | 1e-3 |
| Obj < 1 | 1e-2 or absolute 1e-6 |

Always use **relative error** for non-zero objectives:
```c
rel_err = fabs(actual - expected) / fmax(fabs(expected), 1e-6)
```

### Test Independence

Each test case must be **independent**:
- No shared state between instances
- Don't rely on test execution order
- Each `Cbc_Model` is fresh (no carried-over parameters)

### Debugging Failures

**Local failure:**
```sh
# Run test with verbose output
RUNNING_TEST=CInterfaceTest_my_problem ./CInterfaceTest_my_problem fixtures/

# Run with ASan
./configster --debug --sanitizer=asan --install
LD_LIBRARY_PATH=$MIPSTER_PREFIX-dbg/lib/dbg ./CInterfaceTest_my_problem fixtures/
```

**CI failure:**
```sh
# Download logs
gh run view <run-id> --log > ci.log

# Extract diagnostic output
grep -A 20 "CInterfaceTest_my_problem" ci.log
```

### Test Maintenance

- **Update expected objectives** when algorithmic improvements change optimal values
- **Remove flaky tests** that intermittently timeout
- **Document non-obvious instances** (why this specific structure matters)

---

## Examples

### Minimal Test (Single Instance)

```c
int main() {
  Cbc_Model *m = Cbc_newModel();
  Cbc_readMps(m, "test.mps");
  Cbc_solve(m);
  
  if (!Cbc_isProvenOptimal(m)) return 1;
  
  double obj = Cbc_getObjValue(m);
  if (fabs(obj - 42.0) > 1e-3) return 1;
  
  Cbc_deleteModel(m);
  return 0;
}
```

### Full-Featured Test

See `test/CInterfaceTest_vrp.c` for:
- Multiple test cases with varying difficulty
- Timeout handling
- Solution feasibility validation via `mipster_validate_sol`
- Automated diagnostics on failure
- Structured output for CI parsing

### Synthetic Instance Generation

See `test/gen_vrp_fixtures.py` for:
- Programmatic MPS generation
- Diverse parameter combinations
- Cost matrix generation (Euclidean, Manhattan, asymmetric)
- Demand pattern generation (uniform, skewed, outliers)

---

## References

- **Cut debugging**: `doc/cut-debugging-howto.md`
- **Diagnostic tool**: `mipster_diag` (`docs/utils.html`)
- **Solution validator**: `mipster_validate_sol` (`docs/utils.html`)
- **C API reference**: `src/Cbc_C_Interface.h`
- **Test runner**: `test/run-mipster-tests`

---

## Quick Reference: Test Checklist

When adding a new test:

- [ ] Test fixtures generated and compressed (`.mps.gz`)
- [ ] Test binary written (`CInterfaceTest_<name>.c`)
- [ ] Added to `test/Makefile.am` (`bin_PROGRAMS`, sources, fixtures)
- [ ] Added to `test/run-mipster-tests` (`TESTS` array)
- [ ] Rebuild: `autoreconf --install -I BuildTools && ./configster --opt --install`
- [ ] Local test: `cd test && ./CInterfaceTest_<name> fixtures/`
- [ ] Full suite: `./run-mipster-tests`
- [ ] CI verification: push and check GitHub Actions
- [ ] Documentation: add entry to this guide if introducing new patterns

---

**Questions or issues?** Open an issue at [h-g-s/mipster](https://github.com/h-g-s/mipster/issues).
