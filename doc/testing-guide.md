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
| `test_utils.h` | Shared helpers: `fixture_path()`, `solution_path()`, `validate_all_saved_solutions()`, `mip_obj_tol()`, `build_mps_model()`, `print_perf_summary()` |
| `mipster_validate_sol` | External utility for solution feasibility checking |

#### `test_utils.h` — key helpers

```c
/* Locate a fixture MPS file (honours $MIPSTER_FIXTURE_DIR env var) */
const char *fixture_path(const char *fname);      // e.g. "my_inst.mps.gz"

/* Locate a reference solution file (fixtures/solutions/<base>.sol) */
const char *solution_path(const char *base);      // NULL if file absent

/* Validate every solution in the pool — three checks per solution:
 *   1. Feasibility (bounds, integrality, rows via Cbc_checkFeasibility)
 *   2. Objective consistency: computed c·x must match Cbc_savedSolutionObj
 *   3. Bound: pool obj must not beat certified_opt (within obj_tol)
 * Pass NAN for known_opt to skip checks 2b/3. Returns failure count. */
int validate_all_saved_solutions(Cbc_Model *m, double known_opt,
                                 double obj_tol, const char *tag);

/* Compute effective objective tolerance (ratio gap + abs gap + base_tol) */
double mip_obj_tol(Cbc_Model *m, double certified_opt, double base_tol);

/* Generic mip_diag builder: reads model from MPS path in userdata */
Cbc_Model *build_mps_model(void *userdata);   // pass (void*)path_str
```

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

`mip_diag_wrong_optimal(builder, userdata, certified_obj, timeout)` runs:

1. **Feature sweep** (27 configs): Disables cut families, preprocessing, tolerances, and platform-specific options one-by-one to isolate the cause
2. **Debug cuts pass**: Activates `OsiRowCutDebugger` against a reference solution to catch individually-invalid cuts

### Integration in Tests

```c
/* mip_diag requires a builder function — not a path string directly */
static Cbc_Model *my_builder(void *ud) {
  const char *path = (const char *)ud;
  Cbc_Model *m = Cbc_newModel();
  if (Cbc_readMps(m, path) != 0) { Cbc_deleteModel(m); return NULL; }
  return m;
}

/* On failure, call with builder + path + certified opt + time limit */
if (is_proven_optimal && fabs(obj - certified_opt) > tol) {
  mip_diag_wrong_optimal(my_builder, (void *)path,
                         certified_opt, timeout_sec);
  /* Optional: if a reference solution file exists, also run debugCuts */
  const char *sol = solution_path(tc->name);
  if (sol)
    mip_diag_debug_cuts(my_builder, (void *)path,
                        certified_opt, timeout_sec, sol, NULL, NULL);
}
```

> **Tip:** `test_utils.h` provides `build_mps_model()` which is exactly the builder above. Pass it directly: `mip_diag_wrong_optimal(build_mps_model, (void *)path, ...)`

### Reading Diagnostic Output

From CI logs, look for:

```
[DIAG  8/27] no-mir        → OK  obj=55 certified=55
[DIAG 12/27] no-probing    → OK  obj=55 certified=55
[DIAG 24/27] probing+mir   → WRONG opt=56 (expected 55) — BUG ALSO PRESENT
```

- `OK` = disabling this feature fixed the bug (it is the cause, or at least a co-trigger)
- `WRONG` = same wrong result persists (not the sole cause)
- `LEAD` = found the correct obj without proving optimality (strong suspect)
- `*** LEAD ***` = confirmed correct obj — highest-priority suspect

**If debugCuts shows no `bad row` lines:** The culprit is not an individually-invalid cut. Instead, suspect **LP bound inflation** — floating-point rounding (especially FMA on Apple Silicon/clang) inflates the LP relaxation bound by a few ULPs above the fathoming threshold. Look for `OK` results on `tight-prim-tol`, `scaling-off`, or `tight-cutoff-inc` configs; these pinpoint precision-sensitive code paths.

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

# Generate instances with a Python generator script
python3 gen_my_problem_fixtures.py

# Compress for CI
gzip -c my_instance.mps > fixtures/my_instance.mps.gz

# If the test needs reference solutions for debugCuts, store them in:
#   fixtures/solutions/my_instance.sol
```

### Step 2: Write Test Binary

Create `test/CInterfaceTest_my_problem.c`. Follow the canonical pattern used
by existing tests (`CInterfaceTest_upms.c`, `CInterfaceTest_jssp_fixtures.c`):

```c
#include <math.h>
#include <stdio.h>
#include "Cbc_C_Interface.h"
#include "mip_diag.h"
#include "test_utils.h"

typedef struct {
  const char *name;
  double certified_opt;   /* known optimal; NaN if unknown */
  int max_nodes;          /* MIP primary stop — deterministic across platforms */
  int timeout_sec;        /* loose wall-clock fallback only */
} TestCase;

static const TestCase test_cases[] = {
  /* name                certified_opt  max_nodes  timeout_sec */
  { "my_easy",           100.0,         500000,    60  },
  { "my_medium",         250.5,         500000,    120 },
  { "my_hard",           NAN,           500000,    300 },  /* NAN = any feasible */
};
#define N_TESTS ((int)(sizeof(test_cases)/sizeof(test_cases[0])))

int main(int argc, char **argv)
{
  if (argc < 2) { fprintf(stderr, "Usage: %s <fixture-dir>\n", argv[0]); return 1; }
  int pass = 0, fail = 0;

  for (int i = 0; i < N_TESTS; i++) {
    const TestCase *tc = &test_cases[i];
    char path[512];
    snprintf(path, sizeof(path), "%s/%s.mps.gz", argv[1], tc->name);

    Cbc_Model *m = Cbc_newModel();
    Cbc_setLogLevel(m, 0);
    if (Cbc_readMps(m, path) != 0) {
      printf("  FAIL: cannot read %s\n", tc->name); fail++; continue;
    }

    /* MIP stopping criteria:
     *   - max nodes: PRIMARY deterministic stop (same result on all platforms)
     *   - time limit: loose wall-clock fallback only (never rely on it for
     *     correctness; it exists only to guard against pathological slowdowns)
     * LP tests (CbcSolverLpTest etc.): no limits — LP solves always finish. */
    Cbc_setMaximumNodes(m, tc->max_nodes);
    Cbc_setMaximumSeconds(m, tc->timeout_sec);
    Cbc_solve(m);

    int ok = 1;
    int is_optimal = Cbc_isProvenOptimal(m);
    double obj = Cbc_getObjValue(m);

    /* 1. Validate entire solution pool:
     *    - feasibility (bounds, integrality, rows)
     *    - objective consistency (computed c·x == reported obj)
     *    - bound check (no pool solution may beat certified_opt) */
    ok &= (validate_all_saved_solutions(m, tc->certified_opt,
                                        mip_obj_tol(m, tc->certified_opt, 1e-4),
                                        tc->name) == 0);

    /* 2. If solver claims proven optimality, verify against certified value */
    if (is_optimal && !isnan(tc->certified_opt)) {
      double tol = mip_obj_tol(m, tc->certified_opt, 1e-6);
      if (fabs(obj - tc->certified_opt) > tol) {
        printf("  FAIL: %s proven opt=%.6g expected=%.6g (tol=%.4g)\n",
               tc->name, obj, tc->certified_opt, tol);
        ok = 0;
        mip_diag_wrong_optimal(build_mps_model, (void *)path,
                               tc->certified_opt, tc->timeout_sec);
        const char *sol = solution_path(tc->name);
        if (sol)
          mip_diag_debug_cuts(build_mps_model, (void *)path,
                              tc->certified_opt, tc->timeout_sec,
                              sol, NULL, NULL);
      }
    }

    if (ok) { printf("  PASS %s\n", tc->name); pass++; } else { fail++; }
    Cbc_deleteModel(m);
  }

  printf(fail == 0 ? "=== All %d tests PASSED ===\n"
                   : "=== %d passed, %d FAILED ===\n", pass, fail);
  return (fail == 0) ? 0 : 1;
}
```

### Step 3: Add to Build System

Edit **four** places in `test/Makefile.am`:

```makefile
# 1. Add to bin_PROGRAMS (one long line — append to it)
bin_PROGRAMS = ... CInterfaceTest_my_problem

# 2. Add sources + LDADD at the end of the file
CInterfaceTest_my_problem_SOURCES = CInterfaceTest_my_problem.c
nodist_EXTRA_CInterfaceTest_my_problem_SOURCES = dummy.cpp
CInterfaceTest_my_problem_CPPFLAGS = $(AM_CPPFLAGS) -DFIXTURE_DIR=\"$(FIXTURE_DIR_VAL)\"
CInterfaceTest_my_problem_LDADD = ../src/libmipster.la

# 3. Add fixture files to dist_pkgtestfixture_DATA
dist_pkgtestfixture_DATA = \
  ... \
  fixtures/my_easy.mps.gz \
  fixtures/my_medium.mps.gz \
  fixtures/my_hard.mps.gz

# 4. If you have reference solution files, add them to dist_pkgtestsol_DATA
dist_pkgtestsol_DATA = \
  ... \
  fixtures/solutions/my_easy.sol

# 5. Add to the ctests target (needed for 'make test' without the runner script)
ctests: ... CInterfaceTest_my_problem$(EXEEXT)
```

### Step 4: Add to Test Runner

Edit **two** arrays in `test/run-mipster-tests`:

```bash
# ── Test list ──────────────────────────────────────────────────────────────
TESTS=(
  ...
  CInterfaceTest_my_problem    # ← add here
)

# IMPORTANT: this line MUST immediately follow the closing ) of TESTS.
# It was the cause of a real CI failure when accidentally deleted.
TOTAL=${#TESTS[@]}

# ── If your test takes <fixture-dir> as its first argument ─────────────────
# Add it to FIXTURE_TESTS as well, or it will be run without the fixture path:
FIXTURE_TESTS=(
  ...
  CInterfaceTest_my_problem    # ← add here if your main() takes argv[1]
)
```

> ⚠️ **`TOTAL` gotcha:** `TOTAL=${#TESTS[@]}` must appear right after the `TESTS=(...)` block. If you accidentally delete it while editing the array, all tests will fail with `TOTAL: unbound variable`.

> ⚠️ **`FIXTURE_TESTS` gotcha:** If your binary takes `<fixture-dir>` as its first argument (i.e., `int main(int argc, char **argv)`), it **must** be in `FIXTURE_TESTS` too, or it will be invoked without the fixture path and fail silently.

### Step 5: Rebuild and Test

```sh
# Regenerate Makefile.in from Makefile.am
autoreconf --install -I BuildTools

# Reconfigure and build
./configster --opt --install

# Run your new test directly
cd test && ./CInterfaceTest_my_problem fixtures/

# Run the full suite
./test/run-mipster-tests
```
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

### Stopping Criteria

| Test type | Primary stop | Secondary stop |
|---|---|---|
| **MIP** | `Cbc_setMaximumNodes(m, N)` | `Cbc_setMaximumSeconds(m, T)` (loose fallback) |
| **LP** | none — LP solves always finish | none |

**Why max nodes for MIP:** node count is deterministic across platforms (same branching decisions → same nodes explored). Time limits produce different results on different machines and sanitizer builds.

**Never** use time as the primary MIP stop in CI tests. A reasonable node budget is 500 000 — enough to find and verify optimal for small-to-medium fixtures; accepting "not proven" as a valid outcome when the node limit is hit.

### Solution Validation (all test types)

Every test must validate **all solutions in the pool**, not just the best incumbent:

```c
/* validate_all_saved_solutions checks per solution:
 *   1. feasibility (bounds, integrality, all rows)
 *   2. objective consistency: computed c·x == Cbc_savedSolutionObj()
 *   3. bound: pool obj does not beat certified_opt
 * Returns number of failures. */
ok &= (validate_all_saved_solutions(m, certified_opt,
                                    mip_obj_tol(m, certified_opt, 1e-4),
                                    instance_name) == 0);

/* Additionally: if solver claims proven optimality, check the value */
if (Cbc_isProvenOptimal(m) && !isnan(certified_opt)) {
  double tol = mip_obj_tol(m, certified_opt, 1e-6);
  if (fabs(Cbc_getObjValue(m) - certified_opt) > tol) {
    /* wrong-optimal bug — run mip_diag */
  }
}
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

See `test/CInterfaceTest_upms.c` or `test/CInterfaceTest_jssp_fixtures.c` for the current canonical pattern:
- Node limit as primary termination, time limit as fallback
- `validate_all_saved_solutions()` for full pool validation
- `mip_obj_tol()` for tolerance computation
- `mip_diag_wrong_optimal()` + `mip_diag_debug_cuts()` on failure
- `build_mps_model` builder passed to `mip_diag_*` functions

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

- [ ] Test fixtures generated and compressed (`.mps.gz`) in `test/fixtures/`
- [ ] Reference solutions (if any) in `test/fixtures/solutions/<name>.sol`
- [ ] Test binary written (`CInterfaceTest_<name>.c`) — uses node limit as primary stop, validates full solution pool, calls `mip_diag_wrong_optimal` on failure
- [ ] `test/Makefile.am` updated in **five** places:
  - [ ] `bin_PROGRAMS` — append binary name
  - [ ] Binary sources + `_CPPFLAGS` + `_LDADD` block
  - [ ] `dist_pkgtestfixture_DATA` — list `.mps.gz` files
  - [ ] `dist_pkgtestsol_DATA` — list `.sol` files (if any)
  - [ ] `ctests` target — append `CInterfaceTest_<name>$(EXEEXT)`
- [ ] `test/run-mipster-tests` updated in **two** places:
  - [ ] `TESTS=(...)` array — append binary name
  - [ ] `TOTAL=${#TESTS[@]}` line is intact right after the closing `)` ← easy to accidentally delete
  - [ ] `FIXTURE_TESTS=(...)` array — append name **if binary takes `<fixture-dir>` as first arg**
- [ ] Rebuild: `autoreconf --install -I BuildTools && ./configster --opt --install`
- [ ] Local test: `cd test && ./CInterfaceTest_<name> fixtures/`
- [ ] Full suite: `./test/run-mipster-tests`
- [ ] CI verification: push and check GitHub Actions

---

**Questions or issues?** Open an issue at [h-g-s/mipster](https://github.com/h-g-s/mipster/issues).
