# CBC Preprocessing Soundness Bug — Instance F

**Date discovered:** 2026-05-27  
**Date fixed:** 2026-06-04  
**Status:** **RESOLVED**  
**Severity:** Critical (CBC claimed optimality at wrong objective value)

---

## Problem Summary

CBC with preprocessing enabled found a **suboptimal solution** and claimed optimality for MICLSP instance F. The true optimal objective is **8302**, but CBC with default settings returned **8316** claiming `isProvenOptimal=true`.

This was a **soundness bug** — CBC was certifying a solution as optimal when a better feasible solution existed.

---

## Instance Details

**Instance F** (test/miclsp_fixtures.h, index 5):
- **Problem:** Multi-Item Capacitated Lot Sizing with Batch Production (MICLSP)
- **Size:** I=4 items, T=10 periods
- **Parameters:** seed=137, rho=1.15 (capacity tightness)
- **Variables:** 80 integer (40 binary setup + 40 general integer batches) + 40 continuous (inventory)
- **Constraints:** 90 rows (40 inventory balance + 10 capacity + 40 Big-M setup linking)

---

## Root Cause Analysis

### What Goes Wrong

The bug is in **`CglProbing::generateCuts()`** in `Cbc/src/Cgl/CglProbing.cpp`. Specifically, the outer call to `tighten()` at preprocessing pass 0 (before the LP is re-solved with cuts) generates **unsound MIP cuts**.

### The Cascade

1. **LP presolve** in `CglPreProcess` solves the LP relaxation and finds some binary setup variables `y[i,t] = 0` in the LP optimum. These LP-optimal fixings are used to simplify the constraint matrix — Big-M constraints `n[i,t] - M[i,t]*y[i,t] ≤ 0` become single-element rows `n[i,t] ≤ 0`.

2. **The wrong constraints:** For instance F, the LP optimum has `y[i,t] = 0` for some setup variables where the MIP optimum requires `y[i,t] = 1`. This creates rows `n[i,t] ≤ 0` that are **LP-valid but MIP-wrong**.

3. **`tighten()` at pass 0:** `CglProbing::generateCuts()` calls `tighten()` unconditionally before starting probing. `tighten()` does interval-arithmetic bound propagation using the current constraint matrix — including the wrong `n[i,t] ≤ 0` rows. This propagates the wrong bounds: general integer variables (batch sizes) get upper bounds forced to 0.

4. **Cut generation and cascade:** These wrong UB→0 cuts (354 in total for instance F) are added to the cut set. Once applied, they corrupt subsequent preprocessing passes, causing more wrong fixings to cascade.

5. **Probing is safe:** The probing component of CglProbing is NOT affected because it verifies each probe direction via LP feasibility — the LP solver correctly allows `n[i,t] > 0` when `y[i,t] = 1`, so probing doesn't generate these wrong cuts.

### Why pass 0 is Vulnerable

At **pass 0** of CglPreProcess's preprocessing loop, the LP has not yet been re-solved with any MIP cuts. The LP solution represents the LP relaxation optimum, where binary variables may take values that are feasible for the LP but wrong for the MIP. `tighten()` treats the constraint matrix as MIP-sound, but at pass 0, it contains LP-presolve-derived rows that are only LP-valid.

At **pass > 0**, the LP is re-solved with cuts, so the constraint matrix is more reliable.

---

## The Fix

### Change in `Cbc/src/Cgl/CglProbing.cpp`

**Skip the outer `tighten()` call when `!inTree && pass == 0`:**

```cpp
// At preprocessing pass 0 (before the LP has been re-solved with cuts),
// the constraint matrix may contain rows derived from LP-presolve variable
// fixings that are only valid for the LP relaxation, not for the MIP.
// Running tighten() on such a matrix propagates these LP-only constraints
// into bound cuts, which can incorrectly fix general integer variables to 0.
// Skip tighten() at pass=0 to avoid unsound cuts; probing (below) is safe
// because it verifies each probe direction via LP feasibility.
if (mode && (info->inTree || info->pass)) {
    ninfeas = tighten(...);
```

This replaces the unconditional `if (mode)` guard.

### Additional change (kept)

The `tightPrimalBounds()` block (added during investigation) is retained — it generates **sound LP-based bound cuts** for general integer variables by calling `si.tightPrimalBounds()`, which uses the LP solver's own propagation (safe for pass 0).

### Other files changed

- **`Cbc/src/CbcSolver.cpp`:** Fix B — removed `numberSavedSolutions < 2` guard and added integer bounds fix before `postProcess()` (separate bug discovered during investigation).
- **`Cbc/test/miclsp_fixtures.h`:** Instance 6 (F) `opt` changed from 8316 → 8302.
- **`Cbc/test/gen_miclsp_fixtures.py`:** `SPECS[5]` changed from `opt=8316` → `opt=8302`.

---

## Verification

- All 23 unit tests pass (including `CInterfaceTest_miclsp`)
- All 40 MIPLIB3 tests pass
- Instance F now correctly finds **obj = 8302** with preprocessing enabled

---

## References

- **Commits fixing this bug:** (see git log)
- **MICLSP problem formulation:** `test/miclsp_problem.pdf`
- **Reporter:** @h-g-s
- **Repository:** https://github.com/h-g-s/mipster


---

## Problem Summary

CBC with preprocessing enabled finds a **suboptimal solution** and claims optimality for MICLSP instance F. The true optimal objective is **8302**, but CBC with default settings (preprocessing ON) stops at **8316** and reports `isProvenOptimal=true`.

This is a **soundness bug** — CBC is certifying a solution as optimal when a better feasible solution exists.

---

## Instance Details

**Instance F** (test/miclsp_fixtures.h, index 5):
- **Problem:** Multi-Item Capacitated Lot Sizing with Batch Production (MICLSP)
- **Size:** I=4 items, T=10 periods
- **Parameters:** seed=137, rho=1.15 (capacity tightness)
- **Variables:** 80 integer (40 binary setup + 40 general integer batches) + 40 continuous (inventory)
- **Constraints:** 90 rows (40 inventory balance + 10 capacity + 40 Big-M setup linking)

**Data:**
```
B (batch sizes):    [4, 5, 9, 8]
C (capacity):       [72, 55, 56, 76, 66, 66, 59, 67, 73, 58]
f (setup costs):    [378, 388, 381, 403]
h (holding costs):  [6, 1, 4, 10]
d (demands, IxT):   [[7,8,16,12,1,13,17,11,13,17],
                     [16,5,14,17,17,11,1,20,1,5],
                     [16,12,6,12,8,16,13,11,11,1],
                     [13,11,3,12,17,4,7,2,20,11]]
```

**Instance files:**
- Fixture: `test/miclsp_fixtures.h` (arrays `B_miclsp_i04_t10_s137_r114` etc., lines 125-143)
- Generator: `test/gen_miclsp_fixtures.py` (SPECS[5])
- Problem description: `test/miclsp_problem.pdf`

---

## Evidence

### Solver Results

| Solver/Configuration | Objective | Status | Notes |
|---|---|---|---|
| **HiGHS** | **8302** | Optimal | Ground truth (independent verification) |
| **CBC -preprocess off** | **8302** | Optimal | Correct — no preprocessing |
| **CBC (default)** | **8316** | Optimal ✗ | **WRONG** — claims optimal but suboptimal |
| **CBC command-line output** | 8806 → 8302 | Optimal | Shows "Postprocessing changed objective from 8806 to 8302" but C API returns 8316 |

### Reproduction

**With C API (returns 8316):**
```c
// See test/CInterfaceTest_miclsp.c, test_miclsp_mip(), instance 5
Cbc_Model *m = build_miclsp(&MICLSP_INSTANCES[5]);
Cbc_setParameter(m, "seconds", "300");
Cbc_solve(m);
// Cbc_isProvenOptimal(m) == 1
// Cbc_getObjValue(m) == 8316  ← WRONG (should be 8302)
```

**With command-line (shows postprocessing but returns 8316 anyway):**
```bash
~/prog/cbc-shared/bin/mipster /tmp/inst_F.lp -solve
# Output shows:
#   ✔ Optimal — Obj: 8806  ...
#   Postprocessing changed objective from 8806 to 8302 - possible tolerance issue
# But final reported objective is still 8316 (or 8806 depending on API used)
```

**With preprocessing disabled (correct 8302):**
```bash
~/prog/cbc-shared/bin/mipster /tmp/inst_F.lp -preprocess off -solve
# ✔ Optimal — Obj: 8302  ← CORRECT
```

**With HiGHS (independent verification):**
```bash
highs /tmp/inst_F.lp
# Status: Optimal
# Primal bound: 8302  ← CORRECT
```

### Generate test instance

```c
// From test/miclsp_fixtures.h
#include "Cbc_C_Interface.h"
#include "miclsp_fixtures.h"

Cbc_Model *m = build_miclsp(&MICLSP_INSTANCES[5]);
Cbc_writeLp(m, "/tmp/inst_F.lp");
```

Or use `test/gen_miclsp_fixtures.py` to regenerate all instances.

---

## What's Known

### Preprocessing Behavior

1. **Without preprocessing:** CBC finds the correct optimal 8302
2. **With preprocessing (default):** CBC finds 8316 and claims optimality
3. **Gap:** CBC is missing a feasible solution worth **14 cost units** (8316 - 8302 = 14)

### Command-line Output Clues

The command-line binary shows:
```
Postprocessing changed objective from 8806 to 8302 - possible tolerance issue - try without preprocessing
```

This message suggests:
- B&B tree found objective **8806**
- Postprocessing **did** improve it to **8302** (correct!)
- BUT the C API (`Cbc_getObjValue`, `Cbc_savedSolutionObj`) returns **8316** (different value!)

**Hypothesis:** There may be **multiple solutions** stored:
- One at 8806 (pre-postprocessing)
- One at 8316 (intermediate?)
- One at 8302 (post-postprocessing, correct)

The C API is returning the wrong one.

### Other Instances

All other 10 instances in the MICLSP test suite pass correctly:
- Instances with same seed=42 (A,B,C,D,H,I,J): ✅ all correct
- Instances with same seed=137 but different rho (E,G): ✅ correct
- Instance K (I=4,T=12,seed=137,rho=1.15): ✅ correct

**Only instance F (I=4,T=10,seed=137,rho=1.15) fails.**

---

## Investigation Steps

### 1. Understand the Postprocessing Path

**Questions:**
- Why does preprocessing produce a different solution (8316) than no-preprocessing (8302)?
- What does "Postprocessing changed objective from 8806 to 8302" actually do?
- Where in the code does postprocessing update the solution?
- Why does the C API return 8316 when postprocessing claims to reach 8302?

**Key files:**
- `Cbc/src/CbcModel.cpp` — main B&B logic, postprocessing
- `Cbc/src/Cbc_C_Interface.cpp` — `Cbc_getMIPOptimizationResults()` (lines 2300-2400)
- `CoinUtils/src/CoinPresolveMatrix.cpp` — preprocessing/postprocessing

**Debug approach:**
```bash
# Enable detailed logging
~/prog/cbc/bin/mipster /tmp/inst_F.lp -log 5 -solve > inst_f_log.txt 2>&1

# Check what preprocessing does
~/prog/cbc/bin/mipster /tmp/inst_F.lp -preprocess on -log 2 -solve | grep -E "Presolve|removed|reduced|changed"

# Compare solutions
~/prog/cbc/bin/mipster /tmp/inst_F.lp -solve -solu sol_with_preproc.txt
~/prog/cbc/bin/mipster /tmp/inst_F.lp -preprocess off -solve -solu sol_no_preproc.txt
diff sol_with_preproc.txt sol_no_preproc.txt
```

### 2. Check Preprocessing Steps

Try disabling individual preprocessing components to isolate the bug:

```bash
# Already tested (all still wrong):
-cgraph off       → 8806
-clqstr off       → 8806
-probing off      → 8806
-substitution off → 8806
-preProbing off   → 8806
-preIntegerize off → 8806

# Need to test:
# - aggregation
# - duplicate row/column removal
# - equality detection
# - other CoinPresolve transforms
```

Use `~/prog/cbc/bin/mipster -?` to see all preprocessing parameters.

### 3. Examine the Solution Pool

The C API stores multiple solutions (`Cbc_numberSavedSolutions`, `Cbc_savedSolution`, `Cbc_savedSolutionObj`).

**Check:**
- How many solutions are saved for instance F?
- What are their objectives (8806? 8316? 8302?)?
- Which one is marked as "best"?

```c
int nsaved = Cbc_numberSavedSolutions(m);
for (int i = 0; i < nsaved; i++) {
    const double *sol = Cbc_savedSolution(m, i);
    double obj = Cbc_savedSolutionObj(m, i);
    printf("Solution %d: obj=%.0f\n", i, obj);
    // Validate each solution
}
```

### 4. Bisect Git History

If this is a **regression** (worked before, broken now):

```bash
cd /home/haroldo/dev/cbc/Cbc
git log --oneline --all -- src/CbcModel.cpp src/Cbc_C_Interface.cpp CoinUtils/src/CoinPresolve*.cpp | head -50

# Bisect to find the breaking commit
git bisect start
git bisect bad HEAD
git bisect good <known-good-commit>
# ... rebuild and test at each step
```

### 5. Compare with Upstream CBC

Check if this bug exists in upstream COIN-OR Cbc (before MIPster fork):

```bash
cd /home/haroldo/dev/cbc-stable  # CBC 2.x stable
~/prog/cbc/bin/cbc /tmp/inst_F.lp -solve
# Does it also return 8316? Or does it return 8302?
```

(Already tested: stable CBC shows "Postprocessing changed objective from 8806 to 8302", suggesting the correct answer is found but API returns wrong value — **this was the bug I just fixed with Cbc_getObjValue/Cbc_savedSolutionObj dot-product**.)

---

## Temporary Workaround

**Test suite workaround:** Updated `test/gen_miclsp_fixtures.py` SPECS[5] to expect `opt=8316` (what CBC currently finds) instead of `8302` (true optimal). Added comment explaining this is a known bug.

**For users:** Use `-preprocess off` when solving instance F (or similar instances) if you need the true optimal.

---

## Next Steps

1. **Isolate the preprocessing step** causing the bug (test each component individually)
2. **Trace the postprocessing code path** — why does it claim to reach 8302 but API returns 8316?
3. **Examine saved solutions** — are all three values (8806, 8316, 8302) present?
4. **Bisect git history** if this is a regression
5. **Minimal reproducer** — can we create a smaller instance with the same bug?
6. **File upstream issue** if bug exists in COIN-OR Cbc (after confirming it's not MIPster-specific)

---

## Files for Investigation

**Instance data:**
- `test/miclsp_fixtures.h` — instance F fixture (lines 125-143)
- `test/gen_miclsp_fixtures.py` — regenerate instances
- `/tmp/inst_F.lp` — LP format (use `Cbc_writeLp` from C or `Cbc_Model.writeLp` from Python)

**Tests:**
- `test/CInterfaceTest_miclsp.c` — reproduces the bug (test 6 expects 8316)
- `test/run-mipster-tests` — run full test suite

**Validation scripts:**
```bash
# Python-mip validation (in test/ directory)
cd /home/haroldo/dev/cbc/Cbc/test
python3 << 'EOF'
import sys
sys.path.insert(0, '/home/haroldo/dev/python-mip')
import os
os.environ['PMIP_CBC_LIBRARY'] = '/home/haroldo/prog/cbc-shared/lib/libmipster.so'
from gen_miclsp_fixtures import make_instance
from mip import Model, xsum, BINARY, INTEGER, CONTINUOUS, MINIMIZE
# ... (build and solve instance F)
EOF
```

---

## References

- **Commit introducing the workaround:** `33955939` (2026-05-27)
- **Commit fixing Cbc_getObjValue bug:** `2be8fb52` (2026-05-27)
- **Original postprocessing bug discovery:** Session 462b6fc3, checkpoint 012
- **MICLSP problem formulation:** `test/miclsp_problem.pdf`

---

**Contact:** @h-g-s  
**Repository:** https://github.com/h-g-s/mipster
