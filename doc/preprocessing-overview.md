# MIPster Preprocessing Pipeline — Technical Overview

This document describes the full preprocessing pipeline that MIPster applies
before launching the branch-and-bound search.  It is based on a close reading
of the source code in `Cbc/src/Cgl/CglPreProcess.cpp`, `Cbc/src/Clp/OsiPresolve.cpp`,
`Cbc/src/Cgl/CglProbing.cpp`, and `Cbc/src/CbcSolver.cpp`.

---

## 1. Purpose of Preprocessing

MIP preprocessing transforms the original problem into a smaller, tighter
equivalent that is faster to solve.  It pursues three goals:

1. **Size reduction** — eliminate variables and constraints that are provably
   redundant or fixable.
2. **Bound tightening** — narrow variable domains so the LP relaxation is
   closer to the MIP optimum, reducing the integrality gap.
3. **Strengthening** — add or replace constraints with tighter formulations
   (e.g. substitute inequalities into equalities, identify cliques).

All transformations record enough information to map the optimal solution of
the preprocessed problem back to the original variable space.

---

## 2. Entry Point

Preprocessing is invoked from `CbcSolver::preprocess()` (CbcSolver.cpp,
~line 3409):

```cpp
CglPreProcess process;
// ... configure process, register cut generators ...
OsiSolverInterface *solver2 =
    process.preProcessNonDefault(*saveSolver_, translate[preProcess_],
                                 numberPasses, tunePreProcess_);
```

The returned `solver2` is the (smaller, tighter) model that is handed to
`CbcModel` for branch-and-bound.  If `solver2 == NULL`, preprocessing proved
infeasibility.

### Cut generators registered before preprocessing

Only `CglProbing` is always registered with `process`:

```cpp
CglProbing generator1;
generator1.setUsingObjective(true);
generator1.setMaxPass(3);
generator1.setMaxPassRoot(3);
process.addCutGenerator(&generator1);
```

`CglBKClique` and `CglDuplicateRow` are added conditionally depending on the
problem structure (e.g. whether all-+1 rows exist).

---

## 3. `CglPreProcess::preProcessNonDefault()`

This is the orchestrating function (CglPreProcess.cpp, line 1516).  It runs
a series of phases before a main iterative loop.

### 3.1 Optional Papilo presolve (if compiled in)

When `CBC_USE_PAPILO` is defined and the option is enabled, Papilo's
independent MIP presolve runs first.  Its output is swapped into the Clp
model and the rest of preprocessing continues on the smaller problem.

### 3.2 Mini dual presolve

A lightweight `OsiPresolve::miniPresolvedModel()` pass runs early to apply
cheap dual reductions (e.g. dual fixing of variables with one-sided objective
contributions) before the main presolve.

### 3.3 Cloning and initial setup

`startModel_` is cloned from the input model.  Key bookkeeping arrays are
initialised:

- `originalColumn_[j]` — maps preprocessed column `j` back to the original column.
- `originalRow_[i]` — maps preprocessed row `i` back to the original row.
- `prohibited_[j]` — marks columns that must not be touched (e.g. SOS variables).
- `rowType_[i]` — records special row semantics (objective, SOS rows, etc.).

### 3.4 Initial bound tightening

```
tightenPrimalBounds(*startModel2, false, scBound)
```

`tightenPrimalBounds()` is an `OsiSolverInterface` method that propagates
bound information using the LP solver's own bound-propagation (not the
constraint matrix rewrite used by `tighten()` in CglProbing).  It tightens
row and column bounds based on LP feasibility.

### 3.5 Reduced-cost fixing

```
reducedCostFix(*startModel2)
```

After solving the LP relaxation, variables whose reduced costs guarantee they
will stay at a bound (given the current LP objective) are fixed.  This is
sound because the MIP optimal cannot improve the LP bound.

### 3.6 Clique preprocessing (optional)

When the problem contains all-+1 rows (set-packing structure) or when tuning
enables it, `cliqueIt()` is called.  It uses `CglBKClique` or `CglClique` to:

- Identify cliques (sets of binary variables of which at most one can be 1).
- Replace the original set of covering/packing inequalities with their
  clique-strengthened versions.
- Merge duplicate rows.

The resulting model can be smaller (fewer rows) and tighter.

### 3.7 makeIntegers

`makeIntegers2()` checks whether continuous variables whose LP solution is
always integral (given the structure of the constraints) can be declared
integer, making them available for branching.

### 3.8 Initial LP solve

The strengthened model is solved with the LP simplex before entering the main
pass loop.  If the LP is infeasible, preprocessing returns `NULL` immediately.

---

## 4. Main Iterative Pass Loop

```cpp
for (int iPass = doInitialPresolve;
     iPass < numberSolvers_ && withinTimeLimit;
     iPass++) {
    // VUB analysis
    // OsiPresolve (LP presolve)
    // LP solve
    // reducedCostFix
    // modified()  ← cut generation + strengthening
    // tightenPrimalBounds (if no other changes)
}
```

`numberSolvers_` (the loop bound) equals `numberPasses` from the caller
(default 5 in MIPster).  The special value `99` selects a lightweight
single-pass mode (used by `-preprocess light`).

### 4.1 VUB analysis (Variable Upper Bound)

At the start of each pass, two-element inequality rows of the form
`a·x + b·y ≤ rhs` where one variable appears in only that row are examined.
If the objective drives the singleton variable to its upper bound, the row
can safely be converted to an equality, tightening the formulation.

### 4.2 OsiPresolve (LP presolve)

```cpp
OsiPresolve pinfo;
OsiSolverInterface *presolvedModel =
    pinfo.presolvedModel(*oldModel, feasibilityTolerance,
                         /*keepSolution=*/true, /*maxPasses=*/5,
                         prohibited_, ...);
```

`OsiPresolve` is the LP presolve layer (implemented over `CoinPresolve`).
It applies a sequence of *presolve actions*, each of which reduces the model
size and records the inverse transformation for postsolve.  Key actions:

| Action | What it does |
|---|---|
| **Empty rows/cols** | Remove rows with no nonzeros; fix trivially-bounded columns |
| **Singleton rows** | A row with one variable directly tightens that variable's bounds |
| **Singleton columns** | A column appearing in one row can be substituted out |
| **Fixed variables** | Variables with `LB = UB` are eliminated everywhere |
| **Doubleton rows** | `a·x + b·y = rhs` — substitute one variable, eliminating a row and a column |
| **Implied free variables** | Detect when variable bounds are never active and treat variable as free |
| **Duplicate columns** | Columns with identical constraint coefficients can be merged |
| **Forcing rows** | Rows where all variables are pushed to one bound fix all those variables |
| **Dual fixing** | Variables whose reduced cost sign is determined can be fixed at a bound |

**Important:** LP presolve works with the LP relaxation only.  When it
substitutes a binary variable fixed at 0 by the LP optimum, the resulting
single-element rows are LP-valid but potentially MIP-unsound (the binary may
need to be 1 in the MIP optimum).  This is the root cause of the CglProbing
bug described in §6 below.

After presolve, the model arrays `model_[iPass]` and `presolve_[iPass]` are
saved so postsolve can reconstruct the solution later.

### 4.3 LP solve

The presolved model is solved with primal simplex (`initialSolve()`).  If it
proves infeasible, the outer loop breaks and preprocessing returns `NULL`.

### 4.4 Reduced-cost fixing

Another round of `reducedCostFix()` on the presolved model fixes more
variables that are determined by reduced costs.

### 4.5 `modified()` — cut generation and constraint strengthening

This is the most complex internal function (~2000 lines, line 6431).  It:

1. **Clones** the presolved model to create the next pass's working model.
2. **Sets up probing info** (`CglTreeProbingInfo info`) with `pass = 0`,
   `inTree = false`.
3. **Runs each registered cut generator** (at minimum `CglProbing`) by calling
   `generateCuts()` or `generateCutsAndModify()`.
4. **Applies cuts as row cuts or strengthened rows** — each cut either:
   - tightens an existing row (stored in `whichCut[iRow]`), or
   - is added as a new row.
5. **Counts changes** (`numberChanges`) — used to decide whether to continue
   the outer loop.
6. Does a final LP resolve to get an updated dual solution.

The `constraints` flag (true for all passes except the last) controls whether
row cuts from probing can be added as new constraints, or only used to
strengthen existing rows.

### 4.6 `tightenPrimalBounds` (fallback)

If `modified()` found no changes and no variables were fixed by reduced costs,
`tightenPrimalBounds()` is called once more as a last attempt to find
tightenable bounds.  If free rows are found, the loop extends by one extra
pass.

### 4.7 Loop termination

The loop stops when any of these conditions holds:

- The time limit is exceeded.
- OsiPresolve returns `NULL` (infeasible or unbounded).
- The presolved model has zero rows (problem fully solved).
- Both `numberChanges == 0` and `numberFixed == 0` (no progress).
- `iPass` reaches `numberSolvers_ - 1`.

---

## 5. CglProbing in Preprocessing

`CglProbing::generateCuts()` (CglProbing.cpp) is the main cut generator
called from `modified()`.  During preprocessing it runs with `info.inTree =
false` and `info.pass = 0` (first call from `modified()` in each big pass).

The generator has three sub-components:

### 5.1 `tightPrimalBounds()` block

Runs `OsiSolverInterface::tightPrimalBounds()` — the LP solver's own
propagation on the full LP model.  This is safe because it queries the LP
solver (which uses `getInfinity()` = 1e308 as its infinity threshold) rather
than directly reading the possibly-corrupted constraint matrix.

For each variable whose LP-propagated bound is tighter than its current
model bound, an `OsiColCut` is generated.

This block runs only at `!inTree && !pass` (root, pass 0).

### 5.2 `tighten()` — interval-arithmetic bound propagation

`tighten()` (line 326) propagates bounds through the constraint matrix using
interval arithmetic:

- For each row `a·x ≤ rhs`, computes the minimum and maximum activity
  from current variable bounds.
- Derives tighter bounds: if `max_activity - a[j]*UB[j] > rhs` for some `j`,
  then `UB[j]` must be at most `(rhs - min_activity_without_j) / a[j]`.
- Iterates until no more changes or `maxpass` is reached.

**Critical safety guard (introduced in this session):**

```cpp
if (mode && (info->inTree || info->pass)) {
    ninfeas = tighten(...);
```

`tighten()` is **skipped at pass=0 preprocessing** (`!inTree && pass==0`).
The reason: at this point the constraint matrix may contain single-element
rows like `n[i,t] ≤ 0` introduced by LP presolve when it fixed binary
variables to 0.  These rows are LP-valid but MIP-unsound.  `tighten()` would
propagate them into OsiColCuts that force `UB[n[i,t]] = 0`, which is wrong
when the MIP optimum requires the binary to be 1.

At `pass > 0` or `inTree`, the LP has been re-solved with accumulated cuts
and the constraint matrix is more reliable, so `tighten()` runs normally.

### 5.3 Probing

Probing iterates over binary variables and, for each, does two sub-solves:

1. Fix the binary to 0, propagate bounds, check LP feasibility.
2. Fix the binary to 1, propagate bounds, check LP feasibility.

From the results it extracts:

- **Global cuts** — bounds or row cuts valid regardless of the branching
  direction (both probe directions agree).
- **Implication cuts** — if fixing x=0 forces y=1 (or vice versa), this is
  recorded as an implication.

Probing is safe even at pass=0 because each probe direction is verified
against the LP solver (not just the constraint matrix), so LP-presolve-derived
single-element rows do not cause wrong global cuts.

Root limits: `maxProbeRoot_`, `maxLookRoot_`, `maxElementsRoot_`.
Tree limits: `maxProbe_`, `maxLook_`, `maxElements_`.

---

## 6. The LP-Presolve Constraint Corruption Bug

This section documents a soundness bug found in 2026 and fixed in commit
`19d3e4b8`.

### What happened

For MICLSP instances with Big-M constraints of the form:

```
n[i,t] - M[i,t] * y[i,t] ≤ 0
```

the LP relaxation sometimes found `y[i,t] = 0` optimal.  LP presolve
substituted this into the constraint matrix, creating:

```
n[i,t] ≤ 0   (single-element row)
```

This row is valid for the LP relaxation but wrong for the MIP when
`y[i,t] = 1` is needed in the MIP optimal solution.

### Why it caused wrong cuts

The outer `tighten()` call (before the fix: unconditional `if (mode)`) read
this row and generated an `OsiColCut` setting `UB[n[i,t]] = 0`.  This cut
was sound for the LP but unsound for the MIP.

354 such cuts were generated for instance F (I=4, T=10, seed=137, ρ=1.15).
Once applied, they corrupted subsequent preprocessing passes and the final
MIP solve, causing CBC to certify the wrong objective value 8316 as optimal
(true optimum: 8302).

### The fix

Skip the outer `tighten()` call when `!inTree && pass == 0`:

```cpp
// At preprocessing pass 0 the constraint matrix may contain rows derived
// from LP-presolve variable fixings that are only valid for the LP
// relaxation.  Running tighten() on such a matrix generates unsound MIP
// cuts.  Skip here; probing (below) remains safe because it validates
// each probe direction via LP feasibility.
if (mode && (info->inTree || info->pass)) {
    ninfeas = tighten(...);
```

---

## 7. Postprocessing: Mapping the Solution Back

After branch-and-bound finds an optimal solution to the preprocessed problem,
`CglPreProcess::postProcess()` (line 5601) maps it back to the original
variable space.

The presolve passes are replayed in **reverse order**:

```cpp
for (int iPass = numberSolvers_ - 1; iPass >= 0; iPass--) {
    presolve_[iPass]->postsolve(true);
    // copy solution from model_[iPass] back to model_[iPass-1]
}
```

`OsiPresolve::postsolve()` inverts each presolve action in reverse:
doubleton substitutions are undone, fixed variables are restored, implied
free variables get their values back, etc.

Finally, the solution is written into `originalModel_` (the model the caller
passed to `preProcessNonDefault()`).

`CbcSolver::preprocess()` additionally:
- Fixes integer bounds on the postprocessed solution before calling
  `postProcess()` (rounding LP solution values to the nearest integer where
  appropriate) to avoid subtle mismatches between integer and continuous
  solutions.
- Updates `CbcModel` with the original column indices
  (`model_.setOriginalColumns(process.originalColumns(), ...)`).

---

## 8. Data Flow Summary

```
Original model (OsiSolverInterface)
    │
    ├─ [optional] Papilo presolve
    ├─ Mini dual presolve (OsiPresolve::miniPresolvedModel)
    ├─ tightenPrimalBounds (initial)
    ├─ reducedCostFix (initial)
    ├─ cliqueIt / makeIntegers
    ├─ LP solve
    │
    └─ for iPass in 0..numberPasses-1:
           │
           ├─ VUB analysis (convert inequalities → equalities)
           ├─ OsiPresolve::presolvedModel()
           │      ├─ singleton rows/cols, fixed vars
           │      ├─ doubleton substitution
           │      ├─ implied free, duplicate cols
           │      └─ dual fixing
           ├─ LP solve (presolved model)
           ├─ reducedCostFix (presolved model)
           └─ modified()
                  ├─ CglProbing::generateCuts()
                  │      ├─ tightPrimalBounds() block  [pass=0, !inTree]
                  │      ├─ tighten()  [SKIP at pass=0; run at pass>0 or inTree]
                  │      └─ probing (fix binary, LP check, extract global cuts)
                  ├─ apply cuts / strengthen rows
                  └─ LP resolve
    │
    └─ Preprocessed model → CbcModel (branch-and-bound)
    │
    └─ postProcess() (reverse all OsiPresolve passes, recover solution)
```

---

## 9. Key Parameters

| Parameter | Meaning | Default |
|---|---|---|
| `-preprocess on/off/sos/equal/equalall/strategy/aggregate/light` | Preprocessing mode | `on` |
| `numberPasses` | LP presolve + cut loop iterations | 5 |
| `tunePreProcess_` | Bitmask enabling/disabling sub-features | varies |
| `CglProbing::maxProbeRoot_` | Max binary variables probed at root | 100 |
| `CglProbing::maxLookRoot_` | Max binaries examined at root | 500 |
| `CglProbing::maxPassRoot_` | Tighten passes per call at root | 3 |

---

## 10. Related Source Files

| File | Role |
|---|---|
| `src/Cgl/CglPreProcess.cpp` | Orchestrator — all preprocessing phases |
| `src/Clp/OsiPresolve.cpp` | LP presolve actions and postsolve |
| `src/Cgl/CglProbing.cpp` | Bound propagation (`tighten`) and probing |
| `src/CbcSolver.cpp` | Calls `preProcessNonDefault`, handles postProcess |
| `src/Clp/ClpModel.cpp` | Clp model; normalises UB > 1e27 → COIN_DBL_MAX |
| `src/Osi/OsiSolverInterface.cpp` | `tightPrimalBounds()` implementation |
