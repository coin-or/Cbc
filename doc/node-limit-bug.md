# Node Limit Bug: Integer Preprocessing Bypasses User-Specified Limit

**Status:** Critical bug, unfixed  
**Affects:** All solves with integer preprocessing enabled (default)  
**Workaround:** Use time limits instead of node limits  
**Discovered:** 2026-06-02

## Symptoms

When a user sets a node limit via `-maxNodes` parameter or `Cbc_setMaximumNodes()` C API:
- The solver processes **significantly more nodes** than the limit specifies
- Example: `-maxNodes 1000` → 62,487 nodes processed
- Example: `-maxNodes 1000` → 3,633,351 nodes processed (with restarts)
- The behavior is **non-deterministic** depending on instance characteristics

Time limits work correctly and are not affected by this bug.

## Reproduction

```bash
# Command line
$MIPSTER_PREFIX/bin/mipster test/fixtures/bpc_n30_c180_sd42_geometric_tight.mps.gz \
  -maxNodes 1000 -solve

# Observed output:
# ✔ Stopped (node limit) — BestSol: 10   Bound: 8.45555   Gap: 15.44%   Nodes: 62.5K
# Enumerated nodes:               62487

# C API (similar behavior)
Cbc_Model *m = Cbc_newModel();
Cbc_readMps(m, "problem.mps.gz");
Cbc_setMaximumNodes(m, 1000);
Cbc_solve(m);
printf("Nodes processed: %d\n", Cbc_getNodeCount(m));  // Prints 62487
```

The solver correctly detects it reached the node limit (status message confirms), but the actual node count far exceeds the limit.

## Root Cause Analysis

### Parameter Flow (C API)

The node limit is set correctly through the full chain:

1. **C API:** `Cbc_setMaximumNodes(model, 1000)` sets `model->int_param[INT_PARAM_MAX_NODES] = 1000`
2. **Input queue:** C API converts to command-line parameter: `-maxN=1000`
3. **Parameter parsing:** CbcParameters reads `-maxN=1000` and sets `parameters_[CbcParam::MAXNODES]`
4. **Model synchronization:** `synchronizeModel()` transfers to `cbcModel.intParam_[CbcMaxNumNode] = 1000`
5. **Limit checking:** `isNodeLimitReached()` correctly checks `numberNodes_ >= intParam_[CbcMaxNumNode]`

All parameters are set correctly. The issue is **not** in parameter handling.

### Integer Preprocessing Creates Independent Models

When integer preprocessing is enabled (default), the solve flow is:

```
1. Main CbcModel created by user
   - intParam_[CbcMaxNumNode] = 1000
   - numberNodes_ = 0
   
2. Initial root node processing (heuristics, cuts)
   - numberNodes_ remains 0 (or small)
   
3. integerPresolve() called (CbcModel.cpp:19271)
   - Creates: newModel = new CbcModel(*this)  // Copy constructor
   - Copies: intParam_[CbcMaxNumNode] = 1000  ✓
   - But: newModel->numberNodes_ = 0  ← PROBLEM
   
4. integerPresolveThisModel() (line 19302)
   - Calls: branchAndBound() on newModel (line 19675)
   - newModel processes nodes 0 → 62,487
   - newModel's isNodeLimitReached() checks: numberNodes_ >= 1000
   - But numberNodes_ starts at 0, so limit is effectively ignored
   
5. originalModel() copies results back (line 19619)
   - numberNodes_ = presolvedModel->numberNodes_  // Overwrites with 62,487
```

**Key Issue:** Each model has its **own independent `numberNodes_` counter**, and the preprocessed model's counter starts at 0 regardless of nodes already processed.

### Restarts Compound the Problem

After some nodes, the preprocessed model may decide to restart, creating yet another new model:

```
1. Preprocessed model processes nodes 0 → N
2. Restart decision triggered
3. New model created for restart
   - numberNodes_ = 0 again!
4. Continues processing more nodes
5. Multiple restarts → each adds more nodes
```

This explains instances where millions of nodes are processed with a 1000-node limit.

### Debug Evidence

Added debug output to `isNodeLimitReached()` and `stoppingCriterionReached()`:

```
DEBUG: isNodeLimitReached() nodes=1000 limit=1000 reached=1
DEBUG: stoppingCriterionReached() called, nodes=1000 result=1
DEBUG: isNodeLimitReached() nodes=1000 limit=1000 reached=1
DEBUG: stoppingCriterionReached() called, nodes=1000 result=1
DEBUG: isNodeLimitReached() nodes=1000 limit=1000 reached=1
Enumerated nodes:               62487
```

Observations:
- The check fires correctly at 1000 nodes and returns `true`
- Only **3 checks occur** total (not every 1000 nodes as expected in main loop)
- This confirms the checks stop after the main model hits its limit
- But the solve continues in the preprocessed model with its own counter

### Code Locations

**Integer preprocessing model creation:**
```cpp
// src/CbcModel.cpp:19276
CbcModel *newModel = new CbcModel(*this);  // Copy constructor

// src/CbcModel.cpp:6640 (copy constructor)
memcpy(intParam_, rhs.intParam_, sizeof(intParam_));  // Copies limit ✓

// But numberNodes_ inherits from source or is initialized to 0
```

**Preprocessed model solve:**
```cpp
// src/CbcModel.cpp:19675
branchAndBound();  // Called on preprocessed model
```

**Results copied back:**
```cpp
// src/CbcModel.cpp:19698 (originalModel function)
numberNodes_ = presolvedModel->numberNodes_;  // OVERWRITES instead of ADDING
```

**Restart/subBranchAndBound code paths:**
```cpp
// src/CbcModel.cpp:20306 and 20413
model2->setIntParam(CbcModel::CbcMaxNumNode, maximumNodes);
model2->branchAndBound();
```

The `subBranchAndBound()` function takes a `maximumNodes` parameter and explicitly sets it, but this is for heuristics/special solves, not the main preprocessing path.

## Why Time Limits Work

Time limits are checked differently:
- Single global `CoinCpuTime()` or wall-clock timer
- Checked in `maximumSecondsReached()` 
- Same timer used across all models/restarts
- No per-model counter to reset

## Related: Recent Parameter Refactoring

The codebase has `#define CBC_CLUMSY_CODING` sections (CbcParameters.hpp:66) that handle parameter synchronization:

```cpp
// src/CbcSolver.cpp:9124
#ifdef CBC_CLUMSY_CODING
  parameters.synchronizeModel();  // Transfers parameters to CbcModel
#endif
```

This refactoring was done to centralize parameter handling through CbcParameters instead of direct model setters. However, the node limit bug predates this refactoring - it's an architectural issue with how models are copied, not how parameters are set.

## Attempted Fixes (Did Not Work)

### 1. Set parameter directly on CbcModel before run()

```cpp
// src/Cbc_C_Interface.cpp:2720
cbcSolver.model()->setIntParam(CbcModel::CbcMaxNumNode, 
                               model->int_param[INT_PARAM_MAX_NODES]);
cbcSolver.run(inputQueue, cbc_callb);
```

**Result:** No effect. The preprocessed model still has its own counter.

### 2. Adjust limit in subBranchAndBound

```cpp
// src/CbcModel.cpp:20305
int effectiveMaxNodes = maximumNodes;
if (intParam_[CbcMaxNumNode] != COIN_INT_MAX) {
  int remaining = intParam_[CbcMaxNumNode] - numberNodes_;
  if (remaining < effectiveMaxNodes)
    effectiveMaxNodes = std::max(0, remaining);
}
model2->setIntParam(CbcModel::CbcMaxNumNode, effectiveMaxNodes);
```

**Result:** No effect. This path is for heuristic subproblems, not the main preprocessing flow.

## Proper Fix Options

### Option 1: Global Node Counter (Recommended)

Add a global/shared node counter that persists across all models in a solve:

```cpp
class CbcModel {
  // ...
  int *sharedNodeCounter_;  // Points to shared counter
  int localNodeOffset_;     // This model's starting offset
  
  bool isNodeLimitReached() const {
    if (sharedNodeCounter_)
      return (*sharedNodeCounter_ >= intParam_[CbcMaxNumNode]);
    else
      return numberNodes_ >= intParam_[CbcMaxNumNode];
  }
};
```

When creating preprocessed model:
```cpp
newModel->sharedNodeCounter_ = &globalNodeCounter;
newModel->localNodeOffset_ = globalNodeCounter;
```

**Pros:** Clean, works for all restart scenarios  
**Cons:** Requires threading through model creation, copy constructors, and restart logic

### Option 2: Pass Remaining Budget

When creating a preprocessed or restarted model, calculate and pass the remaining node budget:

```cpp
// In integerPresolve()
int nodesAlreadyProcessed = numberNodes_;
int remainingBudget = intParam_[CbcMaxNumNode] - nodesAlreadyProcessed;
newModel->setIntParam(CbcModel::CbcMaxNumNode, 
                      std::max(0, remainingBudget));
```

**Pros:** Simpler, doesn't require new fields  
**Cons:** Breaks down with multiple restarts; needs to track total across all models

### Option 3: Accumulate Instead of Replace

When copying results back from preprocessed model:

```cpp
// In originalModel() - line 19698
numberNodes_ += presolvedModel->numberNodes_;  // ADD instead of REPLACE
```

**Pros:** Minimal change, correct accounting  
**Cons:** Doesn't actually stop the preprocessed model early; just reports correctly after the fact

### Option 4: Disable Preprocessing When Node Limit Used

Document that node limits don't work with preprocessing and auto-disable it:

```cpp
if (intParam_[CbcMaxNumNode] != COIN_INT_MAX) {
  // Disable integer preprocessing when node limit is set
  // (otherwise limit will be ignored)
  preprocessingMode = 0;
}
```

**Pros:** Simple workaround  
**Cons:** Defeats the purpose of preprocessing; may hurt performance

## Recommended Approach

Implement **Option 1 (Global Node Counter)** as it properly handles all cases:

1. Add `CbcSolveContext` struct to track global solve state:
   ```cpp
   struct CbcSolveContext {
     int globalNodeCounter;
     int maxNodes;
     // ... other global state
   };
   ```

2. Pass context pointer through model creation chain:
   ```cpp
   CbcModel *integerPresolve(bool weak, CbcSolveContext *ctx = nullptr);
   ```

3. Check global counter in `isNodeLimitReached()`:
   ```cpp
   bool CbcModel::isNodeLimitReached() const {
     if (solveContext_)
       return solveContext_->globalNodeCounter >= solveContext_->maxNodes;
     else
       return numberNodes_ >= intParam_[CbcMaxNumNode];
   }
   ```

4. Increment global counter alongside local counter throughout B&B code

## Workaround for Users

Until fixed, users should:

1. **Use time limits instead:** `-sec 300` or `Cbc_setMaximumSeconds(m, 300)`
2. **Disable preprocessing if node limit critical:** `-preprocess off -maxNodes 1000`
3. **Accept approximate node limits:** Understand reported counts may exceed limit by 10-100x

## Test Case

Added to `test/CInterfaceTest_bpc.c`:

```c
// Diagnostic: warn if node limit exceeded by 2x
int nodes_processed = Cbc_getNodeCount(m);
if (nodes_processed > tc->max_nodes * 2) {
  printf("    WARNING: processed %d nodes (limit was %d)\n",
         nodes_processed, tc->max_nodes);
}
```

BPC n=30 test now uses `INT_MAX` node limit and relies on 180s time limit.

## Related Files

- `src/CbcModel.cpp`: Main B&B logic, preprocessing, restarts
  - Lines 19271-19298: `integerPresolve()`
  - Lines 19302-19709: `integerPresolveThisModel()` 
  - Line 19675: preprocessed model `branchAndBound()` call
  - Lines 19619-19703: `originalModel()` result copy-back
  - Line 7665: `isNodeLimitReached()` check
  - Line 19719: `stoppingCriterionReached()` 
  - Lines 20283-20440: `subBranchAndBound()` functions

- `src/CbcParameters.cpp`: Parameter handling
  - Lines 1000-1040: `synchronizeModel()` 
  - Lines 3387-3390: MAXNODES parameter definition

- `src/Cbc_C_Interface.cpp`: C API
  - Lines 4920-4923: `Cbc_setMaximumNodes()` 
  - Lines 2618-2619: Convert to `-maxN` parameter
  - Line 2723: `cbcSolver.run()` processes parameter queue

- `test/CInterfaceTest_bpc.c`: Test demonstrating bug

## References

- Commit `05ffc903`: Documents bug in BPC test
- Commit `83199878`: QAP/BPC MPS fixes (separate issue)
