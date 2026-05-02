# RINS Heuristic — Complete Parameter Reference

> **Scope:** Every parameter, option, and hardcoded constant that influences the
> RINS (Relaxed Induced Neighborhood Search) heuristic in CBC, including the
> `smallBranchAndBound` sub-MIP solver it delegates to.
>
> **Default values** are those in effect when running `cbc problem.mps -solve`
> (i.e. after `doHeuristics()` in `CbcSolverHeuristics.cpp` applies its overrides).

---

## 1. Command-Line Parameter

| CLI flag | Default | Values | Source |
|---|---|---|---|
| `-rins` | `on` | `off`, `on`, `both`, `before`, `often` | `CbcParameters.cpp:665,2516` |

Effect of each mode:

| Mode | Enum | Installs at | fractionSmall | decayFactor |
|---|---|---|---|---|
| `off` | 0 | — | — | — |
| `on` | 1 | B&B phase 1 | 0.5 | 5.0 |
| `both` | 2 | B&B phase 1 or 2 | 0.5 | 5.0 |
| `before` | 3 | B&B phase 2 or 3 | 0.5 | 5.0 |
| `often` | 5 | B&B phase 5 or 6 | 0.6 | 1.5 |

Source: `CbcSolverHeuristics.cpp:1659–1670`

---

## 2. RINS-Specific Parameters

These are set in the `CbcHeuristicRINS` constructor and may be overridden by
`doHeuristics()` or the API.

### 2.1 Frequency Control

| Parameter | Class default | Effective default | Configurable | Source |
|---|---|---|---|---|
| `howOften_` | 100 | 100 | `setHowOften()` | `CbcHeuristicRINS.cpp:30` |
| `decayFactor_` | 0.5 | **5.0** | `setDecayFactor()` | RINS ctor:32, overridden at `CbcSolverHeuristics.cpp:1664` |
| `lastNode_` | −999999 | −999999 | `setLastNode()` | `CbcHeuristicRINS.cpp:31` |

**How frequency works:**

1. When `howOften_ == 100` (the default), a special path applies:
   - Minimum 12-node gap since last run (hardcoded).
   - Forced early-run windows at nodes 41–50 and 91–99 (hardcoded).
2. Starvation guard: if `howOften_ >= 100` and `2 × howOften_` nodes have
   elapsed since last run, force execution regardless of modulo check.
3. Main gate: `(nodeCount % howOften_) == 0` AND
   `(currentPassNumber <= 1 OR currentPassNumber == 999999)`.

**Decay mechanism** (checked every 10 attempts):

```
if (numberTries_ % 10 == 0 && numberSuccesses_ * 3 < numberTries_)
    howOften_ += howOften_ * decayFactor_;
```

With the effective `decayFactor_ = 5.0`, this is extremely aggressive:

| After N failed checks | howOften_ |
|---|---|
| 10 | 100 → 600 |
| 20 | 600 → 3,600 |
| 30 | 3,600 → 21,600 |

RINS effectively disables itself after ~30 unsuccessful attempts.

Source: `CbcHeuristicRINS.cpp:316–317`

### 2.2 Variable Fixing Strategy

| Parameter | Default | Configurable | Source |
|---|---|---|---|
| `shallowDepth_` (repurposed as fixing mode) | 0 | `setShallowDepth()` | `CbcHeuristicRINS.cpp:29` |

Fixing modes:

| Value | Strategy | Description |
|---|---|---|
| 0 | Fix all agreeing | Fix any integer where LP relaxation ≈ incumbent (default) |
| 1 | Fix at lower bound only | Fix only if agreed value equals original lower bound |
| 2 | Fix away from lower bound | Fix only if agreed value differs from original lower bound |
| 3 | Fix at LB and unused | Fix only if at lower bound AND `usedInSolution[col] == 0` |

Source: `CbcHeuristicRINS.cpp:236–258`

### 2.3 Continuous Variable Fixing

| Parameter | Default | Configurable | Source |
|---|---|---|---|
| `stateOfFixing_` | 0 | No (internal state machine) | `CbcHeuristicRINS.cpp:29` |

Continuous fixing is triggered when ALL of:
- `numberContinuous > 2 × numberIntegers` (model is continuous-heavy)
- Either first-time conditions met OR `stateOfFixing_ != 0`

When active, fixes continuous variables at their lower bound sorted by
largest positive reduced cost. The fraction fixed is `1/divisor` of free
columns.

| State | divisor | Meaning |
|---|---|---|
| `stateOfFixing_ == 0` | 4 | Default: fix 25% of eligible continuous vars |
| `stateOfFixing_ > 0` | stateOfFixing_ | Previous success with this divisor |
| `stateOfFixing_ < -1` | −stateOfFixing_ − 1 | Previous failure, retry with this divisor |

Source: `CbcHeuristicRINS.cpp:260–270`

### 2.4 Scheduling (whereFrom)

| Parameter | Default | Configurable | Source |
|---|---|---|---|
| `whereFrom_` | `1 + 8 + 255×256` = 65289 | `setWhereFrom()` | `CbcHeuristicRINS.cpp:36` |

Bit interpretation (low byte):

| Bit | Position | Set? | Meaning |
|---|---|---|---|
| 0 | Before cuts at root | ✓ | RINS runs here |
| 1 | During cuts at root | ✗ | |
| 2 | After cuts at root | ✗ | |
| 3 | After cuts at non-root nodes | ✓ | RINS runs here |
| 4 | During cuts at non-root nodes | ✗ | |

High byte (`255 × 256`): all bits set — also runs when a previous heuristic
just found a solution.

---

## 3. Inherited Base-Class Parameters

From `CbcHeuristic`, set in its constructor and potentially overridden.

| Parameter | Base default | Effective for RINS | Configurable | Source |
|---|---|---|---|---|
| `numberNodes_` | 200 | 200 | `setNumberNodes()` | `CbcHeuristic.cpp:100` |
| `fractionSmall_` | 1.0 | **0.5** | `setFractionSmall()` | Base:102, overridden `CbcSolverHeuristics.cpp:1663` |
| `when_` | 2 | 2 | `setWhen()` | `CbcHeuristic.cpp:99` |
| `switches_` | 0 | 0 (or `-heurOptions` value) | `setSwitches()` | `CbcHeuristic.cpp:107` |
| `feasibilityPumpOptions_` | −1 | −1 (off) | `setFeasibilityPumpOptions()` | `CbcHeuristic.cpp:101` |
| `minDistanceToRun_` | 1 | 1 | `setMinDistanceToRun()` | `CbcHeuristic.cpp:115` |

---

## 4. smallBranchAndBound Parameters

These control the sub-MIP solve that RINS delegates to after fixing variables.

### 4.1 Size Gating

| Parameter | Value | Configurable | Source |
|---|---|---|---|
| `fractionSmall` (local copy) | 0.5 (from `fractionSmall_`) | Via `fractionSmall_` | `CbcHeuristic.cpp:696` |
| Large-model threshold | `2×rows + cols > 40000` | No (hardcoded) | `CbcHeuristic.cpp:742` |
| Large-model reduction | Up to 30% reduction of fractionSmall | No (hardcoded) | `CbcHeuristic.cpp:744` |
| Minimum presolved size for rejection | 300 (`after > 300`) | No (hardcoded) | `CbcHeuristic.cpp:834` |

**Size ratio calculation** (`sizeRatio()`, `CbcHeuristic.cpp:653–688`):

The ratio is `valueNow / valueStart` with a penalty multiplier:

| Row reduction | Multiplier |
|---|---|
| > 20% fewer rows | 1.0× |
| > 10% fewer rows | 1.1× |
| Any fewer rows | 1.5× |
| No reduction | 2.0× |

For long-and-thin models (`rows×10 < cols`), rows are weighted more heavily
(10× or 200× instead of 2×).

If `ratio > fractionSmall`, the sub-MIP is rejected (returns −1). A
**second chance** mechanism then tries fixing binary variables that appeared
most frequently in solutions (`usedInSolution`), with a 30% random acceptance
rate when `maxUsed == 1`.

### 4.2 Self-Decay of fractionSmall_

After each sub-MIP, if total iterations exceeded `100 × (numberNodes + 10)`:

```cpp
fractionSmall_ *= 0.9;  // permanent shrink
```

This is a **permanent** reduction applied to the member variable, affecting
all future calls. Over many sub-MIPs, this can shrink fractionSmall_ to near
zero.

Source: `CbcHeuristic.cpp:1736–1740`

### 4.3 Sub-MIP Solver Configuration

| Parameter | Value | Configurable | Source |
|---|---|---|---|
| Max nodes | `numberNodes_` = 200 | Via `numberNodes_` | `CbcHeuristic.cpp:1232` |
| Iteration limit | `iterationMultiplier × (numberNodes + 10)` = 21,000 | Via event handler | `CbcHeuristic.cpp:1580` |
| `iterationMultiplier` | 100 | Via event handler only | `CbcHeuristic.cpp:699` |
| Max cut passes (root) | `min(20, \|parent value\|)` | No (hardcoded) | `CbcHeuristic.cpp:1244` |
| Max cut passes (tree) | `min(10, parent value)` | No (hardcoded) | `CbcHeuristic.cpp:1245` |
| Hot start iterations | 10 | No (hardcoded) | `CbcHeuristic.cpp:1243` |
| Strategy | `CbcStrategyDefaultSubTree(model_, 1, 5, 1, 0)` | No (hardcoded) | `CbcHeuristic.cpp:1241` |
| Strong branching | Off (`setStrongStrategy(0)`) | No (hardcoded) | `CbcHeuristic.cpp:1215` |
| Fast node depth | −1 (off) | No (hardcoded) | `CbcHeuristic.cpp:1213` |
| Conflict analysis | Off (bit cleared) | No (hardcoded) | `CbcHeuristic.cpp:1237` |

### 4.4 Preprocessing in Sub-MIP

| Phase | What | Options | Source |
|---|---|---|---|
| Fast MILP preprocessing | Singleton tightening + knapsack BT | MILPbt, 100 rounds | `CbcHeuristic.cpp:791` |
| OsiPresolve | LP presolve (size check) | `presolveActions = 1` | `CbcHeuristic.cpp:802` |
| LP solve | `solver->initialSolve()` | — | `CbcHeuristic.cpp:931` |
| CglPreProcess | MIP preprocessing | See below | `CbcHeuristic.cpp:1030` |

CglPreProcess options:

| Condition | Options | Passes |
|---|---|---|
| Normal | `setOptions(16)` — no dupcol | 2 |
| `moreSpecialOptions & 65536` | `setOptions(2+4+8+16)` — no cuts | 2 |
| `moreSpecialOptions2 & 16` | `setOptions(2+4+8+16)` — no cuts | 1 |
| Restart (`numberNodes < 0`) | `setOptions(16)` | 5 |

### 4.5 Event Handler Hook

Before any processing, the event handler can override:

| Field | Default | Can be changed to |
|---|---|---|
| `fractionSmall` | 0.5 | Any value |
| `numberNodes` | 200 | Any value |
| `iterationMultiplier` | 100 | Any value |
| `howOften_` | 100 | Any value |
| `maximumSolutions` | model's value | Any value |

Source: `CbcHeuristic.cpp:700–727`

---

## 5. Hardcoded Constants

| Constant | Value | Location | Purpose |
|---|---|---|---|
| Minimum node gap | 12 | `CbcHeuristicRINS.cpp:178` | Min nodes between runs when `howOften_==100` |
| Early-run window 1 | nodes 41–50 | `CbcHeuristicRINS.cpp:181` | Force run near node 50 |
| Early-run window 2 | nodes 91–99 | `CbcHeuristicRINS.cpp:181` | Force run near node 100 |
| Starvation multiplier | 2× | `CbcHeuristicRINS.cpp:185` | Force run after `2×howOften_` nodes |
| Agreement tolerance | 10 × primalTolerance | `CbcHeuristicRINS.cpp:210` | LP–incumbent agreement check |
| Minimum fix ratio | 20% | `CbcHeuristicRINS.cpp:232` | `5×nFix > numberIntegers` |
| Continuous-heavy threshold | 2× | `CbcHeuristicRINS.cpp:233` | `nContinuous > 2×nIntegers` |
| Continuous fix DJ threshold | 1e-6 | `CbcHeuristicRINS.cpp:253` | Min reduced cost to fix continuous var |
| Continuous fix LB tolerance | 1e-8 | `CbcHeuristicRINS.cpp:247` | At-lower-bound check |
| Default continuous divisor | 4 | `CbcHeuristicRINS.cpp:262` | Fix 25% of eligible continuous vars |
| Decay check interval | 10 tries | `CbcHeuristicRINS.cpp:316` | Check success rate every N tries |
| Decay success threshold | 33% | `CbcHeuristicRINS.cpp:316` | `successes×3 < tries` |
| Large-model threshold | 40,000 | `CbcHeuristic.cpp:742` | `2×rows+cols` |
| Large-model max reduction | 30% | `CbcHeuristic.cpp:744` | Max fractionSmall reduction |
| Size rejection min | 300 | `CbcHeuristic.cpp:834` | `after > 300` to reject |
| Random skip (usedInSolution) | 0.3 | `CbcHeuristic.cpp:862` | 70% chance to skip fixing when maxUsed==1 |
| Iteration multiplier | 100 | `CbcHeuristic.cpp:699` | Sub-MIP iteration limit factor |
| fractionSmall self-decay | 0.9× | `CbcHeuristic.cpp:1740` | Permanent shrink on excess iterations |
| Sub-MIP max cut passes (root) | min(20, parent) | `CbcHeuristic.cpp:1244` | Cut generation limit |
| Sub-MIP max cut passes (tree) | min(10, parent) | `CbcHeuristic.cpp:1245` | Cut generation limit |
| Sub-MIP hot start iters | 10 | `CbcHeuristic.cpp:1243` | Hot start iteration limit |
| CglPreProcess passes (normal) | 2 | `CbcHeuristic.cpp:941` | Preprocessing depth |
| CglPreProcess passes (quick) | 1 | `CbcHeuristic.cpp:944` | Quick preprocessing |
| CglPreProcess passes (restart) | 5 | `CbcHeuristic.cpp:947` | Restart preprocessing |

---

## 6. Tuning Recommendations

### 6.1 High-Impact Parameters (already configurable)

| Parameter | Why tune | Suggestion |
|---|---|---|
| **`decayFactor_`** | At 5.0, RINS kills itself after ~30 failures. For hard instances where RINS needs many attempts, this is too aggressive. | Try 1.0–2.0 for hard instances. The `often` mode already uses 1.5. |
| **`fractionSmall_`** | At 0.5, many sub-MIPs are rejected as "too large". Combined with the permanent 0.9× self-decay, RINS becomes increasingly conservative. | Try 0.6–0.8 for instances where heuristics are important. |
| **`numberNodes_`** | 200 nodes may be too few for hard sub-MIPs, or wasteful for easy ones. | Adaptive: start at 200, increase if sub-MIPs are finding solutions, decrease if not. |
| **`howOften_`** | 100 is reasonable but the special-case logic for `howOften_==100` (early windows, 12-node gap) is fragile. | Consider making the early-run windows configurable. |

### 6.2 Hardcoded Values Worth Exposing as Parameters

| Constant | Current value | Why expose | Suggested parameter |
|---|---|---|---|
| **Minimum fix ratio** | 20% (`5×nFix > nIntegers`) | On some instances, fixing 15% of integers still produces useful sub-MIPs. On others, even 30% is too few. | `-rinsMinFixRatio` (0.0–1.0) |
| **Sub-MIP node limit** | 200 | Already configurable via API but not CLI. This is the single most impactful sub-MIP parameter. | `-rinsNodes` (integer) |
| **Iteration multiplier** | 100 | Controls the hard iteration cap. Only accessible via event handler. | `-rinsIterMult` (integer) |
| **fractionSmall self-decay** | 0.9× | The permanent shrink is a blunt instrument. On long runs, it can make fractionSmall_ approach zero, effectively disabling all sub-MIP heuristics. | `-rinsDecayFraction` (0.0–1.0), or disable entirely |
| **Continuous variable divisor** | 4 | Controls how aggressively continuous variables are fixed. | `-rinsContDivisor` (integer) |

### 6.3 Structural Improvements to Consider

**Decay mechanism redesign.** The current decay has two problems:
1. `decayFactor_ = 5.0` is too aggressive — RINS disables itself quickly.
2. The `fractionSmall_ *= 0.9` permanent shrink compounds across all
   heuristics sharing `smallBranchAndBound`, not just RINS.

A better approach might be:
- Use a **bounded decay** with a floor: `howOften_ = min(howOften_ + howOften_ * decay, maxHowOften)`.
- Replace the permanent `fractionSmall_` shrink with a **per-call** adjustment
  that doesn't persist.
- Add a **recovery mechanism**: if RINS finds a solution after a long gap,
  reset `howOften_` closer to its initial value.

**Fixing strategy selection.** The `shallowDepth_` repurposing (modes 0–3) is
not exposed as a CLI parameter. Mode 1 (fix only at lower bound) could be
valuable for set-packing instances where most variables are 0 in the
incumbent. Mode 3 (fix at LB and unused) could help diversify the search.

**Sub-MIP preprocessing.** The CglPreProcess options in `smallBranchAndBound`
are controlled by `moreSpecialOptions` bits (65536, moreSpecialOptions2 bit 4)
which are inherited from the parent model. There is no way to independently
control sub-MIP preprocessing aggressiveness. A dedicated parameter could
allow:
- Skip CglPreProcess entirely (rely only on fast preprocessing + OsiPresolve)
- Use quick mode (1 pass, no cuts)
- Use normal mode (2 passes)

**Adaptive node budget.** Instead of a fixed 200-node budget, the sub-MIP
could adapt based on:
- How much the sub-MIP was reduced (more reduction → more nodes)
- Whether previous sub-MIPs found solutions (success → more nodes)
- Remaining time budget

---

## 7. Related Parameters (Not RINS-Specific)

| Parameter | CLI | Default | Effect on RINS |
|---|---|---|---|
| `-heurOptions` | integer | 0 | Sets `switches_` on ALL heuristics after setup |
| `-smallBab` | double | 0.5 | Sets `fractionSmall` for FPump only (not RINS) |
| `-fastPreProcessLevel` | keyword | `milpbt` | Controls root fast preprocessing; sub-MIP fast preprocessing is hardcoded to MILPbt with 100 rounds |
