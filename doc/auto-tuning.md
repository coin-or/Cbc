# Problem-class auto-tuning of cut generators

CBC classifies the (preprocessed) model from its `OsiFeatures` and, for
recognised 0-1 structures, switches off cut generators that only waste time and
keeps the ones that work. Nothing is forced on the user: a parameter is changed
only while it is `auto`.

## The `auto` state

Each cut generator parameter below accepts the value `auto` and starts as
`auto`:

`cliqueCuts`, `oddwheelCuts`, `probingCuts`, `gomoryCuts`, `twomirCuts`,
`mixedIntegerRoundingCuts`, `GMICuts`, `knapsackCuts`, `flowCoverCuts`,
`liftAndProjectCuts`, `reduceAndSplitCuts`, `reduce2AndSplitCuts`,
`zeroHalfCuts`.

* `auto` = built-in default, unless `autoTuneCuts` is on (the default) and the
  class below recommends something else.
* Any explicit value, including one equal to the built-in default
  (`-gomoryCuts ifmove`), is a user setting and is **never** overridden.
* `-gomoryCuts auto` goes back to `auto`.
* `-autoTuneCuts off` disables all automatic changes (everything behaves as the
  built-in defaults).

The solver log shows what happened, e.g.

```
AutoTuneCuts: SetPartitioningPacking (100% binary variables, 0-1 coefficients, ...)
  -> Set Partitioning & Packing: active: BKClique, OddWheel, Probing; disabled: ...
  auto-set: gomoryCuts=off ... cliqueCuts=ifmove oddwheelCuts=ifmove probingCuts=ifmove
  user-set (kept): <parameters you set explicitly>
```

## Classes and recommendations

Classes are tested in this order (all require 100% binary variables).

| Class (`-problemClass`) | Rows | Disabled | Enabled (`ifmove`) |
|---|---|---|---|
| `partitioning` | all `sum x = 1` | Gomory, TwoMir, MIR, GMI, RedSplit/2, LandP, Knapsack, Flow, ZeroHalf | Clique, OddWheel, Probing |
| `packing` | all `sum x <= 1` | Knapsack, Flow, LandP, RedSplit/2, GMI | Clique, OddWheel, Gomory, TwoMir, Probing |
| `partpack` | `= 1` / `<= 1`, 0-1 coefficients | as `partitioning` | as `partitioning` |
| `covering` | all `sum x >= 1` | Clique, OddWheel, Knapsack, Flow, LandP, RedSplit/2, GMI | ZeroHalf, Gomory, TwoMir, Probing |
| `mixedspp` / `mostlyspp` | 100% / >= 85% of rows partitioning, packing or covering | Knapsack, Flow, LandP, RedSplit/2 | ZeroHalf, Clique, OddWheel, Gomory, TwoMir, MIR, Probing |
| `knapsack` | >= 85% 0-1 knapsack rows | none | Knapsack, Gomory, TwoMir, Clique, Probing |
| `general01`, `general` | anything else / continuous or general integer variables | no change | no change |

The recommendations come from root-fixture experiments (`mip-root-replay`):
the disabled generators produced no cuts on those structures, while the
enabled ones carry the bound improvement. On a 46-instance benchmark cut time
dropped 4.2x on set partitioning/packing, with equal or better root bounds on
nearly all instances.

## Related parameters

* `-problemClass <class>`: force a class instead of detecting it (`auto`).
* `-cutPassRatio <r>`: multiply the resolved number of root cut passes by `r`
  (e.g. `1.5`, `0.9`).
* `-minDropRatio <r>`: multiply the resolved `minimumDrop` by `r`.

The ratios are relative to whatever value would otherwise be used (including
an `auto` `cutPass`), so they compose with automatic settings instead of
replacing them.

## Extending to other classes

Add an enumerator to `CbcProblemClassifier::ProblemClass`, detect it in
`classify()`, and add a `case` to `getCutTuning()` listing `disableParams` /
`enableParams`. New tunable parameters must be added to
`CbcProblemClassifier::autoCutParams()` so they accept `auto`.

# Heuristic invocation control

How often a tree heuristic runs is decided in one place,
`CbcHeuristicControl` (`src/CbcHeuristicControl.hpp`), the same for RINS, VND,
DINS, combine, combine2 (crossover) and proximity search. For each of them
there are eight parameters, `<h>` being `rins`, `vnd`, `dins`, `combine`,
`combine2` or `proximity`:

| Parameter | Meaning |
|---|---|
| `<h>Schedule` | `auto`, `legacy`, `depth`, `nodes`, `noimprove` (every K nodes and straight after each new incumbent) |
| `<h>Every` | K of the schedule |
| `<h>Decay` | back-off: K is multiplied by `1 + decay` ... |
| `<h>Window` | ... every `window` runs, if fewer than ... |
| `<h>MinSuccess` | ... this fraction of the runs improved the incumbent |
| `<h>MaxFail` | pause after this many failed runs on the same incumbent (resumes when it changes; 0 = never) |
| `<h>TimeShare` | do not start a run while the heuristic has used more than this fraction of the elapsed time (0 = no limit) |
| `<h>Nodes` | node limit of the sub-MIP |

All accept `auto`, which keeps the built-in value (RINS/VND: `noimprove`,
K=100, window 10, min success 1/3, decay 5 for RINS/DINS with `on` and 1.5
otherwise, none for VND). For the heuristics triggered by a new incumbent
(combine, combine2, proximity) a non-legacy schedule means "at least K nodes
since the last run".

Class-specific `auto` (with `-autoTuneCuts on`): on a `PureSetCovering`
problem VND's `auto` is `noimprove`, K=500, decay 0, max failures 3. Screening
on 102 0-1 instances showed VND rarely pays there, and this gave PDI
-0.0088 on that class and kept the optimal count (48). The same setting on
all classes lost `e201600` (a MostlySPP instance that needs frequent VND), so
other classes keep the default above. Explicit parameters override this.

Back-off previously changed an interval that the default `noimprove` schedule
never read, so it had no effect; it now does. With `-log 2` a per-heuristic
summary (runs, improvements, seconds, back-offs, skips, final policy) is
printed at the end, and `CBC_HEUR_TRACE=1` traces every decision.
`Cbc/test/test-heuristic-control` checks that the parameters control the
number of runs. The summary and the controls apply to the model's own
heuristics, so use the default `-threads 0` when measuring them.

To put another heuristic under the same control, add a line to
`CBC_HEUR_CONTROL_LIST` in `CbcParam.hpp`, route its gating through
`shouldRunBySchedule()` (periodic) or `controlAllows()` (incumbent-triggered)
and bracket its run with `CbcHeuristic::RunGuard`.
