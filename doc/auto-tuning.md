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
