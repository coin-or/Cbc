---
inclusion: always
---

# MIPster — Experiments & Benchmarking

## Environment Variables (set in shell profile)

```sh
export MIPSTER_PREFIX="$HOME/prog/cbc"
export MIPSTER_INSTANCES="$HOME/inst"
export MIPSTER_EXPERIMENTS="$HOME/experiments/cbc"
```

## Instance Collections

| Path | Instances | Description |
|---|---|---|
| `$MIPSTER_INSTANCES/miplib/2017+spp` | 358 | **PREFERRED** — MIPLIB 2017 + set packing/partitioning/covering |
| `$MIPSTER_INSTANCES/super/` | 1339 | Large curated superset; grouped by problem class |
| `$MIPSTER_INSTANCES/miplib/2017/` | 251 | MIPLIB 2017 benchmark set |
| `$MIPSTER_INSTANCES/cbc-test-set/` | ~985 | CBC regression/test set (with `.bas` warm-start files) |

All instances are `.mps.gz`. **`solvable_instances.tsv`** lists instances solvable in <5 min — use for smoke tests.

## Running the Solver

```sh
$MIPSTER_PREFIX/bin/mipster problem.mps -solve
```

> ⚠️ **Parameter order matters:** all settings must be placed **before** `-solve` to take effect (except action parameters like `-solu`).

## Parallel Experiment Pattern

```sh
# Adjust -j to available cores (use $(nproc) locally; higher on many-core servers)
ls $MIPSTER_INSTANCES/miplib/2017/*.mps.gz | \
  parallel -j$(nproc) '$MIPSTER_PREFIX/bin/mipster {} -sec 300 -solve > {/.}.log 2>&1'
```

## Experiment Results Format

Results in `$MIPSTER_EXPERIMENTS/<experiment_name>/<condition>/`:
```
STATUS  obj  bound  gap%  time_seconds  nodes  condition_name
```

Pick quick instances by checking solve times in recent results:
```sh
awk '$5 < 10' $MIPSTER_EXPERIMENTS/<experiment_name>/C0_baseline/*.result
```

## Capturing Statistics to CSV

```sh
# Append one row per solve:
$MIPSTER_PREFIX/bin/mipster problem.mps -csvStatistics results.csv -solve -writeStatistics

# Safe for parallel jobs (all append to same file):
ls $MIPSTER_INSTANCES/miplib/2017/*.mps.gz | \
  parallel -j$(nproc) '$MIPSTER_PREFIX/bin/mipster {} -sec 300 -csvStatistics results.csv -solve -writeStatistics'
```

## LP Relaxation & Racing Portfolio Analysis

Scripts: `run_lp_experiments.sh`, `analyze_lp_params.py`, `analyze_racing_portfolios.py`, `make_lp_report.py`

```sh
EXP=$MIPSTER_EXPERIMENTS/lp_relax_YYYY_MM_DD
python3 analyze_lp_params.py --dir $EXP --timelimit 14400
python3 analyze_racing_portfolios.py --dir $EXP --timelimit 14400
python3 make_lp_report.py $EXP --timelimit 14400   # → lp_report.pdf
```

Full workflow: `doc/lp_racing_analysis_howto.md`

## Debugging Invalid Cuts

```sh
# Get a known good solution:
$MIPSTER_PREFIX/bin/mipster problem.mps -preprocess off -solve -solu /tmp/known.sol
# Find cuts that exclude the optimal:
$MIPSTER_PREFIX/bin/mipster problem.mps -debugCuts /tmp/known.sol -solve 2>&1 | grep "bad row"
```

Full guide: `doc/cut-debugging-howto.md`
