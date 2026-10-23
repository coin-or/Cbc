# Analyze parameter tunnings for Solving the LP relaxation of MIPs in MIPster

## Intro

We executed experiments to check how different LP methods (primal / dual) simplex and their parameters (pertv, psi, idiot) can improve the performance of CLP. We are interested in getting configurations which can optimally solve the LPs, producing **correct results**, in **minimimal time**. To be able to compare any result with another, __we convert every result to time__. When wrong results are produced, i.e. infeasible solution or not optimal objective value, we just return a penalty (i.e. 2x the time limit) for the time that that execution took. This way, we can compare every result with another.

## Instance set

The set of instances that we execute experiments is located at  ~/inst/miplib/2017+spp/ . Each insance is store in a .mps.gz file. The features of each instance are available in tabular format in file features.csv in the instance folder. 

### K-Fold Partitions

To enable out-of-sample validation of parameter recommendations, the 380 instances have been split into **10 folds** for k-fold cross-validation. Partition files are stored in:

    ~/inst/miplib/2017+spp/partitions/

Each file `fold_NN.txt` (fold_00.txt … fold_09.txt) contains 38 instances (one name per line, no extension). The split was generated with a fixed random seed (42) by shuffling all instances and assigning them round-robin across the 10 folds. There is no overlap between folds and every instance appears in exactly one fold.

**Workflow for k-fold evaluation:**
- For each fold `k` (0–9): train/select parameters on the other 9 folds, validate on fold `k`
- This gives unbiased estimates of how well a selection strategy generalises to unseen instances
- Report both in-sample and out-of-sample (mean over 10 folds) SGM / solve-rate


## Experiments

A recent experiment is available at 

	~/experiments/cbc/lp_relax_2026_05_15_noblas/
	
	These experiments were executed with a timelimit of 10800 seconds. 
	
The script to analyze results used was:

```
python analyze_lp_params.py --dir /home/haroldo/experiments/cbc/lp_relax_2026_05_15_noblas/ --timelimit 10800
```
	
generating the file lp_analysis.txt in the experiment folder.

since even for the same instance and parameter setting small perturbations can introduce differences, we executed different random seeds for the same setting and instance. Average runtimes per instance and parameter are available in the file lp_avg_times.csv in the experiments folder. Full results (with all different seeds) are available in lp_results.csv . 

The default parameter setting used by cbc is the parameter "cbc_default".


## Analysis Scripts

All k-fold evaluation scripts live in `lp_tune_scripts/` (under the workspace root `/home/haroldo/dev/cbc/`).

| Script | Purpose |
|--------|---------|
| `kfold_best_single_param.py` | K-fold CV for the "best single parameter" strategy: selects the param with the lowest training-set SGM, evaluates on the test fold |
| `kfold_dtree.py` | K-fold CV for feature-based decision trees (depths 1, 2, 3): trains a DTree on each 9-fold subset, traverses the tree for test instances |

Both scripts read `lp_avg_times.csv` from the experiment directory (produced by `analyze_lp_params.py`) and fold files from `~/inst/miplib/2017+spp/partitions/`. They write a text report to the experiment directory and print it to stdout.

**Common options (both scripts):**

| Option | Default | Description |
|--------|---------|-------------|
| `--dir DIR` | *(required)* | Experiment directory containing `lp_avg_times.csv` |
| `--timelimit T` | 10800 | LP time limit in seconds |
| `--penalty-mult M` | 2.0 | Multiplier for penalty time on failures |
| `--shift S` | 1.0 | SGM shift (seconds) |
| `--baseline PARAM` | `cbc_default` | Baseline param for speedup comparison |
| `--partitions DIR` | `~/inst/miplib/2017+spp/partitions` | Fold files directory |
| `--out FILE` | `kfold_*.txt` in exp dir | Output report path |

**`kfold_dtree.py`-only options:**

| Option | Default | Description |
|--------|---------|-------------|
| `--features PATH` | `~/inst/miplib/2017+spp/features.csv` | Instance features CSV |
| `--depths LIST` | `1,2,3` | Comma-separated tree depths to evaluate |
| `--min-leaf N` | 20 | Minimum instances per leaf node |

**Quick re-run on current experiment:**
```sh
cd /home/haroldo/dev/cbc

python3 lp_tune_scripts/kfold_best_single_param.py \
  --dir ~/experiments/cbc/lp_relax_2026_05_15_noblas/ --timelimit 10800

python3 lp_tune_scripts/kfold_dtree.py \
  --dir ~/experiments/cbc/lp_relax_2026_05_15_noblas/ --timelimit 10800
```

Output files written to the experiment directory:
- `kfold_best_single.txt` — single-param k-fold report
- `kfold_dtree.txt` — decision tree k-fold report (all depths)


## Tuning Parameters

We would like to select the **best** parameter setting for solving the LP in these instances — from the options which return correct results (feasible, with the correct objective value), the fastest options.

We are considering two strategies:

### Strategy 1 — Best single parameter

Select one parameter setting that performs best overall. Evaluated with k-fold cross-validation using the script `lp_tune_scripts/kfold_best_single_param.py`:

```
python3 lp_tune_scripts/kfold_best_single_param.py \
  --dir ~/experiments/cbc/lp_relax_2026_05_15_noblas/ --timelimit 10800
```

**Key results** (out-of-sample, 10-fold, baseline = `cbc_default`, SGM shift = 1s):
- Best param selected in 7/10 folds: **`dual_pesteep_psineg1`**
- Pooled test-set SGM: **12.923s** (tree) vs **13.234s** (baseline) → **~2.4% faster**
- In-sample upper bound (all instances): **1.036x faster**
- Gain is real but modest; motivates per-instance selection.

### Strategy 2 — Feature-based decision tree

Use instance features to decide automatically which parameter to apply. Trees are built with the `fbps` greedy algorithm (minimises sum of solve times in each leaf). Evaluated with k-fold cross-validation using `lp_tune_scripts/kfold_dtree.py`:

```
python3 lp_tune_scripts/kfold_dtree.py \
  --dir ~/experiments/cbc/lp_relax_2026_05_15_noblas/ --timelimit 10800
```

**How the decision tree works:**
- `run_lp_dtree.py` (and `kfold_dtree.py`) use the **fbps** library (`fbps/fbps/dtree.py`).
- The tree greedily selects feature splits that minimise the total sum of `avg_wall_seconds` across training instances. Each leaf assigns the param with the lowest total time.
- To apply the tree to a new instance, traverse from the root using feature comparisons (`feat ≤ value` → left, else right) until a leaf is reached.
- `avg_wall_seconds` already encodes **penalty time** (2× timelimit) for failed/wrong-result runs, so the optimisation naturally penalises unreliable params.

**Key results** (out-of-sample, 10-fold, baseline = `cbc_default`, SGM shift = 1s):

| Depth | Pooled test SGM | Baseline SGM | Speedup |
|-------|----------------|--------------|---------|
| 1     | 13.327s        | 13.234s      | **0.993× (slightly slower)** |
| 2     | 12.685s        | 13.234s      | **1.043× faster** |
| 3     | 12.399s        | 13.234s      | **1.067× faster** |

In-sample reference (optimistic, no hold-out): depth-3 achieves 1.128× speedup.  
Depth-2 and depth-3 trees consistently improve over the default; depth-1 does not generalise (too few degrees of freedom).

**Depth-3 in-sample IF-THEN rules (informative for understanding the structure):**
```
IF rhsAvg ≤ 0.103                                              → dual_pesteep_psineg1
IF rhsAvg ≤ 0.103 AND percColsLess4Nz > 17.3                  → dual_pesteep_psineg1_pertv61
IF rhsAvg ≤ 0.103 AND percColsLess4Nz > 17.3 AND percColsLeast4Nz > 76.4  → primal_idiot10
IF rhsAvg > 0.103 AND rInvKnapsack ≤ 630 AND rPercVarBnd ≤ 8.9  → primal_idiot50
IF rhsAvg > 0.103 AND rInvKnapsack ≤ 630 AND rPercVarBnd > 8.9  → dual_pesteep
IF rhsAvg > 0.103 AND rInvKnapsack > 630                       → dual_pesteep_psineg1
```

