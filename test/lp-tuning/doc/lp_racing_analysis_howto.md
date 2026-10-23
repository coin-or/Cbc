# LP Relaxation & Racing Portfolio Analysis — How-To

This document describes the full workflow for running LP relaxation parameter
tuning experiments and analysing the results, including for the **racingLP**
feature in CBC (parallel opportunistic LP solving).

---

## Overview of scripts

| Script | Purpose |
|---|---|
| `run_lp_experiments.sh` | Run LP relaxation solves across params × instances × seeds in parallel |
| `run_one_lp.sh` | Worker called by the above via GNU parallel (not invoked directly) |
| `lp_params.txt` | Parameter configurations to test (one `tag|cbc_params` per line) |
| `analyze_lp_params.py` | Individual param analysis: SGM rankings, head-to-head, errors — produces `lp_avg_times.csv` |
| `analyze_racing_portfolios.py` | **Racing/portfolio analysis**: which K params to race for best expected time |
| `run_lp_dtree.py` | **Decision tree**: builds fbps decision tree from features + lp_avg_times.csv; outputs `lp_dtree.txt` + `lp_dtree.svg` |
| `make_lp_report.py` | **PDF report**: visual summary of all discoveries (6 pages, charts + decision tree + key findings) |
| `summarize_lp_results.py` | Fake-optimal detection and per-instance validation |

---

## Step 1 — Run the experiment

```sh
./run_lp_experiments.sh \
  --bin ~/prog/cbc/bin/cbc \
  --params lp_params.txt \
  --instances /home/haroldo/inst/miplib/2017+spp \
  --timelimit 14400 \
  --seeds 1,2,3 \
  --parallel $(( $(nproc) - 2 )) \
  --outdir /home/haroldo/experiments/cbc/lp_relax_YYYY_MM_DD/
```

Produces per-job `.sol`, `.bas`, `.log`, `.result` files and a consolidated
`lp_results.csv` in the output directory. The experiment is **resumable** —
re-running skips jobs whose `.result` file already exists.

---

## Step 2 — Individual param analysis

```sh
python3 analyze_lp_params.py \
  --dir /home/haroldo/experiments/cbc/lp_relax_YYYY_MM_DD/ \
  --timelimit 14400 \
  --baseline dual_default
```

Outputs:
- `lp_analysis.txt` — full console report (SGM ranking, head-to-head, crashes, timeouts)
- `lp_avg_times.csv` — per-(instance, param) averaged times; **required by the racing script**

### What gets penalised
- **Timeout / killed**: wall time replaced by `timelimit × penalty_mult` (default 2×)
- **Crash** (SIGSEGV/SIGABRT): same penalty
- **Wrong result**: status=OPTIMAL but objective worse than best-known by `> --wrong-tol`
  (default 1e-4 relative). These are also penalised — the solver lied about optimality.
- **TRIVIALLY_OPTIMAL** (STATUS=ERROR but check passes): treated as genuinely solved at
  actual wall time; the solver found the answer without any LP iterations.

Best-known objective per instance is the minimum `obj_from_check` across all runs
that passed the feasibility check (`check_result` starts with `yes`).

---

## Step 3b — Decision tree (feature-based parameter selection)

The fbps decision tree partitions instances by their structural features (nr. of
columns, non-zeros, constraint types, etc.) and recommends the fastest LP parameter
setting for each partition — building the rationale for instance-adaptive racingLP.

```sh
python3 Cbc/test/lp-tuning/run_lp_dtree.py \
  --dir /home/haroldo/experiments/cbc/lp_relax_YYYY_MM_DD/ \
  --timelimit 7200 \
  --depth 3 \
  --min-leaf 20 \
  --features ~/inst/miplib/2017+spp/features.csv
```

Outputs:
- `lp_dtree.txt` — text report: for each leaf the recommended param, instance count, speedup
- `lp_dtree.svg` — visual tree diagram

The decision tree is also embedded in the PDF report (page 5) automatically.
Use `--no-dtree` to skip it, or `--dtree-depth`/`--dtree-min-leaf` to tune it:

```sh
python3 make_lp_report.py $EXP --timelimit 7200 --dtree-depth 4 --dtree-min-leaf 15
python3 make_lp_report.py $EXP --timelimit 7200 --no-dtree  # skip dtree page
```

### Interpreting the decision tree

- **Root split feature**: the single most discriminating instance property.
  A split on `colNzMax ≤ 7` means dense vs sparse column structure drives method choice.
- **Leaf speedup**: how much faster the recommended param is vs the global best param
  (i.e. within that leaf sub-population).
- **Leaf size**: leaves with < `--dtree-min-leaf` instances are not split further;
  small leaves may overfit — check with `--min-leaf 30` if suspicious.
- **Overall speedup**: geometric mean improvement over `dual_default` baseline.

---

## Step 3 — Racing / portfolio analysis

> **racingLP** = CBC feature where K CLP solves run in parallel with different
> parameter settings; whichever finishes first provides the LP optimal, and the
> others are cancelled.  For a portfolio of size K the effective solve time per
> instance is `min(t_1, ..., t_K)`.

```sh
python3 analyze_racing_portfolios.py \
  --dir /home/haroldo/experiments/cbc/lp_relax_YYYY_MM_DD/ \
  --timelimit 14400 \
  --max-k 8 \
  --exhaustive-k 4 \
  --baseline dual_default
```

**Requires `lp_avg_times.csv`** — run `analyze_lp_params.py` first.

Outputs:
- `lp_racing_analysis.txt` — full report (portfolio curve, exhaustive optimal, marginal
  contributions, per-param uniqueness, instance coverage map)
- `lp_racing_curve.csv` — greedy and exhaustive SGM values per K; easy to plot

### Key tables in the report

| Table | What it shows |
|---|---|
| Portfolio construction curve | SGM, solve%, and `ΔStep%` for K=1..max_k; shows where returns diminish |
| Exhaustive optimal portfolios | True global optimum for K ≤ `--exhaustive-k`; greedy quality check |
| Marginal contribution | How many instances each newly added param "wins" in the portfolio |
| Per-param uniqueness | Params that are the **only** solver on some instances (high portfolio value) |
| Racing contributions K=1..5 | How many instances each param covers at each portfolio size |
| Most improved instances | Instances where racing (K=2) helps most vs best single param |
| Hard instances | Instances no param solves within timelimit |

### Interpreting results

- **Big ΔStep% drop at K=2**: typically the biggest gain — a primal and a dual
  variant together cover complementary instance types.
- **ΔStep% < 5% by K=4–5**: diminishing returns; racingLP with 3–4 runners is
  usually sufficient.
- **Greedy quality**: if greedy SGM is within ~2% of exhaustive optimal it's fine to
  trust greedy for larger K where exhaustive is infeasible.
- **Uniqueness table**: params with unique wins are irreplaceable; those with zero
  unique wins can be substituted freely.

### Options

| Option | Default | Description |
|---|---|---|
| `--timelimit` | 14400 | Must match the experiment's `-sec` value |
| `--penalty-mult` | 2.0 | Penalty = timelimit × mult for failures |
| `--shift` | 1.0 | SGM shift in seconds |
| `--max-k` | 8 | Max portfolio size for greedy search |
| `--exhaustive-k` | 4 | Max K for exhaustive optimal (C(28,4) ≈ 20k, fast) |
| `--baseline` | dual_default | Reference single param for comparison columns |

---

## Step 4 — Validation (optional)

```sh
python3 summarize_lp_results.py \
  /home/haroldo/experiments/cbc/lp_relax_YYYY_MM_DD/ \
  --verbose
```

Produces `lp_summary.csv` annotating each run with `OK / FAKE_OPTIMAL / TIMEOUT`.
Useful if you suspect numerical issues with specific param/instance combinations.

---

## Experiment result location

Past LP relaxation experiments: `/home/haroldo/experiments/cbc/`

| Experiment | Notes |
|---|---|
| `lp_relax_2026_05_12_noblas` | 358 instances (miplib/2017+spp), 24 param configs (guess_* removed), seeds 1–3; timelimit 7200s; no OpenBLAS |

---

## Quick re-run recipe

```sh
EXP=/home/haroldo/experiments/cbc/lp_relax_YYYY_MM_DD
TL=7200  # match the --sec value used in the experiment

# 1. Individual analysis (produces lp_avg_times.csv):
python3 Cbc/test/lp-tuning/analyze_lp_params.py \
  --dir $EXP --timelimit $TL --baseline dual_default

# 2. Racing portfolio analysis:
python3 Cbc/test/lp-tuning/analyze_racing_portfolios.py \
  --dir $EXP --timelimit $TL --max-k 8 --exhaustive-k 4

# 3. Decision tree (standalone, also embedded in PDF):
python3 Cbc/test/lp-tuning/run_lp_dtree.py \
  --dir $EXP --timelimit $TL --depth 3 --min-leaf 20 \
  --features ~/inst/miplib/2017+spp/features.csv

# 4. PDF report (6 pages incl. decision tree, requires steps 1-3):
python3 Cbc/test/lp-tuning/make_lp_report.py $EXP --timelimit $TL

# View reports:
less $EXP/lp_analysis.txt
less $EXP/lp_racing_analysis.txt
less $EXP/lp_dtree.txt
xdg-open $EXP/lp_report.pdf
```
