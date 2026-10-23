# Teaching a MIP Solver to Choose Its Own LP Method

*One of the most impactful — and most overlooked — decisions in branch-and-bound is which
LP algorithm to use at each node. We trained a machine-learning model to make that choice
automatically, and embedded it directly into [MIPster](https://github.com/h-g-s/mipster),
our open-source MIP solver forked from [COIN-OR CBC](https://github.com/coin-or/Cbc).*

---

## The Problem

Every LP relaxation in B&B is solved with simplex — either dual or primal, with various
warm-start and perturbation strategies. CBC defaults to dual simplex, which works well on
average, but some instances are dramatically faster with primal simplex plus an
[idiot crash](https://github.com/coin-or/Clp/blob/master/src/ClpSimplex.cpp) warm-start,
or with a steepest-edge pivot rule, or a non-default perturbation value.
The trouble is: knowing *which* method to pick requires understanding the problem structure
before solving it.

## Feature Extraction

We used [OsiFeatures](https://github.com/h-g-s/mipster) — a library of 207 structural
features extracted directly from the LP relaxation matrix — as input to our model.
Features include constraint-type fractions (set packing, partitioning, covering, flow),
variable type fractions (binary, integer, continuous), density, objective coefficient
statistics, and row/column NZ distributions.
Extraction adds negligible overhead (sub-millisecond) even for large instances.

## ML Modelling

We framed the problem as **classification**: given 207 features, predict which LP
parameter configuration yields the fastest solve on this instance. We evaluated:

| Model | k-fold LP speedup vs default |
|---|---|
| Decision tree depth 3 | 1.18× |
| [Random Forest](https://scikit-learn.org/stable/modules/ensemble.html#forests-of-randomized-trees) n=150, depth=8 | **1.33×** |
| [XGBoost](https://xgboost.readthedocs.io/) (classification) | 1.32× |
| [LightGBM](https://lightgbm.readthedocs.io/) (classification) | 1.32× |
| Multi-output regression (predict runtime per param) | 1.28× |

Random Forest won on the k-fold estimate and — crucially — generates compact
embeddable C++ code via [m2cgen](https://github.com/BayesWitnesses/m2cgen)
(~34K lines, no runtime dependencies). XGBoost 3.x is not yet supported by m2cgen 0.10.
LightGBM generates ~280K lines — too large to compile reliably.

## Reducing the Label Space via Win-Profile Clustering

Naively expanding the parameter set (we eventually tested 70 configurations) actually
*hurt* classification accuracy: 70 classes × 380 instances gives roughly 5 examples per
class on average — not enough for shallow trees to learn reliably.

The insight is that many parameters are **redundant specialists**: they win on almost the
same subset of instances, just with slightly different timings. We identified these groups
using [Jaccard distance](https://en.wikipedia.org/wiki/Jaccard_index) on binary
*win-profile vectors* — for each parameter, a bit-vector of whether it was within 5% of
the best time on each instance. Hierarchical
[Ward linkage](https://en.wikipedia.org/wiki/Ward%27s_method) clustering on these vectors
groups parameters that specialise on the same instances. One representative (most cluster
wins, geo-mean as tiebreaker) is selected per cluster.

With **k=12 clusters** and a 5% competitiveness threshold:

| Model | Params / classes | k-fold speedup |
|---|---|---|
| RF depth=8 | 24 (original) | 1.33× |
| RF depth=10 | 64 (extended) | 1.30× |
| RF depth=8 | **12 (clustered)** | **1.36×** |

Fewer, well-separated classes means more training examples per class and better
generalisation — depth=8 is sufficient again, and the speedup exceeds the original set.
The 12 representatives cover all the key solver modes: dual simplex with various
perturbation strengths and pivot strategies, and primal+idiot at different iteration
budgets (10 through 500) for increasingly large sparse problems.

## Important Features

The top predictors of which LP method wins (by Random Forest feature importance and
Pearson correlation with speedup) were:

- **Objective coefficient scale** (`objAvg`, `objStdDev`, `objRatioLSA`) — the single
  strongest signal; large heterogeneous objectives indicate transport/logistics structure.
- **Covering row fraction** (`rPercCovering`) — set covering rows strongly predict
  primal+idiot dominance.
- **Matrix density** — sparse problems (< 1%) overwhelmingly favour primal+idiot crash.
- **Sparse row fraction** (`percRowsLeast256Nz`) — large problems with thin rows.
- **Binary flow rows** (`rFlowBin`) — flow-type structure also pulls toward primal.

## Which Instances Benefit Most?

We benchmarked on the [MIPLIB 2017](https://miplib.zib.de/) + set-packing extension
(380 instances). Speedup on the *LP relaxation* (root node only):

| Tier | N | Typical structure | LP speedup |
|---|---|---|---|
| ★★★ > 5× | 47 | Sparse (0.4%), ≥30% covering rows, 60K+ vars | up to 50× |
| ★★  2–5× | 64 | Sparse (1.8%), 13% covering, large | 2–5× |
| ★   1.1–2× | 199 | Medium density (4%), mixed types | modest |
| ◻  ≤ 1.1× | 63 | Dense (6%+), smaller | no benefit |

The pattern is striking: **large, very sparse, set-covering/packing problems** are where
the model shines. These are exactly the problem classes for which the idiot crash was
designed — it finds a near-feasible primal point in a handful of passes over a sparse
matrix, something dual simplex cannot exploit. Our pet instance, `brazil3`
(density 0.04%, 14K×24K), went from 11.9 s to 2.0 s — a **6× speedup** on the LP alone.

## Integration

The trained forest is exported to a single C++ file with no dependencies.
Activating it is one parameter: `-lpMethod=auto` (or `Cbc_setLpMethod(model, LPM_Auto)`
via the C interface). The solver extracts features, scores 12 configurations in ~0.1 ms,
and proceeds with the winner — fully transparent in the solve log:

```
LP auto: recommended primal_idiot50
✔ LP Optimal — Obj: … Iters: 334 Time: 0.14s
```

---

*All code is in the MIPster repository. The training script lives in
`lp_tune_scripts/build_lp_param_scorer.py`; the generated scorer in
`Cbc/src/CbcLpParamScorer.{hpp,cpp}`.*
