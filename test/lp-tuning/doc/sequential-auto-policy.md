# Sequential root-LP auto policy

`lpMethod auto` uses the following ordered policy when fewer than two LP
threads are available. `lpMethod recommend` explicitly requests this sequential
policy regardless of thread count. With at least two threads, `auto` still
uses the existing parallel racing portfolio; this change does not alter it.
Explicit `dual`, `primal`, `barrier`, or simplex actions bypass the policy.

The first matching rule wins. Features are computed on the strengthened root
LP, before Clp's internal LP presolve. Thresholds retain full discovery precision.

| Rule | Conditions (AND) | Root configuration | Historical support | Big wins |
|---|---|---|---:|---:|
| 1 | `rNzPercRowsCovering > 50.7092105993763`; `colsLeast4Nz > 1056.5`; `percRowsLeast4096Nz <= 0.0445202274795385` | primal, idiot60, pertValue100, psi0 | 18 | 18 |
| 2 | `rNzPercRowsHubImplication > 55.30788084035105`; `percBin > 99.04976740485506`; `cols > 810` | primal, idiot50, pertValue100, psi0 | 15 | 10 |
| 3 | `rNzpercRowsPartPacking > 99.82268606961816`; `rowsLeast8Nz > 1135.5` | primal, idiot60, pertValue100, psi0 | 15 | 9 |
| 4 | `objStdDev <= 0.001552965637549735`; `rNzrowsPartPacking > 2026` | primal, idiot50, pertValue50 | 15 | 7 |
| Otherwise | No match | dual, pertValue72, idiot0, psi0 | 339 | - |

Sprint is disabled in all five configurations. Rule 4 was measured under the
legacy tag `racing_primal_idiot50`, but it is a **single sequential primal
solve**, not racing. Its production tag is `primal_idiot50_pertv50`.
Its omitted psi uses the measured CLI default (-0.5); like psi0, this does
not install Positive Edge pivots. Other settings, including scaling and
presolve, retain their existing behavior. Feature indices use Osi's symbolic
enumerators, not positions in the retired 207-feature RF model. Each rule's
original discovery-time NaN imputations are retained in the scorer.

## Root settings versus reoptimization

These are one-shot initial-LP choices, not a policy for later LP resolves.
The existing capture/restore path restores pivot rules, scaling and
perturbation after the initial solve. Idiot and sprint choices are local solve
options. Then the independent `reoptMethod`, `reoptPertValue`, and
`reoptPresolve` controls govern later resolves. Their defaults and explicit
user overrides are unchanged: root pertValue72/100 must not leak into later
cut, heuristic, or search-node LPs. This also applies through the C interface's
`LPM_Auto`/`LPM_Recommend` and `Cbc_resolve`.

## Measurement campaign

The November 2026 campaign used **403 instances: all 240 MIPLIB 2017 benchmark
instances plus 163 additional instances** from the local `2017+spp` collection,
including covering/packing/partitioning and other structured problems.
This means the complete 240-instance **benchmark set**, not the larger MIPLIB
2017 collection. The checked-in `../ml/auto_policy_results.csv` lists the exact
403 names, the ten policy features, selected rule/configuration, and the
baseline/policy penalized times and validation flags for all three seeds.

The frozen campaign tested 19 configurations x 403 instances x three seeds
(1, 2, 3): **22,971 runs**. It included dual variants, primal/Idiot and sprint
variants, the previous recommender, and strict barrier with native, dense,
and University of Florida/AMD Cholesky plus crossover. No barrier exception
was retained. The final measurements include Clp's November dual-cycling fix.

- Host: Linux aarch64, 128 Neoverse-N1 cores, approximately 250 GiB RAM.
- Optimized static build, `-O3 -march=native -mtune=native -ffp-contract=off`.
- One sequential LP solve per process, BLAS/OpenMP each limited to one thread.
- 32 concurrent jobs, 12G cgroup memory cap per process, swap disabled.
- 10,800-second solver limit, 600-second overtime grace, then process timeout.
- `-initialSolve -rowReductions force -timeMode elapsed` so LP-only experiments
  run the same pre-root strengthening as the ordinary root solve: bound
  propagation, clique/coefficient strengthening, and row reductions.
- Root features captured immediately before method selection, after
  strengthening but before Clp LP presolve. Repeated features were required
  to agree across methods and seeds.
- Success requires OPTIMAL status, exit code zero, and independent primal
  feasibility and dual/optimality validation (`lp_feasible=yes`,
  `lp_optimal=yes`, validation tolerances 1e-6). Strict barrier runs cannot
  silently fall back to another method.
- Every unsuccessful seed costs 21,600 seconds (twice the solver limit).
  Per-instance cost is the arithmetic mean of three seed costs; the shifted
  geometric mean (SGM) uses a one-second shift.

`trdta5581` consistently terminates infeasible during propagation, before
selection, so policy metrics cover the other **402 instances / 1,206 seeds**.

Frozen source revisions:

| Component | Commit |
|---|---|
| CoinUtils | `0073dcdd3e81e055cbf6514fc6b90406a38bdad7` |
| Osi | `f9502c198d12313c9c9fe5019e65effe67919497` |
| Clp | `4d57ccdbba1f57197b1e83d85f17f4e569bf7d24` |
| Cgl | `2c48cb796730a3093da85044a4ba26680f3e777d` |
| Cbc | `7f97b27de6baf84793e6b96bdb00abe4b1f55f7e` |

The snapshot also records source diffs, instrumented feature-capture and
strict-barrier objects, input/binary/library hashes, and commands. The
authoritative raw-results SHA-256 is
`a1dd22611f4b5d631ed4ca2566da9dcf0b4320855fbdccbf4ee86c84b8081f19`.
Later offline analyses explicitly allowed and recorded host libc/libm/loader
updates; they did not rerun or mix solver measurements.

## Four exceptions versus the best single configuration

The comparator is `-psi 0 -pertValue 72 -idiot 0 -sprint 0 -lpMethod dual`,
not the retired RF recommender and not an instance-wise oracle.
It is the best single setting by SGM both on all 402 eligible instances and
on the 339-instance fallback subset. The next-best fallback is dual_pertv58
(residual SGM 2.735 s versus 2.575 s for dual_pertv72).

| Metric | Best single everywhere | Four exceptions + best single fallback |
|---|---:|---:|
| Eligible instances | 402 | 402 |
| SGM, shift 1 s | 4.191659 s | 3.394919 s |
| Sum of per-instance mean penalized time | 101,471.010 s | 62,891.977 s |
| Big improvements | - | 44 |
| Significant regressions | - | 0 |
| Validated seeds | 1,200 / 1,206 | 1,202 / 1,206 |
| Instances validated in all three seeds | 400 | 400 |

This is **19.0% lower SGM** and **38.0% less total penalized time**.
The exceptions cover 63/402 instances (15.7%); 57 are faster and six are
slower. The largest slowdown is 0.463334 seconds. Rule 4 recovers two seeds
of `a2864-99blp`, but one still fails; `supportcase19` still fails all seeds.
There are no newly failed seeds.

| Exception | Incremental mean penalized seconds saved |
|---|---:|
| 1 | 5,789.108 |
| 2 | 8,754.195 |
| 3 | 11,661.418 |
| 4 | 12,374.312 |
| Total | 38,579.033 |

A big improvement requires `(baseline+1)/(alternative+1) >= 2` **and**
at least one second saved. A significant regression requires
`(alternative+1)/(baseline+1) > 1.10` **and** at least one second lost.
The zero significant-regression claim does not mean zero slowdowns.

### Interpretation and limitations

These are **fitted historical results**, assembled from measured configurations,
not fresh timings of the integrated selector and not independent test-set
estimates. The search froze one supported exception at a time and excluded its
covered instances. Candidates required at least 15 fitting instances, five big
wins, at most 25% coverage, no significant fitting regressions and no new failed
seeds. Depth-three conjunctions were searched with a bounded 64-mask beam.

Random tenfold rule-generation results were encouraging for the early rules,
but family-held-out evidence was weaker; the fourth search had four random-fold
wins and one regression. Fold-trained rules differ from these fixed rules.
Residual folds are conditional on prior full-data selections, not independent
validation of the entire ordered policy. Repeated exploration used the same
historical data. In particular, SCP and related-family concentration limits
generalization claims.

A fifth exception was investigated but **not included** after held-out rules
produced large proportional slowdowns. No claim is made here about full-MIP
speedup or uniform behavior on unseen instances.

## Replaying the historical comparison

From a configured build with installed Osi headers:

```sh
python3 test/lp-tuning/ml/replay_auto_policy.py
# Or supply the directory containing OsiFeatures.hpp:
python3 test/lp-tuning/ml/replay_auto_policy.py --include-dir=/prefix/include/coin-or
```

This compiles the actual production scorer in a temporary shared library,
checks its choices against every eligible CSV row, and recomputes the summary
with Python's standard library. It does not solve LPs or fabricate fresh timings.
Unused features are zero-filled; the CSV contains every feature read by this
policy. Boundary, precedence and NaN-imputation tests also run as part of
`lp-relaxation-test`, alongside C API LP solves, racing and independent
reoptimization controls.

To additionally exercise measured examples of all four exceptions through
the C API, comparing their LP objectives with explicit dual/pertValue72 and
checking both inherited and explicitly overridden reoptimization settings:

```sh
make -C test lp-relaxation-test
test/lp-relaxation-test --data-dir=/path/to/mip-sanity-data \
  --auto-instance-dir=/path/to/2017+spp
```

The optional directory must contain `scpk1_core`, `qap10`,
`neos-3226448-wkra`, and `neos-807456` as `.mps.gz` files. It is not required
for the ordinary self-contained LP test suite.

## Integration checks

The compiled policy replay reproduces all 402 eligible selections and the
summary above. Fresh command-line LP solves exercised the four examples above
and the fallback on `p0033`, each compared with explicit dual/pertValue72.
All ten runs passed independent primal-feasibility and LP-optimality checks;
auto and dual objectives agreed within 1e-6 absolute/relative tolerance.

The optional C API examples passed under both `LPM_Auto` (zero threads) and
`LPM_Recommend` (two threads, still sequential). They verify restoration of
the default perturbation 50 after each root solve, then check `Cbc_resolve`
with both inherited settings and explicit primal/perturbation61/presolve-off
overrides. The ordinary LP API/racing/reoptimization suite passed 219/219
reported checks; C API MIP smoke tests passed 9/9.

A separate full-MIP sanity sweep covered 508 instances: 504 validated,
four hard-timeout overtimes, zero wrong-answer failures and zero harness
errors. The pre-change workspace sweep had 505 validated and three overtimes.
The additional overtime was `supportcase12`: with the new fallback its root
LP reached optimality in 58.5 seconds, then the run timed out after entering
the later MIP pipeline. The comparison tool reported 41 performance
regressions and 39 improvements (including gap changes at limited solves).
This is not a clean isolated performance comparison: the PR worktree is based
on newer upstream master including independent reoptimization controls,
whereas the pre-change workspace used the older measurement revision.
Neither sweep establishes an end-to-end MIP speedup. The claimed 19.0%/38.0%
improvements refer only to the frozen root-LP comparison above.
