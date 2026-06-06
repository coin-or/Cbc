# MIPster Development Instructions

> ⚠️ **IMPORTANT — Git commits:** **NEVER commit all files blindly.** Always review `git status` and `git diff --staged` to select only the files relevant to the change. Committing unrelated files (build artifacts, logs, experiment outputs, scratch files) is a very common mistake — stage files explicitly with `git add <file>` and never use `git add .` or `git add -A`.

## Repository Overview

This is a development workspace for **MIPster** (`github.com/h-g-s/mipster`) — an open-source mixed-integer linear programming (MIP) solver written in C++, forked from COIN-OR Cbc. The fork diverges from upstream in several important ways (see *What was removed* below) and is rebranded as MIPster.

The **repository root is the monorepo** — it embeds CoinUtils, Clp, and Cgl sources directly under `src/` and builds everything in a single autotools project.

**Build target:** the repo root. `configster` operates directly on the repo root.

**Component layering (architecture reference only — the build system handles ordering automatically):** `CoinUtils → Clp (includes Osi) → Cgl → Cbc`. This means Cgl may reference Clp/CoinUtils but not Cbc, etc.

> **Note:** Osi sources (`OsiSolverInterface`, `OsiBranchingObject`, `OsiFeatures`, etc.) are compiled as part of the monorepo under `src/OsiCbc/`. They are not a separate project.

### What was removed from upstream Cbc

- **GLPK/GMPL support** — removed from CoinUtils and Cbc
- **Non-Clp LP solvers** — Cplex, DyLP, Mosek, SoPlex, HiGHS, Symphony support removed from Cbc and Cgl
- **GPL dependencies** — readline (interactive mode) and nauty removed
- **Interactive mode** — `mipster` binary no longer enters interactive mode when called without arguments; it prints help and exits
- **Osi project** — merged into the monorepo

### CI & Packaging

GitHub Actions workflows live in `.github/workflows/`. There is a single workflow file:

| Workflow | Platforms | Produces |
|---|---|---|
| `linux-ci.yml` | Linux x86_64, Linux aarch64, macOS x86_64, macOS arm64, Windows | Shared library `.tar.gz` for each platform; `.deb`, `.rpm`, `.pkg.tar.zst` packages on Linux |

On push to `main` and on tags, CI also uploads release assets. The library is named **`libmipster`** and headers install to `include/mipster/`.

> ⚠️ **IMPORTANT — Before tagging a release:** always bump the version in `configure.ac` (`AC_INIT` line) to match the new tag **before** pushing the tag. The PyPI publish workflow validates that the tag version and `configure.ac` version match and will fail otherwise.

## Environment Setup

Set these variables in your shell profile to match your local directory layout:

```sh
export MIPSTER_PREFIX="$HOME/prog/cbc"              # install prefix (opt build)
export MIPSTER_CBC_2X="$HOME/dev/cbc-stable"        # stable Cbc 2.x reference build
export MIPSTER_CBC_MASTER="$HOME/dev/cbc-master"    # upstream Cbc master reference build
export MIPSTER_INSTANCES="$HOME/inst"               # root of benchmark instance collections
export MIPSTER_EXPERIMENTS="$HOME/experiments/cbc"  # root of experiment results
```

Debug builds install to `$MIPSTER_PREFIX-dbg` (ASan) and `$MIPSTER_PREFIX-tsan` (TSan) by default when using `configster`.

## Build Commands

The workflow is always: **configure → build → install** in dependency order.

### Preferred workflow — `configster` (all-in-one)

`configster` is the **single entry point** for configuring, building, and installing
the full stack.
```sh
# Interactive TUI (guides you through all options):
./configster

# Optimised build — configure + build + install in one shot:
./configster --opt --install

# Debug build with AddressSanitizer:
./configster --debug --sanitizer=asan --install

# Debug build with ThreadSanitizer:
./configster --debug --sanitizer=tsan --install

# Valgrind-friendly debug build (no sanitizer, keeps frame pointers):
./configster --debug --sanitizer=valgrind --install

# Configure only (requires a previously installed prefix for .pc files):
./configster --opt

# Custom prefix, jobs, shared libs:
./configster --opt --install --prefix=/opt/cbc --jobs=8 --shared
```

Key options:

| Option | Default | Description |
|---|---|---|
| `--opt` / `--debug` | `--opt` | Build mode |
| `--sanitizer=none\|asan\|tsan\|valgrind` | `none` | Sanitizer (debug only) |
| `--march-native` / `--no-march-native` | on | `-march=native -mtune=native` |
| `--static` / `--shared` / `--both` | `--static` | Library type |
| `--prefix=PATH` | `~/prog/cbc` | Install prefix |
| `--jobs=N` | `nproc` | Parallel make jobs |
| `--install` | off | Build + install after configure |
| `--dry-run` | off | Show commands without running |

> **Important:** When `--install` (or `--build`) is given, `configster` interleaves
> configure → build → install **per project**, so each dependency's `.pc` files are in
> `$PREFIX` before the next project is configured. This is required on a fresh prefix.
> Configure-only mode assumes a prior install already populated `$PREFIX/lib/pkgconfig`.

### Incremental rebuild after source changes

The build system automatically descends into subdirectory libraries (CoinUtils, Clp, Cgl) and rebuilds them when their source files change.

```sh
# From repo root or src/ directory - automatically rebuilds all changed sources:
cd src && make -j$(nproc) && make install

# Or with verbose output:
make V=1 -j$(nproc) && make install
```

The Makefile uses `SUBDIRS = CoinUtils Clp Cgl` so changes to any file in:
- `src/CoinUtils/` → rebuilds CoinUtils library
- `src/Clp/` → rebuilds Clp library  
- `src/Cgl/` → rebuilds Cgl library
- `src/*.cpp` → rebuilds main Cbc sources

All automatically, in dependency order.

**Verify changes took effect:**
```sh
# Check if debug string appears in binary:
strings src/Cgl/.libs/libCgl.a | grep "your debug string"

# Or run the binary:
./src/mipster test/fixtures/gesa3.mps.gz -solve 2>&1 | grep "your debug output"
```

### Adding new source files

When adding new `.cpp`/`.hpp` files to the monorepo, edit `src/Makefile.am`:

1. **Edit `src/Makefile.am`** — add the `.cpp` to the appropriate `lib*_la_SOURCES` and the `.hpp` to `includecoin_HEADERS`.
2. **Regenerate `Makefile.in`** using system autotools:
   ```sh
   autoreconf --install -I BuildTools
   ```
3. **Reconfigure + build + install** with `./configster --opt --install`.

> **System autotools are sufficient.** The repo root has a `BuildTools` symlink (pointing to `../BuildTools`) and `AC_CONFIG_MACRO_DIRS([BuildTools])` so `aclocal` finds all COIN-OR m4 macros automatically. If the symlink is broken (e.g. after a fresh clone), verify that `../BuildTools` exists relative to the repo root.

### Clean
```sh
make clean        # remove build objects (keeps Makefiles)
make distclean    # full clean (removes Makefiles too; requires reconfigure)
```

### Run tests

The monorepo test suite lives in `test/`:

```sh
# Run the monorepo test suite:
cd test && make -j$(nproc) test
```

Individual test binaries: `test/CInterfaceTest`, `test/CbcSolverLpTest`, etc. They can be run directly.

**Writing integration tests?** See `doc/testing-guide.md` for comprehensive guidance on:
- Choosing diverse test instances
- Validating objectives and solution feasibility (C API and external tools)
- Solution pool validation (`Cbc_checkFeasibility`, `Cbc_savedSolution`)
- Cross-solver validation
- Automated debugging with `mip_diag`
- CI integration

### Test Validation Priorities

**Feasibility is paramount** — an infeasible solution is always wrong, regardless of objective or optimality status.

When writing MIP tests, validate in this order:

1. **ALWAYS check feasibility** (highest priority):
   ```c
   const double *bestSol = Cbc_getColSolution(m);
   double maxViolRow = 0.0; int rowIdx = -1;
   double maxViolCol = 0.0; int colIdx = -1;
   if (!Cbc_checkFeasibility(m, bestSol, &maxViolRow, &rowIdx, &maxViolCol, &colIdx)) {
     printf("FAIL: infeasible solution\n");
     return 0;
   }
   // Also check ALL solutions in the pool
   for (int s = 0; s < Cbc_numberSavedSolutions(m); s++) {
     const double *sol = Cbc_savedSolution(m, s);
     if (!Cbc_checkFeasibility(m, sol, ...)) {
       printf("FAIL: pool solution %d infeasible\n", s);
       return 0;
     }
   }
   ```

2. **IF solver claims optimality**, verify objective:
   ```c
   int is_optimal = Cbc_isProvenOptimal(m);
   if (is_optimal) {
     double obj = Cbc_getObjValue(m);
     if (fabs(obj - expected_obj) > tolerance) {
       printf("FAIL: claims optimal but obj=%.6f != expected=%.6f\n", obj, expected_obj);
       return 0;
     }
   }
   ```

3. **IF solver doesn't claim optimality**, accept feasible solution:
   - Don't fail if time/node limit hit
   - Don't check objective (solver may have found suboptimal solution)
   - Feasibility alone is sufficient

**Rationale:**
- An infeasible solution indicates a serious bug (wrong cuts, wrong postsolve, etc.)
- A feasible solution without optimality claim is acceptable (time limit, incomplete search)
- Only when solver claims "proven optimal" must we verify the objective matches expected

**Example test logic:**
```c
int pass = 1;
int is_optimal = Cbc_isProvenOptimal(m);
double obj = Cbc_getObjValue(m);

// 1. Check feasibility (always)
if (!Cbc_checkFeasibility(m, bestSol, ...)) {
  printf("FAIL: infeasible\n");
  pass = 0;
}

// 2. Check objective only if claims optimal
if (is_optimal && fabs(obj - expected) > tol) {
  printf("FAIL: claims optimal but wrong obj\n");
  pass = 0;
}

// 3. Report result
if (pass) {
  if (is_optimal)
    printf("PASS (proven optimal)\n");
  else
    printf("PASS (feasible, not proven optimal)\n");
}
```

## Code Formatting

All C++ source files use **clang-format** with the project's `.clang-format` config (WebKit-based style, 2-space indent, no column limit). Apply formatting with:

```sh
# Format all sources in a project:
./format-all-sources.sh

# Format a single file:
clang-format -i path/to/file.cpp
```

Key style rules: 2-space indent, no tabs, `PointerAlignment: Right`, `BreakBeforeBraces: WebKit` (function braces on new line, control braces same line), no column limit.

## Architecture

### Project Layers

- **CoinUtils** — Utility types and data structures (`CoinPackedMatrix`, `CoinMessage`, `CoinBronKerbosch`, adjacency vectors, clique utilities)
- **Clp** — LP solver (simplex). Also contains the merged Osi layer:
  - `src/OsiCbc/` — Abstract LP solver interface implementation (`OsiSolverInterface`, `OsiBranchingObject`, `OsiFeatures`, etc.)
- **Cgl** — Cut generators. Each cut type is in its own subdirectory under `src/Cgl/` (e.g. `CglBKClique`, `CglGomory`, `CglOddWheel`, `CglCliqueStrengthening`)
- **Cbc** — MIP solver. Orchestrates B&B using Clp as LP relaxation solver and Cgl cut generators

### Key Cbc Source Files

| File | Purpose |
|---|---|
| `CbcModel.cpp/.hpp` | Core B&B model; central class for the solver |
| `CbcSolver.cpp/.hpp` | Command-line driver and parameter handling |
| `CbcParameters.cpp/.hpp` | Parameter definitions |
| `Cbc_C_Interface.cpp/.h` | Public C API |
| `CbcHeuristic*.cpp` | Primal heuristics (FPump, RINS, RENS, Dive variants, etc.) |
| `CbcBranch*.cpp` | Branching objects and strategies |
| `CbcTree.cpp/.hpp` | B&B search tree |
| `CbcInstanceFeatures.cpp/.hpp` | Extracts MIP instance features to CSV (for ML/analysis) |

### Conflict Graph Infrastructure

Cbc has a conflict graph subsystem (recent addition) enabling:
- `CglBKClique` — Bron-Kerbosch clique separator (default, replaces `CglClique`)
- `CglOddWheel` — Odd-wheel cut separator
- `CglCliqueStrengthening` — Preprocessing for set-packing constraints

Controlled via `-cgraph`, `-clqstr`, `-clique`, `-oddwheel` parameters.

## Running the Solver

```sh
# Basic solve:
$MIPSTER_PREFIX/bin/mipster problem.mps -solve

# Debug launch example (see .vscode/launch.json):
./src/mipster inf1.lp -preprocess off -clqstr before -solve -solu a.sol

# Performance profiling:
./scripts/run_perf.sh   # wraps perf record on an MIPLIB instance
```

> ⚠️ **IMPORTANT — Parameters order:** only parameters added **before** `-solve` influence the solver behavior. With a few exceptions (e.g. `-solu` action) all parameter settings should be added _before_ `-solve`.

## Testing & Debugging

- **AddressSanitizer build:** `./configster --debug --sanitizer=asan --install` (installs to `$MIPSTER_PREFIX-dbg`)
- **ThreadSanitizer build:** `./configster --debug --sanitizer=tsan --install` (installs to `$MIPSTER_PREFIX-tsan`)
- `LD_LIBRARY_PATH` must include `.libs/` subdirs of each project when running binaries from the build tree (see `.vscode/launch.json` for the full path list)

### Diagnosing Wrong Results with `mipster_diag`

`mipster_diag` is the **first tool to reach for** when MIPster claims wrong-optimal or
wrong-infeasible. It automates the two most useful debugging passes and produces clean,
actionable output.

**Default flow (fast, ~1–3 min):**
1. Baseline run — confirm the bug reproduces within the time limit.
2. Auto-obtain reference solution (if `-sol` not provided).
3. `debugCuts` pass — activates `OsiRowCutDebugger` and prints:
   - `bad row N lb <= act <= ub` — an invalid row cut excludes the reference solution
   - `nodeBoundProp BAD FIXING (phase,rnd): col N ...` — a bound propagation step wrongly
     fixed a variable outside the reference solution's value

```sh
# Minimal usage: supply the certified optimal value
mipster_diag problem.mps -opt <certified_value>

# Provide a known-good solution to skip auto-obtain
mipster_diag problem.mps -opt <N> -sol /tmp/ref.sol

# Reproduce bugs triggered by specific parameter settings
mipster_diag problem.mps -opt <N> -p nodeBoundProp on -p nodeBoundPropMinDepth 2

# 60 s per run is the default; increase for hard instances
mipster_diag problem.mps -opt <N> -time 300

# Add -scan for the full 13-config feature-disabling sweep (slow)
mipster_diag problem.mps -opt <N> -sol /tmp/ref.sol -scan
```

**Actionable output example** (enlight_hard, wrong-infeasible bug):
```
nodeBoundProp BAD FIXING (propagation, round 0): col 156 (y#9#2) type=binary old=[1,3] new=[1,1] but optimal has 2
* 156 2
BAD 156 1 <= 2 <= 1
```
→ `nodeBoundProp` wrongly fixed column 156 to `[1,1]` even though the optimal value is 2.
The bug is in the bound propagation code, not in any row cut.

If the debugCuts pass shows no `bad row` or `BAD FIXING` lines, add `-scan` to run
the 13-config feature sweep and identify which cut family or feature is responsible.

### Debugging Invalid Cuts (manual)

For finer control beyond `mipster_diag`, use `-debugCuts` directly:

```sh
# Get a known good solution (preprocessing off)
$MIPSTER_PREFIX/bin/mipster problem.mps -preprocess off -solve -solu /tmp/known.sol

# Reproduce the wrong result with cut debugging active
$MIPSTER_PREFIX/bin/mipster problem.mps -debugCuts /tmp/known.sol -solve 2>&1 | grep "bad row"
# Output: "bad row <idx> <lb> <= <activity> <= <ub>" for any invalid cut
```

**Full workflow, cut type codes, tolerance details, and the TwoMIR bug analysis are documented in:**
`doc/cut-debugging-howto.md`

### Investigating Wrong-Optimal Failures from CI

When a CI test claims optimal at a wrong objective value (e.g. `CInterfaceTest_miclsp` prints
`FAIL: claimed optimal obj=X != certified opt=Y`), the test automatically runs two diagnostic
passes defined in `test/mip_diag.h`:

1. **`mip_diag_wrong_optimal`** — re-solves with 13 configurations, each disabling one cut family
   or feature. Runs are time-limited (no node limit), so even bugs that only manifest deep in
   the tree can reproduce.
2. **`mip_diag_debug_cuts`** — re-solves with `OsiRowCutDebugger` active against the reference
   solution file, printing any cut that would exclude it.

#### Reading the diagnostic output from CI logs

Fetch the log for the failing job and grep for `DIAG`:

```sh
gh run view <RUN_ID> --log | grep -E "DIAG|LEAD|bad row"
```

Each diagnostic result is printed on its own labeled line:

```
[DIAG  1/13] no-cgraph    → not proven  obj=8235  bound=7620  gap=7.5%  (time limit)
[DIAG  8/13] no-flow      → not proven  obj=8181  bound=8050  gap=1.6%  (time limit)  *** LEAD ***
[DIAG 12/13] no-twomir    → WRONG       obj=8204  (same wrong result)
```

**Interpretation:**
- `OK` — disabling this feature fixed the bug (it was the cause).
- `WRONG` — same wrong result persists; this feature is not the cause.
- `not proven` — ran out of time without proving optimality; obj value still informative.
- `*** LEAD ***` — this config found the **certified optimal obj** even if not proven; strong
  signal that the disabled feature contains the faulty cut. Prioritise these for follow-up.

#### Follow-up steps

1. **Identify LEADs** — the feature(s) marked `LEAD` are the prime suspects.
2. **Reproduce locally** with `mipster_diag` (pass LEAD's parameter via `-p`):
   ```sh
   mipster_diag problem.mps -opt <certified_value> -p <param> off
   ```
3. **Run debugCuts** to find the exact invalid cut (requires the reference `.sol` file):
   ```sh
   $MIPSTER_PREFIX/bin/mipster problem.mps -debugCuts /path/to/known.sol -solve 2>&1 | grep "bad row"
   ```
4. **Cross-reference** the cut type code from `bad row` output with `doc/cut-debugging-howto.md`
   to identify the generator responsible.

#### Benign CI cancellations

When a newer push arrives while a CI run is in flight, GitHub Actions cancels it with
*"Canceling since a higher priority waiting request exists"* — this shows as a failure but
is **not a real problem**. Always check `gh run list` to confirm whether a more recent
successful run superseded the cancelled one.

### Validating Solution Feasibility with `mipster_validate_sol`

`mipster_validate_sol` checks whether a `.sol` file is truly feasible for a given problem —
useful for catching preprocessing bugs or other solver errors that produce infeasible solutions.

```sh
mipster_validate_sol problem.mps[.gz] solution.sol
```

It checks:
- Variable bounds (`lb ≤ x[i] ≤ ub`)
- Integrality (`|x[i] − round(x[i])| ≤ 1e-5` for integer variables)
- Row activities (`rowLB ≤ Ax ≤ rowUB` for every constraint)
- Objective value (computed `c·x` vs. claimed in `.sol` header)

Exit codes: `0` = feasible, `1` = violations found, `2` = file/usage error.

To validate all solutions from an experiment:
```sh
for sol in /path/to/experiment/*.sol; do
  inst=$(basename "$sol" .sol)
  mipster_validate_sol "$MIPSTER_INSTANCES/miplib/2017+spp/${inst}.mps.gz" "$sol" || \
    echo "INFEASIBLE: $inst"
done
```

Full documentation: `docs/utils.html` (mipster_validate_sol section).

## Hardware & Parallelism

Exploit all cores of the machine where the experiment is running.

- **`GNU parallel` is installed** — use it to dispatch independent solver runs in parallel.
- Each MIPster solve is single-threaded, so parallelism comes from running many instances concurrently.
- Tune `-j<N>` to the number of available cores (`$(nproc)` for local machines; higher values for many-core servers).

### Typical parallel experiment pattern

```sh
# Run instances in parallel (adjust -j to available cores):
ls $MIPSTER_INSTANCES/miplib/2017/*.mps.gz | \
  parallel -j$(nproc) '$MIPSTER_PREFIX/bin/mipster {} -sec 300 -solve > {/.}.log 2>&1'

# With a list file and custom command template:
cat instances.txt | \
  parallel -j$(nproc) '$MIPSTER_PREFIX/bin/mipster $MIPSTER_INSTANCES/super/{} -sec 300 -solve > results/{/.}.log 2>&1'

# Preserve output order, show progress:
parallel --eta -j$(nproc) ...
```

### Builds

`make` already uses all cores via `-j$(nproc)` in the build scripts. No change needed for builds.

## Instance Collections for Experiments

Benchmark collections available under `$MIPSTER_INSTANCES`:

| Path | Instances | Description |
|---|---|---|
| `$MIPSTER_INSTANCES/miplib/2017+spp` | 358 | *PREFERRED INSTANCE SET* — MIPLIB 2017 plus set packing/partitioning/covering instances |
| `$MIPSTER_INSTANCES/super/` | 1339 | Large curated superset (includes miplib); grouped by problem class |
| `$MIPSTER_INSTANCES/miplib/2017/` | 251 | MIPLIB 2017 benchmark set |
| `$MIPSTER_INSTANCES/cbc-test-set/` | ~985 | CBC regression/test set (includes `.bas` warm-start files) |

All instances are compressed `.mps.gz` files. MIPster reads them directly:
```sh
$MIPSTER_PREFIX/bin/mipster $MIPSTER_INSTANCES/miplib/2017+spp/air05.mps.gz -solve
```

Experiment runner scripts in the `scripts/` directory (e.g. `scripts/run_experiments.sh`, `scripts/summarize_results.py`) use these paths.

### Instances for quick tests

The file `solvable_instances.tsv` lists instances solvable in under 5 minutes. Use it to pick quick smoke-test candidates and plan runtimes accordingly.

### Experiment Results Directory

Past experiment results are stored in `$MIPSTER_EXPERIMENTS/`, with one subdirectory per experiment (named by type and date, e.g. `opt_YYYY_MM_DD`, `some_feature_YYYY_MM_DD`). Each experiment subdirectory contains one subdirectory per condition, and each condition contains per-instance `.result` files with the format:

```
STATUS  obj  bound  gap%  time_seconds  nodes  condition_name
```

Example: `SOLVED	26374	26374	0.00%	49.5	326	C0_baseline`

**Use this to pick instances for quick experiments:** before selecting an instance, check a recent experiment's `.result` files to find instances with acceptable solve times.

```sh
# Find instances solved in under 10 seconds in a recent experiment:
awk '$5 < 10' $MIPSTER_EXPERIMENTS/<experiment_name>/C0_baseline/*.result

# Show instance name and time, sorted by time:
for f in $MIPSTER_EXPERIMENTS/<experiment_name>/C0_baseline/*.result; do
  echo "$(basename $f .result) $(awk '{print $5}' $f)"
done | sort -k2 -n | head -20
```

### Capturing per-solve statistics to CSV

Two parameters control this — **both are needed**:

| Parameter | Role |
|---|---|
| `-csvStatistics <file.csv>` | Sets the output file (default: `stats.csv`). Appends one row per solve; writes a header if the file is new/empty. |
| `-writeStatistics [<file.csv>]` | **Action** — triggers the write. If no filename is given, uses the current `-csvStatistics` value. |

Each row captures time, gap, cut generation counts, node counts, and other per-solve metrics.
```sh
# Typical usage — set file then trigger write after solve:
$MIPSTER_PREFIX/bin/mipster problem.mps -csvStatistics results.csv -solve -writeStatistics

# Or pass the filename directly to -writeStatistics:
$MIPSTER_PREFIX/bin/mipster problem.mps -solve -writeStatistics results.csv
```

In parallel experiments, all jobs can safely append to the **same file** because each invocation appends a single row:
```sh
ls $MIPSTER_INSTANCES/miplib/2017/*.mps.gz | \
  parallel -j$(nproc) '$MIPSTER_PREFIX/bin/mipster {} -sec 300 -csvStatistics results.csv -solve -writeStatistics'
```

`CbcInstanceFeatures` extracts OsiFeature values and appends a row to a CSV (used for ML-based parameter tuning research); see `$MIPSTER_INSTANCES/miplib/2017+spp/features.csv` for an example.

## Experiment Scripts

All experiment and analysis scripts live in `scripts/`:

| Script | Purpose |
|---|---|
| `scripts/run_experiments.sh` | Run MIPster on a set of instances, verify objectives, collect logs |
| `scripts/summarize_results.py` | Rich terminal summary of a `run_experiments.sh` output directory |
| `scripts/compare_experiments.py` | Side-by-side PDF comparison of two experiment directories |
| `scripts/compare_multi_experiments.py` | Multi-way comparison (N≥2 experiments) to PDF |
| `scripts/make_report.py` | Multi-page PDF report from a single experiment directory |
| `scripts/generate_features.sh` | Extract `CbcInstanceFeatures` CSV for a folder of instances |
| `scripts/run_perf.sh` | Profile MIPster with `perf record` |

## LP Relaxation & Racing Portfolio Analysis

Scripts for tuning LP relaxation parameters and selecting optimal racing portfolios
for the **racingLP** feature (multiple CLP runs in parallel, first to finish wins):

| Script | Purpose |
|---|---|
| `scripts/run_lp_experiments.sh` | Run LP solves across params × instances × seeds |
| `scripts/analyze_lp_params.py` | Individual param rankings (SGM, head-to-head, errors) — produces `lp_avg_times.csv` |
| `scripts/analyze_racing_portfolios.py` | Portfolio analysis: greedy + exhaustive optimal K-param racing sets |
| `scripts/make_lp_report.py` | **PDF report**: 5-page visual summary (rankings, performance profiles, racing curve, key findings) |
| `scripts/summarize_lp_results.py` | Fake-optimal / wrong-result validation per instance |

**Full workflow documented in:** `doc/lp_racing_analysis_howto.md`

Quick re-run on an existing experiment:
```sh
EXP=$MIPSTER_EXPERIMENTS/lp_relax_YYYY_MM_DD
python3 scripts/analyze_lp_params.py --dir $EXP --timelimit 14400
python3 scripts/analyze_racing_portfolios.py --dir $EXP --timelimit 14400
python3 scripts/make_lp_report.py $EXP --timelimit 14400  # → lp_report.pdf
```

Key correctness note: **wrong results are penalised** — if a param reports OPTIMAL but
its objective exceeds the best-known LP value by `> 1e-4` (relative), it is treated as
a failure (penalty time) in both individual and racing analyses.
