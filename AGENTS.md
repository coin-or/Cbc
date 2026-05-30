# MIPster Development Instructions

> ⚠️ **IMPORTANT — Git commits:** **NEVER commit all files blindly.** Always review `git status` and `git diff --staged` to select only the files relevant to the change. Committing unrelated files (build artifacts, logs, experiment outputs, scratch files) is a very common mistake — stage files explicitly with `git add <file>` and never use `git add .` or `git add -A`.

## Repository Overview

This is a development workspace for **MIPster** (`github.com/h-g-s/mipster`) — an open-source mixed-integer linear programming (MIP) solver written in C++, forked from COIN-OR Cbc. The fork diverges from upstream in several important ways (see *What was removed* below) and is rebranded as MIPster.

The `Cbc/` directory is the main **monorepo** — it embeds CoinUtils, Clp, and Cgl sources directly under `src/` and builds everything in a single autotools project. The sibling directories (`CoinUtils/`, `Clp/`, `Cgl/`) are the original separate subprojects kept for reference and for the stable cbc-2.x build at `$MIPSTER_CBC_2X`.

**Build target:** just `Cbc` (the monorepo). `configster` and `cbc_build.sh` only operate on `Cbc/`.

**Dependency order within the monorepo:** `CoinUtils → Clp (includes Osi) → Cgl → Cbc`

> **Note:** Osi sources (`OsiSolverInterface`, `OsiBranchingObject`, `OsiFeatures`, etc.) are compiled as part of the monorepo under `Cbc/src/Osi/`. They are not a separate project.

### What was removed from upstream Cbc

- **GLPK/GMPL support** — removed from CoinUtils and Cbc
- **Non-Clp LP solvers** — Cplex, DyLP, Mosek, SoPlex, HiGHS, Symphony support removed from Cbc and Cgl
- **GPL dependencies** — readline (interactive mode) and nauty removed
- **Interactive mode** — `mipster` binary no longer enters interactive mode when called without arguments; it prints help and exits
- **Osi project** — merged into the monorepo

### CI & Packaging

GitHub Actions workflows live in `.github/workflows/` and build/package for three platforms:

| Workflow | Produces |
|---|---|
| `linux-ci.yml` | Manylinux (x86-64 + aarch64) shared library `.tar.gz`, then `.deb`, `.rpm`, `.pkg.tar.zst` packages |
| `macos-ci.yml` | macOS `.pkg` installer with headers, library, binary, and man page |
| `windows-ci.yml` | Windows NSIS installer with DLL, headers, and binary |

On push to `main` and on tags, CI also uploads release assets. The library is named **`libmipster`** and headers install to `include/mipster/`.

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

```sh
# Rebuild and reinstall (no reconfigure needed):
./cbc_build.sh --install

# Rebuild only (no install):
./cbc_build.sh

# Rebuild the monorepo directly:
cd Cbc && make V=1 -j$(nproc) && make install
```

> **Note:** Since everything is a single monorepo (`Cbc/`), any source change is picked up by a plain `make` in `Cbc/` — there is no inter-project cascade to worry about. `./cbc_build.sh --install` is always the safest option.

### Adding new source files

When adding new `.cpp`/`.hpp` files to the monorepo, edit `Cbc/src/Makefile.am`:

1. **Edit `Cbc/src/Makefile.am`** — add the `.cpp` to the appropriate `lib*_la_SOURCES` and the `.hpp` to `includecoin_HEADERS`.
2. **Regenerate `Makefile.in`** using system autotools:
   ```sh
   cd Cbc
   autoreconf --install -I BuildTools
   ```
3. **Reconfigure + build + install** with `./configster --opt --install`.

> **System autotools are sufficient.** `Cbc/` has a `BuildTools` symlink and `AC_CONFIG_MACRO_DIRS([BuildTools])` so `aclocal` finds all COIN-OR m4 macros automatically.

### Clean
```sh
./cbc_clean.sh                  # make clean + distclean in each project
./cbc_clean.sh --distclean-only # distclean only (removes Makefiles too)
```

### Run tests

The monorepo test suite lives in `Cbc/test/`:

```sh
# Run the monorepo test suite:
cd Cbc/test && make -j$(nproc) test
```

Individual test binaries: `Cbc/test/CInterfaceTest`, `Cbc/test/cbc_unittest`, etc. They can be run directly.

## Code Formatting

All C++ source files use **clang-format** with the project's `.clang-format` config (WebKit-based style, 2-space indent, no column limit). Apply formatting with:

```sh
# Format all sources in a project:
cd Cbc && ./format-all-sources.sh

# Format a single file:
clang-format -i path/to/file.cpp
```

Key style rules: 2-space indent, no tabs, `PointerAlignment: Right`, `BreakBeforeBraces: WebKit` (function braces on new line, control braces same line), no column limit.

## Architecture

### Project Layers

- **CoinUtils** — Utility types and data structures (`CoinPackedMatrix`, `CoinMessage`, `CoinBronKerbosch`, adjacency vectors, clique utilities)
- **Clp** — LP solver (simplex). Also contains the merged Osi layer:
  - `src/Osi/` — Abstract LP solver interface (`OsiSolverInterface`, `OsiBranchingObject`, `OsiFeatures`, etc.) compiled into `libOsi.la`
  - `src/OsiClp/` — Clp implementation of the Osi interface (`OsiClpSolverInterface`), compiled into `libOsiClp.la`
  - `src/OsiCommonTest/` — Shared unit test helpers (`libOsiCommonTest.la`)
  - `ClpSimplex` (stable) and `AbcSimplex` (experimental parallel version)
- **Cgl** — Cut generators. Each cut type is in its own subdirectory under `Cgl/src/` (e.g. `CglBKClique`, `CglGomory`, `CglOddWheel`, `CglCliqueStrengthening`)
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
./Cbc/src/mipster inf1.lp -preprocess off -clqstr before -solve -solu a.sol

# Performance profiling:
./run_perf.sh   # wraps perf record on an MIPLIB instance
```

> ⚠️ **IMPORTANT — Parameters order:** only parameters added **before** `-solve` influence the solver behavior. With a few exceptions (e.g. `-solu` action) all parameter settings should be added _before_ `-solve`.

## Testing & Debugging

- **AddressSanitizer build:** `./configster --debug --sanitizer=asan --install` (installs to `$MIPSTER_PREFIX-dbg`)
- **ThreadSanitizer build:** `./configster --debug --sanitizer=tsan --install` (installs to `$MIPSTER_PREFIX-tsan`)
- `LD_LIBRARY_PATH` must include `.libs/` subdirs of each project when running binaries from the build tree (see `.vscode/launch.json` for the full path list)

### Debugging Invalid Cuts

If MIPster produces a wrong (suboptimal) proven-optimal, the likely cause is a cut that incorrectly excludes the optimal solution. Use the `-debugCuts <solution.sol>` parameter to detect this automatically via `OsiRowCutDebugger`:

```sh
# Get a known good solution (preprocessing off)
$MIPSTER_PREFIX/bin/mipster problem.mps -preprocess off -solve -solu /tmp/known.sol

# Reproduce the wrong result with cut debugging active
$MIPSTER_PREFIX/bin/mipster problem.mps -debugCuts /tmp/known.sol -solve 2>&1 | grep "bad row"
# Output: "bad row <idx> <lb> <= <activity> <= <ub>" for any invalid cut
```

**Full workflow, cut type codes, tolerance details, and the TwoMIR bug analysis are documented in:**
`doc/cut-debugging-howto.md`

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
