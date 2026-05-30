# Contributing to MIPster

MIPster is an open-source mixed-integer programming (MIP) solver and we welcome contributions of all kinds — bug reports, performance improvements, new cut generators, heuristics, documentation, and experiments.

This guide walks you through setting up a complete development and benchmarking environment from scratch.

---

## Prerequisites

| Tool | Purpose |
|---|---|
| `git` + `gh` (GitHub CLI) | Repository management |
| GCC 7+ or Clang 6+ | C++14 compiler |
| `automake`, `autoconf`, `libtool`, `m4` | Build system |
| `clang-format` | Code formatting |
| `pkg-config` | Dependency discovery |
| `GNU parallel` | Parallel benchmark runs |
| `python3` | Experiment analysis scripts |
| `wget` or `curl` | Instance downloads |

On Ubuntu/Debian:
```sh
sudo apt install git gh build-essential automake autoconf libtool \
     pkg-config clang-format parallel python3 wget
```

On macOS (Homebrew):
```sh
brew install gh automake autoconf libtool pkg-config clang-format \
     parallel python3 wget
```

---

## 1. Clone the Repository

```sh
gh repo clone h-g-s/mipster
cd mipster
```

Or with plain git:
```sh
git clone https://github.com/h-g-s/mipster.git
cd mipster
```

---

## 2. Set Up Environment Variables

MIPster uses environment variables to locate builds, instances, and experiments. Add these to your shell profile.

<details>
<summary><strong>bash / zsh</strong> — add to <code>~/.bashrc</code> or <code>~/.zshrc</code></summary>

```sh
export MIPSTER_PREFIX="$HOME/prog/cbc"              # install prefix (dev build)
export MIPSTER_CBC_2X="$HOME/dev/cbc-stable"        # stable Cbc 2.x reference (optional)
export MIPSTER_CBC_MASTER="$HOME/dev/cbc-master"    # upstream Cbc master reference (optional)
export MIPSTER_INSTANCES="$HOME/inst"               # benchmark instance collections
export MIPSTER_EXPERIMENTS="$HOME/experiments/cbc"  # experiment results

# Add the installed binary to PATH:
export PATH="$MIPSTER_PREFIX/bin:$PATH"
```

Reload: `source ~/.bashrc`
</details>

<details>
<summary><strong>fish</strong> — add to <code>~/.config/fish/config.fish</code></summary>

```fish
set -gx MIPSTER_PREFIX "$HOME/prog/cbc"
set -gx MIPSTER_CBC_2X "$HOME/dev/cbc-stable"
set -gx MIPSTER_CBC_MASTER "$HOME/dev/cbc-master"
set -gx MIPSTER_INSTANCES "$HOME/inst"
set -gx MIPSTER_EXPERIMENTS "$HOME/experiments/cbc"

fish_add_path "$MIPSTER_PREFIX/bin"
```

Reload: `source ~/.config/fish/config.fish`
</details>

<details>
<summary><strong>PowerShell</strong> — add to your <code>$PROFILE</code></summary>

```powershell
$env:MIPSTER_PREFIX    = "$HOME\prog\cbc"
$env:MIPSTER_CBC_2X    = "$HOME\dev\cbc-stable"
$env:MIPSTER_CBC_MASTER = "$HOME\dev\cbc-master"
$env:MIPSTER_INSTANCES = "$HOME\inst"
$env:MIPSTER_EXPERIMENTS = "$HOME\experiments\cbc"

$env:PATH = "$env:MIPSTER_PREFIX\bin;$env:PATH"
```

Reload: `. $PROFILE`
</details>

### Variable reference

| Variable | Default | Description |
|---|---|---|
| `MIPSTER_PREFIX` | `~/prog/cbc` | Install prefix for the dev (optimised) build |
| `MIPSTER_CBC_2X` | `~/dev/cbc-stable` | Stable Cbc 2.x reference build — kept for comparison |
| `MIPSTER_CBC_MASTER` | `~/dev/cbc-master` | Upstream Cbc master reference build — kept for comparison |
| `MIPSTER_INSTANCES` | `~/inst` | Root directory for benchmark instance collections |
| `MIPSTER_EXPERIMENTS` | `~/experiments/cbc` | Root directory for experiment results |

Debug builds install automatically to `$MIPSTER_PREFIX-dbg` (AddressSanitizer) and `$MIPSTER_PREFIX-tsan` (ThreadSanitizer).

---

## 3. Build and Install

```sh
# Configure, build, and install the full stack in one command:
./Cbc/scripts/configster --opt --install

# Verify:
$MIPSTER_PREFIX/bin/mipster --help
```

For a debug build with AddressSanitizer:
```sh
./Cbc/scripts/configster --debug --sanitizer=asan --install
```

See [`AGENTS.md`](AGENTS.md) for the full build reference including incremental rebuilds, adding new source files, and sanitizer options.

---

## 4. Download MIPLIB 2017 Benchmark Instances

[MIPLIB 2017](https://miplib.zib.de/) is the standard benchmark suite for MIP solvers.
MIPster uses it extensively for performance testing and experiments.

```sh
# Create the instances directory:
mkdir -p $MIPSTER_INSTANCES/miplib/2017

# Download the benchmark set (240 instances, ~317 MB):
wget -O /tmp/benchmark.zip https://miplib.zib.de/downloads/benchmark.zip

# Extract:
unzip /tmp/benchmark.zip -d $MIPSTER_INSTANCES/miplib/2017
rm /tmp/benchmark.zip
```

> The full collection (1,000+ instances, ~3.5 GB) is available at  
> `https://miplib.zib.de/downloads/collection.zip` — only needed for large-scale experiments.

### Verify the download

```sh
ls $MIPSTER_INSTANCES/miplib/2017/*.mps.gz | wc -l
# Expected: 240

# Run a quick solve to confirm everything works:
$MIPSTER_PREFIX/bin/mipster \
  $MIPSTER_INSTANCES/miplib/2017/air05.mps.gz -sec 60 -solve
```

### Create the experiments directory

```sh
mkdir -p $MIPSTER_EXPERIMENTS
```

---

## 5. Running Experiments

All experiment scripts live in `scripts/`. They respect the environment variables set above and accept command-line flags to override defaults.

### 5a. Running a benchmark with `run_experiments.sh`

`run_experiments.sh` is the main experiment runner. It dispatches instances in parallel via GNU parallel, collects per-instance `.log` and `.result` files, verifies objectives, and writes a `summary.tsv`.

```sh
# Minimal: run on all MIPLIB 2017 instances, 5-minute limit
scripts/run_experiments.sh \
  --bin $MIPSTER_PREFIX/bin/mipster \
  --instances $MIPSTER_INSTANCES/miplib/2017 \
  --timelimit 300 \
  --parallel $(nproc) \
  --outdir $MIPSTER_EXPERIMENTS/my_first_run

# Quick smoke test on easy instances only:
scripts/run_experiments.sh \
  --bin $MIPSTER_PREFIX/bin/mipster \
  --instances $MIPSTER_INSTANCES/miplib/2017 \
  --timelimit 60 \
  --parallel $(nproc) \
  --outdir $MIPSTER_EXPERIMENTS/smoke_test
```

`solvable_instances.tsv` in the repository root lists instances known to solve in under 5 minutes — use it to filter a fast smoke-test set.

Key options:

| Flag | Description |
|---|---|
| `--bin PATH` | Path to the `mipster` binary |
| `--instances DIR` | Directory of `.mps.gz` instances |
| `--timelimit N` | Per-instance time limit in seconds |
| `--parallel N` | Number of concurrent jobs |
| `--outdir DIR` | Output directory (auto-named if omitted) |
| `--dry-run` | Print jobs without running |

Run `scripts/run_experiments.sh --help` for the full option list.

### 5b. Summarising results

After a run completes, inspect the results:

```sh
EXP=$MIPSTER_EXPERIMENTS/my_first_run

# Rich terminal table:
python3 scripts/summarize_results.py --outdir $EXP

# Plain text (CI-friendly):
python3 scripts/summarize_results.py --outdir $EXP --plain
```

### 5c. Comparing two experiments

After making a code change and running a new experiment, compare against a baseline:

```sh
BASELINE=$MIPSTER_EXPERIMENTS/baseline_run
NEW=$MIPSTER_EXPERIMENTS/my_change_run

# Side-by-side PDF comparison:
python3 scripts/compare_experiments.py $BASELINE $NEW \
  --label-a baseline --label-b my_change \
  comparison.pdf

# Compare more than two experiments:
python3 scripts/compare_multi_experiments.py \
  $BASELINE $NEW $MIPSTER_EXPERIMENTS/another_run \
  -o multi_comparison.pdf
```

### 5d. Generating a full PDF report

```sh
python3 scripts/make_report.py $MIPSTER_EXPERIMENTS/my_first_run report.pdf
```

The report includes performance profiles, per-instance breakdowns, solve rate curves, and a comparison table.

### 5e. Extracting instance features

```sh
# Extract features for all instances in a folder (writes features.csv):
scripts/generate_features.sh $MIPSTER_INSTANCES/miplib/2017

# Control parallelism:
scripts/generate_features.sh $MIPSTER_INSTANCES/miplib/2017 --jobs $(nproc)
```

### 5f. LP relaxation parameter tuning

These scripts are used to tune the `racingLP` feature — running multiple LP parameter configurations in parallel and keeping the first to finish.

```sh
# 1. Run LP experiments across parameter settings × instances × seeds:
scripts/run_lp_experiments.sh \
  --bin $MIPSTER_PREFIX/bin/mipster \
  --instances $MIPSTER_INSTANCES/miplib/2017+spp \
  --timelimit 3600 \
  --outdir $MIPSTER_EXPERIMENTS/lp_run_$(date +%Y_%m_%d)

EXP=$MIPSTER_EXPERIMENTS/lp_run_YYYY_MM_DD

# 2. Analyse individual parameter rankings:
python3 scripts/analyze_lp_params.py --dir $EXP --timelimit 3600

# 3. Find the best racing portfolio of K configurations:
python3 scripts/analyze_racing_portfolios.py --dir $EXP --timelimit 3600

# 4. Generate the full PDF report:
python3 scripts/make_lp_report.py $EXP --timelimit 3600  # → lp_report.pdf
```

Full workflow: `doc/lp_racing_analysis_howto.md`

### 5g. Performance profiling

```sh
# Profile a single instance (output: perf.data):
scripts/run_perf.sh $MIPSTER_INSTANCES/miplib/2017/air05.mps.gz
perf report
```

---

## 6. Code Style

All C++ source files must be formatted with `clang-format` using the project's config:

```sh
cd Cbc && ./format-all-sources.sh
```

Or for a single file:
```sh
clang-format -i Cbc/src/MyNewFile.cpp
```

Style: WebKit-based, 2-space indent, no tabs, no column limit.

---

## 7. Running Tests

```sh
cd Cbc/test && make -j$(nproc) test
```

Run a quick smoke test before submitting a PR:
```sh
./cbc_build.sh --install
cd Cbc/test && make -j$(nproc) test
```

---

## 8. Submitting Changes

1. **Fork** the repository and create a branch from `main`.
2. **Make your changes** — keep commits focused and atomic.
3. **Format code**: `cd Cbc && ./format-all-sources.sh`
4. **Run tests**: `cd Cbc/test && make -j$(nproc) test`
5. **Open a Pull Request** against `main` with a clear description.

### Commit messages

- Use the imperative mood: *"Add clique strengthening for set-covering rows"*
- Reference issues where relevant: *"Fix invalid cut in TwoMIR separator (#42)"*
- **Never commit build artifacts, logs, or experiment outputs** — only source files relevant to the change.

---

## Getting Help

- Open an issue on GitHub for bugs, questions, or feature requests.
- For development discussions, use GitHub Discussions.
