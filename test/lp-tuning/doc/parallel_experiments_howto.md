# Running Parallel Experiments (120 cores)

## TL;DR — The Pattern That Works

We use **GNU parallel** (`/usr/bin/parallel`, version 20231122) for all parallel
experiments. It replaced `xargs -P` and is strictly better for our use case.

> **Note:** before GNU parallel was installed, `/usr/bin/parallel` was the moreutils
> `parallel` — a completely different, far less capable tool. Always verify with
> `parallel --version` that you see `GNU parallel`.

---

## Canonical Script Pair

Every experiment consists of two scripts:

### 1. Worker script  (`run_one_<exp>.sh`)

Called once per job. Arguments come from the jobs file via GNU parallel's
`{1}`, `{2}` etc. positional replacements (or `{}` for a single argument).

**Critical rules (learned from failures):**
- Use **`realpath "$0"`** to find the script's own directory (not `dirname "$0"`, which
  breaks when called via `bash /path/to/script.sh`)
- Use **bash arrays** for option lists — never a plain string (word-splitting bugs):
  ```bash
  OPTS=("-dualpivot" "pesteep" "-psi" "1.0")   # correct
  OPTS="-dualpivot pesteep -psi 1.0"            # WRONG — splits badly on spaces
  ```
- Pass shell variables into `awk` with **`-v`**, never via string substitution:
  ```bash
  awk -v s="$START" -v e="$END" 'BEGIN {printf "%.3f", e-s}'   # correct
  awk 'BEGIN {printf "%.3f", '"$END"'-'"$START"'}'              # WRONG
  ```
- For floating-point comparison use **python3** with `sys.argv` — never shell `bc` or
  awk for toleranced comparisons:
  ```bash
  python3 -c "
  import sys, math
  a, b = float(sys.argv[1]), float(sys.argv[2])
  print('true' if abs(a-b) <= 1e-6*(1+abs(b)) else 'false')
  " "$OBJ" "$REF_OBJ"
  ```
- Use **`printf`** not `echo -e` for TSV output (portable, no escape surprises)
- Workers **may exit non-zero** for genuine setup errors (binary missing, bad args);
  exit 0 means "job was processed" (even if CBC itself failed/was killed). GNU parallel
  records exit codes in the joblog — use `--resume-failed` to retry non-zero jobs.

### 2. Orchestrator script (`run_exp_<exp>.sh`)

**Critical rules:**
- Use **`realpath "$0"`** to resolve WORKER path
- Use **`$HOME`** not `~` in non-interactive scripts (tilde not expanded in all
  contexts)
- Add sanity checks before starting (CBC binary, worker, input files) — fail fast,
  not thousands of jobs later
- Launch with **`nohup … &`** so the experiment survives terminal close:
  ```bash
  nohup bash run_exp_myexp.sh > /tmp/myexp.log 2>&1 &
  echo "PID: $!"
  ```

**GNU parallel invocation:**
```bash
parallel -j 120 --joblog "$OUT_DIR/joblog.tsv" --resume --halt never \
  /bin/bash "$WORKER" {1} {2} "$OUT_DIR" "$CBC_BIN" "$REF_TSV" \
  :::: "$OUT_DIR/jobs.txt"
```

Key flags explained:

| Flag | Purpose |
|---|---|
| `-j 120` | 120 parallel workers |
| `--joblog <file>` | Records seq, timing, exit code for every job |
| `--resume` | Skips jobs already in joblog (by Seq) — built-in idempotency, **replaces manual result-file checks** |
| `--resume-failed` | Alternative: also re-runs jobs that previously exited non-zero |
| `--halt never` | Never abort the run even if workers exit non-zero |
| `--progress` | Live progress to stderr (jobs/sec, ETA) |
| `--eta` | ETA estimate |
| `--colsep '\t'` | Split each input line on tab for `{1}`, `{2}` positional args |
| `::::` | Read job arguments from a file (vs `:::` for inline values or stdin) |

**Multi-column jobs file** (instance + method on same line, tab-separated):
```bash
# Generate jobs.txt with two columns
for inst in "${INSTANCES[@]}"; do
  for method in idiot30 idiot60 dualpesteep idiot80; do
    printf '%s\t%s\n' "$inst" "$method"
  done
done > "$OUT_DIR/jobs.txt"

# Launch — {1}=instance, {2}=method
parallel -j 120 --colsep '\t' --joblog "$OUT_DIR/joblog.tsv" --resume --halt never \
  /bin/bash "$WORKER" {1} {2} "$OUT_DIR" "$CBC_BIN" "$REF_TSV" \
  :::: "$OUT_DIR/jobs.txt"
```

**Note:** `--resume` requires `--joblog`. The joblog must remain at the same path
between the original run and the resume. It identifies jobs by Seq number (position in
the input), so the jobs file must also remain unchanged.

---

## ASan (AddressSanitizer) Experiments

Use the debug build at `~/prog/cbc-dbg/bin/cbc` (built with `-O1 -g -fsanitize=address`).

Set ASAN_OPTIONS **before** launching xargs so every child inherits it:
```bash
export ASAN_OPTIONS="halt_on_error=0:log_path=$OUT_DIR/asan_logs/asan"
mkdir -p "$OUT_DIR/asan_logs"
```
- `halt_on_error=0` — process continues after an ASan error (don't lose other results)
- `log_path=<prefix>` — ASan writes errors to `<prefix>.<pid>` files, separate from
  CBC's stdout/stderr
- CBC stdout+stderr → log file via `> "$LOG" 2>&1` in the worker

Check for ASan errors after the run:
```bash
find "$OUT_DIR/asan_logs" -name 'asan_*' | wc -l   # should be 0
```

---

## Capturing Errors

Both ASan errors and CBC assertion failures must be captured:

| Error type | Where it goes | How to capture |
|---|---|---|
| ASan memory error | `$ASAN_OPTIONS log_path` files | Automatic via ASAN_OPTIONS |
| C assert failure | stderr | Redirect `2>&1` in worker |
| CBC error message | stdout | Redirect `> log 2>&1` in worker |
| xargs meta-errors | orchestrator stdout | Captured in `/tmp/exp.log` |

In the worker, always do:
```bash
timeout $TIMELIMIT /bin/bash -c \
  "$CBC_BIN $INSTANCE_PATH $OPTS -quit > $LOG 2>&1"
```
(`timeout` sends SIGTERM at limit, logs the kill; CBC exits non-zero → worker catches it)

---

## Timing

```bash
START=$(date +%s.%N)
# ... run CBC ...
END=$(date +%s.%N)
ELAPSED=$(awk -v s="$START" -v e="$END" 'BEGIN {printf "%.3f", e-s}')
```

---

## Monitoring a Running Experiment

```bash
OUT_DIR=/path/to/exp_results/myexp_YYYYMMDD_HHMMSS

# Progress
ls "$OUT_DIR/results/" | wc -l   # completed jobs

# Active processes
ps aux | grep 'cbc-dbg' | grep -v grep | wc -l

# ASan errors so far
find "$OUT_DIR/asan_logs" -name 'asan_*' | wc -l

# Match breakdown
awk -F'\t' 'NR>1 {print $6}' "$OUT_DIR/results/"*.tsv | sort | uniq -c

# Mismatches detail
awk -F'\t' 'NR>1 && $6 != "true"' "$OUT_DIR/results/"*.tsv
```

---

## Quick Reference

| Thing | Command/Path |
|---|---|
| Opt CBC binary | `~/prog/cbc/bin/cbc` |
| Debug+ASan CBC binary | `~/prog/cbc-dbg/bin/cbc` |
| Reference LP objectives | `~/inst/super/lp_relaxation_results_gurobi.tsv` |
| Per-instance time limits | `~/inst/super/suggested_lp_time_limits.tsv` |
| Parallelism to use | **120 cores** |
| Orchestrator log | `/tmp/<expname>.log` |
| GNU parallel version | 20231122 (`parallel --version`) |

---

## What NOT to Use

- **`export -f funcname`** + parallel/xargs — function exports do NOT survive into
  child bash processes. Always use a separate `.sh` file for the worker.
- **`~` in scripts** — use `$HOME` instead; tilde is not expanded in all contexts.
- **`dirname "$0"`** alone — unreliable when script is called as `bash /path/to/script.sh`;
  use `dirname "$(realpath "$0")"`.
