#!/usr/bin/env bash
# run_experiments.sh
#
# Run CBC on MIP instances, verify objective values, and collect per-instance
# error logs from ASan / TSan / Valgrind.  Up to --parallel N instances run
# concurrently via GNU parallel.
#
# Optionally builds the cbcbox debug binary (with or without ASan) before
# running, so one command covers build + test.
#
# Output layout ($OUTDIR/):
#   summary.tsv              — one row per instance: status, obj, elapsed, gap, ...
#   <instance>.log           — full CBC + sanitizer stdout/stderr
#   <instance>.sol           — solution file (only present when a feasible solution was found)
#   <instance>.validate.log  — solution validation output (mipster_validate_sol)
#   stats.csv                — CBC statistics for all instances when supported
#   <instance>.vg.log        — valgrind report (when --valgrind is used)
#   <instance>.err           — extracted error snippet (sanitizer / wrong obj / crash / invalid sol)
#   <instance>.zh-constraints.csv
#                           — ZeroHalf best/worst constraint report (--zh-constraint-report)
#   parallel.log             — GNU parallel job log (timing, exit codes)
#   report.txt               — human-readable summary
#   memory_usage.tsv         — memory utilisation sampled every 5 min during the run
#
# Status values in summary.tsv:
#   SOLVED          solved, objective matches an exact reference
#   SOLVED(best)    solved, objective matches the best-known reference
#   SOLVED(improved) solved, objective improves on the best-known reference
#   SOLVED(no_ref)  solved, no reference objective available
#   INVALID_SOL     CBC found a solution but validation (mipster_validate_sol) failed
#   WRONG_OBJ       solved but objective is worse than the reference
#   NO_SOLUTION     CBC exited normally but found no integer feasible solution
#   TIMEOUT       reached --timelimit; CBC exited on its own (objective accepted)
#   OVERTIME      CBC ignored its time limit and was killed after --overtime-grace
#   CRASH         non-zero exit, no sanitizer markers detected
#   ASAN_ERROR    AddressSanitizer error found in log
#   TSAN_RACE     ThreadSanitizer race found in log
#   VG_ERROR      valgrind reported errors (--error-exitcode=42 triggered)
#
# Usage:
#   ./run_experiments.sh [OPTIONS]
#
# Build options (mutually exclusive; sets --bin automatically):
#   --build-asan         Build cbcbox debug+ASan binary, then run tests
#   --build              Build cbcbox debug binary (no sanitizer), then run tests
#   --cbcbox-dir DIR     Path to cbcbox repo (default: /home/haroldo/dev/cbcbox)
#
# Required (unless --build or --build-asan is given):
#   --bin PATH           Path to cbc binary
#
# Optional:
#   --parallel N         Instances to run concurrently (default: 10)
#   --threads N          CBC -threads flag per instance (default: 1)
#   --timelimit T        CBC time limit passed via -seconds (default: 3600)
#   --overtime-grace G   Seconds to wait after timelimit before hard-killing CBC
#                        (default: 900).  Instance is tagged OVERTIME if killed.
#   --outdir DIR         Output directory (default: auto-named under exp_results/)
#   --instances DIR      Directory with .mps.gz files
#                        (default: /home/haroldo/inst/miplib/2017)
#   --valgrind TOOL      Wrap with valgrind TOOL: helgrind | memcheck | drd
#   --vg-supp FILE       Extra --suppressions=FILE passed to valgrind (repeatable)
#   --env KEY=VAL        Extra environment variable for each run (repeatable)
#   --opts FLAG          Extra CBC flag before -solve (repeatable)
#                        e.g. --opts "-rens 1000" --opts "-rounding both"
#   --zh-constraint-report
#                       Save per-instance ZeroHalf best/worst constraint CSV
#                       report to <outdir>/<instance>.zh-constraints.csv
#   --dry-run            Print commands without executing them
#
# Examples:
#   # Build ASan binary, then test 10 instances at a time on MIPLIB 2017:
#   ./run_experiments.sh --build-asan --parallel 10
#
#   # Use a pre-built ASan binary, override instances dir:
#   ASAN_OPTIONS="detect_leaks=1:abort_on_error=1" \
#     ./run_experiments.sh \
#       --bin /home/haroldo/dev/cbcbox/cbc_dist_debug_avx2/bin/cbc \
#       --parallel 10
#
#   # Test specific heuristics (Setting 1: RENS + Rounding):
#   ./run_experiments.sh --build-asan --parallel 10 \
#     --opts "-rens 1000" --opts "-rounding both"
#
#   # Helgrind, 3 CBC threads, 2 instances at a time (very slow):
#   ./run_experiments.sh --bin ~/prog/cbc-hwpar-vg/bin/cbc \
#     --threads 3 --parallel 2 \
#     --valgrind helgrind \
#     --vg-supp "$(pwd)/helgrind_openblas.supp"

set -euo pipefail
# Ignore SIGHUP so GNU parallel workers survive terminal disconnection.
# Child processes inherit SIG_IGN and are therefore also protected.
trap '' HUP

TIMEOUT_CMD="timeout"
GREP_CMD="grep"

if [[ "$(uname -s)" == "Darwin" ]]; then
  if ! command -v timeout &>/dev/null && ! command -v gtimeout &>/dev/null; then
    echo "Error: 'timeout' or 'gtimeout' is required. On macOS, run: brew install coreutils" >&2
    exit 1
  fi
  command -v timeout &>/dev/null || TIMEOUT_CMD="gtimeout"
  command -v ggrep &>/dev/null && GREP_CMD="ggrep"
fi

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(dirname "$SCRIPT_DIR")"

# Flat monorepo — project dir is always the repo root.
_project_dir() {
  echo "$REPO_DIR"
}

_project_config_status() {
  echo "$(_project_dir "$1")/config.status"
}

_configure_invocation_from_status() {
  local cs="$1"
  local sq="'"
  local cfg
  cfg=$(grep "^ac_cs_config=" "$cs" 2>/dev/null \
    | sed "s/^ac_cs_config=${sq}//; s/${sq}$//" \
    | sed "s/${sq}\\\\${sq}${sq}/${sq}/g" || true)
  if [[ -n "$cfg" ]]; then
    printf '%s\n' "$cfg"
    return
  fi
  awk '
    /^[[:space:]]*exec .*\.\/configure / {
      line = $0
      sub(/^[[:space:]]*exec .*\.\/configure /, "", line)
      sub(/ \$ac_configure_extra_args --no-create --no-recursion$/, "", line)
      print line
      exit
    }' "$cs"
}

_subst_from_status() {
  local cs="$1" key="$2"
  awk -v k="$key" '
    $0 ~ "^S\\[\"" k "\"\\]=" {
      line = $0
      sub("^S\\[\"" k "\"\\]=\"", "", line)
      if (line ~ /"\\$/) {
        sub(/"\\$/, "", line)
        if ((getline nxt) > 0) {
          sub(/^"/, "", nxt)
          sub(/"$/, "", nxt)
          line = line nxt
        }
      } else {
        sub(/"$/, "", line)
      }
      print line
      exit
    }
    $0 ~ "^s,@" k "@," {
      line = $0
      sub("^s,@" k "@,", "", line)
      sub(/,;t t$/, "", line)
      print line
      exit
    }' "$cs"
}

# ── Hardcoded reference values (override instance metadata when present) ─────
# These are verified-correct objectives used primarily for race-detection tests.
declare -A EXPECTED_OVERRIDE=(
  [air03]=340160.0
  [air04]=56137.0
  [air05]=26374.0
  [gesa2-o]=25779856.3717
  [j3050_8]=1.0
  [mas76]=40005.054142
  [mzzv11]=-21718.0
  [neos-1281048]=601.0
  [neos-810286]=2877.0
  [nursesched-sprint02]=58.0
  [nw04]=16862.0
  [pk1]=11.0
  [pp08a]=7350.0
  [qiu]=-132.873136947
  [sprint_hidden06_j]=130.0
  [stein45]=30.0
  [trd445c]=-153419.078836
  # MIPster proves 202.185252817 optimal (gap=0); MIPLIB reference 202.190074193 is stale
  [dt_optimization]=202.185252817
)
REL_TOL=1e-6
ABS_TOL=1e-4

# ── Help ─────────────────────────────────────────────────────────────────────
show_help() {
  cat << 'EOF'
Usage: ./run_experiments.sh [OPTIONS]

Run CBC on MIP instances, verify objective values, and collect per-instance
error logs from ASan / TSan / Valgrind.  Up to --parallel N instances run
concurrently via GNU parallel.

Output layout (<outdir>/):
  summary.tsv              — one row per instance: status, obj, elapsed, gap
  <instance>.log           — full CBC + sanitizer stdout/stderr
  <instance>.sol           — solution file (only present when a feasible solution was found)
  <instance>.validate.log  — solution validation output (mipster_validate_sol)
  stats.csv                — CBC statistics for all instances when supported
  <instance>.vg.log        — valgrind report (when --valgrind is used)
  <instance>.err           — extracted error snippet (sanitizer / wrong obj / crash / invalid sol)
  <instance>.zh-constraints.csv
                           — ZeroHalf best/worst constraint report (--zh-constraint-report)
  parallel.log             — GNU parallel job log (timing, exit codes)
  report.txt               — human-readable summary
  memory_usage.tsv         — memory utilisation sampled every 5 min during the run

Status values in summary.tsv:
  SOLVED            solved, objective matches an exact reference
  SOLVED(best)      solved, objective matches the best-known reference
  SOLVED(improved)  solved, objective improves on the best-known reference
  SOLVED(no_ref)    solved, no reference objective available
  SOLVED(inf)       CBC proved infeasible; benchmark agrees (or no ref)
  INVALID_SOL       CBC found a solution but validation (mipster_validate_sol) failed
  WRONG_OBJ         solved but objective is worse than the reference
  INFEASIBLE_WRONG  CBC proved infeasible but a numeric objective was expected
  NO_SOLUTION       CBC exited normally but found no integer feasible solution
  TIMEOUT           reached --timelimit; CBC exited on its own
  TIMEOUT(no_sol)   timed out without finding any feasible solution
  TIMEOUT(gap=X%)   timed out with a feasible solution, gap X%
  OVERTIME          CBC ignored its time limit and was killed after --overtime-grace
  CRASH             non-zero exit, no sanitizer markers detected
  ASAN_ERROR        AddressSanitizer error found in log
  TSAN_RACE         ThreadSanitizer race found in log
  VG_ERROR          valgrind reported errors (--error-exitcode=42 triggered)

Build options (mutually exclusive; sets --bin automatically):
  --build-asan         Build cbcbox debug+ASan binary, then run tests
  --build              Build cbcbox debug binary (no sanitizer), then run tests
  --cbcbox-dir DIR     Path to cbcbox repo (default: /home/haroldo/dev/cbcbox)

Required (unless --build or --build-asan is given):
  --bin PATH           Path to cbc binary

Optional:
  --parallel N         Instances to run concurrently (default: 10)
  --threads N          CBC -threads flag per instance (default: 1)
  --timelimit T        CBC time limit in seconds
                       (default: 3600 normally, 36000 when --valgrind is used)
  --overtime-grace G   Seconds to wait after timelimit before hard-killing CBC
                       (default: 900).  Instance is tagged OVERTIME if killed.
  --outdir DIR         Output directory (default: auto-named under exp_results/)
  --instances DIR      Directory with .mps.gz files
                       (default: /home/haroldo/inst/miplib/2017)
  --valgrind TOOL      Wrap with valgrind TOOL: helgrind | memcheck | drd
  --vg-supp FILE       Extra --suppressions=FILE passed to valgrind (repeatable)
  --env KEY=VAL        Extra environment variable for each run (repeatable)
  --opts FLAG          Extra CBC flag before -solve (repeatable)
                       e.g. --opts "-rens 1000" --opts "-rounding both"
  --zh-constraint-report
                       Save per-instance ZeroHalf best/worst constraint CSV
                       report to <outdir>/<instance>.zh-constraints.csv
  --restart            Resume an interrupted run: skip instances whose .result
                       file already exists in --outdir and only run the rest.
                       Requires --outdir to point at the previous output dir.
  --dry-run            Print commands without executing them
  -h, --help           Show this help message and exit

Examples:
  # Build ASan binary, then test 10 instances at a time on MIPLIB 2017:
  ./run_experiments.sh --build-asan --parallel 10

  # Use a pre-built binary, override instances dir:
  ./run_experiments.sh \
    --bin ~/prog/cbc/bin/cbc \
    --instances /home/haroldo/inst/miplib/2017 \
    --parallel 20

  # Use a pre-built ASan binary with leak detection:
  ASAN_OPTIONS="detect_leaks=1:abort_on_error=1" \
    ./run_experiments.sh \
      --bin /home/haroldo/dev/cbcbox/cbc_dist_debug_avx2/bin/cbc \
      --parallel 10

  # Test specific heuristics (Setting 1: RENS + Rounding):
  ./run_experiments.sh --build-asan --parallel 10 \
    --opts "-rens 1000" --opts "-rounding both"

  # Helgrind, 3 CBC threads, 2 instances at a time (very slow):
  ./run_experiments.sh --bin ~/prog/cbc-hwpar-vg/bin/cbc \
    --threads 3 --parallel 2 \
    --valgrind helgrind \
    --vg-supp "$(pwd)/helgrind_openblas.supp"

  # Dry run to preview the commands that would be executed:
  ./run_experiments.sh --bin ~/prog/cbc/bin/cbc --dry-run

  # Resume an interrupted run (reuses existing .result files):
  ./run_experiments.sh --bin ~/prog/cbc/bin/mipster \
    --outdir /home/haroldo/experiments/cbc/opt_06_05 \
    --restart
EOF
}

# ── Defaults ──────────────────────────────────────────────────────────────────
CBC_BIN=""
PARALLEL=10
THREADS=1
TIMELIMIT=""
OVERTIME_GRACE=900        # seconds to wait beyond timelimit before hard-killing
INSTANCES_DIR="/home/haroldo/inst/miplib/2017"
OUTDIR=""
VG_TOOL=""
VG_SUPPS=()
EXTRA_ENV=()
CBC_EXTRA_OPTS=()
DRY_RUN=0
RESTART=0
BUILD_SANITIZER=""        # "asan" | "none" | "" (no build)
CBCBOX_DIR="/home/haroldo/dev/cbcbox"
ZH_CONSTRAINT_REPORT=0
CBC_STATS_MODE="auto"      # auto | writeStat | none

# ── Parse arguments ───────────────────────────────────────────────────────────
while [[ $# -gt 0 ]]; do
  case "$1" in
    --bin)            CBC_BIN="$2";             shift 2 ;;
    --parallel)       PARALLEL="$2";            shift 2 ;;
    --threads)        THREADS="$2";             shift 2 ;;
    --timelimit)      TIMELIMIT="$2";           shift 2 ;;
    --overtime-grace) OVERTIME_GRACE="$2";      shift 2 ;;
    --outdir)         OUTDIR="$2";              shift 2 ;;
    --instances)      INSTANCES_DIR="$2";       shift 2 ;;
    --valgrind)       VG_TOOL="$2";             shift 2 ;;
    --vg-supp)        VG_SUPPS+=("$2");         shift 2 ;;
    --env)            EXTRA_ENV+=("$2");        shift 2 ;;
    --opts)           CBC_EXTRA_OPTS+=("$2");   shift 2 ;;
    --cbcbox-dir)     CBCBOX_DIR="$2";          shift 2 ;;
    --build-asan)     BUILD_SANITIZER="asan";   shift   ;;
    --build)          BUILD_SANITIZER="none";   shift   ;;
    --zh-constraint-report) ZH_CONSTRAINT_REPORT=1; shift ;;
    --restart)          RESTART=1;                 shift   ;;
    --dry-run)          DRY_RUN=1;                 shift   ;;
    -h|--help)        show_help; exit 0                 ;;
    *) echo "Unknown option: $1" >&2; exit 1  ;;
  esac
done

# ── Build cbcbox debug binary if requested ────────────────────────────────────
if [[ -n "$BUILD_SANITIZER" ]]; then
  [[ -d "$CBCBOX_DIR" ]] || { echo "Error: cbcbox dir not found: $CBCBOX_DIR" >&2; exit 1; }
  build_script="$CBCBOX_DIR/scripts/build_debug.sh"
  [[ -x "$build_script" ]] || { echo "Error: build script not found: $build_script" >&2; exit 1; }

  echo "==> Building cbcbox debug binary in $CBCBOX_DIR ..."
  build_args=()
  [[ "$BUILD_SANITIZER" == "asan" ]] && build_args+=("--asan")
  (cd "$CBCBOX_DIR" && bash "$build_script" "${build_args[@]}")

  # Detect output binary (arch-dependent path from build_debug.sh)
  case "$(uname -m)" in
    x86_64|amd64) built_bin="$CBCBOX_DIR/cbc_dist_debug_avx2/bin/cbc" ;;
    arm64|aarch64) built_bin="$CBCBOX_DIR/cbc_dist_debug/bin/cbc" ;;
    *) echo "Error: unsupported arch $(uname -m)" >&2; exit 1 ;;
  esac
  [[ -x "$built_bin" ]] || { echo "Error: built binary not found: $built_bin" >&2; exit 1; }
  echo "==> Build complete: $built_bin"
  echo ""
  CBC_BIN="$built_bin"
fi

[[ -z "$CBC_BIN" ]] && { echo "Error: --bin is required (or use --build / --build-asan)" >&2; exit 1; }
[[ -x "$CBC_BIN" ]] || { echo "Error: $CBC_BIN is not executable" >&2; exit 1; }

# ── Snapshot binary to a temp dir so rebuilds during the run don't interfere ──
# Exception: system-installed binaries (/usr/bin, /usr/local/bin) are used
# directly — no rebuild can affect them and they rely on sibling variant
# binaries + system library paths that must stay in place.
REAL_CBC_BIN="$(realpath "$CBC_BIN")"
BIN_DIR="$(dirname "$REAL_CBC_BIN")"

if [[ "$BIN_DIR" == /usr/bin || "$BIN_DIR" == /usr/local/bin ]]; then
  TMP_CBC_BIN="$REAL_CBC_BIN"
  echo "==> System binary — using in place: $TMP_CBC_BIN"
else
  EXP_TMPDIR=$(mktemp -d /tmp/cbc_exp_XXXXXXXX)
  trap 'rm -rf "${EXP_TMPDIR:-}"' EXIT
  mkdir -p "$EXP_TMPDIR/bin"

  # Copy launcher + variant binaries
  cp "$REAL_CBC_BIN" "$EXP_TMPDIR/bin/mipster"
  chmod +x "$EXP_TMPDIR/bin/mipster"
  for vbin in "$BIN_DIR"/mipster-generic "$BIN_DIR"/mipster-neon \
              "$BIN_DIR"/mipster-avx2   "$BIN_DIR"/mipster-haswell; do
    [[ -x "$vbin" ]] || continue
    cp "$vbin" "$EXP_TMPDIR/bin/"
    echo "    also copied variant: $(basename "$vbin")"
  done

  # Mirror lib dirs so $ORIGIN/../lib/<variant> RPATHs resolve correctly
  ORIG_LIB_DIR="$(realpath "$BIN_DIR/../lib" 2>/dev/null || true)"
  if [[ -d "$ORIG_LIB_DIR" ]]; then
    for variant in generic neon avx2 haswell; do
      [[ -d "$ORIG_LIB_DIR/$variant" ]] || continue
      mkdir -p "$EXP_TMPDIR/lib"
      cp -rP "$ORIG_LIB_DIR/$variant" "$EXP_TMPDIR/lib/"
      echo "    also copied lib: $variant"
    done
  fi

  TMP_CBC_BIN="$EXP_TMPDIR/bin/mipster"
  echo "==> Binary snapshotted: $TMP_CBC_BIN"
  echo "    (original: $REAL_CBC_BIN)"
fi

ORIGINAL_CBC_BIN="$CBC_BIN"
CBC_BIN="$TMP_CBC_BIN"
echo ""

# ── Detect optional solver features by querying the actual binary ─────────────
if [[ "$CBC_STATS_MODE" == "auto" ]]; then
  if "$CBC_BIN" -? 2>&1 | grep -qi 'writeStat'; then
    CBC_STATS_MODE="writeStat"
  else
    CBC_STATS_MODE="none"
  fi
fi

if [[ "$CBC_STATS_MODE" == "writeStat" ]]; then
  echo "==> Solver stats export: writeStat"
else
  echo "==> Solver stats export: unavailable on this CBC binary; using log parsing only"
fi
echo ""

# ── Resolve time limit default ────────────────────────────────────────────────
# Valgrind slows CBC by 20-100x; use 10 h per instance so slow instances still
# finish.  Non-valgrind runs default to 1 h.
if [[ -z "$TIMELIMIT" ]]; then
  if [[ -n "$VG_TOOL" ]]; then
    TIMELIMIT=36000   # 10 hours
  else
    TIMELIMIT=3600    # 1 hour
  fi
fi
WALLCLOCK_LIMIT=$((TIMELIMIT + OVERTIME_GRACE))  # hard kill deadline

# ── Output directory ──────────────────────────────────────────────────────────
if [[ -z "$OUTDIR" ]]; then
  VARIANT=$(basename "$(dirname "$(dirname "$CBC_BIN")")")  # e.g. cbc-hwpar-tsan
  TS=$(date +%Y%m%d_%H%M%S)
  OUTDIR="exp_results/${VARIANT}_t${THREADS}_p${PARALLEL}_${TS}"
  [[ -n "$VG_TOOL" ]] && OUTDIR="${OUTDIR}_${VG_TOOL}"
fi
mkdir -p "$OUTDIR"

SUMMARY="$OUTDIR/summary.tsv"
REPORT="$OUTDIR/report.txt"

# ── Build expected.sh from instance metadata + hardcoded overrides ────────────
# Workers source this file to get:
#   EXPECTED   — numeric reference objective / best-known value
#   INFEASIBLE — instances known to be infeasible
#   REF_KIND   — "exact" for proven optima, "best" for best-known incumbents
#   REF_SENSE  — "min" / "max" when inferable from the associated bound
EXPECTED_FILE="$OUTDIR/expected.sh"
python3 - "$INSTANCES_DIR" "$EXPECTED_FILE" \
  "${!EXPECTED_OVERRIDE[@]}" \
  "${EXPECTED_OVERRIDE[@]}" << 'PYEOF'
import csv
import math
import os
import sys

instances_dir = sys.argv[1]
out_file = sys.argv[2]

# Hardcoded overrides are passed as two parallel lists: keys then values.
# sys.argv = [script, instances_dir, out, k0, k1, ..., kN, v0, v1, ..., vN]
rest = sys.argv[3:]
n_keys = len(rest) // 2
override_keys = rest[:n_keys]
override_values = rest[n_keys:]
hardcoded = {k: float(v) for k, v in zip(override_keys, override_values)}

expected: dict[str, float] = {}
infeasible: set[str] = set()
ref_kind: dict[str, str] = {}
ref_sense: dict[str, str] = {}


def parse_float(value):
    if value is None:
        return None
    text = str(value).strip()
    if not text:
        return None
    try:
        number = float(text)
    except ValueError:
        return None
    return None if math.isnan(number) else number


def infer_sense(reference, bound):
    if reference is None or bound is None:
        return None
    tol = 1e-9 * max(1.0, abs(reference), abs(bound))
    if bound < reference - tol:
        return "min"
    if bound > reference + tol:
        return "max"
    return None


def first_not_none(*values):
    for value in values:
        if value is not None:
            return value
    return None


def load_bks(path):
    refs = {}
    if not os.path.isfile(path):
        return refs
    with open(path, newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            name = (row.get("instance") or "").strip()
            value = parse_float(row.get("bks"))
            if name and value is not None:
                refs[name] = value
    return refs


def load_solution_metadata(path, bks_lookup):
    if not os.path.isfile(path):
        return
    with open(path, newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            name = (row.get("instance") or "").strip()
            if not name:
                continue

            miplib_status = (row.get("miplib_solu_status") or "").strip().lower()
            super_status = (row.get("super_status") or "").strip().lower()

            if miplib_status == "inf" or super_status == "infeasible":
                infeasible.add(name)
                ref_kind.setdefault(name, "infeasible")
                continue

            reference = first_not_none(
                parse_float(row.get("miplib_solu_value")),
                parse_float(row.get("super_obj")),
                parse_float(row.get("cbc_bks")),
                bks_lookup.get(name),
            )
            if reference is None:
                continue

            exact = miplib_status == "opt" or super_status == "optimal"
            kind = "exact" if exact else "best"

            bound = first_not_none(
                parse_float(row.get("super_bound")),
                parse_float(row.get("cbc_lp")),
                parse_float(row.get("relax_obj")),
            )
            sense = infer_sense(reference, bound)

            expected[name] = reference
            ref_kind[name] = kind
            if sense is not None:
                ref_sense[name] = sense


def load_benchmark_csv(path):
    if not os.path.isfile(path):
        return
    with open(path, newline="") as f:
        reader = csv.reader(f)
        next(reader, None)
        for row in reader:
            if len(row) < 11:
                continue
            name = row[0].strip('"')
            obj = row[10].strip('"')
            if not name:
                continue
            if obj == "Infeasible":
                infeasible.add(name)
                ref_kind.setdefault(name, "infeasible")
                continue
            value = parse_float(obj)
            if value is None:
                continue
            expected.setdefault(name, value)
            ref_kind.setdefault(name, "exact")


bks_lookup = load_bks(os.path.join(instances_dir, "bks.tsv"))
load_solution_metadata(os.path.join(instances_dir, "solution_metadata.csv"), bks_lookup)
load_benchmark_csv(os.path.join(instances_dir, "benchmark.csv"))

# Any bks.tsv entry not already covered by solution_metadata or benchmark is added
# as a best-known (not proven-optimal) reference.
for name, value in bks_lookup.items():
    expected.setdefault(name, value)
    ref_kind.setdefault(name, "best")

# Hardcoded values take priority and are treated as exact references.
for name, value in hardcoded.items():
    expected[name] = value
    ref_kind[name] = "exact"

with open(out_file, "w") as f:
    f.write("declare -A EXPECTED=(\n")
    for key in sorted(expected):
        f.write(f"  [{key}]={repr(expected[key])}\n")
    f.write(")\n")
    f.write("declare -A INFEASIBLE=(\n")
    for key in sorted(infeasible):
        f.write(f"  [{key}]=1\n")
    f.write(")\n")
    f.write("declare -A REF_KIND=(\n")
    for key in sorted(ref_kind):
        f.write(f"  [{key}]={ref_kind[key]}\n")
    f.write(")\n")
    f.write("declare -A REF_SENSE=(\n")
    for key in sorted(ref_sense):
        f.write(f"  [{key}]={ref_sense[key]}\n")
    f.write(")\n")
PYEOF

# ── Serialise VG_SUPPS (array) into a colon-separated string for export ──────
VG_SUPPS_STR=$(IFS=:; echo "${VG_SUPPS[*]+"${VG_SUPPS[*]}"}")

# ── Serialise EXTRA_ENV into newline-separated string for export ──────────────
EXTRA_ENV_STR=$(printf '%s\n' "${EXTRA_ENV[@]+"${EXTRA_ENV[@]}"}")

# ── Serialise CBC_EXTRA_OPTS into newline-separated string for export ─────────
CBC_EXTRA_OPTS_STR=$(printf '%s\n' "${CBC_EXTRA_OPTS[@]+"${CBC_EXTRA_OPTS[@]}"}")

# ── Colour support detection ──────────────────────────────────────────────────
# Enabled only when stdout is a TTY and NO_COLOR is unset (https://no-color.org)
if [[ -t 1 && -z "${NO_COLOR:-}" ]]; then
  USE_COLOR=1
else
  USE_COLOR=0
fi

# ── Export everything run_instance needs ─────────────────────────────────────
export CBC_BIN OUTDIR THREADS TIMELIMIT OVERTIME_GRACE WALLCLOCK_LIMIT TIMEOUT_CMD GREP_CMD
if [[ "$(uname -s)" == "Darwin" ]]; then
  export PARALLEL_SHELL="$BASH"
fi
export VG_TOOL VG_SUPPS_STR EXTRA_ENV_STR CBC_EXTRA_OPTS_STR
export EXPECTED_FILE REL_TOL ABS_TOL DRY_RUN RESTART SCRIPT_DIR ZH_CONSTRAINT_REPORT
export USE_COLOR CBC_STATS_MODE

# ── ASan environment defaults (only when no ASAN_OPTIONS already set) ─────────
# abort_on_error: ensures CBC exits non-zero so we detect the error in exit_code
# detect_leaks:   catches heap leaks (requires lsan to be linked, which the
#                 ASan build includes)
# halt_on_error:  stop after first error to keep logs readable
if [[ "$BUILD_SANITIZER" == "asan" && -z "${ASAN_OPTIONS:-}" ]]; then
  export ASAN_OPTIONS="abort_on_error=1:detect_leaks=1:halt_on_error=1:print_stats=1"
  echo "==> ASAN_OPTIONS set to: $ASAN_OPTIONS"
fi

# ── Helper: check objective within tolerance ──────────────────────────────────
obj_ok() {
  local got="$1" ref="$2"
  awk -v got="$got" -v ref="$ref" \
      -v rtol="$REL_TOL" -v atol="$ABS_TOL" \
    'BEGIN {
       diff = got - ref; if (diff < 0) diff = -diff
       absref = ref < 0 ? -ref : ref
       tol = rtol * absref; if (atol > tol) tol = atol
       exit (diff <= tol) ? 0 : 1
     }'
}
export -f obj_ok

ref_allows_obj() {
  local got="$1" ref="$2" kind="${3:-exact}" sense="${4:-}"
  if [[ "$kind" == "best" ]]; then
    case "$sense" in
      min)
        awk -v got="$got" -v ref="$ref" -v rtol="$REL_TOL" -v atol="$ABS_TOL" '
          BEGIN {
            tol = rtol * ((ref < 0) ? -ref : ref); if (atol > tol) tol = atol
            exit (got <= ref + tol) ? 0 : 1
          }'
        return
        ;;
      max)
        awk -v got="$got" -v ref="$ref" -v rtol="$REL_TOL" -v atol="$ABS_TOL" '
          BEGIN {
            tol = rtol * ((ref < 0) ? -ref : ref); if (atol > tol) tol = atol
            exit (got + tol >= ref) ? 0 : 1
          }'
        return
        ;;
    esac
  fi
  obj_ok "$got" "$ref"
}
export -f ref_allows_obj

gap_vs_ref() {
  local got="$1" ref="$2" kind="${3:-exact}" sense="${4:-}"
  awk -v got="$got" -v ref="$ref" -v kind="$kind" -v sense="$sense" \
      -v rtol="$REL_TOL" -v atol="$ABS_TOL" '
    function abs(x) { return x < 0 ? -x : x }
    BEGIN {
      absref = abs(ref)
      tol = rtol * absref
      if (atol > tol) tol = atol

      if (kind == "best" && sense == "min") diff = got - ref
      else if (kind == "best" && sense == "max") diff = ref - got
      else diff = abs(got - ref)

      if (diff <= tol) {
        print "0.0%"
        exit 0
      }

      if (absref <= 1e-10) {
        print ">100%"
        exit 0
      }

      pct = diff / absref * 100.0
      if (pct > 100.0) print ">100%"
      else printf "%.1f%%\n", pct
    }'
}
export -f gap_vs_ref

# ── Helper: extract sanitizer / crash snippet into .err file ─────────────────
extract_errors() {
  local logfile="$1" errfile="$2" vg_logfile="$3"

  if grep -q "ERROR: AddressSanitizer\|ERROR: LeakSanitizer" "$logfile" 2>/dev/null; then
    { echo "=== AddressSanitizer errors ==="
      grep -A 40 "ERROR: AddressSanitizer\|ERROR: LeakSanitizer" "$logfile" | head -120
    } >> "$errfile"
  fi

  if grep -q "WARNING: ThreadSanitizer\|ERROR: ThreadSanitizer" "$logfile" 2>/dev/null; then
    { echo "=== ThreadSanitizer reports ==="
      grep -A 40 "WARNING: ThreadSanitizer\|ERROR: ThreadSanitizer" "$logfile" | head -200
    } >> "$errfile"
  fi

  if [[ -f "$vg_logfile" ]] && grep -q "ERROR SUMMARY: [^0]" "$vg_logfile" 2>/dev/null; then
    { echo "=== Valgrind ${VG_TOOL} errors ==="
      cat "$vg_logfile"
    } >> "$errfile"
  fi
}
export -f extract_errors

# ── Helper: format elapsed seconds as compact human-readable string ───────────
# Output is always right-aligned to 6 chars: "   45s", " 1m23s", " 1h02m"
fmt_time() {
  local s=$1 t
  if   (( s < 60 ));   then t=$(printf "%ds"      "$s")
  elif (( s < 3600 )); then t=$(printf "%dm%02ds"  "$((s/60))"    "$((s%60))")
  else                      t=$(printf "%dh%02dm"  "$((s/3600))"  "$(( (s%3600)/60 ))")
  fi
  printf "%6s" "$t"
}
export -f fmt_time

# ── Per-instance worker (called by GNU parallel) ──────────────────────────────
run_instance() {
  local mps="$1"
  local name
  name=$(basename "$mps" .mps.gz)

  local logfile="$OUTDIR/${name}.log"
  local errfile="$OUTDIR/${name}.err"
  local vg_logfile="$OUTDIR/${name}.vg.log"
  local resultfile="$OUTDIR/${name}.result"
  local csv_statsfile="$OUTDIR/${name}.stats.csv"
  local zh_reportfile="$OUTDIR/${name}.zh-constraints.csv"
  local solutionfile="$OUTDIR/${name}.sol"

  if [[ "$RESTART" == "1" && -f "$resultfile" ]]; then
    echo "  [skip] $name (result exists)"
    return 0
  fi

  rm -f "$csv_statsfile"

  # Reconstruct VG command from exported scalars
  local vg_cmd=()
  if [[ -n "$VG_TOOL" ]]; then
    vg_cmd=(valgrind "--tool=${VG_TOOL}" "--error-exitcode=42"
                     "--log-file=${vg_logfile}")
    if [[ "$VG_TOOL" == "helgrind" ]]; then
      vg_cmd+=("--history-level=approx" "--conflict-cache-size=1000000")
    elif [[ "$VG_TOOL" == "memcheck" ]]; then
      vg_cmd+=("--leak-check=full" "--show-leak-kinds=definite,indirect"
               "--track-origins=yes")
    fi
    if [[ -n "$VG_SUPPS_STR" ]]; then
      IFS=: read -ra _supps <<< "$VG_SUPPS_STR"
      for s in "${_supps[@]}"; do vg_cmd+=("--suppressions=$s"); done
    fi
  fi

  # Reconstruct EXTRA_ENV
  local env_assignments=()
  if [[ -n "$EXTRA_ENV_STR" ]]; then
    while IFS= read -r kv; do
      [[ -n "$kv" ]] && env_assignments+=("$kv")
    done <<< "$EXTRA_ENV_STR"
  fi
  if [[ "$ZH_CONSTRAINT_REPORT" -eq 1 ]]; then
    rm -f "$zh_reportfile"
    env_assignments+=("CGL_ZH_CONSTRAINT_PROFILE_FILE=${zh_reportfile}")
  fi

  # Reconstruct CBC_EXTRA_OPTS — each entry may contain multiple words
  # (e.g. "-DivingCoefficient both"), so word-split each line.
  local cbc_extra=()
  if [[ -n "$CBC_EXTRA_OPTS_STR" ]]; then
    while IFS= read -r opt; do
      [[ -n "$opt" ]] || continue
      read -ra words <<< "$opt"
      cbc_extra+=("${words[@]}")
    done <<< "$CBC_EXTRA_OPTS_STR"
  fi

  local cmd=("${vg_cmd[@]}"
             "$CBC_BIN" "$mps"
             "-seconds=${TIMELIMIT}")
  [[ "$THREADS" -gt 1 ]] && cmd+=("-threads=${THREADS}")
  [[ ${#cbc_extra[@]} -gt 0 ]] && cmd+=("${cbc_extra[@]}")
  cmd+=("-solve")
  if [[ "$CBC_STATS_MODE" == "writeStat" ]]; then
    cmd+=("-writeStat" "$csv_statsfile")
  fi
  cmd+=("-solu" "$solutionfile" "-quit")

  if [[ "$DRY_RUN" -eq 1 ]]; then
    if [[ ${#env_assignments[@]} -gt 0 ]]; then
      echo "DRY-RUN: env ${env_assignments[*]} ${cmd[*]}"
    else
      echo "DRY-RUN: ${cmd[*]}"
    fi
    return 0
  fi

  local start_epoch=$SECONDS
  set +e
  # Wrap with timeout: send SIGTERM at WALLCLOCK_LIMIT, then SIGKILL after 30s.
  # Exit code 124 means timeout fired (CBC ignored its -seconds limit → OVERTIME).
  if [[ ${#env_assignments[@]} -gt 0 ]]; then
    "$TIMEOUT_CMD" --kill-after=30 "$WALLCLOCK_LIMIT" \
      env "${env_assignments[@]}" "${cmd[@]}" > "$logfile" 2>&1
  else
    "$TIMEOUT_CMD" --kill-after=30 "$WALLCLOCK_LIMIT" \
      "${cmd[@]}" > "$logfile" 2>&1
  fi
  local exit_code=$?
  set -e
  local elapsed=$((SECONDS - start_epoch))

  # ── Extract results from per-instance stats CSV (reliable) or log (fallback) ─
  # The CSV is written by mipster before exit, so it is present for all normal
  # terminations.  It is absent only when the process was killed externally
  # (OVERTIME / crash), in which case we fall back to log-line parsing.
  local obj="" timed_out=0 no_feasible_solution=0 proven_infeasible=0

  if [[ "$CBC_STATS_MODE" == "writeStat" && -f "$csv_statsfile" ]]; then
    # CSV columns (1-indexed): Name,result,integer_feasible,time,sys,elapsed,objective,...
    local _csv_data _csv_result _csv_int_feas _csv_obj
    _csv_data=$(tail -1 "$csv_statsfile")
    _csv_result=$(awk  -F',' '{print $2}' <<< "$_csv_data")
    _csv_int_feas=$(awk -F',' '{print $3}' <<< "$_csv_data")
    _csv_obj=$(awk    -F',' '{print $7}' <<< "$_csv_data")

    [[ "$_csv_int_feas" == "1" ]] && obj="$_csv_obj"
    [[ "$_csv_result" == *"time limit"* ]] && timed_out=1
    if [[ "$_csv_result" == "Problem proven infeasible" \
       || "$_csv_result" == "Linear relaxation infeasible" ]]; then
      proven_infeasible=1
    fi
    [[ "$_csv_int_feas" == "0" && $proven_infeasible -eq 0 ]] && no_feasible_solution=1
  else
    # Fallback: parse from log (process was killed/crashed before writing CSV).
    obj=$("$GREP_CMD" -iP '^optimal\b.*\bobjective value\b' "$logfile" \
          | "$GREP_CMD" -oP '[-\d.eE+]+$' | tail -1 || true)
    [[ -z "$obj" ]] && \
      obj=$("$GREP_CMD" -iP '^objective value[: ]' "$logfile" \
            | "$GREP_CMD" -oP '[-\d.eE+]+$' | tail -1 || true)
    [[ -z "$obj" ]] && \
      obj=$("$GREP_CMD" -iP '^objective:?\s' "$logfile" \
            | "$GREP_CMD" -oP '[-\d.eE+]+$' | tail -1 || true)
    grep -qi "time limit\|stopped on time" "$logfile" 2>/dev/null && timed_out=1
    grep -qi '^No feasible solution found' "$logfile" 2>/dev/null && no_feasible_solution=1
    "$GREP_CMD" -qiP "Result - (Problem proven infeasible|Linear relaxation infeasible)|Problem is infeasible|Fast preprocessing: infeasibility proved" \
      "$logfile" 2>/dev/null && proven_infeasible=1
    # Discard sentinel "no incumbent" value and explicit no-solution flag.
    if [[ -n "$obj" ]] && awk -v x="$obj" 'BEGIN { exit (x >= 1e49 || x <= -1e49) ? 0 : 1 }'; then
      obj=""
    fi
    (( no_feasible_solution )) && obj=""
  fi

  # Load expected values and infeasible set
  # shellcheck source=/dev/null
  source "$EXPECTED_FILE"
  local exp="${EXPECTED[$name]:-}"
  local exp_infeasible="${INFEASIBLE[$name]:-}"
  local ref_kind="${REF_KIND[$name]:-exact}"
  local ref_sense="${REF_SENSE[$name]:-}"

  # Extract dual bound (lower bound for minimisation; present in Result block of both
  # stable and dev formats as "Lower bound:  X").
  local dual
  dual=$("$GREP_CMD" -oP '^Lower bound:\s+\K[-\d.eE+]+' "$logfile" | tail -1 || true)
  # Discard CBC's "no bound yet" sentinel value.
  if [[ -n "$dual" ]] && awk -v x="$dual" 'BEGIN { exit (x >= 1e49 || x <= -1e49) ? 0 : 1 }'; then
    dual=""
  fi
  # For proven-optimal cases CBC sometimes omits the Lower bound line because
  # dual == primal; fill it in from obj so callers always have a dual bound.
  if [[ -z "$dual" && -n "$obj" && $timed_out -eq 0 && $proven_infeasible -eq 0 ]]; then
    dual="$obj"
  fi

  local solution_found=0
  [[ -n "$obj" ]] && solution_found=1

  # ── Validate solution file if present ────────────────────────────────────────
  local validation_failed=0
  local validation_msg=""
  if [[ -f "$solutionfile" ]]; then
    # Derive path to mipster_validate_sol from CBC_BIN (same directory)
    local validate_bin="$(dirname "$CBC_BIN")/mipster_validate_sol"
    if [[ -x "$validate_bin" ]]; then
      local validate_log="$OUTDIR/${name}.validate.log"
      set +e
      "$validate_bin" "$mps" "$solutionfile" > "$validate_log" 2>&1
      local validate_exit=$?
      set -e
      if [[ $validate_exit -eq 1 ]]; then
        # Exit code 1 means validation failed (solution is infeasible)
        validation_failed=1
        validation_msg=$(head -20 "$validate_log" | grep -E "INVALID|violation|Bound|Integrality|Row" || echo "Validation failed (see ${name}.validate.log)")
      fi
      # Exit code 0 = valid, 2 = file error (ignore file errors, they're environment issues)
    fi
  fi

  # Extract primal-dual gap for timeout cases where a feasible solution was found
  local cbc_gap="-"
  if [[ $timed_out -eq 1 && -n "$obj" ]]; then
    local raw_gap
    raw_gap=$("$GREP_CMD" -oP '^Gap:\s+\K[-\d.eE+]+' "$logfile" | tail -1 || true)
    if [[ -n "$raw_gap" ]]; then
      cbc_gap=$(awk -v g="$raw_gap" 'BEGIN {
        pct = g * 100
        if (pct > 100) print ">100%"
        else           printf "%.1f%%\n", pct
      }')
    elif [[ -n "$exp" ]]; then
      cbc_gap=$(gap_vs_ref "$obj" "$exp" "$ref_kind" "$ref_sense")
    fi
  fi

  # Classify status
  local status
  if grep -q "ERROR: AddressSanitizer\|ERROR: LeakSanitizer" "$logfile" 2>/dev/null; then
    status="ASAN_ERROR"
  elif grep -q "WARNING: ThreadSanitizer\|ERROR: ThreadSanitizer" "$logfile" 2>/dev/null; then
    status="TSAN_RACE"
  elif [[ $exit_code -eq 42 ]]; then
    status="VG_ERROR"
  elif [[ $exit_code -eq 124 ]]; then
    status="OVERTIME"
  elif [[ $exit_code -ne 0 ]]; then
    status="CRASH(exit=$exit_code)"
  elif [[ $validation_failed -eq 1 ]]; then
    status="INVALID_SOL"
  elif [[ $proven_infeasible -eq 1 ]]; then
    if [[ -n "$exp_infeasible" || -z "$exp" ]]; then
      status="SOLVED(inf)"      # CBC proved infeasible; benchmark agrees (or no ref)
    else
      status="INFEASIBLE_WRONG"       # CBC proved infeasible but we expect a numeric obj
    fi
  elif [[ -z "$obj" ]]; then
    if [[ $timed_out -eq 1 ]]; then
      status="TIMEOUT(no_sol)"        # timed out without finding any feasible solution
    else
      status="NO_SOLUTION"
    fi
  elif [[ $timed_out -eq 1 ]]; then
    if [[ "$cbc_gap" != "-" ]]; then
      status="TIMEOUT(gap=${cbc_gap})"
    else
      status="TIMEOUT"
    fi
  elif [[ -n "$exp" ]]; then
    if ref_allows_obj "$obj" "$exp" "$ref_kind" "$ref_sense"; then
      if [[ "$ref_kind" == "best" ]]; then
        if obj_ok "$obj" "$exp"; then
          status="SOLVED(best)"
        else
          status="SOLVED(improved)"
        fi
      else
        status="SOLVED"
      fi
    else
      status="WRONG_OBJ"
    fi
  else
    status="SOLVED(no_ref)"
  fi

  # Extract error snippets into .err (only created when non-empty)
  > "$errfile"
  if [[ "$status" == "INVALID_SOL" ]]; then
    { echo "=== Solution validation failed ==="
      echo "  CBC reported: $obj"
      [[ -n "$validation_msg" ]] && echo "$validation_msg"
      echo ""
      echo "  See ${name}.validate.log for full validation output"
      echo ""
    } >> "$errfile"
  fi
  if [[ "$status" == "WRONG_OBJ" ]]; then
    { echo "=== Objective worse than reference ==="
      echo "  got:      $obj"
      echo "  expected: $exp"
      echo "  ref kind: $ref_kind"
      [[ -n "$ref_sense" ]] && echo "  sense:    $ref_sense"
      echo ""
    } >> "$errfile"
  fi
  if [[ "$status" == "OVERTIME" ]]; then
    { echo "=== CBC ignored -seconds time limit ==="
      printf "  Timelimit:  %ds\n" "$TIMELIMIT"
      printf "  Wallclock:  %ds (killed after %ds overtime)\n" "$elapsed" "$OVERTIME_GRACE"
      echo ""
    } >> "$errfile"
  fi
  if [[ "$status" =~ ^(ASAN_ERROR|TSAN_RACE|VG_ERROR|CRASH) ]]; then
    { echo "=== Exit code: $exit_code ==="; echo ""; } >> "$errfile"
  fi
  extract_errors "$logfile" "$errfile" "$vg_logfile"
  [[ ! -s "$errfile" ]] && rm -f "$errfile"

  # Write machine-readable result for summary assembly
  printf "%s\t%s\t%s\t%s\t%d\t%s\t%s\t%d\t%d\t%d\n" \
    "$status" "${obj:--}" "${dual:--}" "${exp:--}" "$elapsed" "$THREADS" "$cbc_gap" \
    "$solution_found" "$proven_infeasible" "$timed_out" > "$resultfile"

  if [[ "$ZH_CONSTRAINT_REPORT" -eq 1 && ! -s "$zh_reportfile" ]]; then
    rm -f "$zh_reportfile"
  fi

  # One-line console output (printed atomically, interleaves safely under parallel)
  #
  # Layout:  ICON  NAME(38)  STATUS(20)  OBJ(16)  wall=TIME(6)  [EXTRA]
  #   ✓  air05                                   SOLVED               26374.0          wall=  45s
  #   ⚠  chromaticindex1024-7                    TIMEOUT(gap=5.2%)    42.0             wall=1h00m
  #   ✗  neos-4532248-waihi                      WRONG_OBJ            1234.5           wall=1m03s  exp=1235.0

  local C_RESET='' C_BGREEN='' C_YELLOW='' C_BYELLOW='' C_BRED=''
  if [[ "${USE_COLOR:-0}" -eq 1 ]]; then
    C_RESET='\033[0m'
    C_BGREEN='\033[1;32m'
    C_YELLOW='\033[0;33m'
    C_BYELLOW='\033[1;33m'
    C_BRED='\033[1;31m'
  fi

  local icon sc
  case "$status" in
    SOLVED|SOLVED\(*)      icon="✓"; sc="$C_BGREEN"  ;;
    TIMEOUT*)              icon="⚠"; sc="$C_YELLOW"  ;;
    NO_SOLUTION|OVERTIME)  icon="⚠"; sc="$C_BYELLOW" ;;
    *)                     icon="✗"; sc="$C_BRED"    ;;
  esac

  # Truncate name to 38 chars
  local display_name="$name"
  (( ${#name} > 38 )) && display_name="${name:0:37}…"

  # Extra context: show expected value for wrong-answer cases
  local extra=""
  [[ "$status" == "INVALID_SOL" ]]      && extra="  validation failed"
  [[ "$status" == "WRONG_OBJ" ]]        && extra="  exp=${exp}"
  [[ "$status" == "INFEASIBLE_WRONG" ]] && extra="  exp=${exp}"

  local padded_status padded_obj time_str
  padded_status=$(printf "%-20s" "$status")
  padded_obj=$(printf "%-16s" "${obj:--}")
  time_str=$(fmt_time "$elapsed")

  printf "  %b%s%b %-38s  %b%s%b  %s  wall=%s%s\n" \
    "$sc" "$icon" "$C_RESET" \
    "$display_name" \
    "$sc" "$padded_status" "$C_RESET" \
    "$padded_obj" \
    "$time_str" \
    "$extra"
}
export -f run_instance

# ── Find instances ────────────────────────────────────────────────────────────
mapfile -t MPS_FILES < <(find "$INSTANCES_DIR" -maxdepth 1 -name "*.mps.gz" | sort)
[[ ${#MPS_FILES[@]} -eq 0 ]] && { echo "No .mps.gz files in $INSTANCES_DIR" >&2; exit 1; }

TOTAL=${#MPS_FILES[@]}

# ── Collect and print build environment + hardware info ───────────────────────
# Reads config.status from each COIN-OR project dir under $SCRIPT_DIR, then
# queries ldd, pkg-config, lscpu, and free(1) to describe exactly how the
# binary was built and on what machine the experiment is running.
print_env_info() {
  local bin="$1"
  local src_dir="$SCRIPT_DIR"

  echo "═══════════════════════════════════════════════════════════════════════"
  echo "  BUILD ENVIRONMENT"
  echo "───────────────────────────────────────────────────────────────────────"
  printf "  %-18s %s\n" "Source dir:" "$src_dir"
  echo ""

  # Configure invocation (best-effort extraction from config.status).
  local cs
  cs="$(_project_config_status Cbc)"
  if [[ -f "$cs" ]]; then
    local cfg
    cfg=$(_configure_invocation_from_status "$cs")
    if [[ -n "$cfg" ]]; then
      printf "  configure:\n"
      # Indent + soft-wrap at 96 chars; continuation lines get 4-space indent
      printf "    %s\n" "$cfg" \
        | fold -s -w 96 | awk 'NR==1{print} NR>1{print "    " $0}'
      echo ""
    fi
  fi

  # Compiler identity (CC / CXX resolved from Cbc's config.status)
  local cbc_cs
  cbc_cs="$(_project_config_status Cbc)"
  if [[ -f "$cbc_cs" ]]; then
    local cfg cc cxx
    cfg=$(_configure_invocation_from_status "$cbc_cs")
    cc=$(printf '%s\n' "$cfg" | "$GREP_CMD" -oP 'CC=[^'"'"' ]+' | head -1 | cut -d= -f2 || true)
    cxx=$(printf '%s\n' "$cfg" | "$GREP_CMD" -oP 'CXX=[^'"'"' ]+' | head -1 | cut -d= -f2 || true)
    echo "  Compiler versions:"
    if [[ -n "$cc" ]] && command -v "$cc" &>/dev/null; then
      printf "    CC  (%s):  %s\n" "$cc"  "$("$cc"  --version 2>&1 | head -1)"
    fi
    if [[ -n "$cxx" ]] && command -v "$cxx" &>/dev/null; then
      printf "    CXX (%s): %s\n" "$cxx" "$("$cxx" --version 2>&1 | head -1)"
    fi

    local add_cxxflags
    add_cxxflags=$(_subst_from_status "$cbc_cs" "ADD_CXXFLAGS")
    [[ -n "$add_cxxflags" ]] && printf "    ADD_CXXFLAGS:      %s\n" "$add_cxxflags"
    echo ""
  fi

  # Binary metadata
  echo "  Binary:"
  printf "    %-10s %s\n" "Path:"  "$bin"
  local bin_mtime bin_size_kb
  bin_mtime=$(stat -c "%y" "$bin" 2>/dev/null | cut -d. -f1 || stat -f "%Sm" -t "%Y-%m-%d %H:%M:%S" "$bin" 2>/dev/null || echo "?")
  bin_size_kb=$(( $(stat -c "%s" "$bin" 2>/dev/null || echo 0) / 1024 ))
  printf "    %-10s %s\n" "Built:"  "$bin_mtime"
  printf "    %-10s %d KB\n" "Size:"   "$bin_size_kb"
  echo ""

  # Linked libraries (ldd on the original installed binary)
  # Capture raw output once; reuse it for the pkg-config version filter below.
  local ldd_raw=""
  ldd_raw=$({ ldd "$bin" 2>/dev/null || otool -L "$bin" 2>/dev/null; } \
    | "$GREP_CMD" -E "(openblas|amd|lapack|blas|nauty|colamd|camd|ccolamd|cholmod|suitesparse|gomp)" \
    || true)
  echo "  Linked libraries (ldd):"
  if [[ -n "$ldd_raw" ]]; then
    awk '{
        lib=$1; path=($3 != "" ? $3 : "(not found)")
        printf "    %-40s  %s\n", lib, path
      }' <<< "$ldd_raw"
  else
    echo "    (none matched / ldd not available)"
  fi
  echo ""

  # Dependency versions via pkg-config — only for libraries actually linked.
  # Map: pkg-config name → grep pattern to test against ldd_raw.
  local -A pkglibs=([openblas]="openblas" [amd]="libamd" [lapack]="liblapack")
  local printed_header=0
  for lib in openblas amd lapack; do
    grep -q "${pkglibs[$lib]}" <<< "$ldd_raw" 2>/dev/null || continue
    local ver
    ver=$(pkg-config --modversion "$lib" 2>/dev/null || echo "(not found)")
    if [[ $printed_header -eq 0 ]]; then
      echo "  Dependency versions (pkg-config):"
      printed_header=1
    fi
    printf "    %-14s %s\n" "$lib" "$ver"
  done
  [[ $printed_header -eq 1 ]] && echo ""

  echo "═══════════════════════════════════════════════════════════════════════"
  echo "  HARDWARE"
  echo "───────────────────────────────────────────────────────────────────────"

  if [[ "$(uname -s)" == "Darwin" ]]; then
    if command -v sysctl &>/dev/null; then
      local arch model cpus cores
      arch=$(sysctl -n hw.machine 2>/dev/null)
      model=$(sysctl -n machdep.cpu.brand_string 2>/dev/null)
      cpus=$(sysctl -n hw.logicalcpu 2>/dev/null)
      cores=$(sysctl -n hw.physicalcpu 2>/dev/null)
      printf "  %-18s %s\n" "Architecture:"  "${arch:-?}"
      printf "  %-18s %s\n" "CPU model:"     "${model:-?}"
      printf "  %-18s %s  (%s core)\n" "CPU count:" "${cpus:-?}" "${cores:-?}"
    else
      echo "  (sysctl not available)"
    fi
  else
    if command -v lscpu &>/dev/null; then
      local lscpu_out
      lscpu_out=$(lscpu 2>/dev/null)
      local arch model cpus sockets cores_per_sock threads_per_core max_mhz min_mhz numa cache_l3
      arch=$(            awk -F': +' '/^Architecture/{print $2}'        <<< "$lscpu_out")
      model=$(           awk -F': +' '/^Model name/{print $2}'          <<< "$lscpu_out")
      cpus=$(            awk -F': +' '/^CPU\(s\):/{print $2; exit}'     <<< "$lscpu_out")
      sockets=$(         awk -F': +' '/^Socket\(s\)/{print $2}'         <<< "$lscpu_out")
      cores_per_sock=$(  awk -F': +' '/^Core\(s\) per socket/{print $2}' <<< "$lscpu_out")
      threads_per_core=$(awk -F': +' '/^Thread\(s\) per core/{print $2}' <<< "$lscpu_out")
      max_mhz=$(         awk -F': +' '/^CPU max MHz/{print $2}'         <<< "$lscpu_out")
      min_mhz=$(         awk -F': +' '/^CPU min MHz/{print $2}'         <<< "$lscpu_out")
      numa=$(            awk -F': +' '/^NUMA node\(s\)/{print $2}'      <<< "$lscpu_out")
      cache_l3=$(        awk -F': +' '/^L3 cache/{print $2}'            <<< "$lscpu_out")
      printf "  %-18s %s\n" "Architecture:"  "${arch:-?}"
      printf "  %-18s %s\n" "CPU model:"     "${model:-?}"
      printf "  %-18s %s  (%s socket × %s core × %s thread)\n" \
        "CPU count:" "${cpus:-?}" "${sockets:-?}" "${cores_per_sock:-?}" "${threads_per_core:-?}"
      if [[ -n "$max_mhz" ]]; then
        printf "  %-18s %s MHz max%s\n" "CPU freq:" "$max_mhz" \
          "${min_mhz:+, $min_mhz MHz min}"
      fi
      [[ -n "$numa"     ]] && printf "  %-18s %s\n" "NUMA nodes:"  "$numa"
      [[ -n "$cache_l3" ]] && printf "  %-18s %s\n" "L3 cache:"    "$cache_l3"
    else
      echo "  (lscpu not available)"
    fi
  fi
  echo ""

  local mem_total mem_used mem_avail
  if [[ "$(uname -s)" == "Darwin" ]]; then
    if command -v sysctl &>/dev/null; then
      local bytes
      bytes=$(sysctl -n hw.memsize 2>/dev/null || echo 0)
      mem_total=$(( bytes / 1048576 ))
      mem_used="?"
      mem_avail="?"
    else
      mem_total="?"
      mem_used="?"
      mem_avail="?"
    fi
  else
    read -r mem_total mem_used mem_avail \
      < <(free -m 2>/dev/null | awk 'NR==2{print $2, $3, $7}' || echo "? ? ?")
  fi
  printf "  %-18s %s MiB total, %s MiB available (%s MiB used)\n" \
    "Memory:" "${mem_total:-?}" "${mem_avail:-?}" "${mem_used:-?}"
  printf "  %-18s %s\n" "Kernel:"   "$(uname -r   2>/dev/null || echo '?')"
  printf "  %-18s %s\n" "Hostname:" "$(hostname    2>/dev/null || echo '?')"
  echo "═══════════════════════════════════════════════════════════════════════"
  echo ""
}

# ── Write experiment_setup.md to the output directory ────────────────────────
# Produces a self-contained Markdown document capturing everything needed to
# identify what was run, how CBC was built, and on what hardware.
write_setup_md() {
  local bin="$1"
  local md="$OUTDIR/experiment_setup.md"
  local src_dir="$SCRIPT_DIR"

  # ── Collect build data ────────────────────────────────────────────────────
  local cbc_cs
  cbc_cs="$(_project_config_status Cbc)"
  local cc="" cxx="" cc_ver="" cxx_ver="" add_cxxflags=""
  local bin_mtime="" bin_size_kb=0

  if [[ -f "$cbc_cs" ]]; then
    local cfg
    cfg=$(_configure_invocation_from_status "$cbc_cs")
    cc=$( printf '%s\n' "$cfg" | "$GREP_CMD" -oP 'CC=[^'"'"' ]+'  | head -1 | cut -d= -f2 || true)
    cxx=$(printf '%s\n' "$cfg" | "$GREP_CMD" -oP 'CXX=[^'"'"' ]+' | head -1 | cut -d= -f2 || true)
    [[ -n "$cc"  ]] && command -v "$cc"  &>/dev/null && cc_ver=$( "$cc"  --version 2>&1 | head -1)
    [[ -n "$cxx" ]] && command -v "$cxx" &>/dev/null && cxx_ver=$("$cxx" --version 2>&1 | head -1)
    add_cxxflags=$(_subst_from_status "$cbc_cs" "ADD_CXXFLAGS")
  fi
  bin_mtime=$(stat -c "%y" "$bin" 2>/dev/null | cut -d. -f1 || stat -f "%Sm" -t "%Y-%m-%d %H:%M:%S" "$bin" 2>/dev/null || echo "?")
  bin_size_kb=$(( $(stat -c "%s" "$bin" 2>/dev/null || echo 0) / 1024 ))

  # Linked libraries (lib name + resolved path, filtered to interesting ones)
  # Capture raw output once; reuse for pkg-config filtering below.
  local ldd_raw_md=""
  ldd_raw_md=$({ ldd "$bin" 2>/dev/null || otool -L "$bin" 2>/dev/null; } \
    | "$GREP_CMD" -E "(openblas|amd|lapack|blas|nauty|colamd|camd|ccolamd|cholmod|suitesparse|gomp)" \
    || true)
  local ldd_out=""
  if [[ -n "$ldd_raw_md" ]]; then
    ldd_out=$(awk '{lib=$1; path=($3!=""?$3:"(not found)"); printf "| `%s` | `%s` |\n", lib, path}' \
      <<< "$ldd_raw_md")
  else
    ldd_out="_(none matched / ldd not available)_"
  fi

  # pkg-config versions — only for libraries actually present in ldd output.
  local -A pkglibs=([openblas]="openblas" [amd]="libamd" [lapack]="liblapack")
  local pkgver_rows=""
  for lib in openblas amd lapack; do
    grep -q "${pkglibs[$lib]}" <<< "$ldd_raw_md" 2>/dev/null || continue
    local ver
    ver=$(pkg-config --modversion "$lib" 2>/dev/null || echo "(not found)")
    pkgver_rows+="| $lib | $ver |"$'\n'
  done

  # Hardware
  local arch="" cpu_model="" cpus="" sockets="" cores_per_sock="" threads_per_core=""
  local max_mhz="" min_mhz="" numa="" cache_l3=""
  local mem_total="" mem_used="" mem_avail=""
  if [[ "$(uname -s)" == "Darwin" ]]; then
    if command -v sysctl &>/dev/null; then
      arch=$(sysctl -n hw.machine 2>/dev/null)
      cpu_model=$(sysctl -n machdep.cpu.brand_string 2>/dev/null)
      cpus=$(sysctl -n hw.logicalcpu 2>/dev/null)
      cores_per_sock=$(sysctl -n hw.physicalcpu 2>/dev/null)
    fi
  else
    if command -v lscpu &>/dev/null; then
      local lscpu_out
      lscpu_out=$(lscpu 2>/dev/null)
      arch=$(            awk -F': +' '/^Architecture/{print $2}'          <<< "$lscpu_out")
      cpu_model=$(       awk -F': +' '/^Model name/{print $2}'            <<< "$lscpu_out")
      cpus=$(            awk -F': +' '/^CPU\(s\):/{print $2; exit}'       <<< "$lscpu_out")
      sockets=$(         awk -F': +' '/^Socket\(s\)/{print $2}'           <<< "$lscpu_out")
      cores_per_sock=$(  awk -F': +' '/^Core\(s\) per socket/{print $2}'  <<< "$lscpu_out")
      threads_per_core=$(awk -F': +' '/^Thread\(s\) per core/{print $2}'  <<< "$lscpu_out")
      max_mhz=$(         awk -F': +' '/^CPU max MHz/{print $2}'           <<< "$lscpu_out")
      min_mhz=$(         awk -F': +' '/^CPU min MHz/{print $2}'           <<< "$lscpu_out")
      numa=$(            awk -F': +' '/^NUMA node\(s\)/{print $2}'        <<< "$lscpu_out")
      cache_l3=$(        awk -F': +' '/^L3 cache/{print $2}'              <<< "$lscpu_out")
    fi
  fi
  
  if [[ "$(uname -s)" == "Darwin" ]]; then
    if command -v sysctl &>/dev/null; then
      local bytes
      bytes=$(sysctl -n hw.memsize 2>/dev/null || echo 0)
      mem_total=$(( bytes / 1048576 ))
      mem_used="?"
      mem_avail="?"
    else
      mem_total="?"
      mem_used="?"
      mem_avail="?"
    fi
  else
    read -r mem_total mem_used mem_avail \
      < <(free -m 2>/dev/null | awk 'NR==2{print $2, $3, $7}' || echo "? ? ?")
  fi

  # ── Write Markdown ────────────────────────────────────────────────────────
  {
    echo "# CBC Experiment Setup"
    echo ""
    echo "**Started:** $(date)"
    echo "**Outdir:** \`$OUTDIR\`"
    echo ""

    echo "## Experiment Parameters"
    echo ""
    echo "| Parameter | Value |"
    echo "|:---|:---|"
    echo "| Binary | \`$bin\` |"
    echo "| Threads per instance | $THREADS |"
    echo "| Concurrent instances | $PARALLEL |"
    echo "| Time limit | ${TIMELIMIT} s (+${OVERTIME_GRACE} s overtime grace) |"
    echo "| Instances dir | \`$INSTANCES_DIR\` ($TOTAL instances) |"
    if [[ ${#CBC_EXTRA_OPTS[@]} -gt 0 ]]; then
      echo "| Extra CBC opts | \`${CBC_EXTRA_OPTS[*]}\` |"
    else
      echo "| Extra CBC opts | _(none)_ |"
    fi
    if [[ -n "$VG_TOOL" ]]; then
      echo "| Valgrind | $VG_TOOL |"
    else
      echo "| Valgrind | _(none)_ |"
    fi
    if [[ -n "${ASAN_OPTIONS:-}" ]]; then
      echo "| ASAN\_OPTIONS | \`$ASAN_OPTIONS\` |"
    else
      echo "| ASAN\_OPTIONS | _(not set)_ |"
    fi
    echo "| ZH constraint report | $([ "$ZH_CONSTRAINT_REPORT" -eq 1 ] && echo yes || echo no) |"
    echo ""

    echo "## Build Environment"
    echo ""
    echo "**Source dir:** \`$src_dir\`"
    echo ""

    echo "### Configure Invocation"
    echo ""
    local cs
    cs="$(_project_config_status Cbc)"
    if [[ -f "$cs" ]]; then
      local cfg
      cfg=$(_configure_invocation_from_status "$cs")
      if [[ -n "$cfg" ]]; then
        echo '```'
        echo "$cfg"
        echo '```'
        echo ""
      fi
    fi

    echo "### Compilers"
    echo ""
    echo "| Role | Binary | Version |"
    echo "|:---|:---|:---|"
    echo "| \`CC\`  | \`${cc:-?}\`  | ${cc_ver:-?} |"
    echo "| \`CXX\` | \`${cxx:-?}\` | ${cxx_ver:-?} |"
    echo ""
    if [[ -n "$add_cxxflags" ]]; then
      echo "**Optimization flags (\`ADD_CXXFLAGS\`):**"
      echo ""
      echo '```'
      echo "$add_cxxflags"
      echo '```'
      echo ""
    fi

    echo "### Binary"
    echo ""
    echo "| | |"
    echo "|:---|:---|"
    echo "| Path  | \`$bin\` |"
    echo "| Built | $bin_mtime |"
    echo "| Size  | ${bin_size_kb} KB |"
    echo ""

    echo "### Linked Libraries"
    echo ""
    echo "| Library | Resolved path |"
    echo "|:---|:---|"
    echo "$ldd_out"
    echo ""

    if [[ -n "$pkgver_rows" ]]; then
      echo "### Dependency Versions"
      echo ""
      echo "| Library | Version |"
      echo "|:---|:---|"
      printf "%s" "$pkgver_rows"
      echo ""
    fi

    echo "## Hardware"
    echo ""
    echo "| Property | Value |"
    echo "|:---|:---|"
    echo "| Hostname | \`$(hostname 2>/dev/null || echo '?')\` |"
    echo "| Kernel | \`$(uname -r 2>/dev/null || echo '?')\` |"
    echo "| Architecture | ${arch:-?} |"
    echo "| CPU model | ${cpu_model:-?} |"
    echo "| CPU count | ${cpus:-?} (${sockets:-?} socket × ${cores_per_sock:-?} core × ${threads_per_core:-?} thread) |"
    if [[ -n "$max_mhz" ]]; then
      echo "| CPU freq | ${max_mhz} MHz max${min_mhz:+, ${min_mhz} MHz min} |"
    fi
    [[ -n "$numa"     ]] && echo "| NUMA nodes | $numa |"
    [[ -n "$cache_l3" ]] && echo "| L3 cache | $cache_l3 |"
    echo "| Memory | ${mem_total:-?} MiB total, ${mem_avail:-?} MiB available (${mem_used:-?} MiB used) |"
    echo ""
  } > "$md"

  echo "==> Experiment setup saved: $md"
}

{
  echo "CBC parallel experiment run"
  echo "Started:        $(date)"
  echo "Binary:         $ORIGINAL_CBC_BIN"
  echo "Binary (snap):  $CBC_BIN"
  echo "Threads:        $THREADS  (per instance)"
  echo "Parallel:       $PARALLEL  (concurrent instances)"
  echo "Timelimit:      ${TIMELIMIT}s  (+${OVERTIME_GRACE}s overtime grace before hard kill)"
  echo "Instances:      $TOTAL  (from $INSTANCES_DIR)"
  echo "Outdir:         $OUTDIR"
  [[ -n "${ASAN_OPTIONS:-}" ]] && echo "ASAN_OPTIONS:   $ASAN_OPTIONS"
  [[ -n "$VG_TOOL" ]] && echo "Valgrind:       $VG_TOOL"
  [[ "$ZH_CONSTRAINT_REPORT" -eq 1 ]] && echo "ZH reports:     enabled"
  echo ""
  print_env_info "$ORIGINAL_CBC_BIN"
} | tee "$REPORT"
write_setup_md "$ORIGINAL_CBC_BIN"

# ── Memory monitor: sample free memory every 5 minutes ───────────────────────
# Writes memory_usage.tsv with columns:
#   timestamp  elapsed_s  total_mb  used_mb  free_mb  available_mb  pct_used
MEM_TSV="$OUTDIR/memory_usage.tsv"
MEM_MONITOR_PID=""
EXP_START_EPOCH=$SECONDS
printf "timestamp\telapsed_s\ttotal_mb\tused_mb\tfree_mb\tavailable_mb\tpct_used\n" > "$MEM_TSV"
(
  trap 'exit 0' TERM INT
  while true; do
    _ts=$(date -Iseconds)
    _elapsed=$(( SECONDS - EXP_START_EPOCH ))
    read -r _total _used _free _avail \
      < <(free -m 2>/dev/null | awk 'NR==2{print $2, $3, $4, $7}') || true
    _pct=$(awk -v u="${_used:-0}" -v t="${_total:-1}" \
      'BEGIN { printf "%.1f", u/t*100 }')
    printf "%s\t%d\t%d\t%d\t%d\t%d\t%s\n" \
      "$_ts" "$_elapsed" "${_total:-0}" "${_used:-0}" \
      "${_free:-0}" "${_avail:-0}" "${_pct:-0}" >> "$MEM_TSV"
    sleep 300 &
    wait $! 2>/dev/null || true
  done
) &
MEM_MONITOR_PID=$!

# ── Run instances via GNU parallel ───────────────────────────────────────────
# setsid puts parallel in a new session so SIGHUP from a closing terminal
# does not reach it (parallel reinstalls its own SIGHUP handler and would
# otherwise stop accepting new jobs mid-run).
if [[ "$(uname -s)" == "Darwin" ]]; then
  if command -v setsid &>/dev/null; then
    SETSID_CMD=(setsid)
  else
    SETSID_CMD=()
  fi
else
  SETSID_CMD=(setsid)
fi

printf '%s\n' "${MPS_FILES[@]}" \
  | "${SETSID_CMD[@]}" parallel \
      --jobs "$PARALLEL" \
      --line-buffer \
      --joblog "$OUTDIR/parallel.log" \
      run_instance \
  || true  # non-zero exit when any worker fails; summary assembly must still run

# ── Merge per-instance stats CSVs into combined stats.csv ────────────────────
if [[ "$CBC_STATS_MODE" == "writeStat" ]]; then
  _stats_merged="$OUTDIR/stats.csv"
  # Build list of per-instance CSV files that exist
  _csv_files=()
  for _mps in "${MPS_FILES[@]}"; do
    _sf="$OUTDIR/$(basename "$_mps" .mps.gz).stats.csv"
    [[ -f "$_sf" ]] && _csv_files+=("$_sf")
  done
  if [[ ${#_csv_files[@]} -gt 0 ]]; then
    # Use Python to merge with proper column unification (handles varying generator columns)
    python3 - "${_csv_files[@]}" "$_stats_merged" <<'PYEOF'
import sys, csv, collections

inputs = sys.argv[1:-1]
output = sys.argv[-1]

# Read all rows, collect all column names in order
all_cols = collections.OrderedDict()  # col -> index (preserves insertion order)
rows = []
for path in inputs:
    with open(path, newline='') as f:
        reader = csv.DictReader(f)
        for col in (reader.fieldnames or []):
            all_cols.setdefault(col, len(all_cols))
        data = list(reader)
        if data:
            rows.append(data[-1])  # take last row (in case of duplicates)

headers = list(all_cols.keys())
with open(output, 'w', newline='') as f:
    writer = csv.DictWriter(f, fieldnames=headers, extrasaction='ignore')
    writer.writeheader()
    for row in rows:
        writer.writerow({col: row.get(col, '0') for col in headers})
PYEOF
    echo "==> Merged per-instance stats → $_stats_merged"
  fi
fi

# ── Stop memory monitor and take one final sample ────────────────────────────
if [[ -n "${MEM_MONITOR_PID:-}" ]]; then
  kill "$MEM_MONITOR_PID" 2>/dev/null || true
  wait "$MEM_MONITOR_PID" 2>/dev/null || true
fi
# Final sample: captures the post-parallel memory state
_ts=$(date -Iseconds)
_elapsed=$(( SECONDS - EXP_START_EPOCH ))
read -r _total _used _free _avail \
  < <(free -m 2>/dev/null | awk 'NR==2{print $2, $3, $4, $7}') || true
_pct=$(awk -v u="${_used:-0}" -v t="${_total:-1}" \
  'BEGIN { printf "%.1f", u/t*100 }')
printf "%s\t%d\t%d\t%d\t%d\t%d\t%s\n" \
  "$_ts" "$_elapsed" "${_total:-0}" "${_used:-0}" \
  "${_free:-0}" "${_avail:-0}" "${_pct:-0}" >> "$MEM_TSV"

# ── Assemble summary from per-instance .result files ─────────────────────────
printf "%-42s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
  "instance" "status" "objective" "dual_bound" "expected" "elapsed_s" "threads" "gap_field" \
  "solution_found" "proven_infeasible" "timed_out" > "$SUMMARY"

SOLVED_COUNT=0; TIMEOUT_COUNT=0; NO_SOLUTION_COUNT=0; FAIL=0; OVERTIME_COUNT=0
FAILED_NAMES=()

# shellcheck source=/dev/null
source "$EXPECTED_FILE"

for mps in "${MPS_FILES[@]}"; do
  name=$(basename "$mps" .mps.gz)
  resultfile="$OUTDIR/${name}.result"
  if [[ -f "$resultfile" ]]; then
    IFS=$'\t' read -r status obj dual exp elapsed thr gap solution_found proven_infeasible timed_out < "$resultfile"
  else
    status="MISSING"; obj="-"; dual="-"; exp="${EXPECTED[$name]:--}"; elapsed=0; thr=$THREADS; gap="-"
    solution_found=0; proven_infeasible=0; timed_out=0
  fi
  printf "%-42s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
    "$name" "$status" "$obj" "${dual:--}" "$exp" "$elapsed" "$thr" "${gap:--}" \
    "$solution_found" "$proven_infeasible" "$timed_out" >> "$SUMMARY"
  if [[ "$status" == SOLVED || "$status" == SOLVED\(* ]]; then
    SOLVED_COUNT=$((SOLVED_COUNT+1))
  elif [[ "$status" == TIMEOUT* ]]; then
    TIMEOUT_COUNT=$((TIMEOUT_COUNT+1))
  elif [[ "$status" == "NO_SOLUTION" ]]; then
    NO_SOLUTION_COUNT=$((NO_SOLUTION_COUNT+1))
  else
    FAIL=$((FAIL+1))
    FAILED_NAMES+=("$name")
    [[ "$status" == "OVERTIME" ]] && OVERTIME_COUNT=$((OVERTIME_COUNT+1))
  fi
done

# ── Compute peak memory utilisation from the TSV ─────────────────────────────
MEM_PEAK_USED=""; MEM_PEAK_TOTAL=""; MEM_PEAK_PCT=""; MEM_PEAK_TS=""
if [[ -f "$MEM_TSV" ]] && [[ $(wc -l < "$MEM_TSV") -gt 1 ]]; then
  read -r MEM_PEAK_USED MEM_PEAK_TOTAL MEM_PEAK_PCT MEM_PEAK_TS < <(
    awk -F'\t' 'NR>1 && $4+0>max { max=$4+0; total=$3+0; pct=$7; ts=$1 }
                END { if (max>0) print max, total, pct, ts }' "$MEM_TSV"
  ) || true
fi

# ── Final report ──────────────────────────────────────────────────────────────
{
  echo ""
  echo "═══════════════════════════════════════════════════"
  printf "  Results: %d/%d solved/proved, %d timed out, %d no-solution, %d failed" \
    "$SOLVED_COUNT" "$TOTAL" "$TIMEOUT_COUNT" "$NO_SOLUTION_COUNT" "$FAIL"
  [[ $OVERTIME_COUNT -gt 0 ]] && printf "  (%d OVERTIME — ignored time limit)" "$OVERTIME_COUNT"
  echo ""
  echo "  Finished: $(date)"
  if [[ -n "$MEM_PEAK_USED" ]]; then
    echo "  Peak memory used:  ${MEM_PEAK_USED} MB / ${MEM_PEAK_TOTAL} MB (${MEM_PEAK_PCT}%)"
    echo "  Peak at:           $MEM_PEAK_TS"
  fi
  echo "═══════════════════════════════════════════════════"

  if [[ ${#FAILED_NAMES[@]} -gt 0 ]]; then
    echo ""
    echo "Failed instances:"
    for n in "${FAILED_NAMES[@]}"; do
      st=$(awk -F'\t' -v n="$n" '$1==n{print $2}' "$SUMMARY")
      printf "  %-42s  %s\n" "$n" "$st"
      [[ -f "$OUTDIR/${n}.err" ]] && printf "    error detail: %s\n" "$OUTDIR/${n}.err"
    done
  fi

  echo ""
  echo "Per-instance logs: $OUTDIR/"
  echo "Summary TSV:       $SUMMARY"
  echo "Solutions:         $OUTDIR/*.sol  (present when feasible solution found)"
  echo "GNU parallel log:  $OUTDIR/parallel.log"
  [[ -f "$MEM_TSV" ]] && echo "Memory usage log:  $MEM_TSV"
  echo "Experiment setup:  $OUTDIR/experiment_setup.md"
} | tee -a "$REPORT"

exit 0
