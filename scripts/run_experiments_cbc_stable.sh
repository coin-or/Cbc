#!/usr/bin/env bash
# run_experiments_cbc_stable.sh
#
# Run CBC stable (COIN-OR Cbc 2.x) on MIP instances and produce a summary file
# compatible with summarize_results.py and compare_experiments.py.
#
# Key differences from run_experiments.sh (MIPster):
#   - Always injects -timeMode elapsed before -seconds (CBC stable defaults to
#     CPU time; we need wall-clock time to be comparable with MIPster results)
#   - No -writeStat support (CBC stable does not have this MIPster extension);
#     all result extraction is done purely from log parsing
#   - No mipster_validate_sol step (not applicable to CBC stable)
#   - No build / ASan / TSan / valgrind options
#   - Log parsing uses precise patterns tied to CBC stable's Result block format:
#       "^Objective value:" for the best incumbent
#       "^Lower bound:"     for the dual bound
#       "^Gap:"             for the primal-dual gap
#       "^Result - Stopped on time" to detect a time-limit stop
#       "^Result - .*infeasible"    to detect proven infeasibility
#
# Output layout ($OUTDIR/):
#   summary.tsv              — one row per instance: status, obj, elapsed, gap
#   <instance>.log           — full solver stdout/stderr
#   <instance>.sol           — solution file (present when feasible solution found)
#   <instance>.err           — error snippet (wrong obj / crash)
#   <instance>.result        — machine-readable single-row result (for summary assembly)
#   parallel.log             — GNU parallel job log (timing, exit codes)
#   report.txt               — human-readable summary
#   memory_usage.tsv         — memory utilisation sampled every 5 min during the run
#   experiment_setup.md      — binary info, hardware, and experiment parameters
#
# Status values in summary.tsv:
#   SOLVED            solved, objective matches exact reference
#   SOLVED(best)      solved, objective matches best-known reference
#   SOLVED(improved)  solved, objective improves on best-known reference
#   SOLVED(no_ref)    solved, no reference available
#   SOLVED(inf)       CBC proved infeasible; benchmark agrees (or no ref)
#   WRONG_OBJ         solved but objective worse than reference
#   INFEASIBLE_WRONG  CBC proved infeasible but a numeric objective was expected
#   NO_SOLUTION       CBC exited normally but found no feasible solution
#   TIMEOUT           timed out
#   TIMEOUT(no_sol)   timed out without finding any feasible solution
#   TIMEOUT(gap=X%)   timed out with feasible solution, gap X%
#   OVERTIME          CBC ignored -seconds limit and was hard-killed
#   CRASH             non-zero exit code, not a known timeout
#
# Usage:
#   ./run_experiments_cbc_stable.sh --bin ~/prog/cbc-stable/bin/cbc [OPTIONS]
#
# Required:
#   --bin PATH        Path to the CBC stable binary
#
# Optional:
#   --parallel N      Instances to run concurrently (default: 10)
#   --timelimit T     Wall-clock time limit in seconds (default: 3600)
#   --overtime-grace G  Seconds to wait beyond timelimit before hard-killing
#                     (default: 900).  Instance is tagged OVERTIME if killed.
#   --outdir DIR      Output directory (default: auto-named under exp_results/)
#   --instances DIR   Directory with .mps.gz files
#                     (default: /home/haroldo/inst/miplib/2017)
#   --opts FLAG       Extra CBC flag inserted before -solve (repeatable)
#                     e.g. --opts "-probing 1" --opts "-gomory 1"
#   --dry-run         Print commands without executing them
#   -h, --help        Show this help and exit
#
# Examples:
#   ./run_experiments_cbc_stable.sh \
#     --bin ~/prog/cbc-stable/bin/cbc \
#     --instances ~/inst/miplib/2017+spp \
#     --timelimit 300 --parallel 20
#
#   # Compare against MIPster on the same instance set:
#   ./run_experiments.sh     --bin ~/prog/cbc/bin/mipster --timelimit 300 --parallel 20
#   ./run_experiments_cbc_stable.sh --bin ~/prog/cbc-stable/bin/cbc \
#                                   --timelimit 300 --parallel 20
#   python3 scripts/compare_experiments.py exp_results/cbc_... exp_results/cbc-stable_...

set -euo pipefail
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

# ── Hardcoded reference values (same set as run_experiments.sh) ───────────────
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
  [dt_optimization]=202.185252817
)
REL_TOL=1e-6
ABS_TOL=1e-4

# ── Help ──────────────────────────────────────────────────────────────────────
show_help() {
  cat << 'EOF'
Usage: ./run_experiments_cbc_stable.sh --bin PATH [OPTIONS]

Run CBC stable (COIN-OR Cbc 2.x) on MIP instances and produce a summary
compatible with summarize_results.py and compare_experiments.py.

Required:
  --bin PATH        Path to the CBC stable binary

Optional:
  --parallel N      Instances to run concurrently (default: 10)
  --timelimit T     Wall-clock time limit in seconds (default: 3600)
  --overtime-grace G  Seconds before hard-killing after timelimit (default: 900)
  --outdir DIR      Output directory (default: auto-named under exp_results/)
  --instances DIR   Directory with .mps.gz files
                    (default: /home/haroldo/inst/miplib/2017)
  --opts FLAG       Extra CBC flag before -solve (repeatable)
  --dry-run         Print commands without executing them
  -h, --help        Show this help and exit

Notes:
  -timeMode elapsed is always injected before -seconds so the limit is measured
  in wall-clock time (CBC stable's default is CPU time).

Examples:
  ./run_experiments_cbc_stable.sh \
    --bin ~/prog/cbc-stable/bin/cbc \
    --instances ~/inst/miplib/2017+spp \
    --timelimit 300 --parallel 20
EOF
}

# ── Defaults ──────────────────────────────────────────────────────────────────
CBC_BIN=""
PARALLEL=10
TIMELIMIT=""
OVERTIME_GRACE=900
INSTANCES_DIR="/home/haroldo/inst/miplib/2017"
OUTDIR=""
CBC_EXTRA_OPTS=()
DRY_RUN=0

# ── Parse arguments ───────────────────────────────────────────────────────────
while [[ $# -gt 0 ]]; do
  case "$1" in
    --bin)            CBC_BIN="$2";             shift 2 ;;
    --parallel)       PARALLEL="$2";            shift 2 ;;
    --timelimit)      TIMELIMIT="$2";           shift 2 ;;
    --overtime-grace) OVERTIME_GRACE="$2";      shift 2 ;;
    --outdir)         OUTDIR="$2";              shift 2 ;;
    --instances)      INSTANCES_DIR="$2";       shift 2 ;;
    --opts)           CBC_EXTRA_OPTS+=("$2");   shift 2 ;;
    --dry-run)        DRY_RUN=1;                shift   ;;
    -h|--help)        show_help; exit 0                 ;;
    *) echo "Unknown option: $1" >&2; exit 1  ;;
  esac
done

[[ -z "$CBC_BIN" ]] && { echo "Error: --bin is required" >&2; exit 1; }
[[ -x "$CBC_BIN" ]] || { echo "Error: $CBC_BIN is not executable" >&2; exit 1; }

# ── Snapshot binary so rebuilds during the run don't interfere ───────────────
EXP_TMPDIR=$(mktemp -d /tmp/cbc_stable_exp_XXXXXXXX)
trap 'rm -rf "${EXP_TMPDIR:-}"' EXIT
mkdir -p "$EXP_TMPDIR/bin"
ORIGINAL_CBC_BIN="$CBC_BIN"
ORIGINAL_BINDIR="$(cd "$(dirname "$CBC_BIN")" && pwd)"
TMP_CBC_BIN="$EXP_TMPDIR/bin/$(basename "$CBC_BIN")"
cp "$CBC_BIN" "$TMP_CBC_BIN"
chmod +x "$TMP_CBC_BIN"
for f in "$ORIGINAL_BINDIR"/*; do
  [[ -f "$f" && -x "$f" ]] && cp "$f" "$EXP_TMPDIR/bin/"
done
if [[ -d "$ORIGINAL_BINDIR/../lib" ]]; then
  mkdir -p "$EXP_TMPDIR/lib"
  cp -rP "$ORIGINAL_BINDIR/../lib/"* "$EXP_TMPDIR/lib/" 2>/dev/null || true
fi
CBC_BIN="$TMP_CBC_BIN"
echo "==> Binary snapshotted: $TMP_CBC_BIN"
echo "    (original: $ORIGINAL_CBC_BIN)"
echo ""

# ── Time limit ────────────────────────────────────────────────────────────────
[[ -z "$TIMELIMIT" ]] && TIMELIMIT=3600
WALLCLOCK_LIMIT=$((TIMELIMIT + OVERTIME_GRACE))

# ── Output directory ──────────────────────────────────────────────────────────
if [[ -z "$OUTDIR" ]]; then
  VARIANT=$(basename "$(dirname "$(dirname "$ORIGINAL_CBC_BIN")")")
  TS=$(date +%Y%m%d_%H%M%S)
  OUTDIR="exp_results/${VARIANT}_p${PARALLEL}_${TS}"
fi
mkdir -p "$OUTDIR"

SUMMARY="$OUTDIR/summary.tsv"
REPORT="$OUTDIR/report.txt"

# ── Build expected.sh from instance metadata + hardcoded overrides ────────────
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

for name, value in bks_lookup.items():
    expected.setdefault(name, value)
    ref_kind.setdefault(name, "best")

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

# ── Serialise CBC_EXTRA_OPTS into newline-separated string for export ─────────
CBC_EXTRA_OPTS_STR=$(printf '%s\n' "${CBC_EXTRA_OPTS[@]+"${CBC_EXTRA_OPTS[@]}"}")

# ── Colour support ────────────────────────────────────────────────────────────
if [[ -t 1 && -z "${NO_COLOR:-}" ]]; then
  USE_COLOR=1
else
  USE_COLOR=0
fi

# ── Export everything run_instance needs ──────────────────────────────────────
export CBC_BIN OUTDIR TIMELIMIT OVERTIME_GRACE WALLCLOCK_LIMIT TIMEOUT_CMD GREP_CMD
if [[ "$(uname -s)" == "Darwin" ]]; then
  export PARALLEL_SHELL="$BASH"
fi
export CBC_EXTRA_OPTS_STR EXPECTED_FILE REL_TOL ABS_TOL DRY_RUN USE_COLOR

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

# ── Helper: format elapsed seconds as compact human-readable string ───────────
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
  local resultfile="$OUTDIR/${name}.result"
  local solutionfile="$OUTDIR/${name}.sol"

  # Reconstruct CBC_EXTRA_OPTS — each entry may contain multiple words.
  local cbc_extra=()
  if [[ -n "$CBC_EXTRA_OPTS_STR" ]]; then
    while IFS= read -r opt; do
      [[ -n "$opt" ]] || continue
      read -ra words <<< "$opt"
      cbc_extra+=("${words[@]}")
    done <<< "$CBC_EXTRA_OPTS_STR"
  fi

  # -timeMode elapsed MUST come before -seconds so the limit is in wall-clock
  # time (CBC stable defaults to CPU time).
  local cmd=("$CBC_BIN" "$mps"
             "-timeMode" "elapsed"
             "-seconds" "$TIMELIMIT")
  [[ ${#cbc_extra[@]} -gt 0 ]] && cmd+=("${cbc_extra[@]}")
  cmd+=("-solve" "-solu" "$solutionfile" "-quit")

  if [[ "$DRY_RUN" -eq 1 ]]; then
    echo "DRY-RUN: ${cmd[*]}"
    return 0
  fi

  local start_epoch=$SECONDS
  set +e
  "$TIMEOUT_CMD" --kill-after=30 "$WALLCLOCK_LIMIT" \
    "${cmd[@]}" > "$logfile" 2>&1
  local exit_code=$?
  set -e
  local elapsed=$((SECONDS - start_epoch))

  # ── Extract results from log ──────────────────────────────────────────────
  # Use the final "Result - ..." line as the canonical status indicator.
  # CBC stable prints this block once per solve (occasionally twice for presolve
  # + main solve); tail -1 always picks the definitive last one.
  local result_line
  result_line=$(grep '^Result - ' "$logfile" | tail -1 || true)

  local proven_infeasible=0
  echo "$result_line" | grep -qi 'infeasible' && proven_infeasible=1

  local timed_out=0
  echo "$result_line" | grep -qi 'stopped on time' && timed_out=1

  # Best incumbent — "Objective value:" appears only in the Result block and
  # only when CBC has found a feasible integer solution.
  local obj
  obj=$(grep '^Objective value:' "$logfile" | tail -1 \
        | "$GREP_CMD" -oP '[-\d.eE+]+$' || true)

  # Discard sentinel "no incumbent" values CBC sometimes emits (±1e50).
  if [[ -n "$obj" ]] && awk -v x="$obj" 'BEGIN { exit (x >= 1e49 || x <= -1e49) ? 0 : 1 }'; then
    obj=""
  fi

  local no_feasible_solution=0
  [[ -z "$obj" && $proven_infeasible -eq 0 ]] && no_feasible_solution=1

  # Dual bound (lower bound for minimisation).
  local dual
  dual=$(grep '^Lower bound:' "$logfile" | tail -1 \
         | "$GREP_CMD" -oP '[-\d.eE+]+$' || true)
  if [[ -n "$dual" ]] && awk -v x="$dual" 'BEGIN { exit (x >= 1e49 || x <= -1e49) ? 0 : 1 }'; then
    dual=""
  fi
  # For proven-optimal cases the dual bound equals the primal; fill from obj.
  if [[ -z "$dual" && -n "$obj" && $timed_out -eq 0 && $proven_infeasible -eq 0 ]]; then
    dual="$obj"
  fi

  local solution_found=0
  [[ -n "$obj" ]] && solution_found=1

  # ── Load reference values ─────────────────────────────────────────────────
  # shellcheck source=/dev/null
  source "$EXPECTED_FILE"
  local exp="${EXPECTED[$name]:-}"
  local exp_infeasible="${INFEASIBLE[$name]:-}"
  local ref_kind="${REF_KIND[$name]:-exact}"
  local ref_sense="${REF_SENSE[$name]:-}"

  # ── Gap for timed-out runs with a feasible solution ───────────────────────
  local cbc_gap="-"
  if [[ $timed_out -eq 1 && -n "$obj" ]]; then
    local raw_gap
    raw_gap=$(grep '^Gap:' "$logfile" | tail -1 \
              | "$GREP_CMD" -oP '[-\d.eE+]+$' || true)
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

  # ── Classify status ───────────────────────────────────────────────────────
  local status
  if [[ $exit_code -eq 124 ]]; then
    status="OVERTIME"
  elif [[ $exit_code -ne 0 ]]; then
    status="CRASH(exit=$exit_code)"
  elif [[ $proven_infeasible -eq 1 ]]; then
    if [[ -n "$exp_infeasible" || -z "$exp" ]]; then
      status="SOLVED(inf)"
    else
      status="INFEASIBLE_WRONG"
    fi
  elif [[ -z "$obj" ]]; then
    if [[ $timed_out -eq 1 ]]; then
      status="TIMEOUT(no_sol)"
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

  # ── Write .err snippet for actionable failures ────────────────────────────
  > "$errfile"
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
  if [[ "$status" =~ ^CRASH ]]; then
    { echo "=== Exit code: $exit_code ==="; echo ""; } >> "$errfile"
  fi
  [[ ! -s "$errfile" ]] && rm -f "$errfile"

  # ── Write machine-readable result for summary assembly ────────────────────
  # Format is identical to run_experiments.sh so compare/summarize scripts work.
  printf "%s\t%s\t%s\t%s\t%d\t%s\t%s\t%d\t%d\t%d\n" \
    "$status" "${obj:--}" "${dual:--}" "${exp:--}" "$elapsed" "1" "$cbc_gap" \
    "$solution_found" "$proven_infeasible" "$timed_out" > "$resultfile"

  # ── One-line console output ───────────────────────────────────────────────
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

  local display_name="$name"
  (( ${#name} > 38 )) && display_name="${name:0:37}…"

  local extra=""
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

# ── Write experiment_setup.md ─────────────────────────────────────────────────
write_setup_md() {
  local bin="$1"
  local md="$OUTDIR/experiment_setup.md"

  local bin_version bin_build_date bin_mtime bin_size_kb
  bin_version=$(  "$bin" 2>&1 </dev/null | grep '^Version:'    | head -1 | awk '{print $2}' || true)
  bin_build_date=$("$bin" 2>&1 </dev/null | grep '^Build Date:' | head -1 | cut -d: -f2- | xargs || true)
  bin_mtime=$(stat -c "%y" "$bin" 2>/dev/null | cut -d. -f1 \
              || stat -f "%Sm" -t "%Y-%m-%d %H:%M:%S" "$bin" 2>/dev/null || echo "?")
  bin_size_kb=$(( $(stat -c "%s" "$bin" 2>/dev/null || echo 0) / 1024 ))

  local arch="" cpu_model="" cpus="" sockets="" cores_per_sock="" threads_per_core=""
  local max_mhz="" numa="" cache_l3="" mem_total="" mem_used="" mem_avail=""
  if [[ "$(uname -s)" == "Darwin" ]]; then
    if command -v sysctl &>/dev/null; then
      arch=$(sysctl -n hw.machine 2>/dev/null)
      cpu_model=$(sysctl -n machdep.cpu.brand_string 2>/dev/null)
      cpus=$(sysctl -n hw.logicalcpu 2>/dev/null)
      cores_per_sock=$(sysctl -n hw.physicalcpu 2>/dev/null)
      local bytes; bytes=$(sysctl -n hw.memsize 2>/dev/null || echo 0)
      mem_total=$(( bytes / 1048576 )); mem_used="?"; mem_avail="?"
    fi
  else
    if command -v lscpu &>/dev/null; then
      local lscpu_out; lscpu_out=$(lscpu 2>/dev/null)
      arch=$(            awk -F': +' '/^Architecture/{print $2}'          <<< "$lscpu_out")
      cpu_model=$(       awk -F': +' '/^Model name/{print $2}'            <<< "$lscpu_out")
      cpus=$(            awk -F': +' '/^CPU\(s\):/{print $2; exit}'       <<< "$lscpu_out")
      sockets=$(         awk -F': +' '/^Socket\(s\)/{print $2}'           <<< "$lscpu_out")
      cores_per_sock=$(  awk -F': +' '/^Core\(s\) per socket/{print $2}'  <<< "$lscpu_out")
      threads_per_core=$(awk -F': +' '/^Thread\(s\) per core/{print $2}'  <<< "$lscpu_out")
      max_mhz=$(         awk -F': +' '/^CPU max MHz/{print $2}'           <<< "$lscpu_out")
      numa=$(            awk -F': +' '/^NUMA node\(s\)/{print $2}'        <<< "$lscpu_out")
      cache_l3=$(        awk -F': +' '/^L3 cache/{print $2}'              <<< "$lscpu_out")
    fi
    read -r mem_total mem_used mem_avail \
      < <(free -m 2>/dev/null | awk 'NR==2{print $2, $3, $7}' || echo "? ? ?")
  fi

  {
    echo "# CBC Stable Experiment Setup"
    echo ""
    echo "**Started:** $(date)"
    echo "**Outdir:** \`$OUTDIR\`"
    echo ""
    echo "## Experiment Parameters"
    echo ""
    echo "| Parameter | Value |"
    echo "|:---|:---|"
    echo "| Binary | \`$bin\` |"
    echo "| CBC version | ${bin_version:-?} |"
    echo "| CBC build date | ${bin_build_date:-?} |"
    echo "| Concurrent instances | $PARALLEL |"
    echo "| Time limit | ${TIMELIMIT} s wall-clock (+${OVERTIME_GRACE} s overtime grace) |"
    echo "| Time mode | elapsed (wall-clock; -timeMode elapsed injected automatically) |"
    echo "| Instances dir | \`$INSTANCES_DIR\` ($TOTAL instances) |"
    if [[ ${#CBC_EXTRA_OPTS[@]} -gt 0 ]]; then
      echo "| Extra CBC opts | \`${CBC_EXTRA_OPTS[*]}\` |"
    else
      echo "| Extra CBC opts | _(none)_ |"
    fi
    echo ""
    echo "## Binary"
    echo ""
    echo "| | |"
    echo "|:---|:---|"
    echo "| Path  | \`$bin\` |"
    echo "| Version | ${bin_version:-?} |"
    echo "| Build date | ${bin_build_date:-?} |"
    echo "| File mtime | $bin_mtime |"
    echo "| Size  | ${bin_size_kb} KB |"
    echo ""
    echo "## Hardware"
    echo ""
    echo "| Property | Value |"
    echo "|:---|:---|"
    echo "| Hostname | \`$(hostname 2>/dev/null || echo '?')\` |"
    echo "| Kernel | \`$(uname -r 2>/dev/null || echo '?')\` |"
    echo "| Architecture | ${arch:-?} |"
    echo "| CPU model | ${cpu_model:-?} |"
    echo "| CPU count | ${cpus:-?} (${sockets:-?} socket × ${cores_per_sock:-?} core × ${threads_per_core:-?} thread) |"
    [[ -n "$max_mhz"  ]] && echo "| CPU freq | ${max_mhz} MHz max |"
    [[ -n "$numa"     ]] && echo "| NUMA nodes | $numa |"
    [[ -n "$cache_l3" ]] && echo "| L3 cache | $cache_l3 |"
    echo "| Memory | ${mem_total:-?} MiB total, ${mem_avail:-?} MiB available (${mem_used:-?} MiB used) |"
    echo ""
  } > "$md"

  echo "==> Experiment setup saved: $md"
}

# ── Print header and setup document ──────────────────────────────────────────
{
  echo "CBC stable parallel experiment run"
  echo "Started:        $(date)"
  echo "Binary:         $ORIGINAL_CBC_BIN"
  echo "Binary (snap):  $CBC_BIN"
  echo "Parallel:       $PARALLEL  (concurrent instances)"
  echo "Timelimit:      ${TIMELIMIT}s wall-clock  (+${OVERTIME_GRACE}s overtime grace)"
  echo "Instances:      $TOTAL  (from $INSTANCES_DIR)"
  echo "Outdir:         $OUTDIR"
  echo ""
} | tee "$REPORT"
write_setup_md "$ORIGINAL_CBC_BIN"

# ── Memory monitor: sample free memory every 5 minutes ───────────────────────
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

# ── Run instances via GNU parallel ────────────────────────────────────────────
if [[ "$(uname -s)" == "Darwin" ]]; then
  command -v setsid &>/dev/null && SETSID_CMD=(setsid) || SETSID_CMD=()
else
  SETSID_CMD=(setsid)
fi

printf '%s\n' "${MPS_FILES[@]}" \
  | "${SETSID_CMD[@]}" parallel \
      --jobs "$PARALLEL" \
      --line-buffer \
      --joblog "$OUTDIR/parallel.log" \
      run_instance \
  || true

# ── Stop memory monitor and take one final sample ────────────────────────────
if [[ -n "${MEM_MONITOR_PID:-}" ]]; then
  kill "$MEM_MONITOR_PID" 2>/dev/null || true
  wait "$MEM_MONITOR_PID" 2>/dev/null || true
fi
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
# Column layout is identical to run_experiments.sh so compare/summarize scripts
# work unchanged on output from both scripts.
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
    IFS=$'\t' read -r status obj dual exp elapsed thr gap solution_found proven_infeasible timed_out \
      < "$resultfile"
  else
    status="MISSING"; obj="-"; dual="-"; exp="${EXPECTED[$name]:--}"; elapsed=0; thr=1; gap="-"
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

# ── Compute peak memory utilisation ──────────────────────────────────────────
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
