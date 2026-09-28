#!/usr/bin/env bash
# run_lp_experiments.sh — LP relaxation parameter tuning experiments
#
# Solves the LP relaxation (-initialSolve) of MILP instances with different
# parameter settings and random seeds.  Produces per-run .sol, .bas, .txt,
# .log files and a consolidated CSV summary.
#
# Features:
#   - Binary snapshot (rebuilds don't affect running experiment)
#   - GNU parallel for concurrent execution
#   - Resumability: skips jobs whose .result file already exists
#   - Hard kill timeout (tolerance beyond LP time limit)
#   - Build/hardware info saved to experiment_setup.md
#
# Usage:
#   ./run_lp_experiments.sh [OPTIONS]
#
# Required:
#   --bin PATH           Path to cbc binary (default: ~/prog/cbc/bin/cbc)
#
# Optional:
#   --params FILE        Parameter settings file (default: ./lp_params.txt)
#   --instances DIR      Directory with .mps.gz files
#                        (default: ~/inst/miplib/2017+spp)
#   --timelimit T        Time limit in seconds via -sec (default: 14400 = 4h)
#   --overtime G         Extra seconds before hard kill (default: 600 = 10min)
#   --seeds S1,S2,...    Comma-separated seeds (default: 123,1234)
#   --parallel N         Concurrent jobs (default: nproc - 2, min 1)
#   --outdir DIR         Experiment directory (default: auto-named)
#   --dry-run            Print job list without executing
#   -h, --help           Show this help

set -euo pipefail
trap '' HUP  # survive terminal disconnection

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# ── Defaults ──────────────────────────────────────────────────────────────────
CBC_BIN="${HOME}/prog/cbc/bin/cbc"
PARAMS_FILE="${SCRIPT_DIR}/lp_params.txt"
INSTANCES_DIR="${HOME}/inst/miplib/2017+spp"
LP_TIMELIMIT=14400    # 4 hours
OVERTIME=600          # 10 minutes grace
SEEDS="123,1234"
PARALLEL=""
OUTDIR=""
DRY_RUN=0

# ── Parse arguments ───────────────────────────────────────────────────────────
show_help() {
  sed -n '2,/^$/{ s/^# \?//; p }' "$0"
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --bin)        CBC_BIN="$2";        shift 2 ;;
    --params)     PARAMS_FILE="$2";    shift 2 ;;
    --instances)  INSTANCES_DIR="$2";  shift 2 ;;
    --timelimit)  LP_TIMELIMIT="$2";   shift 2 ;;
    --overtime)   OVERTIME="$2";       shift 2 ;;
    --seeds)      SEEDS="$2";          shift 2 ;;
    --parallel)   PARALLEL="$2";       shift 2 ;;
    --outdir)     OUTDIR="$2";         shift 2 ;;
    --dry-run)    DRY_RUN=1;           shift   ;;
    -h|--help)    show_help; exit 0            ;;
    *) echo "Unknown option: $1" >&2; exit 1   ;;
  esac
done

# ── Validate inputs ──────────────────────────────────────────────────────────
[[ -x "$CBC_BIN" ]] || { echo "Error: $CBC_BIN not found or not executable" >&2; exit 1; }
[[ -f "$PARAMS_FILE" ]] || { echo "Error: params file not found: $PARAMS_FILE" >&2; exit 1; }
[[ -d "$INSTANCES_DIR" ]] || { echo "Error: instances dir not found: $INSTANCES_DIR" >&2; exit 1; }

# Parse seeds
IFS=',' read -ra SEED_LIST <<< "$SEEDS"

# Default parallel: nproc - 2 (leave headroom), min 1
if [[ -z "$PARALLEL" ]]; then
  NCPU=$(nproc 2>/dev/null || echo 4)
  PARALLEL=$(( NCPU > 3 ? NCPU - 2 : 1 ))
fi

KILL_AFTER=$(( LP_TIMELIMIT + OVERTIME ))

# ── Read parameter settings ──────────────────────────────────────────────────
declare -a PARAM_TAGS=()
declare -a PARAM_ARGS=()
while IFS= read -r line; do
  # Skip comments and blank lines
  [[ "$line" =~ ^[[:space:]]*# ]] && continue
  [[ -z "${line// /}" ]] && continue
  TAG="${line%%|*}"
  ARGS="${line#*|}"
  PARAM_TAGS+=("$TAG")
  PARAM_ARGS+=("$ARGS")
done < "$PARAMS_FILE"

[[ ${#PARAM_TAGS[@]} -eq 0 ]] && { echo "Error: no parameter settings in $PARAMS_FILE" >&2; exit 1; }

# ── Find instances ────────────────────────────────────────────────────────────
mapfile -t INSTANCES < <(find "$INSTANCES_DIR" -maxdepth 1 -name "*.mps.gz" | sort)
[[ ${#INSTANCES[@]} -eq 0 ]] && { echo "Error: no .mps.gz files in $INSTANCES_DIR" >&2; exit 1; }

N_INST=${#INSTANCES[@]}
N_PARAMS=${#PARAM_TAGS[@]}
N_SEEDS=${#SEED_LIST[@]}
N_JOBS=$(( N_INST * N_PARAMS * N_SEEDS ))

# ── Create experiment directory ───────────────────────────────────────────────
if [[ -z "$OUTDIR" ]]; then
  TS=$(date +%Y_%m_%d)
  OUTDIR="${HOME}/experiments/cbc/lp_relaxation_${TS}"
fi
mkdir -p "$OUTDIR"

# ── Snapshot binary ───────────────────────────────────────────────────────────
EXP_TMPDIR=$(mktemp -d /tmp/cbc_lp_exp_XXXXXXXX)
cleanup() { rm -rf "${EXP_TMPDIR:-}"; }
trap cleanup EXIT
SNAP_BIN="$EXP_TMPDIR/cbc"
cp "$CBC_BIN" "$SNAP_BIN"
chmod +x "$SNAP_BIN"

# ── CSV header ────────────────────────────────────────────────────────────────
CSV_FILE="$OUTDIR/lp_results.csv"
if [[ ! -f "$CSV_FILE" ]]; then
  echo "instance,param_tag,seed,status,obj_from_log,obj_from_check,wall_seconds,check_result,exit_code" > "$CSV_FILE"
fi

# ── Write experiment_setup.md ─────────────────────────────────────────────────
SETUP_MD="$OUTDIR/experiment_setup.md"
{
  echo "# LP Relaxation Experiment Setup"
  echo ""
  echo "**Started:** $(date)"
  echo "**Experiment dir:** \`$OUTDIR\`"
  echo ""
  echo "## Parameters"
  echo ""
  echo "| Setting | Value |"
  echo "|:---|:---|"
  echo "| Binary (original) | \`$CBC_BIN\` |"
  echo "| Binary (snapshot) | \`$SNAP_BIN\` |"
  echo "| Params file | \`$PARAMS_FILE\` |"
  echo "| Instances dir | \`$INSTANCES_DIR\` |"
  echo "| Instances | $N_INST |"
  echo "| Parameter settings | $N_PARAMS |"
  echo "| Seeds | ${SEEDS} |"
  echo "| LP time limit (-sec) | ${LP_TIMELIMIT}s |"
  echo "| Hard kill after | ${KILL_AFTER}s (overtime: ${OVERTIME}s) |"
  echo "| Parallel jobs | $PARALLEL |"
  echo "| Total jobs | $N_JOBS |"
  echo ""
  echo "## Parameter Settings"
  echo ""
  echo "| Tag | CBC Parameters |"
  echo "|:---|:---|"
  for i in "${!PARAM_TAGS[@]}"; do
    echo "| ${PARAM_TAGS[$i]} | \`${PARAM_ARGS[$i]}\` |"
  done
  echo ""
  echo "## Build Info"
  echo ""
  echo "\`\`\`"
  "$CBC_BIN" -quit 2>&1 | head -3 || true
  echo "\`\`\`"
  echo ""
  echo "## Hardware"
  echo ""
  echo "\`\`\`"
  echo "Hostname: $(hostname 2>/dev/null || echo '?')"
  echo "Kernel:   $(uname -r 2>/dev/null || echo '?')"
  if command -v lscpu &>/dev/null; then
    lscpu 2>/dev/null | grep -E 'Model name|CPU\(s\):|Thread|Core|Socket|L3 cache'
  fi
  free -h 2>/dev/null | head -2
  echo "\`\`\`"
} > "$SETUP_MD"

# Also copy the params file for reproducibility
cp "$PARAMS_FILE" "$OUTDIR/lp_params.txt"

# ── Print summary ─────────────────────────────────────────────────────────────
echo "═══════════════════════════════════════════════════════════════"
echo "  LP Relaxation Experiment"
echo "───────────────────────────────────────────────────────────────"
echo "  Binary:       $CBC_BIN"
echo "  Snapshot:     $SNAP_BIN"
echo "  Instances:    $N_INST (from $INSTANCES_DIR)"
echo "  Params:       $N_PARAMS settings"
echo "  Seeds:        ${SEEDS}"
echo "  LP timelimit: ${LP_TIMELIMIT}s via -sec (+${OVERTIME}s overtime)"
echo "  Parallel:     $PARALLEL jobs"
echo "  Total jobs:   $N_JOBS"
echo "  Output:       $OUTDIR"
echo "═══════════════════════════════════════════════════════════════"
echo ""

# ── Generate job list ─────────────────────────────────────────────────────────
# Format: instance_path|param_tag|cbc_params|seed
JOBLIST="$EXP_TMPDIR/jobs.txt"
> "$JOBLIST"
for inst in "${INSTANCES[@]}"; do
  for i in "${!PARAM_TAGS[@]}"; do
    for seed in "${SEED_LIST[@]}"; do
      echo "${inst}|${PARAM_TAGS[$i]}|${PARAM_ARGS[$i]}|${seed}" >> "$JOBLIST"
    done
  done
done

# Shuffle job list to interleave instances and reduce peak memory consumption
# (avoids running many parameter settings for the same large instance back-to-back)
shuf "$JOBLIST" -o "$JOBLIST"

# Count how many are already done (resumability)
DONE=0
while IFS='|' read -r inst tag _ seed; do
  iname=$(basename "$inst" .mps.gz)
  [[ -f "$OUTDIR/${iname}_${tag}_s${seed}_fpp.result" ]] && DONE=$((DONE + 1))
done < "$JOBLIST"

if [[ $DONE -gt 0 ]]; then
  echo "  Resuming: $DONE/$N_JOBS jobs already completed, $(( N_JOBS - DONE )) remaining"
  echo ""
fi

if [[ $DRY_RUN -eq 1 ]]; then
  echo "DRY RUN — first 20 jobs (shuffled order):"
  head -20 "$JOBLIST" | while IFS='|' read -r inst tag params seed; do
    echo "  $(basename "$inst" .mps.gz)  tag=$tag  seed=$seed  params='$params'"
  done
  REMAINING=$(( N_JOBS - DONE ))
  echo "  ... ($REMAINING jobs would run)"
  exit 0
fi

# ── Export environment for worker ─────────────────────────────────────────────
export CBC_BIN="$SNAP_BIN"
export EXP_DIR="$OUTDIR"
export LP_TIMELIMIT
export KILL_AFTER
export CSV_FILE

# ── Run via GNU parallel ─────────────────────────────────────────────────────
echo "Starting at $(date) ..."
echo ""

PARALLEL_OPTS=(--jobs "$PARALLEL" --line-buffer --joblog "$OUTDIR/parallel.log")
[[ -t 1 ]] && PARALLEL_OPTS+=(--bar)

parallel "${PARALLEL_OPTS[@]}" "${SCRIPT_DIR}/run_one_lp.sh" < "$JOBLIST"

# ── Final summary ────────────────────────────────────────────────────────────
echo ""
echo "═══════════════════════════════════════════════════════════════"
echo "  Experiment complete at $(date)"
echo "───────────────────────────────────────────────────────────────"

# Count results by status
if [[ -f "$CSV_FILE" ]]; then
  echo "  Results by status:"
  tail -n +2 "$CSV_FILE" | cut -d',' -f4 | sort | uniq -c | sort -rn | while read -r cnt st; do
    printf "    %-20s %d\n" "$st" "$cnt"
  done
fi

TOTAL_RESULTS=$(( $(wc -l < "$CSV_FILE") - 1 ))
echo ""
echo "  Total results:  $TOTAL_RESULTS / $N_JOBS"
echo "  CSV:            $CSV_FILE"
echo "  Parallel log:   $OUTDIR/parallel.log"
echo "  Setup:          $SETUP_MD"
echo "═══════════════════════════════════════════════════════════════"
