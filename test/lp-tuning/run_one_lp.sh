#!/usr/bin/env bash
# run_one_lp.sh — Worker: run one LP relaxation solve job.
# Called by run_lp_experiments.sh via GNU parallel.
#
# Usage: run_one_lp.sh "instance_path|param_tag|cbc_params|seed"
#
# Environment (set by the orchestrator):
#   CBC_BIN       — path to snapshotted cbc binary
#   EXP_DIR       — experiment output directory
#   LP_TIMELIMIT  — LP time limit in seconds (-lpsec)
#   KILL_AFTER    — hard-kill timeout for the process
#   CSV_FILE      — path to results CSV (append, with flock)
#
# Output files per job (in $EXP_DIR):
#   {instance}_{tag}_s{seed}_fpp.sol   — LP solution
#   {instance}_{tag}_s{seed}_fpp.bas   — LP basis
#   {instance}_{tag}_s{seed}_fpp.txt   — checkSolution validation
#   {instance}_{tag}_s{seed}_fpp.log   — full CBC stdout
#   {instance}_{tag}_s{seed}_fpp.error — error output (deleted if empty)
#   {instance}_{tag}_s{seed}_fpp.result — completion marker for resumability
#
# Always exits 0 so parallel never aborts the batch.

set -uo pipefail

JOB="$1"
IFS='|' read -r INST_PATH PARAM_TAG CBC_PARAMS SEED <<< "$JOB"

INAME=$(basename "$INST_PATH" .mps.gz)
PREFIX="${INAME}_${PARAM_TAG}_s${SEED}_fpp"

SOL_FILE="${EXP_DIR}/${PREFIX}.sol"
BAS_FILE="${EXP_DIR}/${PREFIX}.bas"
CHK_FILE="${EXP_DIR}/${PREFIX}.txt"
LOG_FILE="${EXP_DIR}/${PREFIX}.log"
ERR_FILE="${EXP_DIR}/${PREFIX}.error"
RES_FILE="${EXP_DIR}/${PREFIX}.result"

# Skip if already completed (resumability)
[[ -f "$RES_FILE" ]] && exit 0

# Pre-create .bas file (basisOut bug: checks readability before writing)
touch "$BAS_FILE"

# Build command: params before -initialSolve, output actions after
read -ra PARAMS <<< "$CBC_PARAMS"
CMD=("$CBC_BIN" "$INST_PATH"
     -randomSeed "$SEED"
     -sec "$LP_TIMELIMIT")
[[ ${#PARAMS[@]} -gt 0 ]] && CMD+=("${PARAMS[@]}")
CMD+=(-initialSolve
     -writeSolution "$SOL_FILE"
     -basisOut "$BAS_FILE"
     -checkSolution "$CHK_FILE")

# Run with hard timeout
START_NS=$(date +%s%N)
timeout --kill-after=30 "$KILL_AFTER" "${CMD[@]}" > "$LOG_FILE" 2>&1
EXIT_CODE=$?
END_NS=$(date +%s%N)
WALL_MS=$(( (END_NS - START_NS) / 1000000 ))
WALL_S=$(awk -v ms="$WALL_MS" 'BEGIN { printf "%.3f", ms/1000 }')

# Detect errors in log
HAS_ERROR=0
ERROR_MSG=""
if [[ $EXIT_CODE -eq 124 ]]; then
  HAS_ERROR=1
  ERROR_MSG="KILLED: exceeded hard timeout ${KILL_AFTER}s"
elif [[ $EXIT_CODE -ne 0 ]]; then
  HAS_ERROR=1
  ERROR_MSG="EXIT_CODE=$EXIT_CODE"
fi
# Check for CBC error messages in stdout
if grep -qiE 'ERROR|Aborted|SIGABRT|SIGSEGV|infeasible' "$LOG_FILE" 2>/dev/null; then
  if grep -qi 'ERROR' "$LOG_FILE" 2>/dev/null; then
    HAS_ERROR=1
    ERROR_MSG="${ERROR_MSG:+$ERROR_MSG; }$(grep -i 'ERROR' "$LOG_FILE" | head -3 | tr '\n' ' ')"
  fi
fi

# Write .error file (delete if no errors)
if [[ $HAS_ERROR -eq 1 ]]; then
  {
    echo "exit_code=$EXIT_CODE"
    echo "wall_seconds=$WALL_S"
    echo "$ERROR_MSG"
  } > "$ERR_FILE"
else
  rm -f "$ERR_FILE"
fi

# Extract objective from log
OBJ=""
OBJ=$(grep -oP '(?i)optimal.*objective value \K[-+\d.eE]+' "$LOG_FILE" | tail -1 || true)
[[ -z "$OBJ" ]] && OBJ=$(grep -oP '✔ Optimal — Obj: \K[-+\d.eE]+' "$LOG_FILE" | tail -1 || true)

# Extract status from log
STATUS="UNKNOWN"
if [[ $EXIT_CODE -eq 124 ]]; then
  STATUS="TIMEOUT_KILLED"
elif grep -qi '✔ Optimal' "$LOG_FILE" 2>/dev/null; then
  STATUS="OPTIMAL"
elif grep -qi 'infeasible' "$LOG_FILE" 2>/dev/null; then
  STATUS="INFEASIBLE"
elif grep -qi 'stopped on time' "$LOG_FILE" 2>/dev/null; then
  STATUS="TIMEOUT"
else
  STATUS="ERROR"
fi

# Extract check_solution result
CHK_RESULT=""
if [[ -f "$CHK_FILE" ]]; then
  CHK_FEASIBLE=$(awk -F'\t' '$1=="lp_feasible"{print $2}' "$CHK_FILE" || true)
  CHK_PRIMAL=$(awk -F'\t' '$1=="largest_primal_error"{print $2}' "$CHK_FILE" || true)
  CHK_DUAL=$(awk -F'\t' '$1=="largest_dual_error"{print $2}' "$CHK_FILE" || true)
  CHK_OBJ=$(awk -F'\t' '$1=="objective"{print $2}' "$CHK_FILE" || true)
  CHK_RESULT="${CHK_FEASIBLE:-NA};primal=${CHK_PRIMAL:-NA};dual=${CHK_DUAL:-NA}"
fi

# Clean up empty output files
[[ -f "$SOL_FILE" && ! -s "$SOL_FILE" ]] && rm -f "$SOL_FILE"
[[ -f "$BAS_FILE" && ! -s "$BAS_FILE" ]] && rm -f "$BAS_FILE"
[[ -f "$CHK_FILE" && ! -s "$CHK_FILE" ]] && rm -f "$CHK_FILE"

# Write .result marker (resumability)
echo "DONE" > "$RES_FILE"

# Append to CSV (atomic via flock)
CSV_LINE="${INAME},${PARAM_TAG},${SEED},${STATUS},${OBJ:-NA},${CHK_OBJ:-NA},${WALL_S},${CHK_RESULT:-NA},${EXIT_CODE}"
(
  flock -x 200
  echo "$CSV_LINE" >> "$CSV_FILE"
) 200>"${CSV_FILE}.lock"

exit 0
