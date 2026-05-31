#!/usr/bin/env bash
# gen_test_solutions.sh — generate reference optimal solution files for all
# C-interface test fixtures in parallel.
#
# Solutions are written to test/fixtures/solutions/ (created if absent).
# Each .sol file can be passed to -debugCuts to catch cuts that incorrectly
# exclude an integer feasible solution.
#
# Usage (from the repo root):
#   cd test && bash gen_test_solutions.sh [timelimit_seconds]
#   # OR from repo root:
#   bash test/gen_test_solutions.sh [timelimit_seconds]
#
# Default time limit: 1800 seconds (30 minutes).
#
# Requirements:
#   - The test programs must be built (run 'make' in the test directory first).
#   - $MIPSTER_PREFIX must point to a built+installed mipster (for the binary).
#   - GNU parallel is recommended for concurrent MPS fixture solves.

set -euo pipefail

TIMELIMIT="${1:-1800}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TEST_DIR="$SCRIPT_DIR"
OUTDIR="$TEST_DIR/fixtures/solutions"

MIPSTER="${MIPSTER_PREFIX:-$HOME/prog/cbc}/bin/mipster"

echo "==================================================================="
echo "  gen_test_solutions.sh"
echo "  time limit : ${TIMELIMIT}s"
echo "  output dir : ${OUTDIR}"
echo "  mipster    : ${MIPSTER}"
echo "==================================================================="

mkdir -p "$OUTDIR"

# ── 1. MICLSP instances via the C-API generator ────────────────────────────

echo ""
echo "--- MICLSP instances (via gen_miclsp_solutions) ---"
cd "$TEST_DIR"

if [ ! -x "./gen_miclsp_solutions" ]; then
  echo "Building gen_miclsp_solutions..."
  make gen_miclsp_solutions
fi

./gen_miclsp_solutions "$OUTDIR" "$TIMELIMIT"

# ── 2. MPS fixture instances via the mipster binary ────────────────────────

if [ ! -x "$MIPSTER" ]; then
  echo ""
  echo "WARNING: mipster binary not found at $MIPSTER"
  echo "  Set MIPSTER_PREFIX or install mipster to also solve MPS fixtures."
  echo "  Skipping MPS fixture solves."
  exit 0
fi

echo ""
echo "--- MPS fixtures (via mipster binary in parallel) ---"

FIXTURE_DIR="$TEST_DIR/fixtures"

# Collect all MPS fixture files
mapfile -t MPS_FILES < <(ls "$FIXTURE_DIR"/*.mps.gz 2>/dev/null || true)

if [ ${#MPS_FILES[@]} -eq 0 ]; then
  echo "No MPS fixture files found in $FIXTURE_DIR"
else
  echo "Found ${#MPS_FILES[@]} MPS fixture files"

  solve_one() {
    local f="$1"
    local base
    base="$(basename "$f" .mps.gz)"
    local solpath="${OUTDIR}/${base}.sol"
    local logfile
    logfile="$(mktemp /tmp/mipster_XXXXXX.log)"

    if "$MIPSTER" "$f" -preprocess off -sec "$TIMELIMIT" -solve \
        -writeSolution "$solpath" > "$logfile" 2>&1; then
      # "Optimal - objective value X" header is written to the solution file
      if [ -s "$solpath" ] && grep -q "^Optimal" "$solpath"; then
        local obj
        obj=$(grep "^Optimal" "$solpath" | grep -oP '[\d.eE+\-]+$' | head -1)
        echo "  OPTIMAL  ${base}  obj=${obj}"
      elif [ -s "$solpath" ]; then
        echo "  FEASIBLE ${base}  (not proven optimal within ${TIMELIMIT}s)"
        rm -f "$solpath"  # don't keep non-optimal solutions
      else
        echo "  NO SOL   ${base}"
      fi
    else
      echo "  FAILED   ${base}  (mipster exited with error)"
    fi
    rm -f "$logfile"
  }

  export -f solve_one
  export MIPSTER OUTDIR TIMELIMIT

  if command -v parallel &>/dev/null; then
    printf '%s\n' "${MPS_FILES[@]}" | \
      parallel -j"$(nproc)" solve_one
  else
    for f in "${MPS_FILES[@]}"; do
      solve_one "$f"
    done
  fi
fi

echo ""
echo "==================================================================="
echo "  Done.  Solutions written to: $OUTDIR"
ls -1 "$OUTDIR"/*.sol 2>/dev/null | wc -l | xargs printf "  %d .sol files generated\n"
echo "==================================================================="
