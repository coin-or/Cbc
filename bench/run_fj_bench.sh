#!/usr/bin/env bash
# run_fj_bench.sh — Run FJ benchmark on random instances and summarize results.
#
# Usage:
#   ./run_fj_bench.sh [options]
#
# Options:
#   -n N          Number of instances to sample (default: 28)
#   -e N          Effort budget (default: 500000000)
#   -s N          Stall limit (default: 0 = off)
#   -t N          Time limit per instance in seconds (default: 60)
#   --maxsol N    Max solutions per instance (default: 3)
#   --seed N      FJ seed (default: 0)
#   --instances D Instance directory (default: $MIPSTER_INSTANCES/miplib/2017+spp)
#   --all         Run all instances (no sampling)
#   -o FILE       Output TSV file (default: stdout + summary)
#   -j N          Parallel jobs (default: nproc)

set -euo pipefail

NINSTANCES=28
EFFORT=500000000
STALL=0
TIMELIMIT=60
MAXSOL=3
SEED=0
INST_DIR="${MIPSTER_INSTANCES:-$HOME/inst}/miplib/2017+spp"
OUTFILE=""
JOBS=$(nproc)
SAMPLE=1

while [[ $# -gt 0 ]]; do
  case "$1" in
    -n)       NINSTANCES="$2"; shift 2 ;;
    -e)       EFFORT="$2"; shift 2 ;;
    -s)       STALL="$2"; shift 2 ;;
    -t)       TIMELIMIT="$2"; shift 2 ;;
    --maxsol) MAXSOL="$2"; shift 2 ;;
    --seed)   SEED="$2"; shift 2 ;;
    --instances) INST_DIR="$2"; shift 2 ;;
    --all)    SAMPLE=0; shift ;;
    -o)       OUTFILE="$2"; shift 2 ;;
    -j)       JOBS="$2"; shift 2 ;;
    *)        echo "Unknown option: $1" >&2; exit 1 ;;
  esac
done

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
FJ_BENCH="${SCRIPT_DIR}/fj_bench"

if [[ ! -x "$FJ_BENCH" ]]; then
  echo "ERROR: fj_bench not found. Build it first: cd bench && make PREFIX=~/prog/mipster CXX=g++-15" >&2
  exit 1
fi

# Collect instances
if (( SAMPLE )); then
  INSTANCES=$(ls "$INST_DIR"/*.mps.gz | shuf -n "$NINSTANCES")
else
  INSTANCES=$(ls "$INST_DIR"/*.mps.gz)
  NINSTANCES=$(echo "$INSTANCES" | wc -l)
fi

echo "FJ Benchmark: $NINSTANCES instances, effort=$EFFORT, stall=$STALL, timelimit=${TIMELIMIT}s, jobs=$JOBS" >&2

# Run in parallel, collect TSV (fj_bench returns 2 for NOTFOUND, ignore)
RESULTS=$(echo "$INSTANCES" | parallel --will-cite -j"$JOBS" \
  "$FJ_BENCH --effort $EFFORT --stall $STALL --maxsol $MAXSOL --timelimit $TIMELIMIT --seed $SEED --loglevel 0 {} 2>/dev/null || true")

# Save raw TSV if requested
if [[ -n "$OUTFILE" ]]; then
  printf "instance\tstatus\tobj\tlp_time\tfj_time\ttime_first\teffort_first\ttotal_effort\tsolutions\teffort_budget\ttimelimit\tseed\n" > "$OUTFILE"
  echo "$RESULTS" >> "$OUTFILE"
  echo "Raw results saved to: $OUTFILE" >&2
fi

# Summary
FOUND=$(echo "$RESULTS" | grep -c $'^\S\+\tFOUND' || true)
NOTFOUND=$(echo "$RESULTS" | grep -c $'^\S\+\tNOTFOUND' || true)
TOTAL=$((FOUND + NOTFOUND))

echo "" >&2
echo "═══════════════════════════════════════════════════════════" >&2
echo "  FJ Benchmark Summary" >&2
echo "═══════════════════════════════════════════════════════════" >&2
echo "  Instances:    $TOTAL" >&2
echo "  Found:        $FOUND ($(awk "BEGIN{printf \"%.1f\", 100*$FOUND/$TOTAL}")%)" >&2
echo "  Not found:    $NOTFOUND" >&2
echo "" >&2

if (( FOUND > 0 )); then
  echo "  Time to first solution (found instances):" >&2
  echo "$RESULTS" | awk -F'\t' '$2=="FOUND" {print $6}' | sort -n | awk '
    BEGIN { n=0; sum=0 }
    { a[n++]=$1; sum+=$1 }
    END {
      printf "    min:    %.4fs\n", a[0]
      printf "    median: %.4fs\n", a[int(n/2)]
      printf "    mean:   %.4fs\n", sum/n
      printf "    max:    %.4fs\n", a[n-1]
    }' >&2

  echo "" >&2
  echo "  Effort to first solution:" >&2
  echo "$RESULTS" | awk -F'\t' '$2=="FOUND" {print $7}' | sort -n | awk '
    BEGIN { n=0 }
    { a[n++]=$1 }
    END {
      printf "    min:    %.2fM\n", a[0]/1e6
      printf "    median: %.2fM\n", a[int(n/2)]/1e6
      printf "    max:    %.2fM\n", a[n-1]/1e6
    }' >&2
fi

echo "" >&2
echo "  FJ wall-clock time (all instances):" >&2
echo "$RESULTS" | awk -F'\t' '{print $5}' | sort -n | awk '
  BEGIN { n=0; sum=0 }
  { a[n++]=$1; sum+=$1 }
  END {
    printf "    min:    %.2fs\n", a[0]
    printf "    median: %.2fs\n", a[int(n/2)]
    printf "    mean:   %.2fs\n", sum/n
    printf "    max:    %.2fs\n", a[n-1]
    printf "    total:  %.1fs\n", sum
  }' >&2

echo "═══════════════════════════════════════════════════════════" >&2
echo "" >&2

# Print sorted results table
printf "\n%-25s %-10s %12s %8s %8s %10s %10s\n" "Instance" "Status" "Objective" "LP(s)" "FJ(s)" "T1st(s)" "E1st(M)"
printf "%-25s %-10s %12s %8s %8s %10s %10s\n" "-------------------------" "----------" "------------" "--------" "--------" "----------" "----------"
echo "$RESULTS" | sort -t$'\t' -k6 -n | awk -F'\t' '{
  name=$1; status=$2; obj=$3; lp=$4; fj=$5; t1=$6; e1=$7
  if (status=="FOUND") {
    printf "%-25s \033[32m%-10s\033[0m %12s %8.2f %8.2f %10.4f %10.2f\n", name, status, obj, lp, fj, t1, e1/1e6
  } else {
    printf "%-25s \033[31m%-10s\033[0m %12s %8.2f %8.2f %10s %10s\n", name, status, "-", lp, fj, "-", "-"
  }
}'
