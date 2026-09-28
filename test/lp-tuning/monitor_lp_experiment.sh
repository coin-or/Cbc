#!/usr/bin/env bash
# monitor_lp_experiment.sh — status snapshot of a running run_lp_experiments.sh
#
# Usage: monitor_lp_experiment.sh EXP_DIR [--watch SECONDS]
#
# Shows progress, result-status counts, peak-RSS distribution of finished
# jobs, system memory, and every running cbc job (elapsed time, current RSS),
# so the parallelism can be tuned while the experiment runs:
#   echo 48 > EXP_DIR/parallel_jobs    # takes effect as jobs finish

set -uo pipefail

EXP_DIR="${1:?usage: $0 EXP_DIR [--watch SECONDS]}"
EXP_DIR="$(cd "$EXP_DIR" && pwd)"
WATCH=0
[[ "${2:-}" == "--watch" ]] && WATCH="${3:-60}"

snapshot() {
  local csv="$EXP_DIR/lp_results.csv"
  local total
  total=$(grep -oP '^\| Total jobs \| \K[0-9]+' "$EXP_DIR/experiment_setup.md" 2>/dev/null || echo "?")
  local done_n=$(( $(wc -l < "$csv") - 1 ))
  echo "=== $(date '+%F %T')  $EXP_DIR"
  echo "Progress: $done_n / $total jobs    parallel_jobs=$(cat "$EXP_DIR/parallel_jobs" 2>/dev/null || echo ?)"
  echo
  echo "Status counts (by tag):"
  awk -F, 'NR>1 {c[$2" "$4]++} END {for (k in c) print c[k], k}' "$csv" \
    | sort -k2,2 -k1,1nr | awk '{printf "  %-22s %-16s %6d\n", $2, $3, $1}'
  echo
  echo "Peak RSS of finished jobs (MB):"
  awk -F, 'NR>1 && $10 ~ /^[0-9]+$/ {print $10}' "$csv" | sort -n | awk '
    function q(p,  i) { i = int(NR * p); if (i < 1) i = 1; return v[i] }
    {v[NR] = $1}
    END {
      if (NR == 0) { print "  (none yet)"; exit }
      printf "  n=%d  median=%d  p90=%d  p99=%d  max=%d\n", NR, q(0.5), q(0.9), q(0.99), v[NR]
    }'
  local memouts
  memouts=$(awk -F, 'NR>1 && $4=="MEMOUT"' "$csv" | wc -l)
  echo "  MEMOUT jobs: $memouts"
  echo
  echo "System memory:"
  free -g | sed 's/^/  /'
  echo
  echo "Running jobs (sorted by RSS):"
  local n=0 sum_kb=0
  local rows=""
  for pid in $(pgrep -u "$USER" -x cbc); do
    local cmd
    cmd=$(tr '\0' ' ' < "/proc/$pid/cmdline" 2>/dev/null) || continue
    [[ "$cmd" == *"$EXP_DIR/"* ]] || continue
    local rss et job
    rss=$(awk '/^VmRSS/ {print $2}' "/proc/$pid/status" 2>/dev/null)
    et=$(ps -o etimes= -p "$pid" 2>/dev/null | tr -d ' ')
    job=$(grep -oP -- '-writeSolution \S+/\K\S+(?=_fpp\.sol)' <<< "$cmd")
    n=$((n + 1)); sum_kb=$((sum_kb + ${rss:-0}))
    rows+="$(printf '%10d %8d %s' "${rss:-0}" "${et:-0}" "$job")"$'\n'
  done
  printf '  %8s %8s  %s\n' "RSS(MB)" "elapsed" "job"
  sort -rn <<< "$rows" | awk 'NF {printf "  %8d %7ds  %s\n", $1/1024, $2, $3}'
  echo "  running: $n   total RSS: $((sum_kb / 1024 / 1024)) GB"
}

if (( WATCH > 0 )); then
  while true; do clear; snapshot; sleep "$WATCH"; done
else
  snapshot
fi
