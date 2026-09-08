#!/usr/bin/env bash
# fj-tune-sweep.sh — sweep Feasibility Jump / Feasibility Pump configurations
# against a root-fixture collection (see ROOT-FIXTURES.md, mip-root-replay)
# and report, per configuration, how often a root-only replay finds *any*
# incumbent and how good it is relative to bks.tsv.
#
# Each (instance, config) pair is one `mip-root-replay` call -- millisecond-
# scale since the fixture skips pre-processing/root-LP -- so the whole sweep
# (hundreds of instances x several configs) runs in minutes via GNU parallel,
# not hours.
#
# Usage:
#   ./fj-tune-sweep.sh --configs=configs.tsv --outdir=/path/to/results \
#     [--fixture-dir=PATH] [--data-dir=PATH] [--nodes=N] [--sec=N] [--jobs=N]
#
# configs.tsv format: one config per line, tab-separated
#   <name>\t<extra mip-root-replay args, or - for none>
# Blank lines and lines starting with # are skipped.
#
# Output: outdir/<config>.tsv (one row per instance: instance, found,
# best, bound, gap_bks_pct, bbTime) and outdir/summary.tsv (one row per
# config: aggregate stats).
set -uo pipefail

SCRIPT_DIR="$(cd "$(dirname "$(readlink -f "${BASH_SOURCE[0]}")")" && pwd)"
REPLAY="$SCRIPT_DIR/mip-root-replay"

# Worker mode: this script re-invokes itself once per (config, instance) pair
# via GNU parallel. Checked first, before normal option parsing, since the
# lone positional argument in this mode is the instance name, not an option.
if [[ -n "${WORKER_CONFIG_NAME:-}" ]]; then
  name="$WORKER_CONFIG_NAME"
  args="$WORKER_CONFIG_ARGS"
  inst="$1"
  [[ "$args" == "-" ]] && args=""
  # Hard wall-clock backstop: mip-root-replay's own --sec only bounds the
  # B&B/CbcModel time check, which some heuristics (FPump/FJ on very large
  # or numerically hard fixtures) can overrun between polls. `timeout` here
  # guarantees no single instance can stall the whole parallel sweep; give
  # it generous headroom (--sec plus a flat 60s) rather than killing right
  # at the target so we still capture results for borderline-slow runs.
  hard_timeout=$(awk -v s="$SECS_ENV" 'BEGIN { printf "%d", s + 60 }')
  # shellcheck disable=SC2086
  out=$(timeout --kill-after=10 "${hard_timeout}s" "$REPLAY" "$inst" --fixture-dir="$FIXTURE_DIR_ENV" --data-dir="$DATA_DIR_ENV" \
    --nodes="$NODES_ENV" --sec="$SECS_ENV" --quiet $args 2>/dev/null)
  best=$(echo "$out" | sed -n 's/.*best=\([^ ]*\).*/\1/p' | tail -1)
  bound=$(echo "$out" | sed -n 's/.*bound=\([^ ]*\).*/\1/p' | tail -1)
  bbtime=$(echo "$out" | sed -n 's/.*bbTime=\([^ s]*\)s.*/\1/p' | tail -1)
  bks=$(echo "$out" | sed -n 's/.*bks=\([^ ]*\).*/\1/p' | tail -1)
  optimal=$(echo "$out" | sed -n 's/.*optimal=\([01]\).*/\1/p' | tail -1)
  found=0
  [[ -n "$best" && "$best" != "none" ]] && found=1
  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$inst" "$found" "${best:-}" "${bound:-}" "${bks:-}" "${bbtime:-}" "${optimal:-0}"
  exit 0
fi

CONFIGS=""
OUTDIR=""
FIXTURE_DIR="${CBC_ROOT_FIXTURE_DIR:-$HOME/instances/mip-sanity-data/rootFixtures}"
DATA_DIR=""
NODES=1
SECS=30
JOBS="$(nproc)"

for arg in "$@"; do
  case "$arg" in
    --configs=*) CONFIGS="${arg#*=}" ;;
    --outdir=*) OUTDIR="${arg#*=}" ;;
    --fixture-dir=*) FIXTURE_DIR="${arg#*=}" ;;
    --data-dir=*) DATA_DIR="${arg#*=}" ;;
    --nodes=*) NODES="${arg#*=}" ;;
    --sec=*) SECS="${arg#*=}" ;;
    --jobs=*) JOBS="${arg#*=}" ;;
    -h|--help)
      sed -n '2,/^set -uo/p' "$0" | sed '$d' | sed 's/^# \{0,1\}//'
      exit 0
      ;;
    *) echo "Unknown option: $arg" >&2; exit 2 ;;
  esac
done

if [[ -z "$CONFIGS" || -z "$OUTDIR" ]]; then
  echo "Usage: $0 --configs=configs.tsv --outdir=DIR [options]" >&2
  exit 2
fi
if [[ ! -x "$REPLAY" ]]; then
  echo "mip-root-replay not built; run 'make mip-root-replay' in $SCRIPT_DIR first" >&2
  exit 2
fi
if ! command -v parallel >/dev/null 2>&1; then
  echo "GNU parallel is required (apt install parallel)." >&2
  exit 2
fi

mkdir -p "$OUTDIR"

# Instance list: every fixture with a .meta file, by stem basename.
mapfile -t INSTANCES < <(
  for f in "$FIXTURE_DIR"/*.meta; do
    [[ -e "$f" ]] || continue
    b="$(basename "$f" .meta)"
    echo "${b%.root}"
  done | sort
)
echo "Sweeping ${#INSTANCES[@]} fixture instance(s) across $(grep -cvE '^\s*(#|$)' "$CONFIGS") config(s), jobs=$JOBS"

export FIXTURE_DIR_ENV="$FIXTURE_DIR" DATA_DIR_ENV="$DATA_DIR" NODES_ENV="$NODES" SECS_ENV="$SECS"

INSTLIST="$(mktemp)"
printf '%s\n' "${INSTANCES[@]}" > "$INSTLIST"

while IFS=$'\t' read -r cname cargs; do
  [[ -z "$cname" || "$cname" == \#* ]] && continue
  echo "=== $cname ($cargs) ==="
  out_tsv="$OUTDIR/$cname.tsv"
  echo -e "instance\tfound\tbest\tbound\tbks\tbbtime\toptimal" > "$out_tsv"
  WORKER_CONFIG_NAME="$cname" WORKER_CONFIG_ARGS="$cargs" \
    parallel --jobs "$JOBS" -a "$INSTLIST" "$SCRIPT_DIR/$(basename "$0")" >> "$out_tsv"
done < "$CONFIGS"

rm -f "$INSTLIST"

# Aggregate summary across all configs.
python3 "$SCRIPT_DIR/fj-tune-summarize.py" "$OUTDIR" "$CONFIGS" > "$OUTDIR/summary.tsv"
cat "$OUTDIR/summary.tsv"
