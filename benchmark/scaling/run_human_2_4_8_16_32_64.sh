#!/usr/bin/env bash
set -euo pipefail

# Benchmark KLOE compression on 2/4/8/16/32/64/128 human inputs from logenhuman/
# Usage:
#   bash benchmark/scaling/run_human_2_4_8_16_32_64.sh
#
# Optional env vars:
#   BIN=target/release/kloe
#   HUMAN_DIR=logenhuman
#   OUT_BASE=/tmp/kloe_human_scaling_YYYYmmdd_HHMMSS
#   THREADS=32
#   MEMORY_GB=40
#   PARTITION_POWER=10
#   K=31
#   M=7
#   TIG_MODE=simplitig   # one of: simplitig, unitig, matchtig, eulertig
#   BUILD_RELEASE=1      # 1 => cargo build -r before running

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "$ROOT_DIR"

SIZES=(2 4 8 16 32 64 128)
BIN="${BIN:-$ROOT_DIR/target/release/kloe}"
HUMAN_DIR="${HUMAN_DIR:-$ROOT_DIR/logenhuman}"
OUT_BASE="${OUT_BASE:-/tmp/kloe_human_scaling_$(date +%Y%m%d_%H%M%S)}"
THREADS="${THREADS:-32}"
MEMORY_GB="${MEMORY_GB:-40}"
PARTITION_POWER="${PARTITION_POWER:-10}"
K="${K:-31}"
M="${M:-7}"
TIG_MODE="${TIG_MODE:-simplitig}"
BUILD_RELEASE="${BUILD_RELEASE:-1}"

if [[ "$BUILD_RELEASE" == "1" ]]; then
  echo "[build] cargo build -r"
  cargo build -r
fi

if [[ ! -x "$BIN" ]]; then
  echo "ERROR: binary not found/executable: $BIN" >&2
  exit 1
fi

if [[ ! -d "$HUMAN_DIR" ]]; then
  echo "ERROR: human directory not found: $HUMAN_DIR" >&2
  exit 1
fi

if [[ ! -x /usr/bin/time ]]; then
  echo "ERROR: /usr/bin/time not found (required for max RSS / wall time capture)." >&2
  exit 1
fi

mapfile -t GENOMES < <(find "$HUMAN_DIR" -maxdepth 1 -type f -name '*.unitigs.fa.zst' | sort)

if [[ "${#GENOMES[@]}" -lt 128 ]]; then
  echo "ERROR: need at least 128 *.unitigs.fa.zst files in $HUMAN_DIR, found ${#GENOMES[@]}" >&2
  exit 1
fi

mkdir -p "$OUT_BASE/fof" "$OUT_BASE/logs" "$OUT_BASE/out" "$OUT_BASE/tmp"

SUMMARY_TSV="$OUT_BASE/summary.tsv"
printf "n_inputs\texit_code\telapsed_sec\tmax_rss_kb\tarchive_bytes\tfof\tlog\ttime\n" > "$SUMMARY_TSV"

tig_flag=()
case "$TIG_MODE" in
  simplitig) tig_flag=() ;;
  unitig) tig_flag=(--unitig) ;;
  matchtig) tig_flag=(--matchtig) ;;
  eulertig) tig_flag=(--eulertig) ;;
  *)
    echo "ERROR: TIG_MODE must be one of: simplitig|unitig|matchtig|eulertig" >&2
    exit 1
    ;;
esac

to_seconds() {
  local raw="$1"
  awk -v t="$raw" 'BEGIN {
    n=split(t,a,":");
    if (n==3) printf "%.3f", (a[1]*3600 + a[2]*60 + a[3]);
    else if (n==2) printf "%.3f", (a[1]*60 + a[2]);
    else printf "%.3f", a[1];
  }'
}

echo "[info] output base: $OUT_BASE"
echo "[info] human dir:   $HUMAN_DIR"
echo "[info] binary:      $BIN"
echo "[info] mode:        $TIG_MODE"
echo "[info] threads:     $THREADS"
echo "[info] memory GB:   $MEMORY_GB"

echo "[info] selected first 128 genomes (sorted):"
printf '  %s\n' "${GENOMES[@]:0:128}"

for n in "${SIZES[@]}"; do
  fof="$OUT_BASE/fof/fof_${n}.txt"
  printf '%s\n' "${GENOMES[@]:0:$n}" > "$fof"

  run_out="$OUT_BASE/out/n${n}"
  run_tmp="$OUT_BASE/tmp/n${n}"
  log="$OUT_BASE/logs/run_${n}.log"
  tlog="$OUT_BASE/logs/run_${n}_time.txt"

  rm -rf "$run_out" "$run_tmp"
  mkdir -p "$run_out" "$run_tmp"

  cmd=(
    "$BIN" compress
    -i "$fof"
    -o "$run_out/"
    -d "$run_tmp"
    -t "$THREADS"
    -r "$MEMORY_GB"
    -P "$PARTITION_POWER"
    -k "$K"
    -m "$M"
    "${tig_flag[@]}"
  )

  echo "[run n=$n] ${cmd[*]}"
  set +e
  /usr/bin/time -v "${cmd[@]}" >"$log" 2>"$tlog"
  exit_code=$?
  set -e

  wall_raw="$(awk -F': ' '/Elapsed \(wall clock\) time/ {print $2}' "$tlog" | tail -n1)"
  rss_kb="$(awk -F': ' '/Maximum resident set size \(kbytes\)/ {print $2}' "$tlog" | tail -n1)"

  if [[ -n "$wall_raw" ]]; then
    elapsed_sec="$(to_seconds "$wall_raw")"
  else
    elapsed_sec="NA"
  fi

  if [[ -z "$rss_kb" ]]; then
    rss_kb="NA"
  fi

  if [[ -d "$run_out" ]]; then
    archive_bytes="$(du -sb "$run_out" | awk '{print $1}')"
  else
    archive_bytes="0"
  fi

  printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
    "$n" "$exit_code" "$elapsed_sec" "$rss_kb" "$archive_bytes" \
    "$fof" "$log" "$tlog" >> "$SUMMARY_TSV"

  echo "[done n=$n] exit=$exit_code elapsed_s=$elapsed_sec max_rss_kb=$rss_kb archive_bytes=$archive_bytes"
done

echo ""
echo "Benchmark complete"
echo "Summary: $SUMMARY_TSV"
