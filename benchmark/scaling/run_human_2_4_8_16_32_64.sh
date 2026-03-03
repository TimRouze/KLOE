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
#   SIZES="2 4 8 16 32 64 128"
#   THREADS=32
#   MEMORY_GB=40
#   HARD_MEM_GB=0      # 0 disables hard limit; enforced via ulimit -v
#   SYSTEMD_MEMORY_MAX=40G   # if set, run each job under systemd-run --user --scope -p MemoryMax=...
#   K=31
#   TIG_MODE=simplitig   # one of: simplitig, unitig, matchtig, eulertig
#   BUILD_RELEASE=1      # 1 => cargo build -r before running
#   CLEAN_RUN_ARTIFACTS=1 # 1 => remove per-run out/tmp after metrics are captured

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "$ROOT_DIR"

SIZES_RAW="${SIZES:-2 4 8 16 32 64 128}"
read -r -a SIZES <<<"$SIZES_RAW"
BIN="${BIN:-$ROOT_DIR/target/release/kloe}"
HUMAN_DIR="${HUMAN_DIR:-$ROOT_DIR/loganhuman}"
OUT_BASE="${OUT_BASE:-/tmp/kloe_human_scaling_$(date +%Y%m%d_%H%M%S)}"
THREADS="${THREADS:-32}"
MEMORY_GB="${MEMORY_GB:-40}"
HARD_MEM_GB="${HARD_MEM_GB:-0}"
SYSTEMD_MEMORY_MAX="${SYSTEMD_MEMORY_MAX:-}"
K="${K:-31}"
TIG_MODE="${TIG_MODE:-simplitig}"
BUILD_RELEASE="${BUILD_RELEASE:-1}"
CLEAN_RUN_ARTIFACTS="${CLEAN_RUN_ARTIFACTS:-1}"

CURRENT_RUN_OUT=""
CURRENT_RUN_TMP=""

cleanup_current_run() {
  if [[ -n "$CURRENT_RUN_OUT" ]]; then
    rm -rf "$CURRENT_RUN_OUT"
  fi
  if [[ -n "$CURRENT_RUN_TMP" ]]; then
    rm -rf "$CURRENT_RUN_TMP"
  fi
  CURRENT_RUN_OUT=""
  CURRENT_RUN_TMP=""
}

cleanup_on_exit() {
  cleanup_current_run
}

trap cleanup_on_exit EXIT INT TERM

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
GENOME_GLOB="*.unitigs.fa.zst"
if [[ "${#GENOMES[@]}" -lt 128 ]]; then
  mapfile -t GENOMES < <(find "$HUMAN_DIR" -maxdepth 1 -type f -name '*.u.fa.zst' | sort)
  GENOME_GLOB="*.u.fa.zst"
fi

if [[ "${#GENOMES[@]}" -lt 128 ]]; then
  echo "ERROR: need at least 128 human files matching '*.unitigs.fa.zst' or '*.u.fa.zst' in $HUMAN_DIR, found ${#GENOMES[@]}" >&2
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
echo "[info] hard mem GB: $HARD_MEM_GB"
echo "[info] systemd mem: ${SYSTEMD_MEMORY_MAX:-disabled}"
echo "[info] sizes:       ${SIZES[*]}"
echo "[info] cleanup:     CLEAN_RUN_ARTIFACTS=$CLEAN_RUN_ARTIFACTS"
echo "[info] input glob:  $GENOME_GLOB"

echo "[info] selected first 128 genomes (sorted):"
printf '  %s\n' "${GENOMES[@]:0:128}"

for n in "${SIZES[@]}"; do
  fof="$OUT_BASE/fof/fof_${n}.txt"
  printf '%s\n' "${GENOMES[@]:0:$n}" > "$fof"

  run_out="$OUT_BASE/out/n${n}"
  run_tmp="$OUT_BASE/tmp/n${n}"
  log="$OUT_BASE/logs/run_${n}.log"
  tlog="$OUT_BASE/logs/run_${n}_time.txt"
  CURRENT_RUN_OUT="$run_out"
  CURRENT_RUN_TMP="$run_tmp"

  rm -rf "$run_out" "$run_tmp"
  mkdir -p "$run_out" "$run_tmp"

  cmd=(
    "$BIN" compress
    -i "$fof"
    -o "$run_out/"
    -d "$run_tmp"
    -t "$THREADS"
    -r "$MEMORY_GB"
    -k "$K"
    "${tig_flag[@]}"
  )

  echo "[run n=$n] ${cmd[*]}"
  set +e
  if [[ -n "$SYSTEMD_MEMORY_MAX" ]]; then
    if ! command -v systemd-run >/dev/null 2>&1; then
      echo "ERROR: SYSTEMD_MEMORY_MAX is set but systemd-run is unavailable." >&2
      exit 1
    fi
    /usr/bin/time -v systemd-run --user --scope -p "MemoryMax=$SYSTEMD_MEMORY_MAX" "${cmd[@]}" >"$log" 2>"$tlog"
  elif [[ "$HARD_MEM_GB" -gt 0 ]]; then
    hard_mem_kb=$((HARD_MEM_GB * 1024 * 1024))
    (
      ulimit -Sv "$hard_mem_kb"
      /usr/bin/time -v "${cmd[@]}"
    ) >"$log" 2>"$tlog"
  else
    /usr/bin/time -v "${cmd[@]}" >"$log" 2>"$tlog"
  fi
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

  if [[ "$CLEAN_RUN_ARTIFACTS" == "1" ]]; then
    cleanup_current_run
  else
    CURRENT_RUN_OUT=""
    CURRENT_RUN_TMP=""
  fi
done

if [[ "$CLEAN_RUN_ARTIFACTS" == "1" ]]; then
  rmdir "$OUT_BASE/out" "$OUT_BASE/tmp" 2>/dev/null || true
fi

echo ""
echo "Benchmark complete"
echo "Summary: $SUMMARY_TSV"
