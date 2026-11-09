#!/usr/bin/env bash
set -euo pipefail

# Run Salmonella targeted decompression + k-mer content checks.
# Default sizes are 2,4,8,16 and automatically include 32/64 when enough genomes exist.
# No parameters required.
#
# Optional env overrides:
#   BIN=/abs/path/to/kloe
#   SALMON_DIR=/abs/path/to/salmonella128
#   SIZES="2 4 8 16 32 64"   # optional explicit override
#   THREADS=32
#   MEMORY_GB=40
#   K=31
#   CHECK_JOBS=8
#   DECOMP_THREADS=1
#   OUT_BASE=/tmp/kloe_salmonella_checks_YYYYmmdd_HHMMSS
#   CLEAN_RUN_ARTIFACTS=1

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "$ROOT_DIR"

BIN="${BIN:-$ROOT_DIR/target/release/kloe}"
SALMON_DIR="${SALMON_DIR:-/home/nadine/Code/parsou/salmonella128}"
THREADS="${THREADS:-32}"
MEMORY_GB="${MEMORY_GB:-40}"
K="${K:-31}"
CPU_COUNT="$(getconf _NPROCESSORS_ONLN 2>/dev/null || echo 8)"
CHECK_JOBS="${CHECK_JOBS:-$CPU_COUNT}"
DECOMP_THREADS="${DECOMP_THREADS:-1}"
OUT_BASE="${OUT_BASE:-/tmp/kloe_salmonella_checks_$(date +%Y%m%d_%H%M%S)}"
CLEAN_RUN_ARTIFACTS="${CLEAN_RUN_ARTIFACTS:-1}"

if [[ ! -x "$BIN" ]]; then
  echo "[build] cargo build -r"
  cargo build -r
fi

if [[ ! -x "$BIN" ]]; then
  echo "ERROR: binary not found/executable: $BIN" >&2
  exit 1
fi
if [[ ! -d "$SALMON_DIR" ]]; then
  echo "ERROR: salmonella directory not found: $SALMON_DIR" >&2
  exit 1
fi
if ! command -v python3 >/dev/null 2>&1; then
  echo "ERROR: python3 is required." >&2
  exit 1
fi

mapfile -t INPUTS < <(find "$SALMON_DIR" -maxdepth 1 -type f -name '*.fna.gz' | sort)
if [[ "${#INPUTS[@]}" -lt 16 ]]; then
  echo "ERROR: need at least 16 '*.fna.gz' genomes in $SALMON_DIR, found ${#INPUTS[@]}" >&2
  exit 1
fi

if [[ -n "${SIZES:-}" ]]; then
  SIZES_RAW="$SIZES"
  read -r -a SIZES <<<"$SIZES_RAW"
else
  SIZES=(2 4 8 16)
  if [[ "${#INPUTS[@]}" -ge 32 ]]; then
    SIZES+=(32)
  fi
  if [[ "${#INPUTS[@]}" -ge 64 ]]; then
    SIZES+=(64)
  fi
fi

max_n=0
for n in "${SIZES[@]}"; do
  if [[ "$n" -gt "$max_n" ]]; then
    max_n="$n"
  fi
done
if [[ "$max_n" -gt "${#INPUTS[@]}" ]]; then
  echo "ERROR: requested max size n=$max_n but only ${#INPUTS[@]} '*.fna.gz' genomes found in $SALMON_DIR" >&2
  exit 1
fi

mkdir -p "$OUT_BASE"
SUMMARY="$OUT_BASE/summary.tsv"
printf "n\tstatus\tchecks_ok\tchecks_total\trun_dir\n" > "$SUMMARY"

VERIFY_PY="$OUT_BASE/verify_pair.py"
cat >"$VERIFY_PY" <<'PY'
#!/usr/bin/env python3
import gzip
import sys

def open_auto(path):
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "rt")

def revcomp(seq):
    trans = str.maketrans("ACGTacgt", "TGCAtgca")
    return seq.translate(trans)[::-1]

def canonical(kmer):
    rc = revcomp(kmer)
    k = kmer.upper()
    r = rc.upper()
    return k if k <= r else r

def kmers_from_fasta(path, k):
    out = set()
    seq = []
    with open_auto(path) as fh:
        for raw in fh:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if seq:
                    s = "".join(seq)
                    for i in range(0, len(s) - k + 1):
                        km = s[i:i+k]
                        if all(c in "ACGTacgt" for c in km):
                            out.add(canonical(km))
                    seq = []
            else:
                seq.append(line)
        if seq:
            s = "".join(seq)
            for i in range(0, len(s) - k + 1):
                km = s[i:i+k]
                if all(c in "ACGTacgt" for c in km):
                    out.add(canonical(km))
    return out

def main():
    if len(sys.argv) != 4:
        print("usage: verify_pair.py <orig_fasta_or_gz> <dump_fasta> <k>", file=sys.stderr)
        return 2
    orig, dump, kraw = sys.argv[1], sys.argv[2], sys.argv[3]
    k = int(kraw)
    a = kmers_from_fasta(orig, k)
    b = kmers_from_fasta(dump, k)
    only_a = len(a - b)
    only_b = len(b - a)
    status = "OK" if only_a == 0 and only_b == 0 else "FAIL"
    print(f"{status}\torig={len(a)}\tdump={len(b)}\tonly_orig={only_a}\tonly_dump={only_b}")
    return 0 if status == "OK" else 1

if __name__ == "__main__":
    raise SystemExit(main())
PY
chmod +x "$VERIFY_PY"

WORKER_SH="$OUT_BASE/worker.sh"
cat >"$WORKER_SH" <<'SH'
#!/usr/bin/env bash
set -euo pipefail
line="$1"
IFS=$'\t' read -r idx src qfile outdir status_file log_file <<<"$line"

printf '%s\n' "$src" >"$qfile"

"$BIN" decompress \
  -c "$ARCHIVE_DIR/" \
  -o "$outdir/" \
  -Q "$qfile" \
  -t "$DECOMP_THREADS" \
  -r "$MEMORY_GB" \
  >"$log_file" 2>&1

base="$(basename "$src")"
stem="${base%.gz}"
dump="$outdir/Dump_${stem}.fa"
if [[ ! -f "$dump" ]]; then
  printf "FAIL\tidx=%s\treason=missing_dump\tpath=%s\n" "$idx" "$dump" >"$status_file"
  exit 1
fi

if verify_out="$(python3 "$VERIFY_PY" "$src" "$dump" "$K" 2>&1)"; then
  printf "OK\tidx=%s\t%s\n" "$idx" "$verify_out" >"$status_file"
else
  printf "FAIL\tidx=%s\t%s\n" "$idx" "$verify_out" >"$status_file"
  exit 1
fi

if [[ "$CLEAN_RUN_ARTIFACTS" == "1" ]]; then
  rm -f "$qfile"
  rm -rf "$outdir"
fi
SH
chmod +x "$WORKER_SH"

echo "[info] output base: $OUT_BASE"
echo "[info] salmonella dir: $SALMON_DIR"
echo "[info] sizes: ${SIZES[*]}"
echo "[info] threads=$THREADS memory_gb=$MEMORY_GB k=$K"
echo "[info] check_jobs=$CHECK_JOBS decomp_threads=$DECOMP_THREADS clean=$CLEAN_RUN_ARTIFACTS"

for n in "${SIZES[@]}"; do
  run_dir="$OUT_BASE/n${n}"
  archive_dir="$run_dir/archive"
  tmp_dir="$run_dir/tmp"
  logs_dir="$run_dir/logs"
  q_dir="$run_dir/queries"
  d_dir="$run_dir/decompressed"
  s_dir="$run_dir/status"
  mkdir -p "$archive_dir" "$tmp_dir" "$logs_dir" "$q_dir" "$d_dir" "$s_dir"

  fof="$run_dir/fof_${n}.txt"
  printf '%s\n' "${INPUTS[@]:0:$n}" >"$fof"

  echo "[n=$n] compress start"
  "$BIN" compress \
    -i "$fof" \
    -o "$archive_dir/" \
    -d "$tmp_dir" \
    -t "$THREADS" \
    -r "$MEMORY_GB" \
    -k "$K" \
    >"$logs_dir/compress.log" 2>&1
  echo "[n=$n] compress done"

  tasks="$run_dir/tasks.tsv"
  : >"$tasks"
  for ((i = 1; i <= n; i++)); do
    src="${INPUTS[$((i - 1))]}"
    qfile="$q_dir/q_$(printf '%02d' "$i").txt"
    outdir="$d_dir/q_$(printf '%02d' "$i")"
    status_file="$s_dir/q_$(printf '%02d' "$i").status"
    log_file="$logs_dir/decomp_q_$(printf '%02d' "$i").log"
    printf "%s\t%s\t%s\t%s\t%s\t%s\n" "$i" "$src" "$qfile" "$outdir" "$status_file" "$log_file" >>"$tasks"
  done

  echo "[n=$n] targeted decompress + k-mer checks in parallel"
  export BIN ARCHIVE_DIR="$archive_dir" DECOMP_THREADS MEMORY_GB VERIFY_PY K CLEAN_RUN_ARTIFACTS
  if ! xargs -P "$CHECK_JOBS" -d $'\n' -I{} "$WORKER_SH" "{}" <"$tasks"; then
    ok_count=0
    for status_path in "$s_dir"/*.status; do
      [[ -f "$status_path" ]] || continue
      if grep -q '^OK' "$status_path"; then
        ok_count=$((ok_count + 1))
      fi
    done
    printf "%s\tFAIL\t%s\t%s\t%s\n" "$n" "$ok_count" "$n" "$run_dir" >>"$SUMMARY"
    echo "[n=$n] FAILED ($ok_count/$n checks passed). See $run_dir"
    exit 1
  fi

  ok_count=0
  for status_path in "$s_dir"/*.status; do
    [[ -f "$status_path" ]] || continue
    if grep -q '^OK' "$status_path"; then
      ok_count=$((ok_count + 1))
    fi
  done
  if [[ "$ok_count" -ne "$n" ]]; then
    printf "%s\tFAIL\t%s\t%s\t%s\n" "$n" "$ok_count" "$n" "$run_dir" >>"$SUMMARY"
    echo "[n=$n] FAILED ($ok_count/$n checks passed)."
    exit 1
  fi

  printf "%s\tOK\t%s\t%s\t%s\n" "$n" "$ok_count" "$n" "$run_dir" >>"$SUMMARY"
  echo "[n=$n] OK ($ok_count/$n checks passed)"

  if [[ "$CLEAN_RUN_ARTIFACTS" == "1" ]]; then
    rm -rf "$archive_dir" "$tmp_dir"
  fi
done

echo ""
echo "All Salmonella checks passed for sizes: ${SIZES[*]}."
echo "Summary: $SUMMARY"
