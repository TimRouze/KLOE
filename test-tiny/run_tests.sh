#!/usr/bin/env bash
# Round-trip test for KLOE: build, compress, decompress (full + targeted), verify k-mers.
# Usage: ./test-tiny/run_tests.sh [--target-dir DIR]
#
# Runs from the repo root. Exits non-zero on first failure.

set -euo pipefail

REPO_DIR="$(cd "$(dirname "$0")/.." && pwd)"
TARGET_DIR=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        --target-dir) TARGET_DIR="$2"; shift 2 ;;
        *) echo "Unknown option: $1"; exit 1 ;;
    esac
done

if [[ -z "$TARGET_DIR" ]]; then
    TARGET_DIR=$(mktemp -d /tmp/kloe_test_target.XXXXXX)
    CLEANUP_TARGET=1
else
    CLEANUP_TARGET=0
fi

WORK_DIR=$(mktemp -d /tmp/kloe_test_work.XXXXXX)
trap 'rm -rf "$WORK_DIR"; [[ $CLEANUP_TARGET -eq 1 ]] && rm -rf "$TARGET_DIR"' EXIT

BINARY="$TARGET_DIR/release/kloe"
COMPRESS_OUT="$WORK_DIR/compressed"
DECOMP_ALL="$WORK_DIR/decomp_all"
DECOMP_WANTED="$WORK_DIR/decomp_wanted"

echo "=== Building KLOE ==="
RUSTFLAGS="-C target-cpu=native" cargo build --release --manifest-path "$REPO_DIR/Cargo.toml" --target-dir "$TARGET_DIR" 2>&1
echo

echo "=== Compress (test-tiny, 3 files) ==="
mkdir -p "$COMPRESS_OUT"
"$BINARY" compress -i "$REPO_DIR/test-tiny/fof.txt" -o "$COMPRESS_OUT/" -t 1 2>&1
echo

echo "=== Decompress all ==="
mkdir -p "$DECOMP_ALL"
"$BINARY" decompress -c "$COMPRESS_OUT/" -o "$DECOMP_ALL/" 2>&1
# Check all 3 output files exist
for f in sample1.fa sample2.fa sample3.fa; do
    if [[ ! -f "$DECOMP_ALL/Dump_$f" ]]; then
        echo "FAIL: expected $DECOMP_ALL/Dump_$f not found"
        exit 1
    fi
done
echo

echo "=== Verify k-mers (full decompression) ==="
python3 "$REPO_DIR/test-tiny/verify_kmers.py" "$REPO_DIR/test-tiny" "$DECOMP_ALL"
echo

echo "=== Decompress targeted (sample1 only) ==="
WANTED="$WORK_DIR/wanted.txt"
echo "test-tiny/sample1.fa" > "$WANTED"
mkdir -p "$DECOMP_WANTED"
"$BINARY" decompress -c "$COMPRESS_OUT/" -o "$DECOMP_WANTED/" -Q "$WANTED" 2>&1
# Check only sample1 was decompressed
if [[ ! -f "$DECOMP_WANTED/Dump_sample1.fa" ]]; then
    echo "FAIL: expected $DECOMP_WANTED/Dump_sample1.fa not found"
    exit 1
fi
NFILES=$(find "$DECOMP_WANTED" -name 'Dump_*.fa' | wc -l)
if [[ "$NFILES" -ne 1 ]]; then
    echo "FAIL: targeted decompression produced $NFILES files, expected 1"
    exit 1
fi
echo

echo "=== Verify k-mers (targeted decompression) ==="
python3 "$REPO_DIR/test-tiny/verify_kmers.py" "$REPO_DIR/test-tiny" "$DECOMP_WANTED" sample1.fa
echo

echo "=== ALL TESTS PASSED ==="
