#!/usr/bin/env python3
"""Verify canonical 31-mers match between original and decompressed FASTA files.

Usage:
    python verify_kmers.py <original_dir> <decompressed_dir> [sample1.fa sample2.fa ...]

If no sample names are given, all *.fa files in original_dir are checked.
Decompressed files are expected to be named Dump_<samplename> in decompressed_dir.
"""

import glob
import os
import sys

K = 31

COMPLEMENT = str.maketrans("ACGTacgt", "TGCAtgca")


def reverse_complement(seq):
    return seq.translate(COMPLEMENT)[::-1]


def canonical(kmer):
    rc = reverse_complement(kmer)
    return min(kmer.upper(), rc.upper())


def extract_canonical_kmers(fasta_path):
    """Extract all canonical K-mers from a FASTA file."""
    kmers = set()
    seq = []
    with open(fasta_path) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if seq:
                    full_seq = "".join(seq)
                    for i in range(len(full_seq) - K + 1):
                        kmer = full_seq[i : i + K]
                        if all(c in "ACGTacgt" for c in kmer):
                            kmers.add(canonical(kmer))
                    seq = []
            else:
                seq.append(line)
        if seq:
            full_seq = "".join(seq)
            for i in range(len(full_seq) - K + 1):
                kmer = full_seq[i : i + K]
                if all(c in "ACGTacgt" for c in kmer):
                    kmers.add(canonical(kmer))
    return kmers


def main():
    if len(sys.argv) < 3:
        print(__doc__.strip())
        sys.exit(1)

    orig_dir = sys.argv[1]
    decomp_dir = sys.argv[2]
    samples = sys.argv[3:]

    if not samples:
        samples = sorted(
            os.path.basename(p) for p in glob.glob(os.path.join(orig_dir, "*.fa"))
        )

    all_ok = True
    for sample in samples:
        orig_path = os.path.join(orig_dir, sample)
        decomp_name = "Dump_" + sample
        decomp_path = os.path.join(decomp_dir, decomp_name)

        if not os.path.exists(orig_path):
            print(f"SKIP {sample}: original not found at {orig_path}")
            continue
        if not os.path.exists(decomp_path):
            print(f"SKIP {sample}: decompressed not found at {decomp_path}")
            continue

        orig_kmers = extract_canonical_kmers(orig_path)
        decomp_kmers = extract_canonical_kmers(decomp_path)

        match = orig_kmers == decomp_kmers

        print(f"--- {sample} ---")
        print(f"  Original:     {len(orig_kmers)} canonical {K}-mers")
        print(f"  Decompressed: {len(decomp_kmers)} canonical {K}-mers")

        if match:
            print(f"  MATCH")
        else:
            all_ok = False
            only_orig = orig_kmers - decomp_kmers
            only_decomp = decomp_kmers - orig_kmers
            print(f"  MISMATCH")
            print(f"    Only in original:     {len(only_orig)}")
            print(f"    Only in decompressed: {len(only_decomp)}")
            if only_orig and len(only_orig) <= 5:
                for km in sorted(only_orig):
                    print(f"      orig-only: {km}")
            if only_decomp and len(only_decomp) <= 5:
                for km in sorted(only_decomp):
                    print(f"      decomp-only: {km}")
        print()

    if all_ok:
        print("ALL PAIRS MATCH -- verification passed.")
    else:
        print("MISMATCH DETECTED -- verification FAILED.")
        sys.exit(1)


if __name__ == "__main__":
    main()
