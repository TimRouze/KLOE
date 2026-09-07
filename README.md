# kloe

## Compilation
Download source code from Github

```sh
git clone  https://github.com/TimRouze/KLOE
cd KLOE
```

## Usage example
This projects uses [simd-minimizers](https://github.com/rust-seq/simd-minimizers?tab=readme-ov-file) which requires AVX2 or NEON instruction sets, which, on x64, requires either target-cpu=native or target-cpu=x86-64-v3. See [this README](https://github.com/ragnargrootkoerkamp/ensure_simd) for details.

```sh
cargo build -r
./target/release/kloe compress -i path/to/file/of/file -o output/path -d temporary/folder -t 12 -r 32
```
Builds are configured to always use optimized profiles and native CPU flags from `.cargo/config.toml`.

This creates a KLOE v4 archive containing the k-mer content of every file in the
input file-of-files. New archives use a dataset-major membership index; they do
not store the inverse CID-to-dataset relation a second time.

### KLOE v4 archive layout

- `manifest.kloe`: format version, k and minimizer sizes, dataset/CID counts,
  and abundance configuration.
- `filenames_id.txt`: dataset names in dataset-ID order.
- `filenames.idx`: the combined dataset directory. It contains the `N+1`
  dataset-posting boundaries, ordered name offsets/lengths, and a collision-safe
  open-addressed name hash table at a maximum 75% load. Hash slots contain only
  a fingerprint and dataset ID; collisions are verified against the complete
  name.
- `dataset_to_cid.bin`: one independently addressable zstd frame per dataset.
  Sparse datasets use sorted CID deltas encoded as unsigned varints. Dense
  datasets compare that representation with a CID bitmap and with a bitmap XOR
  a single archive-wide majority-membership bitmap; KLOE stores the smallest.
  The bitmap codecs are considered only when scanning the CID universe costs at
  most eight times the number of returned CIDs. In abundance mode, CID deltas
  and abundance deviations use separate zstd columns inside the dataset frame.
- `abundance_base.bin` and `abundance_base.idx` (abundance mode only): one
  block-compressed, directly indexed base abundance per CID. Dataset postings
  store deviations from this base, avoiding duplicate absolute codes.
- `positions_kloe.bin`: delta-varint encoded boundaries into the logical tig
  stream.
- `positions_kloe.idx`: an absolute checkpoint every 256 boundaries. At most
  255 small deltas are decoded to locate an arbitrary CID.
- `bucket_sizes.txt`: variable-length tig sizes, grouped into independently
  zstd-compressed metadata blocks.
- `bucket_sizes.idx`: direct routing from a CID to its compressed size block.
- `tigs_kloe.fa`: 2-bit DNA in independently zstd-compressed 8 MiB logical
  blocks.
- `tigs_kloe.idx`: fixed-size descriptors mapping logical tig offsets directly
  to compressed blocks.

The small `.idx` files deliberately remain uncompressed: they are fixed-width
random-access routing structures. The bulk data stays delta/varint encoded and
zstd compressed. This avoids the previous permanent duplication of the full
dataset/CID relation while retaining direct access.

For decompression, run:
```sh
# WHOLE ARCHIVE DECOMPRESSION
./target/release/kloe decompress -o Output/path/for/decompressed/data -c path/to/compressed/archive/directory -r 16
# TARGETED DECOMPRESSION
./target/release/kloe decompress -o Output/path/for/decompressed/data -c path/to/compressed/archive/directory -Q TARGET/FILES/LIST -r 16
# UNION OF THE COLORS LISTED IN TARGET/FILES/LIST
./target/release/kloe decompress -o Output/path/for/decompressed/data -c path/to/compressed/archive/directory -Q TARGET/FILES/LIST --color-set-operation union -r 16
# INTERSECTION OF THE COLORS LISTED IN TARGET/FILES/LIST, AS UNITIGS
./target/release/kloe decompress -o Output/path/for/decompressed/data -c path/to/compressed/archive/directory -Q TARGET/FILES/LIST --color-set-operation intersection --unitig -r 16
```

To merge two existing archives, run:
```sh
./target/release/kloe merge -c path/to/archive_A --compressed-dir-b path/to/archive_B -o path/to/merged/archive -r 16
```

### Compression parameters
When running kloe in compression mode, add "compress" before any other parameter.

#### Threads used -t
This parameter defines the number of threads used by kloe to construct the archive.
The default value is 1 thread.

#### out-dir -o
The output directory where the compressed archive should be written.
Default is current directory

#### temp_dir -d
A temporary directory to write temporary files used during compression workflow.
Default is current directory

#### k_size -k
K size. Default is 31.

#### minimizer_size -m
Minimizer size. Default is 7.

#### unitigs
If this flag is set, the archive will contain monochromatic unitigs.

#### matchtigs
If this flag is set, the archive will contain monochromatic matchtigs.

#### eulertigs
If this flag is set, the archive will contain monochromatic eulertigs.

By default, if none of the above flags are set, compression uses ggcat monochromatic unitigs mode.
There can only be one flag set at once or 0. If several flags are set the tool will not run and raise an error.

#### Embedded ggcat structured-output backend (hard switch)
KLOE embeds `ggcat` in-process (no external `ggcat` executable required at runtime).

Compression is driven by a forked ggcat structured-sequence backend. Sequence and
color-run data are streamed directly into bounded KLOE binary chunks; no intermediate
graph FASTA is formatted, stored, or parsed.

The streamed data is written directly in the KLOE v4 archive layout described
above. Legacy archives remain readable.

The vendored ggcat fork is pinned to commit:
`fe6a633e64f60cd7266951d73c1def5cc023fa96`

`-r/--memory` is the total compression memory budget in GB. KLOE derives bounded
post-processing windows, graph-record chunks, group batches, and transpose blocks
from this value. Embedded ggcat receives 75% of the budget and uses
disk-backed intermediate storage; the remainder is reserved for KLOE and I/O.

Dataset-to-color indexes use a bounded external transpose into dataset-major,
CID-monotonic posting lists. Graph positions are also generated as
external-memory streams, so their bulk payloads do not accumulate in RAM as
dataset count grows.

### Targeted-query path and complexity

For each requested name, KLOE probes `filenames.idx` on disk and verifies the
stored name. It then opens only that dataset's zstd posting frame. For each CID
in the posting list, the size, logical tig range, and compressed DNA blocks are
located through their indexes; no complete archive metadata file is loaded or
scanned. Multiple requested datasets are processed in bounded batches and their
sorted postings are merged so a shared CID is decoded once per batch.

For one requested dataset `S_j`, the expected work is
`O(query-name bytes + |S_j|)` and memory is bounded independently of the archive
dataset count. The indexes can introduce only fixed format-level read
amplification: at most one 256-entry position checkpoint block, one bounded size
block, and the 8 MiB tig blocks intersecting the requested ranges. Each bitmap
codec is restricted to datasets containing at least one eighth of the CID
universe, so its full-universe scan is bounded by a constant multiple of the
returned posting count. Because every
visited CID contributes sequence to `S_j`, traversal and emitted-data work are
output-sensitive. The output file itself necessarily costs `Theta(|S_j|)` to
write.


## Archive decompression
When running kloe in decompression mode, add "decompress" before any other parameter.

#### compressed-dir -c
Input directory for decompression, the directory where the compressed kloe archive is saved.

#### Wanted files -Q
For targeted decompression, a list of files the user wants to decompress from the archive.

#### out-dir -o
Directory in which the decompressed files should be written to.

#### output compression (`--output-compression`)
Decompressed FASTA is zstd-compressed by default. Select one of:

- `zstd` (default), producing `Dump_*.fa.zst`
- `gz`, producing `Dump_*.fa.gz`
- `xz`, producing `Dump_*.fa.xz`
- `fasta` (aliases: `uncompressed`, `fa`), producing uncompressed `Dump_*.fa`

#### color-set union and intersection (`--color-set-operation`)
Use `-Q/--wanted-files` to provide the archived dataset names (one per line),
then select `--color-set-operation union` or `--color-set-operation intersection`.
KLOE emits each qualifying color-class sequence once in `Dump_union.*` or
`Dump_intersection.*`. Every requested color must exist in the archive. This
operation is available for indexed KLOE archives and cannot be combined with
`--abundance`.

#### ggcat rebuild after decompression (`--ggcat-rebuild`)
After writing `Dump_*` files, run embedded ggcat on those dumps and produce a
rebuilt compacted output in the output directory. Supplying `--unitig`,
`--matchtig`, or `--eulertig` automatically enables this rebuild; those flags
therefore describe the actual decompressed output rather than being ignored.
`--ggcat-rebuild` alone rebuilds simplitigs (implemented via ggcat Pathtigs).
The rebuilt filename is `rebuilt_unitigs.*`, `rebuilt_matchtigs.*`,
`rebuilt_eulertigs.*`, or `rebuilt_simplitigs.*`, with the suffix selected by
`--output-compression`.

`-r/--memory` controls the rebuild budget (GB) passed to embedded ggcat.

## Archive merge
When running kloe in merge mode, add "merge" before any other parameter.

#### compressed-dir -c
Input directory for archive A.

#### compressed-dir-b
Input directory for archive B.

#### out-dir -o
Output directory where the merged archive is written.

#### memory -r
The default merge rebuilds the combined colored graph with embedded GGCAT. It removes
sequence redundancy shared between the two input archives. `-r/--memory` is the
global memory budget for this operation.

#### structural-merge
Pass `--structural-merge` only when a container-level concatenation is desired.
This bounded-memory path preserves the input SPSS records but does not remove
redundancy shared between the archives. It cannot be combined with a tig output mode.

#### output tig mode flags
Merge honors the same output mode flags as compression:
- default: simplitigs produced by GGCAT recompaction
- `--unitig`
- `--matchtig`
- `--eulertig`

Output tig modes apply to the default GGCAT merge.
