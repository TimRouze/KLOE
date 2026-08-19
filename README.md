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

This will create a compressed KLOE archive with the k-mer content of every files in the input file of file.
The archive is composed of 5 files:
- tigs_kloe.fa
- id_to_color_id.txt.zst
- positions_kloe.txt.zst
- bucket_sizes.txt
- filenames_id.txt

For decompression, run:
```sh
# WHOLE ARCHIVE DECOMPRESSION
./target/release/kloe decompress -o Output/path/for/decompressed/data -c path/to/compressed/archive/directory
# TARGETED DECOMPRESSION
./target/release/kloe decompress -o Output/path/for/decompressed/data -c path/to/compressed/archive/directory -Q TARGET/FILES/LIST
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

The on-disk KLOE archive format is unchanged (`tigs_kloe.fa`, `bucket_sizes.txt`, `positions_kloe.bin`, `id_to_color_id.txt.zst`, `filenames_id.txt`).

The vendored ggcat fork is pinned to commit:
`fe6a633e64f60cd7266951d73c1def5cc023fa96`

`-r/--memory` is the total compression memory budget in GB. KLOE derives bounded
post-processing windows, graph-record chunks, group batches, and transpose blocks
from this value. Embedded ggcat receives 75% of the budget and uses
disk-backed intermediate storage; the remainder is reserved for KLOE and I/O.

Dataset-to-color indexes use a linear, CID-monotonic external transpose rather than
comparison sorting. Graph positions are also generated as external-memory streams,
so their bulk payloads do not accumulate in RAM as dataset count grows.


## Archive decompression
When running kloe in compression mode, add "decompress" before any other parameter.

#### compressed-dir -c
Input directory for decompression, the directory where the compressed kloe archive is saved.

#### Wanted files -Q
For targeted decompression, a list of files the user wants to decompress from the archive.

#### out-dir -o
Directory in which the decompressed files should be written to.

#### ggcat rebuild after decompression (`--ggcat-rebuild`)
After writing `Dump_*.fa` files, run embedded ggcat on those dumps and produce a rebuilt compacted output in the output directory:
- `rebuilt_unitigs.fa` when `--unitig` is set
- `rebuilt_matchtigs.fa` when `--matchtig` is set
- `rebuilt_eulertigs.fa` when `--eulertig` is set
- `rebuilt_simplitigs.fa` by default (implemented via ggcat Pathtigs)

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
The default merge is structural: packed SPSS bytes and compact metadata are streamed
from archive A followed by archive B. It does not construct a graph, decode tigs, or
materialize CID memberships. Its memory use is bounded by I/O and one compressed
dataset-metadata record rather than by the number of k-mers or color sets. `-r` is
retained for CLI consistency but is only used by merge recompaction.

#### recompact-merge
Pass `--recompact-merge` to rebuild the combined colored graph with embedded GGCAT.
This can remove sequence redundancy shared by the two input archives, but costs
substantially more time and memory. The `-r/--memory` value is the global budget for
this path. The input archives are still streamed without per-dataset FASTA dumps.

#### output tig mode flags
Merge honors the same output mode flags as compression:
- default structural merge: preserve each input archive's SPSS
- `--unitig`
- `--matchtig`
- `--eulertig`

Selecting an output tig mode enables GGCAT recompaction because a structural merge
cannot change the input SPSS representation.
