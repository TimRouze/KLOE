# kloe
## Goal
KLOE is a compressor for unitigs file such as those present on the [Logan](https://github.com/IndexThePlanet/Logan) database. It takes as input a file of file containing the unitigs filepaths to be compressed.
The main interest of this tool is to allow efficient targeted decompression (decompressing a subpart of what is present in the compressed archive). The tool is a proof of concept aiming to show the interest of such decompression and how it can be used to compress immense databases while authorising the use of the data without needing the computing ressources to handle the entirety of it. 



## Installation

If you do not have rust installed on your computer, follow this [link](https://rustup.rs/) and follow the instructions.

Then, clone this repository and go in the 'KLOE' folder created.:

```sh
git clone  https://github.com/TimRouze/KLOE
cd KLOE
```

Correctness checks for the main binary:

```sh
cargo test --bin kloe -- --test-threads=1
cargo clippy --bin kloe
```

HTTP tests use small local fixtures rather than downloading Logan datasets.

## Usage example
This projects uses [simd-minimizers](https://github.com/rust-seq/simd-minimizers?tab=readme-ov-file) which requires AVX2 or NEON instruction sets, which, on x64, requires either target-cpu=native or target-cpu=x86-64-v3. See [this README](https://github.com/ragnargrootkoerkamp/ensure_simd) for details.

```sh
cargo build --release --bin kloe
./target/release/kloe compress -i path/to/file/of/file -o output/path -d temporary/folder -t 12 -r 32
```
Builds are configured to always use optimized profiles and native CPU flags from `.cargo/config.toml`.

This creates a KLOE archive containing the k-mer content of the available files
in the input file-of-files.

Use a fresh directory for each archive, KLOE does not automatically replace or merge existing archives.

### Stream Logan unitigs from accession names

Use `--logan-accessions accessions.txt` instead of `-i/--input-list` to fetch
Logan unitigs directly from the public bucket. The list should contain one SRA run
accession per line.

For example, create an accession list and choose whether to retain abundances:

```sh
cat > accessions.txt <<'EOF'
ERR10036117
ERR10036119
EOF

# Store the k-mer sets without abundance:
./target/release/kloe compress --logan-accessions accessions.txt \
  -o archive/ -d temporary/ -t 32 -r 32 --logan-downloads 4

# Store approximate abundance (8 bits, maximum 50000):
./target/release/kloe compress --logan-accessions accessions.txt \
  -o archive-abundance/ -d temporary/ -t 32 -r 32 \
  --abundance --abundance-bits 8 --abundance-max 50000 --logan-downloads 4

# Store rounded integer abundance exactly up to the maximum (16 bits):
./target/release/kloe compress --logan-accessions accessions.txt \
  -o archive-exact/ -d temporary/ -t 32 -r 32 \
  --abundance --abundance-bits 16 --abundance-max 50000 --logan-downloads 4
```

`--logan-downloads` limits simultaneous downloads (default 4)

`--logan-retries N` sets the number of additional attempts (default 2) Use `--logan-strict` to stop on
the first invalid, duplicate or failed accession.


To restore abundance-annotated FASTA from either abundance archive, use
`--stored-tigs --abundance`; the precision is read from the archive:

```sh
./target/release/kloe decompress -c archive-abundance/ -o restored/ -r 16 \
  --stored-tigs --abundance
```

Output is Zstandard-compressed FASTA by default. Add
`--output-compression fasta` for plain FASTA.
Abundances are rounded and capped at `--abundance-max` before encoding.

### KLOE archive layout

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

## Compression parameters
When running kloe in compression mode, add "compress" before any other parameter.

#### Threads used -t
Number of threads used by kloe to construct the archive.
The default is the number of available logical CPUs.

#### Max RAM -r
Max RAM allocated for compression.

#### out-dir -o
The output directory where the compressed archive should be written.
Default is current directory

#### temp_dir -d
A temporary directory to write temporary files used during compression workflow.
By default, graph construction uses the operating system's temporary directory.

#### k_size -k
K size. Default is 31.

#### minimizer_size -m
Minimizer size. Default is 12 and must be smaller than k.

#### abundance filter during compression (`--abundance-threshold`)
`--abundance-threshold N` keeps only input FASTA records whose header abundance
is **strictly greater than N**.

```sh
# Filter local inputs without storing abundance metadata:
./target/release/kloe compress -i inputs.fof -o filtered/ -d temporary/ \
  -t 12 -r 32 --abundance-threshold 10

# Filter streamed Logan inputs and preserve abundance:
./target/release/kloe compress --logan-accessions accessions.txt \
  -o filtered-abundance/ -d temporary/ -t 12 -r 32 \
  --abundance-threshold 10 --abundance --abundance-bits 8 --logan-downloads 4
```

#### archive SPSS mode (default: unitigs)
By default, KLOE uses the unitig construction mode of its embedded GGCAT
backend. Select at most one representation flag:

- `--unitig`: unitigs (default).
- `--simplitig`: GGCAT fast simplitigs
- `--matchtig`: GGCAT greedy matchtigs.
- `--eulertig`: GGCAT fast eulertigs.

## Archive decompression
When running kloe in decompression mode, add "decompress" before any other parameter.

#### compressed-dir -c
Input directory for decompression, the directory where the compressed kloe archive is saved.

#### Wanted files -Q
For targeted decompression, a list of files the user wants to decompress from the archive.

#### out-dir -o
Directory in which the decompressed files should be written to.

#### output compression (`--output-compression`)
Decompressed FASTA is zstd-compressed by default. Compression is independent
of the selected graph representation. Select one of:

- `zstd` (default)
- `gz`
- `xz`
- `fasta`

#### color-set operations (`--color-set-operation`)
Use `-Q/--wanted-files` to provide the archived dataset names (one per line),
then select `--color-set-operation union`, `intersection`, `difference`, or
`symmetric-difference`. Difference keeps k-mers from the **first** listed dataset
that are absent from all subsequent datasets. With two datasets, symmetric
difference keeps k-mers present in exactly one; with more datasets it keeps k-mers
present in an odd number of listed datasets.


#### Abundance encoding (`--abundance`)
During compression, choose `--abundance-bits B` (2–32, default **8**) and
`--abundance-max M` (positive integer, default **50,000**; alias `--max-abundance`).

#### decompressed representation (default: unitigs)
Normal decompression rebuilds **each requested dataset independently**, producing
one output file per dataset. Select at most one of `--unitig` (default),
`--simplitig`, `--matchtig`, or `--eulertig`. These options do not combine datasets;
only an explicit color-set query produces a combined result.

## Archive merge
When running kloe in merge mode, add "merge" before any other parameter.

#### compressed-dir -c
Input directory for archive A.

#### compressed-dir-b
Input directory for archive B.

#### out-dir -o
Output directory where the merged archive is written.
It must differ from both input directories, including aliases through symlinks,
and must not already contain archive files.

#### k_size -k
Merge infers k from v4 archive manifests and requires both inputs to use the
same k. An explicit `-k` must match that recorded value. For inputs without a
manifest, provide their original k explicitly; KLOE cannot infer a legacy
archive's k from its sequences.
Graph merge currently supports k values from 1 through 32.

#### memory -r
The default merge rebuilds the combined colored graph with embedded GGCAT. It removes
sequence redundancy shared between the two input archives. `-r/--memory` is the
global memory budget for this operation.

#### structural-merge
Pass `--structural-merge` only when a container-level concatenation is desired.
This bounded-memory path preserves the input SPSS records but does not remove
redundancy shared between the archives. It cannot be combined with a tig output mode.
It supports older compact archives containing `id_to_color_id.txt.zst`; current
v4 archives require the default graph merge.

#### output tig mode flags
Merge honors the same output mode flags as compression:

- default / `--unitig`: unitigs.
- `--simplitig`: embedded GGCAT simplitigs (experimental; not yet benchmarked).
- `--matchtig`: greedy matchtigs (experimental; not yet benchmarked).
- `--eulertig`: fast eulertigs (experimental; not yet benchmarked).

Output tig modes apply to the default GGCAT merge.
