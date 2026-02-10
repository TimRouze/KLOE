# kloe

## Compilation
Download source code from Github

```sh
git clone  https://github.com/TimRouze/KLOE
cd KLOE
```

## Usage example
This projects uses [simd-minimizers](https://github.com/rust-seq/simd-minimizers?tab=readme-ov-file) which requires AVX2 or NEON instruction sets, which, on x64, requires either target-cpu=native or target-cpu=x86-64-v3. See [this README](https://github.com/ragnargrootkoerkamp/ensure_simd) for details.

This projects uses [simd-minimizers](https://github.com/rust-seq/simd-minimizers?tab=readme-ov-file) which requires AVX2 or NEON instruction sets, which, on x64, requires either target-cpu=native or target-cpu=x86-64-v3. See [this README](https://github.com/ragnargrootkoerkamp/ensure_simd) for details.

```sh
RUSTFLAGS="-C target-cpu=native" cargo build -r
./target/release/kloe compress -i path/to/file/of/file -o output/path -d temporary/folder -t 12
```
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
./target/release/kloe decompress -o Output/path/for/decompressed/data -c path/to/compressed/archive/directory
# TARGETED DECOMPRESSION
./target/release/kloe decompress -o Output/path/for/decompressed/data -c path/to/compressed/archive/directory -Q TARGET/FILES/LIST
./target/release/kloe decompress -o Output/path/for/decompressed/data -c path/to/compressed/archive/directory -Q TARGET/FILES/LIST
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

By default, if none of the above flags are set, the archive will contain monochromatic simplitigs.
There can only be one flag set at once or 0. If several flags are set the tool will not run and raise an error.


## Archive decompression
When running kloe in compression mode, add "decompress" before any other parameter.

#### compressed-dir -c
Input directory for decompression, the directory where the compressed kloe archive is saved.
#### compressed-dir -c
Input directory for decompression, the directory where the compressed kloe archive is saved.

#### Wanted files -Q
For targeted decompression, a list of files the user wants to decompress from the archive.

#### out-dir -o
Directory in which the decompressed files should be written to.
#### out-dir -o
Directory in which the decompressed files should be written to.
