# kloe

## Compilation
Download source code from Github

```sh
git clone  https://github.com/TimRouze/KLOE
cd KLOE
```

## Usage example
```sh
I="PATH/TO/FOF" K=31 cargo build -r
./target/release/kloe compress -t 12
```
This will create a compressed KLOE archive with the k-mer content of every files in the input file of file.
The archive is composed of 3 files:
- multicolor.kloe
- omnicolor.kloe
- multicolor_bucket_files.txt.zst

For decompression, Kloe currently needs to have a specific fof in the same directory as the other 3. To create this file, simply use this command:
```sh
awk '{print $0":"NR-1}' PATH/TO/FOF > filename_to_color.txt
```
Then simply run:
```sh
# WHOLE ARCHIVE DECOMPRESSION
./target/release/kloe decompress --omnicolor-file omnicolor.kloe --multicolor-file multicolor.kloe
# TARGETED DECOMPRESSION
./target/release/kloe decompress --omnicolor-file omnicolor.kloe --multicolor-file multicolor.kloe --wanted-files TARGET/FILES/LIST
```
### Compile time parameters
There are two parameters that should be given during build.

#### Input File of File I=
The input file of fasta files.

#### K-mer size -k
length of the k-mers used.

### Compression parameters
When running kloe in compression mode, add "compress" before any other parameter.

#### Threads used -t
This parameter defines the number of threads used by kloe to construct the archive.
The default value is 1 thread.

#### out-dir -o
The output directory where the compressed archive should be written.
Default is current directory

## Archive decompression
When running kloe in compression mode, add "decompress" before any other parameter.

#### Omnicolored file --omnicolor-file
The file containing omnicolored monochromatigs

#### Multicolored file --multicolor-file
The file containing multicolored monochromatigs

#### Input directory -i
The directory where kloe should fetch the interfacing file.

#### Wanted files -Q
For targeted decompression, a list of files the user wants to decompress from the archive.

