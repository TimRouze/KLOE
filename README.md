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
With this version, it is necessary to build the project again each time you change the input file of file and the k-mer size typing this command (replacing the path by the one of your file of file): 
```sh
I="PATH/TO/FOF" K=31 cargo build -r
```
Then, the executable will be in the following folder:
```sh
target/release/
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
## Compile time parameters
There are two parameters that should be given during build.

#### Input File of File I=
The input file of fasta files.

#### K-mer size -k
length of the k-mers used.

## Compression parameters
When running kloe in compression mode, add "compress" before any other parameter.

#### Threads used -t
This parameter defines the number of threads used by kloe to construct the archive.
The default value is 1 thread.

#### out-dir -o
The output directory where the compressed archive should be written.
Default is current directory

## Decompression parameters
When running kloe in compression mode, add "decompress" before any other parameter.

#### Omnicolored file --omnicolor-file
The file containing omnicolored monochromatigs

#### Multicolored file --multicolor-file
The file containing multicolored monochromatigs

#### Input directory -i
The directory where kloe should fetch the interfacing file.

#### Wanted files -Q
For targeted decompression, a list of files the user wants to decompress from the archive.

