# Examples

If you want to try to compress data with Kloe, first download the two human unitigs files listed in the file of file "2_humans.txt"
Or if you want to download your own test data, create a file of file with one accession per line and then run replacing "2_humans.txt" with your own fof name:
```sh
cd Example
awk '{system("wget https://s3.amazonaws.com/logan-pub/u/"$1"/"$1".unitigs.fa.zst")}' 2_humans.txt
readlink -e *.zst > fof_example.txt
cd ..
```

## Compress data
Then, to compress data using Kloe, type in these command lines.
```sh
I="Example/" K=31 cargo build -r
./target/release/kloe compress -o "Example/" -t 20
```

## Decompress
To decompress data from a Kloe archive, type this command to decompress everything:
```sh
./target/release/kloe decompress --omnicolor-file omnicolor.kloe --multicolor-file multicolor.kloe
```
To decompress only a subpart of the archive, first create a file of file with the filenames.
```sh
head -n 1 Example/fof_example.txt > Example/wanted.txt
./target/release/kloe decompress --omnicolor-file omnicolor.kloe --multicolor-file multicolor.kloe --wanted-files Example/wanted.txt
```

