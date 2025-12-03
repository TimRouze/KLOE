# Experiments for RECOMB-SEQ 2025 submission

## File of files

For experiments, we used 2 datasets.
- A collection of 1024 human unitigs downloaded from [Logan](https://github.com/IndexThePlanet/Logan). [Human unitigs accession file](fof_logan_unitigs.txt)
- A collection of 7318 Salmonella genomes downloaded from the NCBI. [Salmonellas accession file](fof_7318_salmo.txt)

## Command lines
### KLOE
```sh
I="PATH/TO/FOF" K=31 cargo build -r --target-dir "PATH/TO/BUILD/DIRECTORY"
\time ./PATH/TO/BUILD/DIRECTORY/kloe compress -o "OUTPATH" -t 20 2> "PATH/TO/TIME/FILE"
```

### GGCAT
```sh
\time ggcat build -l "PATH/TO/FOF" -k 31 -o "OUTPATH" -c -t PATH/TO/TEMP/DIRECTORY -m 500 -s 1 2> "PATH/TO/TIME/FILE"
\time zstd "ggcat_output" -f -o "ggcat_ouput.zst" 2> "PATH/TO/TIME/FILE"
```

### XZ -1
```sh
XZ_OPT=-1 \time tar --create --xz --file="OUTPATH" --files-from=<(echo "PATH/TO/FOF") 2> "PATH/TO/TIME/FILE"
```

### XZ -3
```sh
\time tar --create --xz --file="OUTPATH" --files-from=<(echo "PATH/TO/FOF") 2> "PATH/TO/TIME/FILE"
```

### ZSTD -1
```sh
\time tar -c --files-from=<(echo "PATH/TO/FOF") | zstd -1 -o "OUTPATH" 2> "PATH/TO/TIME/FILE"
```

### ZSTD -6
```sh
\time tar -c --files-from=<(echo "PATH/TO/FOF") | zstd -o "OUTPATH" 2> "PATH/TO/TIME/FILE"
```
