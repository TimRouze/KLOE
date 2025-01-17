USAGE:
```sh
#COMPRESS
I=../../fof_6_salmonellas.txt K=31 cargo r -r -- compress -t 12

#DECOMPRESS W/ QUERY
I=../../fof_6_salmonellas.txt K=31 cargo r -r -- decompress --omnicolor-file omnicolor.kloe --multicolor-file multicolor.kloe  -Q wanted_salmo.txt -t 12

#DECOMPRESS WHOLE ARCHIVE
I=../../fof_6_salmonellas.txt K=31 cargo r -r -- decompress --omnicolor-file omnicolor.kloe --multicolor-file multicolor.kloe-t 12
```

## TODO
- Change bitvec par vec<u8> (liste d'entier)
- taille simplitigs -> taille -K + 1
- FIX K > 31
- DELETE UTILITY FILES AFTER DECOMPRESSION?
- STORE ABUNDANCES
