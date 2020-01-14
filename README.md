# exon_aligneR

Given the exon sequences of two candidate orthologues create an alignment of
the genes as a set of exons in order to identify orthologous exons. The
algorithm considers exons as residues to be aligned using a modified
Needleman-Wunsch algorithm.

This also includes a simpler Needleman-Wunsch implementation that takes an
arbitrary substitution matrix. This can be used to align introns by inserting
a special character at exon boundaries and assigning a high mismatch
penalty. This is probably a better means of identifying orthologous introns,
though exactly what the penalties should be is not completely clear.

The alignments will be implemented in C.

### Compilation

```sh
R CMD SHLIB exon_aligneR.c needleman_wunsch.c

## or for more speed ?
MAKEFLAGS="CFLAGS=-O3" R CMD SHLIB exon_aligneR.c needleman_wunsch.c
```

### Loading of library
From R:
```R
dyn.load(".../exon_aligneR/src/exon_aligneR.so")
### where ... is the path to exon_alignR
```


