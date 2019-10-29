dyn.load("src/exon_aligneR.so")

## then we read in a file of sequences
exons.l <- readLines("exons_PTHR13342")
id.l <- grep(">", exons.l)
s.e <- c(id.l - 1, length(exons.l))[-1]

exons <- lapply( 1:length(id.l), function(i){
    exons.l[(id.l[i]+1):s.e[i]]
})

names(exons) <- sapply( exons.l[id.l], function(x){ strsplit(x, split="\t")[[1]][1] })

### to compare two set of exons:

.Call( "align_exons", exons[[1]], exons[[2]], c(4.0, -4.0, -8.0, -1.0), 2.0)
.Call( "align_exons", exons[[1]], exons[[2]], c(4.0, -4.0, -8.0, -0.5), 2.0)

alignments <- lapply( 1:length(exons), function(i){ .Call( "align_exons", exons[[1]], exons[[i]], c(4.0, -4.0, -8.0, -0.5), 2.0) })
    


