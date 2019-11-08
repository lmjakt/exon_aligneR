dyn.load("src/exon_aligneR.so")
source("functions.R")

## then we read in a file of sequences
exons.l <- readLines("exons_PTHR13342")
id.l <- grep(">", exons.l)
s.e <- c(id.l - 1, length(exons.l))[-1]

exons <- lapply( 1:length(id.l), function(i){
    exons.l[(id.l[i]+1):s.e[i]]
})

names(exons) <- sapply( exons.l[id.l], function(x){ strsplit(x, split="\t")[[1]][1] })


transcripts <- lapply(exons, function(x){
    list('l'=unname(sapply(x, nchar)), 's'=paste(x, collapse=''))
})
    

### to compare two set of exons:

.Call( "align_exons", transcripts[[1]]$s, transcripts[[2]]$s, transcripts[[1]]$l, transcripts[[2]]$l, c(4.0, -4.0, -8.0, -0.5), 0.5)
.Call( "align_exons", transcripts[[1]]$s, transcripts[[4]]$s, transcripts[[1]]$l, transcripts[[4]]$l, c(4.0, -4.0, -8.0, -0.5), 0.5)

.Call( "align_exons", exons[[1]], exons[[2]], c(4.0, -4.0, -8.0, -1.0), 0.5)
.Call( "align_exons", exons[[1]], exons[[2]], c(4.0, -4.0, -8.0, -0.5), 0.5)
.Call( "align_exons", exons[[2]], exons[[1]], c(4.0, -4.0, -8.0, -1.0), 0.5)
.Call( "align_exons", exons[[1]], exons[[1]], c(4.0, -4.0, -8.0, -1.0), 0.5)

## this is aligning 4 exons to 5
.Call( "align_exons", exons[[1]], exons[[4]], c(4.0, -4.0, -12.0, -0.5), 0.5)
.Call( "align_exons", exons[[4]], exons[[1]], c(4.0, -4.0, -12.0, -0.5), 0.5)
## this gives us:
## [[1]]
##      [,1] [,2] [,3] [,4]
## [1,] -126   -9  -58 -345
## [2,]  -47  -10  -41 -301
## [3,]  -89  192  -23 -282
## [4,]  -22  -27  272 -168
## [5,]  -91    8  -28 -232

## [[2]]
##      [,1] [,2] [,3] [,4] [,5]
## [1,]    0    0    0    0    0
## [2,]    0  -70   -9  -58  -70
## [3,]    0  -47  -80  -50 -184
## [4,]    0  -89  145  -91 -308
## [5,]    0  -22  -91  417 -259
## [6,]    0  -91  -14  295  185

## [[3]]
##      [,1] [,2] [,3] [,4] [,5]
## [1,]    0    2    2    2    2
## [2,]    1    1    3    3    1
## [3,]    1    3    3    3    1
## [4,]    1    3    3    2    1
## [5,]    1    3    1    3    3
## [6,]    1    3    3    1    3



##alignments <- lapply( 1:length(exons), function(i){ .Call("align_exons", exons[[1]], exons[[i]], c(4.0, -4.0, -8.0, -0.5), 1.0 ) })
alignments.2 <- lapply( 1:length(transcripts), function(i){ align.exons( transcripts[[1]], transcripts[[i]], c(4.0, -4.0, -8.0, -0.5), 1.0 ) })
alignments.3 <- lapply( 1:length(transcripts), function(i){ align.exons( transcripts[[1]], transcripts[[i]], c(4.0, -4.0, -8.0, -0.5), 1.0, TRUE ) })

alignments.4 <- vector(mode='list', length = (17^2 -17 )/2 )

system.time(
    for(i in 1:100)
        alignments <- lapply( 1:length(exons), function(i){ align.exons( transcripts[[1]], transcripts[[i]], c(4.0, -4.0, -8.0, -0.5), 1.0 ) })
)

a2 <- align.exons( transcripts[[4]], transcripts[[9]], c(4.0, -4.0, -8.0, -0.5), 0.1 )

plot.exon.alignments( alignments[[4]] )


