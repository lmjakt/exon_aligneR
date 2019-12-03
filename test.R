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
exon.dbs <- unname(sapply( exons.l[id.l], function(x){ strsplit(x, split="\t")[[1]][4] }))

## to set suitable attributes:
transcripts <- lapply( 1:length(exons), function(i){
    list('id'=sub("^>", "", names(exons)[i]), 'db'=exon.dbs[i], 'sp'=sub("([^_]+)_([^_]+)_.+$", "\\1 \\2", exon.dbs[i]),
         'l'=unname(sapply(exons[[i]], nchar)), 's'=paste(exons[[i]], collapse='') )
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


l <- length(transcripts)
alignments.4 <- vector(mode='list', length = (l^2 -l )/2 )
o <- 0
for(i in 1:(length(transcripts)-1)){
    for(j in (i+1):length(transcripts)){
        o <- o + 1
        alignments.4[[o]] <- c('i'=i, 'j'=j, align.exons( transcripts[[i]], transcripts[[j]], c(4.0, -4.0, -8.0, -0.5), 0.25 ))
    }
}

alignments.5 <- vector(mode='list', length = (l^2 -l )/2 )
o <- 0
for(i in 1:(length(transcripts)-1)){
    for(j in (i+1):length(transcripts)){
        o <- o + 1
        alignments.5[[o]] <- c('i'=i, 'j'=j, align.exons( transcripts[[i]], transcripts[[j]], c(4.0, -4.0, -8.0, -0.5), 0.25, FALSE, TRUE))
    }
}

alignments.6 <- vector(mode='list', length = (l^2 -l )/2 )
o <- 0
for(i in 1:(length(transcripts)-1)){
    for(j in (i+1):length(transcripts)){
        o <- o + 1
        alignments.6[[o]] <- c('i'=i, 'j'=j, align.exons( transcripts[[i]], transcripts[[j]], c(4.0, -4.0, -8.0, -0.5), 0.25, TRUE, TRUE))
    }
}

transcripts.2 <- read.transcripts("exons_PTHR10044_SF93")
l.2 <- length(transcripts.2)
alignments.7 <- vector(mode='list', length = (l.2^2 -l.2 )/2 )
o <- 0
for(i in 1:(length(transcripts.2)-1)){
    for(j in (i+1):length(transcripts.2)){
        o <- o + 1
        alignments.7[[o]] <- c('i'=i, 'j'=j, align.exons( transcripts.2[[i]], transcripts.2[[j]], c(4.0, -4.0, -8.0, -0.5), 0.25, TRUE, TRUE))
    }
}


for(i in 1:length(alignments.6)){
    plot.exon.alignments( alignments.6[[i]], g.height=6, h1=1, i.sep=1.25, w.radius=4 )
    inpt <- readline("next: ")
}

pdf("exon_alignments_example_1.pdf", width=12, height=8)
plot.exon.alignments( alignments.6[[3]], g.height=6, h1=1, i.sep=1.25, w.radius=4 )
dev.off()

for(i in 1:length(alignments.7)){
    plot.exon.alignments( alignments.7[[i]], g.height=6, h1=1, i.sep=1.25, w.radius=4 )
    inpt <- readline("next: ")
}


## let us make a lookup matrix
lom <- matrix(0, nrow=length(transcripts), ncol=length(transcripts))
for(i in 1:length(alignments.4))
    lom[ alignments.4[[i]]$i, alignments.4[[i]]$j ] = i
## the second is always more expensive.

## then let us consider how I can consider the equivalence of exon alignments
## That is for transcripts A, B, C, .. N,
## e is: which( alignments.4[[ lom[1,t2 ]]]$g.align[, 1] == 0) 
## A1 -> B1   alignments.4[[ lom[1,2 ]]$g.align[ e, 2 ]
## A1 -> C2   alignments.4[[ lom[1,3 ]]$g.align[ e, 2 ]
## A1 -> D1   alignments.4[[ lom[1,4 ]]$g.align[ e, 2 ] 
## then the following should also be true:
## e is: which( alignments.4[[ lom[2,t2 ]]]$g.align[, 1] == B1) 
## B1 -> C2   alignments.4[[ lom[2,3 ]]$g.align[ e, 2 ]
## B1 -> D1   alignments.4[[ lom[2,4 ]]$g.align[ e, 2 ]
##
## C2 -> D1
## 
## So can we simply count this doing a simple loop?

## for exon 1 of A (the first exon we can thus say something like)


## this does seem to give me reasonable tables from which I can look at the overall agreement..
## but it is limited to whatever exon positions I have in the first transcript and that is a bit
## silly. But I don't see an elegant way around this, other than 
g1.ex <- lapply(0:3, function(s){
    seed.e <- s
    g1 <- matrix(-1, nrow=length(transcripts), ncol=length(transcripts))
    for(i in 1:(length(transcripts)-1)){
        for(j in (i+1):length(transcripts)){
            o <- lom[i, j]
            e <- which( alignments.4[[o]]$g.align[,1] == seed.e )
            if(length(e) == 1)
                g1[i, j] <- alignments.4[[o]]$g.align[e,2]
        }
        seed.e <- g1[i,i+1]
    }
    g1
})

g1.ex.l <- lapply(0:3, function(s){
    seed.e <- s
    g1 <- matrix(-1, nrow=length(transcripts), ncol=length(transcripts))
    for(i in 1:(length(transcripts)-1)){
        for(j in (i+1):length(transcripts)){
            o <- lom[i, j]
            e <- which( alignments.5[[o]]$g.align[,1] == seed.e )
            if(length(e) == 1)
                g1[i, j] <- alignments.5[[o]]$g.align[e,2]
        }
        seed.e <- g1[i,i+1]
    }
    g1
})



                 
system.time(
    for(i in 1:100)
        alignments <- lapply( 1:length(exons), function(i){ align.exons( transcripts[[1]], transcripts[[i]], c(4.0, -4.0, -8.0, -0.5), 1.0 ) })
)

a2 <- align.exons( transcripts[[4]], transcripts[[9]], c(4.0, -4.0, -8.0, -0.5), 0.1 )

plot.exon.alignments( alignments[[4]] )


## we also have a new function that takes a modified sequence and something else.
## with introns marked with
exon.i <- sapply( exons, function(x){
    paste( toupper(x), collapse='I' ) })

names(exon.i) <- sub("([^_]+)_([^_]+)_.+$", "\\1 \\2", exon.dbs )

## allow from ! to `
##
al.size <- as.integer(1 + 96 - 33)
al.offset <- as.integer(33)
sub.matrix <- make.sub.matrix(al.offset=al.offset, al.size=al.size,
                              letters=c('A', 'C', 'T', 'G', 'I'), match=c(4,4,4,4, 30) )

dyn.load("src/exon_aligneR.so")

tmp <- .Call("align_seqs", exon.i[1], exon.i[14], al.offset, al.size, sub.matrix, as.integer(c(-10, -1)), TRUE, "I" )

## the following fails because I called the R_registerRoutines
tmp <- .Call("align_seqs", exon.i[1], exon.i[14], al.offset, al.size, sub.matrix, as.integer(c(-10, -1)), TRUE)

source("~/R/general_functions.R")
image( tmp[[1]], col=hsvScale(1:255, sat=1:255/255) )

cols <- c(rgb(0.9, 0.9, 0.9), rgb(1, 0, 0),
          rep(rgb(0.5, 0.5, 0.5), 4) )
names(cols) <- c('-', 'I', 'A', 'C', 'G', 'T')

plot.new()
plot.window(xlim=c(-1, nchar(tmp[[3]][1])), ylim=c(0, 5))
x2 <- 1:nchar(tmp[[3]][1])
x1 <- x2-1
al.a <- strsplit(tmp[[3]][1], '')[[1]]
al.b <- strsplit(tmp[[3]][2], '')[[1]]
rect( x1, 1, x2, 2, col=cols[al.a], border=NA )
rect( x1, 2.5, x2, 3.5, col=cols[al.b], border=NA )

## that seems to work, but I still have a memory corruption that gets hit once in a while.. 
system.time(
    for(i in 1:100){
 ##       print(paste("i: ", i))
        tmp <- .Call("align_seqs", exon.i[1], exon.i[14], al.offset, al.size, sub.matrix, as.integer(c(-10, -1)), TRUE, "I" )
    }
)

l <- length(transcripts)
alignments.4 <- vector(mode='list', length = (l^2 -l )/2 )

l <- length(exon.i)
aligns <- vector(mode='list', length=(l^2-l)/2 )
k <- 0
for(i in 1:(l-1)){
    for(j in (i+1):l){
        k <- k + 1
        aligns[[k]] <- c('i'=i, 'j'=j, align.seqs( exon.i[i], exon.i[j], al.offset, al.size, sub.matrix, as.integer(c(-10, -1)), TRUE, "I" ) )
##        aligns[[k]] <- c('i'=i, 'j'=j, .Call("align_seqs", exon.i[i], exon.i[j], al.offset, al.size, sub.matrix, as.integer(c(-10, -1)), TRUE, "I" ) )
    }
}

for(i in 1:length(aligns)){
    sp1 <- names(exon.i)[ aligns[[i]]$i ]
    sp2 <- names(exon.i)[ aligns[[i]]$j ]
    plot.new()
    plot.window(xlim=c(-1, nchar(aligns[[i]][[5]][1])), ylim=c(0, 5))
    x2 <- 1:nchar(aligns[[i]][[5]][1])
    x1 <- x2-1
    al.a <- strsplit(aligns[[i]][[5]][1], '')[[1]]
    al.b <- strsplit(aligns[[i]][[5]][2], '')[[1]]
    rect( x1, 1, x2, 2, col=cols[al.a], border=NA )
    rect( x1, 2.5, x2, 3.5, col=cols[al.b], border=NA )
    mtext(paste(sp1, "vs", sp2, aligns[[i]]$score), cex=2)
    inpt <- readline('next')
}

i <- 1
plot.new()
plot.window( xlim=c(-1, nchar(aligns[[i]]$seq[1])), ylim=c(0,5) )
draw.aligns( aligns[[i]], 3, 1, 1, cols, sp.a=names(exon.i)[ aligns[[i]]$i ], sp.b=names(exon.i)[ aligns[[i]]$j ], sim.pos=c(1,1) )

plot.new()
plot.window( xlim=c(-1, nchar(aligns[[i]]$seq[1])), ylim=c(0,5) )
draw.aligns( aligns[[i]], 3, 1, 1, cols, sp.a=names(exon.i)[ aligns[[i]]$i ], sp.b=names(exon.i)[ aligns[[i]]$j ], sim.pos=c(1,2) )

plot.new()
plot.window( xlim=c(-1, nchar(aligns[[i]]$seq[1])), ylim=c(0,5) )
draw.aligns( aligns[[i]], 3, 1, 1, cols, sp.a=names(exon.i)[ aligns[[i]]$i ], sp.b=names(exon.i)[ aligns[[i]]$j ], sim.pos=c(2,1) )

i <- which( sapply( 1:length(aligns), function(i){ aligns[[i]]$i }) == 1 )
max.x <- max( sapply( aligns[i], function(x){ nchar(x$seq[1]) }) )
plot.new()
g.sep <- 6
top <- g.sep * length(i) + 1
plot.window( xlim=c(0, max.x), ylim=c(0, top) )
top <- top - 1
for(j in 1:length(i)){
    y1 <- top - (j-1) * g.sep - 1
    y2 <- y1 - 2
    draw.aligns( aligns[[i[j]]], y1, y2, 1, cols, sp.a=names(exon.i)[ aligns[[i[j]]]$i ], sp.b=names(exon.i)[ aligns[[i[j]]]$j ], sim.pos=c(1,1) )
}

for(sp.i in 1:length(exon.i)){
    i <- which( sapply( 1:length(aligns), function(i){ aligns[[i]]$i == sp.i || aligns[[i]]$j == sp.i }) )
    max.x <- max( sapply( aligns[i], function(x){ nchar(x$seq[1]) }) )
    plot.new()
    g.sep <- 6
    top <- g.sep * length(i) + 1
    plot.window( xlim=c(0, max.x), ylim=c(0, top) )
    top <- top - 1
    text( max.x/2, top, names(exon.i[sp.i]), pos=3, cex=2 )
    for(j in 1:length(i)){
        y1 <- top - (j-1) * g.sep - 1
        y2 <- y1 - 2
        draw.aligns( aligns[[i[j]]], y1, y2, 1, cols, sp.a=names(exon.i)[ aligns[[i[j]]]$i ], sp.b=names(exon.i)[ aligns[[i[j]]]$j ], sim.pos=c(1,1) )
    }
    inpt <- readline("next: ")
}
