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

## try the reverse complement
tmp.rc <- .Call("rev_complement", exons[[1]])

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

## lets reorder by rough phylogeny..
exon.i <- exon.i[ c('tetraodon nigroviridis', 'takifugu rubripes', 'electrophorus electricus', 'mastacembelus armatus',
                    'oryzias latipes', 'amphilophus citrinellus', 'neolamprologus brichardi',
                    'poecilia formosa', 'oreochromis niloticus', 'danio rerio',
                    'lepisosteus oculatus', 'callorhinchus milii',
                    'anolis carolinensis', 'gallus gallus', 'ursus americanus', 'mus musculus',
                    'homo sapiens')]

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

## try the stats function..
tmp <- .Call("nucl_align_stats", aligns[[1]]$seq)
tmp <- align.seqs.stats( aligns[[1]]$seq )

aligns.stats <- t(sapply( aligns, function(x){ align.seqs.stats( x$seq )}))


## we can then get the distances from these
align.stats.jc <- apply( aligns.stats, 1, jukes.cantor )
align.stats.k2 <- apply( aligns.stats, 1, kimura.two  )
align.stats.jcg <- apply( aligns.stats, 1, jukes.cantor.indel  )

## we can make a little matrix out of that..
jc.m <- matrix(0, nrow=length(exon.i), ncol=length(exon.i))
score.m <- jc.m
k2.m <- jc.m
jcg.m <- jc.m
for(i in 1:length(aligns)){
    jc.m[ aligns[[i]]$i, aligns[[i]]$j ] <- align.stats.jc[i]
    jc.m[ aligns[[i]]$j, aligns[[i]]$i ] <- align.stats.jc[i]
    ##
    k2.m[ aligns[[i]]$i, aligns[[i]]$j ] <- align.stats.k2[i]
    k2.m[ aligns[[i]]$j, aligns[[i]]$i ] <- align.stats.k2[i]
    ##
    jcg.m[ aligns[[i]]$i, aligns[[i]]$j ] <- align.stats.jcg[i]
    jcg.m[ aligns[[i]]$j, aligns[[i]]$i ] <- align.stats.jcg[i]
    score.m[ aligns[[i]]$i, aligns[[i]]$j ] <- aligns[[i]]$score
}

par(mfrow=c(1,2))
par(mar=c(10.1, 10.1, 4.1, 2.1))
image(1:length(exon.i), 1:length(exon.i), jc.m, col=hsvScale(1:255, sat=1:255/255 ), axes=FALSE, xlab=NA, ylab=NA, main='distance')
axis(1, at=1:length(exon.i), labels=names(exon.i), las=2, cex.axis=0.8)
axis(2, at=1:length(exon.i), labels=names(exon.i), las=2, cex.axis=0.8)
image(1:length(exon.i), 1:length(exon.i), score.m, col=hsvScale(1:255, sat=1:255/255 ), axes=FALSE, xlab=NA, ylab=NA, main='score' )
axis(1, at=1:length(exon.i), labels=names(exon.i), las=2, cex.axis=0.8)
axis(2, at=1:length(exon.i), labels=names(exon.i), las=2, cex.axis=0.8)


par(mar=c(10.1, 10.1, 4.1, 2.1))
image(1:length(exon.i), 1:length(exon.i), jcg.m, col=hsvScale(1:255, sat=1:255/255 ), axes=FALSE, xlab=NA, ylab=NA, main='distance')
axis(1, at=1:length(exon.i), labels=names(exon.i), las=2, cex.axis=0.8)
axis(2, at=1:length(exon.i), labels=names(exon.i), las=2, cex.axis=0.8)
image(1:length(exon.i), 1:length(exon.i), score.m, col=hsvScale(1:255, sat=1:255/255 ), axes=FALSE, xlab=NA, ylab=NA, main='score' )
axis(1, at=1:length(exon.i), labels=names(exon.i), las=2, cex.axis=0.8)
axis(2, at=1:length(exon.i), labels=names(exon.i), las=2, cex.axis=0.8)

par(mar=c(10.1, 10.1, 4.1, 2.1))
image(1:length(exon.i), 1:length(exon.i), k2.m, col=hsvScale(1:255, sat=1:255/255 ), axes=FALSE, xlab=NA, ylab=NA, main='distance')
axis(1, at=1:length(exon.i), labels=names(exon.i), las=2, cex.axis=0.8)
axis(2, at=1:length(exon.i), labels=names(exon.i), las=2, cex.axis=0.8)
image(1:length(exon.i), 1:length(exon.i), score.m, col=hsvScale(1:255, sat=1:255/255 ), axes=FALSE, xlab=NA, ylab=NA, main='score' )
axis(1, at=1:length(exon.i), labels=names(exon.i), las=2, cex.axis=0.8)
axis(2, at=1:length(exon.i), labels=names(exon.i), las=2, cex.axis=0.8)

## do we get a reasonable tree from these distances ?
require('ape')

dimnames(jc.m) <- list( names(exon.i), names(exon.i) )
dimnames(jcg.m) <- list( names(exon.i), names(exon.i) )
dimnames(k2.m) <- list( names(exon.i), names(exon.i) )

jc.nj <- nj(jc.m)
jcg.nj <- nj(jcg.m)
k2.nj <- nj(k2.m)

par(mar=c(5.1, 4.1, 4.1, 2.1))
par(mfrow=c(1,3))
plot( root(jc.nj, outgroup='callorhinchus milii'), cex=1.25)
plot( root(jcg.nj, outgroup='callorhinchus milii'), cex=1.25)
plot( root(k2.nj, outgroup='callorhinchus milii'), cex=1.25)

par(mar=c(5.1, 4.1, 4.1, 2.1))
par(mfrow=c(1,3))
plot( root(jc.nj, outgroup='callorhinchus milii'), cex=1.25, 'radial')
plot( root(jcg.nj, outgroup='callorhinchus milii'), cex=1.25, 'radial')
plot( root(k2.nj, outgroup='callorhinchus milii'), cex=1.25, 'radial')


for(i in 1:length(aligns)){
    sp1 <- names(exon.i)[ aligns[[i]]$i ]
    sp2 <- names(exon.i)[ aligns[[i]]$j ]
    plot.new()
    plot.window(xlim=c(-1, nchar(aligns[[i]]$seq[1])), ylim=c(0, 5))
    draw.aligns(aligns[[i]], 3, 1, 1, cols, sp1, sp2)
    print( aligns[[i]]$stats )
    axis(1)
    inpt <- readline('next: ')
    if(inpt == 'q')
        break
}


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


## let us test the mulithreaded version
dyn.load("src/exon_aligneR.so")
source('functions.R')

system.time(
    aligns.2 <- align.seqs.mt( exon.i[1], exon.i[2:length(exon.i)], al.offset, al.size, sub.matrix, gap=as.integer(c(-10, -1)),
                         tgaps.free=TRUE, sp.char="I", thread.n=6)
)
## 0.053 seconds

aligns.2.stats <- t( sapply(aligns.2, align.mt.to.stats) )

## to confirm that we get the correct data set..
par(mfrow=c(2,1))
for(i in 1:length(aligns.2)){
    plot.new()
    sp1 <- names(aligns.2)[1]
    sp2 <- names(aligns.2)[i+1]
    plot.window(xlim=c(-1, nchar(aligns.2[[i]]$seq[1])), ylim=c(0, 5))
    draw.aligns(aligns.2[[i]], 3, 1, 1, cols, sp1, sp2)
    axis(1)
    print( aligns.2[[i]]$stats )
    plot.new()
    plot.window(xlim=c(-1, nchar(aligns[[i]]$seq[1])), ylim=c(0, 5))
    draw.aligns(aligns[[i]], 3, 1, 1, cols, sp1, sp2)
    print( aligns[[i]]$stats )
    axis(1)
    inpt <- readline('next: ')
    if(inpt == 'q')
        break
}



## it seems that this works. Let us see how threaded it is..
## I have 6 cores and 12 threads on this machine. Let's see how it works with this
### This is with default compilation parameters
t.perf <- sapply(1:12, function(t){
    system.time(
        tmp <- align.seqs.mt( exon.i[1], exon.i[2:length(exon.i)], al.offset, al.size, sub.matrix, gap=as.integer(c(-10, -1)),
                         tgaps.free=TRUE, sp.char="I", thread.n=t)
    )
})

plot(1:ncol(t.perf), t.perf['elapsed',], type='b')
plot(1:ncol(t.perf), t.perf['elapsed',1] / t.perf['elapsed',], type='b')
abline(0,1, lty=2, col='red')
plot(1:ncol(t.perf), 1:ncol(t.perf) * t.perf['elapsed',] / t.perf['elapsed',1], type='b')




## lets make a bigger set of sequences to do this against and see if we get a better
## throughput...
b.seq <- rep(exon.i, 10)

t.perf.2 <- sapply(1:12, function(t){
    system.time(
        tmp <- align.seqs.mt( exon.i[1], b.seq, al.offset, al.size, sub.matrix, gap=as.integer(c(-10, -1)),
                         tgaps.free=TRUE, sp.char="I", thread.n=t)
    )
})


plot(1:ncol(t.perf.2), t.perf.2['elapsed',], type='b')
plot(1:ncol(t.perf.2), t.perf.2['elapsed',1] / t.perf.2['elapsed',], type='b')
abline(0,1, lty=2, col='red')
plot(1:ncol(t.perf.2), 1:ncol(t.perf.2) * t.perf.2['elapsed',] / t.perf.2['elapsed',1], type='b')

## that does give us rather better scaling. Up to a maximum of 5x speedup on 6 cores.
## Which in essence means that I can do this later on..

par(mfrow=c(1,2))
plot(1:ncol(t.perf), t.perf['elapsed',1] / t.perf['elapsed',], type='b', xlab='threads', ylab='relative performance')
abline(0,1,col='red', lty=2)
plot(1:ncol(t.perf.2), t.perf.2['elapsed',1] / t.perf.2['elapsed',], type='b', xlab='threads', ylab='relative performance')
abline(0,1,col='red', lty=2)

## after compiling with
## MAKEFLAGS="CFLAGS=-03" R CMD SHLIB exon_aligneR.c needleman_wunsch.c
dyn.load("src/exon_aligneR.so")
t.perf.2.03 <- sapply(1:12, function(t){
    system.time(
        tmp <- align.seqs.mt( exon.i[1], b.seq, al.offset, al.size, sub.matrix, gap=as.integer(c(-10, -1)),
                         tgaps.free=TRUE, sp.char="I", thread.n=t)
    )
})

t.perf.2['elapsed',] / t.perf.2.03['elapsed',]

## that makes no differece whatsoever..
## with which_i_max defined as static inline
dyn.load("src/exon_aligneR.so")
t.perf.2.o3.i <- sapply(1:12, function(t){
    system.time(
        tmp <- align.seqs.mt( exon.i[1], b.seq, al.offset, al.size, sub.matrix, gap=as.integer(c(-10, -1)),
                         tgaps.free=TRUE, sp.char="I", thread.n=t)
    )
})

t.perf.2.o3.i['elapsed',] / t.perf.2['elapsed',]
## this makes a marginal difference.. takes 10% less time.
## [1] 0.9490085 0.9883041 0.9140625 0.9385666 0.9268293 0.8914027 0.9323671
## [8] 0.9086538 0.9064039 0.8632075 0.8838384 0.8871795

## what abot with O4 ?
dyn.load("src/exon_aligneR.so")
t.perf.2.o4.i <- sapply(1:12, function(t){
    system.time(
        tmp <- align.seqs.mt( exon.i[1], b.seq, al.offset, al.size, sub.matrix, gap=as.integer(c(-10, -1)),
                         tgaps.free=TRUE, sp.char="I", thread.n=t)
    )
})

t.perf.2.o4.i['elapsed',] / t.perf.2['elapsed',]
## [1] 0.8866856 1.0097466 0.9348958 0.9419795 0.9268293 0.8959276 0.9420290
## [8] 0.8942308 0.9162562 0.8443396 0.8939394 0.8769231

## again, this makes no real difference in performance.. 

## declare some variables outside of the looping.. 
dyn.load("src/exon_aligneR.so")
t.perf.3.o4.i <- sapply(1:12, function(t){
    system.time(
        tmp <- align.seqs.mt( exon.i[1], b.seq, al.offset, al.size, sub.matrix, gap=as.integer(c(-10, -1)),
                         tgaps.free=TRUE, sp.char="I", thread.n=t)
    )
})

t.perf.3.o4.i['elapsed',] / t.perf.2['elapsed',]

## as expected that makes bugger all difference. The compiler is smart enough for
## that not to make any difference.

## revert to local declaration
dyn.load("src/exon_aligneR.so")
t.perf.4.o4.i <- sapply(1:12, function(t){
    system.time(
        tmp <- align.seqs.mt( exon.i[1], b.seq, al.offset, al.size, sub.matrix, gap=as.integer(c(-10, -1)),
                         tgaps.free=TRUE, sp.char="I", thread.n=t)
    )
})

t.perf.4.o4.i['elapsed',] / t.perf.2['elapsed',]
## and back to where we started.. With 6 threads we take 11% less time

## remove the use of which_max_i by having a couple of local conditions..
dyn.load("src/exon_aligneR.so")
t.perf.5.o4 <- sapply(1:12, function(t){
    system.time(
        tmp <- align.seqs.mt( exon.i[1], b.seq, al.offset, al.size, sub.matrix, gap=as.integer(c(-10, -1)),
                         tgaps.free=TRUE, sp.char="I", thread.n=t)
    )
})

t.perf.5.o4['elapsed',] / t.perf.2['elapsed',]
## OK, that is faster by about 15 to 20%. 4 days instead of 5.. 

## avoid setting the full score and table to 0 using memset..
dyn.load("src/exon_aligneR.so")
t.perf.6.o4 <- sapply(1:12, function(t){
    system.time(
        tmp <- align.seqs.mt( exon.i[1], b.seq, al.offset, al.size, sub.matrix, gap=as.integer(c(-10, -1)),
                         tgaps.free=TRUE, sp.char="I", thread.n=t)
    )
})

t.perf.6.o4['elapsed',] / t.perf.2['elapsed',]
## [1] 0.7686497 0.8440546 0.8385417 0.8191126 0.8008130 0.7782805 0.8115942
## [8] 0.7884615 0.7931034 0.7358491 0.7676768 0.7846154

## that gives quite variable results, and in general has more effect for a single
## thread than for others. Up to 75%.


## Lets try the new smith waterman alignment

sw.aligns <- .Call("sw_aligns", exon.i[1], exon.i[2], sub.matrix$offset, sub.matrix$size,
                   sub.matrix$sm, as.integer(c(-10, -1)), 10L, 10L )

plot.new()
plot.window( xlim=c(0,nchar(exon.i[2])), ylim=c(0, nchar(exon.i[1])), asp=1 )
rect(0, 0, nchar(exon.i[2]), nchar(exon.i[1]))
segments( sw.aligns[[3]][,3], sw.aligns[[3]][,1], sw.aligns[[3]][,4], sw.aligns[[3]][,2], col=hsvScale(sw.aligns[[3]][,5]), lwd=3 )
#abline( v=sw.aligns[[3]][,3] )
abline( h=sw.aligns[[3]][,1] )
abline( h=sw.aligns[[3]][,2] )
## so that seems to be ok, but the 

sw.aligns <- .Call("sw_aligns", exon.i[1], exon.i[16], al.offset, al.size, sub.matrix, as.integer(c(-10, -1)), 6L, 20L )


system.time(
    for(i in 1:100){
        sw.aligns <- .Call("sw_aligns", exon.i[1], exon.i[16], al.offset, al.size, sub.matrix, as.integer(c(-10, -1)), 6L, 20L )
    }
)
##  user  system elapsed 
## 1.194   0.052   1.247 

nchar( exon.i[1] ) ## 625
nchar( exon.i[2] ) ## 337
## so these are short and the 100 alignments per second is not particularly great..
## but we know that this is not a fast function, because of lots of extra stuff
## being done. 

plot.new()
plot.window( xlim=c(0,nchar(exon.i[16])), ylim=c(0, nchar(exon.i[1])), asp=1 )
rect(0, 0, nchar(exon.i[16]), nchar(exon.i[1]))
segments( sw.aligns[[3]][,3], sw.aligns[[3]][,1], sw.aligns[[3]][,4], sw.aligns[[3]][,2], col=hsvScale(sw.aligns[[3]][,5]), lwd=3 )
#abline( v=sw.aligns[[3]][,3] )
##abline( h=sw.aligns[[3]][,2] )

sw.aligns <- .Call("sw_aligns", substr(exon.i[1], 1, 20), substr(exon.i[16], 1, 20),
                   al.offset, al.size, sub.matrix, as.integer(c(-10, -1)), FALSE )

sw.aligns <- local.aligns( exon.i[1], exon.i[16], sub.matrix$offset, sub.matrix$size, sub.matrix$sm,
                           c(-10, -1), 20L, 40L)

plot.new()
plot.window( xlim=c(0,nchar(exon.i[16])), ylim=c(0, nchar(exon.i[1])), asp=1 )
rect(0, 0, nchar(exon.i[16]), nchar(exon.i[1]))
segments( sw.aligns[[3]][,3], sw.aligns[[3]][,1], sw.aligns[[3]][,4], sw.aligns[[3]][,2], col=hsvScale(sw.aligns[[3]][,5]), lwd=3 )


### try the local.score function.

seq <- c("ACTAGAAIAA-----ACAT",
         "AC-GGAAI--AA---ACAT")

local.score( seq, 3L, -8L, sub.matrix )

plot(local.score( seq, 3L, -8L, sub.matrix ), type='b')
