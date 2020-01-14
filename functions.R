## functions specific to exon alignments

## A horribly slow way to read transcripts.
read.transcripts <- function(fname){
    exons.l <- readLines(fname)
    id.l <- grep(">", exons.l)
    s.e <- c(id.l - 1, length(exons.l))[-1]
    exons <- lapply( 1:length(id.l), function(i){
        exons.l[(id.l[i]+1):s.e[i]]
    })
    names(exons) <- sapply( exons.l[id.l], function(x){ strsplit(x, split="\t")[[1]][1] })
    exon.dbs <- unname(sapply( exons.l[id.l], function(x){ strsplit(x, split="\t")[[1]][4] }))

    lapply( 1:length(exons), function(i){
        list('id'=sub("^>", "", names(exons)[i]), 'db'=exon.dbs[i], 'sp'=sub("([^_]+)_([^_]+)_.+$", "\\1 \\2", exon.dbs[i]),
             'l'=unname(sapply(exons[[i]], nchar)), 's'=paste(exons[[i]], collapse='') )
    })
}

read.exons <- function(fname, ex.spacer='I'){
    exons.l <- readLines(fname)
    id.l <- grep(">", exons.l)
    s.e <- c(id.l - 1, length(exons.l))[-1]
    exons <- lapply( 1:length(id.l), function(i){
        exons.l[(id.l[i]+1):s.e[i]]
    })
    names(exons) <- sapply( exons.l[id.l], function(x){ strsplit(x, split="\t")[[1]][1] })
    exon.dbs <- unname(sapply( exons.l[id.l], function(x){ strsplit(x, split="\t")[[1]][4] }))

    lapply( 1:length(exons), function(i){
        list('id'=sub("^>", "", names(exons)[i]), 'db'=exon.dbs[i], 'sp'=sub("([^_]+)_([^_]+)_.+$", "\\1 \\2", exon.dbs[i]),
             'l'=unname(sapply(exons[[i]], nchar)), e=exons[[i]], 's'=paste(toupper(exons[[i]]), collapse=ex.spacer) )
    })
}


## and and b are lists containing
## $s   a transcript sequence
## $l   the lengths of the exons making up $s
## penalties a vector giving:
## match, mismatch, gap insertion, gap extension
## gap.multiplier : a value used to calculate gap penalties
## for exon alignments
align.exons <- function(a, b, penalties, gap.multiplier, local=FALSE, local.gene=FALSE){
    if(!is.character(a$s) || !is.character(b$s))
        stop("a and b should be character vectors")
    if(!is.numeric(penalties) || !is.numeric(gap.multiplier))
        stop("penalties and the gap.multiplier should both be numeric vectors")
    if(length(penalties) != 4)
        stop("penalties should have 4 entries: match, mismatch, insertion, extension")
    if(length(gap.multiplier) != 1)
        stop("the gap.multiplier should be a single value")
    if(!is.logical(local))
        stop("local should be a logical vector indicating whether the alignment is local")
    al <- .Call("align_exons", a$s, b$s,
                as.integer(a$l), as.integer(b$l),
                as.double(penalties), as.double(gap.multiplier), local, local.gene )
##    names(al) <- c("exon.scores", "scores", "pointers", "align", "exons.a", "exons.b")
    names(al) <- c("exon.scores", "scores", "pointers", "g.align", "a", "b")
    al <- c('a.id'=a$id, 'b.id'=b$id, 'a.sp'=a$sp, 'b.sp'=b$sp, al)
    al
}



make.sub.matrix <- function( al.offset=33, al.size=64, letters=c('A', 'C', 'T', 'G', 'I'),
                            match = c(4, 4, 4, 4, 30), mismatch=-match ){
    m <- matrix(0, nrow=al.size, ncol=al.size)
    for(i in 1:length(letters)){
        for(j in 1:length(letters)){
            row <- 1 + utf8ToInt(letters[i]) - al.offset
            column <- 1 + utf8ToInt(letters[j]) - al.offset
            m[row,column] <- ifelse(i==j, match[i], min(mismatch[c(i,j)]) )
        }
    }
    m <- matrix(as.integer(m), nrow=nrow(m))
    list('offset'=as.integer(al.offset), 'size'=as.integer(al.size), 'sm'=m)
}

align.seqs <- function(a, b, al.offset, al.size,
                       sub.matrix=make.sub.matrix(al.offset=al.offset, al.size=al.size),
                       gap=as.integer(c(-10, -1)), tgaps.free=TRUE, sp.char="I"){
    al <- .Call("align_seqs", a, b, al.offset, al.size, sub.matrix, gap, tgaps.free, sp.char)
    names(al) <- c('scores', 'ptr', 'seq', 'a.pos', 'b.pos', 'al.i', 'stats')
    names(al$stats) <- c('length', 'al_n', 'a_gapi', 'b_gapi', 'a_gap', 'b_gap',
                         'a_gap_l', 'a_gap_r', 'b_gap_l', 'b_gap_r', 'match', 'mismatch')
    c(al, 'score'=al$scores[ nrow(al$scores), ncol(al$scores) ])
}

align.seqs.stats <- function(seqs){
    if(!is.character(seqs) || length(seqs) < 2)
        stop("Provide a character vector of two nucleic acid sequences");
    stats <- .Call("nucl_align_stats", seqs);
    names(stats) <- c('left_gaps', 'right_gaps', 'gaps', 'gap.i', 'ident', 'transition', 'transversion',
                      'A', 'C', 'G', 'T')
    stats <- c(stats, 'al.l'=sum( stats[c('ident', 'transition', 'transversion')] ) )
    class(stats) <- 'nuc.align.stats'
    stats
}

align.seqs.mt <- function(a, b, al.offset, al.size,
                          sub.matrix=make.sub.matrix(al.offset=al.offset, al.size=al.size),
                          gap=as.integer(c(-10, -1)), tgaps.free=TRUE, sp.char="I", thread.n=1){
    al <- .Call("align_seqs_mt", a, b, al.offset, al.size, sub.matrix, gap, tgaps.free, sp.char, as.integer(thread.n))
    names(al) <- names(b)
    for(i in 1:length(al)){
        names(al[[i]]) <- c('stats', 'seq', 'a.pos', 'b.pos', 'al.i');
        names(al[[i]]$stats) <- c('score', 'length', 'al_n', 'a_gapi', 'b_gapi', 'a_gap', 'b_gap',
                                  'a_gap_l', 'a_gap_r', 'b_gap_l', 'b_gap_r', 'match', 'mismatch',
                                  'transition', 'transversion', 'A', 'C', 'G', 'T')
    }
    al
}

align.mt.to.stats <- function( al ){
    as <- al$stats
    stats <- c('left_gaps' = as['a_gap_l'] + as['b_gap_l'],
               'right_gaps' = as['a_gap_r'] + as['b_gap_r'],
               'gaps' = as['a_gap'] + as['b_gap'],
               'gap.i' = as['a_gapi'] + as['b_gapi'],
               'ident' = as['match'], 'transition' = as['transition'], 'transversion'=as['transversion'],
               'A'=as['A'], 'C'=as['C'], 'G'=as['G'], 'T'=as['T'], 'al.l'=as['al_n'])
    names(stats) <- c('left_gaps', 'right_gaps', 'gaps', 'gap.i', 'ident', 'transition', 'transversion',
                      'A', 'C', 'G', 'T', 'al.l')
    class(stats) <- 'nuc.align.stats'
    stats
}

### phylogenetic distance functions. These all take an object of nuc.align.stats as defined above
### I got these equations from :
### http://bioinformatics.psb.ugent.be/downloads/psb/Userman/treecon_distance.html
###
### but I should get the original equations somewhere

## combine a list of stats to a single stat. Simple summing
merge.stats <- function( stats ){
    ss <- NULL
    b <- 1
    while(is.null(ss)){
        ss <- stats[[b]]
        b <- b + 1
    }
    for(i in b:length(stats)){
        if(!is.null(stats[[i]]))
            ss <- ss + stats[[i]]
    }
    ss
}
    

jukes.cantor <- function( stats ){
    (-3/4) * log( 1 - (4/3) * (stats['transition'] + stats['transversion']) / stats['al.l'] )
}

jukes.cantor.indel <- function(stats){
    G <- stats['gap.i'] - (as.logical(stats['left_gaps']) + as.logical(stats['right_gaps']))
    Fu <- (stats['transition'] + stats['transversion']) / stats['al.l']
    I <- stats['ident']
    T <- G + Fu + I
    (-3/4) * log( 1 - (4/3) * (Fu / (1 + Fu)) ) * (1 - G/T) + G/T
}

kimura.two <- function( stats ){
    (-1/2) * log( (1 - 2*stats['transition']/stats['al.l'] - stats['transversion']/stats['al.l']) * sqrt(1-2*stats['transversion']/stats['al.l']) )
}

tajima.nei <- function( stats ){
    l <- stats$al.l - stats$unknown
    b <- 1 - ( (stats['A']/l)^2 + (stats['C']/l)^2 + (stats['G']/l)^2 + (stats['T']/l)^2 )
    -b * log( 1 - (1/b) * (stats$transition + stats$transversion) / l )
}



## al should be a list of 5 components.
## with names:
## exon.scores, scores, pointers, align, exons.a, exons.b
plot.exon.alignments <- function(al, i.sep=1.25, g.height=6, h1=1, w.radius=4, w.sd=w.radius/2, sim.h=h1*0.75,
                                 ex.score.h=1){
    ## we set up two representations; on the top a map of the different
    ## exons simply showing where they align.
    a.seq <- strsplit( al$a, '' )
    b.seq <- strsplit( al$b, '' )
    ex.lengths <- unname( sapply( a.seq, length ))
    sep <- i.sep + h1
    y1 <- 1
    y2 <- y1 + sep
    
    plot.new()
    plot.window(xlim=c(0, sum(ex.lengths)), ylim = c(0, g.height * 1.25 + length(a.seq) * h1 * g.height + ex.score.h * nrow(al$exon.scores)))

    usr <- par("usr")
    top <- usr[4]
    
    ## fill in where we have aligned residues
    a.seq.c <- unlist( a.seq )
    b.seq.c <- unlist( b.seq )
    a.col <- ifelse( a.seq.c == '-', 'white', 'black' )
    b.col <- ifelse( b.seq.c == '-', 'white', 'black' )
    x2 <- 1:length( a.seq.c )
    x1 <- x2 - 1
    rect( x1, top - y1, x2, top - (y1+h1), col=a.col, border=NA )
    rect( x1, top - y2, x2, top - (y2+h1), col=b.col, border=NA )
    krn <- dnorm( -w.radius:w.radius, sd=w.sd )
    r.sim <- sapply(1:length(a.seq.c), function(i){
        b <- ifelse( i > w.radius, i - w.radius, 1 )
        e <- ifelse( i + w.radius <= length(a.seq.c), i + w.radius, length(a.seq.c) )
        k.b <- 1 + w.radius - (i-b);
        k.e <- 1 + w.radius + (e-i)
        sum( as.numeric(a.seq.c[b:e] == b.seq.c[b:e]) * krn[k.b:k.e] ) / sum(krn[k.b:k.e])
    })
    y3 <- (h1 + sim.h + y1 + y2) / 2
    lines( (x1+x2)/2, top - y3 + r.sim * sim.h )
    segments( (x1+x2)/2, top-y3, (x1+x2)/2, top - y3 + ifelse( a.seq.c == b.seq.c, r.sim * sim.h, 0 ) )
    segments( x1[1], top - y3, x2[length(x2)], top - y3 )
    
    ## To see the boundaries of the exons we do:
    x1 <- c(0, cumsum( ex.lengths ))[ -(length(ex.lengths)+1) ]
    x2 <- x1 + ex.lengths
    rect( x1, top - y1, x2, top - (y1+h1), border='red', lwd=2 )
    rect( x1, top - y2, x2, top - (y2+h1), border='red', lwd=2 )

    ex.y1 <- g.height * 1.25 + (1:length(a.seq)-1) * h1 * g.height
    ex.y2 <- ex.y1 + sep

    cols <- c('white', 'black', 'red', 'green', 'blue')
    names(cols) <- c('-', 'A', 'C', 'G', 'T')
    x.m <- sum( ex.lengths ) / max( ex.lengths )
    for(i in 1:length(ex.lengths)){
        x2 <- x.m * 1:ex.lengths[i]
        x1 <- x2 - x.m
        rect( x1, top-ex.y1[i], x2, top-(ex.y1[i] + h1), col=cols[ a.seq[[i]] ], border=NA )
        rect( x1, top-ex.y2[i], x2, top-(ex.y2[i] + h1), col=cols[ b.seq[[i]] ], border=NA )
        r.sim <- sapply(1:length(a.seq[[i]]), function(j){
            b <- ifelse( j > w.radius, j - w.radius, 1 )
            e <- ifelse( j + w.radius <= length(a.seq[[i]]), j + w.radius, length(a.seq[[i]]) )
            k.b <- 1 + w.radius - (j-b);
            k.e <- 1 + w.radius + (e-j)
            sum( as.numeric(a.seq[[i]][b:e] == b.seq[[i]][b:e]) * krn[k.b:k.e] ) / sum(krn[k.b:k.e])
        })
        x <- (x1 + x2)/2
        y <- (ex.y1[i] + ex.y2[i] + h1 + sim.h) / 2
        lines( x, top - y + r.sim * sim.h )
        segments( x, top-y, x, top - y + ifelse( a.seq[[i]] == b.seq[[i]], r.sim * sim.h, 0 ) )
        segments(x[1], top - y, x[length(x)], top-y)
        text(0, top - ex.y1[i] - h1/2, al$g.align[i,1], pos=2)
        text(0, top - ex.y2[i] - h1/2, al$g.align[i,2], pos=2)
    }
    ## from the bottom we can simply do:
    y <- rev(rep(1:nrow(al$exon.scores), ncol(al$exon.scores)) - 1)
    x <- matrix(1:ncol(al$exon.scores), nrow=nrow(al$exon.scores), ncol=ncol(al$exon.scores), byrow=TRUE ) - 1
    asp <- with(par(), (pin[1]/pin[2]) / (diff(usr[1:2])/diff(usr[3:4])) )
    v <- (al$exon.scores - min(al$exon.scores)) / diff(range(al$exon.scores))
    rect(ex.score.h * x/asp, y*ex.score.h, (x+1)*ex.score.h/asp, (1 + y)*ex.score.h, col=hsv(v, 1, v), border=NA)
    
    mtext( paste(al$a.sp, "vs", al$b.sp), side=3, line=2, cex=2 )
    mtext( paste(al$a.id, ", ", al$b.id, sep=""), side=3, line=0, cex=1 )
}

## cols must be a named vector, where the names are the residues present in the aligned sequences
## i.e. generally, ATCG, -, and other special residues (eg. introns)
## sim.pos determines where to draw the similarity plot.
##     (1,1) => draw underneath y1
##     (1,2) => draw above y1
##     (2,1) => draw underneath y2
##     (2,2) => draw above y2
draw.aligns <- function(al, y1, y2, h1, cols, sp.a=NULL, sp.b=NULL, w.radius=4, w.sd=w.radius/2, sim.h=h1*0.75, sim.sep=h1*0.125,
                        sim.pos=c(1,1), border=cols){
    krn <- dnorm( -w.radius:w.radius, sd=w.sd )
    a.seq <- strsplit( al$seq[1], '' )[[1]]
    b.seq <- strsplit( al$seq[2], '' )[[1]]
    sim <- sapply( 1:length(a.seq), function(i){
        b <- ifelse( i > w.radius, i - w.radius, 1 )
        e <- ifelse( i + w.radius <= length(a.seq), i + w.radius, length(a.seq) )
        k.b <- 1 + w.radius - (i-b)
        k.e <- 1 + w.radius + (e-i)
        sum( as.numeric(a.seq[b:e] == b.seq[b:e]) * krn[k.b:k.e] ) / sum(krn[k.b:k.e])
    })
    ## then draw our rectangles using the colours specified
    x2 <- 1:length(a.seq)
    x1 <- x2 - 1
    rect(x1, y1, x2, y1+h1, col=cols[ a.seq ], border=border[a.seq])
    rect(x1, y2, x2, y2+h1, col=cols[ b.seq ], border=border[b.seq])
    ## then we plot the similarity beneath y1, and at 
    ## sim always has a maximum of 1.
    sim.y.base <- ifelse( sim.pos[1] == 1, y1, y2 )
    sim.y <- sim.y.base + ifelse( sim.pos[2] == 1,  -(sim.h + sim.sep), (h1 + sim.sep) )
    lines( (x1+x2)/2, sim.y + sim * sim.h )
    segments( 0.5, sim.y, length(a.seq)-0.5, sim.y )
    segments( (x1+x2)/2, sim.y, (x1+x2)/2, sim.y + ifelse( a.seq == b.seq, sim * sim.h, 0 ) )
    if(!is.null(sp.a) || !is.null(sp.b))
        text( 0, y2, paste(sp.a, "vs", sp.b), adj=c(0,1.2) )
}

## al are two aligned sequences that ought to be the same size
## though we will make use of the first one
align.print <- function( al, w ){
    nc <- max(nchar(al[1:2]))
    beg <- seq(1, nc, w)
    for(b in beg){
        s1 <- substr(al[1], b, b+w)
        s2 <- substr(al[2], b, b+w)
        s1.c <- strsplit(s1, '')[[1]]
        s2.c <- strsplit(s2, '')[[1]]
        id.line <- paste( ifelse( s1.c == s2.c, "|", " " ), collapse='' )
        cat( sprintf("% 4d  %s\n", b, substr(al[1], b, b+w)) )
        cat( sprintf("      %s\n", id.line) )
        cat( sprintf("% 4d  %s\n", b, substr(al[2], b, b+w)) )
        cat("\n")
    }
}


