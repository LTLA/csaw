simsam <- function(f.out, chr, pos, strands, chromosomes, mapq=199, is.dup=NULL, names=NULL, is.paired=FALSE, len=10, cigar=NULL) 
# Constructs a simulated SAM file, compresses it to BAM and indexes it.
# Accepts a set of names, positions and strands, and fills in the rest if not supplied.
# Also handles creation of paired-end data with correct mate information.
{
    samFile <- paste0(f.out, ".sam")
    out <- file(samFile, open="w")

    # Adding SAM header.
    for (rname in names(chromosomes)) {
        write(c("@SQ", paste("SN:", rname, sep=""), paste("LN:", chromosomes[[rname]], sep="")), ncolumns=3, file=out, sep="\t");
    }

    # Setting up values in the CIGAR, names and flags.
    if (is.null(cigar)) { 
        cigar <- paste0(len, "M")
    }
    seq <- strrep("N", len)
    qual <- strrep(".", len)

    if (is.null(names)) { 
        names <- sprintf("x%i", seq_along(pos))
    }

    flags<-ifelse(strands, 0, 16)  
    if (!is.null(is.dup)) { 
        flags <- flags + ifelse(is.dup, 1024, 0) 
    }

    .expander <- function(x)  rep(x, length.out=length(names)) 
    stuff <- data.frame(QNAME=names, FLAG=flags, CHR=chr, POS=pos, MAPQ=.expander(mapq), CIGAR=.expander(cigar), stringsAsFactors=FALSE)

    # Setting a whole host of paired end information in the relevant fields.
    if (is.paired) {
        is.first <- !duplicated(stuff$QNAME)
        first.dex <- which(is.first)
        second.dex <- which(!is.first)

        first.m <- second.dex[match(stuff$QNAME[first.dex], stuff$QNAME[second.dex])]
        second.m <- first.dex[match(stuff$QNAME[second.dex], stuff$QNAME[first.dex])]
        counterpart <- integer(nrow(stuff))
        counterpart[first.dex] <- first.m
        counterpart[second.dex] <- second.m

        mate.strand <- bitwAnd(stuff$FLAG[counterpart], 0x10)==0L
        mate.chr <- stuff$CHR[counterpart]
        mate.pos <- stuff$POS[counterpart]

        reflags <- stuff$FLAG + 1
        reflags <- reflags + ifelse(is.first, 64, 128)
        reflags <- reflags + ifelse(mate.strand, 0, 32)  

        mate.chr.size <- chromosomes[mate.chr] + 1L # Keep above mate.chr reassignment!
        mate.chr <- ifelse(mate.chr==stuff$CHR, "=", mate.chr)
        my.cigar <- GenomicAlignments::cigarWidthAlongReferenceSpace(stuff$CIGAR)
        mate.cigar <- GenomicAlignments::cigarWidthAlongReferenceSpace(stuff$CIGAR[counterpart])

        isize <- pmin(pmax(mate.pos + mate.cigar, pos + my.cigar), mate.chr.size) - pmin(mate.pos, pos)
        isize[mate.pos < pos] <- isize[mate.pos < pos]*-1
        isize[mate.chr!="="] <- 0

        stuff$FLAG <- reflags
        stuff <- cbind(stuff, data.frame(MATECHR=mate.chr, MATEPOS=mate.pos, ISIZE=isize))
    } else {
        stuff <- cbind(stuff, data.frame(MATECHR=.expander("*"), MATEPOS=.expander(0), ISIZE=.expander(0)))
    }

    stuff <- cbind(stuff, data.frame(SEQ=.expander(seq), QUAL=.expander(qual), stringsAsFactors=FALSE))
    write.table(file=out, stuff, row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")
    close(out)

    # We sort and compress to BAM as well.
    require(Rsamtools)
    tempName <- paste(f.out, "_temp", sep="")
    tempName <- asBam(samFile, destination=tempName, overwrite=TRUE, indexDestination=FALSE)
    newName <- sortBam(tempName, destination=f.out)
    indexBam(newName)
    unlink(tempName)
    return(newName)
}

regenSE <- function(nreads, chromos, outfname, ...) 
# Convenient single-end simulations, auto-filling read positions and alignment statistics.
# We generate some complicated CIGAR strings with soft clips for some variety.
{
    chr <- sample(length(chromos), nreads, replace=TRUE)
    pos <- integer(nreads)
    str <- logical(nreads)

    for (i in seq_along(chromos)) {
        current <- chr==i
        pos[current] <- round(runif(sum(current), 1, chromos[i]))
        str[current] <- rbinom(sum(current), 1, 0.5)==1L
    }

    isdup <- rbinom(nreads, 1, 0.8)==0L
    mapq <- round(runif(nreads, 0, 40))

    len <- 20
    left.clip <- round(runif(nreads, 0, len/4))
    right.clip <- round(runif(nreads, 0, len/4))
    cigar <- sprintf("%iM", len - left.clip - right.clip)
    cigar[right.clip > 0] <- sprintf("%s%sS", cigar[right.clip > 0], right.clip[right.clip > 0])
    cigar[left.clip > 0] <- sprintf("%sS%s", left.clip[left.clip > 0], cigar[left.clip > 0])

    simsam(outfname, names(chromos)[chr], pos, str, chromos, is.dup=isdup, mapq=mapq, cigar=cigar, len=len, ...)
}

regenPE <- function(npairs, chromos, outfname) 
# Convenient paired-end simulations, auto-filling read positions and alignment statistics.
# The number of reads is simply double the number of pairs involved. Note the use of seq_len
# to provide consistent behaviour when npairs = 0.
{

    qnames <- sprintf("READ%i", c(sample(seq_len(npairs)), sample(seq_len(npairs))))
    regenSE(npairs * 2L, chromos, outfname, is.paired=TRUE, names=qnames)
}

makeDiscard <- function(ndisc, sizeof, chromos) 
# Construct a set of blacklisted regions to use in read extraction.
{
    chosen <- sample(length(chromos), ndisc, replace=TRUE)
    chosen.pos <- runif(ndisc, 1, chromos[chosen]-sizeof)
    reduce(GRanges(names(chromos)[chosen], IRanges(chosen.pos, chosen.pos+sizeof)))
}

generateWindows <- function(chrs, nwin, winsize) 
# Construct genomic windows of specified number and size for general use.
{
    allregs<-GRanges()
    for (x in names(chrs)) {
        max.step<-floor(chrs[[x]]/nwin)
        stopifnot(max.step >= 1)
        pos<-cumsum(round(runif(nwin, 1, max.step)))
        suppressWarnings(allregs<-c(allregs, GRanges(x, IRanges(pos, 
            pmin(chrs[[x]], pos+winsize-1L)))))
    }
    total.n<-nwin*length(chrs)
    return(allregs)
}
