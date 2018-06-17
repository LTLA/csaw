# This tests the basic internal functions on which all others are built.
# library(testthat); library(csaw); source("test-basic.R")

source("simsam.R")

######################################################################
# Single-end read extraction.

library(Rsamtools)
SREF <- function(bam, param) {
    out <- scanBam(bam)[[1]]
    
    if (!is.na(param$minq)) {
        keep <- out$mapq >= param$minq
        out <- lapply(out, "[", keep)
    }

    if (param$dedup) {
        keep <- bitwAnd(out$flag, 0x400)==0
        out <- lapply(out, "[", keep)
    }

    if (!is.na(param$forward)) {
        if (param$forward) {
            keep <- bitwAnd(out$flag, 0x10)==0
        } else {
            keep <- bitwAnd(out$flag, 0x10)!=0
        }
        out <- lapply(out, "[", keep)
    }

    out$rwidth <- GenomicAlignments::cigarWidthAlongReferenceSpace(out$cigar)

    if (length(param$discard)) {
        pretend <- GRanges(out$rname, IRanges(out$pos, out$pos + out$rwidth - 1L))
        keep <- !overlapsAny(pretend, param$discard, type="within")
        out <- lapply(out, "[", keep)
    }

    return(out)
}

set.seed(1000)
test_that(".extractSE works correctly in the single-end case", {
    chromos <- c(chrA=1000, chrB=3000, chrC=200)
    genome <- GRanges(names(chromos), IRanges(1L, chromos))
    discarder <- makeDiscard(10, 100, chromos)

    for (rparam in list(readParam(), 
                        readParam(minq=10),
                        readParam(dedup=TRUE),
                        readParam(forward=FALSE),
                        readParam(forward=TRUE),
                        readParam(discard=discarder),
                        readParam(minq=10, dedup=TRUE, discard=discarder)
                        )) {

        for (nreads in c(100, 1000, 10000)) {
            obam <- regenSE(nreads, chromos, outfname=tempfile())
            ref <- SREF(obam, rparam)
            
            for (idx in seq_along(genome)) {
                element <- genome[idx]
                out <- csaw:::.extractSE(obam, element, rparam)

                cur.chr <- as.logical(seqnames(element)==ref$rname)
                expect_identical(out$forward$pos, ref$pos[cur.chr & ref$strand=="+"])
                expect_identical(out$forward$qwidth, ref$rwidth[cur.chr & ref$strand=="+"])
                expect_identical(out$reverse$pos, ref$pos[cur.chr & ref$strand=="-"])
                expect_identical(out$reverse$qwidth, ref$rwidth[cur.chr & ref$strand=="-"])
            }
        }
    }

    # Testing what happens with no reads.
    obam <- regenSE(0L, chromos, outfname=tempfile())
    for (idx in seq_along(genome)) {
        element <- genome[idx]
        out <- csaw:::.extractSE(obam, element, readParam())
        expect_identical(out$forward$pos, integer(0))
        expect_identical(out$forward$qwidth, integer(0))
        expect_identical(out$reverse$pos, integer(0))
        expect_identical(out$reverse$qwidth, integer(0))
    }

    expect_error(csaw:::.extractSE(obam, element, readParam(forward=NULL)), "strand extraction")
})

######################################################################
# Paired-end read extraction.

PREF <- function(bam, param) {
    raw <- SREF(bam, param)

    fkeep <- bitwAnd(raw$flag, 0x40)!=0L
    first <- lapply(raw, "[", fkeep)
    skeep <- bitwAnd(raw$flag, 0x80)!=0L
    second <- lapply(raw, "[", skeep)

    if (param$pe=="first") {
        return(first)
    } else if (param$pe=="second") {
        return(second)
    } 

    # Matching names.
    inboth <- intersect(first$qname, second$qname)
    m1 <- match(inboth, first$qname)
    first <- lapply(first, "[", m1)
    m2 <- match(inboth, second$qname)
    second <- lapply(second, "[", m2)

    # Constructing read pairs.
    lower <- ifelse(first$strand=="+", first$pos, second$pos)
    upper <- ifelse(first$strand=="+", second$pos + second$rwidth, first$pos + first$rwidth)
    size <- upper - lower 

    is.valid <- first$rname==second$rname & first$strand!=second$strand & size > 0L & size <= param$max.frag
    return(list(chr=first$rname[is.valid], pos=lower[is.valid], size=size[is.valid]))
}

set.seed(1001)
test_that(".extractPE works correctly in the paired-end case", {
    chromos <- c(chrA=1000, chrB=3000)
    genome <- GRanges(names(chromos), IRanges(1L, chromos))
    discarder <- makeDiscard(10, 100, chromos)

    for (rparam in list(readParam(pe="both"), 
                        readParam(pe="both", minq=10),
                        readParam(pe="both", dedup=TRUE),
                        readParam(pe="both", discard=discarder),
                        readParam(pe="both", max.frag=200),
                        readParam(pe="both", max.frag=1e8),
                        readParam(pe="both", minq=10, dedup=TRUE, discard=discarder)
                        )) {

        for (npairs in c(100, 1000, 10000)) {
            obam <- regenPE(npairs, chromos, outfname=tempfile())
            ref <- PREF(obam, rparam)
            
            for (idx in seq_along(genome)) {
                element <- genome[idx]
                out <- csaw:::.extractPE(obam, element, rparam)

                cur.chr <- as.logical(seqnames(element)==ref$chr)
                ref.pos <- ref$pos[cur.chr]
                ref.size <- ref$size[cur.chr]

                o1 <- order(out$pos, out$size)
                o2 <- order(ref.pos, ref.size)
                expect_identical(out$pos[o1], ref.pos[o2])
                expect_identical(out$size[o1], ref.size[o2])
            }
        }
    }

    # Testing what happens with no reads.
    obam <- regenPE(0L, chromos, outfname=tempfile())
    for (idx in seq_along(genome)) {
        element <- genome[idx]
        out <- csaw:::.extractPE(obam, element, readParam())
        expect_identical(out$pos, integer(0))
        expect_identical(out$size, integer(0))
    }

    # Strand can only be NA.
    expect_error(csaw:::.extractPE(obam, element, readParam(forward=NULL)), "cannot specify read strand")
    expect_error(csaw:::.extractPE(obam, element, readParam(forward=TRUE)), "cannot specify read strand")
    expect_error(csaw:::.extractPE(obam, element, readParam(forward=FALSE)), "cannot specify read strand")
})

set.seed(1002)
test_that(".extractSE works correctly with only first or second reads", {
    chromos <- c(chrA=1000, chrB=3000)
    genome <- GRanges(names(chromos), IRanges(1L, chromos))
    discarder <- makeDiscard(10, 100, chromos)

    for (rparam in list(readParam(pe="first", dedup=TRUE),
                        readParam(pe="second", minq=10)
                        )) {

        for (npairs in c(100, 1000, 10000)) {
            obam <- regenPE(npairs, chromos, outfname=tempfile())
            ref <- PREF(obam, rparam)
            
            for (idx in seq_along(genome)) {
                element <- genome[idx]
                out <- csaw:::.extractSE(obam, element, rparam)

                cur.chr <- as.logical(seqnames(element)==ref$rname)
                expect_identical(out$forward$pos, ref$pos[cur.chr & ref$strand=="+"])
                expect_identical(out$forward$qwidth, ref$rwidth[cur.chr & ref$strand=="+"])
                expect_identical(out$reverse$pos, ref$pos[cur.chr & ref$strand=="-"])
                expect_identical(out$reverse$qwidth, ref$rwidth[cur.chr & ref$strand=="-"])
            }
        }
    }
})

######################################################################
# Chromosome length and end coercion checks.

set.seed(1003)
test_that(".activeChrs works correctly", {
    chromos <- c(chrA=1000L, chrB=3000L)
    obam1 <- regenSE(100, chromos, outfname=tempfile())
    obam2 <- regenPE(100, chromos, outfname=tempfile())
    obam3 <- regenPE(100, c(chromos, chrC=200), outfname=tempfile())
    obam4 <- regenPE(100, chromos * 2, outfname=tempfile())

    # In the simplest case... 
    expect_identical(csaw:::.activeChrs(obam1, NULL), chromos)
    expect_identical(csaw:::.activeChrs(c(obam1, obam2), NULL), chromos)
    expect_identical(csaw:::.activeChrs(c(obam1, obam2), "chrA"), chromos["chrA"])
    expect_identical(csaw:::.activeChrs(c(obam1, obam2), "chrB"), chromos["chrB"])
    expect_identical(csaw:::.activeChrs(c(obam1, obam2), names(chromos)), chromos)

    # Triggering various warnings.
    expect_warning(out <- csaw:::.activeChrs(c(obam1, obam3), NULL), "not identical")
    expect_identical(out, chromos)
    expect_warning(out <- csaw:::.activeChrs(c(obam1, obam4), NULL), "not identical")
    expect_identical(out, chromos)
})

set.seed(1004)
test_that(".coerceFragments works correctly", {
    starts <- as.integer(round(runif(100, 1, 1000)))
    widths <- as.integer(round(runif(100, 1, 100)))
    ends <- starts + widths

    # Chrlen and positivity are enforced.
    chrlen <- 1000L
    out <- csaw:::.coerceFragments(starts, ends, final=NA, chrlen=chrlen)
    expect_identical(out$end, pmin(ends, chrlen))
    expect_identical(out$start, starts)

    out <- csaw:::.coerceFragments(-starts, ends, final=NA, chrlen=1e5)
    expect_true(all(out$start==1L))
    expect_identical(out$end, ends)

    # Final width is enforced.    
    final <- 50
    out <- csaw:::.coerceFragments(starts, ends, final=final, chrlen=1e6)
    ref <- resize(IRanges(starts, ends), width=final, fix="center")
    expect_identical(out$start, pmax(1L, start(ref)))
    expect_identical(out$end, end(ref))

    chrlen <- 500L
    final <- 21 # odd width
    out <- csaw:::.coerceFragments(starts, ends, final=final, chrlen=chrlen)
    ref <- resize(IRanges(starts, ends), width=final, fix="center")
    expect_identical(out$start, pmin(pmax(1L, start(ref)), chrlen))
    expect_identical(out$end, pmin(end(ref), chrlen))

    # Works correcty on empty inputs.
    out <- csaw:::.coerceFragments(integer(0), integer(0), final=NA, chrlen=1e6)
    expect_identical(out$start, integer(0))
    expect_identical(out$end, integer(0))

    out <- csaw:::.coerceFragments(integer(0), integer(0), final=50L, chrlen=1e6)
    expect_identical(out$start, integer(0))
    expect_identical(out$end, integer(0))
})

######################################################################
# Directional read extension.

test_that(".collateExt works correctly", {
    nbam <- 4

    # Simple number.
    ext <- 100
    out <- csaw:::.collateExt(nbam, ext)
    expect_identical(out$ext, rep(as.integer(ext), nbam))
    expect_identical(out$final, NA_integer_)

    variety <- c(100, 40, 20, 80)
    expect_error(csaw:::.collateExt(nbam, variety), "must be an integer scalar")
    expect_error(csaw:::.collateExt(nbam, -ext), "positive")

    # A list of length 2.
    final <- 50L
    out <- csaw:::.collateExt(nbam, list(variety, final))
    expect_identical(out$ext, as.integer(variety))
    expect_identical(out$final, final)

    expect_error(csaw:::.collateExt(nbam, list(ext)), "must be a list of length 2")
    expect_error(csaw:::.collateExt(nbam, list(ext, final)), "not consistent with the number of libraries")
    expect_error(csaw:::.collateExt(nbam, list(variety, -final)), "positive")
})

test_that(".extendSE works correctly", {
    # Simulating reads.
    fstarts <- as.integer(round(runif(100, 1, 1000)))
    fwidths <- as.integer(round(runif(100, 1, 100)))
    rstarts <- as.integer(round(runif(50, 1, 1000)))
    rwidths <- as.integer(round(runif(50, 1, 100)))
    reads <- list(forward=list(pos=fstarts, qwidth=fwidths), reverse=list(pos=rstarts, qwidth=rwidths))
    
    # Running it through the extender with no extension.
    chrlen <- 800L
    out <- csaw:::.extendSE(reads, ext=NA, final=NA, chrlen=chrlen)

    new.fstarts <- pmin(fstarts, chrlen)
    new.fends <- pmin(fstarts + fwidths - 1L, chrlen)
    new.rstarts <- pmin(rstarts, chrlen)
    new.rends <- pmin(rstarts + rwidths - 1L, chrlen)

    expect_identical(out$start, c(new.fstarts, new.rstarts))
    expect_identical(out$end, c(new.fends, new.rends))

    # Running it through the extender.
    ext <- 80L
    chrlen <- 800L
    out <- csaw:::.extendSE(reads, ext=ext, final=NA, chrlen=chrlen)

    new.fstarts <- pmin(fstarts, chrlen)
    new.fends <- pmin(fstarts + ext - 1L, chrlen)
    new.rends <- rstarts + rwidths 
    new.rstarts <- pmin(pmax(1L, new.rends - ext), chrlen)
    new.rends <- pmin(new.rends - 1L, chrlen)

    expect_identical(out$start, c(new.fstarts, new.rstarts))
    expect_identical(out$end, c(new.fends, new.rends))
    
    # Testing it with a final setting.
    final <- 50L
    ext <- 80L
    chrlen <- 800L
    out <- csaw:::.extendSE(reads, ext=ext, final=final, chrlen=chrlen)

    new.fwidth <- rep(ext, length(fstarts))
    new.rends <- rstarts + rwidths 
    new.rstarts <- new.rends - ext
    new.rwidth <- rep(ext, length(rstarts))
    new.reads <- list(forward=list(pos=fstarts, qwidth=new.fwidth),
                      reverse=list(pos=new.rstarts, qwidth=new.rwidth))

    ref <- csaw:::.extendSE(new.reads, ext=NA, final=final, chrlen=chrlen)
    expect_identical(out, ref)
})

######################################################################
# Other miscellaneous functions.

test_that("strand decisions are made correctly", {
    expect_error(csaw:::.decideStrand(readParam(forward=NULL)), "unspecified strandedness")
    expect_identical(csaw:::.decideStrand(readParam(forward=TRUE)), "+")
    expect_identical(csaw:::.decideStrand(readParam(forward=FALSE)), "-")
    expect_identical(csaw:::.decideStrand(readParam(forward=NA)), "*")
})

test_that("colData formatter works correctly", {
    bf <- "whee.bam"
    total <- 1011021
    ext <- list(ext=5L)
    all.extras <- list(list(c(1, 10), c(2, 20)))

    out <- csaw:::.formatColData(bf, total, ext.data=ext, all.extras=all.extras, readParam())
    expect_identical(out$bam.files, bf)
    expect_identical(out$totals, total)
    expect_identical(out$ext, ext$ext)
    expect_equal(out$rlen, round(weighted.mean(1:2, 1:2*10)))

    # Trying with multiple libraries.
    bf <- c("whee.bam", "yay.bam")
    total <- c(11021, 341234)
    ext <- list(ext=c(5L, 10L))
    all.extras <- list(list(c(1, 10), c(2, 20)), list(c(3, 30), c(4, 40)))

    out <- csaw:::.formatColData(bf, total, ext.data=ext, all.extras=all.extras, readParam())
    expect_identical(out$bam.files, bf)
    expect_identical(out$totals, total)
    expect_identical(out$ext, ext$ext)
    expect_equal(out$rlen[1], round(weighted.mean(1:2, 1:2*10)))
    expect_equal(out$rlen[2], round(weighted.mean(3:4, 3:4*10)))
})

set.seed(1005)
test_that("automated range conversion works correctly", {
    nrows <- 200
	ncols <- 6
    counts <- matrix(runif(nrows * ncols, 1, 1e4), nrows)
    rowRanges <- GRanges(rep(c("chr1", "chr2"), c(50, 150)),
                         IRanges(floor(runif(200, 1e5, 1e6)), width=100),
                         strand=sample(c("+", "-"), 200, TRUE))
    rse <- SummarizedExperiment(list(counts=counts), rowRanges=rowRanges)
    
    expect_identical(rowRanges(rse), csaw:::.toGRanges(rse))
    expect_identical(rowRanges(rse), csaw:::.toGRanges(rowRanges(rse)))
    expect_error(csaw:::.toGRanges("whee"), "must be a RangedSummarizedExperiment")
})

