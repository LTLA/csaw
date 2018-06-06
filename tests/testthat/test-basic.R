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
    out <- csaw:::.coerceFragments(starts, ends, final=NA, chrlen=1000)
    expect_identical(out$end, pmin(ends, 1000L))
    expect_identical(out$start, starts)

    out <- csaw:::.coerceFragments(-starts, ends, final=NA, chrlen=1e5)
    expect_true(all(out$start==1L))
    expect_identical(out$end, ends)

    # Final width is enforced.    
    out <- csaw:::.coerceFragments(starts, ends, final=50, chrlen=1e6)
    ref <- resize(IRanges(starts, ends), width=50, fix="center")
    expect_identical(out$start, pmax(1L, start(ref)))
    expect_identical(out$end, end(ref))

    out <- csaw:::.coerceFragments(starts, ends, final=20, chrlen=500)
    ref <- resize(IRanges(starts, ends), width=20, fix="center")
    expect_identical(out$start, pmin(pmax(1L, start(ref)), 500L))
    expect_identical(out$end, pmin(end(ref), 500L))
})

######################################################################

