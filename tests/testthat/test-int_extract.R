# This tests the basic internal read extraction functions.
# library(testthat); library(csaw); source("test-int_extract.R")

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

    # Strand should be specified.
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

    # Strand specification will always fail in PE mode.
    expect_error(readParam(pe="both", forward=NULL), "strand-specific extraction")
    expect_error(readParam(pe="both", forward=TRUE), "strand-specific extraction")
    expect_error(readParam(pe="both", forward=FALSE), "strand-specific extraction")
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
