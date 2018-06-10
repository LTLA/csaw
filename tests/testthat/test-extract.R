# This tests the extractReads function.
# library(csaw); library(testthat); source("test-extract.R")

source("simsam.R")

set.seed(2000)
test_that("extractReads works correctly without extension", {
    chromos <- c(chrA=1000, chrB=3000, chrC=2000)
    discarder <- makeDiscard(10, 100, chromos)
    genome <- GRanges(names(chromos), IRanges(1L, chromos))
    nreads <- 1000L

    for (rparam in list(readParam(),
                        readParam(minq=10),
                        readParam(dedup=TRUE),
                        readParam(forward=FALSE),
                        readParam(forward=TRUE),
                        readParam(discard=discarder),
                        readParam(minq=10, dedup=TRUE, discard=discarder)
                        )) {
        obam <- regenSE(nreads, chromos, outfname=tempfile())
        extractor <- generateWindows(chromos, 1, 500)  # one per chromosome, in order.

        for (idx in seq_along(extractor)) {
            element <- extractor[idx]
            out <- extractReads(obam, element, param=rparam)
            expect_true(all(seqnames(out)==names(chromos)[idx]))

            ref <- csaw:::.extractSE(obam, genome[idx], rparam)
            ref.forward <- IRanges(ref$forward$pos, pmin(ref$forward$pos + ref$forward$qwidth - 1L, chromos[idx]))
            ref.forward <- ref.forward[overlapsAny(ref.forward, ranges(element))]
            ref.reverse <- IRanges(ref$reverse$pos, pmin(ref$reverse$pos + ref$reverse$qwidth - 1L, chromos[idx]))
            ref.reverse <- ref.reverse[overlapsAny(ref.reverse, ranges(element))]

            forward <- out[strand(out)=="+"]
            expect_identical(start(forward), start(ref.forward))
            expect_identical(width(forward), width(ref.forward))

            reverse <- out[strand(out)=="-"]
            expect_identical(start(reverse), start(ref.reverse))
            expect_identical(width(reverse), width(ref.reverse))
        }
    }

    # Testing what happens with no reads.
    obam <- regenSE(0L, chromos, outfname=tempfile())
    extractor <- generateWindows(chromos, 1, 500)
    for (idx in seq_along(genome)) {
        element <- extractor[idx]
        out <- extractReads(obam, element, param=readParam())
        expect_identical(length(out), 0L)
    }
})

set.seed(2001)
test_that("extractReads works correctly with extension", {
    chromos <- c(chrA=1000, chrB=3000, chrC=2000)
    discarder <- makeDiscard(10, 100, chromos)
    genome <- GRanges(names(chromos), IRanges(1L, chromos))
    nreads <- 100L

    for (ext.strat in c(20L, 50L, 100L, 200L)) {
        obam <- regenSE(nreads, chromos, outfname=tempfile())
        extractor <- generateWindows(chromos, 1, 500)

        for (idx in seq_along(extractor)) {
            element <- extractor[idx]
            out <- extractReads(obam, element, ext=ext.strat, param=readParam())

            ref <- extractReads(obam, genome[idx], ext=NA, param=readParam())
            suppressWarnings(ref <- resize(ref, fix="start", width=ext.strat, ignore.strand=FALSE))
            ref <- trim(ref)
            ref <- ref[overlapsAny(ref, element)]

            expect_identical(ref, out)
        }
    }
})

set.seed(2002)
test_that("extractReads works correctly in the paired-end case", {
    chromos <- c(chrA=1000, chrB=3000)
    genome <- GRanges(names(chromos), IRanges(1L, chromos))
    discarder <- makeDiscard(10, 100, chromos)
    npairs <- 10000

    for (rparam in list(readParam(pe="both"),
                        readParam(pe="both", minq=10),
                        readParam(pe="both", dedup=TRUE),
                        readParam(pe="both", discard=discarder),
                        readParam(pe="both", max.frag=200),
                        readParam(pe="both", max.frag=1e8),
                        readParam(pe="both", minq=10, dedup=TRUE, discard=discarder)
                        )) {

        obam <- regenPE(npairs, chromos, outfname=tempfile())
        extractor <- generateWindows(chromos, 1, 500)

        for (idx in seq_along(genome)) {
            element <- extractor[idx]
            out <- extractReads(obam, element, param=rparam)
            expect_true(all(seqnames(out)==names(chromos)[idx]))

            ref <- csaw:::.extractPE(obam, genome[idx], rparam)
            ref <- IRanges(ref$pos, pmin(ref$pos + ref$size - 1L, chromos[idx]))
            ref <- ref[overlapsAny(ref, ranges(element))]
            expect_identical(ref, ranges(out))
        }
    }

    # Testing what happens with no reads.
    obam <- regenPE(0L, chromos, outfname=tempfile())
    for (idx in seq_along(genome)) {
        element <- extractor[idx]
        out <- extractReads(obam, element, param=rparam)
        expect_identical(length(out), 0L)
    }
})

set.seed(2003)
test_that("extractReads works correctly to return the reads themselves", {
    chromos <- c(chrA=1000, chrB=3000, chrC=2000)
    discarder <- makeDiscard(10, 100, chromos)
    genome <- GRanges(names(chromos), IRanges(1L, chromos))
    rparam <- readParam(pe="both")

    npairs <- 10000L
    obam <- regenPE(npairs, chromos, outfname=tempfile())
    extractor <- generateWindows(chromos, 1, 500)
    for (idx in seq_along(genome)) {
        element <- extractor[idx]
        ref <- extractReads(obam, element, param=rparam)
        out <- extractReads(obam, element, param=rparam, as.reads=TRUE)
        expect_identical(start(ref), start(out$forward))
        expect_identical(end(ref), end(out$reverse))
    }

    # Directly returning reads in pe="first" or "second" mode.
    obam <- regenPE(npairs, chromos, outfname=tempfile())
    extractor <- generateWindows(chromos, 1, 500)
    for (idx in seq_along(genome)) {
        ref <- extractReads(obam, extractor[idx], param=readParam())
        out1 <- extractReads(obam, extractor[idx], param=reform(rparam, pe="first"))
        out2 <- extractReads(obam, extractor[idx], param=reform(rparam, pe="second"))
        combined <- c(out1, out2)
        expect_identical(sort(ref), sort(combined))
    }
})