# This script tests the windowCounts function.
# library(csaw); library(testthat); source("setup.R"); source("test-window.R")

expected_ranges <- function(width, shift, spacing, bam.files, param) {
    spacing <- as.integer(spacing)
    width <- as.integer(width)
    shift <- as.integer(shift)

    chrs <- scanBamHeader(bam.files[1])[[1]][[1]]
    if (length(param$restrict)) { 
        chrs <- chrs[names(chrs) %in% param$restrict] 
    }

    output <- vector("list", length(chrs))
    names(output) <- names(chrs)
    for (x in names(chrs)) {
        multiples <- ceiling(chrs[[x]]/spacing)
        all.starts <- 0:multiples*spacing+1L-shift
        all.ends <- all.starts+width-1L
        all.starts <- pmax(all.starts, 1L)
        all.ends <- pmin(all.ends, chrs[[x]])

        keep <- all.starts <= chrs[[x]] & all.ends > 0L
        gr <- GRanges(x, IRanges(all.starts[keep], all.ends[keep]))
        keep <- !GenomicRanges::duplicated(gr) 
        output[[x]] <- gr[keep]
    }

    output <- suppressWarnings(do.call(c, unname(output)))
    if (!is.na(param$forward)) {
        strand(output) <- ifelse(param$forward, "+", "-")
    }
    return(output)
}

CHECKFUN <- function(bam.files, param, fraglen=200, width=100, shift=0, spacing=50, filter=10) { 
    x <- windowCounts(bam.files, param=param, ext=fraglen, width=width, shift=shift, spacing=spacing, filter=-1)

    # Checking that we get the expected ranges.
    expected <- expected_ranges(width=width, shift=shift, spacing=spacing, bam.files=bam.files, param=param)
    expect_true(all(rowRanges(x)==expected))

    # Checking that we get the expected counts.
    y <- regionCounts(bam.files, param=param, ext=fraglen, regions=rowRanges(x))
    expect_identical(assay(y), assay(x))
    expect_identical(y$totals, x$totals)

    # Checking that filtering behaves as expected.
    z <- windowCounts(bam.files, param=param, ext=fraglen, width=width, shift=shift, spacing=spacing, filter=filter)
    keep <- rowSums(assay(x)) >= filter
    expect_equal(z, x[keep,])

    return(x)
}

tempdir <- tempfile()
dir.create(tempdir)
chromos<-c(chrA=1000, chrB=5000)

set.seed(40001)
test_that("windowCounts works with SE data", {
    for (rparam in list(readParam(), 
                        readParam(minq=10),
                        readParam(dedup=TRUE),
                        readParam(forward=FALSE),
                        readParam(forward=TRUE),
                        readParam(discard=makeDiscard(10, 100, chromos))
                        )) {
        bam.files <-c(regenSE(1000, chromos, file.path(tempdir, "A")), 
                      regenSE(2000, chromos, file.path(tempdir, "B")))
        CHECKFUN(bam.files, param=rparam)
    }

    # Checking extension strategies.
    rparam <- readParam(minq=10)
    for (ext.param in list(NA, 11, 111, 1111)) {
        bam.files <-c(regenSE(1000, chromos, file.path(tempdir, "A")), 
                      regenSE(1500, chromos, file.path(tempdir, "B")))
        CHECKFUN(bam.files, param=rparam, fraglen=ext.param)
    }

    # Checking restriction strategies.
    bam.files <-c(regenSE(600, chromos, file.path(tempdir, "A")), 
                  regenSE(500, chromos, file.path(tempdir, "B")))
    y1 <- CHECKFUN(bam.files, param=readParam(), fraglen=200, width=100, spacing=50, shift=0)
    y2 <- windowCounts(bam.files, param=readParam(restrict=names(chromos)), ext=200, width=100, spacing=50, shift=0)
    expect_identical(assay(y1), assay(y2))
    expect_identical(y1$totals, y2$totals)

    yA <- windowCounts(bam.files, param=readParam(restrict="chrA"), ext=200, width=100, spacing=50, shift=0)
    expect_true(all(seqnames(rowRanges(yA))=="chrA"))

    yB <- windowCounts(bam.files, param=readParam(restrict="chrB"), ext=200, width=100, spacing=50, shift=0)
    expect_true(all(seqnames(rowRanges(yB))=="chrB"))

    expect_equivalent(rbind(assay(yA), assay(yB)), assay(y1))
    expect_identical(suppressWarnings(c(rowRanges(yA), rowRanges(yB))), rowRanges(y1))
    expect_equal(yA$totals + yB$totals, y1$totals)
})

set.seed(40002)
test_that("windowCounts works with PE data", {
    for (rparam in list(readParam(pe="both"), 
                        readParam(pe="both", max.frag=1e8),
                        readParam(pe="both", max.frag=200))) {
        bam.files <-c(regenPE(10000, chromos, file.path(tempdir, "A")), 
                      regenPE(11000, chromos, file.path(tempdir, "B")))
        CHECKFUN(bam.files, param=rparam)
    }

    # Checking alternative counting from PE data.
    for (rparam in list(readParam(pe="first"), 
                        readParam(pe="second"))) {
        bam.files <-c(regenPE(20000, chromos, file.path(tempdir, "A")), 
                      regenPE(10000, chromos, file.path(tempdir, "B")))
        CHECKFUN(bam.files, param=rparam, fraglen=100)
        CHECKFUN(bam.files, param=rparam, fraglen=NA)
    }
})

set.seed(400021)
test_that("windowCounts works with irregular width, shift and spacing", {
    rparam <- readParam(minq=10)

    for (width in list(11, 1111)) {
        bam.files <-c(regenSE(1000, chromos, file.path(tempdir, "A")), 
                      regenSE(1500, chromos, file.path(tempdir, "B")))
        CHECKFUN(bam.files, param=rparam, width=width)
    }

    for (shift in list(11, 22)) {
        bam.files <-c(regenSE(1000, chromos, file.path(tempdir, "A")), 
                      regenSE(1500, chromos, file.path(tempdir, "B")))
        CHECKFUN(bam.files, param=rparam, shift=shift)
    }

    for (spacing in list(22, 111)) {
        bam.files <-c(regenSE(1000, chromos, file.path(tempdir, "A")), 
                      regenSE(1500, chromos, file.path(tempdir, "B")))
        CHECKFUN(bam.files, param=rparam, spacing=spacing)
    }

    # Checking that spacing, shift and width interact properly.
    bam.files <-c(regenSE(1000, chromos, file.path(tempdir, "A")), 
                  regenSE(1500, chromos, file.path(tempdir, "B")))
    CHECKFUN(bam.files, param=rparam, shift=49, spacing=50, width=100)
    CHECKFUN(bam.files, param=rparam, shift=50, spacing=51, width=100)
    CHECKFUN(bam.files, param=rparam, shift=49, spacing=50, width=50)
})

set.seed(40003)
test_that("windowCounts works with binning", {
    # SE case:        
    bam.files <-c(regenSE(1000, chromos, file.path(tempdir, "A")), 
                  regenSE(2000, chromos, file.path(tempdir, "B")))
    x <- windowCounts(bam.files, bin=TRUE, width=100)
    y <- windowCounts(bam.files, bin=FALSE, width=100, ext=1, spacing=100, shift=0, filter=1)
    expect_identical(assay(x), assay(y))
    expect_equal(colSums(assay(x)), x$totals)

    z <- windowCounts(bam.files, bin=TRUE, width=max(chromos)) # expect one range per chromosome.
    expect_equal(colSums(assay(z)), z$totals)
    expect_true(all(rowRanges(z)==GRanges(names(chromos), IRanges(1, chromos)))) 

    # PE case (harder to compare exactly, as the midpoints of the fragments are used):
    bam.files <-c(regenPE(10000, chromos, file.path(tempdir, "A")), 
                  regenPE(20000, chromos, file.path(tempdir, "B")))
    rparam <- readParam(pe="both")
    x <- windowCounts(bam.files, bin=TRUE, width=100, param=rparam, filter=0)
    expect_equal(colSums(assay(x)), x$totals)

    y <- windowCounts(bam.files, bin=TRUE, width=200, param=rparam)
    nesting <- findOverlaps(rowRanges(x), rowRanges(y), type="within") 
    expect_identical(anyDuplicated(queryHits(nesting)), 0L) # fully nested.
    aggregator <- aggregate(assay(x) ~ subjectHits(nesting), FUN=sum)
    summed.y <- as.matrix(aggregator[,2:3])
    expect_equivalent(summed.y, assay(y))

    # PE case with very small chromosomes, check for correct end-of-chromosome behaviour.
    bam.files <-c(regenPE(10000, chromos/10, file.path(tempdir, "A")), 
                  regenPE(20000, chromos/10, file.path(tempdir, "B")))
    x <- windowCounts(bam.files, bin=TRUE, width=100, param=rparam, filter=0)
    expect_equal(colSums(assay(x)), x$totals)
})

test_that("windowCounts gives the same result for BamFiles", {
    bam.files <-c(regenSE(1000, chromos, file.path(tempdir, "A")), 
                  regenSE(2000, chromos, file.path(tempdir, "B")))

    x1 <- windowCounts(bam.files, bin=TRUE, width=100)
    x2 <- windowCounts(lapply(bam.files, BamFile), bin=TRUE, width=100)
    expect_identical(x1, x2)
})

set.seed(40004)
test_that("windowCounts behaves correctly with silly inputs", {
    # Doesn't fail with no counts.
    obam <- regenSE(0L, chromos, outfname=tempfile())
    out <- windowCounts(obam)
    expect_true(all(assay(out)==0))
    expect_identical(out$totals, 0L)

    # Throws with invalid arguments.
    expect_error(windowCounts(obam, shift=-1), "must be a non-negative integer")
    expect_error(windowCounts(obam, spacing=0), "must be a positive integer")
    expect_error(windowCounts(obam, width=0), "must be a positive integer")
    expect_error(windowCounts(obam, shift=50, spacing=50), "shift must be less than the spacing")
})
