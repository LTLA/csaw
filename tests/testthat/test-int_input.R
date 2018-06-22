# This tests the basic internal functions on which all others are built.
# library(testthat); library(csaw); source("test-int_input.R")

source("simsam.R")

set.seed(320000)
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

set.seed(320001)
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

