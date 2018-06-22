# This tests the basic internal output functions.
# library(testthat); library(csaw); source("test-int_output.R")

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
