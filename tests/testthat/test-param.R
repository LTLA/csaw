# This tests functions related to the readParam() class.
# library(testthat); library(csaw); source("test-param.R")

test_that("readParam constructor works as expected", {
    expect_s4_class(X <- readParam(), "readParam")
    expect_identical(X$dedup, FALSE)
    expect_identical(X$max.frag, 500L)
    expect_identical(X$minq, NA_integer_)
    expect_identical(X$forward, NA)
    expect_identical(X$restrict, character(0))
    
    expect_s4_class(X <- readParam(dedup=TRUE, max.frag=200L, minq=10, forward=TRUE, restrict="chrB"), "readParam")
    expect_identical(X$dedup, TRUE)
    expect_identical(X$max.frag, 200L)
    expect_identical(X$minq, 10L)
    expect_identical(X$forward, TRUE)
    expect_identical(X$restrict, "chrB")
})

test_that("readParam reform works as expected", {
    X <- readParam()
    Y <- reform(X, minq=10)
    expect_identical(Y$minq, 10L)

    Y <- reform(X, dedup=TRUE)
    expect_identical(Y$dedup, TRUE)

    Y <- reform(X, max.frag=250L)
    expect_identical(Y$max.frag, 250L)

    Y <- reform(X, forward=FALSE)
    expect_identical(Y$forward, FALSE)

    Y <- reform(X, restrict=c("chrA", "chrB"))
    expect_identical(Y$restrict, c("chrA", "chrB"))
})

test_that("readParam discard setup works as expected", {
    gr <- GRanges(c("chrA", "chrB"), IRanges(1:2, 10:11))

    X <- readParam(discard=gr)
    expect_identical(X$processed.discard$chrA$pos, c(1L, 11L))
    expect_identical(X$processed.discard$chrA$id, c(0L, 0L))
    expect_identical(X$processed.discard$chrB$pos, c(2L, 12L))
    expect_identical(X$processed.discard$chrB$id, c(0L, 0L))

    Y <- readParam()
    expect_identical(unname(Y$processed.discard), list())
    Y <- reform(Y, discard=gr)
    expect_identical(X$processed.discard, Y$processed.discard)
})
