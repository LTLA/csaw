# Various tests for expression-computing functions.
# library(csaw); library(testthat); source("test-exprs.R")

set.seed(100)
test_that("scaledAverage works as expected", {
    se <- SummarizedExperiment(list(counts=matrix(rnbinom(1000, mu=100, size=10), ncol=10, nrow=100)))
    se$totals <- runif(ncol(se), 1e6, 2e6)

    y <- asDGEList(se)
    stopifnot(identical(y$samples$lib.size, se$totals))
    ref <- aveLogCPM(y)
    out <- scaledAverage(se, scale=1)
    expect_equal(ref, out)

    se$norm.factors <- runif(10, 0.5, 1.5)
    y <- asDGEList(se)
    stopifnot(identical(y$samples$norm.factors, se$norm.factors))
    ref <- aveLogCPM(y)
    out <- scaledAverage(se, scale=1)
    expect_equal(ref, out)
    
    assay(se, "offset") <- matrix(rnorm(1000), ncol=10, nrow=100) # Neither function should care about the offset.
    y <- asDGEList(se)
    stopifnot(!is.null(y$offset))    
    ref <- aveLogCPM(y)
    out <- scaledAverage(se, scale=1)
    expect_equal(ref, out)
    
    ref <- aveLogCPM(y, dispersion=0.1)
    out <- scaledAverage(se, scale=1, dispersion=0.1)
    expect_equal(ref, out)
})

test_that("scaledAverage works correctly with scaling", {
    se <- SummarizedExperiment(list(counts=matrix(rnbinom(100, mu=100, size=10), ncol=1, nrow=100)))
    se$totals <- runif(ncol(se), 1e6, 2e6)
    ref <- scaledAverage(se, scale=1)
   
    # Downscaling after multiplication of counts should yield the same result 
    # (for a single library, at least!)
    se2 <- se
    assay(se2) <- assay(se) * 2
    out2 <- scaledAverage(se2, scale=2)
    expect_equal(ref, out2)
    
    scalar <- runif(nrow(se), 1, 5)
    se1 <- se
    assay(se1) <- scalar * assay(se)
    out1 <- scaledAverage(se1, scale=scalar)
    expect_equal(ref, out1)
    
    # Checking proper behaviour with invalid inputs.
    stopifnot(all(is.na(scaledAverage(se, scale=-1))))
    stopifnot(all(is.infinite(scaledAverage(se, scale=0))))
    flucscale <- rep(1, nrow(y))
    flucscale[1] <- 0
    flucscale[2] <- -1
    
    out <- scaledAverage(se, scale=flucscale)
    ref <- aveLogCPM(asDGEList(se))
    expect_true(is.infinite(out[1]))
    expect_true(is.na(out[2]))
    expect_equal(out[-(1:2)], ref[-(1:2)])
})

