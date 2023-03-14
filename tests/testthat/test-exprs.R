# Various tests for expression-computing functions.
# library(csaw); library(testthat); source("test-exprs.R")

library(edgeR)

set.seed(100)
test_that("asDGEList works as expected", {
    se <- SummarizedExperiment(list(counts=matrix(rnbinom(1000, mu=100, size=10), ncol=10, nrow=100)))
    se$totals <- runif(ncol(se), 1e6, 2e6)

    y <- asDGEList(se)
    expect_equivalent(y$counts, assay(se))
    expect_identical(y$samples$lib.size, se$totals)
    expect_identical(y$samples$norm.factors, rep(1, ncol(se)))
    expect_identical(y$offset, NULL)

    se$norm.factors <- runif(10, 0.5, 1.5)
    y <- asDGEList(se)
    expect_equivalent(y$counts, assay(se))
    expect_identical(y$samples$lib.size, se$totals)
    expect_identical(y$samples$norm.factors, se$norm.factors)
    expect_identical(y$offset, NULL)

    assay(se, "offset") <- matrix(rnorm(1000), ncol=10, nrow=100) 
    y <- asDGEList(se)
    expect_equivalent(y$counts, assay(se))
    expect_identical(y$samples$lib.size, se$totals)
    expect_identical(y$samples$norm.factors, se$norm.factors)
    expect_equivalent(y$offset, scaleOffset(y, assay(se, "offset"))$offset)

#    # Works on an empty DGEList.
#    emp <- SummarizedExperiment(list(counts=matrix(0,0,10)))
#    emp$totals <- 1
#    y2 <- asDGEList(emp)
#    expect_identical(nrow(y2), 0L)
})

set.seed(1000)
test_that("scaledAverage works as expected", {
    se <- SummarizedExperiment(list(counts=matrix(rnbinom(1000, mu=100, size=10), ncol=10, nrow=100)))
    se$totals <- runif(ncol(se), 1e6, 2e6)

    y <- asDGEList(se)
    ref <- aveLogCPM(y)
    out <- scaledAverage(se, scale=1)
    expect_equal(ref, out)

    se$norm.factors <- runif(10, 0.5, 1.5)
    y <- asDGEList(se)
    ref <- aveLogCPM(y)
    out <- scaledAverage(se, scale=1)
    expect_equal(ref, out)
    
    assay(se, "offset") <- matrix(rnorm(1000), ncol=10, nrow=100) # Neither function should care about the offset.
    y <- asDGEList(se)
    ref <- aveLogCPM(y)
    out <- scaledAverage(se, scale=1)
    expect_equal(ref, out)
    
    ref <- aveLogCPM(y, dispersion=0.1)
    out <- scaledAverage(se, scale=1, dispersion=0.1)
    expect_equal(ref, out)

    # Trying with empty inputs.
    emp <- SummarizedExperiment(list(counts=matrix(0,0,10)))
    emp$totals <- 10
    expect_identical(scaledAverage(emp), numeric(0))
    expect_identical(scaledAverage(emp, scale=numeric(0)), numeric(0))
})

set.seed(1001)
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
    flucscale <- rep(1, nrow(se))
    flucscale[1] <- 0
    flucscale[2] <- -1
    
    out <- scaledAverage(se, scale=flucscale)
    ref <- aveLogCPM(asDGEList(se))
    expect_true(is.infinite(out[1]))
    expect_true(is.na(out[2]))
    expect_equal(out[-(1:2)], ref[-(1:2)])
})

set.seed(1002)
test_that("calculateCPM works correctly with library sizes", {
    se <- SummarizedExperiment(list(counts=matrix(rnbinom(1000, mu=100, size=10), ncol=10, nrow=100)))
    se$totals <- runif(ncol(se), 1e6, 2e6)
    
    y <- asDGEList(se)
    ref <- cpm(y, prior.count=2, log=TRUE)
    out <- calculateCPM(se, prior.count=2, log=TRUE)
    dimnames(out) <- dimnames(ref) <- NULL
    expect_equal(ref, out)

    # Without log-transformation.
    ref1 <- cpm(y, prior.count=2, log=FALSE)
    out1 <- calculateCPM(se, prior.count=2, log=FALSE)
    dimnames(ref1) <- dimnames(out1) <- NULL
    expect_equal(ref1, out1)

    # With normalization factors.
    se$norm.factors <- runif(ncol(se))
    y <- asDGEList(se)
    ref2 <- cpm(y, prior.count=2, log=TRUE)
    out2 <- calculateCPM(se, prior.count=2, log=TRUE)
    dimnames(out2) <- dimnames(ref2) <- NULL
    expect_equal(ref2, out2)

    # Disobey the normalization factors.
    out3 <- calculateCPM(se, prior.count=2, log=TRUE, use.norm.factors=FALSE)
    dimnames(out3) <- NULL
    expect_equal(out, out3)
    expect_false(isTRUE(all.equal(out2, out3)))
})

set.seed(1003)
test_that("calculateCPM works correctly with offsets", {
    se <- SummarizedExperiment(list(counts=matrix(rnbinom(1000, mu=100, size=10), ncol=10, nrow=100)))
    se$totals <- runif(ncol(se), 100, 200)
    ref <- calculateCPM(se, prior.count=2, log=FALSE, use.offset=FALSE)
    expect_equivalent(ref, cpm(asDGEList(se)))

    # Responds to the library sizes. 
    centered.off <- log(se$totals)
    centered.off <- centered.off - mean(centered.off)
    assay(se, "offset") <- matrix(centered.off, nrow=nrow(se), ncol=ncol(se), byrow=TRUE)
    out <- calculateCPM(se, prior.count=2, log=FALSE, use.offset=TRUE)
    ref <- calculateCPM(se, prior.count=2, log=FALSE, use.offset=FALSE)
    expect_equal(out, ref)

    # Responds correctly with log-transformation. 
    ref <- calculateCPM(se, prior.count=2, log=TRUE, use.offset=FALSE)
    out <- calculateCPM(se, prior.count=2, log=TRUE, use.offset=TRUE)
    expect_equal(out, ref)

    # Responds to _different_ offsets per gene.
    se1 <- se
    assay(se1, "offset") <- matrix(runif(ncol(se)), nrow=nrow(se), ncol=ncol(se), byrow=TRUE)
    se2 <- se
    assay(se2, "offset") <- matrix(rnorm(ncol(se)), nrow=nrow(se), ncol=ncol(se), byrow=TRUE)
    out <- calculateCPM(rbind(se1, se2))
    ref <- rbind(calculateCPM(se1), calculateCPM(se2))
    expect_equivalent(out, ref)

    # Behaves correctly with no genes.
    expect_equivalent(calculateCPM(se[0,]), matrix(0,0,ncol(se)))
    expect_equivalent(calculateCPM(se[,0]), matrix(0,nrow(se),0))
})

