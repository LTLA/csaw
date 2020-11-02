# This tests the correctness of the normOffsets functions.
# library(csaw); library(testthat); source("test-norm.R")

set.seed(1000)
data <- SummarizedExperiment(list(counts=matrix(rpois(10000, lambda=10), ncol=10)))
data$totals <- rpois(10, lambda=10000)
library(edgeR)
ref <- calcNormFactors(DGEList(assay(data), lib.size=data$totals), doWeighting=FALSE, method="TMM")$samples$norm.factors

test_that("normFactors works as expected", {
    nf <- normFactors(data, se.out=FALSE)
    expect_identical(nf, ref)
    
    data2 <- normFactors(data, se.out=TRUE)
    expect_identical(data2$norm.factors, ref)
    expect_identical(data2$totals, data$totals)

    # Works with DGEList inputs.
    nf2 <- normFactors(asDGEList(data), se.out=FALSE)
    expect_equivalent(nf, nf2)

    y <- normFactors(asDGEList(data))
    expect_equivalent(y$samples$norm.factors, nf2)

    mat <- normFactors(assay(data))
    expect_type(attr(mat, "norm.factors"), "double")

    # Behaves when empty.
    # 1 is correct, as calcNormFactors() just diverts to that.
    expect_identical(normFactors(data[0,], se.out=FALSE), rep(1, 10)) 
    expect_identical(normFactors(data[,0], se.out=FALSE), numeric(0))
})

test_that("normFactors works with another se.out", {
    data3 <- data

    # Checking that overwriting se.out works.
    assay(data3) <- assay(data3)*runif(nrow(data3))
    data3b <- normFactors(data, se.out=data3)
    expect_identical(assay(data3b), assay(data3))
    expect_identical(data3b$norm.factors, ref)

    # Works with an output DGEList.
    data3c <- normFactors(data, se.out=asDGEList(data3))
    expect_equivalent(data3b$norm.factors, data3c$samples$norm.factors)

    # Throwing an error when the library sizes are different.
    data3$totals <- rpois(10, lambda=10000)
    expect_error(normFactors(data, se.out=data3), "library sizes")
    expect_error(normFactors(data, se.out=data[,1]), "library sizes")
})

set.seed(1001)
test_that("normFactors yields the correct output from undersampling", {
    n <- 1000L
    mu1 <- rep(10, n)
    mu2 <- mu1
    mu2[1:100] <- 100
    undersamp <- sum(mu1)/sum(mu2)
    mu2 <- mu2*undersamp
    counts <- cbind(rnbinom(n, mu=mu1, size=20), rnbinom(n, mu=mu2, size=20))
    
    data <- SummarizedExperiment(list(counts=counts))
    data$totals <- c(sum(mu1), sum(mu2))
    truth <- sqrt(c(1/undersamp, undersamp)) 
    expect_true(all(abs(normFactors(data)$norm.factors - truth) < 0.1))

    # Testing what happens with weighting, after adding some high-abundance DB regions. 
    # The errors should be larger than the unweighted version.
    n <- 100000
    blah <- matrix(rnbinom(2*n, mu=10, size=20), ncol=2)
    tospike <- 10000
    blah[1:tospike,1] <- rnbinom(tospike, mu=1000, size=20)
    blah[1:tospike,2] <- rnbinom(tospike, mu=2000, size=20)
    full.lib.size <- colSums(blah)
    
    true.value <- 1/full.lib.size
    true.value <- true.value/exp(mean(log(true.value)))
    
    data <- SummarizedExperiment(list(counts=blah))
    data$totals <- full.lib.size
    expect_true(all(abs(normFactors(data)$norm.factors - true.value) < 
        abs(normFactors(data, weighted=TRUE)$norm.factors - true.value)))
})

set.seed(1000)
means <- 2^runif(1000)
data <- SummarizedExperiment(list(counts=matrix(rpois(10000, lambda=means), ncol=10)))
data$totals <- rpois(10, lambda=10000)

test_that("normOffsets works correctly", {
    offs <- normOffsets(data, se.out=FALSE)
    data <- normOffsets(data, se.out=TRUE)
    expect_equal(offs, assay(data, "offset"))

    # Behaves with a DGEList input.
    offs2 <- normOffsets(asDGEList(data))
    expect_equivalent(offs2$offset, offs)

    mat <- normOffsets(assay(data))
    expect_equal(dim(attr(mat, "offset")), dim(data))

    # Behaves when empty.
    expect_identical(dim(normOffsets(data[0,], se.out=FALSE)), c(0L, 10L))
})

test_that("normOffsets works with another se.out", {
    data <- normOffsets(data)

    # Checking that overwriting se.out works.
    shuffler <- sample(nrow(data))
    data2 <- data[shuffler,]
    data2 <- normOffsets(data, se.out=data2)
    expect_equal(assay(data2, "offset"), assay(data, "offset")[shuffler,])
 
    # Reference calculation, after subtracting the reference 'ab' from the observed values.
    # This represents the theoretical fit to the M-values computed against the reference.
    lib.sizes <- data$totals
    mat <- assay(data)
    
    cont.cor <- 0.5
    cont.cor.scaled <- cont.cor * lib.sizes/mean(lib.sizes)
    ab <- aveLogCPM(mat, lib.size=lib.sizes, prior.count=cont.cor)/log2(exp(1))
    
    ref <- matrix(0, nrow(mat), ncol(mat), byrow=TRUE, dimnames=dimnames(mat))
    for (x in seq_len(ncol(mat))) {
        fit <- loessFit(log(mat[,x]+cont.cor.scaled[x]) - ab, ab) # explicit subtraction this time.
        ref[,x] <- fit$fitted 
    }
    ref <- scaleOffset(lib.sizes, ref)

    expect_equal(ref, assay(data, "offset"))

    # Breaks when the library sizes are different.
    data3 <- data
    data3$totals <- rpois(10, lambda=10000)
    expect_error(normOffsets(data, se.out=data3), "library sizes")
    expect_error(normOffsets(data, se.out=data[,1]), "library sizes")
})
