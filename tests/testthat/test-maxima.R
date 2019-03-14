# This script tests the findMaxima function.
# library(testthat); library(csaw); source("test-maxima.R")

set.seed(7000)
test_that("findMaxima works correctly", {
    chromos=c(chrA=1000, chrB=2000, chrC=500)
    for (nregs in c(20, 100, 500)) {
        for (winsize in c(10, 100)) {
            for (range in c(50, 100, 200)) {

            	reg.data <- generateWindows(chromos, nregs, winsize)
                reg.data <- reg.data[sample(length(reg.data))]
                metric <- runif(length(reg.data), -1, 1)
                
                # Finding our own maxima.
                new.reg <- reg.data
                start(new.reg) <- start(new.reg) - range
                end(new.reg) <- end(new.reg) + range
                olap <- findOverlaps(reg.data, new.reg)
                above <- metric[queryHits(olap)] >= metric[subjectHits(olap)]
                check.max <- unname(vapply(split(above, queryHits(olap)), all, FUN.VALUE=TRUE))

                # Getting max.
                is.max <- findMaxima(reg.data, range=range, metric=metric)
                expect_identical(is.max, check.max)
            }
        }
    }
})

set.seed(7001)
test_that("findMaxima is correctly strand-aware", {
    chromos=c(chrA=1020, chrB=2001, chrC=2500)
    for (nregs in c(20, 100, 500)) {
        for (winsize in c(10, 100)) {
            for (range in c(50, 100, 200)) {
            
                reg <- generateWindows(chromos, nregs, winsize)
                strand(reg) <- sample(c("+", "-", "*"), length(reg), replace=TRUE)
                metric <- rnorm(length(reg))
                combo <- findMaxima(reg, range=range, metric=metric, ignore.strand=FALSE)

                # Running separately on each strand, and checking that the boundaries are the same.
                ref <- logical(length(combo))
                is.forward <- as.logical(strand(reg)=="+")
                ref[is.forward] <- findMaxima(reg[is.forward], range=range, metric=metric[is.forward])
                is.reverse <- as.logical(strand(reg)=="-")
                ref[is.reverse] <- findMaxima(reg[is.reverse], range=range, metric=metric[is.reverse])
                is.unstrand <- as.logical(strand(reg)=="*")
                ref[is.unstrand] <- findMaxima(reg[is.unstrand], range=range, metric=metric[is.unstrand])

                expect_identical(combo, ref)
            }
        }
    }
})

set.seed(7002)
test_that("findMaxima behaves with silly inputs", {
    chromos=c(chrA=100, chrB=200, chrC=500)
    reg <- generateWindows(chromos, nwin=10, winsize=10)
    expect_identical(findMaxima(reg[0], range=10, metric=numeric(0)), logical(0))
    expect_error(findMaxima(reg, metric=numeric(0), range=10), "per region")
    expect_error(findMaxima(reg, metric=rep(NA, length(reg)), range=10), "missing values")
})
