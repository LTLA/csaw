# This tests the various wrapper functions.
# library(testthat); library(csaw); source("setup.R"); source("test-wrappers.R")

set.seed(1919000)
test_that("mergeResults works correctly", {
    chromos <- c(chrA=10000, chrB=5000, chrC=2000)
    reg <- generateWindows(chromos, nwin=100, winsize=50)
    ns <- length(reg)
    tab <- data.frame(logFC=rnorm(ns), PValue=rbeta(ns, 1, 3), logCPM=rnorm(ns))

    # Yields all output fields, correctly named.
    merged <- mergeResults(reg, tab, tol=100)
    expect_s4_class(merged$regions, "GRanges")
    expect_s4_class(merged$combined, "DataFrame")
    expect_s4_class(merged$best, "DataFrame")
    expect_type(metadata(merged)$ids, "integer")

    # Same results in the metadata.
    reg2 <- reg
    mcols(reg2) <- tab
    expect_identical(merged, mergeResults(reg2, tol=100))

    # Responds to other arguments
    merged2 <- mergeResults(reg, tab, tol=1)
    expect_false(identical(merged, merged2))

    merged2 <- mergeResults(reg, tab, tol=100, merge.args=list(max.width=200))
    expect_false(identical(merged, merged2))

    merged2 <- mergeResults(reg, tab, tol=100, combine.args=list(fc.col="logCPM"))
    expect_false(identical(merged, merged2))

    merged2 <- mergeResults(reg, tab, tol=100, get.best=FALSE)
    expect_false(identical(merged, merged2))

    # Behaves with empty arguments.
    expect_identical(nrow(mergeResults(reg[0], tab[0,], tol=10)), 0L)
})

set.seed(1919001)
test_that("overlapResults works correctly", {
    chromos <- c(chrA=10000, chrB=5000, chrC=2000)
    reg <- generateWindows(chromos, nwin=100, winsize=50)
    ref <- generateWindows(chromos, nwin=5, winsize=50)
    ns <- length(reg)
    tab <- data.frame(logFC=rnorm(ns), PValue=rbeta(ns, 1, 3), logCPM=rnorm(ns))

    # Yields all output fields, correctly named.
    merged <- overlapResults(reg, tab, regions=ref)
    expect_identical(nrow(merged), length(ref))

    expect_s4_class(merged$regions, "GRanges")
    expect_s4_class(merged$combined, "DataFrame")
    expect_s4_class(merged$best, "DataFrame")
    expect_s4_class(metadata(merged)$overlaps, "Hits")

    # Same results in the metadata.
    reg2 <- reg
    mcols(reg2) <- tab
    expect_identical(merged, overlapResults(reg2, regions=ref))

    # Responds to other arguments
    merged2 <- overlapResults(reg, tab, regions=ref, overlap.args=list(type="within"))
    expect_false(identical(merged, merged2))

    merged2 <- overlapResults(reg, tab, regions=ref, combine.args=list(fc.col="logCPM"))
    expect_false(identical(merged, merged2))

    merged2 <- overlapResults(reg, tab, regions=ref, get.best=FALSE)
    expect_false(identical(merged, merged2))

    # Behaves with empty arguments.
    expect_identical(nrow(overlapResults(reg[0], tab[0,], ref)), length(ref))
})

set.seed(1919002)
test_that("mergeResultsList works correctly", {
    data.list <- result.list <- list()
    chromos <- c(chrA=10000, chrB=5000, chrC=2000)
    sizes <- c(10, 20, 50)
    for (s in seq_along(sizes)) {
        data.list[[s]] <- generateWindows(chromos, winsize=sizes[s], nwin=20)
        ns <- length(data.list[[s]])
        result.list[[s]] <- data.frame(logFC=rnorm(ns), PValue=rbeta(ns, 1, 3), logCPM=rnorm(ns))
    }

    merged <- mergeResultsList(data.list, result.list, tol=100)
    expect_s4_class(merged$regions, "GRanges")
    expect_s4_class(merged$combined, "DataFrame")
    expect_s4_class(merged$best, "DataFrame")

    expect_type(metadata(merged)$ids, "integer")
    expect_type(metadata(merged)$weights, "double")
    expect_identical(length(metadata(merged)$ranges), sum(lengths(data.list)))
    expect_identical(metadata(merged)$tab, do.call(rbind, result.list))

    # Same results in the metadata.
    data.list2 <- data.list
    for (s in seq_along(data.list2)) {
        mcols(data.list2[[s]]) <- result.list[[s]]
    }
    merged2 <- mergeResultsList(data.list2, tol=100)
    stripped <- merged
    metadata(stripped) <- metadata(merged2) <- list()
    expect_identical(stripped, merged2)

    # Responds to other arguments.
    merged2 <- mergeResultsList(data.list, result.list, tol=1)
    expect_false(identical(merged, merged2))

    merged2 <- mergeResultsList(data.list, result.list, tol=100, merge.args=list(max.width=200))
    expect_false(identical(merged, merged2))

    merged2 <- mergeResultsList(data.list, result.list, tol=100, combine.args=list(fc.col="logCPM"))
    expect_false(identical(merged, merged2))

    merged2 <- mergeResultsList(data.list, result.list, tol=100, get.best=FALSE)
    expect_false(identical(merged, merged2))

    merged2 <- mergeResultsList(data.list, result.list, tol=100, equiweight=FALSE)
    expect_false(identical(merged, merged2))

    # Behaves with single arguments.
    merged3 <- mergeResultsList(data.list[1], result.list[1], tol=100)
    refout <- mergeResults(data.list[[1]], result.list[[1]], tol=100)
    expect_equal(merged3$combined, refout$combined)
    expect_equal(merged3$best, refout$best)
    expect_true(all(merged3$regions==refout$regions))
})

set.seed(1919003)
test_that("overlapResultsList works correctly", {
    data.list <- result.list <- list()
    chromos <- c(chrA=10000, chrB=5000, chrC=2000)
    ref <- generateWindows(chromos, nwin=5, winsize=50)

    sizes <- c(10, 20, 50)
    for (s in seq_along(sizes)) {
        data.list[[s]] <- generateWindows(chromos, winsize=sizes[s], nwin=20)
        ns <- length(data.list[[s]])
        result.list[[s]] <- data.frame(logFC=rnorm(ns), PValue=rbeta(ns, 1, 3), logCPM=rnorm(ns))
    }

    merged <- overlapResultsList(data.list, result.list, region=ref)
    expect_s4_class(merged$regions, "GRanges")
    expect_s4_class(merged$combined, "DataFrame")
    expect_s4_class(merged$best, "DataFrame")

    expect_type(metadata(merged)$weights, "double")
    expect_identical(length(metadata(merged)$ranges), sum(lengths(data.list)))
    expect_identical(metadata(merged)$tab, do.call(rbind, result.list))

    # Same results in the metadata.
    data.list2 <- data.list
    for (s in seq_along(data.list2)) {
        mcols(data.list2[[s]]) <- result.list[[s]]
    }
    merged2 <- overlapResultsList(data.list2, region=ref)

    stripped <- merged
    metadata(stripped) <- metadata(merged2) <- list()
    expect_identical(stripped, merged2)

    # Responds to other arguments.
    merged2 <- overlapResultsList(data.list, result.list, regions=ref, overlap.args=list(type="within"))
    expect_false(identical(merged, merged2))

    merged2 <- overlapResultsList(data.list, result.list, regions=ref, combine.args=list(fc.col="logCPM"))
    expect_false(identical(merged, merged2))

    merged2 <- overlapResultsList(data.list, result.list, regions=ref, get.best=FALSE)
    expect_false(identical(merged, merged2))

    merged2 <- overlapResultsList(data.list, result.list, regions=ref, equiweight=FALSE)
    expect_false(identical(merged, merged2))

    # Behaves with single arguments.
    merged3 <- overlapResultsList(data.list[1], result.list[1], regions=ref)
    refout <- overlapResults(data.list[[1]], result.list[[1]], regions=ref)
    expect_equal(merged3$combined, refout$combined)
    expect_equal(merged3$best, refout$best)
    expect_true(all(merged3$regions==refout$regions))
})
