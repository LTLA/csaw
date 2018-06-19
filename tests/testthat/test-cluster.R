# This tests the post-hoc clustering methods.
# library(csaw); library(testthat); source("test-cluster.R")

compCFDR <- function(ids, threshold, weights) {
    if (is.null(weights)) { 
        obs.sizes <- table(ids) 
    } else {
        obs.sizes <- sapply(split(weights, ids), FUN=sum)
    }
    obs.sizes <- sort(obs.sizes)
    num.fp <- sum(cumsum(obs.sizes) <= sum(obs.sizes) * threshold)
    num.fp/length(obs.sizes)
}

set.seed(100)
test_that('clusterFDR works correctly', {
    for (mu in c(5, 10, 20)) {
        for (size in c(1, 10, 20)) {
            ids <- rnbinom(100, mu=mu, size=size)

            for (threshold in c(0.05, 0.1)) { 
                out <- compCFDR(ids, threshold, NULL)
                test.fdr <- clusterFDR(ids, threshold)
                expect_identical(out, test.fdr)

                w <- runif(100)
                out <- compCFDR(ids, threshold, w)
                test.fdr <- clusterFDR(ids, threshold, weight=w)
                expect_identical(out, test.fdr)
            }
        }
    }

    # Silly input checks.
    expect_identical(clusterFDR(integer(0), 0.05), 0)
    expect_error(clusterFDR(integer(0), 0.05, weight=1), "must be the same")

    expect_identical(clusterFDR(runif(100), 0), 0) # threshold of zero => FDR of zero.
    expect_identical(clusterFDR(1, numeric(0)), 0)
})

set.seed(101)
test_that("weighted p-value calculations are correct", {
    # Beta-distributed, with and without weights.
    pvals <- rbeta(1000, 1, 20)
    expect_equal(p.adjust(pvals, method="BH"), csaw:::.weightedFDR(pvals, rep(1, length(pvals))))

    weight <- sample(5, length(pvals), replace=TRUE)
    exp.p <- rep(pvals, weight)
    exp.bh <- p.adjust(exp.p, method="BH")
    expect_equal(exp.bh[cumsum(weight)], csaw:::.weightedFDR(pvals, weight))

    # Uniformly-distributed, with and without weights.
    pvals <- runif(1000)
    expect_equal(p.adjust(pvals, method="BH"), csaw:::.weightedFDR(pvals, rep(1, length(pvals))))

    weight <- sample(5, length(pvals), replace=TRUE)
    exp.p <- rep(pvals, weight)
    exp.bh <- p.adjust(exp.p, method="BH")
    expect_equal(exp.bh[cumsum(weight)], csaw:::.weightedFDR(pvals, weight))

    # Other bits and pieces.
    pvals <- rep(1, 1000)
    expect_equal(p.adjust(pvals, method="BH"), csaw:::.weightedFDR(pvals, rep(1, length(pvals))))

    expect_equal(numeric(0), csaw:::.weightedFDR(numeric(0), numeric(0)))
})

set.seed(102)
test_that("controlClusterFDR works as expected", {
    nfalse <- 1000
    ntrue <- 100

    for (nsites in c(5, 10, 20)) {
        for (target in c(0.01, 0.05, 0.1)) {
            # Setting up situations with small and large clusters, where the latter are always detected (strong true DB).
            # This provides the maximal chance that the window- and cluster-level FDRs are different.
            ids <- c(seq_len(nfalse), nfalse + sample(nsites, ntrue, replace=TRUE))
            p <- c(runif(nfalse), numeric(ntrue))
            FUN <- function(is.sig) ids[is.sig]

            out <- controlClusterFDR(target=target, adjp=p, FUN=FUN)
            expect_true(out$FDR <= target)
            expect_true(out$threshold <= target)
            expect_identical(clusterFDR(FUN(p <= out$threshold), out$threshold), out$FDR)

            # Exceeding the window-level threshold should generally increase the cluster-level FDR
            # (note, possible failures here as this depends on resolution).
            for (up in c(1.05, 1.1, 1.5, 2)) {
                a.bit.up <- out$threshold*up
                expect_true(clusterFDR(FUN(p <= a.bit.up), a.bit.up) > target || a.bit.up > target)
            }
        }
    }

    # Caps at the target.
    expect_identical(controlClusterFDR(target=0.05, adjp=0, FUN=function(is.sig) { 1 })$threshold, 0.05)
    expect_identical(controlClusterFDR(target=0.1, adjp=0, FUN=function(is.sig) { 1 })$threshold, 0.1)
    expect_identical(controlClusterFDR(target=0.01, adjp=0, FUN=function(is.sig) { 1 })$threshold, 0.01)

    # Checking correct behaviour for empty inputs.
    out <- controlClusterFDR(target=0.05, adjp=numeric(0), FUN=FUN)
    expect_identical(out$threshold, 0.05)
    expect_identical(out$FDR, 0)
})

##################################################
# Reference results for the consolidation function; also implicitly tests clusterWindows.
# There's not much point doing exact checks, because we'd just be re-implementing most of it.

checkResults <- function(data.list, result.list, target, pval.col="PValue", ..., true.pos) {
    out <- consolidateClusters(data.list, result.list, pval.col=pval.col, target=target, ...)

    # Checking that the clustering is fine.
    all.ids <- unlist(out$id)
    ref <- splitAsList(do.call(c, data.list), all.ids)
    names(ref) <- NULL
    expect_identical(unlist(range(ref)), out$region)

    # Checking that the right windows were chosen.
    all.ps <- unlist(lapply(result.list, FUN=function(x) { x[,pval.col] }))
    was.sig <- !is.na(all.ids)
    if (any(was.sig) && any(!was.sig)) { 
        expect_true(max(all.ps[was.sig]) < min(all.ps[!was.sig])) 
    }

    # Comparing the observed and estimated FDRs.
    np <- out$region[!overlapsAny(out$region, true.pos),]
    expect_true(length(np)/length(out$region) <= out$FDR * 1.1) # A bit fragile, hence the tolerance.
    expect_true(out$FDR <= target)
    return(NULL)
}

set.seed(103)
test_that("consolidateClusters works as expected", {
    windows <- GRanges("chrA", IRanges(1:1000, 1:1000))
    test.p <- runif(1000)
    test.p[rep(1:2, 100) + rep(0:99, each=2) * 10] <- 0 
    
    true.pos <- windows[test.p==0]
    checkResults(list(windows), list(data.frame(PValue=test.p)), tol=0, target=0.05, true.pos=true.pos)
    checkResults(list(windows), list(data.frame(PValue=test.p)), tol=10, target=0.05, true.pos=true.pos)
    
    checkResults(list(windows, windows[1:10]), list(data.frame(PValue=test.p), data.frame(PValue=test.p[1:10])), tol=0, target=0.05, true.pos=true.pos) # Multiple entries
    checkResults(list(windows, windows[1:10]), list(data.frame(PValue=test.p), data.frame(PValue=test.p[1:10])), equiweight=FALSE, tol=0, target=0.05, true.pos=true.pos)
    
    # Smaller number of true positive regions
    set.seed(50)
    test.p <- runif(1000)
    test.p[rep(1:2, 50) + rep(0:49, each=2) * 10] <- 0  
    
    true.pos <- windows[test.p==0]
    checkResults(list(windows), list(data.frame(PValue=test.p)), tol=0, target=0.05, true.pos=true.pos)
    checkResults(list(windows), list(data.frame(PValue=test.p)), tol=5, target=0.05, true.pos=true.pos)
    checkResults(list(windows), list(data.frame(PValue=test.p)), tol=5, target=0.1, true.pos=true.pos)
    checkResults(list(windows), list(data.frame(whee=test.p)), tol=2, pval.col="whee", target=0.05, true.pos=true.pos)
 
    # Adding sign information.
    signs <- ifelse(rbinom(1000, 1, 0.5)!=0L, 1, -1)
    checkResults(list(windows, windows[1:10]), list(data.frame(PValue=test.p, logFC=signs), data.frame(PValue=test.p[1:10], logFC=signs[1:10])), 
                 tol=0, target=0.05, true.pos=true.pos)
    checkResults(list(windows, windows[1:10]), list(data.frame(PValue=test.p, logFC=signs), data.frame(PValue=test.p[1:10], logFC=signs[1:10])), 
                 tol=0, fc.col="logFC", target=0.05, true.pos=true.pos)
    
    # Checking behaviour when empty.
    out <- consolidateClusters(list(windows[0]), list(data.frame(PValue=numeric(0))), tol=0, target=0.05)
    expect_identical(out$id, list())
    expect_identical(out$FDR, 0)
    expect_s4_class(out$region, "GRanges")
    expect_identical(length(out$region), 0L)
})

