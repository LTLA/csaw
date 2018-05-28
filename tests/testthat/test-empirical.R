# This tests that the empirical FDR performs its calculations correctly.
# library(csaw); library(testthat); source("test-empirical.R")

autogen <- function(true.clust, false.clust, nwindows) {
    n.clusters <- true.clust + false.clust
    merged.ids <- sample(n.clusters, nwindows, replace=TRUE)
    tab <- data.frame(logFC=abs(rnorm(nwindows)), logCPM=runif(nwindows, -2, 1), PValue=rbeta(nwindows, 1, 100))

    # Adding true nulls.
    is.false <- merged.ids <= false.clust
    tab$PValue[is.false] <- runif(sum(is.false))
    tab$logFC[is.false] <- rnorm(sum(is.false))
    return(list(tab=tab, id=merged.ids))
}

set.seed(90000)
test_that("empiricalFDR works correctly with vanilla input", {
    for (true.clust in c(10, 100)) {
        for (false.clust in c(50, 500)) {
            clust.out <- autogen(true.clust, false.clust, 1000)
            tab <- clust.out$tab
            merged.ids <- clust.out$id

            out <- empiricalFDR(merged.ids, tab)
            expect_identical(rownames(out), as.character(sort(unique(merged.ids))))
        
            # Checking calculations for the p-values. 
            new.p <- tab$PValue/2
            new.p[tab$logFC < 0] <- 1 - new.p[tab$logFC < 0]
            tab2 <- tab
            tab2$PValue <- new.p
            ref <- combineTests(merged.ids, tab2)
            expect_equal(out$PValue, ref$PValue)
        
            alt <- empiricalFDR(merged.ids, tab, neg.down=FALSE)
            tab2 <- tab
            tab2$PValue <- 1-new.p
            ref <- combineTests(merged.ids, tab2)
            expect_equal(alt$PValue, ref$PValue)
        
            # Checking calculations for the FDR.
            emp.fdr <- findInterval(out$PValue, sort(alt$PValue))/rank(out$PValue, ties.method="max")
            emp.fdr <- pmin(emp.fdr, 1)
            o <- order(out$PValue, decreasing=TRUE)
            emp.fdr[o] <- cummin(emp.fdr[o])
            expect_equal(emp.fdr, out$FDR)
        }
    }
})

set.seed(90001)
test_that("empiricalFDR controls the FDR correctly", {
    for (ntrue in c(1000, 2000)) {
        for (nfalse in c(2000, 5000)) {
            overall <- numeric(20)
            for (it in seq_along(overall)) {
                clust.out <- autogen(1000, 5000, 10000)
                out <- empiricalFDR(clust.out$id, clust.out$tab)
    
                is.sig <- out$FDR <= 0.05
                is.false <- as.integer(rownames(out)) <= nfalse 
                overall[it] <- sum(is.sig & is.false)/sum(is.sig)
            }
            expect_true(mean(overall) < 0.05  * 1.1) # tolerance.
        }
    }
})

set.seed(90002)
test_that("empiricalFDR works correctly with weighted input", {
    for (true.clust in c(10, 100)) {
        for (false.clust in c(50, 500)) {
            clust.out <- autogen(true.clust, false.clust, 1000)
            tab <- clust.out$tab
            merged.ids <- clust.out$id

            weight <- runif(length(merged.ids))
            out <- empiricalFDR(merged.ids, tab, weight=weight)

            out2 <- empiricalFDR(merged.ids, tab)
            expect_identical(rownames(out), rownames(out2))
            expect_identical(out$logFC.up, out2$logFC.up)
            expect_identical(out$logFC.down, out2$logFC.down)

            # Checking calculations for the p-values. 
            new.p <- tab$PValue/2
            new.p[tab$logFC < 0] <- 1 - new.p[tab$logFC < 0]
            tab2 <- tab
            tab2$PValue <- new.p
            ref <- combineTests(merged.ids, tab2, weight=weight)
            expect_equal(out$PValue, ref$PValue)
         }
    }
})

set.seed(90003)
test_that("empiricalFDR works correctly with alternative options", {
    clust.out <- autogen(100, 500, 1000)
    tab <- clust.out$tab
    merged.ids <- clust.out$id
    out <- empiricalFDR(merged.ids, tab)
    
    # Checking that we get the same result with character input.
    out2 <- empiricalFDR(as.character(merged.ids), tab)
    out2 <- out2[rownames(out),]
    expect_identical(out, out2)

    # Checking that we get the same result with different column headings.
    tab2 <- tab
    colnames(tab2) <- c("whee", "blah", "yay")
    expect_error(empiricalFDR(merged.ids, tab2))
    expect_error(empiricalFDR(merged.ids, tab2))

    out3 <- empiricalFDR(merged.ids, tab2, fc.col="whee", pval.col="yay")
    expect_equal(out3$yay, out$PValue)
    expect_equal(out3$whee.up, out$logFC.up)
    expect_equal(out3$whee.down, out$logFC.down)

    # Works correctly with empty input.
    out <- empiricalFDR(integer(0), data.frame(PValue=numeric(0), logCPM=numeric(0), logFC=numeric(0)), weight=numeric(0))
    expect_identical(nrow(out), 0L)
})
