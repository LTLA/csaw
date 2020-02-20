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
            tabup <- combineTests(merged.ids, tab2)
            expect_equal(out$PValue, tabup$PValue)
        
            alt <- empiricalFDR(merged.ids, tab, neg.down=FALSE)
            tab2 <- tab
            tab2$PValue <- 1-new.p
            tabdown <- combineTests(merged.ids, tab2)
            expect_equal(alt$PValue, tabdown$PValue)

            # Checking calculations for the FDR.
            emp.fdr <- findInterval(out$PValue, sort(alt$PValue))/rank(out$PValue, ties.method="max")
            emp.fdr <- pmin(emp.fdr, 1)
            o <- order(out$PValue, decreasing=TRUE)
            emp.fdr[o] <- cummin(emp.fdr[o])
            expect_equal(emp.fdr, out$FDR)

            # Checking directionality.
            expect_identical(out[,1:2], tabup[,1:2])
            expect_identical(out[,3], tabdown[,3])
            expect_identical(out$rep.test, tabup$rep.test)
            expect_identical(out$rep.logFC, tabup$rep.logFC)
            expect_identical(alt$rep.test, tabdown$rep.test)
            expect_identical(alt$rep.logFC, tabdown$rep.logFC)
        }
    }
})

set.seed(90001)
test_that("empiricalFDR controls the FDR correctly", {
    for (ntrue in c(500, 1000, 2000, 5000)) {
        for (nfalse in c(500, 1000, 2000, 5000)) {
            original <- autogen(ntrue, nfalse, 10000)
            clust.out <- original

            # Forcing the p-value distribution for true nulls to be the same for pos/neg logFCs.
            # This ensures that the observed FDR == empirical FDR.
            fakes <- clust.out$id <= nfalse
            clust.out$id <- c(clust.out$id[fakes], clust.out$id[fakes], clust.out$id[!fakes])

            retab1 <- retab2 <- clust.out$tab[fakes,]
            retab1$logFC <- abs(retab1$logFC)
            retab2$logFC <- -abs(retab2$logFC)
            clust.out$tab <- rbind(retab1, retab2, clust.out$tab[!fakes,])

            # Empirical FDR MUST be below the threshold.
            out <- empiricalFDR(clust.out$id, clust.out$tab)
            for (target in c(0.01, 0.05, 0.1)) { 
                is.sig <- out$FDR <= target
                is.false <- as.integer(rownames(out)) <= nfalse 
                expect_true(!any(is.sig) || sum(is.sig & is.false)/sum(is.sig) <= target)
            }

            # Same as above, but with new IDs for the true null with pos/neg log-fold changes.
            clust.out$id <- c(original$id[fakes], original$id[fakes] + nfalse, original$id[!fakes] + nfalse)
            out <- empiricalFDR(clust.out$id, clust.out$tab)
            for (target in c(0.01, 0.05, 0.1)) { 
                is.sig <- out$FDR <= target
                is.false <- as.integer(rownames(out)) <= nfalse * 2L
                expect_true(!any(is.sig) || sum(is.sig & is.false)/sum(is.sig) <= target)
            }
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
