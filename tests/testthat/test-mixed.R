# Tests the mixedClusters function.
# library(csaw); library(testthat); source("test-mixed.R")

autogen <- function(n.clusters, total.n) {
    ids <- sample(n.clusters, total.n, replace=TRUE)
    tab <- data.frame(logFC=runif(total.n, -1, 1), logCPM=runif(total.n, -2, 1), PValue=rbeta(total.n, 1, 10))
    return(list(id=ids, table=tab))
}

set.seed(60000)
test_that("mixedClusters works as expected on vanilla inputs", {
    for (ntests in c(10, 50, 100)) {
        test <- autogen(20, ntests)
        ids <- test$id
        tab <- test$table
        tabmix <- mixedClusters(ids, tab)

        uptab <- tab
        uptab$PValue <- ifelse(tab$logFC > 0, tab$PValue/2, 1-tab$PValue/2)
        tabup <- combineTests(ids, uptab)

        downtab <- tab
        downtab$PValue <- ifelse(tab$logFC < 0, tab$PValue/2, 1-tab$PValue/2)
        tabdown <- combineTests(ids, downtab)

        expect_identical(tabmix[,1:3], tabup[,1:3])
        expect_identical(tabmix$PValue, pmax(tabup$PValue, tabdown$PValue))
        expect_identical(tabmix$FDR, p.adjust(tabmix$PValue, method="BH"))
    }
})

set.seed(60001)
test_that("mixedClusters works as expected with weights", {
    for (ntests in c(10, 50, 100)) {
        test <- autogen(20, ntests)
        ids <- test$id
        tab <- test$table
        w <- runif(ntests)
        tabmix <- mixedClusters(ids, tab, weight=w)

        uptab <- tab
        uptab$PValue <- ifelse(tab$logFC > 0, tab$PValue/2, 1-tab$PValue/2)
        tabup <- combineTests(ids, uptab, weight=w)

        downtab <- tab
        downtab$PValue <- ifelse(tab$logFC < 0, tab$PValue/2, 1-tab$PValue/2)
        tabdown <- combineTests(ids, downtab, weight=w)

        expect_identical(tabmix[,1:3], tabup[,1:3])
        expect_identical(tabmix$PValue, pmax(tabup$PValue, tabdown$PValue))
        expect_identical(tabmix$FDR, p.adjust(tabmix$PValue, method="BH"))
    }
})

set.seed(60002)
test_that("mixedClusters works as expected with alternative options", {
    test <- autogen(20, 100)
    ids <- test$id
    tab <- test$table
    ref <- mixedClusters(ids, tab)

    retab <- tab
    retab$whee <- tab$logFC
    retab$yay <- tab$PValue
    retab$logFC <- NULL
    retab$PValue <- NULL
    out <- mixedClusters(ids, retab, fc.col="whee", pval.col="yay")
 
    expect_identical(out$whee.up, ref$logFC.up)
    expect_identical(out$whee.down, ref$logFC.down)
    expect_identical(out$yay, ref$PValue)
})

set.seed(60003) 
test_that("mixedClusters behaves correctly with silly inputs", {
    test <- autogen(20, 100)
    ids <- test$id
    tab <- test$table
    ref <- mixedClusters(ids[0], tab[0,])
    expect_identical(dim(ref), c(0L, 5L))

    expect_error(mixedClusters(ids, tab[0,]), "not TRUE")
    expect_error(mixedClusters(ids, tab, weight=0), "not TRUE")
})
