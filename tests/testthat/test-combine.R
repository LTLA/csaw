# This tests the combining power of the combineTests function.
# library(csaw); library(testthat); source("test-combine.R")

autogen <- function(n.clusters, total.n) {
	ids <- sample(n.clusters, total.n, replace=TRUE)
	tab <- data.frame(logFC=runif(total.n, -1, 1), logCPM=runif(total.n, -2, 1), PValue=rbeta(total.n, 1, 10))
    return(list(id=ids, table=tab))
}

set.seed(50000)
test_that("combineTests works as expected on vanilla inputs", {
    # Always testing differing numbers of windows relative to the number of clusters,
    # to explore the gamut of behaviours when there are few to many windows per cluster.
    for (ntests in c(10, 50, 100, 200)) {
        test <- autogen(20, ntests)
        ids <- test$id
        tab <- test$table
        tabcom <- combineTests(ids, tab)

        # Checking numbers
        by.id <- table(ids)
        expect_identical(as.integer(by.id), tabcom$nWindows)
        expect_identical(names(by.id), rownames(tabcom))

        test.up <- table(ids[tab$logFC > 0.5])
        reference <- as.integer(test.up)[match(rownames(tabcom), names(test.up))]
        reference <- ifelse(is.na(reference), 0L, reference) 
        expect_identical(reference, tabcom$logFC.up)

        test.down <- table(ids[tab$logFC < -0.5])
        reference <- as.integer(test.down)[match(rownames(tabcom), names(test.down))]
        reference <- ifelse(is.na(reference), 0L, reference) 
        expect_identical(reference, tabcom$logFC.down)

        # Checking Simes.
        checker <- split(data.frame(PValue=tab$PValue), ids)
        osimes <- sapply(checker, FUN=function(x) { min(p.adjust(x$PValue, method="BH")) })
        osimes <- unname(osimes)
        expect_equal(osimes, tabcom$PValue)
        expect_equal(p.adjust(osimes, method="BH"), tabcom$FDR)

        # Checking if we get the same results after reversing the ids (ensures internal re-ordering is active).
        re.o <- rev(seq_along(ids))
        out2 <- combineTests(ids[re.o], tab[re.o,])
        expect_equal(tabcom, out2)

        # Checking we get the same results with a character vector (though ordering might be different).
        out3 <- combineTests(as.character(ids), tab)
        out3 <- out3[rownames(tabcom),]
        expect_equal(tabcom, out3)
    }
})

set.seed(50001)
test_that("combineTests works with alternative options", {
    for (ntests in c(10, 50, 100, 200)) {
        test <- autogen(20, ntests)
        ids <- test$id
        tab <- test$table
        
        w <- runif(length(ids), 1, 5)
        tabcom <- combineTests(ids, tab, weight=w)
        uweight <- combineTests(ids, tab)
        expect_identical(tabcom[,c("nWindows", "logFC.up", "logFC.down")],
            uweight[,c("nWindows", "logFC.up", "logFC.down")])

        # Checking weighted Simes.
        checker <- split(data.frame(PValue=tab$PValue, weight=w), ids)
        osimes <- sapply(checker, FUN=function(x) {
            o<-order(x$PValue)
            min(x$PValue[o]/cumsum(x$weight[o])) * sum(x$weight)
        })
        osimes <- unname(osimes)
        expect_equal(osimes, tabcom$PValue)
        expect_equal(p.adjust(osimes, method="BH"), tabcom$FDR)
        
        # Weights respond correctly to re-ordering.
        re.o <- rev(seq_along(ids))
        out2 <- combineTests(ids[re.o], tab[re.o,], weight=w[re.o])
        expect_equal(tabcom, out2)

        # Adding some tests if there's multiple log-FC's in 'tab'.
        retab <- tab
        is.fc <- which(colnames(retab)=="logFC")
        colnames(retab)[is.fc] <- "logFC.1"
        retab$logFC.2 <- -retab$logFC.1

        ref <- combineTests(ids, tab)
        out <- combineTests(ids, retab)
        expect_identical(ref$logFC.up, out$logFC.1.up)
        expect_identical(ref$logFC.down, out$logFC.1.down)
        expect_identical(ref$logFC.up, out$logFC.2.down)
        expect_identical(ref$logFC.down, out$logFC.2.up)
        expect_identical(ref$PValue, out$PValue)

        retab <- tab
        colnames(retab) <- c("whee", "blah", "yay")
        out4 <- combineTests(ids, retab, pval.col="yay", fc.col="whee")
        expect_equal(ref$logFC.up, out4$whee.up)
        expect_equal(ref$logFC.down, out4$whee.down)
        expect_equal(ref$PValue, out4$yay)
    }
})

set.seed(50002)
test_that("combineTests direction inference works as expected", {
    for (nclusters in c(10, 20, 50, 100)) { 
        test <- autogen(nclusters, 200)
        ids <- test$id
        tab <- test$table
        tabcom <- combineTests(ids, tab)

        # Checking inferred direction by comparing to behaviour when 
        # the p-value for all tests in one direction are set to 1.
        going.up <- tab$logFC > 0
        tab.up <- tab
        tab.up$PValue[!going.up] <- 1
        out.up <- combineTests(ids, tab.up)

        tab.down <- tab
        tab.down$PValue[going.up] <- 1
        out.down <- combineTests(ids, tab.down)
        
        direction <- rep("mixed", nrow(tabcom))
        tol <- 1e-6
        up.same <- out.up$PValue/tabcom$PValue - 1 <= tol  # No need to use abs(), up/down cannot be lower.
        down.same <- out.down$PValue/tabcom$PValue - 1 <= tol
        direction[up.same & !down.same] <- "up"
        direction[!up.same & down.same] <- "down"

        expect_identical(direction, tabcom$direction)

        # Also works when weights get involved.
        w <- runif(length(ids), 1, 2)
        tabcom <- combineTests(ids, tab, weight=w)
        out.up <- combineTests(ids, tab.up, weight=w)
        out.down <- combineTests(ids, tab.down, weight=w)

        direction <- rep("mixed", nrow(tabcom))
        tol <- 1e-6
        up.same <- out.up$PValue/tabcom$PValue - 1 <= tol  # No need to use abs(), up/down cannot be lower.
        down.same <- out.down$PValue/tabcom$PValue - 1 <= tol
        direction[up.same & !down.same] <- "up"
        direction[!up.same & down.same] <- "down"

        expect_identical(direction, tabcom$direction)
    }
})

set.seed(50003)
test_that("combineTests handles edge cases correctly", {
    # Checking that unsaturation does not compromise the results.
    test <- autogen(100, 50)
    tab <- test$table 
    ids <- test$id
    expect_true(any(tabulate(ids)==0L)) # saturated.

    out <- combineTests(ids, tab)
    asfac <- factor(ids)
    out2 <- combineTests(as.integer(asfac), tab)
    expect_identical(rownames(out), levels(asfac)[as.integer(rownames(out2))])
    rownames(out) <- rownames(out2)
    expect_equal(out, out2)

    # Checking what happens if some proportion of the IDs become NA.
    for (prop in c(0.2, 0.5, 0.8)) {
        test <- autogen(20, 100)
        ids <- test$id
        tab <- test$table

        invalid <- sample(length(ids), length(ids) * prop)
        na.ids <- ids
        na.ids[invalid] <- NA_integer_
        out.na <- combineTests(na.ids, tab)
        out.ref <- combineTests(na.ids[-invalid], tab[-invalid,])
        expect_equal(out.na, out.ref) 
    }

    # Checking for sane behaviour when no log-fold changes are supplied.
    out3 <- combineTests(ids, tab, fc.col=integer(0))
    ref <- combineTests(ids, tab)
    ref$logFC.up <- ref$logFC.down <- ref$direction <- NULL
    expect_identical(ref, out3)

    # Checking for sane behaviour when no IDs are supplied.
    emp <- combineTests(ids[0], tab[0,])
    expect_identical(nrow(emp), 0L)
    expect_error(combineTests(ids[0], tab))
    expect_error(combineTests(ids, tab[0,]))
    expect_error(combineTests(ids, tab, weight=numeric(0)))
})
