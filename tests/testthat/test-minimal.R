# This tests the combining power of the minimalTests function.
# library(csaw); library(testthat); source("test-minimal.R")

autogen <- function(n.clusters, total.n) {
	ids <- sample(n.clusters, total.n, replace=TRUE)
	tab <- data.frame(logFC=runif(total.n, -1, 1), logCPM=runif(total.n, -2, 1), PValue=rbeta(total.n, 1, 10))
    return(list(id=ids, table=tab))
}

SELECTOR <- function(p) p[min(max(3, ceiling(0.4*length(p))), length(p))]

set.seed(50000)
test_that("minimalTests works as expected on vanilla inputs", {
    # Always testing differing numbers of windows relative to the number of clusters,
    # to explore the gamut of behaviours when there are few to many windows per cluster.
    for (ntests in c(10, 50, 100, 200)) {
        test <- autogen(20, ntests)
        ids <- test$id
        tab <- test$table
        tabcom <- minimalTests(ids, tab)

        # Checking numbers.
        by.id <- table(ids)
        expect_identical(as.integer(by.id), tabcom$num.tests)
        expect_identical(names(by.id), rownames(tabcom))

        p.by.id <- split(tab$PValue, ids)
        p.by.id <- lapply(p.by.id, p.adjust, method="holm")
        ncounter <- function(p, dir) { sum(p[dir] <= 0.05) }

        up.by.id <- split(tab$logFC > 0, ids)
        reference <- mapply(p=p.by.id, dir=up.by.id, FUN=ncounter, USE.NAMES=FALSE)
        expect_identical(reference, tabcom$num.up.logFC)

        down.by.id <- split(tab$logFC < 0, ids)
        reference <- mapply(p=p.by.id, dir=down.by.id, FUN=ncounter, USE.NAMES=FALSE)
        expect_identical(reference, tabcom$num.down.logFC)

        # Checking the representative.
        expect_identical(tab$logFC[tabcom$rep.test], tabcom$rep.logFC)

        # Checking Holm.
        checker <- split(data.frame(PValue=tab$PValue), ids)
        oholm <- sapply(checker, FUN=function(x) { 
            p <- p.adjust(x$PValue, method="holm")
            p <- sort(p)
            SELECTOR(p) 
        })
        oholm <- unname(oholm)
        expect_equal(oholm, tabcom$PValue)
        expect_equal(p.adjust(oholm, method="BH"), tabcom$FDR)

        # Checking if we get the same results after reversing the ids (ensures internal re-ordering is active).
        re.o <- rev(seq_along(ids))
        out2 <- minimalTests(ids[re.o], tab[re.o,])
        out2$rep.test <- re.o[out2$rep.test]
        expect_equal(tabcom, out2)

        # Checking we get the same results with a character vector (though ordering might be different).
        out3 <- minimalTests(as.character(ids), tab)
        out3 <- out3[rownames(tabcom),]
        expect_equal(tabcom, out3)
    }
})

set.seed(50001)
test_that("minimalTests works with alternative options", {
    for (ntests in c(10, 50, 100, 200)) {
        test <- autogen(20, ntests)
        ids <- test$id
        tab <- test$table
        
        w <- runif(length(ids), 1, 5)
        tabcom <- minimalTests(ids, tab, weight=w)
        uweight <- minimalTests(ids, tab)
        expect_identical(tabcom$num.tests, uweight$num.tests)

        # Checking weighted Holm.
        checker <- split(data.frame(PValue=tab$PValue, weight=w), ids)
        oholm <- sapply(checker, FUN=function(x) {
            S <- x$PValue/x$weight 
            o <- order(S)
            S <- S[o] * rev(cumsum(rev(x$weight[o])))
            S <- cummax(S)
            SELECTOR(S)
        })
        oholm <- unname(oholm)
        expect_equal(oholm, tabcom$PValue)
        expect_equal(p.adjust(oholm, method="BH"), tabcom$FDR)

        # Weights respond correctly to re-ordering.
        re.o <- rev(seq_along(ids))
        out2 <- minimalTests(ids[re.o], tab[re.o,], weight=w[re.o])
        out2$rep.test <- re.o[out2$rep.test]
        expect_equal(tabcom, out2)

        # Adding some tests if there's multiple log-FC's in 'tab'.
        retab <- tab
        is.fc <- which(colnames(retab)=="logFC")
        colnames(retab)[is.fc] <- "logFC.1"
        retab$logFC.2 <- -retab$logFC.1

        ref <- minimalTests(ids, tab)
        out <- minimalTests(ids, retab)
        expect_identical(ref$num.up.logFC, out$num.up.logFC.1)
        expect_identical(ref$num.down.logFC, out$num.down.logFC.1)
        expect_identical(ref$num.up.logFC, out$num.down.logFC.2)
        expect_identical(ref$num.down.logFC, out$num.up.logFC.2)
        expect_identical(ref$PValue, out$PValue)

        retab <- tab
        colnames(retab) <- c("whee", "blah", "yay")
        out4 <- minimalTests(ids, retab, pval.col="yay", fc.col="whee")
        expect_equal(ref$num.up.logFC, out4$num.up.whee)
        expect_equal(ref$num.down.logFC, out4$num.down.whee)
        expect_equal(ref$PValue, out4$yay)
    }
})

set.seed(50002)
test_that("minimalTests direction inference works as expected", {
    for (nclusters in c(10, 20, 50, 100)) { 
        test <- autogen(nclusters, 200)
        ids <- test$id
        tab <- test$table
        tabcom <- minimalTests(ids, tab)

        # Checking inferred direction by comparing to behaviour when 
        # the p-value for all tests in one direction are set to 1.
        going.up <- tab$logFC > 0
        tab.up <- tab
        tab.up$PValue[!going.up] <- 1
        out.up <- minimalTests(ids, tab.up)

        tab.down <- tab
        tab.down$PValue[going.up] <- 1
        out.down <- minimalTests(ids, tab.down)

        tol <- 1e-6
        up.same <- out.up$PValue/tabcom$PValue - 1 <= tol  # No need to use abs(), up/down cannot be lower.
        down.same <- out.down$PValue/tabcom$PValue - 1 <= tol

        # Due to ties in Holm procedure, the 'setting p-values to unity' interpretation above is not exact;
        # clusters reported as 'mixed' may actually not change when setting p-values to unity in one direction.
        # But any clusters reported as 'up' MUST be 'up.same' and '!down.same'; opposite is true for reported 'down'.
        # Similarly, if you are '!up.same' and '!down.same', you MUST be reported as 'mixed'.
        expect_identical(up.same & !down.same & tabcom$direction=="up", tabcom$direction=="up")
        expect_identical(!up.same & down.same & tabcom$direction=="down", tabcom$direction=="down")
        expect_identical(!up.same & !down.same & tabcom$direction=="mixed", !up.same & !down.same)

        # Also works when weights get involved.
        w <- runif(length(ids), 1, 2)
        tabcom <- minimalTests(ids, tab, weight=w)
        out.up <- minimalTests(ids, tab.up, weight=w)
        out.down <- minimalTests(ids, tab.down, weight=w)

        tol <- 1e-6
        up.same <- out.up$PValue/tabcom$PValue - 1 <= tol  # No need to use abs(), up/down cannot be lower.
        down.same <- out.down$PValue/tabcom$PValue - 1 <= tol

        expect_identical(up.same & !down.same & tabcom$direction=="up", tabcom$direction=="up")
        expect_identical(!up.same & down.same & tabcom$direction=="down", tabcom$direction=="down")
        expect_identical(!up.same & !down.same & tabcom$direction=="mixed", !up.same & !down.same)
    }
})

set.seed(50003)
test_that("minimalTests handles edge cases correctly", {
    # Checking that unsaturation does not compromise the results.
    test <- autogen(100, 50)
    tab <- test$table 
    ids <- test$id
    expect_true(any(tabulate(ids)==0L)) # saturated.

    out <- minimalTests(ids, tab)
    asfac <- factor(ids)
    out2 <- minimalTests(as.integer(asfac), tab)
    expect_identical(rownames(out), levels(asfac)[as.integer(rownames(out2))])
    rownames(out) <- rownames(out2)
    expect_equal(out, out2)

    # Handles extreme minimal requests.
    out <- minimalTests(ids, tab, min.sig.n=0, min.sig.prop=0)
    best <- getBestTest(ids, tab)
    expect_identical(out, best)

    out <- minimalTests(ids, tab, min.sig.n=100, min.sig.prop=1)
    by.gene <- split(tab$PValue, ids)
    max.p <- vapply(by.gene, FUN=function(p) max(p.adjust(p, "holm")), 0)
    expect_identical(out$PValue, unname(max.p))

    # Checking what happens if some proportion of the IDs become NA.
    for (prop in c(0.2, 0.5, 0.8)) {
        test <- autogen(20, 100)
        ids <- test$id
        tab <- test$table

        invalid <- sample(length(ids), length(ids) * prop)
        na.ids <- ids
        na.ids[invalid] <- NA_integer_
        out.na <- minimalTests(na.ids, tab)
        out.ref <- minimalTests(na.ids[-invalid], tab[-invalid,])
        out.ref$rep.test <- which(!is.na(na.ids))[out.ref$rep.test]
        expect_equal(out.na, out.ref) 
    }

    # Checking for sane behaviour when no log-fold changes are supplied.
    out3 <- minimalTests(ids, tab, fc.col=integer(0))
    ref <- minimalTests(ids, tab)
    ref$num.up.logFC <- ref$num.down.logFC <- ref$direction <- ref$rep.logFC <- NULL
    expect_identical(ref, out3)

    # Checking for sane behaviour when no IDs are supplied.
    emp <- minimalTests(ids[0], tab[0,])
    expect_identical(nrow(emp), 0L)
    expect_error(minimalTests(ids[0], tab))
    expect_error(minimalTests(ids, tab[0,]))
    expect_error(minimalTests(ids, tab, weight=numeric(0)))
})
