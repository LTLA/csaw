# This tests the correctness of the getBestTest function.
# library(csaw); library(testthat); source("test-best.R")

autogen <- function(nids, alpha, beta) {
	n <- 1000L
    ids <- round(runif(n, 1, nids))
	tab <- data.frame(logFC=rnorm(n), logCPM=rnorm(n), PValue=rbeta(n, alpha, beta))
    return(list(id=ids, table=tab))
}

set.seed(40000)
test_that("getBestTest works as expected on vanilla input", {
    test <- autogen(100, 1, 1)
    tab <- test$table
    ids <- test$id
	best <- getBestTest(ids, tab)

    # Checking we get the same result after Bonferroni correction.
	ref <- aggregate(tab$PValue ~ ids, FUN=function(x) { min(1, x*length(x)) }, data=NULL)	
	xref <- aggregate(seq_along(ids) ~ ids, FUN=function(x) { x[which.min(tab$PValue[x])] }, data=NULL)	
    expect_equal(best$PValue, ref[,2])
    expect_equal(best$best, xref[,2])
    expect_identical(rownames(best), as.character(sort(unique(ids))))

    # After shuffling things around, to check that internal re-orderings have no effect.
    re.o <- rev(seq_along(ids))
    out2 <- getBestTest(ids[re.o], tab[re.o,])
    expect_equal(best$PValue, out2$PValue)
    expect_identical(rownames(best), rownames(out2))
    expect_identical(out2$best, nrow(tab) - best$best + 1L) # because the ordering of 'tab' changes.

    # Testing what happens with a character vector as input.
    out3 <- getBestTest(as.character(ids), tab)
    expect_identical(sort(rownames(best)), sort(rownames(out3)))
    expect_identical(best, out3[rownames(best),])
})

set.seed(40001)
test_that("getBestTest works with alternative options", {
    test <- autogen(100, 1, 10)
    tab <- test$table
    ids <- test$id

	# Effectively frequency weights.
	w <- runif(length(ids), 1, 10)
	best <- getBestTest(ids, tab, weight=w)
	ref <- aggregate(seq_along(ids) ~ ids, FUN=function(x) { min(1, tab$PValue[x]/w[x]*sum(w[x])) }, data=NULL)	
	xref <- aggregate(seq_along(ids)~ ids, FUN=function(x) { x[which.min(tab$PValue[x]/w[x])] }, data=NULL)	
    expect_equal(best$PValue, ref[,2])
    expect_equal(best$best, xref[,2])

    # Weights also respond to re-ordering.
    re.o <- rev(seq_along(ids))
    out2 <- getBestTest(ids[re.o], tab[re.o,], weight=w[re.o])
    expect_equal(best$PValue, out2$PValue)
    expect_identical(rownames(best), rownames(out2))
    expect_identical(out2$best, nrow(tab) - best$best + 1L) # because the ordering of 'tab' changes.

	# Now, searching for the max log-CPM.
	mostab <- getBestTest(ids, tab, by.pval=FALSE)
    ref <- aggregate(seq_along(ids) ~ ids, FUN=function(x) { x[which.max(tab$logCPM[x])] }, data=NULL)
    expect_identical(ref[,2], mostab$best)

    # Changing the column headings.
    ref <- getBestTest(ids, tab)
    retab <- tab
    colnames(retab) <- c("whee", "blah", "yay")
    expect_error(getBestTest(ids, retab), "failed to find")

    out <- getBestTest(ids, retab, pval.col="yay")
    expect_equal(out$yay, ref$PValue)
    expect_identical(out$best, ref$best)

    out2 <- getBestTest(ids, retab, pval.col=3)
    expect_equal(out2, out)

    expect_error(getBestTest(ids, retab, pval.col=c("yay", "yay")), "multiple")
})


set.seed(40001)
test_that("getBestTest handles edge cases correctly", {
    # Checking that unsaturation does not compromise the results.
    test <- autogen(2000, 1, 2)
    tab <- test$table
    ids <- test$id
    expect_true(any(tabulate(ids)==0L)) # saturated.
    
    out <- getBestTest(ids, tab)
    asfac <- factor(ids)
    out2 <- getBestTest(as.integer(asfac), tab)
    expect_identical(rownames(out), levels(asfac)[as.integer(rownames(out2))])
    rownames(out) <- rownames(out2)
    expect_equal(out, out2)

    # Checking what happens if the first id becomes NA.
    test <- autogen(20, 1, 1)
    tab <- test$table
    ids <- test$id
 
    na.ids <- ids 
    na.ids[1] <- NA_integer_
    out.na <- getBestTest(na.ids, tab, weight=w)
    out.ref <- getBestTest(na.ids[-1], tab[-1,], weight=w[-1])
    expect_equal(out.na$PValue, out.ref$PValue)
    expect_identical(rownames(out.na), rownames(out.ref))
    expect_identical(out.na$best, out.ref$best + 1L)  # because of the missing first row.

    # Checking for sane behaviour when no IDs are supplied.
    emp <- getBestTest(ids[0], tab[0,])
    expect_identical(nrow(emp), 0L)
    expect_error(getBestTest(ids[0], tab))
    expect_error(getBestTest(ids, tab[0,]))
    expect_error(getBestTest(ids, tab, weight=numeric(0)))
})

