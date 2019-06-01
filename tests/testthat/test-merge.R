# This tests the mergeWindows function, separately from the combineTests function.
# library(csaw); library(testthat); source("test-merge.R")

nativemerge <- function(reg, tol, sign=NULL) {
    n <- length(reg)
    o <- order(reg)
    reg <- reg[o]
    
    increment <- integer(n)
    last.end <- end(reg)
    by.chr <- split(seq_len(n), as.character(seqnames(reg)))
    for (x in by.chr) {
        increment[x[1]] <- 1L
        last.end[x] <- cummax(last.end[x])
    }
    to.next <- c(0L, start(reg)[-1] - last.end[-n] - 1L)
    increment[to.next > tol] <- 1L

    # If a sign is supplied...
    if (!is.null(sign)) { 
        posfc <- sign[o]
        altered.sign <- c(TRUE, posfc[-1]!=posfc[-n])
        increment[altered.sign] <- 1L
    }

    merged.ids <- cumsum(increment)
    merged.ids[o] <- merged.ids
    return(merged.ids)
}

set.seed(900001)
test_that("mergeWindows works correctly in the basic case", {
    chromos <- c(chrA=10000, chrB=5000, chrC=2000) 
    for (nwin in c(50, 100, 200)) {
        for (winsize in c(1, 10, 100)) {
            for (tol in c(10, 50, 200)) {

                reg <- generateWindows(chromos, nwin=nwin, winsize=winsize)
                merged.ids <- nativemerge(reg, tol)
                out <- mergeWindows(reg, tol=tol)
                expect_identical(merged.ids, out$id)

                # Checking the reported coordinates of each region.
                ostarts <- aggregate(start(reg) ~ merged.ids, FUN=min, data=NULL)
                expect_identical(ostarts[,2], start(out$regions))
                oends <- aggregate(end(reg) ~ merged.ids, FUN=max, data=NULL)
                expect_identical(oends[,2], end(out$regions))
                expect_identical(seqnames(reg), seqnames(out$regions[merged.ids]))

                expect_true(all(strand(out$regions)=="*"))
            }
        }
    }
})

set.seed(900002)
test_that("mergeWindows works correctly with signed non-nested windows", {
    chromos <- c(chrA=10000, chrB=5000, chrC=2000) 
    for (nwin in c(50, 100, 200)) {
        for (winsize in c(1, 10, 100)) {
            for (tol in c(10, 50, 200)) {

                # Note that these are fixed width, so open-ended nesting cannot occur.
                reg <- generateWindows(chromos, nwin=nwin, winsize=winsize)
                posfc <- rbinom(length(reg), 1, 0.5)==1L
                merged.ids <- nativemerge(reg, tol, sign=posfc)
                out <- mergeWindows(reg, tol=tol, sign=posfc)
                expect_identical(merged.ids, out$id)

                # Checking the reported value of each region.
                ostarts <- aggregate(start(reg)~ merged.ids, FUN=min, data=NULL)
                expect_identical(ostarts[,2], start(out$regions))
                oends <- aggregate(end(reg)~merged.ids, FUN=max, data=NULL)
                expect_identical(oends[,2], end(out$regions))
                expect_identical(seqnames(reg), seqnames(out$regions[merged.ids]))

                expect_true(all(strand(out$regions)=="*"))
            }
        }
    }
})

set.seed(900003)
test_that("mergeWindows works correctly with nested signed windows", {
    # Should be okay, start point equality
    gr <- GRanges("chrA", IRanges(c(1,1,1), c(10, 30, 50))) 
    x <- mergeWindows(gr, tol=10, sign=c(TRUE, TRUE, TRUE))
    expect_identical(x$id, c(1L, 1L, 1L))
    x <- mergeWindows(gr, tol=10, sign=c(TRUE, FALSE, TRUE))
    expect_identical(x$id, 1:3)
    
    # Should be okay, end point equality
    gr <- GRanges("chrA", IRanges(c(10, 20, 40), c(200, 200, 200))) 
    x <- mergeWindows(gr, tol=10, sign=c(TRUE, TRUE, TRUE))
    expect_identical(x$id, c(1L, 1L, 1L))
    x <- mergeWindows(gr, tol=10, sign=c(TRUE, FALSE, TRUE))
    expect_identical(x$id, 1:3)

    # Nested and of different sign; no merging should happen.
    gr <- GRanges("chrA", IRanges(c(1, 3, 50), c(200, 100, 80))) 
    x <- mergeWindows(gr, tol=10, sign=c(TRUE, FALSE, TRUE))
    expect_identical(x$id, 1:3)
   
    # A more complex nesting; no merging should happen.
    gr2 <- GRanges("chrA", IRanges(c(1, 3, 50, 90), c(200, 100, 80, 1000))) 
    x <- mergeWindows(gr2, tol=10, sign=c(TRUE, FALSE, TRUE, FALSE))
    expect_identical(x$id, 1:4)
})

splitter <- function(ids, regions, max.width) {
    by.id <- split(seq_along(regions), ids)
    new.ids <- integer(length(ids))
    last <- 0L

    for (i in seq_along(by.id)) {
        chosen <- by.id[[i]]
        cur.start <- start(regions)[chosen]
        cur.end <- end(regions)[chosen]

        total.width <- max(cur.end) - min(cur.start) + 1L
        MULT <- ceiling(total.width/max.width)
        sub.width <- total.width/MULT
        
        subcluster <- (cur.start + cur.end)/2 - min(cur.start)
        subcluster <- as.integer(floor(subcluster/sub.width))
        subcluster <- as.integer(factor(subcluster)) 

        new.ids[chosen] <- subcluster + last
        last <- last + max(subcluster)
    }
    
    return(new.ids)
}

set.seed(900004)
test_that("mergeWindows responds to the maximum width", {
    chromos <- c(chrA=10000, chrB=5000, chrC=2000) 
    for (winsize in c(1, 10, 100)) {
        for (tol in c(10, 50, 200)) {
            for (maxwidth in c(10000, 1000, 100)) {

                # Without signage:
                reg <- generateWindows(chromos, nwin=100, winsize=winsize)
                free <- mergeWindows(reg, tol=tol)
                restricted <- mergeWindows(reg, tol=tol, max.width=maxwidth)
                expect_identical(restricted$id, splitter(free$id, reg, maxwidth))
                
                # With signage:
                posfc <- rbinom(length(reg), 1, 0.5)==1L
                free <- mergeWindows(reg, tol=tol, sign=posfc)
                restricted <- mergeWindows(reg, tol=tol, sign=posfc, max.width=maxwidth)
                expect_identical(restricted$id, splitter(free$id, reg, maxwidth))
            }
        }
    }
})

set.seed(900005)
test_that("mergeWindows works correctly for stranded input", {
    chromos <- c(chrA=2000, chrB=5000, chrC=1000) 
    for (nwin in c(50, 100, 200)) {
        for (winsize in c(1, 10, 100)) {
            for (tol in c(10, 50, 200)) {

                reg <- generateWindows(chromos, nwin=nwin, winsize=winsize)
                strand(reg) <- sample(c("+", "-", "*"), length(reg), replace=TRUE)
                combo <- mergeWindows(reg, tol=tol, ignore.strand=FALSE)

                # Checking that each set of IDs only has one strand, and it is the expected one.
                by.id <- split(strand(reg), combo$id)
                expect_true(all(lengths(lapply(by.id, FUN=runValue))==1))
                expect_identical(strand(combo$regions)[combo$id], strand(reg))

                # Running separately on each strand, and checking that the boundaries are the same.
                is.forward <- as.logical(strand(reg)=="+")
                forward <- mergeWindows(reg[is.forward], tol=tol)
                is.reverse <- as.logical(strand(reg)=="-")
                reverse <- mergeWindows(reg[is.reverse], tol=tol)
                is.unstrand <- as.logical(strand(reg)=="*")
                unstrand <- mergeWindows(reg[is.unstrand], tol=tol)

                strand(forward$regions) <- "+"
                strand(reverse$regions) <- "-"
                strand(unstrand$regions) <- "*"
                expect_identical(c(forward$regions, reverse$regions, unstrand$regions), combo$regions)

                final.out <- integer(length(reg))
                final.out[is.forward] <- forward$id
                final.out[is.reverse] <- reverse$id+length(forward$regions) 
                final.out[is.unstrand] <- unstrand$id+length(forward$regions)+length(reverse$regions)
                expect_identical(final.out, combo$id)
            }
        }
    }
})

test_that("mergeWindows works correctly for silly inputs", {
    # empty inputs.
    out <- mergeWindows(GRanges(), tol=10)
    expect_identical(out$id, integer(0))
    expect_identical(out$regions, GRanges())

    out <- mergeWindows(GRanges(), tol=10, max.width=1000)
    expect_identical(out$id, integer(0))
    expect_identical(out$regions, GRanges())

    # mismatch inputs.
    expect_error(mergeWindows(GRanges("chrA:1-1000"), tol=10, sign=logical(5)), "must be the same")
})
