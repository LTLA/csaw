# This tests the various functions related to consolidation.
# library(csaw); library(testthat); source("setup.R"); source("test-cons.R")

chromos <- c(chrA=12332, chrB=34892)

set.seed(200000)
test_that("mergeWindowsList works correctly", {
    for (nwin in c(5, 50, 500)) {
        for (sizes in list(c(10, 50), c(50, 50), c(100, 50))) {
            all.win <- list()
            for (s in seq_along(sizes)) {
                all.win[[s]] <- generateWindows(chromos, winsize=sizes[s], nwin=nwin)
            }

            combined.win <- do.call(c, all.win)
            ref <- mergeWindows(combined.win, tol=100)

            cons <- mergeWindowsList(all.win, tol=100)
            expect_identical(cons$ids, ref$ids)
            expect_identical(cons$regions, ref$regions)
            expect_true(all(combined.win==cons$ranges))

            # Works with GRanges objects.
            cons2 <- mergeWindowsList(all.win, tol=100)
            expect_identical(cons, cons2)

            # Weights are computed correctly.
            origins <- rep(seq_along(all.win), lengths(all.win))
            f <- paste0(cons$id, ".", origins)
            summed <- by(cons$weights, INDICES=f, FUN=sum)
            expect_equal(as.numeric(summed), rep(1, length(unique(f))))

            # Responds to the sign.
            signs <- rbinom(length(combined.win), 1, 0.5)==1L
            cons3 <- mergeWindowsList(all.win, tol=100, signs=signs)
            ref2 <- mergeWindows(combined.win, tol=100, signs=signs)
            expect_identical(cons3$ids, ref2$ids)
            expect_identical(cons3$regions, ref2$regions)
        }
    }

     # Responds correctly to empty inputs.
     cons.E <- mergeWindowsList(list(all.win[[1]], all.win[[2]][0]), tol=100)
     ref <- mergeWindows(all.win[[1]], tol=100)
     expect_identical(cons.E$ids, ref$ids)
     expect_identical(cons.E$regions, ref$regions)

     cons.E <- mergeWindowsList(list(all.win[[1]][0], all.win[[2]]), tol=100)
     ref <- mergeWindows(all.win[[2]], tol=100)
     expect_identical(cons.E$ids, ref$ids)
     expect_identical(cons.E$regions, ref$regions)

     cons.E <- mergeWindowsList(list(all.win[[1]][0], all.win[[2]][0]), tol=100)
     expect_identical(cons.E$ids, integer(0))
     expect_identical(length(cons.E$regions), 0L)
     expect_identical(cons.E$weights, numeric(0))
})

set.seed(200001)
test_that("findOverlapsList works correctly", {
    for (nwin in c(5, 50, 500)) {
        for (sizes in list(c(10, 50), c(50, 50), c(100, 50))) {
            all.win <- list()
            regions <- generateWindows(chromos, winsize=1000, nwin=5)
            for (s in seq_along(sizes)) {
                all.win[[s]] <- generateWindows(chromos, winsize=sizes[s], nwin=nwin)
            }

            combined.win <- do.call(c, all.win)
            ref <- findOverlaps(regions, combined.win)

            cons <- findOverlapsList(all.win, regions)
            expect_identical(cons$overlaps, ref)
            expect_true(all(combined.win==cons$ranges))

            # Works with GRanges objects.
            cons2 <- findOverlapsList(all.win, regions)
            expect_identical(cons, cons2)

            # Checking that we respond to arguments.
            ref <- findOverlaps(regions, combined.win, type="within")
            cons2 <- findOverlapsList(all.win, regions, type="within")
            expect_identical(cons2$overlaps, ref)
            expect_identical(cons$ranges, cons2$ranges)

            # Weights are computed correctly.
            origins <- rep(seq_along(all.win), lengths(all.win))
            f <- paste0(queryHits(cons$overlaps), ".", origins[subjectHits(cons$overlaps)])
            summed <- by(cons$weights, INDICES=f, FUN=sum)
            expect_equal(as.numeric(summed), rep(1, length(unique(f))))
        }
    }

    # Responds correctly to empty inputs.
    cons.E <- findOverlapsList(list(all.win[[1]], all.win[[2]][0]), regions)
    expect_identical(cons.E$overlaps, findOverlaps(regions, all.win[[1]]))

    cons.E <- findOverlapsList(list(all.win[[1]][0], all.win[[2]]), regions)
    expect_identical(cons.E$overlaps, findOverlaps(regions, all.win[[2]]))

    cons.E <- findOverlapsList(list(all.win[[1]][0], all.win[[2]][0]), regions)
    expect_identical(length(cons.E$overlaps), 0L)
    expect_identical(cons.E$weights, numeric(0))
})
