# This tests the various functions related to consolidation.
# library(csaw); library(testthat); source("setup.R"); source("test-cons.R")

chromos <- c(chrA=12332, chrB=34892)

set.seed(200000)
test_that("mergeWindowsList works correctly", {
    for (nwin in c(5, 50, 500)) {
        for (sizes in list(c(10, 50), c(50, 50), c(100, 50))) {
	        data.list <- list()
        	for (s in seq_along(sizes)) {
        		windows <- generateWindows(chromos, winsize=sizes[s], nwin=nwin)
        		counts <- matrix(0, ncol=1, nrow=length(windows))
        		data.list[[s]] <- SummarizedExperiment(list(counts=counts), rowRanges=windows)
        	}

            all.win <- lapply(data.list, rowRanges)
            combined.win <- do.call(c, all.win)
            ref <- mergeWindows(combined.win, tol=100)

        	cons <- mergeWindowsList(data.list, tol=100)
            expect_identical(cons$id, ref$id)
            expect_identical(cons$merged, ref$region)
            expect_true(all(combined.win==cons$ranges))

            # Works with GRanges objects.
            cons2 <- mergeWindowsList(all.win, tol=100)
            expect_identical(cons, cons2)

            # Weights are computed correctly.
            origins <- rep(seq_along(all.win), lengths(all.win))
            f <- paste0(cons$id, ".", origins)
            summed <- by(cons$weight, INDICES=f, FUN=sum)
            expect_equal(as.numeric(summed), rep(1, length(unique(f))))

            # Responds to the sign.
            sign.list <- list()
            for (s in seq_along(sizes)) {
                sign.list[[s]] <- rbinom(length(data.list[[s]]), 1, 0.5)==1
            }
            cons3 <- mergeWindowsList(data.list, tol=100, sign.list=sign.list)
            ref2 <- mergeWindows(combined.win, tol=100, sign=unlist(sign.list))
            expect_identical(cons3$id, ref2$id)
            expect_identical(cons3$merged, ref2$region)
        }
    }

     # Responds correctly to empty inputs.
     cons.E <- mergeWindowsList(list(all.win[[1]], all.win[[2]][0]), tol=100)
     ref <- mergeWindows(all.win[[1]], tol=100)
     expect_identical(cons.E$id, ref$id)
     expect_identical(cons.E$merged, ref$region)

     cons.E <- mergeWindowsList(list(all.win[[1]][0], all.win[[2]]), tol=100)
     ref <- mergeWindows(all.win[[2]], tol=100)
     expect_identical(cons.E$id, ref$id)
     expect_identical(cons.E$merged, ref$region)

     cons.E <- mergeWindowsList(list(all.win[[1]][0], all.win[[2]][0]), tol=100)
     expect_identical(cons.E$id, integer(0))
     expect_identical(length(cons.E$region), 0L)
     expect_identical(cons.E$weight, numeric(0))

     # Throws an error when expected.
     expect_error(mergeWindowsList(all.win, tol=100, sign.list=list()), "are not identical")
     expect_error(mergeWindowsList(all.win, tol=100, sign.list=vector("list", length(all.win))), "are not identical")
})

set.seed(200001)
test_that("findOverlapsList works correctly", {
    for (nwin in c(5, 50, 500)) {
        for (sizes in list(c(10, 50), c(50, 50), c(100, 50))) {
	    	data.list <- list()
	        regions <- generateWindows(chromos, winsize=1000, nwin=5)
            for (s in seq_along(sizes)) {
	        	windows <- generateWindows(chromos, winsize=sizes[s], nwin=nwin)
	        	counts <- matrix(0, ncol=1, nrow=length(windows))
	        	data.list[[s]] <- SummarizedExperiment(list(counts=counts), rowRanges=windows)
	        }

            all.win <- lapply(data.list, rowRanges)
            combined.win <- do.call(c, all.win)
            ref <- findOverlaps(regions, combined.win)

	        cons <- findOverlapsList(data.list, ref=regions)
            expect_identical(cons$olap, ref)
            expect_true(all(combined.win==cons$ranges))

            # Works with GRanges objects.
            cons2 <- findOverlapsList(all.win, ref=regions)
            expect_identical(cons, cons2)

            # Checking that we respond to arguments.
            ref <- findOverlaps(regions, combined.win, type="within")
            cons2 <- findOverlapsList(data.list, ref=regions, type="within")
            expect_identical(cons2$olap, ref)
            expect_identical(cons$ranges, cons2$ranges)

            # Weights are computed correctly.
            origins <- rep(seq_along(all.win), lengths(all.win))
            f <- paste0(queryHits(cons$olap), ".", origins[subjectHits(cons$olap)])
            summed <- by(cons$weight, INDICES=f, FUN=sum)
            expect_identical(as.numeric(summed), rep(1, length(unique(f))))
        }
    }

    # Responds correctly to empty inputs.
    cons.E <- findOverlapsList(list(all.win[[1]], all.win[[2]][0]), ref=regions)
    expect_identical(cons.E$olap, findOverlaps(regions, all.win[[1]]))

    cons.E <- findOverlapsList(list(all.win[[1]][0], all.win[[2]]), ref=regions)
    expect_identical(cons.E$olap, findOverlaps(regions, all.win[[2]]))

    cons.E <- findOverlapsList(list(all.win[[1]][0], all.win[[2]][0]), ref=regions)
    expect_identical(length(cons.E$olap), 0L)
    expect_identical(cons.E$weight, numeric(0))
})

set.seed(200002)
test_that("consolidateTests works correctly", {
    for (sizes in list(c(10, 50), c(50, 50), c(100, 50))) {
        id.list <- result.list <- weight.list <- list()
        for (s in seq_along(sizes)) {
            ns <- sizes[s]
            result.list[[s]] <- data.frame(logFC=rnorm(ns), PValue=rbeta(ns, 1, 3), logCPM=rnorm(ns))
            id.list[[s]] <- sample(20, ns, replace=TRUE)
            weight.list[[s]] <- runif(ns)
        }

        # Combining things for easier calls.
        combo.id <- unlist(id.list)
        combo.res <- do.call(rbind, result.list)
        combo.weight <- unlist(weight.list)

        # Running consolidation.
        cons <- consolidateTests(id.list, result.list, weight.list=weight.list)
        expect_identical(cons, combineTests(combo.id, combo.res, weight=combo.weight))

        cons <- consolidateTests(id.list, result.list, weight.list=NULL)
        expect_identical(cons, combineTests(combo.id, combo.res))

        expect_warning(cons2 <- consolidateTests(id.list, result.list), "should be specified")
        expect_identical(cons2, cons)

        # Alternative function works, with proper re-indexing if necessary.
        ref <- getBestTest(combo.id, combo.res, weight=combo.weight)
        cons <- consolidateTests(id.list, result.list, weight.list=weight.list, FUN=getBestTest, reindex=NULL)
        expect_identical(ref, cons)

        cons <- consolidateTests(id.list, result.list, weight.list=weight.list, FUN=getBestTest)
        expect_identical(ref$best, cons$best$row + c(0L, as.integer(sizes[1]))[cons$best$origin])
        cons$best <- ref$best <- NULL
        expect_identical(ref, cons)

        cons <- consolidateTests(id.list, result.list, weight.list=NULL, FUN=getBestTest)
        ref <- getBestTest(combo.id, combo.res)
        expect_identical(ref$best, cons$best$row + c(0L, as.integer(sizes[1]))[cons$best$origin])
        cons$best <- ref$best <- NULL
        expect_identical(ref, cons)

        # Argument specification works.
        cons <- consolidateTests(id.list, result.list, weight.list=NULL, fc.col="logCPM")
        expect_identical(cons, combineTests(combo.id, combo.res, fc.col="logCPM"))
    }

    # Testing all the error conditions.
    expect_error(consolidateTests(id.list[1], result.list), "should be the same")
    expect_error(consolidateTests(list(integer(0), integer(0)), result.list), "should match")
    expect_error(consolidateTests(id.list[1], result.list, weight.list=weight.list[1]), "should be the same")
    expect_error(consolidateTests(id.list, result.list, weight.list=list(integer(0), integer(0))), "should match")
})

set.seed(200003)
test_that("consolidateOverlaps works correctly", {
   for (np in c(5, 50, 500)) {
        for (sizes in list(c(10, 50), c(50, 50), c(100, 50))) {
        	olap.list <- result.list <- weight.list <- list()
        	for (s in seq_along(sizes)) {
                ns <- sizes[s]
        		result.list[[s]] <- data.frame(logFC=rnorm(ns), PValue=rbeta(ns, 1, 3), logCPM=rnorm(ns))
        		olap.list[[s]] <- Hits(sample(20, np, replace=TRUE), sample(ns, np, replace=TRUE),
                        nLnode=25, # deliberately chosne to be >20, i.e., some regions have no overlapping windows.
                        nRnode=ns, sort.by.query=TRUE)
                weight.list[[s]] <- runif(np)
            }

            # Combining the Hits.
            last <- 0
            all.S <- vector("list", length(sizes))
            for (s in seq_along(sizes)) {
                all.S[[s]] <- subjectHits(olap.list[[s]]) + last
                last <- last + nRnode(olap.list[[s]])
            }
            combo.olap <- Hits(unlist(lapply(olap.list, queryHits)), unlist(all.S),
                nLnode=nLnode(olap.list[[s]]), nRnode=last, sort.by.query=TRUE)
            combo.res <- do.call(rbind, result.list)
            combo.weight <- unlist(weight.list)

            # Running consolidation.
            cons <- consolidateOverlaps(olap.list, result.list, weight.list=weight.list)
            expect_identical(cons, combineOverlaps(combo.olap, combo.res, o.weight=combo.weight))

            cons <- consolidateOverlaps(olap.list, result.list, weight.list=NULL)
            expect_identical(cons, combineOverlaps(combo.olap, combo.res))

            expect_warning(cons2 <- consolidateOverlaps(olap.list, result.list), "should be specified")
            expect_identical(cons2, cons)

            # Alternative function works, with proper re-indexing if necessary.
            ref <- getBestOverlaps(combo.olap, combo.res, o.weight=combo.weight)
            cons <- consolidateOverlaps(olap.list, result.list, weight.list=weight.list, FUN=getBestOverlaps, reindex=NULL)
            expect_identical(ref, cons)
            
            cons <- consolidateOverlaps(olap.list, result.list, weight.list=weight.list, FUN=getBestOverlaps)
            expect_identical(ref$best, cons$best$row + c(0L, as.integer(sizes))[cons$best$origin])
            cons$best <- ref$best <- NULL
            expect_identical(ref, cons)
            
            ref <- getBestOverlaps(combo.olap, combo.res)
            cons <- consolidateOverlaps(olap.list, result.list, weight.list=NULL, FUN=getBestOverlaps)
            expect_identical(ref$best, cons$best$row + c(0L, as.integer(sizes))[cons$best$origin])
            cons$best <- ref$best <- NULL
            expect_identical(ref, cons)

            # Argument specification works.
            cons <- consolidateOverlaps(olap.list, result.list, weight.list=NULL, fc.col="logCPM")
            expect_identical(cons, combineOverlaps(combo.olap, combo.res, fc.col="logCPM"))
        }
   }

    # Testing all the error conditions.
    expect_error(consolidateOverlaps(olap.list[1], result.list), "should be the same")
    expect_error(consolidateOverlaps(olap.list, list(result.list[[1]][0,], result.list[[2]])), "different length")
    expect_error(consolidateOverlaps(olap.list[1], result.list, weight.list=weight.list[1]), "should be the same")
    expect_error(consolidateOverlaps(olap.list, result.list, weight.list=list(integer(0), integer(0))), "should match")
})
