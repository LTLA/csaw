# This tests the consolidateSizes function. As that function is basically written in fairly easy R, the
# test function below doesn't really do much except repeat the function itself.
# library(csaw); library(testthat); source("test-cons.R")

source("simsam.R")
chromos <- c(chrA=12332, chrB=34892)

set.seed(200000)
test_that("consolidateWindows works correctly with mergeWindows", {
    for (nwin in c(5, 50, 500)) {
        for (sizes in list(c(10, 50), c(50, 50), c(100, 50))) {
	        data.list <- list()
        	for (s in seq_along(sizes)) {
        		windows <- generateWindows(chromos, winsize=sizes[s], nwin=nwin)
        		counts <- matrix(0, ncol=1, nrow=length(windows))
        		data.list[[s]] <- SummarizedExperiment(list(counts=counts), rowRanges=windows)
        	}

        	# Running consolidation and checking the output.
        	expect_warning(cons <- consolidateWindows(data.list), "default")
            all.win <- lapply(data.list, rowRanges)
            ref <- mergeWindows(do.call(c, all.win), tol=100)
            expect_identical(unlist(cons$id), ref$id)
            expect_identical(cons$region, ref$region)

            cons2 <- consolidateWindows(data.list, merge.args=list(tol=50))
            ref2 <- mergeWindows(do.call(c, all.win), tol=50)
            expect_identical(unlist(cons2$id), ref2$id)
            expect_identical(cons2$region, ref2$region)

            expect_warning(cons3 <- consolidateWindows(all.win), "default") # Works with GRanges objects.
            expect_identical(cons, cons3)
            expect_warning(cons3 <- consolidateWindows(all.win, merge.args=list(tol=100)), NA) # shuts up the warning.
            expect_identical(cons, cons3)

            # Weights are computed correctly.
            for (i in seq_along(cons$weight)) {
                curid <- cons$id[[i]]
                summed <- by(cons$weight[[i]], INDICES=curid, FUN=sum)
                expect_identical(as.numeric(summed), rep(1, length(unique(curid))))
            }

            cons4 <- consolidateWindows(all.win, equiweight=FALSE, merge.args=list(tol=100))
            expect_identical(cons[c("id", "region")], cons4[c("id", "region")])
            expect_identical(cons4$weight, NULL)

            # Responds to the sign.
            sign.list <- list()
            for (s in seq_along(sizes)) {
                sign.list[[s]] <- rbinom(length(data.list[[s]]), 1, 0.5)==1
            }
            cons5 <- consolidateWindows(data.list, merge.args=list(tol=100), sign.list=sign.list)
            ref3 <- mergeWindows(do.call(c, all.win), tol=100, sign=unlist(sign.list))
            expect_identical(unlist(cons5$id), ref3$id)
            expect_identical(cons5$region, ref3$region)
        }
    }

     # Responds correctly to empty inputs.
     cons.E <- consolidateWindows(list(all.win[[1]], all.win[[2]][0]), merge.args=list(tol=100))
     ref <- mergeWindows(all.win[[1]], tol=100)
     expect_identical(cons.E$id[[1]], ref$id)
     expect_identical(cons.E$id[[2]], integer(0))
     expect_identical(cons.E$region, ref$region)

     cons.E <- consolidateWindows(list(all.win[[1]][0], all.win[[2]]), merge.args=list(tol=100))
     ref <- mergeWindows(all.win[[2]], tol=100)
     expect_identical(cons.E$id[[1]], integer(0))
     expect_identical(cons.E$id[[2]], ref$id)
     expect_identical(cons.E$region, ref$region)

     cons.E <- consolidateWindows(list(all.win[[1]][0], all.win[[2]][0]), merge.args=list(tol=100))
     expect_identical(lengths(cons.E$id), integer(2))
     expect_identical(length(cons.E$region), 0L)
     expect_identical(lengths(cons.E$weight), integer(2))

     # Throws an error when expected.
     expect_error(consolidateWindows(all.win, merge.args=list(tol=100), sign.list=list()), "are not identical")
     expect_error(consolidateWindows(all.win, merge.args=list(tol=100), sign.list=vector("list", length(all.win))), "are not identical")
})

# Defining a function that unpacks and sorts the query/subject hits.
unpack <- function(cons.olap) {
    comp.Q <- unlist(lapply(cons.olap, queryHits))
    o <- order(comp.Q)
    comp.Q <- comp.Q[o]

    comp.S <- lapply(cons.olap, subjectHits)
    for (x in 2:length(comp.S)) {
        comp.S[[x]] <- comp.S[[x]] + nRnode(cons.olap[[x-1]])
    }
    comp.S <- unlist(comp.S)[o]
    return(list(Q=comp.Q, S=comp.S))
}

set.seed(200001)
test_that("consolidateWindows works correctly with findOverlaps", {
    for (nwin in c(5, 50, 500)) {
        for (sizes in list(c(10, 50), c(50, 50), c(100, 50))) {
	    	data.list <- list()
	        regions <- generateWindows(chromos, winsize=1000, nwin=5)
            for (s in seq_along(sizes)) {
	        	windows <- generateWindows(chromos, winsize=sizes[s], nwin=nwin)
	        	counts <- matrix(0, ncol=1, nrow=length(windows))
	        	data.list[[s]] <- SummarizedExperiment(list(counts=counts), rowRanges=windows)
	        }

            # Running consolidation and checking the output.
	        cons <- consolidateWindows(data.list, region=regions)
            unpacked <- unpack(cons$olap)

            all.win <- lapply(data.list, rowRanges)
            ref <- findOverlaps(regions, do.call(c, all.win))
            expect_identical(unpacked$Q, queryHits(ref))
            expect_identical(unpacked$S, subjectHits(ref))

            # Checking that we respond to arguments.
            cons2 <- consolidateWindows(data.list, region=regions, overlap.args=list(type="within"))
            unpacked2 <- unpack(cons2$olap)

            ref <- findOverlaps(regions, do.call(c, all.win), type="within")
            expect_identical(unpacked2$Q, queryHits(ref))
            expect_identical(unpacked2$S, subjectHits(ref))

            cons3 <- consolidateWindows(all.win, region=regions) # Works with GRanges objects.
            expect_identical(cons, cons3)

            # Weights are computed correctly.
            for (i in seq_along(cons$weight)) {
                curid <- queryHits(cons$olap[[i]])
                summed <- by(cons$weight[[i]], INDICES=curid, FUN=sum)
                expect_identical(as.numeric(summed), rep(1, length(unique(curid))))
            }

            cons4 <- consolidateWindows(all.win, region=regions, equiweight=FALSE)
            expect_identical(cons[c("olap", "region")], cons4[c("olap", "region")])
            expect_identical(cons4$weight, NULL)
        }
    }

    # Responds correctly to empty inputs.
    cons.E <- consolidateWindows(list(all.win[[1]], all.win[[2]][0]), region=regions)
    expect_identical(cons.E$olap[[1]], findOverlaps(regions, all.win[[1]]))
    expect_identical(length(cons.E$olap[[2]]), 0L)

    cons.E <- consolidateWindows(list(all.win[[1]][0], all.win[[2]]), region=regions)
    expect_identical(cons.E$olap[[2]], findOverlaps(regions, all.win[[2]]))
    expect_identical(length(cons.E$olap[[1]]), 0L)

    cons.E <- consolidateWindows(list(all.win[[1]][0], all.win[[2]][0]), region=regions)
    expect_identical(lengths(cons.E$olap), integer(2))
    expect_identical(length(cons.E$region), 0L)
    expect_identical(lengths(cons.E$weight), integer(2))
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
        expect_identical(ref$best, cons$best$row + c(0L, as.integer(sizes))[cons$best$origin])
        cons$best <- ref$best <- NULL
        expect_identical(ref, cons)

        cons <- consolidateTests(id.list, result.list, weight.list=NULL, FUN=getBestTest)
        ref <- getBestTest(combo.id, combo.res)
        expect_identical(ref$best, cons$best$row + c(0L, as.integer(sizes))[cons$best$origin])
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
