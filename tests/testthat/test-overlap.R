# This checks the overlap summarization functions, relative to the expected values.
# library(testthat); library(csaw); source("setup.R"); source("test-overlap.R")

chromos <- c(A=1000, B=2000)

set.seed(130000)
test_that("combineOverlaps works correctly", {
    for (nreg in c(2, 10)) {
        for (nwin in c(1, 10, 100)) {
            regions <- generateWindows(chromos, nreg, 500)
            windows <- generateWindows(chromos, nwin, 50)

            olap <- findOverlaps(regions, windows)
            ns <- length(windows)
            tab <- data.frame(logFC=rnorm(ns), PValue=rbeta(ns, 1, 3), logCPM=rnorm(ns))

            # Straight-up comparison to combineTests, after discarding all NA's.
            output <- combineOverlaps(olap, tab)
            refstats <- combineTests(queryHits(olap), tab[subjectHits(olap),])
        	refstats$rep.test <- subjectHits(olap)[refstats$rep.test]
            expect_identical(output[!is.na(output$PValue),], refstats)

        	# Testing with weights.
            test.weight <- runif(ns)
            output <- combineOverlaps(olap, tab, i.weight=test.weight)
            refstats <- combineTests(queryHits(olap), tab[subjectHits(olap),], weight=test.weight[subjectHits(olap)])
        	refstats$rep.test <- subjectHits(olap)[refstats$rep.test]
            expect_identical(output[!is.na(output$PValue),], refstats)

        	# More weight testing, where o.weight is constructed from the weight for each i.weight.
            output2 <- combineOverlaps(olap, tab, o.weight=test.weight[subjectHits(olap)])
            expect_identical(output, output2)
        }
    }

    # Testing with empty inputs.
    out <- combineOverlaps(Hits(), data.frame(logFC=numeric(0), PValue=numeric(0), logCPM=numeric(0)))
    expect_identical(nrow(out), 0L)
    expect_identical(out$PValue, numeric(0)) 
})

set.seed(130001)
test_that("getBestOverlaps works correctly", {
    for (nreg in c(2, 10)) {
        for (nwin in c(1, 10, 100)) {
            regions <- generateWindows(chromos, nreg, 500)
            windows <- generateWindows(chromos, nwin, 50)
            
        	olap <- findOverlaps(regions, windows)
        	ns <- length(windows)
        	tab <- data.frame(logFC=rnorm(ns), PValue=rbeta(ns, 1, 3), logCPM=rnorm(ns))
        
        	output <- getBestOverlaps(olap, tab)
        	refstats <- getBestTest(queryHits(olap), tab[subjectHits(olap),])
        	refstats$rep.test <- subjectHits(olap)[refstats$rep.test]
        	expect_identical(output[!is.na(output$PValue),], refstats)
        
        	# Testing with weights.
        	test.weight <- runif(ns)
        	output <- getBestOverlaps(olap, tab, i.weight=test.weight)
        	refstats <- getBestTest(queryHits(olap), tab[subjectHits(olap),], weight=test.weight[subjectHits(olap)])
        	refstats$rep.test <- subjectHits(olap)[refstats$rep.test]
        	expect_identical(output[!is.na(output$PValue),], refstats)
        
        	# More weight testing.
        	output2 <- getBestOverlaps(olap, tab, o.weight=test.weight[subjectHits(olap)])
        	expect_identical(output, output2)
        }
    }

    # Testing with empty inputs.
    out <- getBestOverlaps(Hits(), data.frame(logFC=numeric(0), PValue=numeric(0), logCPM=numeric(0)))
    expect_identical(nrow(out), 0L)
    expect_identical(out$PValue, numeric(0)) 
})

set.seed(130002)
test_that("empiricalOverlaps works correctly", {
    for (nreg in c(2, 10)) {
        for (nwin in c(1, 10, 100)) {
            regions <- generateWindows(chromos, nreg, 500)
            windows <- generateWindows(chromos, nwin, 50)
            
        	olap <- findOverlaps(regions, windows)
        	ns <- length(windows)
        	tab <- data.frame(logFC=rnorm(ns), PValue=rbeta(ns, 1, 3), logCPM=rnorm(ns))
        
        	# Straight-up comparison to empiricalFDR, after discarding all NA's.
        	output <- empiricalOverlaps(olap, tab)
        	refstats <- empiricalFDR(queryHits(olap), tab[subjectHits(olap),])
        	refstats$rep.test <- subjectHits(olap)[refstats$rep.test]
        	expect_identical(output[!is.na(output$PValue),], refstats)
        
        	# Testing with weights.
        	test.weight <- runif(ns)
        	output <- empiricalOverlaps(olap, tab, i.weight=test.weight)
        	refstats <- empiricalFDR(queryHits(olap), tab[subjectHits(olap),], weight=test.weight[subjectHits(olap)])
        	refstats$rep.test <- subjectHits(olap)[refstats$rep.test]
        	expect_identical(output[!is.na(output$PValue),], refstats)
        
        	# More weight testing, where o.weight is constructed from the weight for each i.weight.
        	output2 <- empiricalOverlaps(olap, tab, o.weight=test.weight[subjectHits(olap)])
        	expect_identical(output, output2)
        }
    }

    # Testing with empty inputs.
    out <- empiricalOverlaps(Hits(), data.frame(logFC=numeric(0), PValue=numeric(0), logCPM=numeric(0)))
    expect_identical(nrow(out), 0L)
    expect_identical(out$PValue, numeric(0)) 
})

set.seed(1300021)
test_that("mixedOverlaps works correctly", {
    for (nreg in c(2, 10)) {
        for (nwin in c(1, 10, 100)) {
            regions <- generateWindows(chromos, nreg, 500)
            windows <- generateWindows(chromos, nwin, 50)

            olap <- findOverlaps(regions, windows)
            ns <- length(windows)
            tab <- data.frame(logFC=rnorm(ns), PValue=rbeta(ns, 1, 3), logCPM=rnorm(ns))

            # Straight-up comparison to combineTests, after discarding all NA's.
            output <- mixedOverlaps(olap, tab)
            refstats <- mixedClusters(queryHits(olap), tab[subjectHits(olap),])
        	refstats$rep.up.test <- subjectHits(olap)[refstats$rep.up.test]
        	refstats$rep.down.test <- subjectHits(olap)[refstats$rep.down.test]
            expect_identical(output[!is.na(output$PValue),], refstats)

        	# Testing with weights.
            test.weight <- runif(ns)
            output <- mixedOverlaps(olap, tab, i.weight=test.weight)
            refstats <- mixedClusters(queryHits(olap), tab[subjectHits(olap),], weight=test.weight[subjectHits(olap)])
        	refstats$rep.up.test <- subjectHits(olap)[refstats$rep.up.test]
        	refstats$rep.down.test <- subjectHits(olap)[refstats$rep.down.test]
            expect_identical(output[!is.na(output$PValue),], refstats)

        	# More weight testing, where o.weight is constructed from the weight for each i.weight.
            output2 <- mixedOverlaps(olap, tab, o.weight=test.weight[subjectHits(olap)])
            expect_identical(output, output2)
        }
    }

    # Testing with empty inputs.
    out <- mixedOverlaps(Hits(), data.frame(logFC=numeric(0), PValue=numeric(0), logCPM=numeric(0)))
    expect_identical(nrow(out), 0L)
    expect_identical(out$PValue, numeric(0)) 
})

set.seed(130003)
test_that("summitOverlaps works correctly", {
    for (nreg in c(2, 10)) {
        for (nwin in c(1, 10, 100)) {
            regions <- generateWindows(chromos, nreg, 500)
            windows <- generateWindows(chromos, nwin, 50)

        	olap <- findOverlaps(regions, windows)
        	ns <- length(windows)
        	tab <- data.frame(logFC=rnorm(ns), PValue=rbeta(ns, 1, 3), logCPM=rnorm(ns))
        	output <- getBestOverlaps(olap, tab)
        
        	# Checking summit calls.
        	re.weight <- summitOverlaps(olap, output$rep.test)
        	best.win <- output$rep.test[queryHits(olap)]
        	is.summit <- !is.na(best.win) & best.win==subjectHits(olap)
        	re.weight2a <- summitOverlaps(olap, o.summit=is.summit)
        	re.weight2b <- summitOverlaps(olap, o.summit=which(is.summit))
        	expect_identical(re.weight, re.weight2a)
            expect_identical(re.weight, re.weight2b)
        
        	isummits <- rbinom(ns, 1, 0.1)==1L
        	re.weight3 <- summitOverlaps(olap, o.summit=isummits[subjectHits(olap)])
        	re.weight4 <- summitOverlaps(olap, i.summit=isummits)
        	expect_identical(re.weight3, re.weight4)
        
        	# Checking the core upweightSummit machinery itself.
        	by.region <- split(is.summit, queryHits(olap))	
        	nu.weight <- lapply(by.region, FUN=function(x) {
        		N <- length(x)
        		output <- rep(1, N)
        		output[x] <- N/sum(x)
        		output	
        	})
            if (length(nu.weight)) {
                expect_identical(re.weight, unlist(nu.weight, use.names=FALSE))
            } else {
                expect_identical(re.weight, numeric(0))
            }
        }
    }

    # Testing with empty inputs.
    out <- summitOverlaps(Hits(), data.frame(logFC=numeric(0), PValue=numeric(0), logCPM=numeric(0)))
    expect_identical(out, numeric(0))
})
