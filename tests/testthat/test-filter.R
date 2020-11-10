# This tests the sensibility of the filterWindows*() functions. In particular,
# we want to make sure that the filter is calculated properly, despite the 
# manipulations of width and prior count.
# library(csaw); library(testthat); source("test-filter.R")

library(edgeR)

test_that("global filtering works correctly with equal widths", {
    windowed <- SummarizedExperiment(assays=SimpleList(counts=matrix(10, 1, 1)),
    	rowRanges=GRanges("chrA", IRanges(1, 1)), colData=DataFrame(totals=1e6, ext=100),
    	metadata=list(final.ext=NA))
   
    # Filter statistic should be zero, as effective length after extension is the same as the bin width.
    binned <- SummarizedExperiment(assays=SimpleList(counts=matrix(10, nrow=10, 1)),
    	rowRanges=GRanges("chrA", IRanges(0:9*100+1, 1:10*100), seqinfo=Seqinfo("chrA", 1000)), 
    	colData=DataFrame(totals=1e6, ext=1), metadata=list(final.ext=NA, spacing=100))
    out <- filterWindowsGlobal(windowed, binned)
    expect_equivalent(out$filter, 0) 
    
    # Effective length after extension is exactly 10-fold less than the bin width.
    binned <- SummarizedExperiment(assays=SimpleList(counts=matrix(100, nrow=10, 1)),
    	rowRanges=GRanges("chrA", IRanges(0:9*1000+1, 1:10*1000), seqinfo=Seqinfo("chrA", 10000)), 
    	colData=DataFrame(totals=1e6, ext=1), metadata=list(final.ext=NA, spacing=1000))
    out <- filterWindowsGlobal(windowed, binned) 
    expect_equivalent(out$filter, 0)

    # Continues to work when the spacing for the bins is not specified.
    metadata(binned)$spacing <- NULL
    out2 <- filterWindowsGlobal(windowed, binned)
    expect_identical(out, out2)

    # Testing out different prior counts.
    out <- filterWindowsGlobal(windowed, binned, prior.count=3.5) 
    expect_equivalent(out$filter, 0)
    out <- filterWindowsGlobal(windowed, binned, prior.count=0.5) 
    expect_equivalent(out$filter, 0)
    
    # Testing what happens when the median is below the number of recorded bins.
    zeroed <- windowed
    assay(zeroed)[1] <- 0
    binned <- SummarizedExperiment(assays=SimpleList(counts=matrix(10, nrow=10, 1)),
    	rowRanges=GRanges("chrA", IRanges(0:9*100+1, 1:10*100), seqinfo=Seqinfo("chrA", 10000)), 
    	colData=DataFrame(totals=1e6, ext=1), metadata=list(final.ext=NA, spacing=100))
    out <- filterWindowsGlobal(zeroed, binned)
    expect_equivalent(out$filter, 0) # Background estimate is also computed from all-zeroes.
    
    # Testing what happens when the median is within the recorded bins, but not quite the median of them.
    binned <- SummarizedExperiment(assays=SimpleList(counts=matrix(10, nrow=10, 1)),
    	rowRanges=GRanges("chrA", IRanges(0:9*100+1, 1:10*100), seqinfo=Seqinfo("chrA", 1100)), 
    	colData=DataFrame(totals=1e6, ext=1), metadata=list(final.ext=NA, spacing=100))
    out <- filterWindowsGlobal(windowed, binned)
    expect_equivalent(out$filter, 0)

    # Testing what happens when you don't specify the background.
    seqinfo(rowRanges(zeroed)) <- Seqinfo("chrA", 100)
    metadata(zeroed)$spacing <- 10
    out <- filterWindowsGlobal(zeroed) 
    expect_equivalent(out$filter, 0) # Should be zero, as both median and count are based on zero's.

    win2 <- windowed
    seqinfo(rowRanges(win2)) <- Seqinfo("chrA", 100)
    metadata(win2)$spacing <- 10
    out <- filterWindowsGlobal(win2)
    expect_equal(out$abundances, aveLogCPM(asDGEList(win2)))
    expect_equal(out$filter, out$abundances - aveLogCPM(DGEList(matrix(0, 1, ncol(win2)), lib.size=win2$totals)))
    
    metadata(win2)$spacing <- 100
    out2 <- filterWindowsGlobal(win2)
    expect_equivalent(out2$filter, 0)
    expect_equivalent(out2$abundances, out$abundances)

    # Works correctly on empty inputs.
    emp <- filterWindowsGlobal(zeroed[0,])
    expect_equal(emp$filter, numeric(0))
    emp <- filterWindowsGlobal(zeroed[0,], binned)
    expect_equal(emp$filter, numeric(0))
    emp <- filterWindowsGlobal(zeroed, binned[0,])
    expect_equal(emp$filter, 0)
})

test_that("global filtering works correctly with unequal widths", {
    intervals <- SummarizedExperiment(assays=SimpleList(counts=matrix(1:100, 100, 1)),
        rowRanges=GRanges("chrA", IRanges(1, width=1:100)), 
        colData=DataFrame(totals=1e6))
   
    binned <- SummarizedExperiment(assays=SimpleList(counts=matrix(100, nrow=10, 1)),
    	rowRanges=GRanges("chrA", IRanges(0:9*100+1, 1:10*100), seqinfo=Seqinfo("chrA", 1000)), 
    	colData=DataFrame(totals=1e6, ext=1), metadata=list(final.ext=NA, spacing=100))

    # All intervals have counts that scale with width, so all filter stats should be zero. 
    out <- filterWindowsGlobal(intervals[1:10,], binned)
    expect_true(all(abs(out$filter) < 1e-10))

    out <- filterWindowsGlobal(intervals, binned)
    expect_true(all(abs(out$filter) < 1e-5)) # a bit of give required, because the interpolation is approximate.
})

test_that("local filtering works correctly", {
    windowed <- SummarizedExperiment(assays=SimpleList(counts=matrix(10, 1, 1)),
    	rowRanges=GRanges("chrA", IRanges(1, 1)), colData=DataFrame(totals=1e6, ext=100),
    	metadata=list(final.ext=NA))

    # Should be zero, as the count/width for the window is subtracted from the background bin.
    binned <- SummarizedExperiment(assays=SimpleList(counts=matrix(20, 1, 1)),
    	rowRanges=GRanges("chrA", IRanges(1, 200)), colData=DataFrame(totals=1e6, ext=1),
    	metadata=list(final.ext=NA))
    out <- filterWindowsLocal(windowed, binned)
    expect_equivalent(out$filter, 0)

    # After subtraction, the background is still 10-times wider, but again this should be zero.
    binned <- SummarizedExperiment(assays=SimpleList(counts=matrix(110, 1, 1)),
    	rowRanges=GRanges("chrA", IRanges(1, 1100)), colData=DataFrame(totals=1e6, ext=1),
    	metadata=list(final.ext=NA))
    out <- filterWindowsLocal(windowed, binned)
    expect_equivalent(out$filter, 0)
  
    # Testing out different prior counts.
    out <- filterWindowsLocal(windowed, binned, prior.count=3.5) 
    expect_equivalent(out$filter, 0)
    out <- filterWindowsLocal(windowed, binned, prior.count=0.5) 
    expect_equivalent(out$filter, 0)
 
    # Another case of subtraction, this time with a different extension for the background.
    binned <- SummarizedExperiment(assays=SimpleList(counts=matrix(110, 1, 1)),
    	rowRanges=GRanges("chrA", IRanges(1, 1001)), colData=DataFrame(totals=1e6, ext=100),
    	metadata=list(final.ext=NA))
    out <- filterWindowsLocal(windowed, binned)
    expect_equivalent(out$filter, 0)
    
    # More subtraction.
    binned <- SummarizedExperiment(assays=SimpleList(counts=matrix(100, 1, 1)),
    	rowRanges=GRanges("chrA", IRanges(1, 1000)), colData=DataFrame(totals=1e6, ext=1),
    	metadata=list(final.ext=1))
    out <- filterWindowsLocal(windowed, binned)
    expect_equivalent(out$filter, 0)

    # Works correctly on empty inputs.
    emp <- filterWindowsLocal(windowed[0,], binned[0,])
    expect_equal(emp$filter, numeric(0))
    expect_error(filterWindowsLocal(windowed, binned[0,]), "same length")
})

test_that("control-based filtering works correctly", {
    windowed <- SummarizedExperiment(assays=SimpleList(counts=matrix(10, 1, 1)),
    	rowRanges=GRanges("chrA", IRanges(1, 1)), colData=DataFrame(totals=1e6, ext=100),
    	metadata=list(final.ext=NA))

    countered <- windowed
    expect_warning(out <- filterWindowsControl(windowed, countered), "not specified")
    expect_equivalent(out$filter, 0)
    expect_warning(out <- filterWindowsControl(windowed, countered, prior.count=5), "not specified")
    expect_equivalent(out$filter, 0)

    # Also seeing what happens when the library size of the control changes.
    countered2 <- countered
    assay(countered2)[1] <- 20
    countered2$totals <- 2e6
    expect_warning(out <- filterWindowsControl(windowed, countered2))

    # With normalization; in this case, a trivial scaling.
    binned.chip <- SummarizedExperiment(assays=SimpleList(counts=matrix(100, 1, 1)),
        rowRanges=GRanges("chrA", IRanges(1, 1000)), colData=DataFrame(totals=1e6, ext=1),
        metadata=list(final.ext=NA))
    binned.con <- SummarizedExperiment(assays=SimpleList(counts=matrix(100, 1, 1)),
        rowRanges=GRanges("chrA", IRanges(1, 1000)), colData=DataFrame(totals=1e6, ext=1),
        metadata=list(final.ext=NA))

    sc.info <- scaleControlFilter(binned.chip, binned.con)
    expect_equal(sc.info$scale, 1)
    expect_identical(sc.info$data.totals, windowed$totals)
    expect_identical(sc.info$back.totals, countered$totals)
    
    out <- filterWindowsControl(windowed, countered, scale.info=sc.info)
    expect_equivalent(out$filter, 0)

    # More effortful normalization, assuming undersampling in control.
    countered2 <- countered
    assay(countered2)[1] <- 5 
    binned.con2 <- binned.con
    assay(binned.con2)[1] <- 50

    sc.info <- scaleControlFilter(binned.chip, binned.con2)
    expect_equal(sc.info$scale, 50/100)
    expect_identical(sc.info$data.totals, windowed$totals)
    expect_identical(sc.info$back.totals, countered2$totals)

    out <- filterWindowsControl(windowed, countered2, scale.info=sc.info)
    expect_equivalent(out$filter, 0)
    out <- filterWindowsControl(windowed, countered2, scale.info=sc.info, prior.count=5)
    expect_equivalent(out$filter, 0)

    # Checking that it is unhappy when the regions are of different size.
    countered2 <- SummarizedExperiment(assays=SimpleList(counts=matrix(20, 1, 1)),
    	rowRanges=GRanges("chrA", IRanges(1, 200)), colData=DataFrame(totals=1e6, ext=1),
    	metadata=list(final.ext=NA))
    expect_warning(out <- filterWindowsControl(windowed, countered2))
    expect_equal(out$filter, 0)

    # Works correctly on empty inputs.
    expect_warning(emp <- filterWindowsControl(windowed[0,], countered[0,]))
    expect_equal(emp$filter, numeric(0))
    expect_error(filterWindowsControl(windowed, countered[0,]), "same length")
})

test_that("proportional filtering works as expected", {
    multi.win <- SummarizedExperiment(assays=list(counts=matrix(11:20, 10, 1)),
    	rowRanges=GRanges("chrA", IRanges(1:10, 1:10), seqinfo=Seqinfo("chrA", 1000)), 
    	colData=DataFrame(totals=1e6, ext=100), metadata=list(final.ext=NA, spacing=1))

    # works if not all windows are available, assuming the lost windows are lower abundance.
    out <- filterWindowsProportion(multi.win)$filter
    expect_equal(tail(out,1), 1)
    expect_true(all(diff(out) > 0))
    expect_true(all(out > 0.99))

    # Still works upon changes to the chromosome size.
    seqinfo(rowRanges(multi.win)) <- Seqinfo("chrA", 10)
    out <- filterWindowsProportion(multi.win)$filter
    expect_equal(out, 1:10/10)
    expect_true(all(diff(out) > 0))

    # Works correctly on empty inputs.
    emp <- filterWindowsProportion(multi.win[0,])
    expect_equal(emp$filter, numeric(0))
})
