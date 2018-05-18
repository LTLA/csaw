# This tests the sensibility of the filterWindows() function. In particular,
# we want to make sure that the filter is calculated properly, despite the 
# manipulations of width and prior count.

suppressWarnings(suppressPackageStartupMessages(require(csaw)))
suppressWarnings(suppressPackageStartupMessages(require(edgeR)))

comp <- function(..., tol=1e-8) {
	out <- filterWindows(...)
	stopifnot(all(abs(out$filter) < tol))
	out
}

####################################################################################################
# Checking that scaledAverage mimics aveLogCPM with prior.count=1.

set.seed(100)
se <- SummarizedExperiment(list(counts=matrix(rnbinom(1000, mu=100, size=10), ncol=10, nrow=100)))
se$totals <- colSums(assay(se))

y <- asDGEList(se)
stopifnot(identical(y$samples$lib.size, se$totals))
ref <- aveLogCPM(y)
out <- scaledAverage(se, scale=1)
stopifnot(length(ref)==length(out))
stopifnot(all(abs(ref-out) < 1e-10))
head(out)

se$norm.factors <- runif(10, 0.5, 1.5)
y <- asDGEList(se)
stopifnot(identical(y$samples$norm.factors, se$norm.factors))
ref <- aveLogCPM(y)
out <- scaledAverage(se, scale=1)
stopifnot(length(ref)==length(out))
stopifnot(all(abs(ref-out) < 1e-10))
head(out)

assay(se, "offset") <- matrix(rnorm(1000), ncol=10, nrow=100) # Neither function should care about the offset.
y <- asDGEList(se)
stopifnot(!is.null(y$offset))    
ref <- aveLogCPM(y)
out <- scaledAverage(se, scale=1)
stopifnot(length(ref)==length(out))
stopifnot(all(abs(ref-out) < 1e-10))
head(out)

ref <- aveLogCPM(y, dispersion=0.1)
out <- scaledAverage(se, scale=1, dispersion=0.1)
stopifnot(length(ref)==length(out))
stopifnot(all(abs(ref-out) < 1e-10))
head(out)

# Checking proper behaviour of the scaling.

se <- SummarizedExperiment(list(counts=matrix(rnbinom(100, mu=100, size=10), ncol=1, nrow=100)))
se$totals <- colSums(assay(se))    
ref <- scaledAverage(se, scale=1)

se2 <- se
assay(se2) <- assay(se) * 2
out2 <- scaledAverage(se2, scale=2)
stopifnot(all(abs(ref - out2) < 1e-10))
head(out2)

scalar <- runif(nrow(se), 1, 5)
se1 <- se
assay(se1) <- scalar * assay(se)
out1 <- scaledAverage(se1, scale=scalar)
stopifnot(all(abs(ref - out1) < 1e-10))
head(out1)

# Checking proper behaviour with invalid inputs.

stopifnot(all(is.na(scaledAverage(se, scale=-1))))
stopifnot(all(is.infinite(scaledAverage(se, scale=0))))
flucscale <- rep(1, nrow(y))
flucscale[1] <- 0
flucscale[2] <- -1

out <- scaledAverage(se, scale=flucscale)
ref <- aveLogCPM(asDGEList(se))
stopifnot(is.infinite(out[1]))
stopifnot(is.na(out[2]))
stopifnot(all(abs(out-ref)[-c(1,2)] < 1e-10))
head(out)

####################################################################################################
# Matching up global filtering.

windowed <- SummarizedExperiment(assays=SimpleList(counts=matrix(10, 1, 1)),
	rowRanges=GRanges("chrA", IRanges(1, 1)), colData=DataFrame(totals=1e6, ext=100),
	metadata=list(final.ext=NA))

binned <- SummarizedExperiment(assays=SimpleList(counts=matrix(10, nrow=10, 1)),
	rowRanges=GRanges("chrA", IRanges(0:9*100+1, 1:10*100), seqinfo=Seqinfo("chrA", 1000)), 
	colData=DataFrame(totals=1e6, ext=1), metadata=list(final.ext=NA, spacing=100))
comp(windowed, binned, type="global")

binned <- SummarizedExperiment(assays=SimpleList(counts=matrix(100, nrow=10, 1)),
	rowRanges=GRanges("chrA", IRanges(0:9*1000+1, 1:10*1000), seqinfo=Seqinfo("chrA", 10000)), 
	colData=DataFrame(totals=1e6, ext=1), metadata=list(final.ext=NA, spacing=1000))
comp(windowed, binned, type="global") 

# Testing what happens when the median is below the number of recorded bins.

zeroed <- windowed
assay(zeroed)[1] <- 0
binned <- SummarizedExperiment(assays=SimpleList(counts=matrix(10, nrow=10, 1)),
	rowRanges=GRanges("chrA", IRanges(0:9*100+1, 1:10*100), seqinfo=Seqinfo("chrA", 10000)), 
	colData=DataFrame(totals=1e6, ext=1), metadata=list(final.ext=NA, spacing=100))
comp(zeroed, binned, type="global")

# Testing what happens when the median is within the recorded bins, but not quite the median of them.

binned <- SummarizedExperiment(assays=SimpleList(counts=matrix(10, nrow=10, 1)),
	rowRanges=GRanges("chrA", IRanges(0:9*100+1, 1:10*100), seqinfo=Seqinfo("chrA", 1100)), 
	colData=DataFrame(totals=1e6, ext=1), metadata=list(final.ext=NA, spacing=100))
comp(windowed, binned, type="global")

# Testing what happens when you don't specify the bin.

seqinfo(rowRanges(zeroed)) <- Seqinfo("chrA", 100)
metadata(zeroed)$spacing <- 10
comp(zeroed, type="global") # Should be zero, as both median and count are based on zero's.

win2 <- windowed
seqinfo(rowRanges(win2)) <- Seqinfo("chrA", 100)
metadata(win2)$spacing <- 10
comp(win2, type="global", tol=Inf) # Should be large, as median is based on zero's.

metadata(win2)$spacing <- 100
comp(win2, type="global") # Should be zero, as median is now based on the window itself (only window in genome).

####################################################################################################
# Matching up local filtering (both count and width subtracted from the background).

binned <- SummarizedExperiment(assays=SimpleList(counts=matrix(20, 1, 1)),
	rowRanges=GRanges("chrA", IRanges(1, 200)), colData=DataFrame(totals=1e6, ext=1),
	metadata=list(final.ext=NA))
comp(windowed, binned, type="local")

binned <- SummarizedExperiment(assays=SimpleList(counts=matrix(110, 1, 1)),
	rowRanges=GRanges("chrA", IRanges(1, 1100)), colData=DataFrame(totals=1e6, ext=1),
	metadata=list(final.ext=NA))
comp(windowed, binned, type="local")

binned <- SummarizedExperiment(assays=SimpleList(counts=matrix(110, 1, 1)),
	rowRanges=GRanges("chrA", IRanges(1, 1001)), colData=DataFrame(totals=1e6, ext=100),
	metadata=list(final.ext=NA))
comp(windowed, binned, type="local")

binned <- SummarizedExperiment(assays=SimpleList(counts=matrix(100, 1, 1)),
	rowRanges=GRanges("chrA", IRanges(1, 1000)), colData=DataFrame(totals=1e6, ext=1),
	metadata=list(final.ext=1))
comp(windowed, binned, type="local")

####################################################################################################
# Matching up control filtering.

countered <- windowed
suppressWarnings(comp(windowed, countered, type="control"))
suppressWarnings(comp(windowed, countered, type="control", prior.count=5))

# With normalization.

binned.chip <- SummarizedExperiment(assays=SimpleList(counts=matrix(100, 1, 1)),
	rowRanges=GRanges("chrA", IRanges(1, 1000)), colData=DataFrame(totals=1e6, ext=1),
	metadata=list(final.ext=NA))
binned.con <- SummarizedExperiment(assays=SimpleList(counts=matrix(100, 1, 1)),
	rowRanges=GRanges("chrA", IRanges(1, 1000)), colData=DataFrame(totals=1e6, ext=1),
	metadata=list(final.ext=NA))

sc.info <- scaleControlFilter(binned.chip, binned.con)
stopifnot(isTRUE(all.equal(sc.info$scale, 1)))
stopifnot(identical(sc.info$data.totals, windowed$totals))
stopifnot(identical(sc.info$back.totals, countered$totals))

comp(windowed, countered, type="control", scale.info=sc.info)$filter

# More effortful normalization.

assay(countered)[1] <- 5 # assuming undersampling in control.
assay(binned.con)[1] <- 50

sc.info <- scaleControlFilter(binned.chip, binned.con)
stopifnot(isTRUE(all.equal(sc.info$scale, 50/100)))
stopifnot(identical(sc.info$data.totals, windowed$totals))
stopifnot(identical(sc.info$back.totals, countered$totals))

comp(windowed, countered, type="control", scale.info=sc.info)
comp(windowed, countered, type="control", scale.info=sc.info, prior.count=5)

# Also seeing what happens when the library size of the control changes.

assay(countered)[1] <- 20
countered$totals <- 2e6
suppressWarnings(comp(windowed, countered, type="control"))

####################################################################################################
# Matching up proportional changes.

multi.win <- SummarizedExperiment(assays=SimpleList(counts=matrix(11:20, 10, 1)),
	rowRanges=GRanges("chrA", IRanges(1:10, 1:10), seqinfo=Seqinfo("chrA", 1000)), 
	colData=DataFrame(totals=1e6, ext=100), metadata=list(final.ext=NA, spacing=1))
out <- filterWindows(multi.win, type="proportion")$filter
stopifnot(tail(out,1)==1)
stopifnot(diff(out) > 0)
out

seqinfo(rowRanges(multi.win)) <- Seqinfo("chrA", 10)
out <- filterWindows(multi.win, type="proportion")$filter
stopifnot(tail(out,1)==1)
stopifnot(diff(out) > 0)
out

####################################################################################################
# End.

