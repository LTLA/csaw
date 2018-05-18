# This tests the correctness of the normOffsets functions.

suppressPackageStartupMessages(library(csaw))
suppressPackageStartupMessages(library(edgeR))

###################################################

# First we test for scaling normalization.

set.seed(1000)
data <- SummarizedExperiment(list(counts=matrix(rpois(10000, lambda=10), ncol=10)))
data$totals <- rpois(10, lambda=10000)

nf <- normFactors(data, se.out=FALSE)
ref <- calcNormFactors(DGEList(assay(data), lib.size=data$totals), doWeighting=FALSE)$samples$norm.factors
stopifnot(identical(nf, ref))
ref

data2 <- normFactors(data, se.out=TRUE)
stopifnot(identical(data2$norm.factors, ref))
stopifnot(identical(data2$totals, data$totals))

# Checking that overwriting se.out works.

data3 <- data
assay(data3) <- assay(data3)*runif(nrow(data3))
data3b <- normFactors(data, se.out=data3)
stopifnot(identical(assay(data3b), assay(data3)))
stopifnot(identical(data3b$norm.factors, ref))

data3$totals <- rpois(10, lambda=10000)
try(normFactors(data, se.out=data3)) # should chuck an error.

###################################################

# Additional fixed-value checks; testing what happens when we add some undersampling.

n <- 1000L
mu1 <- rep(10, n)
mu2 <- mu1
mu2[1:100] <- 100
undersamp <- sum(mu1)/sum(mu2)
mu2 <- mu2*undersamp
counts <- cbind(rnbinom(n, mu=mu1, size=20), rnbinom(n, mu=mu2, size=20))

actual.lib.size <- c(sum(mu1), sum(mu2))
normFactors(counts, lib.sizes=actual.lib.size)
normFactors(counts, logratioTrim=0.4, lib.sizes=actual.lib.size)
normFactors(counts, sumTrim=0.3, lib.size=actual.lib.size)
sqrt(c(1/undersamp, undersamp)) # True values, looks pretty good.

# Testing what happens with weighting, after adding some high-abundance DB regions. 

n <- 100000
blah <- matrix(rnbinom(2*n, mu=10, size=20), ncol=2)
tospike <- 10000
blah[1:tospike,1] <- rnbinom(tospike, mu=1000, size=20)
blah[1:tospike,2] <- rnbinom(tospike, mu=2000, size=20)
full.lib.size <- colSums(blah)

true.value <- 1/full.lib.size
true.value <- true.value/exp(mean(log(true.value)))
true.value

normFactors(blah, lib.sizes=full.lib.size)
normFactors(blah, weighted=TRUE, lib.sizes=full.lib.size) # less accurate.

###################################################

# Now we throw in some tests for loess normalization.

set.seed(1000)
means <- 2^runif(1000)
data <- SummarizedExperiment(list(counts=matrix(rpois(10000, lambda=means), ncol=10)))
data$totals <- rpois(10, lambda=10000)

offs <- normOffsets(data, type="loess", se.out=FALSE)
data <- normOffsets(data, type="loess", se.out=TRUE)
stopifnot(identical(dim(offs), dim(assay(data, "offset"))))
stopifnot(all(abs(offs-assay(data, "offset")) < 1e-8))
head(offs)

try(normOffsets(data, type="loess", se.out=data))

# Reference calculation, after subtracting the reference 'ab' from the observed values.

lib.sizes <- data$totals
mat <- assay(data)

cont.cor <- 0.5
cont.cor.scaled <- cont.cor * lib.sizes/mean(lib.sizes)
ab <- aveLogCPM(mat, lib.size=lib.sizes, prior.count=cont.cor)/log2(exp(1))

ref <- matrix(0, nrow(mat), ncol(mat), byrow=TRUE)
for (x in seq_len(ncol(mat))) {
    fit <- loessFit(log(mat[,x]+cont.cor.scaled[x]) - ab, ab) # explicit subtraction this time.
    ref[,x] <- fit$fitted 
}
ref <- ref-rowMeans(ref)

stopifnot(isTRUE(all.equal(ref, offs)))

###################################################
# End.
