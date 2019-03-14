# Testing the correlateReads function.
# library(csaw); library(testthat); source("test-ccf.R")

library(Rsamtools)
CHECKFUN <- function(bam.files, max.dist, param, cross=FALSE) {
	chromos <- scanBamHeader(bam.files[1])[[1]][[1]]
	out <- integer(max.dist + 1)
	total <- 0L

	for (chr in names(chromos)) {
		clen <- chromos[[chr]]
        chromosome <- GRanges(chr, IRanges(1, clen))
        
		all.f <- all.r <- vector("list", length(bam.files))
		for (b in seq_along(bam.files)) {
            if (param$pe!="both") {
                reads <- extractReads(bam.files[b], chromosome, param=param)
                is.forward <- as.logical(strand(reads)=="+")
                all.f[[b]] <- start(reads)[is.forward]
                all.r[[b]] <- pmin(end(reads)[!is.forward], clen) + 1L
            } else {
                reads <- extractReads(bam.files[b], chromosome, param=param, as.reads=TRUE)
                all.f[[b]] <- start(reads$forward)
                all.r[[b]] <- pmin(end(reads$reverse), clen) + 1L
            }
		}

		f <- tabulate(unlist(all.f), clen+1L)
		r <- tabulate(unlist(all.r), clen+1L)

		# Autocorrelations, if not cross-correlations, so we just fuse them together.
		if (!cross) {
			comb <- f+r
			r <- f <- comb
		}
      
		fnreads <- sum(lengths(all.f))
		rnreads <- sum(lengths(all.r))
		nreads <- fnreads + rnreads
        if (nreads==0L) {
            next
        }

        # Ignoring chromosomes without forward or reverse reads, if we want 'cross'.
        if (cross && (fnreads==0L || rnreads==0L)) { 
            next 
        }

		out <- out + nreads * vapply(0:max.dist, FUN=function(i){ 
			if (i>=length(f)-1L || i>=length(r)-1L) { return(0) }
			fr <- f[1:(length(f)-i)]
			rr <- r[(i+1):length(r)]
			if (sd(fr)==0 || sd(rr)==0) { return(0) }
			cor(fr, rr)
		}, FUN.VALUE=0)
		total <- total + nreads
	}

    # Avoid problems when total=0.
	out/pmax(1, total)
}

set.seed(800000)
test_that("correlateReads works correctly with single BAM files", {
    n <- 1000
    for (rparam in list(readParam(), readParam(dedup=TRUE, minq=10))) {
        for (nreads in c(100, 1000, 10000)) {
            for (chromos in list(c(chrA=10000), c(chrA=1000, chrB=2000))) { 
                bam <- regenSE(nreads, chromos, outfname=tempfile())
                out <- CHECKFUN(bam, max.dist=n, cross=TRUE, param=rparam)
                out2 <- correlateReads(bam, max.dist=n, cross=TRUE, param=rparam)
                expect_equal(out, out2)
            }
        }
	}
})

set.seed(800001)
test_that("correlateReads works correctly with autocorrelations", {
    n <- 1000
    rparam <- readParam()
    for (nreads in c(100, 1000, 10000)) {
        for (chromos in list(c(chrA=10000), c(chrA=1000, chrB=2000))) {
            bam <- regenSE(nreads, chromos, outfname=tempfile())
            out <- CHECKFUN(bam, max.dist=n, cross=FALSE, param=rparam)
            out2 <- correlateReads(bam, max.dist=n, cross=FALSE, param=rparam)
            expect_equal(out, out2)
        }
    }
})

set.seed(800002)
test_that("correlateReads works correctly with multiple BAM files", {
    chromos <- c(chrA=1000, chrB=2000)
    nreads <- 10000
    n <- 1000
    rparam <- readParam()

    bam.files <- c(regenSE(nreads, chromos, outfname=tempfile()),
                   regenSE(nreads, chromos, outfname=tempfile()))

    # Cross-correlations.        
    out <- CHECKFUN(bam.files, max.dist=n, cross=TRUE, param=rparam)
    out2 <- correlateReads(bam.files, max.dist=n, cross=TRUE, param=rparam)
    expect_equal(out, out2)

    # Auto-correlations.        
    out <- CHECKFUN(bam.files, max.dist=n, cross=FALSE, param=rparam)
    out2 <- correlateReads(bam.files, max.dist=n, cross=FALSE, param=rparam)
    expect_equal(out, out2)
})

set.seed(8000021)
test_that("correlateReads works correctly with paired-end BAM files", {
    chromos <- c(chrA=1000, chrB=2000)
    nreads <- 10000
    n <- 1000
    rparam <- readParam(pe="both")

    bam.files <- c(regenPE(nreads, chromos, outfname=tempfile()),
                   regenPE(nreads, chromos, outfname=tempfile()))

    # Cross-correlations.        
    out <- CHECKFUN(bam.files, max.dist=n, cross=TRUE, param=rparam)
    out2 <- correlateReads(bam.files, max.dist=n, cross=TRUE, param=rparam)
    expect_equal(out, out2)

    # Auto-correlations.        
    out <- CHECKFUN(bam.files, max.dist=n, cross=FALSE, param=rparam)
    out2 <- correlateReads(bam.files, max.dist=n, cross=FALSE, param=rparam)
    expect_equal(out, out2)
})

set.seed(800003)
test_that("correlateReads works correctly when 'max.dist' is longer than the chromosome", {
    rparam <- readParam()
    n <- 1000
    nreads <- 10000

    # Testing shorter than 'n'.
    chromos <- c(chrA=500, chrB=100, chrC=800) 
    bam <- regenSE(nreads, chromos, outfname=tempfile())
    out <- correlateReads(bam, max.dist=n, cross=TRUE, param=rparam)
    out2 <- CHECKFUN(bam, max.dist=n, cross=TRUE, param=rparam)
    expect_equal(out, out2)
    expect_true(all(out[(max(chromos)+1):length(out)]==0))

    # Testing on a tiny chromosome.
    chromos <- c(chrA=1)
    bam2 <- regenSE(1000, chromos, outfname=tempfile())
    out <- correlateReads(bam2, max.dist=n, cross=TRUE, param=rparam)
    out2 <- CHECKFUN(bam2, max.dist=n, cross=TRUE, param=rparam)
    expect_equal(out, out2)
    expect_equal(out[1], -1) # all forwards at 1, all reverse at 1-past-end, i.e., 2.
    expect_true(all(out[-1]==0))
})

set.seed(800003)
test_that("correlateReads works correctly with silly inputs", {
    chromos <- c(chrA=1000, chrB=2000)
    rparam <- readParam()
    n <- 1000

    # Empty files.
    bam0 <- regenSE(0L, chromos, outfname=tempfile())
    out <- correlateReads(bam0, max.dist=n, cross=TRUE, param=rparam)
    expect_equal(out, numeric(n+1L))

    # Reads on only one strand, which should trigger cross=TRUE checks.
    bam1 <- regenSE(1000, chromos, outfname=tempfile())

    out <- correlateReads(bam1, max.dist=n, cross=TRUE, param=readParam(forward=TRUE))
    expect_equal(out, numeric(n+1L))
    out <- correlateReads(bam1, max.dist=n, cross=FALSE, param=readParam(forward=TRUE))
    out2 <- CHECKFUN(bam1, max.dist=n, cross=FALSE, param=readParam(forward=TRUE))
    expect_equal(out, out2)

    out <- correlateReads(bam1, max.dist=n, cross=TRUE, param=readParam(forward=FALSE))
    expect_equal(out, numeric(n+1L))
    out <- correlateReads(bam1, max.dist=n, cross=FALSE, param=readParam(forward=FALSE))
    out2 <- CHECKFUN(bam1, max.dist=n, cross=FALSE, param=readParam(forward=FALSE))
    expect_equal(out, out2)

    # Negative or zero width.
    expect_error(correlateReads(bam1, max.dist=0), "positive")
})
