# This checks the checkBimodality function.
# library(csaw); library(testthat); source("setup.R"); source("test-bimodal.R")

CHECKFUN <- function(bam.files, regions, width, param, prior.count) {
    width <- rep(width, length.out=length(bam.files))
	output <- rep(NA, length(regions))
	chr.sizes <- Rsamtools::scanBamHeader(bam.files[1])[[1]][[1]]
    if (length(param$restrict) > 0) {
        chr.sizes <- chr.sizes[names(chr.sizes) %in% param$restrict]
    }

	for (chr in names(chr.sizes)) { 
		full.chr <- GRanges(chr, IRanges(1, chr.sizes[[chr]]))
		left.forward <- right.forward <- left.reverse <- right.reverse <- 0

		for (b in seq_along(bam.files)) { 
            if (param$pe!="both") { 
                all.reads <- extractReads(bam.files[b], full.chr, param=param)
                is.forward <- strand(all.reads)=="+"
                forward.reads <- all.reads[is.forward]
                reverse.reads <- all.reads[!is.forward]
            } else {
                all.reads <- extractReads(bam.files[b], full.chr, param=param, as.reads=TRUE)
                forward.reads <- all.reads$forward
                reverse.reads <- all.reads$reverse
            }

            # Computing coverage, somewhat counterintuitively; 'left.forward', for example, is extended to 
            # the right, as it must find all overlaps to the left of the region of interest.
			start.F <- start(forward.reads)
			end.F <- end(forward.reads)
			left.forward <- left.forward + coverage(IRanges(start.F, start.F+width[b]-1), width=end(full.chr))
			right.forward <- right.forward + coverage(IRanges(end.F-width[b]+1, end.F), width=end(full.chr))

			end.R <- end(reverse.reads)
			start.R <- start(reverse.reads)
			left.reverse <- left.reverse + coverage(IRanges(start.R, start.R+width[b]-1), width=end(full.chr))
			right.reverse <- right.reverse + coverage(IRanges(end.R-width[b]+1, end.R), width=end(full.chr))
		}

		for (r in BiocGenerics::which(seqnames(regions)==chr)) { # using BiocGenerics:: as testthat is unhappy during CHECK.
			cur.reg <- regions[r]

			# Computing the bimodality statistic.
			lf.set <- as.integer(left.forward[start(cur.reg):end(cur.reg)])
			rf.set <- as.integer(right.forward[start(cur.reg):end(cur.reg)])
			lr.set <- as.integer(left.reverse[start(cur.reg):end(cur.reg)])
			rr.set <- as.integer(right.reverse[start(cur.reg):end(cur.reg)])

  			bistat <- max(pmin((lf.set+prior.count)/(lr.set+prior.count), (rr.set+prior.count)/(rf.set+prior.count)))
			output[r] <- bistat
		}
	}

    return(output)
}


chromos<-c(chrA=5000, chrB=5000, chrC=8000)
tmpdir <- tempfile()
dir.create(tmpdir)

set.seed(39000)
test_that("checkBimodality works correctly", {
    for (param in list(readParam(), 
                       readParam(minq=10),
                       readParam(discard=makeDiscard(10, 200, chromos)),
                       readParam(restrict=c("chrC", "chrB")))) {
        bam.files <-c(regenSE(1000, chromos, file.path(tmpdir, "A")), regenSE(1000, chromos, file.path(tmpdir, "B")))
        my.ranges <- generateWindows(chromos, 10, 10)
        out <- checkBimodality(bam.files, my.ranges, param=param, width=100, prior.count=2)
        ref <- CHECKFUN(bam.files, my.ranges, param=param, width=100, prior.count=2)
        expect_equal(out, ref)
    }

    # Responds to different width and prior counts.
    for (width in c(20, 100, 500)) {
        for (prior in c(1,5)) { 
            bam.files <-c(regenSE(1000, chromos, file.path(tmpdir, "A")), regenSE(1000, chromos, file.path(tmpdir, "B")))
            my.ranges <- generateWindows(chromos, 10, 10)
            out <- checkBimodality(bam.files, my.ranges, param=param, width=width, prior.count=prior)
            ref <- CHECKFUN(bam.files, my.ranges, param=param, width=width, prior.count=prior)
            expect_equal(out, ref)
        }
    }
})

set.seed(39001)
test_that("checkBimodality responds to variable width and window sizes", {
    bam.files <-c(regenSE(1000, chromos, file.path(tmpdir, "A")), regenSE(1000, chromos, file.path(tmpdir, "B")))
    param <- readParam()

    for (lower in c(10, 20, 50)) { 
        my.ranges <- generateWindows(chromos, 20, round(runif(20, lower, lower*10)))
        out <- checkBimodality(bam.files, my.ranges, param=param, width=100, prior.count=2)
        ref <- CHECKFUN(bam.files, my.ranges, param=param, width=100, prior.count=2)
        expect_equal(out, ref)
    }

    out <- checkBimodality(bam.files, my.ranges, param=param, width=list(c(100, 200), NA), prior.count=2)
    ref <- CHECKFUN(bam.files, my.ranges, param=param, width=c(100, 200), prior.count=2)
    expect_equal(out, ref)    
})

set.seed(39002)
test_that("checkBimodality handles paired-end reads correctly", {
    bam.files <-c(regenPE(10000, chromos, file.path(tmpdir, "A")), regenPE(10000, chromos, file.path(tmpdir, "B")))
    param <- readParam(pe="both")
    my.ranges <- generateWindows(chromos, 10, 10)

    out <- checkBimodality(bam.files, my.ranges, param=param, width=100, prior.count=2)
    ref <- CHECKFUN(bam.files, my.ranges, param=param, width=100, prior.count=2)
    expect_equal(out, ref)    
})


set.seed(39003)
test_that("checkBimodality fails correctly on empty inputs", {
    bam.files <-c(regenSE(0, chromos, file.path(tmpdir, "A")), regenSE(0, chromos, file.path(tmpdir, "B")))
    my.ranges <- generateWindows(chromos, 20, 20)

    out <- checkBimodality(bam.files, my.ranges[0])
    expect_identical(out, numeric(0))

    out <- checkBimodality(bam.files, my.ranges[1:10]) # prior counts come into play.
    expect_identical(out, rep(1, 10L))
})
