# This tests the profileSites command, to ensure that it's actually giving proper advice.
# library(csaw); library(testthat); source("test-profile.R")

# Running the reference analysis.
CHECKFUN <- function(bam, regions, param, ext, range, use.strand=FALSE, match.strand=FALSE) {
    chromos <- scanBamHeader(bam)[[1]][[1]]

    totally <- totally.reverse <- list()
    for (chr in names(chromos)) {
        cur.chr <- GRanges(chr, IRanges(1, chromos[[chr]]))

        if (!match.strand) { 
            out <- extractReads(bam, cur.chr, param=param, ext=ext)
            totally[[chr]] <- coverage(ranges(out), width=chromos[[chr]]) 
        } else {
            if (param$pe!="both") {
                out <- extractReads(bam, cur.chr, param=param, ext=ext)
                rev.read <- strand(out)=="-"
                totally[[chr]] <- coverage(ranges(out)[!rev.read], width=chromos[[chr]]) 
                totally.reverse[[chr]] <- coverage(ranges(out)[rev.read], width=chromos[[chr]]) 
            } else {
                out <- extractReads(bam, cur.chr, param=param, ext=ext, as.reads=TRUE)
                totally[[chr]] <- coverage(ranges(out$forward), width=chromos[[chr]]) 
                totally.reverse[[chr]] <- coverage(ranges(out$reverse), width=chromos[[chr]]) 
            }
        }
    } 

    relevant.start <- start(regions) - range 
    relevant.end <- start(regions) + range
    if (use.strand) { 
        reverse <- as.logical(strand(regions)=="-")
        relevant.start[reverse] <- end(regions[reverse]) + range # Automatic reversal for reverse-stranded regions.
        relevant.end[reverse] <- end(regions[reverse]) - range
    }

    totes.prof <- matrix(0, nrow=length(regions), ncol=range*2+1)
    for (x in seq_along(regions)) {
        curchr <- as.character(seqnames(regions[x]))
        relevant <- relevant.start[x]:relevant.end[x]
        valid <- relevant > 0L & relevant <= chromos[[curchr]]
                
        # Using reverse coverage if match.strand is TRUE.
        if (match.strand && reverse[x]) { 
            chosen <- totally.reverse            
        } else {
            chosen <- totally
        }

        cur.prof <- as.integer(chosen[[curchr]][relevant[valid]])
        totes.prof[x,valid] <- cur.prof
    }

    colnames(totes.prof) <- (-range):range
    return(totes.prof)
}


nreads <- 5000
chromos <- c(chrA=10000, chrB=5000)
outfname <- tempfile()

set.seed(120000)
test_that("profileSites works correctly with different read extraction settings", {
    for (rparam in list(readParam(), readParam(dedup=TRUE, minq=10))) {
        windows <- generateWindows(chrs=chromos, winsize=10, nwin=20)
        bam <- regenSE(nreads, chromos, outfname)
        
        all.profiles <- profileSites(bam, windows, ext=100, range=100, average=FALSE, param=rparam)
        ref <- CHECKFUN(bam, windows, ext=100, range=100, param=rparam)
        expect_equal(ref, all.profiles)
    }

    # And also on paired-end data.
    rparam <- readParam(pe="both")
    windows <- generateWindows(chrs=chromos, winsize=10, nwin=20)
    bam <- regenPE(nreads, chromos, outfname)
        
    all.profiles <- profileSites(bam, windows, ext=100, range=100, average=FALSE, param=rparam)
    ref <- CHECKFUN(bam, windows, ext=100, range=100, param=rparam)
    expect_equal(ref, all.profiles)
})

set.seed(120001)
test_that("profileSites works correctly with different width and extension", {
    rparam <- readParam()
    for (range in c(20, 100, 500)) {
        for (ext in c(20, 100, 500)) {
            windows <- generateWindows(chrs=chromos, winsize=10, nwin=20)
            bam <- regenSE(nreads, chromos, outfname)

            all.profiles <- profileSites(bam, windows, ext=ext, range=range, average=FALSE, param=rparam)
            ref <- CHECKFUN(bam, windows, ext=ext, range=range, param=rparam)
            expect_equal(ref, all.profiles)
        }
    }
})

set.seed(120002)
test_that("profileSites works correctly with strand-specific extraction", {
    rparam <- readParam()

    for (range in c(20, 100, 500)) {
        windows <- generateWindows(chrs=chromos, winsize=10, nwin=20)
        strand(windows) <- sample(c("+", "-", "*"), length(windows), replace=TRUE)
        bam <- regenSE(nreads, chromos, outfname)
    
        all.profiles <- profileSites(bam, windows, ext=100, range=range, average=FALSE, param=reform(rparam, forward=NULL), strand="use")
        ref <- CHECKFUN(bam, windows, ext=100, range=range, param=rparam, use.strand=TRUE)
        expect_equal(ref, all.profiles)
    
        all.profiles <- profileSites(bam, windows, ext=100, range=range, average=FALSE, param=reform(rparam, forward=NULL), strand="match")
        ref <- CHECKFUN(bam, windows, ext=100, range=range, param=rparam, use.strand=TRUE, match.strand=TRUE)
        expect_equal(ref, all.profiles)
    }

    expect_error(profileSites(bam, windows, ext=100, range=100, average=FALSE, param=rparam, strand="match"), "set forward=NULL")
})

set.seed(120003)
test_that("profileSites correctly averages profiles together", {
    rparam <- readParam()
    windows <- generateWindows(chrs=chromos, winsize=10, nwin=20)
    bam <- regenSE(nreads, chromos, outfname)
    
    all.profiles <- profileSites(bam, windows, ext=100, range=100, average=FALSE, param=rparam)

    none <- profileSites(bam, windows, ext=100, range=100, param=rparam)
    ref.none <- colMeans(all.profiles)
    expect_equal(ref.none, none)
    
    total <- profileSites(bam, windows, ext=100, range=100, param=rparam, normalize="total")
    ref.total <- colMeans(all.profiles/rowSums(all.profiles))
    expect_equal(ref.total, total)

    max <- profileSites(bam, windows, ext=100, range=100, param=rparam, normalize="max")
    ref.max <- colMeans(all.profiles/apply(all.profiles, 1, max))
    expect_equal(ref.max, max)
})

set.seed(120004)
test_that("profileSites correctly handles silly inputs", {
    chromos <- c(chrA=1000, chrB=2000)
    windows <- generateWindows(chrs=chromos, winsize=10, nwin=20)
    rparam <- readParam()
    n <- 1000
    
    # Empty files.
    bam0 <- regenSE(0L, chromos, outfname=tempfile())
    out <- profileSites(bam0, windows, ext=100, range=100, param=rparam, average=FALSE)
    expect_identical(ncol(out), 201L)
    expect_identical(nrow(out), length(windows))
    expect_true(all(out==0))

    # Empty inputs.
    out <- profileSites(bam0, windows[0], ext=100, range=100, param=rparam, average=FALSE)
    expect_identical(ncol(out), 201L)
    expect_identical(nrow(out), 0L)

    out <- profileSites(bam0, windows[0], ext=100, range=100, param=rparam, average=TRUE)
    expect_equivalent(out, rep(NA_real_, 201L))

    # Negative or zero width.
    expect_error(out <- profileSites(bam0, windows, ext=100, range=-100, param=rparam), "positive")
    expect_error(out <- profileSites(bam0, windows, ext=-100, range=100, param=rparam), "positive")
})
