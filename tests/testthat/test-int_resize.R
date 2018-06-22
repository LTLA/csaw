# This tests the internal resizing functions.
# library(testthat); library(csaw); source("test-int_resize.R")

set.seed(230000)
test_that(".coerceFragments works correctly", {
    starts <- as.integer(round(runif(100, 1, 1000)))
    widths <- as.integer(round(runif(100, 1, 100)))
    ends <- starts + widths

    # Chrlen and positivity are enforced.
    chrlen <- 1000L
    out <- csaw:::.coerceFragments(starts, ends, final=NA, chrlen=chrlen)
    expect_identical(out$end, pmin(ends, chrlen))
    expect_identical(out$start, starts)

    out <- csaw:::.coerceFragments(-starts, ends, final=NA, chrlen=1e5)
    expect_true(all(out$start==1L))
    expect_identical(out$end, ends)

    # Final width is enforced.    
    final <- 50
    out <- csaw:::.coerceFragments(starts, ends, final=final, chrlen=1e6)
    ref <- resize(IRanges(starts, ends), width=final, fix="center")
    expect_identical(out$start, pmax(1L, start(ref)))
    expect_identical(out$end, end(ref))

    chrlen <- 500L
    final <- 21 # odd width
    out <- csaw:::.coerceFragments(starts, ends, final=final, chrlen=chrlen)
    ref <- resize(IRanges(starts, ends), width=final, fix="center")
    expect_identical(out$start, pmin(pmax(1L, start(ref)), chrlen))
    expect_identical(out$end, pmin(end(ref), chrlen))

    # Works correcty on empty inputs.
    out <- csaw:::.coerceFragments(integer(0), integer(0), final=NA, chrlen=1e6)
    expect_identical(out$start, integer(0))
    expect_identical(out$end, integer(0))

    out <- csaw:::.coerceFragments(integer(0), integer(0), final=50L, chrlen=1e6)
    expect_identical(out$start, integer(0))
    expect_identical(out$end, integer(0))
})

test_that(".collateExt works correctly", {
    nbam <- 4

    # Simple number.
    ext <- 100
    out <- csaw:::.collateExt(nbam, ext)
    expect_identical(out$ext, rep(as.integer(ext), nbam))
    expect_identical(out$final, NA_integer_)

    variety <- c(100, 40, 20, 80)
    expect_error(csaw:::.collateExt(nbam, variety), "must be an integer scalar")
    expect_error(csaw:::.collateExt(nbam, -ext), "positive")

    # A list of length 2.
    final <- 50L
    out <- csaw:::.collateExt(nbam, list(variety, final))
    expect_identical(out$ext, as.integer(variety))
    expect_identical(out$final, final)

    expect_error(csaw:::.collateExt(nbam, list(ext)), "must be a list of length 2")
    expect_error(csaw:::.collateExt(nbam, list(ext, final)), "not consistent with the number of libraries")
    expect_error(csaw:::.collateExt(nbam, list(variety, -final)), "positive")
})

set.seed(230001)
test_that(".extendSE works correctly", {
    # Simulating reads.
    fstarts <- as.integer(round(runif(100, 1, 1000)))
    fwidths <- as.integer(round(runif(100, 1, 100)))
    rstarts <- as.integer(round(runif(50, 1, 1000)))
    rwidths <- as.integer(round(runif(50, 1, 100)))
    reads <- list(forward=list(pos=fstarts, qwidth=fwidths), reverse=list(pos=rstarts, qwidth=rwidths))
    
    # Running it through the extender with no extension.
    chrlen <- 800L
    out <- csaw:::.extendSE(reads, ext=NA, final=NA, chrlen=chrlen)

    new.fstarts <- pmin(fstarts, chrlen)
    new.fends <- pmin(fstarts + fwidths - 1L, chrlen)
    new.rstarts <- pmin(rstarts, chrlen)
    new.rends <- pmin(rstarts + rwidths - 1L, chrlen)

    expect_identical(out$start, c(new.fstarts, new.rstarts))
    expect_identical(out$end, c(new.fends, new.rends))

    # Running it through the extender.
    ext <- 80L
    chrlen <- 800L
    out <- csaw:::.extendSE(reads, ext=ext, final=NA, chrlen=chrlen)

    new.fstarts <- pmin(fstarts, chrlen)
    new.fends <- pmin(fstarts + ext - 1L, chrlen)
    new.rends <- rstarts + rwidths 
    new.rstarts <- pmin(pmax(1L, new.rends - ext), chrlen)
    new.rends <- pmin(new.rends - 1L, chrlen)

    expect_identical(out$start, c(new.fstarts, new.rstarts))
    expect_identical(out$end, c(new.fends, new.rends))
    
    # Testing it with a final setting.
    final <- 50L
    ext <- 80L
    chrlen <- 800L
    out <- csaw:::.extendSE(reads, ext=ext, final=final, chrlen=chrlen)

    new.fwidth <- rep(ext, length(fstarts))
    new.rends <- rstarts + rwidths 
    new.rstarts <- new.rends - ext
    new.rwidth <- rep(ext, length(rstarts))
    new.reads <- list(forward=list(pos=fstarts, qwidth=new.fwidth),
                      reverse=list(pos=new.rstarts, qwidth=new.rwidth))

    ref <- csaw:::.extendSE(new.reads, ext=NA, final=final, chrlen=chrlen)
    expect_identical(out, ref)
})

