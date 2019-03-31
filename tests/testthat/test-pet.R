# This script tests the getPESizes function.
# library(csaw); library(testthat); source("setup.R"); source("test-pet.R")

tempdir <- tempfile()
dir.create(tempdir)
chromos <- c(chrA=1000, chrB=5000)

CHECKFUN <- function(npairs, nsingles, chromosomes, param) {
    SE <- regenSE(nsingles, chromosomes, file.path(tempdir, "se"))
    PE <- regenPE(npairs, chromosomes, file.path(tempdir, "pe"))
    out <- mergeBam(c(SE, PE), file.path(tempdir, "combined"), indexDestination=TRUE, overwrite=TRUE)
    stuff <- getPESizes(out, param=param)

    chrs.to.use <- names(chromosomes)
    if (length(param$restrict)) {
        chrs.to.use <- intersect(chrs.to.use, param$restrict)
    }

    # Checking sizes.
    param0 <- reform(param, max.frag=Inf)
    ref3 <- lapply(chrs.to.use, FUN=function(chr) { 
        width(extractReads(out, GRanges(chr, IRanges(1, chromosomes[[chr]])), param=param0))
    })
    expect_identical(sort(stuff$sizes), sort(unlist(ref3)))

    # Checking failure statistics.
    everything <- do.call(DataFrame, scanBam(PE)[[1]])
    if (!is.na(param$minq)) everything <- everything[everything$mapq >= param$minq,]    
    if (param$dedup) everything <- everything[bitwAnd(everything$flag, 0x400)==0,]

    everything <- everything[order(everything$qname),]
    first <- everything[bitwAnd(everything$flag, 0x40)!=0,]
    second <- everything[bitwAnd(everything$flag, 0x80)!=0,]

    paired <- intersect(first$qname, second$qname)
    first <- first[match(paired, first$qname),]
    second <- second[match(paired, second$qname),]

    if (length(param$restrict)) {
        # Keep all pairs with at least one read on the restricted chromosomes.
        keep <- first$rname %in% chrs.to.use | second$rname %in% chrs.to.use
        first <- first[keep,]
        second <- second[keep,]
    }

    inter.chr <- first$rname!=second$rname
    expect_identical(stuff$diagnostics[["inter.chr"]], sum(inter.chr)) 
    first <- first[!inter.chr,]
    second <- second[!inter.chr,]

    wrong.strand <- (first$strand==second$strand | 
        (first$strand=="+" & first$pos > end(.make_gr(second, chromos))) | 
        (second$strand=="+" & second$pos > end(.make_gr(first, chromos))))
    expect_identical(stuff$diagnostics[["unoriented"]], sum(wrong.strand))
    first <- first[!wrong.strand,]
    second <- second[!wrong.strand,]

    if (length(param$discard)) {
        ref <- getPESizes(out, param=reform(param, discard=GRanges()))
        expect_identical(stuff$diagnostics[["discarded"]], length(ref$sizes) - length(stuff$sizes))
    } else {
        expect_identical(stuff$diagnostics[["discarded"]], 0L)
    }

    invisible(NULL)
}

.make_gr <- function(df, chromos) {
    maxed <- chromos[df$rname]
    alen <- GenomicAlignments::cigarWidthAlongReferenceSpace(df$cigar)
    GRanges(df$rname, IRanges(df$pos, pmin(maxed, df$pos + alen - 1L)))
}

set.seed(70001)
test_that("getPESizes calculates the diagnostics correctly", {
    chromos <- c(chrA=1000, chrB=3000)
    genome <- GRanges(names(chromos), IRanges(1L, chromos))
    discarder <- makeDiscard(10, 100, chromos)

    for (rparam in list(readParam(pe="both"), 
                        readParam(pe="both", minq=10),
                        readParam(pe="both", dedup=TRUE),
                        readParam(pe="both", discard=discarder),
                        readParam(pe="both", minq=10, dedup=TRUE, discard=discarder),
                        readParam(pe="both", restrict="chrB")
                        )) {
        for (npairs in c(1000L, 5000L)) {
            CHECKFUN(npairs=npairs, nsingles=20L, chromosomes=chromos, param=rparam)
        }
    }
            
    # Checking that it still does something sensible with no singles.
    set.seed(12)
    CHECKFUN(npairs=1000L, nsingles=0L, chromosomes=chromos, param=readParam(pe="both", minq=10))
})

set.seed(70002)
test_that("getPESizes responds to silly inputs correctly", {
    obam <- regenPE(0L, chromos, outfname=tempfile())
    out <- getPESizes(obam, readParam(pe="both"))
    expect_identical(out$sizes, integer(0))
    expect_true(all(out$diagnostics==0L))

    expect_error(out <- getPESizes(obam, readParam()), "paired-end inputs required")
})
