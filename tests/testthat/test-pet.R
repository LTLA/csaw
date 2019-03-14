# This script tests the getPESizes function.
# library(csaw); library(testthat); source("test-pet.R")

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

    # Checking singletons. 
    ref1 <- lapply(chrs.to.use, FUN=function(chr) { 
        extractReads(SE, GRanges(chr, IRanges(1, chromosomes[[chr]])), param=reform(param, pe="none")) 
    })
    se.nreads <- sum(sapply(ref1, length))
    expect_identical(stuff$diagnostics[["single"]], se.nreads)

    # Checking mapped reads.
    ref2 <- lapply(chrs.to.use, FUN=function(chr) { 
        extractReads(PE, GRanges(chr, IRanges(1, chromosomes[[chr]])), param=reform(param, pe="none")) 
    })
    pe.nreads <- sum(sapply(ref2, length))
    expect_identical(stuff$diagnostics[["mapped.reads"]], se.nreads + pe.nreads)

    # Checking total reads.
    expect_identical(stuff$diagnostics[["total.reads"]], npairs*2L+nsingles)

    # Checking sizes.
    ref3 <- lapply(chrs.to.use, FUN=function(chr) { 
        width(extractReads(out, GRanges(chr, IRanges(1, chromosomes[[chr]])), param=reform(param, max.frag=1e8)))
    })
    expect_identical(sort(stuff$sizes), sort(unlist(ref3)))

    # Checking interchromosomal reads and all that. 
    everything <- do.call(DataFrame, scanBam(PE)[[1]])
    if (!is.na(param$minq)) everything <- everything[everything$mapq >= param$minq,]    
    if (param$dedup) everything <- everything[bitwAnd(everything$flag, 0x400)==0,]
    if (length(param$restrict)) everything <- everything[everything$rname %in% param$restrict,]
    if (length(param$discard)) {
        gr <- .make_gr(everything, chromos)
        everything <- everything[!overlapsAny(gr, param$discard, type="within"),]
    }

    everything <- everything[order(everything$qname),]
    first <- everything[bitwAnd(everything$flag, 0x40)!=0,]
    second <- everything[bitwAnd(everything$flag, 0x80)!=0,]
    paired <- intersect(first$qname, second$qname)
    expect_identical(stuff$diagnostics[["mate.unmapped"]], nrow(first) + nrow(second) - length(paired) * 2L) 

    first <- first[match(paired, first$qname),]
    second <- second[match(paired, second$qname),]
    expect_identical(stuff$diagnostics[["inter.chr"]], sum(first$rname!=second$rname)) 

    wrong.strand <- first$rname==second$rname & 
            (first$strand==second$strand | 
            (first$strand=="+" & first$pos > end(.make_gr(second, chromos))) | 
            (second$strand=="+" & second$pos > end(.make_gr(first, chromos))))
    expect_identical(stuff$diagnostics[["unoriented"]], sum(wrong.strand))

    expect_identical(stuff$diagnostics[["mapped.reads"]], 
            (length(stuff$sizes) + sum(stuff$diagnostics[c("inter.chr", "unoriented")]))*2L + sum(stuff$diagnostics[c("mate.unmapped", "single")]))
    return(NULL)
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
    set.seed(12); CHECKFUN(npairs=1000L, nsingles=0L, chromosomes=chromos, param=readParam(pe="both", minq=10))
})

set.seed(70002)
test_that("getPESizes responds to silly inputs correctly", {
    obam <- regenPE(0L, chromos, outfname=tempfile())
    out <- getPESizes(obam, readParam(pe="both"))
    expect_identical(out$sizes, integer(0))
    expect_true(all(out$diagnostics==0L))

    expect_error(out <- getPESizes(obam, readParam()), "paired-end inputs required")
})

