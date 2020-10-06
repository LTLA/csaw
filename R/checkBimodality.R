#' @export
#' @importFrom GenomeInfoDb seqnames
#' @importFrom BiocGenerics start end order
#' @importFrom S4Vectors split
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
#' @importFrom BiocParallel bpmapply bpisup bpstart bpstop SerialParam
checkBimodality <- function(bam.files, regions, width=100, param=readParam(), prior.count=2, invert=FALSE, BPPARAM=SerialParam()) 
# Gives the maximum strand bimodality score for any base pair in each region. 
# The idea is to try to distinguish between genuine TF binding sites and read stacks or other artifacts.
#
# written by Aaron Lun
# created 1 May 2015
{
    bam.files <- .make_BamFiles(bam.files)
    nbam <- length(bam.files)
    ext.data <- .collateExt(nbam, width) 
    invert <- as.logical(invert)

    totals <- integer(nbam)
    nx <- length(regions)
    out.scores <- rep(NA_real_, nx)
    indices <- split(seq_len(nx), seqnames(regions))
    prior.count <- as.double(prior.count)

    if (!bpisup(BPPARAM)) {
        bpstart(BPPARAM)
        on.exit(bpstop(BPPARAM))
    }

    extracted.chrs <- .activeChrs(bam.files, param$restrict)
    for (chr in names(extracted.chrs)) {
        chosen <- indices[[chr]]
        if (length(chosen)==0L) { next } 
        outlen <- extracted.chrs[[chr]]
        where <- GRanges(chr, IRanges(1, outlen))

        # Pulling out read data.
        collected <- bpmapply(FUN=.check_bimodality, bam.file=bam.files, init.ext=ext.data$ext, 
             MoreArgs=list(where=where, param=param, final.ext=ext.data$final, outlen=outlen),
             BPPARAM=BPPARAM, SIMPLIFY=FALSE)

        # Checking region order.
        rstarts <- start(regions)[chosen]
        rends <- end(regions)[chosen]
        ro <- order(rstarts)

        # Computing bimodality scores.
        out <- .Call(cxx_check_bimodality, collected, rstarts[ro], rends[ro], prior.count, invert)
        out.scores[chosen][ro] <- out 
    }

    return(out.scores)
}

.check_bimodality <- function(bam.file, where, param, init.ext, final.ext, outlen)
{
    if (param$pe=="both") {
        reads <- .extractPE(bam.file, where=where, param=param, with.reads=TRUE)
    } else {
        reads <- .extractSE(bam.file, where=where, param=param)
    }

    # Computing what would happen if we extended one way and the other.
    Fstandard <- .extendSEdir(reads$forward, ext=init.ext, final=final.ext, chrlen=outlen, forward=TRUE)
    Fflipped <- .extendSEdir(reads$forward, ext=init.ext, final=final.ext, chrlen=outlen, forward=FALSE)
    Rstandard <- .extendSEdir(reads$reverse, ext=init.ext, final=final.ext, chrlen=outlen, forward=FALSE)
    Rflipped <- .extendSEdir(reads$reverse, ext=init.ext, final=final.ext, chrlen=outlen, forward=TRUE)
    
    # Standard extension for originally forward reads goes to (2), flipped extension go to (1) as they'll be at an earlier position.
    # Opposite is true for reverse reads; standard extension goes to (1), and flipped extension goes to (2).
    earlier <- mapply(c, Fflipped, Rstandard, SIMPLIFY=FALSE)
    later <- mapply(c, Fstandard, Rflipped, SIMPLIFY=FALSE)
    start1 <- earlier$start
    end1 <- earlier$end
    start2 <- later$start
    end2 <- later$end
    original.strand <- rep(c(1L, 0L), c(length(reads$forward$pos), length(reads$reverse$pos)))

    # Sorting, as required.
    o <- order(start1)
    if (any(start1 > start2)) { 
        stop("extension of flipped alignment should not be before the original alignment")
    }
    list(start1[o], end1[o], start2[o], end2[o], original.strand[o])
}

