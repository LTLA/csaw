#' @importFrom GenomicAlignments readGAlignments
#' @importFrom BiocGenerics start width 
#' @importFrom IRanges overlapsAny
#' @importFrom GenomeInfoDb seqnames
.extractSE <- function(bam.file, where, param) 
# Extracts single-end read data from a BAM file with removal of unmapped,
# duplicate and poorly mapped/non-unique reads. We also discard reads in the
# specified discard regions. 
#
# written by Aaron Lun
# created 8 December 2013
{
    sbp <- .generate_sbp(where, param)
    reads <- readGAlignments(bam.file, param=sbp)

    is.forward <- strand(reads)=="+"
    forwards <- reads[is.forward]
    reverses <- reads[!is.forward]

    if (any(seqnames(param$discard) == as.character(seqnames(where)))) {
        forwards <- forwards[!overlapsAny(forwards, param$discard, type="within")]
        reverses <- reverses[!overlapsAny(reverses, param$discard, type="within")]
    }

    list(
        forward=list(pos=start(forwards), qwidth=width(forwards)),
        reverse=list(pos=start(reverses), qwidth=width(reverses))
    )
}

#' @importFrom Rsamtools ScanBamParam scanBamFlag
.generate_sbp <- function(where, param) {
    flags <- list(isUnmappedQuery=FALSE,
        isSecondaryAlignment=FALSE,
        isSupplementaryAlignment=FALSE)

    if (param$pe=="first" || param$pe=="second") {
        flags$isFirstMateRead <- param$pe=="first"
        flags$isSecondMateRead <- param$pe=="second"
    } else if (param$pe=="second") {
        flags$isPaired <- TRUE
        flags$hasUnmappedMate <- FALSE
    }

    if (length(param$forward)==0L) { 
        stop("read strand extraction must be specified") 
    } else if (!is.na(param$forward)) {
        flags$isMinusStrand <- !param$forward
    }

    if (param$dedup) {
        flags$isDuplicate <- FALSE
    }

    ScanBamParam(flag=do.call(scanBamFlag, flags), mapqFilter=param$minq, which=where)
}

#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges overlapsAny IRanges
#' @importFrom BiocGenerics start end strand
#' @importFrom GenomeInfoDb seqnames
#' @importFrom GenomicAlignments readGAlignmentPairs first last
.extractPE <- function(bam.file, where, param, with.reads=FALSE, diagnostics=FALSE)
# A function to extract PE data for a particular chromosome. 
# 
# written by Aaron Lun
# created 8 December 2013
{
    sbp <- .generate_sbp(where, param)
    reads <- readGAlignmentPairs(bam.file, param=sbp)
    first.read <- first(reads)
    second.read <- last(reads)

    # Checking that reads in a pair are on the same chromosome and different strands.
    cur.chr <- as.character(seqnames(where))
    same.chr1 <- seqnames(first.read) == cur.chr
    same.chr2 <- seqnames(second.read) == cur.chr
    same.chr <- same.chr1 & same.chr2

    if (diagnostics) {
        inter.set <- table(seqnames(first.read)[!same.chr1]) + table(seqnames(second.read)[!same.chr2])
    }

    oriented <- strand(first.read) != strand(second.read)
    keep <- oriented & same.chr
    first.read <- first.read[keep]
    second.read <- second.read[keep]

    # Reorganizing the objects in terms of forward and reverse reads.
    first.forward <- strand(first.read)=="+"
    second.forward <- strand(second.read)=="+"
    forward.read <- c(first.read[first.forward], second.read[second.forward])
    reverse.read <- c(second.read[first.forward], first.read[second.forward])

    # Insert size from forward start to reverse end.
    frag.starts <- start(forward.read)
    frag.sizes <- end(reverse.read) - frag.starts + 1L
    inward <- frag.sizes > 0L

    keep <- inward & frag.sizes <= param$max.frag 
    frag.starts <- frag.starts[keep]
    frag.sizes <- frag.sizes[keep]

    # Discarding *fragments* in blacklist regions.
    if (any(seqnames(param$discard) == as.character(seqnames(where))) && length(frag.sizes)) {
        cur.ranges <- GRanges(cur.chr, IRanges(frag.starts, width=frag.sizes))
        blacklisted <- overlapsAny(cur.ranges, param$discard, type="within")
        frag.starts <- frag.starts[!blacklisted]
        frag.sizes <- frag.sizes[!blacklisted]
        keep[keep] <- !blacklisted
    } else {
        blacklisted <- NULL
    }

    output <- list(pos=frag.starts, size=frag.sizes)

    if (with.reads) {
        forward.read <- forward.read[keep]
        reverse.read <- reverse.read[keep]
        output$forward <- list(pos=start(forward.read), qwidth=width(forward.read))
        output$reverse <- list(pos=start(reverse.read), qwidth=width(reverse.read))
    }

    if (diagnostics) {
        output$diagnostics <- list(
            inter.chr=inter.set,
            unoriented=sum(same.chr & !oriented) + sum(!inward), 
            discarded=sum(blacklisted)
        )
    }

    output
}
