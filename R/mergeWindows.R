#' @export
#' @importFrom BiocGenerics strand strand<- start end
#' @importFrom S4Vectors runValue
#' @importFrom GenomeInfoDb seqnames Seqinfo seqlevels seqlengths
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
mergeWindows <- function(ranges, tol, signs=NULL, max.width=NULL, ignore.strand=TRUE)
# This function merges the adjacent windows if they lie within 'tol' of each other,
# Any abundance filtering should be done beforehand. Negative values of tol refer
# to a minimum overlap. A value of zero means that the windows must be adjacent
# (i.e. non-overlapping and contiguous).
# 
# written by Aaron Lun
# created 30 July 2013
{
    ranges <- .toGRanges(ranges)
    strs <- strand(ranges)
    if (!ignore.strand && length(runValue(strs))!=1) {
        # Strand-specific clustering.
        forward <- as.logical(strs=="+")
        reverse <- as.logical(strs=="-")
        neither <- as.logical(strs=="*")

        ids <- integer(length(ranges))
        regs <- vector("list", 3)
        collected <- 0L

        if (any(forward)) { 
            out <- Recall(ranges=ranges[forward], tol=tol, signs=signs[forward], max.width=max.width, ignore.strand=TRUE) 
            ids[forward] <- out$id
            strand(out$region) <- "+"
            regs[[1]] <- out$region
            collected <- collected + length(out$region)
        }

        if (any(reverse)) { 
            out <- Recall(ranges=ranges[reverse], tol=tol, signs=signs[reverse], max.width=max.width, ignore.strand=TRUE) 
            ids[reverse] <- out$id + collected
            strand(out$region) <- "-"
            regs[[2]] <- out$region
            collected <- collected + length(out$region)
        }

        if (any(neither)) { 
            out <- Recall(ranges=ranges[neither], tol=tol, signs=signs[neither], max.width=max.width, ignore.strand=TRUE) 
            ids[neither] <- out$id + collected
            regs[[3]] <- out$region
        }

        return(list(ids=ids, regions=do.call(c, regs)))
    }

    tol <- as.integer(tol)
    max.width <- as.integer(max.width)

    chrs <- as.integer(seqnames(ranges))
    starts <- start(ranges)
    ends <- end(ranges)
    o <- order(chrs, starts, ends)

    if (is.null(signs)) { 
        signs <- logical(length(ranges)) 
    } else {
        if (length(signs)!=length(ranges)) {
            stop("lengths of 'signs' and 'ids' must be the same")
        }
        signs <- signs[o]
    }

    # Running the merge. Indices correspond with positions in 'clustered'.
    out <- .Call(cxx_merge_windows, chrs[o], starts[o], ends[o], signs, tol, max.width)
    out[[1]][o] <- out[[1]]
    clustered <- GRanges(levels(seqnames(ranges))[out[[2]]], IRanges(out[[3]], out[[4]]),
        seqinfo=Seqinfo(seqlevels(ranges), seqlengths(ranges)))

    list(ids=out[[1]], regions=clustered)
}
