#' @export
clusterWindowsList <- function(ranges.list, tab.list, equiweight=TRUE, ...) 
# This does the same as clusterWindows, but for results from many different analyses
# (ostensibly with different window sizes).
#
# written by Aaron Lun
# created 8 January 2016
{
    nset <- length(ranges.list)
    set.it.vec <- seq_len(nset)
    if (nset!=length(tab.list)) { stop("ranges.list must have same length as tab.list") }
    
    for (x in set.it.vec) {
        ranges.list[[x]] <- .toGRanges(ranges.list[[x]])
        currows <- length(ranges.list[[x]])
        ntab <- nrow(tab.list[[x]])
        if (currows!=ntab) { stop("corresponding entries of data and tab.lists must have same number of entries") }
    }
    
    # Merging everyone together.
    all.data <- do.call(c, ranges.list)
    all.result <- do.call(rbind, tab.list)
    
    # Computing weights based on number of windows; each analysis contributes same effective number of tests.
    # Not quite the same as consolidateWindows()'s weighting, but that can't be helped as the adjusted p-values
    # are computed before clustering.
    if (equiweight) { 
        weights <- rep(1/lengths(ranges.list), lengths(ranges.list))
    } else {
        weights <- NULL
    }   

    out <- clusterWindows(all.data, all.result, weights=weights, ...)
    c(list(ranges=all.data), out)
}

#' @export
consolidateClusters <- function(data.list, result.list, equiweight=TRUE, ...) {
    .Deprecated("clusterWindowsList")
    out <- clusterWindowsList(ranges.list=data.list, tab.list=result.list, equiweight=equiweight, ...)
    groupings <- rep(seq_along(data.list), lengths(data.list)) 
    out$ids <- split(out$ids, groupings)
    names(out$ids) <- names(data.list)
    out
}

