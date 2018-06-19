#' @export
consolidateClusters <- function(data.list, result.list, equiweight=TRUE, ...) 
# This does the same as clusterWindows, but for results from many different analyses
# (ostensibly with different window sizes).
#
# written by Aaron Lun
# created 8 January 2016
{
    nset <- length(data.list)
    set.it.vec <- seq_len(nset)
    if (nset!=length(result.list)) { stop("data list must have same length as result list") }
    
    for (x in set.it.vec) {
        data.list[[x]] <- .toGRanges(data.list[[x]])
        currows <- length(data.list[[x]])
        ntab <- nrow(result.list[[x]])
        if (currows!=ntab) { stop("corresponding entries of data and result lists must have same number of entries") }
    }
    
    # Merging everyone together.
    all.data <- do.call(c, data.list)
    all.result <- do.call(rbind, result.list)
    groupings <- rep(seq_along(data.list), lengths(data.list)) 
    
    # Computing weights based on number of windows; each analysis contributes same effective number of tests.
    # Not quite the same as consolidateWindows()'s weighting, but that can't be helped as the adjusted p-values are computed before clustering.
    if (equiweight) { 
        weights <- rep(1/lengths(data.list), lengths(data.list))
    } else {
        weights <- NULL
    }   

    out <- clusterWindows(all.data, all.result, weight=weights, ...)
    out$id <- split(out$id, groupings)
    names(out$id) <- names(data.list)
    return(out)
}

