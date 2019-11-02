#' @export
clusterWindowsList <- function(ranges.list, tab.list, equiweight=TRUE, ...) 
# This does the same as clusterWindows, but for results from many different analyses
# (ostensibly with different window sizes).
#
# written by Aaron Lun
# created 8 January 2016
{
    .verify_ranges_tabs(ranges.list, tab.list)
    ranges.list <- lapply(ranges.list, .toGRanges)
    
    # Merging everyone together.
    all.data <- do.call(c, ranges.list)
    all.data$origin <- rep(seq_along(ranges.list), lengths(ranges.list))
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
