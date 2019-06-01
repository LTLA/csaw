#' @export
mergeWindowsList <- function(ranges.list, tol=100, ...) {
	all.ranges <- do.call(c, lapply(ranges.list, .toGRanges))
    mcols(all.ranges) <- NULL
	merged <- mergeWindows(all.ranges, tol=tol, ...)
    all.ranges$origin <- rep(seq_along(ranges.list), lengths(ranges.list))

    # Computing weights inversely proportional to the number of windows of each width in each cluster.
    last <- 0L
    all.weights <- numeric(length(all.ranges))
    for (i in seq_along(ranges.list)) {
        currows <- last + seq_along(ranges.list[[i]])
        curid <- merged$id[currows]
        all.weights[currows] <- (1/tabulate(curid))[curid]	
        last <- last + length(ranges.list[[i]])
    }

	list(ranges=all.ranges, ids=merged$id, regions=merged$regions, weights=all.weights)
}
