#' @export
#' @importFrom IRanges findOverlaps
#' @importFrom S4Vectors queryHits subjectHits
findOverlapsList <- function(x, ref, ...) {
	all.ranges <- do.call(c, lapply(x, .toGRanges))
    mcols(all.ranges) <- NULL
    all.ranges$origin <- rep(seq_along(x), lengths(x))
    olaps <- findOverlaps(query=ref, subject=all.ranges, ...)

    # Computing weights inversely proportional to the number of windows of each width in each region.
    by.origin <- split(seq_along(olaps), all.ranges$origin[subjectHits(olaps)])
    rel.weights <- numeric(length(olaps))
    for (i in seq_along(by.origin)) {
        indices <- by.origin[[i]]
        curid <- queryHits(olaps)[indices]
        rel.weights[indices] <- (1/tabulate(curid))[curid]	
    }

    list(ranges=all.ranges, olap=olaps, weight=rel.weights)
}

