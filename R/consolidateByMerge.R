#' @export
consolidateByMerge <- function(x, tol=100, sign.list=NULL, ...) {
	if (!is.null(sign.list)) {
		slen <- lengths(sign.list)
		dlen <- lengths(x)
		if (length(slen)!=length(dlen) || any(slen!=dlen)) { 
			stop("vector lengths of 'sign.list' and 'x' are not identical")
		}
		signs <- unlist(sign.list)
	} else {
        signs <- NULL
    }

	all.ranges <- do.call(c, lapply(x, .toGRanges))
    mcols(all.ranges) <- NULL
	merged <- mergeWindows(all.ranges, sign=signs, tol=tol, ...)
    all.ranges$origin <- rep(seq_along(x), lengths(x))
    all.ranges$id <- merged$id

    # Computing weights inversely proportional to the number of windows of each width in each cluster.
    last <- 0L
    all.weights <- numeric(length(all.ranges))
    for (i in seq_along(x)) {
        currows <- last + seq_along(x[[i]])
        curid <- merged$id[currows]
        all.weights[currows] <- (1/tabulate(curid))[curid]	
        last <- last + length(x[[i]])
    }
    all.ranges$weight <- all.weights

	list(ranges=all.ranges, merged=merged$region)
}
