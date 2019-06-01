#' @importFrom S4Vectors queryLength queryHits subjectHits
.overlapStats <- function(overlaps, tab, o.weights=NULL, i.weights=NULL, FUN, ...) { 
	region.dex <- queryHits(overlaps)
	win.dex <- subjectHits(overlaps)

	# Setting up weights.
	if (is.null(o.weights)) { 
		if (!is.null(i.weights)) {
			o.weights <- i.weights[win.dex]
		}
	}

    output <- FUN(region.dex, tab[win.dex,], weights=o.weights, ...)

	N <- queryLength(overlaps)
    expand.vec <- rep(NA_integer_, N)
    row.dex <- as.integer(rownames(output))
    if (any(N <= 0L | row.dex > N)) { 
        stop("IDs are not within [1, nregions]") 
    }

    expand.vec[row.dex] <- seq_along(row.dex)
    output <- output[expand.vec,,drop=FALSE]
    rownames(output) <- as.character(seq_len(N))
    output
}

#' @export
combineOverlaps <- function(overlaps, tab, o.weights=NULL, i.weights=NULL, ...) 
# Wrapper around combineTests for Hits from findOverlaps,
# when windows are overlapped with regions.
#
# written by Aaron Lun
# created 25 March 2015
# last modified 26 March 2015
{
	.overlapStats(overlaps, tab, o.weights=o.weights, i.weights=i.weights, FUN=combineTests, ...)
}

#' @export
#' @importFrom S4Vectors subjectHits
getBestOverlaps <- function(overlaps, tab, o.weights=NULL, i.weights=NULL, ...) 
# Wrapper around getBestTest for Hits from findOverlaps,
# when windows are overlapped with regions.
#
# written by Aaron Lun
# created 25 March 2015
# last modified 26 March 2015
{
	output <- .overlapStats(overlaps, tab, o.weights=o.weights, i.weights=i.weights, FUN=getBestTest, ...)
    output$best <- subjectHits(overlaps)[output$best]
    output
}

#' @export
empiricalOverlaps <- function(overlaps, tab, o.weights=NULL, i.weights=NULL, ...) 
# Wrapper around empiricalFDR for Hits from findOverlaps,
# when windows are overlapped with regions
#
# written by Aaron Lun
# created 7 January 2017
{
    .overlapStats(overlaps, tab, o.weights=o.weights, i.weights=i.weights, FUN=empiricalFDR, ...)
}

#' @export
mixedOverlaps <- function(overlaps, tab, o.weights=NULL, i.weights=NULL, ...) 
# Wrapper around mixedClusters for Hits from findOverlaps,
# when windows are overlapped with regions
#
# written by Aaron Lun
# created 7 January 2017
{
    .overlapStats(overlaps, tab, o.weights=o.weights, i.weights=i.weights, FUN=mixedClusters, ...)
}

#' @export
#' @importFrom S4Vectors queryHits subjectHits 
summitOverlaps <- function(overlaps, region.best, o.summit=NULL, i.summit=NULL)
# Wrapper around upweightSummits for Hits from findOverlaps.
#
# written by Aaron Lun
# created 25 March 2015
# last modified 26 March 2015
{
	region.dex <- queryHits(overlaps)
	win.dex <- subjectHits(overlaps)

	if (!missing(region.best)) { 
		summit.dex <- region.best[region.dex]
		summits <- !is.na(summit.dex) & win.dex==summit.dex

	} else if (!is.null(o.summit)) {
		if (is.integer(o.summit)) { 
			out <- logical(length(overlaps))
			out[o.summit] <- TRUE
			o.summit <- out
		} else {
			stopifnot(length(o.summit)==length(overlaps))
		}	
		summits <- o.summit

 	} else if (!is.null(i.summit)) { 	
		if (is.integer(i.summit)) { 
			out <- logical(max(win.dex, i.summit))
			out[i.summit] <- TRUE
			i.summit <- out
		}
		summits <- i.summit[win.dex]

	} else {
		stop("either region.best, i.summit or o.summit must be specified")
	}

	upweightSummit(region.dex, summits)
}



