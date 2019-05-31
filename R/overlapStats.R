#' @importFrom S4Vectors queryLength queryHits subjectHits
.overlapStats <- function(olap, tab, o.weight=NULL, i.weight=NULL, FUN, ...) { 
	region.dex <- queryHits(olap)
	win.dex <- subjectHits(olap)

	# Setting up weights.
	if (is.null(o.weight)) { 
		if (!is.null(i.weight)) {
			o.weight <- i.weight[win.dex]
		}
	}

    output <- FUN(region.dex, tab[win.dex,], weight=o.weight, ...)

	N <- queryLength(olap)
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
combineOverlaps <- function(olap, tab, o.weight=NULL, i.weight=NULL, ...) 
# Wrapper around combineTests for Hits from findOverlaps,
# when windows are overlapped with regions.
#
# written by Aaron Lun
# created 25 March 2015
# last modified 26 March 2015
{
	.overlapStats(olap, tab, o.weight=o.weight, i.weight=i.weight, FUN=combineTests, ...)
}

#' @export
getBestOverlaps <- function(olap, tab, o.weight=NULL, i.weight=NULL, ...) 
# Wrapper around getBestTest for Hits from findOverlaps,
# when windows are overlapped with regions.
#
# written by Aaron Lun
# created 25 March 2015
# last modified 26 March 2015
{
	.overlapStats(olap, tab, o.weight=o.weight, i.weight=i.weight, FUN=getBestTest, ...)
}

#' @export
empiricalOverlaps <- function(olap, tab, o.weight=NULL, i.weight=NULL, ...) 
# Wrapper around empiricalFDR for Hits from findOverlaps,
# when windows are overlapped with regions
#
# written by Aaron Lun
# created 7 January 2017
{
    .overlapStats(olap, tab, o.weight=o.weight, i.weight=i.weight, FUN=empiricalFDR, ...)
}

#' @export
mixedOverlaps <- function(olap, tab, o.weight=NULL, i.weight=NULL, ...) 
# Wrapper around mixedClusters for Hits from findOverlaps,
# when windows are overlapped with regions
#
# written by Aaron Lun
# created 7 January 2017
{
    .overlapStats(olap, tab, o.weight=o.weight, i.weight=i.weight, FUN=mixedClusters, ...)
}

#' @export
#' @importFrom S4Vectors queryHits subjectHits 
summitOverlaps <- function(olap, region.best, o.summit=NULL, i.summit=NULL)
# Wrapper around upweightSummits for Hits from findOverlaps.
#
# written by Aaron Lun
# created 25 March 2015
# last modified 26 March 2015
{
	region.dex <- queryHits(olap)
	win.dex <- subjectHits(olap)

	if (!missing(region.best)) { 
		summit.dex <- region.best[region.dex]
		summits <- !is.na(summit.dex) & win.dex==summit.dex

	} else if (!is.null(o.summit)) {
		if (is.integer(o.summit)) { 
			out <- logical(length(olap))
			out[o.summit] <- TRUE
			o.summit <- out
		} else {
			stopifnot(length(o.summit)==length(olap))
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



