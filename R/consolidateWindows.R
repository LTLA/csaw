#' @export
#' @importFrom IRanges findOverlaps
#' @importFrom S4Vectors queryHits
consolidateWindows <- function(data.list, equiweight=TRUE, 
    merge.args=list(), region=NULL, overlap.args=list()) 
# Consolidates results for multiple window sizes into a result for the
# genomic region over which those windows are tiled. Returns the combined
# results, as well as ID vectors for cross-referencing and inspection.
#
# written by Aaron Lun
# created 26 February 2015
{
	nset <- length(data.list)
	for (x in seq_len(nset)) {
        data.list[[x]] <- .toGRanges(data.list[[x]])
	}
    rel.weights <- NULL

	if (is.null(region)) { 
        merge.call <- do.call(call, c("mergeWindows", merge.args))
        merge.call <- match.call(mergeWindows, merge.call)
        merge.args <- as.list(merge.call)[-1]
        if (is.null(merge.args$tol)) { 
            merge.args$tol <- 100
            warning("'tol' for 'mergeWindows' set to a default of 100 bp")
        }

		all.ranges <- do.call(c, data.list)
		merged <- do.call(mergeWindows, c(merge.args, regions=all.ranges)) 

		# Formatting for nice output.
		final.ids <- vector("list", nset)
		last <- 0L
		for (x in seq_len(nset)) { 
            currows <- length(data.list[[x]])
			final.ids[[x]] <- merged$id[last + seq_len(currows)]
			last <- last + currows
		}

        if (equiweight) {
            # Computing weights inversely proportional to the number of windows of each width in each cluster.
    		rel.weights <- vector("list", nset)
            for (x in seq_len(nset)) {
                curid <- final.ids[[x]]
                rel.weights[[x]] <- (1/tabulate(curid))[curid]	
            }
            names(rel.weights) <- names(data.list)
		}
        
        names(final.ids) <- names(data.list)
        return(list(id=final.ids, weight=rel.weights, region=merged$region))

	} else {
		final.olaps <- vector("list", nset)
		for (x in seq_len(nset)) {
			final.olaps[[x]] <- do.call(findOverlaps, c(query=region, subject=data.list[[x]], overlap.args))
		}

        if (equiweight) {
            # Computing weights inversely proportional to the number of windows of each width in each region.
    		rel.weights <- vector("list", nset)
            for (x in seq_len(nset)) {
                curlap <- final.olaps[[x]]
                curid <- queryHits(curlap)
                rel.weights[[x]] <- (1/tabulate(curid))[curid]	
            }
            names(rel.weights) <- names(data.list)
        }
	
        names(final.olaps) <- names(data.list)
        return(list(olap=final.olaps, weight=rel.weights))
	}
}

consolidateSizes <- function(data.list, result.list, equiweight=TRUE, 
    merge.args=list(), region=NULL, overlap.args=list(), ...) 
# Deprecated function that calls consolidateWindows.    
{
    .Deprecated(new="consolidateWindows")
    out <- consolidateWindows(data.list, equiweight=equiweight, 
        merge.args=merge.args, region=region, overlap.args=overlap.args)
    
    if (!is.null(region)) {
        tabcom <- consolidateTests(olap.list=out$olap, result.list=result.list,
            weight.list=out$weight, ...) 
    } else {
        tabcom <- consolidateTests(id.list=out$id, result.list=result.list,
            weight.list=out$weight, ...) 
    }
   
    out$table <- tabcom
    out$weight <- NULL
    return(out)
}

