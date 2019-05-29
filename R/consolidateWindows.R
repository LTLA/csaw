#' @export
#' @importFrom IRanges findOverlaps
#' @importFrom S4Vectors queryHits
consolidateWindows <- function(data.list, equiweight=TRUE, merge.args=list(), region=NULL, overlap.args=list(), sign.list=NULL) 
# Consolidates results for multiple window sizes into a result for the
# genomic region over which those windows are tiled. Returns the combined
# results, as well as ID vectors for cross-referencing and inspection.
#
# written by Aaron Lun
# created 26 February 2015
{
	if (is.null(region)) { 
        .Deprecated(new="consolidateByMerge")
        output <- do.call(consolidateByMerge, c(list(data.list, sign.list=sign.list), merge.args))

        weight <- if (equiweight) {
            split(output$ranges$weight, output$ranges$origin)
        } else {
            NULL
        }

        list(
            id=split(output$ranges$id, output$ranges$origin), 
            weight=weight, 
            region=output$merged
        )

	} else {
        .Deprecated(new="consolidateByOverlap")
        output <- do.call(consolidateByOverlap, c(list(data.list, equiweight=equiweight, region=region), overlap.args))

        ranges.index <- unlist(lapply(data.list, seq_along))
        subject.origin <- output$ranges$origin[subjectHits(output$olap)]

        by.ranges <- split(output$olap, subject.origin)
        for (i in seq_along(by.ranges)) {
            current <- by.ranges[[i]]
            by.ranges[[i]] <- Hits(queryHits(current), ranges.index[subjectHits(current)],
                nLnode=nLnode(current), nRnode=nRnode(current), sort.by.query=TRUE)
        }

        weight <- if (equiweight) {
            split(output$olap$weight, subject.origin)
        } else {
            NULL
        }

        list(olap=by.ranges, weight=weight)
	}
}
