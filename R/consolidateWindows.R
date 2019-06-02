#' @export
#' @importFrom IRanges findOverlaps
#' @importFrom S4Vectors queryHits
consolidateWindows <- function(data.list, equiweight=TRUE, merge.args=list(), regions=NULL, overlap.args=list(), sign.list=NULL) 
# Consolidates results for multiple window sizes into a result for the
# genomic region over which those windows are tiled. Returns the combined
# results, as well as ID vectors for cross-referencing and inspection.
#
# written by Aaron Lun
# created 26 February 2015
{
	if (is.null(regions)) { 
        .Deprecated(new="mergeWindowsList")
        output <- do.call(mergeWindowsList, c(list(data.list, signs=unlist(sign.list)), merge.args))

        weight <- if (equiweight) {
            split(output$weights, output$ranges$origin)
        } else {
            NULL
        }

        list(
            id=split(output$ids, output$ranges$origin), 
            weight=weight, 
            region=output$regions
        )

	} else {
        .Deprecated(new="findOverlapsList")
        output <- do.call(findOverlapsList, c(list(data.list, regions=regions), overlap.args))

        ranges.index <- unlist(lapply(data.list, seq_along))
        subject.origin <- output$ranges$origin[subjectHits(output$overlaps)]

        by.ranges <- vector("list", length(data.list))
        for (i in seq_along(by.ranges)) {
            current <- output$overlaps[subject.origin==i]
            by.ranges[[i]] <- Hits(queryHits(current), ranges.index[subjectHits(current)],
                nLnode=nLnode(current), nRnode=length(data.list[[i]]), sort.by.query=TRUE)
        }

        weight <- if (equiweight) {
            split(output$weights, subject.origin)
        } else {
            NULL
        }

        list(olap=by.ranges, weight=weight)
	}
}
