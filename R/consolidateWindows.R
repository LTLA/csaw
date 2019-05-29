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
        output <- do.call(consolidateByMerge, c(list(data.list, equiweight=equiweight, sign.list=sign.list), merge.args))
        list(
            id=split(output$original$id, output$original$origin), 
            weight=split(output$original$weight, output$original$origin),
            region=output$region
        )

	} else {
        .Deprecated(new="consolidateByOverlap")
        output <- do.call(consolidateByOverlap, c(list(data.list, equiweight=equiweight, region=region), overlap.args))

        original.index <- unlist(lapply(data.list, seq_along))
        subject.origin <- output$original$origin[subjectHits(output$olap)]

        by.original <- split(output$olap, subject.origin)
        for (i in seq_along(by.original)) {
            current <- by.original[[i]]
            by.original[[i]] <- Hits(queryHits(current), original.index[subjectHits(current)],
                nLnode=nLnode(current), nRnode=nRnode(current), sort.by.query=TRUE)
        }

        list(
            olap=by.original,
            weight=split(output$weight, subject.origin)
        )
	}
}
