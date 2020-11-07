#' @export
#' @importFrom S4Vectors DataFrame metadata<- I
#' @importFrom IRanges findOverlaps
overlapResults <- function(ranges, tab=mcols(ranges), regions, get.best=TRUE, overlap.args=list(), combine.args=list(), best.args=list()) {
    olap <- do.call(findOverlaps, c(list(query=regions, subject=ranges), overlap.args))
    combined <- do.call(combineOverlaps, c(list(olap, tab), combine.args))
    output <- DataFrame(regions=I(regions), combined=I(combined))
    if (get.best) {
        output$best <- do.call(getBestOverlaps, c(list(olap, tab), best.args))
    }
    metadata(output)$overlaps <- olap
    output
}
