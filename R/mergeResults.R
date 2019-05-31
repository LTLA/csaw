#' @export
#' @importFrom S4Vectors DataFrame
mergeResults <- function(regions, tab, tol, get.best=TRUE, merge.args=list(), combine.args=list(), best.args=list()) {
    merged <- do.call(mergeWindows, c(list(regions, tol=tol), merge.args))
    combined <- do.call(combineTests, c(list(merged$id, tab), combine.args))
    output <- DataFrame(merged=I(merged$region), combined=I(combined))
    if (get.best) {
        output$best <- do.call(getBestTest, c(list(merged$id, tab), best.args))
    }
    output
}

#' @export
#' @importFrom S4Vectors DataFrame
#' @importFrom IRanges findOverlaps
overlapResults <- function(regions, tab, ref, get.best=TRUE, overlap.args=list(), combine.args=list(), best.args=list()) {
    olap <- do.call(findOverlaps, c(list(query=ref, subject=regions), overlap.args))
    combined <- do.call(combineOverlaps, c(list(olap, tab), combine.args))
    output <- DataFrame(ref=I(ref), combined=I(combined))
    if (get.best) {
        output$best <- do.call(getBestOverlaps, c(list(olap, tab), best.args))
    }
    output
}
