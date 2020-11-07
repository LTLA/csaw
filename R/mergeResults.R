#' @export
#' @importFrom S4Vectors DataFrame metadata<- mcols I
mergeResults <- function(ranges, tab=mcols(ranges), tol, get.best=TRUE, merge.args=list(), combine.args=list(), best.args=list()) {
    merged <- do.call(mergeWindows, c(list(ranges, tol=tol), merge.args))
    combined <- do.call(combineTests, c(list(merged$ids, tab), combine.args))
    output <- DataFrame(regions=I(merged$regions), combined=I(combined))
    if (get.best) {
        output$best <- do.call(getBestTest, c(list(merged$id, tab), best.args))
    }
    metadata(output)$ids <- merged$ids
    output
}
