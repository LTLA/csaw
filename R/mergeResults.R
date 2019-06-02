#' @export
#' @importFrom S4Vectors DataFrame
mergeResults <- function(ranges, tab, tol, get.best=TRUE, merge.args=list(), combine.args=list(), best.args=list()) {
    merged <- do.call(mergeWindows, c(list(ranges, tol=tol), merge.args))
    combined <- do.call(combineTests, c(list(merged$id, tab), combine.args))
    output <- DataFrame(regions=I(merged$regions), combined=I(combined))
    if (get.best) {
        output$best <- do.call(getBestTest, c(list(merged$id, tab), best.args))
    }
    output
}
