#' @export
#' @importFrom S4Vectors DataFrame metadata<- I
#' @importFrom IRanges findOverlaps
overlapResultsList <- function(ranges.list, tab.list=NULL, regions, equiweight=TRUE, get.best=TRUE, 
    overlap.args=list(), combine.args=list(), best.args=list()) 
{
    if (is.null(tab.list)) {
        tab.list <- lapply(ranges.list, mcols)
    }

    .verify_ranges_tabs(ranges.list, tab.list)
    olap.out <- do.call(findOverlapsList, c(list(ranges.list, regions), overlap.args))
    tab <- do.call(rbind, tab.list)

    combine.args$overlaps <- best.args$overlaps <- olap.out$overlaps
    combine.args$tab <- best.args$tab <- tab
    if (equiweight) {
        combine.args$o.weights <- best.args$o.weights <- olap.out$weights
    }

    combined <- do.call(combineOverlaps, combine.args)
    output <- DataFrame(regions=I(regions), combined=I(combined))
    if (get.best) {
        output$best <- do.call(getBestOverlaps, best.args)
    }

    metadata(output) <- olap.out
    metadata(output)$tab <- tab
    output
}
