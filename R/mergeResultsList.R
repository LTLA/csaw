#' @export
#' @importFrom S4Vectors DataFrame metadata<- mcols I
mergeResultsList <- function(ranges.list, tab.list=NULL, tol, equiweight=TRUE, get.best=TRUE, 
    merge.args=list(), combine.args=list(), best.args=list()) 
{
    if (is.null(tab.list)) {
        tab.list <- lapply(ranges.list, mcols)
    }

    .verify_ranges_tabs(ranges.list, tab.list)
    m.out <- do.call(mergeWindowsList, c(list(ranges.list, tol=tol), merge.args))
    tab <- do.call(rbind, tab.list)

    combine.args$ids <- best.args$ids <- m.out$ids 
    combine.args$tab <- best.args$tab <- tab
    if (equiweight) {
        combine.args$weights <- best.args$weights <- m.out$weights
    }

    combined <- do.call(combineTests, combine.args)
    output <- DataFrame(regions=I(m.out$regions), combined=I(combined))
    if (get.best) {
        output$best <- do.call(getBestTest, best.args)
    }

    metadata(output) <- m.out[c("ranges", "ids", "weights")]
    metadata(output)$tab <- tab
    output
}

.verify_ranges_tabs <- function(ranges.list, tab.list) {
    all.len <- vapply(ranges.list, length, FUN.VALUE=0L)
    all.nrow <- vapply(tab.list, nrow, FUN.VALUE=0L)
    if (!identical(length(all.len), length(all.nrow))) {
        stop("'ranges.list' and 'tab.list' should have equal length")
    }
    if (!identical(all.len, all.nrow)) {
        stop("elements of 'ranges.list' and 'tab.list' should have equal length")
    }
}
