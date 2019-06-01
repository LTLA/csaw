#' @export
#' @importFrom S4Vectors DataFrame metadata<-
mergeResultsList <- function(ranges.list, tab.list, tol, equiweight=TRUE, get.best=TRUE, 
    merge.args=list(), combine.args=list(), best.args=list()) 
{
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

