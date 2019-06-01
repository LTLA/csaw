#' @export
clusterWindows <- function(ranges, tab, target, pval.col=NULL, fc.col=NA, tol, ..., weights=NULL, grid.length=21, iterations=4)
# This does a search for the clusters based on DB windows. 
# It aims to achieve a cluster-level FDR of 'target'.
#
# written by Aaron Lun
# created 8 January 2016
{
    ranges <- .toGRanges(ranges)
    if (nrow(tab)!=length(ranges)) { 
        stop("number of ranges is not consistent with entries in 'tab'") 
    }
    if (missing(tol)) {
        tol <- 100
        warning("'tol' for 'mergeWindows' set to a default of 100 bp")
    }
    if (missing(target)) { 
        target <- 0.05
        warning("unspecified 'target' for the cluster-level FDR set to 0.05")
    }

    # Computing a frequency-weighted adjusted p-value.
    pval.col <- .getPValCol(pval.col, tab)
    if (is.null(weights)) { 
        weights <- rep(1, nrow(tab)) 
    }
    adjp <- .weightedFDR(tab[,pval.col], weights)

    # Getting the sign.
    if (is.na(fc.col)) { 
        sign <- NULL
    } else {
        sign <- tab[,fc.col] > 0
    }

    # Controlling the cluster-level FDR
    FUN <- function(sig) { 
        mergeWindows(ranges[sig], tol=tol, signs=sign[sig], ...) 
    }
    out <- controlClusterFDR(target=target, adjp=adjp, 
        FUN=function(sig) { FUN(sig)$ids }, 
        weights=weights, grid.length=grid.length, iterations=iterations)
    sig <- adjp <= out$threshold
    clusters <- FUN(sig)

    # Cleaning up
    full.ids <- rep(NA_integer_, nrow(tab))
    full.ids[sig] <- clusters$ids
    clusters$ids <- full.ids
    clusters$FDR <- out$FDR 
    clusters
}

.weightedFDR <- function(p, w) {
    if (length(p)!=length(w)) { stop("weights and p-value vector are not of same length") }
    o <- order(p)
    p <- p[o]
    w <- w[o]
    adjp <- numeric(length(o))
    adjp[o] <- rev(cummin(rev(sum(w)*p/cumsum(w))))
    pmin(adjp, 1)
}
