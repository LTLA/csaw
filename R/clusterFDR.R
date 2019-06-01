#' @export
clusterFDR <- function(ids, threshold, weights=NULL)
# This computes an informal estimate of the cluster-level FDR,
# given the cluster IDs for all significant windows. The idea
# is to allow clustering of significant windows to explicitly
# identify the differential subinterval in complex regions.
#
# written by Aaron Lun
# created 13 April 2015
{
    ids <- as.integer(ids)
    o <- order(ids)
	ids <- ids[o]

    if (is.null(weights)) { 
        weights <- rep(1, length(ids)) 
    } else {
        if (length(weights)!=length(ids)) { 
            stop("lengths of 'weights' and 'ids' must be the same")
        }
        weights <- as.double(weights)
        weights <- weights[o]
    }

	num.fp <- sum(weights) * threshold
	cluster.sizes <- tapply(weights, INDEX=ids, FUN=sum)
	num.fp.cluster <- sum(cumsum(sort(cluster.sizes)) <= num.fp)

    if (length(cluster.sizes)) { 
    	return(num.fp.cluster/length(cluster.sizes))
    } else {
        return(0)
    }
}

#' @export
controlClusterFDR <- function(target, adjp, FUN, ..., weights=NULL, grid.length=21, iterations=4)
# Identifies the window-level FDR threshold that is required to 
# control the cluster-level threshold at 'target', given the 
# window-level adjusted p-values and the clustering function FUN.
#
# written by Aaron Lun
# created 5 January 2016
{
    if (is.null(weights)) { 
        weights <- rep(1, length(adjp)) 
    } 

    # Doesn't make sense to have window-level FDR > cluster-level FDR. 
    # We'll be conservative and only search grid points that are lower.
    upper <- target
    lower <- 0

    # Using an iterative grid search, as this tends to be most robust for a discrete and discontinuous function.
    for (it in seq_len(iterations)) { 
        grid <- seq(lower, upper, length.out=grid.length)

        fdrs <- integer(length(grid))
        for (tx in seq_along(grid)) { 
            threshold <- grid[tx]
            is.sig <- adjp <= threshold
            fdrs[tx] <- clusterFDR(FUN(is.sig, ...), threshold, weights=weights[is.sig])
        }

        # Picking the highest grid point above which the FDR > target.
        # There must be at least one below the target, as FDR = 0 at threshold=0.
        controlled <- fdrs <= target
        chosen <- max(which(controlled))
        if (chosen==grid.length) {
            break
        } 
        lower <- grid[chosen]
        upper <- grid[chosen+1L]
    }

    return(list(threshold=grid[chosen], FDR=fdrs[chosen]))
}

