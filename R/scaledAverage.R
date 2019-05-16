#' @export
#' @importFrom methods is
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
#' @importFrom SummarizedExperiment assay
#' @importFrom edgeR aveLogCPM aveLogCPM.DGEList addPriorCount mglmOneGroup
scaledAverage <- function(y, scale=1, prior.count=NULL, dispersion=NULL, assay.id="counts")
# This computes the scaled average abundance, with some finesse to deal with
# the interaction between scaling and the prior count. The `scale` factor
# represents the downscaling factor for the abundances, so the prior count
# has to be upscaled during the calculation itself.
#
# written by Aaron Lun
# created 5 November 2014
{
    counts <- assay(y, i=assay.id, withDimnames=FALSE)
    common.disp <- NULL
    lib.size <- y$totals
    if (!is.null(y$norm.factors)) {
        lib.size <- lib.size * y$norm.factors
    }

    # Setting values.
    if (is.null(dispersion)) { 
        dispersion <- common.disp
        if (is.null(dispersion)) dispersion <- 0.05
    }
	if (is.null(prior.count)) { prior.count <- formals(aveLogCPM.DGEList)$prior.count }

    # Checking validity of scale
    is.zero <- scale==0
    is.neg <- scale < 0
    failed <- is.zero | is.neg
	if (any(failed)) { 
        scale[failed] <- 1
    }

    # Computing the prior to add, scaling it, and then adding it.
    # This way ensures that the offsets are constant regardless of 'scale'.
    empty <- matrix(0, nrow(counts), ncol(counts))
    ap <- addPriorCount(empty, lib.size=lib.size, prior.count=prior.count)
    ap$y <- ap$y * scale
    ap$y <- ap$y + counts

    # Computing the average abundances in ave-logCPM.
	ave <- mglmOneGroup(y=ap$y, offset=ap$offset, dispersion=dispersion)
    ave <- (ave - log(scale) + log(1e6))/log(2)

    # Replacing values with invalid scale.
    if (any(failed)) {
        if (length(failed)==1) { 
            ave[] <- ifelse(is.zero, -Inf, NA_real_)
        } else {
            ave[is.neg] <- NA_real_
            ave[is.zero] <- -Inf
        }
    } 
    return(ave)
}
