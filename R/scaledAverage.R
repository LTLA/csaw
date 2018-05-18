scaledAverage <- function(y, scale=1, prior.count=NULL, dispersion=NULL, assay.id="counts")
# This computes the scaled average abundance, with some finesse to deal with
# the interaction between scaling and the prior count. The `scale` factor
# represents the downscaling factor for the abundances, so the prior count
# has to be upscaled during the calculation itself.
#
# written by Aaron Lun
# created 5 November 2014
# last modified 3 March 2017
{
    if (is(y, "DGEList")) {
        .Deprecated(msg="DGEList inputs to scaledAverage are deprecated.\nUse SummarizedExperiment inputs instead.")
        counts <- y$counts
        common.disp <- y$common.dispersion 
        lib.size <- y$samples$lib.size * y$samples$norm.factors
    } else if (is(y, "SummarizedExperiment")) {
        counts <- assay(y, i=assay.id, withDimnames=FALSE)
        common.disp <- NULL
        lib.size <- y$totals
        if (!is.null(y$norm.factors)) {
            lib.size <- lib.size * y$norm.factors
        }
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

getWidths <- function(data) 
# This computes the effective width of the data in the
# RangedSummarizedExperiment object. This is done by accounting for the effect
# of read extension; or, for paired end data, the median fragment length.
#
# written by Aaron Lun
# created 5 November 2014
# last modified 17 December 2015
{
	flen <- metadata(data)$final.ext

	if (is.na(flen)) { 
		flen <- data$ext
		is.missing <- is.na(flen)
		if (any(is.missing)) { 
			if (is.null(data$rlen) || any(is.na(data$rlen[is.missing]))) { 
				stop("need to specify read lengths in 'data$rlen'")
			}
			flen[is.missing] <- data$rlen[is.missing]
		}
		flen <- as.integer(mean(flen))
	}

	width(rowRanges(data)) + flen - 1L
}

