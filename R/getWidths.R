#' @export
#' @importFrom BiocGenerics width
#' @importFrom S4Vectors metadata
#' @importFrom SummarizedExperiment rowRanges
getWidths <- function(data) 
# This computes the effective width of the data in the
# RangedSummarizedExperiment object. This is done by accounting for the effect
# of read extension; or, for paired end data, the median fragment length.
#
# written by Aaron Lun
# created 5 November 2014
{
	flen <- metadata(data)$final.ext

    if (is.null(flen)) {
        flen <- 1L
    } else if (is.na(flen)) { 
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

