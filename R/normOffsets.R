setGeneric("normOffsets", function(object, ...) standardGeneric("normOffsets"))

setMethod("normOffsets", "matrix", function(object, lib.sizes=NULL, type=c("scaling", "loess"), weighted=FALSE, ...) 
# This provides a wrapper to perform TMM normalization with non-standard
# library sizes (e.g. due to filtering) and weighting turned off.
# Alternatively, it can do a form a fast loess-like normalization which uses
# the average count as the covariate, rather than the typical A-value-based
# shenanigans. This avoids instability at low abundances.
#
# written by Aaron Lun
# created 19 November 2013
# last modified 29 August 2015
{
	if (is.null(lib.sizes)) { 
		lib.sizes <- colSums(object) 
	}

	type <- match.arg(type)
	if (type=="scaling") { 
		y <- DGEList(object, lib.size=lib.sizes)
		y <- calcNormFactors(y, doWeighting=weighted, ...)
		return(y$samples$norm.factors)

	} else if (type=="loess") { 
		# Scaled corrections squeeze offsets towards relative log-library sizes.
		# Constant value of `cont.cor' would squeeze them towards zero.
		cont.cor <- 0.5
		cont.cor.scaled <- cont.cor * lib.sizes/mean(lib.sizes)
		
		# Using it as a prior.count for abundance ensures linearity with log-object.
		ab <- aveLogCPM(object, lib.size=lib.sizes, prior.count=cont.cor)/log2(exp(1))

		offs <- matrix(0, nrow(object), ncol(object), byrow=TRUE)
		for (x in seq_len(ncol(object))) {
			fit <- loessFit(log(object[,x]+cont.cor.scaled[x]), ab, ...)
			offs[,x] <- fit$fitted 
		}
		offs <- offs-rowMeans(offs)
		return(offs)
	}
})

setMethod("normOffsets", "SummarizedExperiment", function(object, assay=1, type="scaling", ..., se.out=TRUE) {
    if (is.null(object$totals)) { 
        stop("missing 'totals' from SummarizedExperiment")
    } 
    lib.sizes <- object$totals 
	
    out <- normOffsets(assay(object, assay), lib.sizes=lib.sizes, type=type, ...)
    
    if (!is.logical(se.out)) { 
        if (type!="scaling") {
            stop("alternative output object not supported for loess normalization")
        }
        if (!identical(se.out$totals, lib.sizes)) {
            stop("library sizes of 'se.out' and 'object' are not identical")
        }
        object <- se.out
        se.out <- TRUE
    }

    if (!se.out) {
        return(out)
    } else {
        if (type=="scaling") {
            object$norm.factors <- out
        } else {
            assay(object, "offset") <- out
        }
        return(object)
    }
})
