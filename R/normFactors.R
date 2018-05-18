setGeneric("normFactors", function(object, ...) standardGeneric("normFactors"))

setMethod("normFactors", "matrix", function(object, lib.sizes=NULL, weighted=FALSE, ...) 
# This provides a wrapper to perform TMM normalization with non-standard
# library sizes (e.g. due to filtering) and weighting turned off.
# Alternatively, it can do a form a fast loess-like normalization which uses
# the average count as the covariate, rather than the typical A-value-based
# shenanigans. This avoids instability at low abundances.
#
# written by Aaron Lun
# created 18 May 2018
{
	if (is.null(lib.sizes)) { 
		lib.sizes <- colSums(object) 
	}
	calcNormFactors(object, lib.size=lib.sizes, doWeighting=weighted, ...)
})

setMethod("normFactors", "SummarizedExperiment", function(object, assay.id="counts", ..., se.out=TRUE) {
    if (is.null(object$totals)) { 
        stop("missing 'totals' from SummarizedExperiment")
    } 
    lib.sizes <- object$totals 
    out <- normFactors(assay(object, i=assay.id, withDimnames=FALSE), lib.sizes=lib.sizes, ...)
   
	# Choosing to put these values in a different location, if requested. 
    if (!is.logical(se.out)) { 
        if (!identical(se.out$totals, lib.sizes)) {
            stop("library sizes of 'se.out' and 'object' are not identical")
        }
        object <- se.out
        se.out <- TRUE
    }

    if (!se.out) {
        return(out)
    } else {
        object$norm.factors <- out
        return(object)
    }
})
