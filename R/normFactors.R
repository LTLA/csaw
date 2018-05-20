normFactors <- function(object, weighted=FALSE, ..., assay.id="counts", se.out=TRUE) 
# This provides a wrapper to perform TMM normalization with non-standard
# library sizes (e.g. due to filtering) and weighting turned off.
# Alternatively, it can do a form a fast loess-like normalization which uses
# the average count as the covariate, rather than the typical A-value-based
# shenanigans. This avoids instability at low abundances.
#
# written by Aaron Lun
# created 18 May 2018
{
    lib.sizes <- object$totals
    if (is.null(lib.sizes)) {
        stop("library sizes not present in 'object$totals'")
    }

	out <- calcNormFactors(assay(object, i=assay.id, withDimnames=FALSE), 
        lib.size=lib.sizes, doWeighting=weighted, ...)

	# Choosing to put these values in a different location, if requested. 
    if (!is.logical(se.out)) { 
        if (ncol(object)!=ncol(se.out)) {
            stop("number of libraries differs between 'se.out' and 'object'")
        }

        if (!identical(se.out$totals, lib.sizes)) {
            stop("library sizes of 'se.out' and 'object' are not identical")
        }

        object <- se.out
        se.out <- TRUE
    }

    if (!se.out) {
        return(out)
    } 
    object$norm.factors <- out
    return(object)
}
