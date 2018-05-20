#' @export
setGeneric("asDGEList", function(object, ...) { standardGeneric("asDGEList") })

#' @export
#' @importFrom edgeR DGEList scaleOffset
#' @importFrom SummarizedExperiment assay assayNames
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
setMethod("asDGEList", "SummarizedExperiment", function(object, lib.sizes, norm.factors, assay.id="counts", ...) 
# This defines a wrapper function to convert a SummarizedExperiment class
# object into a DGEList object for input into edgeR.
#
# written by Aaron Lun
# created 2 September 2014
{
	all.args <- list(...)
	if (missing(lib.sizes)) { 
		if (is.null(object$totals)) { warning("library sizes not found in 'totals', setting to NULL") }
		lib.sizes <- object$totals
	}
	all.args$lib.size <- lib.sizes

	if (missing(norm.factors)) {
		if (!is.null(object$norm.factors)) { 
			all.args$norm.factors <- object$norm.factors 
		}
	} else {
		all.args$norm.factors <- norm.factors 
	}

    all.args$counts <- assay(object, i=assay.id, withDimnames=FALSE) 
	y <- do.call(DGEList, all.args)

    if ("offset" %in% assayNames(object)) {
        offset <- assay(object, i="offset", withDimnames=FALSE)
        y <- scaleOffset(y, offset)
    }
	return(y)
})
