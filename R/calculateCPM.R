#' @export
#' @importFrom edgeR scaleOffset.default cpm addPriorCount
#' @importFrom SummarizedExperiment assay
calculateCPM <- function(object, use.norm.factors=TRUE, use.offsets=FALSE, 
        log=TRUE, prior.count=1, assay.id="counts")
# This is a convenience wrapper to compute CPMs, avoiding the
# need to convert to a DGEList and then call edgeR functions.
# It also provides a more convenient way to compute offset-adjusted functions.
#
# written by Aaron Lun
# created 18 May 2018
{
    mat <- assay(object, i=assay.id)
    lib.size <- object$totals
    if (is.null(object$totals)) { 
        stop("missing 'totals' from SummarizedExperiment 'colData'")
    } 
    
    if (!use.offsets) {
        # Computing effective library sizes.
        if (!is.null(object$norm.factors) && use.norm.factors) {
            lib.size <- lib.size * object$norm.factors
        }

        return(cpm(mat, lib.size=lib.size, prior.count=prior.count, log=log))
    }

    # Computing a offset-adjusted CPM matrix.
    off.mat <- assay(object, i="offset", withDimnames=FALSE)
    new.offset <- scaleOffset.default(lib.size, off.mat)
    if (!log) {
        cpm.out <- mat/exp(new.offset) * 1e6
    } else {
        ap <- addPriorCount(mat, offset=new.offset, prior.count=prior.count)
        cpm.out <- log2(ap$y) - as.matrix(ap$offset)/log(2) + log2(1e6)
    }

    dimnames(cpm.out) <- dimnames(mat)
    return(cpm.out)
}
