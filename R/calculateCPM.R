calculateCPM <- function(object, use.norm.factors=TRUE, use.offsets=FALSE, prior.count=1, log=TRUE, assay.id="counts")
# This is a convenience wrapper to compute CPMs, avoiding the
# need to convert to a DGEList and then call edgeR functions.
# It also provides a more convenient way to compute offset-adjusted functions.
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
    if (!log) {
        new.offset <- scaleOffset(lib.size, off.mat)
        cpm.out <- mat/exp(new.offset - log(1e6))
        
    } else {
        ap <- addPriorCount(mat, offset=off.mat, prior.count=prior.count)
        new.offset <- scaleOffset(lib.size, as.matrix(ap$offset))/log(2)
        cpm.out <- log2(ap$y) - new.offset + log2(1e6)
    }

    dimnames(cpm.out) <- dimnames(mat)
    return(cpm.out)
}
