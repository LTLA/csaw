#' Scaling normalization across libraries
#' 
#' Calculate normalization factors using count data from multiple libraries.
#' 
#' @param object A \linkS4class{SummarizedExperiment} object containing a count matrix and library sizes in the \code{totals} field of the \code{\link{colData}}.
#'
#' Alternatively, a \link{DGEList} object containing a count matrix in \code{object$counts} and library sizes in \code{object$samples$lib.size}.
#'
#' Alternatively, an ordinary matrix containing counts.
#' @param method Deprecated argument, ignored.
#' @param weighted A logical scalar indicating whether precision weights should be used for TMM normalization.
#' @param ... Other arguments to be passed to \code{\link{calcNormFactors}}.
#' @param assay.id An integer scalar or string specifying the assay values to use for normalization.
#' @param se.out A logical scalar indicating whether or not a SummarizedExperiment object should be returned.
#'
#' Alternatively, a SummarizedExperiment or DGEList object in which normalization factors are to be stored.
#' 
#' @details
#' This function uses the trimmed mean of M-values (TMM) method to remove composition biases, typically in background regions of the genome.
#' The key difference from standard TMM is that precision weighting is turned off by default so as to avoid upweighting high-abundance regions. 
#' These are more likely to be bound and thus more likely to be differentially bound. 
#' Assigning excessive weight to such regions will defeat the purpose of trimming when normalizing the coverage of background regions.
#' 
#' % Large changes in the normalization factor estimates with changes in the prior suggest that the counts are too low i.e. not
#' % enough new information in the dataset. This can be overcome by (obviously) increasing the counts. For example, binning
#' % can be performed with a larger bin size in \code{windowCounts} to obtain proportionally larger counts.
#' 
#' The normalization factors are always computed from \code{object}.
#' However, if \code{se.out} is a (different) SummarizedExperiment object, these factors are stored in \code{se.out} and the modified object.
#' This is useful when \code{se.out} contains counts for windows, but the normalization factors are computed using larger bins in \code{object}.
#' The same logic applies when \code{se.out} is a (different) DGEList object.
#' 
#' Note that an error is raised if the library sizes in \code{se.out} are not identical to \code{object$totals}.
#' This is because the normalization factors are only comparable when the library sizes are the same.
#' Consistent library sizes can be achieved by using the same \code{\link{readParam}} object in \code{\link{windowCounts}} and related functions.
#' 
#' @return
#' If \code{se.out=FALSE}, a numeric vector containing the relative normalization factors for each library.
#'
#' If \code{se.out=TRUE}, the same vector is stored in the \code{norm.factors} field of \code{mcols(object)} (if \code{object} is a SummarizedExperiment)
#' or \code{object$samples} (if \code{object} is a DGEList) and the modified \code{object} is returned.
#' 
#' If \code{se.out} is a separate SummarizedExperiment or DGEList object, 
#' the normalization factors are stored inside \code{se.out} and the modified object is returned.
#' 
#' @author Aaron Lun
#' 
#' @references
#' Robinson MD, Oshlack A (2010). 
#' A scaling normalization method for differential expression analysis of RNA-seq data. 
#' \emph{Genome Biology} 11, R25.
#' 
#' @examples
#' counts <- matrix(rnbinom(400, mu=10, size=20), ncol=4)
#' data <- SummarizedExperiment(list(counts=counts))
#' data$totals <- colSums(counts)
#' 
#' # TMM normalization.
#' normFactors(data)
#' 
#' @seealso
#' \code{\link{calcNormFactors}}, for the base method.
#'
#' \code{\link{normOffsets}}, for the trended normalization strategy.
#' 
#' @keywords normalization
#' @export
#' @importFrom edgeR calcNormFactors
normFactors <- function(object, method=NULL, weighted=FALSE, ..., assay.id="counts", se.out=TRUE) {
    if (!is.null(method)) {
        .Deprecated("'method=' in 'normFactors' is deprecated and ignored.")
    }

    info <- .fetch_norm_details(object, assay.id=assay.id)
    x <- info$x
    lib.sizes <- info$lib.sizes

	out <- calcNormFactors(x, lib.size=lib.sizes, doWeighting=weighted, method=method, ...)

	# Choosing to put these values in a different location, if requested. 
    if (!is.logical(se.out)) { 
        other <- .fetch_norm_details(se.out, assay.id=assay.id, msg="se.out")
        if (!identical(other$lib.sizes, lib.sizes)) {
            stop("library sizes of 'se.out' and 'object' are not identical")
        }
        object <- se.out
        se.out <- TRUE
    }

    if (!se.out) {
        out
    } else {
        if (is(object, "DGEList")) {
            object$samples$norm.factors <- out
        } else if (is(object, "SummarizedExperiment")) {
            object$norm.factors <- out
        } else {
            attr(object, "norm.factors") <- out
        }
        object
    }
}

#' @importFrom edgeR DGEList
#' @importFrom SummarizedExperiment assay
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
#' @importFrom Matrix colSums
.fetch_norm_details <- function(object, assay.id, msg="object") {
    if (is(object, "DGEList")) {
        x <- object$counts
        lib.sizes <- object$samples$lib.size
    } else if (is(object, "SummarizedExperiment")) {
        x <- assay(object, i=assay.id, withDimnames=FALSE)
        lib.sizes <- object$totals
    } else {
        x <- object
        lib.sizes <- colSums(x)
    }
    if (is.null(lib.sizes)) {
        stop(sprintf("library sizes not present in '%s'", msg))
    }
    list(x=x, lib.sizes=lib.sizes)
}
