#' Normalize trended biases across libraries
#' 
#' Calculate normalization offsets by performing a loess fit to count data from multiple libraries.
#' 
#' @inheritParams normFactors
#' @param ... Other arguments to be passed to \code{\link{loessFit}}. 
#' 
#' @details
#' This function performs non-linear normalization similar to the fast loess algorithm in \code{\link{normalizeCyclicLoess}}. 
#' The aims is to account for mean dependencies in the efficiency biases between libraries.
#' For each sample, a lowess curve is fitted to the log-counts against the log-average count. 
#' The fitted value for each genomic window is used as an offset in a generalized linear model for that feature and sample. 
#' The use of the average count provides more stability than the average log-count when low counts are present for differentially bound regions.
#' 
#' The trend fits are always computed from \code{object}.
#' However, if \code{se.out} is a (different) SummarizedExperiment or DGEList object, 
#' the trend fits will be used to compute offsets for each entry in \code{se.out} using spline interpolation.
#' This is useful when \code{se.out} contains counts for windows in an endogenous genome, but the trend fits are computed using spike-in chromatin regions.
#' 
#' An error is raised if the library sizes in \code{se.out$totals} are not identical to \code{object$totals}.
#' This is because the average abundances are only comparable when the library sizes are the same.
#' Consistent library sizes can be achieved by using the same \code{\link{readParam}} object in \code{\link{windowCounts}} and related functions.
#' 
#' @return
#' If \code{se.out=FALSE}, a numeric matrix of dimensions equal to \code{object}, containing the offset for each observation.
#' These offsets have already been scaled to be comparable in magnitude to the log-library sizes.
#'
#' If \code{se.out=TRUE}, the same matrix is stored in the \code{offset} assay of \code{object} (if it is a SummarizedExperiment)
#' or \code{object$offset} (if a DGEList) and the modified \code{object} is returned.
#' 
#' If \code{se.out} is a separate SummarizedExperiment or DGEList object, the offset matrix instead has dimensions equal to \code{se.out}.
#' This matrix is stored inside \code{se.out} and the modified object is returned.
#' 
#' @author Aaron Lun
#' 
#' @references
#' Ballman KV, Grill DE, Oberg AL, Therneau TM (2004). 
#' Faster cyclic loess: normalizing RNA arrays via linear models. 
#' \emph{Bioinformatics} 20, 2778-86.
#' 
#' @examples
#' counts <- matrix(rnbinom(400, mu=10, size=20), ncol=4)
#' data <- SummarizedExperiment(list(counts=counts))
#' data$totals <- colSums(counts)
#' 
#' # TMM normalization.
#' normFactors(data)
#' 
#' # Using loess-based normalization, instead.
#' offsets <- normOffsets(data)
#' head(offsets)
#' offsets <- normOffsets(data, span=0.4)
#' offsets <- normOffsets(data, iterations=1)
#' 
#' @seealso
#' \code{\link{loessFit}}, for the fitting algorithm.
#'
#' \code{\link{normalizeCyclicLoess}}, for the original inspiration for this method.
#' 
#' @keywords normalization
#' 
#' @export
#' @importFrom limma loessFit
#' @importFrom stats spline
#' @importFrom edgeR aveLogCPM scaleOffset
#' @importFrom methods is
#' @importFrom SummarizedExperiment assay assay<-
normOffsets <- function(object, ..., assay.id="counts", se.out=TRUE) {
    info <- .fetch_norm_details(object, assay.id=assay.id)
    mat <- info$x
    lib.sizes <- info$lib.sizes

    nlibs <- ncol(mat)
    nwin <- nrow(mat)

    # Scaled continuity corrections squeeze offsets towards relative log-library sizes.
    # Constant value of `cont.cor' would squeeze them towards zero, which could be misleading.
    cont.cor <- 0.5
    cont.cor.scaled <- cont.cor * lib.sizes/mean(lib.sizes)
		
    # Using it as a prior.count for abundance ensures linearity with log-counts
    ab <- aveLogCPM(mat, lib.size=lib.sizes, prior.count=cont.cor)/log2(exp(1))

    # Computing the output mean.
    mean.out <- NULL
    if (!is.logical(se.out)) {
        other <- .fetch_norm_details(se.out, assay.id=assay.id, msg="se.out")
        if (!identical(other$lib.sizes, lib.sizes)) {
            stop("library sizes of 'se.out' and 'object' are not identical")
        }

        mat2 <- other$x
        mean.out <- aveLogCPM(mat2, lib.size=lib.sizes, prior.count=cont.cor)/log2(exp(1))

        object <- se.out
        se.out <- TRUE
        nwin <- nrow(object)
    }
    
    # Computing the offsets as log-differences from the average (implicit if the same average is used across libraries).
    offs <- matrix(0, nwin, nlibs, byrow=TRUE, dimnames=dimnames(object))
    for (x in seq_len(nlibs)) {
    	fit <- loessFit(log(mat[,x]+cont.cor.scaled[x]), ab, ...)
        if (is.null(mean.out)) {
            offs[,x] <- fit$fitted 
        } else {
            offs[,x] <- spline(x=ab, y=fit$fitted, xout=mean.out)$y
        }
    }

    offs <- scaleOffset(lib.sizes, offs)

    # Deciding what to return.
    if (!se.out) { 
        offs
    } else {
        if (is(object, "DGEList")) {
            object$offset <- offs
        } else if (is(object, "SummarizedExperiment")) {
            assay(object, "offset") <- offs
        } else {
            attr(object, "offset") <- offs
        }
        object
    }
}

