#' @importFrom limma loessFit
#' @importFrom stats spline
#' @importFrom edgeR aveLogCPM
#' @importFrom methods is
#' @importFrom SummarizedExperiment assay assay<-
normOffsets <- function(object, lib.sizes=NULL, type=c("scaling", "loess"), ..., assay.id="counts", se.out=TRUE) 
# Perform a fast loess normalization which uses the average count as the covariate, 
# rather than the typical A-value-based methods to avoid instability at low abundances.
#
# written by Aaron Lun
# created 19 November 2013
{
	type <- match.arg(type)
	if (type=="scaling") { 
        .Deprecated(msg='type="scaling" is deprecated.\nUse normFactors() instead.')
		return(normFactors(object, lib.size=lib.size, ..., assay.id=assay.id, se.out=se.out))
	} 

    mat <- assay(object, i=assay.id, withDimnames=FALSE)
    if (is.null(lib.sizes)) {
        lib.sizes <- object$totals
        if (is.null(lib.sizes)) {
            stop("library sizes not present in 'object$totals'")
        }
    }
    nlibs <- length(lib.sizes)
    nwin <- nrow(mat)

    # Scaled corrections squeeze offsets towards relative log-library sizes.
    # Constant value of `cont.cor' would squeeze them towards zero, which could be misleading.
    cont.cor <- 0.5
    cont.cor.scaled <- cont.cor * lib.sizes/mean(lib.sizes)
		
    # Using it as a prior.count for abundance ensures linearity with log-counts
    ab <- aveLogCPM(mat, lib.size=lib.sizes, prior.count=cont.cor)/log2(exp(1))

    # Computing the output mean.
    mean.out <- NULL
    if (is(se.out, "SummarizedExperiment")) { 
        object <- se.out
        se.out <- TRUE
        
        nwin <- nrow(object)
        if (ncol(object)!=nlibs) {
            stop("number of libraries differs between 'se.out' and 'object'")
        }

        se.lib.sizes <- object$totals
        if (is.null(se.lib.sizes)) {
            stop("library sizes not present in 'se.out$totals'")
        }
        if (!identical(se.lib.sizes, lib.sizes)) {
            stop("library sizes should be identical between 'se.out' and 'object'")
        }

        mat2 <- assay(object, i=assay.id, withDimnames=FALSE)
        mean.out <- aveLogCPM(mat2, lib.size=lib.sizes, prior.count=cont.cor)/log2(exp(1))
    }
    
    # Computing the offsets as log-differences from the average (implicit if the same average is used across libraries).
    offs <- matrix(0, nwin, nlibs, byrow=TRUE)
    for (x in seq_len(nlibs)) {
    	fit <- loessFit(log(mat[,x]+cont.cor.scaled[x]), ab, ...)

        if (is.null(mean.out)) {
            offs[,x] <- fit$fitted 
        } else {
            offs[,x] <- spline(x=ab, y=fit$fitted, xout=mean.out)$y
        }
    }
    offs <- offs - rowMeans(offs)

    # Deciding what to return.
    if (!se.out) { 
        return(offs)
    } else {
        assay(object, "offset") <- offs
        return(object)
    }
}

