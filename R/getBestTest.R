#' Get the best test in a cluster
#' 
#' Find the test with the greatest significance or the highest abundance in each cluster.
#' 
#' @inheritParams combineTests
#' @param by.pval Logical scalar indicating whether the best test should be selected on the basis of the smallest p-value.
#' If \code{FALSE}, the best test is defined as that with the highest abundance.
#' @param cpm.col An integer scalar or string specifying the column of \code{tab} containing the log-CPM values.
#' Defaults to \code{"logCPM"}.
#' 
#' @inherit combineTests return
#'
#' @details
#' Each cluster is defined as a set of tests with the same value of \code{ids} (any \code{NA} values are ignored).
#' If \code{by.pval=TRUE}, this function identifies the test with the lowest p-value as that with the strongest evidence against the null in each cluster.
#' The p-value of the chosen test is adjusted using the (Holm-)Bonferroni correction, based on the total number of tests in the parent cluster. 
#' This is necessary to obtain strong control of the family-wise error rate such that the best test can be taken from each cluster for further consideration.
#' 
#' % i.e. The configuration, in this case, is taking the best test.
#' 
#' The importance of each window in each cluster can be adjusted by supplying different relative \code{weights} values. 
#' Each weight is interpreted as a different threshold for each test in the cluster using the weighted Holm procedure. 
#' Larger weights correspond to lower thresholds, i.e., less evidence is needed to reject the null for tests deemed to be more important. 
#' This may be useful for upweighting particular tests such as those for windows containing a motif for the TF of interest.
#' 
#' Note the difference between this function and \code{\link{combineTests}}. 
#' The latter presents evidence for any rejections within a cluster. 
#' This function specifies the exact location of the rejection in the cluster, which may be more useful in some cases but at the cost of conservativeness. 
#' In both cases, clustering procedures such as \code{\link{mergeWindows}} can be used to identify the cluster.
#' 
#' % The vagueness of combineTests may be good enough in most applications (i.e. wanting to get a location
#' % to look at the genomic context, or in instances where differential binding is obvious). If error control
#' % at specific locations is needed, then getBestTests is probably more appropriate..
#' 
#' If \code{by.pval=FALSE}, the best test is defined as that with the highest log-CPM value. 
#' This should be independent of the p-value so no adjustment is necessary. Weights are not applied here. 
#' This mode may be useful when abundance is correlated to rejection under the alternative hypothesis, e.g., picking high-abundance regions that are more likely to contain peaks.
#' 
#' To obtain \code{ids}, a simple clustering approach for genomic windows is implemented in \code{\link{mergeWindows}}.
#' However, anything can be used so long as it is independent of the p-values and does not compromise type I error control, e.g., promoters, gene bodies, independently called peaks. 
#' Any tests with \code{NA} values for \code{ids} will be ignored.
#'
#' @examples
#' ids <- round(runif(100, 1, 10))
#' tab <- data.frame(logFC=rnorm(100), logCPM=rnorm(100), PValue=rbeta(100, 1, 2))
#' best <- getBestTest(ids, tab)
#' head(best)
#' 
#' best <- getBestTest(ids, tab, cpm.col="logCPM", pval.col="PValue")
#' head(best)
#' 
#' # With window weighting.
#' w <- round(runif(100, 1, 5))
#' best <- getBestTest(ids, tab, weight=w)
#' head(best)
#' 
#' # By logCPM.
#' best <- getBestTest(ids, tab, by.pval=FALSE)
#' head(best)
#' 
#' best <- getBestTest(ids, tab, by.pval=FALSE, cpm.col=2, pval.col=3)
#' head(best)
#' 
#' @seealso
#' \code{\link{combineTests}} and \code{\link{minimalTests}}, for other methods for obtaining cluster-level p-values.
#'
#' \code{\link{mergeWindows}}, to generate \code{ids}.
#'
#' \code{\link{glmQLFTest}}, for one method of generating \code{tab}.
#' 
#' @author Aaron Lun
#' 
#' @keywords testing
#'
#' @export
#' @importFrom S4Vectors splitAsList
getBestTest <- function(ids, tab, by.pval=TRUE, weights=NULL, pval.col=NULL, fc.col=NULL, fc.threshold=0.05, cpm.col=NULL) {
    if (!by.pval) {
        if (is.null(cpm.col)) { 
            cpm.col <- "logCPM" 
        }
        if (length(cpm.col)!=1L) {
            stop("absent or multiple logCPM columns are not supported")
        }

        cur.col <- tab[,cpm.col]
        by.id <- splitAsList(seq_along(ids), ids)
        best <- unlist(lapply(by.id, function(x) x[which.max(cur.col[x])]))
        is.best <- seq_along(ids) %in% best

        is.pval <- .getPValCol(pval.col, tab)
        tab[!is.best,is.pval] <- NA_real_
    }

    minimalTests(ids, tab, weights=weights, pval.col=pval.col, fc.col=fc.col, 
        fc.threshold=fc.threshold, min.sig.n=0, min.sig.prop=0)
}

# You can reduce the conservativeness of the Bonferroni method by using windows that 
# don't overlap. You can also just use combineTests with smaller clusters if you're getting
# lost within each cluster.
#    More complex approaches would attempt to estimate the correlation based on genome-wide
# data. I think the best approach would be to compute the 1st order statistic for p-values
# in each cluster of a given size across the dataset. You can then fit the observed statistics
# robustly (i.e., get rid of very low p-values) to a B(1, n) distribution to determine the
# effective number of independent tests 'n'. You can then apply that to the Bonferroni 
# correction for a cluster of that size.
