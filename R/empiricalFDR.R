#' Control the empirical FDR
#' 
#' Control the empirical FDR across clusters for comparisons to negative controls, 
#' based on tests that are significant in the \dQuote{wrong} direction.
#' 
#' @inheritParams mixedTests
#' @param neg.down A logical scalar indicating if negative log-fold changes correspond to the \dQuote{wrong} direction.
#' 
#' @details
#' Some experiments involve comparisons to negative controls where there should be no signal/binding.
#' In such case, genuine differences should only occur in one direction, i.e., up in the non-control samples.
#' Thus, the number of significant tests that change in the wrong direction can be used as an estimate of the number of false positives.
#' 
#' This function converts two-sided p-values in \code{tab[,pval.col]} to one-sided counterparts in the wrong direction.
#' It combines the one-sided p-values for each cluster using \code{\link{combineTests}}.
#' The number of significant clusters at some p-value threshold represents the estimated number of false positive clusters.
#' 
#' The same approach is applied for one-sided p-values in the right direction, where the number of detected clusters at the threshold represents the total number of discoveries.
#' Dividing the number of false positives by the number of discoveries yields the empirical FDR at each p-value threshold.
#' Monotonicity is enforced (i.e., the empirical FDR only decreases with decreasing p-value) as is the fact that the empirical FDR must be below unity.
#' 
#' The p-values specified in \code{pval.col} are assumed to be originally computed from some two-sided test,
#' where the distribution of p-values is the same regardless of the direction of the log-fold change (under both the null and alternative hypothesis).
#' This rules out p-values computed from ANODEV where multiple contrasts are tested at once;
#' or from methods that yield asymmetric p-value distributions, e.g., GLM-based TREAT.
#' 
#' @section Caution:
#' Control of the empirical FDR is best used for very noisy data sets where the BH method is not adequate.
#' The BH method only protects against statistical false positives under the null hypothesis that the log-fold change is zero.
#' However, the empirical FDR also protects against experimental false positives, caused by non-specific binding that yields uninteresting (but statistically significant) DB.
#' 
#' The downside is that the empirical FDR calculation relies on the availability of a good estimate of the number of false positives.
#' It also assumes that the distribution of p-values is the same for non-specific binding events in both directions 
#' (i.e., known events with negative log-FCs and unknown events among those with positive log-FCs).
#' Even if the log-fold changes are symmetric around zero, this does not mean that the p-value distributions will be the same,
#' due to differences in library size and number between control and ChIP samples.
#' 
#' In summary, the BH method in \code{\link{combineTests}} is more statistically rigorous and should be preferred for routine analyses.
#' 
#' @return
#' A \linkS4class{DataFrame} with one row per cluster and various fields:
#' \itemize{
#' \item An integer field \code{num.tests}, specifying the total number of tests in each cluster.
#' \item Two integer fields \code{num.up.*} and \code{num.down.*} for each log-FC column in \code{tab}, containing the number of tests with log-FCs significantly greater or less than 0, respectively.
#' See \code{?"\link{cluster-direction}"} for more details.
#' \item A numeric field containing the cluster-level p-value. 
#' If \code{pval.col=NULL}, this column is named \code{PValue}, otherwise its name is set to \code{colnames(tab[,pval.col])}.
#' \item A numeric field \code{FDR}, containing the empirical FDR corresponding to that cluster's p-value.
#' \item A character field \code{direction} (if \code{fc.col} is of length 1), specifying the dominant direction of change for tests in each cluster.
#' See \code{?"\link{cluster-direction}"} for more details.
#' \item One integer field \code{rep.test} containing the row index (for \code{tab}) of a representative test for each cluster.
#' See \code{?"\link{cluster-direction}"} for more details.
#' \item One numeric field \code{rep.*} for each log-FC column in \code{tab}, containing a representative log-fold change for the differential tests in the cluster.
#' See \code{?"\link{cluster-direction}"} for more details.
#' }
#' Each row is named according to the ID of the corresponding cluster.
#' 
#' @author
#' Aaron Lun
#' 
#' @seealso
#' \code{\link{combineTests}}, used to combine the p-values in each direction.
#' 
#' @references
#' Zhang Y et al. (2008). 
#' Model-based Analysis of ChIP-Seq (MACS). 
#' \emph{Genome Biol.} 9, R137.
#'
#' @examples
#' ids <- round(runif(100, 1, 10))
#' tab <- data.frame(logFC=rnorm(100), logCPM=rnorm(100), PValue=rbeta(100, 1, 2))
#' empirical <- empiricalFDR(ids, tab)
#' head(empirical)
#'
#' @export
empiricalFDR <- function(ids, tab, weights=NULL, pval.col=NULL, fc.col=NULL, fc.threshold=0.05, neg.down=TRUE) {
    fc.col <- .parseFCcol(fc.col, tab, multiple=FALSE)
    pval.col <- .getPValCol(pval.col, tab)
    com.out <- .get_two_one_sided_results(ids, tab, pval.col=pval.col, fc.col=fc.col,
        weights=weights, fc.threshold=fc.threshold)

    if (neg.down) { 
        right <- "up"
        wrong <- "down"
    } else {
        right <- "down"
        wrong <- "up"
    }

    right.com <- com.out[[right]]
    wrong.com <- com.out[[wrong]]
    all.down <- sprintf("num.%s.%s", wrong, colnames(tab)[fc.col])
    right.com[,all.down] <- wrong.com[,all.down]

    # Computing the empirical FDR.
    pval.colname <- colnames(tab)[pval.col]
    right.comp <- right.com[,pval.colname]
    o <- order(right.comp)
    right.comp <- right.comp[o]
    empirical <- findInterval(right.comp, sort(wrong.com[,pval.colname]))/seq_along(right.comp)
    
    empirical <- pmin(1, empirical)
    empirical <- rev(cummin(rev(empirical)))
    empirical[o] <- empirical
    right.com$FDR <- empirical

    # Mopping up.
    right.com$direction <- rep(right, nrow(right.com))
    right.com
}

.make_one_sided <- function(tab, pval.col, fc.col) {
    cur.fc <- tab[,fc.col]
    going.up <- cur.fc > 0
    pval <- tab[,pval.col]
    
    # Calculating each set fresh, to avoid numeric 
    # imprecision from repeated "1-" operations
    up.p <- pval/2
    up.p[!going.up] <- 1 - up.p[!going.up]
    down.p <- pval/2
    down.p[going.up] <- 1 - down.p[going.up]

    list(up=up.p, down=down.p)
}
