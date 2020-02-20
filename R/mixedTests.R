#' Tests for mixed DB clusters
#' 
#' Intersects two one-sided tests to determine if a cluster contains tests with changes in both directions.
#' 
#' @inheritParams combineTests
#' @param fc.col An integer or string specifying the single column of \code{tab} containing the log-fold change.
#' @param ... Further arguments to pass to \code{\link{mixedTests}}.
#' 
#' @return
#' A \linkS4class{DataFrame} with one row per cluster and various fields:
#' \itemize{
#' \item An integer field \code{num.tests}, specifying the total number of tests in each cluster.
#' \item Two integer fields \code{num.up.*} and \code{num.down.*} for each log-FC column in \code{tab}, containing the number of tests with log-FCs significantly greater or less than 0, respectively.
#' See \code{?"\link{cluster-direction}"} for more details.
#' \item A numeric field containing the cluster-level p-value. 
#' If \code{pval.col=NULL}, this column is named \code{PValue}, otherwise its name is set to \code{colnames(tab[,pval.col])}.
#' \item A numeric field \code{FDR}, containing the BH-adjusted cluster-level p-value.
#' \item A character field \code{direction}, set to \code{"mixed"} for all clusters. 
#' See \code{?"\link{cluster-direction}"} for more details.
#' \item Two integer fields \code{rep.up.test} and \code{rep.down.test}, containing the row index (for \code{tab}) of representative tests with positive and negative sign, respectively, for each cluster.
#' See \code{?"\link{cluster-direction}"} for more details.
#' \item One numeric field \code{rep.up.*} and \code{rep.down.*} for each log-FC column in \code{tab}, containing log-fold changes for the representative tests in the cluster.
#' See \code{?"\link{cluster-direction}"} for more details.
#' }
#' Each row is named according to the ID of the corresponding cluster.
#'
#' @details
#' This function converts two-sided p-values to one-sided counterparts for each direction of log-fold change.
#' For each direction, the corresponding one-sided p-values are combined by \code{\link{combineTests}} to yield a combined p-value for each cluster.
#' Each cluster is associated with two combined p-values (one in each direction), which are intersected using the Berger's intersection-union test (IUT).
#' 
#' The IUT p-value provides evidence against the null hypothesis that either direction is not significant.
#' In short, a low p-value is only possible if there are significant changes in both directions.
#' This formally identifies genomic regions containing complex DB events, i.e., where depletion in one subinterval of the bound/enriched region is accompanied by increasing binding in another subinterval. 
#' Examples include swaps in adjacent TF binding locations between conditions or shifts in histone mark patterns in bidirectional promoters.
#' 
#' We expect that the p-values in \code{pval.col} are two-sided and independent of the sign of the log-fold change under the null hypothesis.
#' This is true for likelihood ratio tests but may not be true for others (e.g., from \code{\link{glmTreat}}), so caution is required when supplying values in \code{tab}.
#' 
#' @seealso
#' \code{\link{combineTests}}, for a more general-purpose method of combining tests.
#' 
#' @author Aaron Lun
#' 
#' @examples
#' ids <- round(runif(100, 1, 10))
#' tab <- data.frame(logFC=rnorm(100), logCPM=rnorm(100), PValue=rbeta(100, 1, 2))
#' mixed <- mixedTests(ids, tab)
#' head(mixed)
#' 
#' @references
#' Berger RL and Hsu JC (1996). 
#' Bioequivalence trials, intersection-union tests and equivalence confidence sets.
#' \emph{Statist. Sci.} 11, 283-319.
#' 
#' @export
#' @importFrom stats p.adjust
mixedTests <- function(ids, tab, weights=NULL, pval.col=NULL, fc.col=NULL, fc.threshold=0.05) {
    fc.col <- .parseFCcol(fc.col, tab, multiple=FALSE)
    pval.col <- .getPValCol(pval.col, tab)

    com.out <- .get_two_one_sided_results(ids, tab, pval.col=pval.col, fc.col=fc.col,
        weights=weights, fc.threshold=fc.threshold)
    up.com <- com.out$up
    down.com <- com.out$down

    # Taking the IUT p-value.
    pval.colname <- colnames(tab)[pval.col]
    up.com[,pval.colname] <- pmax(up.com[,pval.colname], down.com[,pval.colname])
    up.com$FDR <- p.adjust(up.com[,pval.colname], method="BH")
    up.com$direction <- rep("mixed", nrow(up.com))

    # Replacing the down counts.
    all.down <- sprintf("num.down.%s", colnames(tab)[fc.col])
    up.com[,all.down] <- down.com[,all.down]

    # Adding representative up/down statistics.
    rep.up <- grep("^rep\\.", colnames(up.com))
    colnames(up.com)[rep.up] <- sub("^rep\\.", "rep.up.", colnames(up.com)[rep.up])
    rep.down <- grep("^rep\\.", colnames(down.com))
    colnames(down.com)[rep.down] <- sub("^rep\\.", "rep.down.", colnames(down.com)[rep.down])

    cbind(up.com, down.com[,rep.down])
}

.get_two_one_sided_results <- function(ids, tab, pval.col, fc.col, ...) {
    all.p <- .make_one_sided(tab, pval.col=pval.col, fc.col=fc.col)

    # Combining the one-sided p-values.
    up.tab <- tab
    up.tab[,pval.col] <- all.p$up
    up.com <- combineTests(ids, up.tab, pval.col=pval.col, fc.col=fc.col, ...)

    # Repeating in the other direction.
    down.tab <- tab
    down.tab[,pval.col] <- all.p$down
    down.com <- combineTests(ids, down.tab, pval.col=pval.col, fc.col=fc.col, ...)

    list(up=up.com, down=down.com)
}

#' @export
#' @rdname mixedTests
mixedClusters <- function(...) {
    .Deprecated(new="mixedTests")
    mixedTests(...)
}
