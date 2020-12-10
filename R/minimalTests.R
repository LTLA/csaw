#' Require rejection of a minimal number of tests
#'
#' Compute a p-value for each cluster based around the rejection of a minimal number or proportion of tests from that cluster.
#' 
#' @inheritParams combineTests
#' @param min.sig.n Integer scalar containing the minimum number of significant barcodes when \code{method="holm-min"}.
#' @param min.sig.prop Numeric scalar containing the minimum proportion of significant barcodes when \code{method="holm-min"}.
#'
#' @inherit combineTests return
#'
#' @details
#' All tests with the same value of \code{ids} are used to define a single cluster.
#' For each cluster, this function applies the Holm-Bonferroni correction to the p-values from all of its tests.
#' It then chooses the \eqn{x}th-smallest adjusted p-value as the cluster-level p-value,
#' where \eqn{x} is defined from the larger of \code{min.sig.n} and the product of \code{min.sig.prop} and the number of tests.
#' (If \eqn{x} is larger than the total number of tests, the largest per-test p-value is used instead.)
#'
#' % More formally, the cluster-level null hypothesis is that the per-test null hypothesis is false for fewer than \eqn{x} tests.
#' % We reject the cluster-level null at a nominal threshold \eqn{T} if we can reject \eqn{x} or more tests while controlling the FWER at \eqn{T}.
#' % This is done using the Holm-Bonferroni method to obtain strong FWER control over the rejected tests in the cluster.
#' % It is then straightforward to invert this logic to obtain a cluster-level p-value from the Holm-adjusted per-test p-values.
#'
#' Here, a cluster can only achieve a low p-value if at least \eqn{x} tests also have low p-values.
#' This favors clusters that exhibit consistent changes across all tests,
#' which is useful for detecting, e.g., systematic increases in binding across a broad genomic region spanning many windows.
#' By comparison, \code{\link{combineTests}} will detect a strong change in a small subinterval of a large region,
#' which may not be of interest in some circumstances.
#'
#' The importance of each test within a cluster can be adjusted by supplying different relative \code{weights} values. 
#' This may be useful for downweighting low-confidence tests, e.g., those in repeat regions. 
#' In the weighted Holm procedure, weights are used to downscale the per-test p-values,
#' effectively adjusting the distribution of per-test errors that contribute to family-wise errors.
#' Note that these weights have no effect between clusters.
#'
#' To obtain \code{ids}, a simple clustering approach for genomic windows is implemented in \code{\link{mergeWindows}}.
#' However, anything can be used so long as it is independent of the p-values and does not compromise type I error control, e.g., promoters, gene bodies, independently called peaks. 
#' Any tests with \code{NA} values for \code{ids} will be ignored.
#'
#' @author Aaron Lun
#'
#' @examples
#' ids <- round(runif(100, 1, 10))
#' tab <- data.frame(logFC=rnorm(100), logCPM=rnorm(100), PValue=rbeta(100, 1, 2))
#' minimal <- minimalTests(ids, tab)
#' head(minimal)
#'
#' @seealso
#' \code{\link{groupedHolmMin}}, which does the heavy lifting.
#'
#' \code{\link{combineTests}} and \code{\link{getBestTest}}, for another method of combining p-values for each cluster.
#'
#' \code{\link{mergeWindows}}, for one method of generating \code{ids}.
#' 
#' \code{\link{glmQLFTest}}, for one method of generating \code{tab}.
#' 
#' @references
#' Holm S (1979).
#' A simple sequentially rejective multiple test procedure.
#' \emph{Scand. J. Stat.} 6, 65-70.
#'
#' @export
#' @importFrom metapod groupedHolmMin
minimalTests <- function(ids, tab, min.sig.n=3, min.sig.prop=0.4, weights=NULL, pval.col=NULL, fc.col=NULL, fc.threshold=0.05) {
    .general_test_combiner(ids=ids, tab=tab, weights=weights, 
        pval.col=pval.col, fc.col=fc.col, fc.threshold=fc.threshold,
        FUN=function(...) groupedHolmMin(..., min.n=min.sig.n, min.prop=min.sig.prop),
        count.correct="holm"
    )
}
