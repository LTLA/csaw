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
#' For each cluster, this function applies the Holm-Bonferroni correction to all of its tests.
#' It then chooses the \eqn{x}th-smallest adjusted p-value as the cluster-level p-value,
#' where \eqn{x} is defined from the larger of \code{min.sig.n} and the product of \code{min.sig.prop} and the number of tests.
#' (If \eqn{x} is larger than the total number of tests, the largest per-test p-value is used instead.)
#'
#' More formally, the cluster-level null hypothesis is that the per-test null hypothesis is false for fewer than \eqn{x} tests.
#' We reject the cluster-level null at a nominal threshold \eqn{T} if we can reject \eqn{x} or more tests while controlling the FWER at \eqn{T}.
#' This is done using the Holm-Bonferroni method to obtain strong FWER control over the rejected tests in the cluster.
#' It is then straightforward to invert this logic to obtain a cluster-level p-value from the Holm-adjusted per-test p-values.
#'
#' The idea is that a cluster can only achieve a low p-value if at least \eqn{x} tests also have low p-values.
#' This favors clusters that exhibit consistent changes across all tests,
#' which is useful for detecting, e.g., systematic increases in binding across a broad genomic region spanning many windows.
#'
#' The importance of each test within a cluster can be adjusted by supplying different relative \code{weights} values. 
#' This may be useful for downweighting low-confidence tests, e.g., those in repeat regions. 
#' In Holm's procedure, weights are interpreted as scaling factors on the nominal threshold for each test to adjust the distribution of errors.
#' Note that these weights have no effect between clusters and will not be used to adjust the computed FDR.
#' 
#' @author Aaron Lun
#'
#' @examples
#' ids <- round(runif(100, 1, 10))
#' tab <- data.frame(logFC=rnorm(100), logCPM=rnorm(100), PValue=rbeta(100, 1, 2))
#' minimal <- minimalTests(ids, tab)
#' head(minimal)
#'
#' @references
#' Benjamini Y and Hochberg Y (1997). 
#' Multiple hypotheses testing with weights. 
#' \emph{Scand. J. Stat.} 24, 407-418.
#' 
#' @export
#' @importFrom stats p.adjust
#' @importFrom S4Vectors DataFrame
minimalTests <- function(ids, tab, min.sig.n=3, min.sig.prop=0.4, weights=NULL, pval.col=NULL, fc.col=NULL, fc.threshold=0.05) {
    .general_test_combiner(ids=ids, tab=tab, weights=weights, 
        pval.col=pval.col, fc.col=fc.col, fc.threshold=fc.threshold,
        FUN=function(...) {
            .Call(cxx_compute_cluster_holm, ..., min.sig.n, min.sig.prop)
        }
    )
}
