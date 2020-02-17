#' Require rejection of a minimal number of tests
#'
#' Compute a p-value for each cluster based around the rejection of a minimal number or proportion of tests from that cluster.
#' 
#' @inheritParams combineTests
#' @param min.sig.n Integer scalar containing the minimum number of significant barcodes when \code{method="holm-min"}.
#' @param min.sig.prop Numeric scalar containing the minimum proportion of significant barcodes when \code{method="holm-min"}.
#'
#' @return A \linkS4class{DataFrame} containing:
#' \itemize{
#' \item An integer field \code{NumTests}, specifying the total number of tests in each cluster.
#' \item One or more numeric fields containing the log-fold changes of the test with the \eqn{x}th-smallest p-value in each cluster.
#' \item A numeric field containing the combined p-value. 
#' If \code{pval.col=NULL}, this column is named \code{PValue}, otherwise its name is set to \code{colnames(tab[,pval.col])}.
#' \item A numeric field \code{FDR}, containing the q-value corresponding to the combined p-value.
#' }
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
#' @author Aaron Lun
#'
#' @examples
#' ids <- round(runif(100, 1, 10))
#' tab <- data.frame(logFC=rnorm(100), logCPM=rnorm(100), PValue=rbeta(100, 1, 2))
#' minimal <- minimalTests(ids, tab)
#' head(minimal)
#'
#' @export
#' @importFrom stats p.adjust
#' @importFrom S4Vectors DataFrame
minimalTests <- function(ids, tab, min.sig.n=3, min.sig.prop=0.4, weight=NULL, pval.col=NULL, fc.col=NULL) {
    is.pval <- .getPValCol(pval.col, tab)
    fc.col <- .parseFCcol(fc.col, tab) 

    # Counting the number of tests.
    by.id <- split(seq_along(ids), ids)
    output <- DataFrame(NumTests=lengths(by.id), row.names=names(by.id))

    # Computing the Holm minimal p-value.
    barcode.p <- numeric(nrow(output))
    chosen <- integer(nrow(output))
    all.p <- tab[,is.pval]

    for (i in seq_along(by.id)) {
        current <- by.id[[i]]
        p <- p.adjust(all.p[current], method="holm")

        n <- min(max(min.sig.n, ceiling(min.sig.prop * length(p))), length(p))
        o <- order(p)
        idx <- o[n]

        barcode.p[i] <- p[idx]
        chosen[i] <- current[idx]
    }

    # Throwing in statistics for the representative window.
    output <- cbind(output, tab[chosen,fc.col])

    output[,is.pval] <- barcode.p
    output$FDR <- p.adjust(output[[is.pval]], method="BH")
    output
}

