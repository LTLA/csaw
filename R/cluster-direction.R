#' Reporting cluster-level direction in \pkg{csaw}
#'
#' An overview of the strategies used to obtain cluster-level summaries of the direction of change,
#' based on the directionality information of individual tests.
#' This is relevant to all functions that aggregate per-test statistics into a per-cluster summary,
#' e.g., \code{\link{combineTests}}, \code{\link{minimalTests}}.
#' It assumes that there are zero, one or many columns of log-fold changes in the data.frame of per-test statistics,
#' typically specified using a \code{fc.cols} argument.
#'
#' @section Counting the per-test directions:
#' For each cluster, we will report the number of tests that are up (positive log-fold change) or down (negative) for each column in \code{fc.col}.
#' This provide some indication of whether binding increases or decreases - or both - across tests in the cluster.
#' If a cluster contains non-negligble numbers of both up and down tests, this indicates that there may be a complex difference within that cluster. 
#' Similarly, complex differences may be present if the total number of tests is larger than the number of tests in either category (i.e., change is not consistent across the cluster).
#'
#' To count up/down tests, we apply a multiple testing correction to the p-values \emph{within} each cluster.
#' Only the tests with adjusted p-values no greater than \code{fc.threshold} are counted as being up or down.
#' We can interpret this as a test of conditional significance; assuming that the cluster is interesting (i.e., contains at least one true positive), what is the distribution of the signs of the changes within that cluster?
#' Note that this procedure has no bearing on the p-value reported for the cluster itself.
#'
#' The nature of the per-test correction within each cluster depends on the p-value-combining procedure implemented in each function.
#' In most cases, there is a per-test correction that naturally accompanies the per-cluster p-value:
#' \itemize{
#' \item For \code{\link{combineTests}}, the Benjamini-Hochberg correction is used.
#' \item For \code{\link{minimalTests}}, the Holm correction is used. 
#' \item For \code{\link{getBestTest}} with \code{by.pval=TRUE}, the Holm correction is used.
#' We could also use the Bonferroni correction here but Holm is universally more powerful so we use that instead.
#' \item For \code{\link{getBestTest}} with \code{by.pval=FALSE}, 
#' all tests bar the one with the highest abundance are simply removed.
#' This mimics the application of an independent filter.
#' No correction is applied as only one test remains.
#' \item For \code{\link{mixedTests}}, the Benjamini-Hochberg correction is used, 
#' given that this function just calls \code{\link{combineTests}} on the one-sided p-values in each direction.
#' Here, though, the number of up tests is obtained using the one-sided p-values for a positive change;
#' similarly, the number of down tests is obtained using the one-sided p-values for a negative change.
#' }
#' 
#' @section Determining the cluster-level direction:
#' When only one log-fold change column is specified, we will try to determine which direction contributes to the combined p-value.
#' This is done by considering whether the cluster-level p-value would increase if all tests in one direction were assigned p-values of unity.
#' If there is an increase, then tests changing in that direction must contribute to the combined p-value calculations. 
#' In this manner, clusters are labelled based on whether their combined p-values are driven by tests with only positive (\code{"up"}) or negative log-fold changes (\code{"down"}) or both (\code{"mixed"}).
#' 
#' The label for each cluster is stored as the \code{"direction"} field in the returned data frame.
#' However, keep in mind that the label only describes the direction of change among the most significant tests in the cluster.
#' Clusters with complex differences may still be labelled as changing in only one direction, if the tests changing in one direction have much lower p-values than the tests changing in the other direction (even if both sets of p-values are significant).
#' More rigorous checks for mixed changes should be performed with \code{\link{mixedTests}}.
#'
#' Speaking of which, the \code{"direction"} field for \code{\link{mixedTests}} is simply set to \code{"mixed"} for all clusters.
#' This reflects the fact that the reported p-value represents the evidence for mixed directionality in this function;
#' indeed, the field itself is simply reported for consistency, given that we already know we are looking for mixed clusters! 
#'
#' @section Representative log-fold changes:
#' For each combining procedure, it is usually possible to identify a representative test for the entire cluster.
#' The index of this test is reported in the output as the \code{"rep.test"} field along with its log-fold changes.
#' \itemize{
#' \item For \code{\link{combineTests}}, the test with the lowest BH-adjusted p-value before enforcing monotonicity is used.
#' \item For \code{\link{minimalTests}}, the test with the \eqn{x}th-smallest p-value is used (see the documentation for details).
#' \item For \code{\link{getBestTest}}, the test with the lowest p-value is used if \code{by.pval=TRUE}.
#' Otherwise, the test with the highest abundance is used.
#' \item For \code{\link{mixedTests}}, two representative tests are reported in each direction.
#' The representative test in each direction is defined using \code{\link{combineTests}} as described above.
#' }
#' Note that the chosen test is only representative in the sense that its p-value is especially important for computing the cluster-level p-value.
#' This is usually sufficient to provide proxy log-fold changes for clusters with simple differences,
#' but is obviously inadequate for representing more complex changes.
#' 
#' @author Aaron Lun
#'
#' @name cluster-direction
NULL

