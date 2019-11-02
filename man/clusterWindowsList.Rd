\name{clusterWindowsList}
\alias{clusterWindowsList}

\title{Consolidate DB clusters}
\description{Consolidate DB results from multiple analyses with cluster-level FDR control.}

\usage{
clusterWindowsList(ranges.list, tab.list, equiweight=TRUE, ...)
}

\arguments{
\item{ranges.list}{A list of \linkS4class{GRanges} or \linkS4class{RangedSummarizedExperiment} objects,
usually containing windows of varying sizes from separate calls to \code{\link{windowCounts}}.}
\item{tab.list}{A list of data.frames of differential binding results, 
usually from separate analyses at differing window sizes.
Each should contain one row per interval for the corresponding entry of \code{ranges.list}.}
\item{equiweight}{a logical scalar indicating whether equal weighting from each analysis should be enforced}
\item{...}{arguments to be passed to \code{\link{clusterWindows}}}
}

\value{
A named list is returned containing:
\describe{
\item{\code{ranges}:}{A \linkS4class{GRanges} object containing the concatenated intervals from all elements of \code{x}.
The element-wise metadata of this object contains the integer field \code{origin}, 
an integer field specifying the index of \code{x} from which each interval was obtained.}
\item{\code{ids}:}{An integer vectors indicating the cluster ID for each window in \code{ranges}.
Non-significant windows that were not assigned to a cluster have IDs of \code{NA}, as described in \code{?\link{clusterWindows}}.}
\item{\code{regions}:}{A \linkS4class{GRanges} object containing the coordinates for each cluster.}
\item{\code{FDR}:}{A numeric scalar containing the cluster-level FDR estimate.}
}
}

\details{
This function consolidates DB results from multiple analyses, typically involving different window sizes.
The aim is to provide comprehensive detection of DB at a range of spatial resolutions.
Significant windows from each analysis are identified and used for clustering with \code{\link{clusterWindows}}.
This represents the post-hoc counterpart to \code{\link{mergeResultsList}}.

Some effort is required to equalize the contribution of the results from each analysis.
This is done by setting \code{equiweight=TRUE}, 
where the weight of each window is inversely proportional to the number of windows from that analysis.
These weights are used as frequency weights for window-level FDR control (to identify DB windows prior to clustering).
Otherwise, the final results would be dominated by large number of small windows.

Users can cluster by the sign of log-fold changes, to obtain clusters of DB windows of the same sign.
However, note that nested windows with opposite signs may give unintuitive results - see \code{\link{mergeWindows}} for details.
}

\seealso{
\code{\link{clusterWindows}}, the equivalent function for a single \linkS4class{GRanges} input.

\code{\link{mergeResultsList}}, for a more rigorous approach to clustering windows.
}

\author{
Aaron Lun
}

\examples{
# Making up some GRanges.
low <- GRanges("chrA", IRanges(runif(100, 1, 1000), width=5))
med <- GRanges("chrA", IRanges(runif(40, 1, 1000), width=10))
high <- GRanges("chrA", IRanges(runif(10, 1, 1000), width=20))

# Making up some DB results.
dbl <- data.frame(logFC=rnorm(length(low)), PValue=rbeta(length(low), 1, 20))
dbm <- data.frame(logFC=rnorm(length(med)), PValue=rbeta(length(med), 1, 20))
dbh <- data.frame(logFC=rnorm(length(high)), PValue=rbeta(length(high), 1, 20))
result.list <- list(dbl, dbm, dbh)

# Consolidating.
cons <- clusterWindowsList(list(low, med, high), result.list, tol=20)
cons$region
cons$id
cons$FDR

# Without weights.
cons <- clusterWindowsList(list(low, med, high), result.list, tol=20, 
    equiweight=FALSE)
cons$FDR

# Using the signs.
cons <- clusterWindowsList(list(low, med, high), result.list, tol=20, 
    fc.col="logFC")
cons$FDR
}

\keyword{clustering}
