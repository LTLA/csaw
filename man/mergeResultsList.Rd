\name{ResultList wrappers}
\alias{mergeResultsList}
\alias{overlapResultsList}

\title{Wrappers to combine multiple analyses}

\description{Convenience wrappers around lower-level functions to easily obtain region-level results 
from multiple sets of windows, each with their own statistics.}

\usage{
mergeResultsList(ranges.list, tab.list=NULL, tol, equiweight=TRUE, 
    get.best=TRUE, merge.args=list(), combine.args=list(), 
    best.args=list())

overlapResultsList(ranges.list, tab.list=NULL, regions, equiweight=TRUE, 
    get.best=TRUE, overlap.args=list(), combine.args=list(), 
    best.args=list())
}

\arguments{
\item{ranges.list}{A list of \linkS4class{GRanges} or \linkS4class{RangedSummarizedExperiment} objects,
usually containing windows of varying sizes from separate calls to \code{\link{windowCounts}}.}
\item{tab.list}{A list of data.frames of differential binding results, usually from separate analyses at differing window sizes.
Each should contain one row per interval for the corresponding entry of \code{ranges.list}.
If \code{NULL}, defaults to a list containing the \code{\link{mcols}} of each entry of \code{ranges.list}.}
\item{tol}{A numeric scalar to be passed to \code{\link{mergeWindows}}.}
\item{equiweight}{A logical scalar indicating whether windows should be weighted
so that each analysis contributes equally to the combined p-values.}
\item{get.best}{A logical scalar indicating whether \code{\link{getBestTest}} or \code{\link{getBestOverlaps}} should be run.}
\item{regions}{A \linkS4class{GRanges} class containing reference regions to overlap with the windows.}
\item{merge.args}{Further arguments to pass to \code{\link{mergeWindows}}.}
\item{overlap.args}{Further arguments to pass to \code{\link{findOverlaps}}.}
\item{combine.args}{Further arguments to pass to \code{\link{combineTests}} or \code{\link{combineOverlaps}}.}
\item{best.args}{Further arguments to pass to \code{\link{getBestTest}} or \code{\link{getBestOverlaps}}.}
}

\value{
Both functions will return a \linkS4class{DataFrame} containing
\code{regions}, a \linkS4class{GRanges} containing intervals of clustered windows (for \code{mergeResultsList}) 
or the input \code{regions} (for \code{overlapResultsList});
\code{combined}, the DataFrame containing combined statistics for each region;
and \code{best}, the DataFrame containing statistics for the best window in each region.

The \code{metadata} of the DataFrame will contain
\code{ranges}, a \linkS4class{GRanges} object containing the concatenated intervals from all elements of \code{ranges.list};
\code{tab}, a data.frame containing the combined results from all \code{tab.list};
and \code{weights}, a numeric vector of weights per window (for \code{mergeResultsList}) 
or per overlap (for \code{overlapsResultList}).

For \code{mergeResultsList}, the metadata will also contain \code{ids}, 
an integer vector specifying the cluster ID in the output \code{regions} for each interval in \code{ranges}.

For \code{overlapResultsList}, the metadata will also contain \code{overlaps},
a \linkS4class{Hits} object specifying the overlaps between the input \code{regions} (query) and \code{ranges}.
}

\details{
\code{mergeResultsList} calls \code{\link{mergeWindowsList}} to identify clusters of windows corresponding to regions.
These clusters are used in \code{\link{combineTests}} and (optionally) \code{\link{getBestTest}} to obtain per-region statistics.

\code{overlapResultsList} calls \code{\link{findOverlapsList}} to identify all windows overlapping each reference region.
These overlaps are used in \code{\link{combineOverlaps}} and (optionally) \code{\link{getBestOverlaps}},
where each set of overlapping windows for a reference region are used to define the latter's statistics.

Weighting is performed as described in \code{\link{mergeWindowsList}} and \code{\link{findOverlapsList}}.
If \code{equiweight=FALSE}, weights are reported but are not used in the p-value calculations.
}

\author{
Aaron Lun
}

\seealso{
\code{\link{mergeWindowsList}}, \code{\link{combineTests}} and \code{\link{getBestTest}}, 
for the functions that actually do the work in \code{mergeResultsList}.

\code{\link{findOverlapsList}}, \code{\link{combineOverlaps}} and \code{\link{getBestOverlaps}}, 
for the functions that actually do the work in \code{overlapResultsList}.
}

\examples{
# Making up some data.
tab1 <- data.frame(logFC=rnorm(20), logCPM=rnorm(20), PValue=rbeta(20, 1, 2))
gr1 <- GRanges("chrA", IRanges(sample(1000, 20), width=40))

tab2 <- data.frame(logFC=rnorm(20), logCPM=rnorm(20), PValue=rbeta(20, 1, 2))
gr2 <- GRanges("chrA", IRanges(sample(1000, 20), width=40))

# Merging windows.
(merge.out <- mergeResultsList(list(gr1, gr2), list(tab1, tab2), tol=1))
merge.out$regions
merge.out$combined

# Overlapping windows with a reference.
ref <- GRanges("chrA", IRanges(c(10, 500, 900), width=20))
(olap.out <- overlapResultsList(list(gr1, gr2), list(tab1, tab2), regions=ref))
olap.out$regions
olap.out$combined
}

