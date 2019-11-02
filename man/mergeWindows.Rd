\name{mergeWindows}
\alias{mergeWindows}

\title{Merge windows into clusters}

\description{Uses a simple single-linkage approach to merge adjacent or overlapping windows into clusters.}

\usage{
mergeWindows(ranges, tol, signs=NULL, max.width=NULL, ignore.strand=TRUE)
}

\arguments{
	\item{ranges}{A \linkS4class{GRanges} or \linkS4class{RangedSummarizedExperiment} object containing window coordinates.}
	\item{tol}{A numeric scalar specifying the maximum distance between adjacent windows.}
	\item{signs}{A logical vector specifying whether each window has a positive log-FC.}
	\item{max.width}{A numeric scalar specifying the maximum size of merged intervals.}
	\item{ignore.strand}{A logical scalar indicating whether to consider the strandedness of \code{ranges}.}
}

\value{
A list containing \code{ids}, an integer vector containing the cluster ID for each window; 
and \code{regions}, a GRanges object containing the start/stop coordinates for each cluster of windows.
}

\details{
Windows in \code{ranges} are merged if the gap between the end of one window and the start of the next is no greater than \code{tol}. 
Adjacent windows can then be chained together to build a cluster of windows across the linear genome. 
A value of zero for \code{tol} means that the windows must be contiguous whereas negative values specify minimum overlaps.

Specification of \code{max.width} prevents the formation of excessively large clusters when many adjacent regions are present. 
Any cluster that is wider than \code{max.width} is split into multiple subclusters of (roughly) equal size.
Specifically, the cluster interval is partitioned into the smallest number of equally-sized subintervals where each subinterval is smaller than \code{max.width}.
Windows are then assigned to each subinterval based on the location of the window midpoints.
Suggested values range from 2000 to 10000 bp, but no limits are placed on the maximum size if it is \code{NULL}. 

The tolerance should reflect the minimum distance at which two regions of enrichment are considered separate. 
If two windows are more than \code{tol} apart, they \emph{will} be placed into separate clusters. 
In contrast, the \code{max.width} value reflects the maximum distance at which two windows can be considered part of the same region. 

If \code{ignore.strand=FALSE}, the entries in \code{ranges} are split into their separate strands.
The function is run separately on the entries for each strand, and the results collated.
The \code{region} returned in the output will be stranded to reflect the strand of the contributing input regions.
This may be useful for strand-specific applications.

Note that, in the output, the cluster ID reported in \code{id} corresponds to the index of the cluster coordinates in \code{region}.
}

\section{Splitting clusters by sign}{
If \code{sign!=NULL}, windows are only merged if they have the same sign of the log-FC and are not separated by intervening windows with opposite log-FC values. 
This can occasionally be useful to ensure consistent changes when summarizing adjacent DB regions of opposing sign.
However, it is \emph{not} recommended for routine clustering in differential analyses as the resulting clusters will not be independent of the p-value.

To illustrate, consider any number of non-DB sites, some of which will have large log-fold change by chance.
Sites with large log-fold changes will be more likely to form large clusters when \code{signs} is specified, as the overlapping windows are more likely to be consistent in their sign.
In contrast, sites with log-fold changes close to zero are likely to form smaller clusters as the overlapping windows will oscillate around a log-fold change of zero.
At best, this results in conservativeness in the correction, as smaller p-values are grouped together while larger p-values form more (smaller) clusters.
At worst, this results in anticonservativeness if further filtering is applied to remove smaller clusters with few windows.

Also, if any nested regions are present with opposing sign, sign-aware clustering may become rather unintuitive.
Imagine a chain of overlapping windows with positive log-fold changes, and in a window in the middle of this chain, a single window with a negative log-fold change is nested.
The chain would ordinarily form a single cluster, but this is broken by the negative log-FC window.
Thus, two clusters form (before and after the negative window - three clusters, if one includes the negative window itself) despite complete overlaps between all clusters.
}

\examples{
x <- round(runif(10, 100, 1000))
gr <- GRanges(rep("chrA", 10), IRanges(x, x+40))
mergeWindows(gr, 1)
mergeWindows(gr, 10)
mergeWindows(gr, 100)
mergeWindows(gr, 100, sign=rep(c(TRUE, FALSE), 5))
}

\seealso{
\code{\link{combineTests}}, \code{\link{windowCounts}}
}

\author{Aaron Lun}

\keyword{clustering}
