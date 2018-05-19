\name{consolidateTests}
\alias{consolidateTests}
\alias{consolidateOverlaps}

\title{Consolidate window statistics}
\description{Consolidate DB statistics from analyses using multiple window sizes.}

\usage{
consolidateTests(id.list, result.list, weight.list, 
    FUN=combineTests, ...)

consolidateOverlaps(olap.list, result.list, weight.list, 
    FUN=combineOverlaps, ...)
}

\arguments{
\item{id.list}{A list of integer vectors specifying the identity of the cluster to which each window in each analysis belongs.
Typically produced by \code{\link{consolidateWindows}} with \code{region=NULL}.}
\item{olap.list}{A list of Hits objects specifying the overlaps between windows of each analysis and a common set of pre-defined regions.
Typically produced by \code{\link{consolidateWindows}} with \code{region!=NULL}.}
\item{result.list}{A list of data.frames containing the DB test results for each analysis.}
\item{weight.list}{A list of numeric vectors specifying the weight of each window in \code{id.list} or overlapping window in \code{olap.list}.}
\item{FUN}{A function specifying how the statistics should be consolidated.}
\item{...}{Further arguments to pass to \code{FUN}.}
}

\value{
A data.frame of results depending on \code{FUN}, see \code{\link{combineTests}} and \code{\link{getBestTest}} for examples.
Each row corresponds to a cluster for \code{consolidateTests} or to a pre-defined region for \code{consolidateOverlaps}.
}

\details{
By default, these functions will use \code{\link{combineTests}} or \code{\link{combineOverlaps}} to consolidate statistics.
This will yield a single combined p-value for each cluster or region.
Users can obtain the best windows in each cluster/region by setting \code{FUN=getBestTest} or \code{FUN=getBestOverlaps} instead.

If \code{weight.list} is not specified, it is set to \code{NULL} and a warning is raised.
Users should generally pass the weights produced by \code{\link{consolidateWindows}} with \code{equiweight=TRUE},
to ensure that all window sizes contribute equally to the final result.
}

\seealso{
\code{\link{consolidateWindows}},
\code{\link{combineTests}},
\code{\link{combineOverlaps}}
}

\author{
Aaron Lun
}

\examples{
bamFiles <- system.file("exdata", c("rep1.bam", "rep2.bam"), package="csaw")
low <- windowCounts(bamFiles, width=1, filter=1)
med <- windowCounts(bamFiles, width=100, filter=1)
high <- windowCounts(bamFiles, width=500, filter=1)

# Mocking up some pretend DB results.
dbl <- data.frame(logFC=rnorm(nrow(low)), PValue=runif(nrow(low)), logCPM=0)
dbm <- data.frame(logFC=rnorm(nrow(med)), PValue=runif(nrow(med)), logCPM=0)
dbh <- data.frame(logFC=rnorm(nrow(high)), PValue=runif(nrow(high)), logCPM=0)

# Consolidating.
cons <- consolidateWindows(list(low, med, high), 
	merge.args=list(tol=100, max.width=300))
comb <- consolidateTests(cons$id, result.list=list(dbl, dbm, dbh),
    weight.list=cons$weight)
head(comb)

# Trying with a custom region.
of.interest <- GRanges(c('chrA', 'chrA', 'chrB', 'chrC'), 
    IRanges(c(1, 500, 100, 1000), c(200, 1000, 700, 1500)))
cons <- consolidateWindows(list(low, med, high), region=of.interest)
comb <- consolidateOverlaps(cons$olap, result.list=list(dbl, dbm, dbh),
    weight.list=cons$weight)
head(comb)
}

\keyword{clustering}
