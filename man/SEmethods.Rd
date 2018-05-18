\name{SEmethods}
\alias{asDGEList}
\alias{asDGEList,SummarizedExperiment-method}

\title{SummarizedExperiment to DGEList}

\description{Converting a SummarizedExperiment object to a DGEList object for analysis with \pkg{edgeR}.}

\usage{
\S4method{asDGEList}{SummarizedExperiment}(object, lib.sizes, norm.factors, assay.id="counts", ...)
}

\arguments{
  \item{object}{A SummarizedExperiment object or its derived classes, like that produced by \code{\link{windowCounts}}.}
  \item{lib.sizes}{An (optional) integer vector of library sizes.}
  \item{norm.factors}{An (optional) numeric vector of normalization factors.}
  \item{assay.id}{A string or integer scalar indicating which assay in \code{object} contains the count matrix.}
  \item{...}{Other arguments to be passed to \code{\link{DGEList}}.}
}

\details{
Counts are extracted from specified \code{assay} matrix in the SummarizedExperiment object and used to construct a DGEList object via \code{\link{DGEList}}.
If not specified in \code{lib.sizes}, library sizes are taken from the \code{totals} field in the column data of \code{object}.
Warnings will be generated if this field is not present.

If \code{norm.factors} is not specified, \code{asDGEList} will attempt to extract normalization factors from \code{object$norm.factors}.
If this is not available, factors will be set to the default (all unity).
If \code{assays(object)$offset} is present, this will be assigned to the \code{offset} field of the output \code{DGEList} object.
}

\value{
A DGEList object is returned containing counts and normalization information.
}

\author{Aaron Lun}

\examples{
bamFiles <- system.file("exdata", c("rep1.bam", "rep2.bam"), package="csaw")
data <- windowCounts(bamFiles, width=100, filter=1)

asDGEList(data)
asDGEList(data, lib.sizes=c(10, 100))
asDGEList(data, norm.factors=c(1.11, 2.23), group=c("a", "b"))
}

\seealso{
\code{\link{DGEList}}
}

\keyword{normalization}
