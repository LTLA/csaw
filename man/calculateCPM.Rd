\name{calculateCPM}
\alias{calculateCPM}

\title{Calculate CPM}
\description{Calculate counts-per-million (CPM) values for each feature.}

\usage{
calculateCPM(object, use.norm.factors=TRUE, use.offsets=FALSE, 
    log=TRUE, prior.count=1, assay.id="counts")
}

\arguments{
\item{object}{A SummarizedExperiment object containing a count matrix, as produced by \code{\link{windowCounts}}.}
\item{use.norm.factors}{A logical scalar indicating whether to use normalization factors, if available.}
\item{use.offsets}{A logical scalar indicating whether to use offsets, if available.}
\item{log}{A logical scalar indicating whether log2-transformed CPM values should be returned.}
\item{prior.count}{A numeric scalar specifying the prior count to add when \code{log=TRUE}.}
\item{assay.id}{A string or integer scalar indicating which assay of \code{y} contains the counts.}
}

\value{
A matrix of the same dimensions as \code{object}, containing (log-)transformed CPM values for each feature in each sample.
}

\details{
CPMs are calculated in the standard manner when \code{log=FALSE}, \code{use.offsets=FALSE} and \code{use.norm.factors=FALSE}.

When \code{log=TRUE}, a library size-adjusted prior count is added to both the counts and the library sizes, see \code{\link[edgeR]{cpm}} for details.

When \code{use.norm.factors=TRUE}, the effective library size is used for computing CPMs, provided that normalization factors are available in \code{object}.
This is defined as the product of the library size in \code{object$totals} and the normalization factor in \code{object$norm.factors}.

If \code{use.offsets=TRUE}, the offsets are converted into effective library sizes using \code{\link{scaleOffset}}.
If \code{log=TRUE}, this is done after addition of a prior count to both the counts and library sizes, see \code{\link{addPriorCount}} for details.
}

\author{
Aaron Lun
}

\seealso{
\code{\link[edgeR]{cpm}},
\code{\link{scaleOffset}},
\code{\link{addPriorCount}}
}

\examples{
bamFiles <- system.file("exdata", c("rep1.bam", "rep2.bam"), package="csaw")
data1 <- windowCounts(bamFiles, width=50, filter=1)
head(calculateCPM(data1))

data1$norm.factors <- c(0.5, 1.5)
head(calculateCPM(data1))

head(calculateCPM(data1, log=FALSE))

# With offsets:
assay(data1, "offset") <- matrix(rnorm(nrow(data1)*ncol(data1)), 
    nrow=nrow(data1), ncol=ncol(data1))
head(calculateCPM(data1, use.offsets=TRUE))

head(calculateCPM(data1, use.offsets=FALSE))
}
