\name{normOffsets}
\alias{normOffsets}
\alias{normOffsets,matrix-method}
\alias{normOffsets,SummarizedExperiment-method}

\alias{normalize}
\alias{normalize,SummarizedExperiment-method}

\title{Normalize counts between libraries}

\description{Calculate normalization factors or offsets using count data from multiple libraries.}

\usage{
\S4method{normOffsets}{matrix}(object, lib.sizes=NULL, type=c("scaling", "loess"), 
    weighted=FALSE, ...)

\S4method{normOffsets}{SummarizedExperiment}(object, assay=1, type="scaling", ..., se.out=NULL)

\S4method{normalize}{SummarizedExperiment}(object, ...)
}

\arguments{
  \item{object}{A matrix of numeric counts with one column per library and one row per genomic feature (e.g., window).}
  \item{lib.sizes}{A numeric vector specifying the total number of reads per library.
This is extracted from \code{object$totals} for SummarizedExperiment inputs, or computed as \code{colSums(object)} for matrices.}
  \item{type}{A character string indicating what type of normalization is to be performed.}
  \item{weighted}{A logical scalar indicating whether precision weights should be used for TMM normalization.}
  \item{...}{Other arguments to be passed to \code{\link{calcNormFactors}} for \code{type="scaling"}, or \code{\link{loessFit}} for \code{type="loess"}.}
  \item{assay}{An integer scalar or string specifying the assay values to use for normalization.}
  \item{se.out}{A logical scalar indicating whether or not a SummarizedExperiment object should be returned.
Alternatively, a SummarizedExperiment object in which normalization factors are to be stored.}
}

\details{
If \code{type="scaling"}, this function provides a convenience wrapper for the \code{\link{calcNormFactors}} function in the \pkg{edgeR} package.
Specifically, it uses the trimmed mean of M-values (TMM) method to perform normalization. 
Precision weighting is turned off by default so as to avoid upweighting high-abundance regions. 
These are more likely to be bound and thus more likely to be differentially bound. 
Assigning excessive weight to such regions will defeat the purpose of trimming when normalizing the coverage of background regions.

% We can't min/maximize some quantity with regards to the estimated normalization factors, and then choose a reference
% from that, as that would bias the resulting estimates (or at least, interact with existing bias unpredictably).

% Large changes in the normalization factor estimates with changes in the prior suggest that the counts are too low i.e. not
% enough new information in the dataset. This can be overcome by (obviously) increasing the counts. For example, binning
% can be performed with a larger bin size in \code{windowCounts} to obtain proportionally larger counts.

If \code{type="loess"}, this function performs non-linear normalization similar to the fast loess algorithm in \code{\link{normalizeCyclicLoess}}. 
For each sample, a lowess curve is fitted to the log-counts against the log-average count. 
The fitted value for each genomic window is used as an offset in a generalized linear model for that feature and sample. 
The use of the average count provides more stability than the average log-count when low counts are present for differentially bound regions.

If a SummarizedExperiment object is supplied, the values to be used for normalization will be extracted according to the specified \code{assay} field.
Furthermore, \code{se.out=TRUE} means that the SummarizedExperiment method will return a modified version of \code{object} containing normalization information.
Normalization factors are stored in the \code{"norm.factors"} field in the \code{mcols}, while the offset matrix is stored in the \code{"offset"} field in the \code{assays}.
\code{se.out=NULL} will soon be deprecated and defaults to \code{FALSE}.

\code{normalize,SummarizedExperiment-method} has the same behaviour as \code{normOffsets} with \code{se.out=TRUE}.
However, this will be deprecated to avoid clashes with other \code{normalize} methods for SummarizedExperiment objects.
}

\section{Additional details for SummarizedExperiment inputs}{
If \code{se.out} is a SummarizedExperiment object and \code{type="scaling"}, 
the function will calculate the normalization factors from \code{object} but return them in a modified version of \code{se.out}.
This is useful when \code{se.out} contains counts for windows, but the normalization factors are computed using larger bins in \code{object}.

Note that the normalization factors can only be interpreted with respect to the library sizes used to calculate them.
As such, the function will throw an error if the library sizes in \code{se.out$totals} are not identical to \code{object$totals}.
Consistent library sizes can be achieved by using the same \code{\link{readParam}} object in \code{\link{windowCounts}} and related functions.

The use of a SummarizedExperiment object in \code{se.out} with \code{type="loess"} is not yet supported.
}

\value{
For \code{type="scaling"}, a numeric vector containing the relative normalization factors for each library is returned.

For \code{type="loess"}, a numeric matrix of the same dimensions as \code{counts}, containing the log-based offsets for use in GLM fitting.

If \code{se.out=TRUE}, a SummarizedExperiment is returned that contains the computed normalization factors/offsets but is otherwise identical to \code{object}.
}

\author{Aaron Lun}

\references{
Robinson MD, Oshlack A (2010). A scaling normalization method for differential
expression analysis of RNA-seq data. \emph{Genome Biology} 11, R25.

Ballman KV, Grill DE, Oberg AL, Therneau TM (2004). Faster cyclic loess:
normalizing RNA arrays via linear models. \emph{Bioinformatics} 20, 2778-86.
}

\examples{
# A trivial example
counts <- matrix(rnbinom(400, mu=10, size=20), ncol=4)
normOffsets(counts)
normOffsets(counts, lib.sizes=rep(400, 4))

# Using loess-based normalization, instead.
offsets <- normOffsets(counts, type="loess", lib.size=full.lib.size)
head(offsets)
offsets <- normOffsets(counts, type="loess", span=0.4, lib.size=full.lib.size)
offsets <- normOffsets(counts, type="loess", iterations=1, lib.size=full.lib.size)

# Same for SummariedExperiment objects. 
bamFiles <- system.file("exdata", c("rep1.bam", "rep2.bam"), package="csaw")
data <- windowCounts(bamFiles, width=100, filter=1)

normOffsets(data, se.out=FALSE)
normOffsets(data, se.out=TRUE)

another.data <- windowCounts(bamFiles, width=10)
normOffsets(data, se.out=another.data)

normOffsets(data, type="loess", se.out=TRUE)
head(normOffsets(data, type="loess", se.out=FALSE))
}

\seealso{
\code{\link{calcNormFactors}},
\code{\link{loessFit}},
\code{\link{normalizeCyclicLoess}}
}

\keyword{normalization}
