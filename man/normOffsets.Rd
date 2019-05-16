\name{normOffsets}
\alias{normOffsets}
\alias{normFactors}

\title{Normalize counts across libraries}

\description{Calculate normalization factors or offsets using count data from multiple libraries.}

\usage{
normOffsets(object, ..., assay.id="counts", se.out=TRUE)

normFactors(object, method="TMM", weighted=FALSE, ..., 
    assay.id="counts", se.out=TRUE)
}

\arguments{
  \item{object}{A SummarizedExperiment object containing a count matrix.}
  \item{method}{String specifying the type of scaling normalization method to use.}
  \item{weighted}{A logical scalar indicating whether precision weights should be used for TMM normalization.}
  \item{...}{Other arguments to be passed to \code{\link{calcNormFactors}} for \code{type="scaling"}, or \code{\link{loessFit}} for \code{type="loess"}.}
  \item{assay.id}{An integer scalar or string specifying the assay values to use for normalization.}
  \item{se.out}{A logical scalar indicating whether or not a SummarizedExperiment object should be returned.
Alternatively, a SummarizedExperiment object in which normalization factors are to be stored.}
}

\details{
The \code{normFactors} function provides a convenience wrapper for the \code{\link{calcNormFactors}} function in the \pkg{edgeR} package.
This uses the trimmed mean of M-values (TMM) method to remove composition biases, typically in background regions of the genome.
Precision weighting is turned off by default so as to avoid upweighting high-abundance regions. 
These are more likely to be bound and thus more likely to be differentially bound. 
Assigning excessive weight to such regions will defeat the purpose of trimming when normalizing the coverage of background regions.

% Large changes in the normalization factor estimates with changes in the prior suggest that the counts are too low i.e. not
% enough new information in the dataset. This can be overcome by (obviously) increasing the counts. For example, binning
% can be performed with a larger bin size in \code{windowCounts} to obtain proportionally larger counts.

The \code{normOffsets} function performs non-linear normalization similar to the fast loess algorithm in \code{\link{normalizeCyclicLoess}}. 
This aims to account for mean dependencies in the efficiency biases between libraries.
For each sample, a lowess curve is fitted to the log-counts against the log-average count. 
The fitted value for each genomic window is used as an offset in a generalized linear model for that feature and sample. 
The use of the average count provides more stability than the average log-count when low counts are present for differentially bound regions.

Both functions expect SummarizedExperiment objects as input.
The count matrix to be used for normalization will be extracted according to the specified \code{assay.id} field.
Library sizes are extracted from \code{object$totals}.
}

\section{Different SummarizedExperiment outputs}{
For \code{normFactors}, the normalization factors are always computed from \code{object}.
However, if \code{se.out} is a (different) SummarizedExperiment object, these factors are stored and returned in \code{se.out}.
This is useful when \code{se.out} contains counts for windows, but the normalization factors are computed using larger bins in \code{object}.

For \code{normOffsets}, the trend fits are always computed from \code{object}.
However, if \code{se.out} is a (different) SummarizedExperiment object, the trend fits will be used to compute offsets for each entry in \code{se.out} using spline interpolation.
This is useful when \code{se.out} contains counts for windows in an endogenous genome, but the trend fits are computed using spike-in chromatin regions.

In both functions, an error is raised if the library sizes in \code{se.out$totals} are not identical to \code{object$totals}.
This is because the normalization factors (for \code{normFactors}) and average abundances (for \code{normOffsets}) are only comparable when the library sizes are the same.
Consistent library sizes can be achieved by using the same \code{\link{readParam}} object in \code{\link{windowCounts}} and related functions.
}

\value{
For \code{normFactors}, a numeric vector containing the relative normalization factors for each library is computed.
This is returned directly if \code{se.out=FALSE}, otherwise it is stored in the \code{norm.factors} field of the \code{mcols} of the output object.

For \code{normOffsets} with \code{type="loess"}, a numeric matrix of the same dimensions as \code{counts} is computed, containing the log-based offsets for use in GLM fitting.
This is returned directly if \code{se.out=FALSE}, otherwise it is stored in the \code{offsets} assay of the output object.

If \code{se.out=TRUE}, a SummarizedExperiment is returned that contains the computed normalization factors/offsets but is otherwise identical to \code{object}.
If \code{se.out} is a SummarizedExperiment object, the normalization factors and offsets will be stored in an object that is otherwise identical to \code{se.out}.
}

\author{Aaron Lun}

\references{
Robinson MD, Oshlack A (2010). 
A scaling normalization method for differential expression analysis of RNA-seq data. 
\emph{Genome Biology} 11, R25.

Ballman KV, Grill DE, Oberg AL, Therneau TM (2004). 
Faster cyclic loess: normalizing RNA arrays via linear models. 
\emph{Bioinformatics} 20, 2778-86.
}

\examples{
counts <- matrix(rnbinom(400, mu=10, size=20), ncol=4)
data <- SummarizedExperiment(list(counts=counts))
data$totals <- colSums(counts)

# TMM normalization.
normFactors(data)

# Using loess-based normalization, instead.
offsets <- normOffsets(data)
head(offsets)
offsets <- normOffsets(data, span=0.4)
offsets <- normOffsets(data, iterations=1)
}

\seealso{
\code{\link{calcNormFactors}},
\code{\link{loessFit}},
\code{\link{normalizeCyclicLoess}}
}

\keyword{normalization}
