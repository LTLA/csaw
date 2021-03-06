\name{regionCounts}
\alias{regionCounts}

\title{Count reads overlapping each region}

\description{Count the number of extended reads overlapping pre-specified regions}

\usage{
regionCounts(bam.files, regions, ext=100, param=readParam(),
    BPPARAM=SerialParam())
}

\arguments{
\item{bam.files}{A character vector containing paths to sorted and indexed BAM files.
Alternatively, a list of \linkS4class{BamFile} objects.}
\item{regions}{A \linkS4class{GenomicRanges} object containing the regions over which reads are to be counted.}
\item{ext}{An integer scalar or list describing the average length of the sequenced fragment in each library, see \code{?\link{windowCounts}}.}
\item{param}{A \linkS4class{readParam} object containing read extraction parameters, or a list of such objects (one for each BAM file).}
\item{BPPARAM}{A \linkS4class{BiocParallelParam} specifying how parallelization is to be performed across files.}
}

\value{
A \linkS4class{RangedSummarizedExperiment} object is returned containing one integer matrix.
Each entry of the matrix contains the count for each library (column) at each region (row).
The coordinates of each region are stored as the \code{rowRanges}.
The total number of reads, read length and extension length used in each library are stored in the \code{colData}.
Other parameters (e.g., \code{param}) are stored in the \code{metadata}.
}

\details{
This function simply provides a wrapper around \code{\link{countOverlaps}} for read counting into specified regions. 
It is provided so as to allow for counting with awareness of the other parameters, e.g., \code{ext}, \code{pe}.
This allows users to coordinate region-based counts with those from \code{\link{windowCounts}}. 
Checking that the output \code{totals} are the same between the two calls is strongly recommended.

Note that the strandedness of \code{regions} will not be considered when computing overlaps.
In other words, both forward and reverse reads will be counted into each region, regardless of the strandedness of that region.
This can be altered by setting the \code{forward} slot in the \code{param} object to only count reads from one strand or the other.
The strandedness of the output \code{rowRanges} will depend on the strand(s) from which reads were counted.

See \code{\link{windowCounts}} for more details on read extension.
}

\author{Aaron Lun}

\seealso{
\code{\link{countOverlaps}},
\code{\link{windowCounts}},
\code{\link{readParam}}
}

\examples{
# A low filter is only used here as the examples have very few reads.
bamFiles <- system.file("exdata", c("rep1.bam", "rep2.bam"), package="csaw")
incoming <- GRanges(c('chrA', 'chrA', 'chrB', 'chrC'), 
    IRanges(c(1, 500, 100, 1000), c(200, 1000, 700, 1500)))
regionCounts(bamFiles, regions=incoming)
regionCounts(bamFiles, regions=incoming, param=readParam(restrict="chrB"))

# Loading PE data.
bamFile <- system.file("exdata", "pet.bam", package="csaw")
regionCounts(bamFile, regions=incoming, param=readParam(pe="both"))
regionCounts(bamFile, regions=incoming, param=readParam(max.frag=100, 
	pe="first", restrict="chrA"))
}

\keyword{counting}
