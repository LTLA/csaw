\name{getPESizes}
\alias{getPESizes}

\title{Compute fragment lengths for paired-end tags}
 
\description{Compute the length of the sequenced fragment for each read pair in paired-end tag (PE) data.}

\usage{
getPESizes(bam.file, param=readParam(pe="both"))
}

\arguments{
	\item{bam.file}{A character string containing the file path to a sorted and indexed BAM file.}
    \item{param}{A \linkS4class{readParam} object containing read extraction parameters.}
}

\value{
A list containing:
\describe{
\item{\code{sizes}:}{An integer vector of fragment lengths for all valid read pairs in the library.}
\item{\code{diagnostics}:}{An integer vector containing the number of inter-chromosomal read pairs; 
    number of improperly orientated read pairs; and number of read pairs in discarded regions.}
}
}

\details{
This function computes a number of diagnostics for paired-end data in a supplied BAM file.
The aims is to provide an indication of the quality of library preparation and sequencing.

The output \code{sizes} will only contain fragment sizes for all correctly oriented pairs on the same chromosome.
Users are referred to \code{?\link{readParam}} for what \pkg{csaw} considers to be a proper pair.
Statistics are reported for the number of reads that fail prior to this point, for various reasons, e.g.,
inter-chromosomal, incorrectly oriented or discarded in blacklist regions.

Note that filtering by \code{max.frag} is not performed here.
The point of this function is to allow users to choose a threshold.
}

\author{Aaron Lun}

\seealso{
\code{\link{readParam}}
}

\examples{
bamFile <- system.file("exdata", "pet.bam", package="csaw")
out <- getPESizes(bamFile, param=readParam(pe="both"))
out <- getPESizes(bamFile, param=readParam(pe="both", restrict="chrA"))
out <- getPESizes(bamFile, param=readParam(pe="both", discard=GRanges("chrA", IRanges(1, 50))))
}

\keyword{diagnostics}
