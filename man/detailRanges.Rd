\name{detailRanges}
\alias{detailRanges}

\title{Add annotation to ranges}
\description{Add detailed exon-based annotation to specified genomic regions.}

\usage{
detailRanges(incoming, txdb, orgdb, dist=5000, promoter=c(3000, 1000), 
    key.field="ENTREZID", name.field="SYMBOL", ignore.strand=TRUE)
}

\arguments{
\item{incoming}{A \linkS4class{GRanges} object containing the ranges to be annotated.}
\item{txdb}{A TxDb object for the genome of interest.}
\item{orgdb}{An OrgDb object for the genome of interest.}
\item{dist}{An integer scalar specifying the flanking distance to annotate.}
\item{promoter}{An integer vector of length 2, where first and second values define the promoter as some distance upstream and downstream from the TSS, respectively.}
\item{key.field}{A character scalar specifying the key type in \code{orgdb} corresponding to the gene IDs in \code{txdb}.}
\item{name.field}{A character scalar specifying the column from \code{orgdb} to use as the gene name.}
\item{ignore.strand}{A logical scalar indicating whether strandedness in \code{incoming} should be ignored.}
}

\value{
If \code{incoming} is not provided, a \linkS4class{GRanges} object will be returned containing ranges for the exons, promoters and gene bodies. 
Gene keys (e.g., Entrez IDs) are povided as row names.
Gene symbols and feature types are stored as metadata.

If \code{incoming} is a \linkS4class{GRanges} object, a list will be returned with \code{overlap}, \code{left} and \code{right} elements. 
Each element is a character vector of length equal to the number of ranges in \code{incoming}.
Each non-empty string records the gene symbol, the overlapped exons and the strand. 
For \code{left} and \code{right}, the gap between the range and the annotated feature is also included.
}

\details{
This function adds annotations to a given set of genomic regions in the form of compact character strings specifying the features overlapping and flanking each region.
The aim is to determine the genic context of empirically identified regions, for some basic biological interpretation of binding/marking in those regions. 
All neighboring genes within a specified range are reported, rather than just the closest gene to the region.
If a region in \code{incoming} is stranded and \code{ignore.strand=FALSE}, annotated features will only be reported if they lie on the same strand as that region.

If \code{incoming} is missing, then the annotation will be provided directly to the user in the form of a \linkS4class{GRanges} object. 
This may be more useful when further work on the annotation is required.
Features are labelled as exons (\code{"E"}), promoters (\code{"P"}) or gene bodies (\code{"G"}).
Overlaps to introns can be identified by finding those regions that overlap with gene bodies but not with any of the corresponding exons.

The default settings for \code{key.field} and \code{name.field} will work for human and mouse genomes, but may not work for other organisms.
The \code{key.field} should specify the key type in the \code{orgdb} object that corresponds to the gene IDs of the \code{txdb} object.
For example, in S. cerevisiae, \code{key.field} is set to \code{"ORF"} to match the gene IDs in the corresponding TxDb object,
while \code{name.field} is set to \code{"GENENAME"} to obtain the gene symbols.
}

\section{Explanation of fields}{
For annotated features overlapping a region, the character string in the \code{overlap} output vector will be of the form \code{GENE:STRAND:TYPE}. 
\code{GENE} is the gene symbol by default, but reverts to the key (default Entrez ID) if no symbol is defined.
\code{STRAND} is simply the strand of the gene, either \code{"+"} or \code{"-"}.
The \code{TYPE} indicates the feature types that are overlapped - exon (\code{"E"}), promoter (\code{"P"}) and/or intron (\code{"I"}).
Note that intron overlaps are only reported if the region does not overlap an exon directly.

For annotated features flanking the region within a distance of \code{dist}, the \code{TYPE} is instead the distance to the feature.
This represents the gap between the edge of the region and the closest exon for that gene.
Flanking promoters are not reported, as it is more informative to report the distance to the exon directly;
and flanking an intron should be impossible without overlapping an exon directly (and thus should not be reported, see above).
Note that exons directly overlapping the supplied region are not considered for flanking annotation, as the distance would be negative.

The strand information is often useful in conjunction with the left/right flanking features.
For example, if an exon for a negative-strand gene is to the left, the current region must be upstream of that exon.
Conversely, if the exon for a positive-strand gene is to the left, the region must be downstream.
The opposite applies for features to the right of the current region.
}

\author{Aaron Lun}

\examples{ 
library(org.Mm.eg.db)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)

current <- readRDS(system.file("exdata", "exrange.rds", package="csaw"))
output <- detailRanges(current, orgdb=org.Mm.eg.db,
    txdb=TxDb.Mmusculus.UCSC.mm10.knownGene) 
head(output$overlap)
head(output$right)
head(output$left)
}

\keyword{annotation}
