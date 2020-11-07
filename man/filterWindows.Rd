\name{filterWindows}
\alias{filterWindowsProportion}
\alias{filterWindowsGlobal}
\alias{filterWindowsLocal}
\alias{filterWindowsControl}
\alias{scaleControlFilter}

\title{Filtering methods for RangedSummarizedExperiment objects}
\description{Convenience functions to compute filter statistics for windows, 
based on proportions or using enrichment over background.}

\usage{
filterWindowsGlobal(data, background, assay.data="counts", 
    assay.back="counts", prior.count=2, grid.pts=21)

filterWindowsLocal(data, background, assay.data="counts", 
    assay.back="counts", prior.count=2)

filterWindowsControl(data, background, assay.data="counts", 
    assay.back="counts", prior.count=2, scale.info=NULL)

filterWindowsProportion(data, assay.data="counts", prior.count=2)

scaleControlFilter(data.bin, back.bin, assay.data="counts", 
    assay.back="counts")
}

\arguments{
\item{data}{A \linkS4class{RangedSummarizedExperiment} object containing window-level counts.}
\item{background}{A \linkS4class{RangedSummarizedExperiment} object to be used for estimating background enrichment.
\itemize{
    \item For \code{filterWindowsGlobal}, this should contain counts for large contiguous bins across the genome,
    for the same samples used to construct \code{data}.
    \item For \code{filterWindowsLocal}, this should contain counts for regions in which \code{rowRanges(data)} is nested,
    for the same samples used to construct \code{data}.
    \item For \code{filterWindowsControl}, this should contain count for the same regions as \code{rowRanges(data)}
    for control samples.
}}
\item{assay.data}{A string or integer scalar specifying the assay containing window/bin counts in \code{data}.}
\item{assay.back}{A string or integer scalar specifying the assay containing window/bin counts in \code{background}.}
\item{prior.count}{A numeric scalar, specifying the prior count to use in \code{\link{scaledAverage}}.}
\item{scale.info}{A list containing the output of \code{scaleControlFilter}, 
i.e., a normalization factor and library sizes for ChIP and control samples.}
\item{data.bin}{A RangedSummarizedExperiment containing bin-level counts for ChIP libraries.}
\item{back.bin}{A RangedSummarizedExperiment containing bin-level counts for control libraries.}
\item{grid.pts}{An integer scalar specifying the number of grid points to use for interpolation when \code{data} contains variable-width intervals,
e.g., for peaks or other regions instead of windows.}
}

\value{
All filtering functions return a named list containing:
\itemize{
    \item \code{abundances}, a numeric vector containing the average abundance of each row in \code{data}.
    \item \code{filter}, a numeric vector containing the filter statistic for the given \code{type} for each row.
        The definition of this filter statistic will vary across the different methods.
    \item \code{back.abundances}, a numeric vector containing the background abundance for each row in \code{data}.
        Not reported for \code{filterWindowsProportion}.
}

For \code{scaleControlFilter}, a named list is returned containing:
\itemize{
    \item \code{scale}, a numeric scalar containing the scaling factor for multiplying the control counts.
    \item \code{data.totals}, a numeric vector containing the library sizes for \code{data}.
    \item \code{back.totals}, anumeric vector containing the library sizes for \code{background}.
}
}

\details{
The aim of these functions is to compute a filter statistic for each window,
according to various abundance-based definitions that are discussed below.
Windows can then be filtered to retain those with large filter statistics.
This selects for high-abundance windows that are more likely to contain genuine binding sites
and thus are more interesting for downstream (differential binding) analyses.
}

\section{Proportion-based filtering}{
\code{filterWindowsProportion} supposes that a certain percentage of the genome (by length) is genuinely bound.
The filter statistic is defined as the ratio of the rank in abundance to the total number of windows.
Rank is in ascending order, i.e., higher abundance windows have higher ratios.
Windows are retained that have rank ratios above a threshold, e.g., 0.99 if 1\% of the genome is assumed to be bound.

The definition of the rank is dependent on the total number of windows in the genome.
However, empty windows or bins are automatically discarded in \code{\link{windowCounts}} 
(exacerbated if \code{filter} is set above unity).
This will result in underestimation of the rank or overestimation of the global background.
To avoid this, the total number of windows is inferred from the spacing.
}

\section{Global background filtering}{
\code{filterWindowsGlobal} uses the median average abundance across the genome as a global estimate of the background abundance.
This assumes that \code{background} contains unfiltered counts for large (2 - 10 kbp) genomic bins, 
from which the background abundance can be stably computed.
The filter statistic for each window is defined as the difference between the window abundance and the global background,
after adjusting for the differences in widths between windows and bins.

Similarly to \code{filterWindowsProportion}, the total number of bins is inferred from their width in \code{background}.
This avoids overestimating the global background if some filtering has already been applied to the bins.
% For example, windowCounts will automatically remove bins with no counts (e.g., telomeres, centromeres).

The calculation is fairly straightforward for window-level counts where all \code{\link{rowRanges}(data)} are of the same width. 
If the widths are variable (e.g., because \code{data} contains counts for peaks, genes or other irregular features), 
the adjustment for differences in width needs to performed separately for each unique width.
If there are more than \code{grid.pts} unique widths, we expedite this process by computing the adjustment for \code{grid.pts} widths
and interpolating to obtain the adjusted background at each width.

If \code{background} is not supplied, the background abundance is directly computed from entries in \code{data}.
This assumes that \code{data} contains windows for most of the regions in the genome,
and that the coverage is sufficiently high across most windows to obtain a stable background estimate.
}

\section{Local background filtering}{
\code{filterWindowsLocal} compares the abundance of each window to the flanking regions.
This considers each window to contain the entirety of a binding event, 
where any coverage of the surrounding regions is treated as background enrichment.
It is analogous to the behaviour of peak-calling methods and accounts for local fluctuations in background,
e.g., due to differences in mappability, sequenceability or accessibility.

We assume that each region in \code{data} is nested within each corresponding region of \code{background}.
The counts of each row in \code{data} are then subtracted from those of the corresponding row in \code{background}.
The average abundance of the remaining counts is computed and used as the background abundance.
The filter statistic is defined by subtracting the background abundance from the corresponding window abundance for each row.

To generate \code{background}, we suggest using \code{\link{regionCounts}} on a \code{\link{resize}}d \code{rowRanges(data)} -
see Examples.
}

\section{Control-based filtering}{ 
In \code{filterWindowsControl}, we assume that \code{data} contains window-level counts for ChIP samples,
while \code{background} contains counts for the same windows in the control samples.
(However, it is also possible to supply nested regions as described in \code{filterWindowsLocal},
where each interval in \code{background} includes the flanking regions around the corresponding entry in \code{data}.)
For each window in \code{data}, 
the background abundance is defined as the average abundance of the corresponding row in \code{background}.
The filter statistic is defined as the difference between the window's average abundance and its background abundance.

Composition biases are likely to be present between ChIP and control samples, 
where increased binding at some loci reduces coverage of other loci in the ChIP samples.
This incorrectly results in smaller filter statistics for the latter loci, as the fold-change over the input is reduced.
To correct for this, a normalization factor between ChIP and control samples can be computed with \code{scaleControlFilter}
and passed to \code{filterWindowsControl} using the \code{scale.info} argument.
A warning is raised if \code{scale.info=NULL}. 

To use \code{scaleControlFilter}, users should supply two \linkS4class{RangedSummarizedExperiment} objects, 
each containing the counts for large (~10 kbp) bins in the ChIP and control samples.
The difference in the average abundance between the two objects is computed for each bin.
The median of the differences across all bins is used as a normalization factor to correct the filter statistics for each window.
The assumption is that most bins represent background regions, 
such that a systematic difference in abundance between ChIP and control should represent the composition bias.

\code{scaleControlFilter} will also store the library sizes for each object in its output.
This is used to check that \code{data} and \code{background} have the same library sizes.
Otherwise, the normalization factor computed with bin-level counts cannot be correctly applied to the window-level counts.
}

\section{Details on filter statistic calculations}{
When computing the filter statistic in background-based methods
the abundances of bins/regions in \code{background} must be rescaled for comparison to those of smaller windows 
- see \code{\link{getWidths}} and \code{\link{scaledAverage}} for more details.
In particular, the effective width of the window is often larger than the \code{width} used in \code{\link{windowCounts}}
due to the counting of fragments rather than reads.
The fragment length is extracted from \code{data$ext} and \code{background$ext}, 
though users will need to set \code{data$rlen} or \code{background$rlen} for unextended reads (i.e., \code{ext=NA}).

The \code{prior.count} protects against inflated log-fold increases when the background counts are near zero.
A low prior is sufficient if \code{background} has large counts, which is usually the case for wide regions.
Otherwise, if the set of windows with large filter statistics are dominated by low counts,
\code{prior.count} should be increased to a larger value like 5.
}

\seealso{
\code{\link{windowCounts}},
\code{\link{aveLogCPM}},
\code{\link{getWidths}},
\code{\link{scaledAverage}}
}

\examples{
bamFiles <- system.file("exdata", c("rep1.bam", "rep2.bam"), package="csaw")
data <- windowCounts(bamFiles, filter=1)

# Proportion-based (keeping top 1%)
stats <- filterWindowsProportion(data)
head(stats$filter)
keep <- stats$filter > 0.99 
new.data <- data[keep,]

# Global background-based (keeping fold-change above 3).
background <- windowCounts(bamFiles, bin=TRUE, width=300)
stats <- filterWindowsGlobal(data, background)
head(stats$filter)
keep <- stats$filter > log2(3)

# Local background-based.
locality <- regionCounts(bamFiles, resize(rowRanges(data), fix="center", 300))
stats <- filterWindowsLocal(data, locality)
head(stats$filter)
keep <- stats$filter > log2(3)

# Control-based, with binning for normalization (pretend rep2.bam is a control).
binned <- windowCounts(bamFiles, width=10000, bin=TRUE)
chip.bin <- binned[,1]
con.bin <- binned[,2]
scinfo <- scaleControlFilter(chip.bin, con.bin)

chip.data <- data[,1]
con.data <- data[,2]
stats <- filterWindowsControl(chip.data, con.data,
    prior.count=5, scale.info=scinfo)

head(stats$filter)
keep <- stats$filter > log2(3)
}

\keyword{filtering}
