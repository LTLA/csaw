\name{csawnews}
\title{csaw News}
\encoding{UTF-8}

\section{Version 1.26.0}{\itemize{
\item \code{filterWindowsGlobal()} finally behaves correctly for variable-width \code{data=}.

\item \code{normFactors()} and \code{normOffsets()} accept DGEList inputs in their \code{object=} and \code{se.out=} arguments.

\item \code{mergeResults()} and friends now default to taking \code{tab=} from the \code{mcols()} of the inputs.
}}

\section{Version 1.24.0}{\itemize{
\item Accept a list of BamFiles in \code{bam.files=} for all functions.
}}

\section{Version 1.22.0}{\itemize{
\item Removed deprecated functionality in consolidateTests() and related functions.
\item Modified the method used to count the number of up/down windows in combineTests().
Reported a representative test in the output along with its log-fold changes.
\item Renamed mixedClusters() to mixedTests() with deprecation of the former.
\item Added minimalTests() for cluster-level summarization with a minimum proportion of significant tests.
\item Standardized the output of getBestTest(), empiricalFDR(), mixedTests() and minimalTests() to that of combineTests().
In particular, getBestTest() now returns the index of its best window in the rep.test field.
}}

\section{Version 1.20.0}{\itemize{
\item Removed deprecated functionality in normOffsets(), readParam(), scaledAverage().

\item Added mergeResults(), overlapResults() wrapper functions to simplify getting region-level results.

\item Added mergeWindowsList(), findOverlapsList() functions for consolidating windows from multiple analyses.
Deprecated consolidateWindows().

\item Added mergeResultsList(), overlapResultsList() wrapper functions to obtain consolidated region-level results.
Deprecated consolidateTests(), consolidateOverlaps().

\item Renamed regions= to ranges= in mergeWindows() for consistency. 

\item Added clusterWindowsList() to replace consolidateClusters().

\item Added filterWindowsGlobal(), filterWindowsLocal(), filterWindowsProportion() and filterWindowsControl().
Deprecated filterWindows().

\item Minor renaming of scaleControlInfo() arguments, added assay.data= and assay.back= arguments.
}}

\section{Version 1.18.0}{\itemize{
\item Deprecated the type= argument in normOffsets().

\item Removed the deprecated consolidateSizes() function.

\item Deprecated the use of BPPARAM= in readParam().
Added the BPPARAM= argument to all relevant functions.

\item Modified csawUsersGuide() to point to the workflow package.
}}

\section{Version 1.16.0}{\itemize{
\item Added normFactors() function to avoid confusion when normOffsets() returns factors.

\item Deprecated type="scaling" option in normOffsets().

\item Added calculateCPM() function for convenient calculation of (log-)CPMs.

\item Split up consolidateSizes() function into consolidateWindows(), consolidateTests() and consolidateOverlaps().
Deprecated consolidateSizes() itself.

\item Switched output of combineTests() and getBestTest() and related functions to a DataFrame.

\item Modified mergeWindows() behaviour with specified sign=, for dealing with nested windows of opposing sign.

\item Altered controlClusterFDR() to take the largest adjusted p-value threshold that yields a cluster-level FDR below target=.

\item Simplified detailRanges() output so that it no longer returns an arbitrary exon number.
}}

\section{Version 1.14.0}{\itemize{
\item Fully removed support for paramList objects. 

\item Removed support for normalize(), modified default option for se.out= in normOffsets().
}}

\section{Version 1.12.0}{\itemize{
\item Removed support for paramList objects. 

\item Added option for normOffsets() to return SummarizedExperiment objects containing normalization data.

\item Deprecated normalize() to avoid S4 method clashes. 

\item Moved scaling prior to control-based filtering into a new function, scaleControlFilter(), for greater modularity.

\item Updated user's guide.
}}

\section{Version 1.10.0}{\itemize{
\item 
Added calculation of dominant directionality in combineTests().
Fixed out-of-array indexing bug in the C++ code. 

\item
Supported factor input for ids argument in combineTests(), getBestTest().

\item
Added the empiricalFDR(), empiricalOverlaps() functions for controlling the empirical FDR.

\item
Added the mixedClusters(), mixedOverlaps() functions for testing for mixed clusters.

\item 
Ensured that window-level FDR threshold chosen by controlClusterFDR() is not above the cluster-level FDR.

\item
Minor fix to scaledAverage() to avoid slightly inaccurate results.
Also, zero or negative scale factors now return -Inf and NA, respectively.

\item
Switched to new scaleOffset() function for adding offsets in asDGEList().
Added option to specify the assay to be used.

\item
Added multi-TSS support in detailRanges().

\item 
Modified paired-end machinery in windowCounts(), getPESizes() to be more accommodating of overruns.

\item
Ignored secondary and supplementary alignments in all functions.

\item
Added options to specify assay in SE objects in filterWindows().

\item
Replaced weighting with normalization options in profileSites().

\item 
Updated user's guide.
}}

\section{Version 1.8.0}{\itemize{
\item Added protection against NA values in filterWindows().

\item Deprecated the use of parameter lists in any param= arguments.

\item Tightened up allowable values of ext= arguments in various functions.

\item Added the BPPARAM slot in the readParam class to store BiocParallelParam objects.

\item Added support for parallelization in windowCounts(), regionCounts() and others.

\item Updated documentation, user's guide.
}}

\section{Version 1.6.0}{\itemize{
\item 
Restored normalize() as a S4 method returning a RangedSummarizedExperiment object.

\item 
Modified asDGEList() to use any available normalization data in the input object.

\item 
Generalized S4 methods to apply on SummarizedExperiment objects.

\item
Removed the rescue.ext option for PE handling, to maintain consistent totals calculations.

\item
Removed the fast.pe option for PE data handling, in favour of improved default processing.

\item
Removed dumpPE(), which is not required without the fast.pe option.

\item
Removed makeExtVector() in favour of list/DataFrame specification.

\item
windowCounts() and regionCounts() now compute and store the mean PE size and read length.

\item
Minor fix in correlateReads() for end-of-chromosome behaviour.

\item
Modified checkBimodality() so that the width argument behaves like ext in windowCounts().

\item
extractReads() with as.reads=TRUE for PE data now returns a GRangesList.

\item
Added the controlClusterFDR(), clusterWindows() and consolidateClusters() functions to automate control of the cluster-level FDR.

\item 
Added protection against NA values in the cluster IDs for combineTests(), getBestTest(), upweightSummits().

\item
All read extraction methods are now CIGAR-aware and will ignore soft-clipped parts of the alignment. 
}}

\section{Version 1.4.0}{\itemize{
\item
Added clusterFDR() function to compute the FDR for clusters of DB windows.

\item 
Added checkBimodality() function to compute bimodality scores for regions.

\item
Modified normalize(), asDGEList() to allow manual specification of library sizes.

\item
Switched from normalizeCounts(), normalize() to S4 method normOffsets().

\item 
Modified default parameter specification in strandedCounts(), to avoid errors.

\item
Switched to warning from error when a restricted chromosome is specified in extractReads().

\item 
Modified extractReads() interface with improved support for extended read and paired read extraction.

\item 
Added normalization options to filterWindows() when using control samples.

\item
Fixed bug for proportional filtering in filterWindows().

\item
Allowed correlateReads() to accept paired-end specification when extracting data.

\item
Added maximizeCcf() function to estimate the average fragment length.

\item 
Added support for strand-specific overlapping in detailRanges().

\item
Increased the fidelity of retained information in dumped BAM file from dumpPE().

\item
Modified strand specification arguments for profileSites(), allowed reporting of individual profiles.

\item
Removed param= specification from wwhm().

\item 
Switched to RangedSummarizedExperiment conventions for all relevant functions.

\item
Switched to mapqFilter for scanBam() when filtering on mapping quality.

\item
Added tests for previously untested functions.

\item
Slight updates to documentation, user's guide.
}}

\section{Version 1.2.0}{\itemize{
\item 
Added getWidths(), scaledAverage() and filterWindows(), to facilitate comparison of abundances during filtering.

\item 
Added findMaxima() to identify locally maximal windows from range data.

\item 
Added profileSites() to examine the coverage profile around specified regions, with wwhm() to guess the ideal window size from the profile.

\item
Changed default window width in windowCounts() to 50 bp, default filter to a fixed count of 10. Also, filter=0 is honored when bin=TRUE.

\item
Switched from the depracated rowData to rowRanges for all manipulations of SummarizedExperiment.

\item 
Changed all instances of `pet' to `pe' in read parameter specification, and renamed getPETSizes() to getPESizes().

\item 
Removed the redundant rescue.pairs parameter in readParam().

\item
Added fast.pe parameter in readParam(), for fast paired-end data extraction. Added dumpPE() to pre-process paired-end BAM files for fast downstream extraction.

\item 
Added support for custom column specification in getBestTest(), combineTests().

\item
Switched from reporting average log-FC to numbers of up/down windows in combineTests().

\item
Allowed getBestTest() to return all fields associated with the best window in the output table.

\item
Added upweightSummits() to compute weights favouring high-abundance windows.

\item 
Added combineOverlaps(), getBestOverlaps() and summitOverlaps() wrapper functions for processing of Hits objects.

\item 
Added consolidateSizes(), to consolidate DB results from multiple window sizes.

\item
Added support for custom key/name specification in detailRanges() for non-human/mouse systems.

\item 
Added support for strand-specific read extraction in readParam(), strand-specific counting via strandedCounts().

\item 
Added strand-awareness to mergeWindows(). Added protection against stranded input regions in extractReads(), detailRanges().

\item
Changed algorithm for splitting of large peaks in mergeWindows().

\item
Stored counting parameters in exptData for windowCounts(), regionCounts().

\item
Fixed small inaccuracies with continuity correction addition in normalizeCounts() for NB-loess.

\item 
Switched to fragment midpoint for binning of paired-end data in windowCounts().

\item 
Added support for lists of library-specific readParam objects in windowCounts(), regionCounts(), correlateReads().

\item 
Added makeExtVector(), to support variable read extension lengths between libraries in windowCounts(), regionCounts().

\item
Added support for read extension within extractReads().

\item 
Updated the user's guide to reflect new and modified functions.

\item 
Added sra2bam.sh in inst/doc to reproducibly generate BAM files to run UG examples.

\item 
Cleaned up code in inst/tests for modified functions, added new tests for new functions.
}}

\section{Version 1.0.0}{\itemize{
\item 
New package \code{csaw}, for de novo detection of differential binding in ChIP-seq data.
}}
