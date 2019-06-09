#' @export
filterWindows <- function(data, background, type="global", assay.data="counts", assay.back="counts",
                          prior.count=2, scale.info=NULL)
# This is a function for proportion- or background-based filtering of a
# RangedSummarizedExperiment object. For the former, it computes the relative
# ranks that can be used to determine the proportion of highest-abundance
# windows to keep. For the latter, it returns the enrichment term between data
# and background.
#
# written by Aaron Lun
# created 18 February 2015	
{
	type <- match.arg(type, c("global", "local", "control", "proportion"))

	if (type=="proportion") {
        .Deprecated("filterWindowsProportion")
        filterWindowsProportion(data, assay.data=assay.data, prior.count=prior.count)
	} else if (type=="global") {
        .Deprecated("filterWindowsGlobal")
        filterWindowsGlobal(data, background, assay.data=assay.data, assay.back=assay.back, prior.count=prior.count)
    } else if (type=="local") {
        .Deprecated("filterWindowsLocal")
        filterWindowsLocal(data, background, assay.data=assay.data, assay.back=assay.back, prior.count=prior.count)
    } else {
        .Deprecated("filterWindowsControl")
        filterWindowsControl(data, background, assay.data=assay.data, assay.back=assay.back, 
            prior.count=prior.count, scale.info=scale.info)
	}
}

#' @export
filterWindowsProportion <- function(data, assay.data="counts", prior.count=2) 
# Defines the filter statistic for each window as the relative rank.
{
    abundances <- scaledAverage(data, assay.id=assay.data, scale=1, prior.count=prior.count)
    genome.windows <- .getWindowNum(data)
    ranked.win <- rank(abundances)
    nwin <- nrow(data)

    if (genome.windows <= nwin) {
        relative.rank <- ranked.win/nwin
    } else {
        relative.rank <- 1 + (ranked.win - nwin)/genome.windows
    }
    list(abundances=abundances, filter=relative.rank)
}

#' @export
filterWindowsGlobal <- function(data, background, assay.data="counts", assay.back="counts", prior.count=2) 
# Defines the filter statistic for each window as the log-fold change
# above a global estimate for the background enrichment.
{
    abundances <- scaledAverage(data, assay.id=assay.data, scale=1, prior.count=prior.count)

    if (missing(background)) {
        filter.stat <- abundances  - .getGlobalBg(data, abundances, prior.count)
        return(list(abundances=abundances, filter=filter.stat))
    } 

    # We don't adjust each abundance for the size of the individual regions, 
    # as that would be a bit too much like FPKM and make strong assumptions 
    # about the uniformity of coverage with respect to region width. We only 
    # use the widths to adjust the background abundances for comparison.
    bwidth <- getWidths(background)
    dwidth <- getWidths(data)

    .checkLibSizes(data, background)
    relative.width <- median(bwidth)/median(dwidth)
    if (is.na(relative.width)) { 
        relative.width <- 1 
    }

    bg.ab <- scaledAverage(background, assay.id=assay.back, scale=relative.width, prior.count=prior.count)
    filter.stat <- abundances - .getGlobalBg(background, bg.ab, prior.count)

    list(abundances=abundances, back.abundances=bg.ab, filter=filter.stat)
}

#' @export
filterWindowsLocal <- function(data, background, assay.data="counts", assay.back="counts", prior.count=2)
# Defines the filter statistic for each window as the log-fold change
# above a local estimate for the background enrichment, based on 
# flanking regions.
{
    .checkNestedRanges(data, background)
    .checkLibSizes(data, background)

    bwidth <- getWidths(background)
    dwidth <- getWidths(data)
    relative.width <- (bwidth  - dwidth)/dwidth
    assay(background, i=assay.back) <- (assay(background, i=assay.back, withDimnames=FALSE) 
        - assay(data, assay=assay.data, withDimnames=FALSE))

    # Some protection for negative widths (counts should be zero, so only the prior gets involved in bg.ab).
    subzero <- relative.width <= 0
    if (any(subzero)) { 
        relative.width[subzero] <- 1
        bg.y$counts[subzero,] <- 0L
    }	

    abundances <- scaledAverage(data, assay.id=assay.data, scale=1, prior.count=prior.count)
    bg.ab <- scaledAverage(background, assay.id=assay.back, scale=relative.width, prior.count=prior.count)
    filter.stat <- abundances - bg.ab
    list(abundances=abundances, back.abundances=bg.ab, filter=filter.stat)
}

#' @export
filterWindowsControl <- function(data, background, assay.data="counts", assay.back="counts",
    prior.count=2, scale.info=NULL)
# Defines the filter statistic for each window as the log-fold change
# above a local estimate for the background enrichment, based on 
# control libraries at the same position.
{
    if (!is.null(scale.info)) {
        if (!identical(scale.info$data.totals, data$totals)) {
            stop("'data$totals' are not the same as those used for scaling")
        } 
        if (!identical(scale.info$back.totals, background$totals)) { 
            stop("'back$totals' are not the same as those used for scaling")
        }
        background$totals <- background$totals * scale.info$scale 
    } else {
        warning("normalization factor not specified for composition bias")
    }
    .checkNestedRanges(data, background)

    bwidth <- getWidths(background)
    dwidth <- getWidths(data)
    relative.width <- bwidth/dwidth
    lib.adjust <- prior.count * mean(background$totals)/mean(data$totals) # Account for library size differences.

    abundances <- scaledAverage(data, assay.id=assay.data, scale=1, prior.count=prior.count)
    bg.ab <- scaledAverage(background, assay.id=assay.back, scale=relative.width, prior.count=lib.adjust)
    filter.stat <- abundances - bg.ab 

    list(abundances=abundances, back.abundances=bg.ab, filter=filter.stat)
}

##########################
### Internal functions ###
##########################

.checkLibSizes <- function(data, background) {
	if (!identical(data$totals, background$totals)) { 
		stop("'data$totals' and 'background$totals' should be identical")
	}
	return(NULL)
}

#' @importFrom SummarizedExperiment assay assay<- rowRanges
#' @importFrom BiocGenerics width
#' @importFrom IRanges pintersect
.checkNestedRanges <- function(data, background) {
    if (!identical(nrow(data), nrow(background))) { 
        stop("'data' and 'background' should be of the same length") 
    }	
    if (!identical(width(pintersect(rowRanges(data), rowRanges(background))), width(rowRanges(data)))) {
        stop("'rowRanges(data)' are not nested within 'rowRanges(background)'")
    }
    NULL
}

#' @importFrom GenomeInfoDb seqlengths
#' @importFrom SummarizedExperiment rowRanges
#' @importFrom S4Vectors metadata
.getWindowNum <- function(data) 
# Get the total number of windows, to account for those not 
# reported in windowCounts (for empty windows/those lost by filter > 1).
{
	spacing <- metadata(data)$spacing
	if (is.null(spacing)) { stop("failed to find spacing for windows") }
	sum(ceiling(seqlengths(rowRanges(data))/spacing)) 
}

#' @importFrom edgeR aveLogCPM
#' @importFrom stats quantile
.getGlobalBg <- function(data, ab, prior.count)
# Getting the quantile of those windows that were seen, corresponding to 
# the median of all windows in the genome. Assumes that all lost windows
# have lower abundances than those that are seen, which should be the
# case for binned data (where all those lost have zero counts).
{
	prop.seen <- length(ab)/.getWindowNum(data)
 	if (prop.seen > 1) { return(median(ab)) }
	if (prop.seen < 0.5) { return(aveLogCPM(rbind(integer(ncol(data))), lib.size=data$totals, prior.count=prior.count)) }
	quantile(ab, probs=1 - 0.5/prop.seen) 
}

###############################
### Miscellaneous functions ###
###############################

#' @export
#' @importFrom stats median
scaleControlFilter <- function(data.bin, back.bin, assay.data="counts", assay.back="counts")
# Computes the normalization factor due to composition bias
# between the ChIP and background samples.
{
    adjusted <- filterWindowsControl(data.bin, back.bin, prior.count=0, 
        assay.data=assay.data, assay.back=assay.back,
        scale.info=list(scale=1, data.totals=data.bin$totals, back.totals=back.bin$totals))
    nf <- 2^-median(adjusted$filter, na.rm=TRUE) # protect against NA's from all-zero. 
    list(scale=nf, data.totals=data.bin$totals, back.totals=back.bin$totals)
}
