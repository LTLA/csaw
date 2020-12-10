#' @importFrom Rsamtools BamFile
#' @importClassesFrom Rsamtools BamFile
.make_BamFile <- function(bam.file) {
    if (!is(bam.file, "BamFile")) {
        BamFile(bam.file)
    } else {
        bam.file
    }
}

.make_BamFiles <- function(bam.files) {
    lapply(bam.files, .make_BamFile)
}

#' @importFrom Rsamtools scanBamHeader
.activeChrs <- function(bam.files, restrict) 
# Processes the incoming data; checks that bam headers are all correct, truncates the list according to 'restrict'.
# 
# written by Aaron Lun
# created 12 December 2014
{ 
    originals <- NULL
    bam.files <- .make_BamFiles(bam.files)

    for (bam in bam.files) {
        chrs <- scanBamHeader(bam)$targets
        chrs <- chrs[order(names(chrs))]

        if (is.null(originals)) { 
            originals <- chrs 
        } else if (!identical(originals, chrs)) { 
            warning("chromosomes are not identical between BAM files")
            pairing <- match(names(originals), names(chrs))
            originals <- pmin(originals[!is.na(pairing)], chrs[pairing[!is.na(pairing)]])
        }
    }
    if (length(restrict)) { 
        originals <- originals[names(originals) %in% restrict] 
    }
    return(originals)
}

#' @importClassesFrom SummarizedExperiment RangedSummarizedExperiment
#' @importClassesFrom GenomicRanges GenomicRanges
#' @importFrom methods is
.toGRanges <- function(x) 
# Converts the input to a GRanges, if it wasn't already one. 
{
    if (is(x, "RangedSummarizedExperiment")) { 
        x <- rowRanges(x)
    } else if (!is(x, "GenomicRanges")) {
        stop("'x' must be a RangedSummarizedExperiment or GRanges object")
    }
    return(x)
}

.parseFCcol <- function(fc.col, tab, multiple=TRUE) 
# Checks 'tab' for a column containing log-fold changes.
# Checks for any "logFC.*"-named column if a fc.col=NULL.
{
    if (is.null(fc.col)) { 
        fc.col <- grep("logFC", colnames(tab))
    } else if (is.character(fc.col)) { 
        fc.col <- match(fc.col, colnames(tab)) 
        if (any(is.na(fc.col))) { stop("failed to match logFC column names") }
    }
    if (!multiple && length(fc.col)!=1L) {
        stop("exactly one column should be specified by 'fc.col'")
    }
    as.integer(fc.col)
}

.getPValCol <- function(pval.col, tab) 
# Checks 'tab' for a column containing p-values.
# Checks for any "PValue"-named column if pval.col=NULL.
{
    if (length(pval.col)>1L) { 
        stop("multiple p-value columns are not supported")
    }
	if (is.null(pval.col)) { 
		pval.col <- "PValue"
    }
    if (is.character(pval.col)) {
        pval.col <- which(colnames(tab)==pval.col)
    } else {
        pval.col <- as.integer(pval.col) # coerce to integer, just in case.
    }
    if (length(pval.col)==0) { 
        stop("failed to find any p-value column")
    }
    return(pval.col)
}
