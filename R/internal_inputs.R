#' @importFrom Rsamtools scanBamHeader
.activeChrs <- function(bam.files, restrict) 
# Processes the incoming data; checks that bam headers are all correct,
# truncates the list according to 'restrict'.
# 
# written by Aaron Lun
# created 12 December 2014
{ 
    originals <- NULL
    for (bam in bam.files) {
        chrs <- scanBamHeader(bam)[[1]][[1]]
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
# Converts the input to a GRanges, if it wasn't before.
{
    if (is(x, "RangedSummarizedExperiment")) { 
        x <- rowRanges(x)
    } else if (!is(x, "GenomicRanges")) {
        stop("'x' must be a RangedSummarizedExperiment or GRanges object")
    }
    return(x)
}
