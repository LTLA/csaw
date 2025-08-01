#' @export
#' @importFrom BiocGenerics start end strand start<-
#' @importFrom S4Vectors split
#' @importFrom BiocParallel bpmapply bpisup bpstart bpstop SerialParam
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
#' @importFrom Seqinfo seqnames
profileSites <- function(bam.files, regions, param=readParam(), range=5000, ext=100, 
    average=TRUE, normalize="none", strand=c("ignore", "use", "match"), BPPARAM=SerialParam()) 
# Computes the coverage profile around putative binding sites. The 5' edge of the
# binding site is identified by counting reads into a window of size `width`, on the left and
# right of a given position, and determining if the right/left ratio is greater than 5. It then
# records the coverage of the resulting bases, up to `range`.
#
# written by Aaron Lun
# created 2 July 2014
{
    average <- as.logical(average)

    # A bit of work for strand-specificity.
    strand <- match.arg(strand)
    use.strand <- (strand!="ignore")
    match.strand <- (strand=="match")
    if (match.strand && length(param$forward)) { 
        stop("set forward=NULL in param for strand-specific profiling")  
    }
    if (use.strand) { 
        reverse <- strand(regions)=="-"
        if (any(reverse)) {
            reverse <- as.logical(reverse)
            rregs <- regions[reverse]
            start(rregs) <- end(rregs) # Using the 5' end of the reverse-stranded region.

            rprof <- Recall(bam.files=bam.files, regions=rregs, 
                param=reform(param, forward=ifelse(match.strand, FALSE, NA)), 
                range=range, ext=ext, average=average, normalize=normalize,
                strand="ignore", BPPARAM=BPPARAM) 

            if (any(!reverse)) { 
                fprof <- Recall(bam.files=bam.files, regions=regions[!reverse], 
                    param=reform(param, forward=ifelse(match.strand, TRUE, NA)), 
                    range=range, ext=ext, average=average, normalize=normalize, 
                    strand="ignore", BPPARAM=BPPARAM)
            } else { 
                fprof <- 0 
            }
            
            if (average) { 
                prop.rstr <- sum(reverse)/length(reverse) # Weighting by the number of regions with each strand.
                return(fprof * (1-prop.rstr) + rev(rprof) * prop.rstr) # Flipping the profile for reverse-strand.
            } else {
                total.len <- ncol(rprof)
                final.mat <- matrix(0L, length(regions), total.len)
                final.mat[reverse,] <- rprof[,rev(seq_len(total.len))]
                final.mat[!reverse,] <- fprof
                colnames(final.mat) <- colnames(fprof)
                return(final.mat)
            }
        }
    }

    # Determining the normalization type.
    norm.types <- c("none", "total", "max") 
    normalize <- match.arg(normalize, norm.types)
    norm.type <- match(normalize, norm.types)

    # Setting up.
    bam.files <- .make_BamFiles(bam.files)
    nbam <- length(bam.files)
    ext.data <- .collateExt(nbam, ext)
    range <- as.integer(range)
    if (range <= 0L) { 
        stop("range should be positive") 
    }

    if (average) { 
        total.profile <- numeric(range*2 + 1)
    } else {
        total.profile <- matrix(0L, length(regions), range*2 + 1)
    }
    indices <- split(seq_along(regions), seqnames(regions))

    if (!bpisup(BPPARAM)) {
        bpstart(BPPARAM)
        on.exit(bpstop(BPPARAM))
    }

    # Running through the chromosomes.
    extracted.chrs <- .activeChrs(bam.files, param$restrict)
    for (i in seq_along(extracted.chrs)) {
        chr <- names(extracted.chrs)[i]
        chosen <- indices[[chr]]
        if (!length(chosen)) { 
            next 
        }

        outlen <- extracted.chrs[i]
        where <- GRanges(chr, IRanges(1L, outlen))

        # Reading in the reads for the current chromosome for all the BAM files.
        bp.out <- bpmapply(FUN=.profile_sites, bam.file=bam.files, init.ext=ext.data$ext, 
            MoreArgs=list(where=where, param=param, final.ext=ext.data$final, outlen=outlen),
            BPPARAM=BPPARAM, SIMPLIFY=FALSE)

        starts <- lapply(bp.out, "[[", "starts")
        ends <- lapply(bp.out, "[[", "ends")

        # Pulling out the regions.
        all.starts <- start(regions)[chosen]
        os <- order(all.starts)
        all.starts <- all.starts[os]

        # We call the C++ functions to aggregate profiles.
        starts <- unlist(starts)
        ends <- unlist(ends)
        if (!length(starts)) { next }

        cur.profile <- .Call(cxx_get_profile, starts, ends, all.starts, range, average, norm.type) 
        if (average) { 
            total.profile <- total.profile + cur.profile
        } else {
            cur.profile <- t(cur.profile)
            cur.profile[os,] <- cur.profile
            total.profile[chosen,] <- cur.profile
        }
    }

    # Cleaning up and returning the profiles. 
    if (average) { 
        total.profile <- total.profile/length(regions)
        names(total.profile) <- (-range):range
    } else {
        colnames(total.profile) <- (-range):range
    }
    return(total.profile)
}

.profile_sites <- function(bam.file, where, param, init.ext, final.ext, outlen) {
    if (param$pe!="both") {
        reads <- .extractSE(bam.file, where=where, param=param)
        extended <- .extendSE(reads, ext=init.ext, final=final.ext, chrlen=outlen)
        start.pos <- extended$start
        end.pos <- extended$end
    } else {
        out <- .extractPE(bam.file, where=where, param=param)
        checked <- .coerceFragments(out$pos, out$pos+out$size-1L, final=final.ext, chrlen=outlen)
        start.pos <- checked$start
        end.pos <- checked$end
    }

    list(starts=start.pos, ends=end.pos)
}

#' @export
#' @importFrom stats median
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @importFrom S4Vectors DataFrame
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
wwhm <- function(profile, regions, ext=100, proportion=0.5, rlen=NULL)
# This function computes the window width at half its maximum. This uses
# the output of profileSites to get the full width of the peak; it then
# subtracts twice the extension length to obtain the window width. 
# 
# written by Aaron Lun
# created 2 March 2015
# last modified 23 July 2015
{
    if (proportion <= 0 | proportion >= 1) { stop("proportion should be between 0 and 1") }
    is.max <- which.max(profile)
    if (length(is.max)!=1L) { stop("profile cannot be empty or all-NA") }
    cutoff <- proportion * profile[is.max]
    above.max <- profile >= cutoff

    # Getting the width of the peak at half-max.
    out <- rle(above.max)
    ends <- cumsum(out$lengths)
    starts <- c(1L, ends[-length(ends)]+1L)
    chosen <- findInterval(is.max, starts)
    chosen.start <- starts[chosen]
    chosen.end <- ends[chosen]
    if (chosen.end==length(profile) || chosen.start==1L) {
        warning("width at specified proportion exceeds length of profile")
    }
    peak.width <- chosen.end - chosen.start + 1L

    # Getting the median size of the regions.
    if (!missing(regions)) {
        ref.size <- median(width(regions))
    } else {
        warning("regions not supplied, assuming width of 1 bp")
        ref.size <- 1L
    }

    # To get the average extension length across libraries, via getWidths.
    # Using a range of width 1 bp, so extension length is directly returned.
    # Setting start above 1, to future-proof against potential issues with extending before chromosome start.
    nlibs <- length(ext)
    ext.data <- .collateExt(nlibs, ext)
    dummy.data <- SummarizedExperiment(colData=DataFrame(ext=ext.data$ext), 
        metadata=list(final.ext=ext.data$final),
        rowRanges=GRanges("chrA", IRanges(start=100000, width=1)))
    if (!is.null(rlen)) { dummy.data$rlen <- rlen }
    ext.len <- getWidths(dummy.data)

    # Computing the window size. Add 2 to ensure one base overlaps the 
    # most extreme fragments on both sides. Subtract ref.size-1 as random
    # distribution of summits within maximal windows widens the peak.
    max(peak.width - ext.len*2L + 2L - ref.size + 1L, 1L)
}
