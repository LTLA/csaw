#' @export
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @importFrom GenomeInfoDb seqlevels<- seqlengths<-
#' @importFrom BiocParallel bpmapply bpisup bpstart bpstop SerialParam
windowCounts <- function(bam.files, spacing=50, width=spacing, ext=100, shift=0, filter=10, bin=FALSE, 
     param=readParam(), BPPARAM=SerialParam())
# Gets counts from BAM files at each position of the sliding window. 
# Appliesa gentle filter to remove the bulk of window positions with low counts.
# Returns a RangedSummarizedExperiment object with counts and genomic coordinates.
# 
# written by Aaron Lun
# created 5 April 2012
{
    shift <- as.integer(shift)
    width <- as.integer(width)
    spacing <- as.integer(spacing)

    if (bin) { # A convenience flag for binning.
        spacing <- width
        ext <- 1L
        filter <- min(1, filter)
    }

    if (shift < 0L) { stop("shift must be a non-negative integer") }
    if (width <= 0L) { stop("width must be a positive integer") }
    if (spacing <= 0L) { stop("spacing must be a positive integer") }
    if (shift >= spacing) { stop("shift must be less than the spacing") } # avoid redundant windows, see POINT 1 below.
    at.start <- .is_pt_at_start(shift, width) # see POINT 2

    bam.files <- .make_BamFiles(bam.files)
    nbam <- length(bam.files)
    ext.data <- .collateExt(nbam, ext)

    if (!bpisup(BPPARAM)) {
        bpstart(BPPARAM)
        on.exit(bpstop(BPPARAM))
    }

    # Initializing various collectable containers. These are non-empty so the class 
    # will be right even if no chromosomes are around and the loop is empty.
    totals <- integer(nbam)
    extracted.chrs <- .activeChrs(bam.files, param$restrict)
    nchrs <- length(extracted.chrs)

    all.out <- rep(list(matrix(0L, ncol=nbam, nrow=0)), nchrs)
    all.regions <- rep(list(GRanges()), nchrs)
    all.lengths <- rep(list(vector("list", nchrs)), nbam)

    for (i in seq_len(nchrs)) { 
        chr <- names(extracted.chrs)[i]
        outlen <- extracted.chrs[i]
        where <- GRanges(chr, IRanges(1, outlen))
        total.pts <- .get_total_pts(outlen, shift, spacing, at.start) # see POINT 3

        # Parallelized loading.
        bp.out <- bpmapply(FUN=.window_counts, bam.file=bam.files, init.ext=ext.data$ext, 
            MoreArgs=list(where=where, param=param, 
                final.ext=ext.data$final, outlen=outlen, bin=bin, 
                shift=shift, width=width, spacing=spacing, 
                total.pts=total.pts, at.start=at.start),
            BPPARAM=BPPARAM, SIMPLIFY=FALSE)

        outcome <- matrix(0L, total.pts, nbam)
        for (bf in seq_along(bp.out)) {
            outcome[,bf] <- bp.out[[bf]]$counts
            totals[bf] <- totals[bf] + bp.out[[bf]]$totals
            all.lengths[[bf]][[i]] <- bp.out[[bf]]$lengths
        }

        # Filtering on row sums (for memory efficiency). 
        keep <- rowSums(outcome)>=filter 
        if (!any(keep)) { 
            next 
        } else if (!all(keep)) { 
            outcome <- outcome[keep,,drop=FALSE] 
        }
        all.out[[i]] <- outcome

        # Defining genomic coordinates.
        center <- (which(keep) - 1L) * spacing + 1L + ifelse(at.start, 0L, spacing) 
        reg.start <- center - shift
        reg.end <- pmin(outlen, reg.start + width - 1L)
        reg.start <- pmax(1L, reg.start)
        all.regions[[i]] <- GRanges(chr, IRanges(reg.start, reg.end))
    }

    # Generating the remaining GRanges for output (suppressing numerous warnings).
    all.regions <- suppressWarnings(do.call(c, all.regions))
    seqlevels(all.regions) <- names(extracted.chrs)
    seqlengths(all.regions) <- extracted.chrs
    strand(all.regions) <- .decideStrand(param)

    SummarizedExperiment(
        assays=SimpleList(counts=do.call(rbind, all.out)), 
        rowRanges=all.regions, 
        colData=.formatColData(bam.files, totals, ext.data, all.lengths, param),
        metadata=list(
            spacing=spacing, width=width, shift=shift, bin=bin, 
            param=param, 
            final.ext=ifelse(bin, 1L, ext.data$final) # For getWidths with paired-end binning.
        )
    ) 
}

# POINT 1: redundant windows are avoided by enforcing 'shift < spacing'. 
# A window that is shifted to cover the whole chromosome cannot be spaced 
# to a new position where it still covers the whole chromosome. The next 
# one must be inside the chromosome.

# POINT 2: Need to account for the possible loss of a window from the front when
# the shift is non-zero, because the corresponding window is wholly outside the
# chromosome (i.e., shifted so that the width of the window is before position 1).
.is_pt_at_start <- function(shift, width) {
    shift < width 
}
 
# POINT 3: Accounting for the possible gain of a centrepoint from the back when
# shift is non-zero, i.e., does the shift bring the next window start under outlen?
# i.e., given 1L - shift + spacing * X <= outlen, solve for the largest integer X.
.get_total_pts <- function(outlen, shift, spacing, at.start) {
    as.integer(floor((outlen + shift - 1L)/spacing)) + at.start # if the extra point exists at the start.
}

.window_counts <- function(bam.file, where, param, 
        init.ext, final.ext, outlen, bin, 
        shift, width, spacing, 
        total.pts, at.start) 
# Internal function to ensure csaw namespace carries into bplapply.
{
    if (param$pe!="both") {
        reads <- .extractSE(bam.file, where=where, param=param)
        extended <- .extendSE(reads, ext=init.ext, final=final.ext, chrlen=outlen)
        frag.start <- extended$start
        frag.end <- extended$end

        rlengths <- cbind(c(mean(reads$forward$qwidth), mean(reads$reverse$qwidth)),
                       c(length(reads$forward$qwidth), length(reads$reverse$qwidth)))
    } else {
        out <- .extractPE(bam.file, where=where, param=param)
        if (bin) { 
            # Only want to record each pair once in a bin, so forcing it to only use the midpoint.
            mid <- as.integer(out$pos + out$size/2)
            mid <- pmin(mid, outlen)
            frag.end <- frag.start <- mid
        } else {
            checked <- .coerceFragments(out$pos, out$pos+out$size-1L, final=final.ext, chrlen=outlen)
            frag.start <- checked$start
            frag.end <- checked$end
        }
        rlengths <- c(mean(out$size), length(out$size))
    }

    # Windows are parametrized in terms of extension from the window start.
    # Non-unity window width corresponds to extension to the right.
    # Shifting causes it to be extended to the left, at the expense of the right.
    right <- width - shift - 1L
    left <- shift

    # We make it even simpler by extending the reads, which means that we don't 
    # have to actually change the window starts. The *start* of each read must be 
    # extended by 'right' and the *end* of each read must be extended by 'left'. 
    # (This mirrors the extension of the windows.)
    frag.start <- frag.start - right
    frag.end <- frag.end + left    
    
    # We pull out counts at the specified spacing. We do have to keep track of 
    # whether or not we want to use the first point, though.
    out <- .Call(cxx_get_rle_counts, frag.start, frag.end, total.pts, spacing, at.start)
    return(list(counts=out, totals=length(frag.start), lengths=rlengths))
}

