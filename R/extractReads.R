#' @export
#' @importFrom GenomicRanges GRanges GRangesList
#' @importFrom IRanges IRanges overlapsAny
#' @importFrom GenomeInfoDb Seqinfo seqnames
#' @importFrom Rsamtools scanBamHeader
#' @importFrom S4Vectors Rle
#' @importFrom BiocGenerics strand strand<- start end
extractReads <- function(bam.file, region, ext=NA, param=readParam(), as.reads=FALSE)
# Exactly as specified. Takes a region and plots it in bimodal mode, with
# options for duplicate removal, mapping quality enhancement, colorization,
# and PE manipulation.
#
# written by Aaron Lun
# created 1 September 2014
{
    if (length(region)!=1L) {
        stop("exactly one range is required for read extraction")
    }
    if (as.logical(strand(region)!="*")) {
        warning("strandedness of region will be ignored, use param$forward instead")
        strand(region) <- "*"
    }
    ext.data <- .collateExt(1, ext)

    bam.file <- .make_BamFile(bam.file)
    chrs <- scanBamHeader(bam.file)$targets
    cur.chr <- as.character(seqnames(region)[1])
    if (length(param$restrict) && ! cur.chr %in% param$restrict) {
        warning("current chromosome not in restricted subset")
    }
    if (! cur.chr %in% names(chrs)) {
        stop("cannot find current chromosome in the BAM file header")
    }
    max.len <- chrs[[cur.chr]]
    sqi <- Seqinfo(cur.chr, max.len)

    # Kicking out the region by 'max.frag' or 'ext' to guarantee capture of all participants.
    if (param$pe=="both") {
        keepers <- c(param$max.frag, ext.data$final)
    } else {
        keepers <- c(ext.data$ext, ext.data$final)
    }
    keepers <- keepers[!is.na(keepers)]
    if (length(keepers)) {
        max.ext <- max(keepers)
    } else {
        max.ext <- 0L
    }
    actual.region <- GRanges(cur.chr, IRanges(max(1L, start(region)-max.ext), min(max.len, end(region)+max.ext)))

    if (param$pe!="both") {
        read.data <- .extractSE(bam.file, where=actual.region, param=param)
        forward.reads <- .extendSEdir(read.data$forward, ext=ext.data$ext[1], final=ext.data$final, chrlen=max.len, forward=TRUE)
        reverse.reads <- .extendSEdir(read.data$reverse, ext=ext.data$ext[1], final=ext.data$final, chrlen=max.len, forward=FALSE)

        f.chr <- Rle(cur.chr, length(forward.reads$start))
        f.reads <- GRanges(f.chr, IRanges(forward.reads$start, forward.reads$end), strand=Rle("+", length(f.chr)), seqinfo=sqi)
        r.chr <- Rle(cur.chr, length(reverse.reads$start))
        r.reads <- GRanges(r.chr, IRanges(reverse.reads$start, reverse.reads$end), strand=Rle("-", length(r.chr)), seqinfo=sqi)

        fkeep <- overlapsAny(f.reads, region)
        f.reads <- f.reads[fkeep]
        rkeep <- overlapsAny(r.reads, region)
        r.reads <- r.reads[rkeep]
        return(c(f.reads, r.reads))
    } else {
        frag.data <- .extractPE(bam.file, where=actual.region, param=param, with.reads=as.reads)
        cur.frags <- .coerceFragments(frag.data$pos, frag.data$pos + frag.data$size - 1L, final=ext.data$final, chrlen=max.len)

        # Filtering to retain those fragments that actually overlap.
        add.chr <- Rle(cur.chr, length(cur.frags$start))
        of.interest <- GRanges(add.chr, IRanges(cur.frags$start, cur.frags$end), seqinfo=sqi)
        keep <- overlapsAny(of.interest, region)
        if (!as.reads) { 
            return(of.interest[keep]) 
        }

        # Reporting the individual reads, if requested (but only for the *fragments* that overlap the region).
        forward <- GRanges(add.chr, IRanges(frag.data$forward$pos, pmin(frag.data$forward$pos+frag.data$forward$qwidth-1L, max.len)), 
            seqinfo=sqi, strand=Rle("+", length(add.chr)))
        reverse <- GRanges(add.chr, IRanges(frag.data$reverse$pos, pmin(frag.data$reverse$pos+frag.data$reverse$qwidth-1L, max.len)), 
            seqinfo=sqi, strand=Rle("-", length(add.chr)))
        out <- GRangesList(forward[keep], reverse[keep])
        names(out) <- c("forward", "reverse")
        return(out)
    }

    # Returning an empty set, otherwise.
    return(GRanges(seqinfo=sqi))
}

