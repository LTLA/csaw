#' @importFrom GenomeInfoDb seqnames
#' @importFrom BiocGenerics start end
.extractSE <- function(bam.file, where, param) 
# Extracts single-end read data from a BAM file with removal of unmapped,
# duplicate and poorly mapped/non-unique reads. We also discard reads in the
# specified discard regions. 
#
# written by Aaron Lun
# created 8 December 2013
{
    cur.chr <- as.character(seqnames(where)) 
    bam.file <- path.expand(bam.file)
    bam.index <- paste0(bam.file, ".bai")

    if (length(param$forward)==0L) { 
        stop("read strand extraction must be specified") 
    }

    if (param$pe=="first") {
        use.first <- TRUE
    } else if (param$pe=="second") { 
        use.first <- FALSE
    } else {
        use.first <- NA
    }

    out <- .Call(cxx_extract_single_data, bam.file, bam.index, 
        cur.chr, start(where), end(where), 
        param$minq, param$dedup, param$forward, use.first)

    names(out) <- c("forward", "reverse")
    names(out$forward) <- names(out$reverse) <- c("pos", "qwidth")

    disc <- param$discard
    if (length(disc)) {
        for (i in seq_along(out)) {
            curdata <- out[[i]]
            keep <- .filter_blacklist(cur.chr, curdata$pos, curdata$qwidth, disc)
            out[[i]] <- list_subset(curdata, keep)
        }
    }

    out
}

#' @importFrom IRanges overlapsAny IRanges
#' @importFrom GenomicRanges GRanges
.filter_blacklist <- function(chr, pos, size, disc) {
    if (!length(pos)) { return(logical(0)) } # avoid non-parallel errors with 'chr'.
    gr <- GRanges(chr, IRanges(pos, width=size))
    !overlapsAny(gr, disc, type="within")
}

#' @importFrom GenomicsRanges GRanges
#' @importFrom IRanges overlapsAny IRanges
#' @importFrom GenomeInfoDb seqnames
#' @importFrom BiocGenerics start end
.extractPE <- function(bam.file, where, param, with.reads=FALSE, diagnostics=FALSE)
# A function to extract PE data for a particular chromosome. Synchronisation
# is expected.  We avoid sorting by name  as it'd mean we have to process the
# entire genome at once (can't go chromosome-by-chromosome).  This probably
# results in increased memory usage across the board, and it doesn't fit in
# well with the rest of the pipelines which assume coordinate sorting.
# 
# written by Aaron Lun
# created 8 December 2013
{
    cur.chr <- as.character(seqnames(where)) 
    bam.file <- path.expand(bam.file)
    bam.index <- paste0(bam.file, ".bai")

    if (!identical(param$forward, NA)) { 
        stop("cannot specify read strand when 'pe=\"both\"'") 
    }

    out <- .Call(cxx_extract_pair_data, bam.file, bam.index, 
        cur.chr, start(where), end(where), 
        param$minq, param$dedup, 
        diagnostics)

    # Filtering out *fragments* by region.
    disc <- param$discard
    if (length(disc)) {
        starts <- out[[1]][[1]]
        sizes <- out[[2]][[1]] + out[[2]][[2]] - starts
        keep <- .filter_blacklist(cur.chr, starts, sizes, disc)
        out <- list_subset(out, keep)
    }

    if (diagnostics) {
        names(out) <- c("forward", "reverse", "total", "single", "unoriented", "one.unmapped", "inter.chr")
        return(out)
    }

    # Filtering out by size.
    left.pos <- out[[1]][[1]]
    left.len <- out[[1]][[2]]
    right.pos <- out[[2]][[1]]
    right.len <- out[[2]][[2]]
    all.sizes <- right.pos + right.len - left.pos

    keep <- all.sizes <= param$max.frag 
    output <- list(pos=left.pos[keep], size=all.sizes[keep])
    if (with.reads) {
        output$forward <- list(pos=left.pos[keep], qwidth=left.len[keep])
        output$reverse <- list(pos=right.pos[keep], qwidth=right.len[keep])
    }

    output
}

list_subset <- function(x, keep) {
    if (!is.list(x)) {
        return(x[keep])
    } 
    for (i in seq_along(x)) {
        x[[i]] <- list_subset(x[[i]], keep)
    }
    x
}

