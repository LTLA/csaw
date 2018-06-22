#' @importFrom GenomeInfoDb seqnames
#' @importFrom BiocGenerics start end
.extractSE <- function(bam.file, where, param) 
# Extracts single-end read data from a BAM file with removal of unmapped,
# duplicate and poorly mapped/non-unique reads. We also discard reads in the
# specified discard regions. In such cases, the offending reads must be wholly
# within the repeat region.  We use the real alignment width, just in case we
# have very long reads in the alignment that are heavily soft-clipped (i.e., they
# should be reported as within but the read length will put them out).
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

    out <- .Call(cxx_extract_single_data, bam.file, bam.index, cur.chr,
            start(where), end(where), param$minq, param$dedup, param$forward, use.first) 
    names(out) <- c("forward", "reverse")
    names(out$forward) <- names(out$reverse) <- c("pos", "qwidth")

    # Filtering by discard regions. Using alignment width so long reads can escape repeats.
    for (i in names(out)) { 
        current <- out[[i]]
        keep <- .discardReads(cur.chr, current[[1]], current[[2]], param$discard)
        current <- lapply(current, "[", keep)
        out[[i]] <- current
    }
    return(out)
}

#' @importFrom GenomeInfoDb seqnames
#' @importFrom IRanges IRanges overlapsAny ranges
.discardReads <- function(chr, pos, alen, discard) {
    if (!length(pos)) { 
        return(logical(0)) 
    }
    relevant <- seqnames(discard)==chr
    if (any(relevant)) { 
        keep <- !overlapsAny(IRanges(pos, pos+alen-1L), ranges(discard)[relevant], type="within")
    } else {
        keep <- !logical(length(pos))
    }
    return(keep)
}

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
    out <- .Call(cxx_extract_pair_data, bam.file, bam.index, cur.chr,
            start(where), end(where), param$minq, param$dedup, diagnostics)

    if (diagnostics) {
        names(out) <- c("forward", "reverse", "total", "single", "ufirst", "usecond", "one.mapped", "ifirst", "isecond")
        return(out)
    }
    left.pos <- out[[1]][[1]]
    left.len <- out[[1]][[2]]
    right.pos <- out[[2]][[1]]
    right.len <- out[[2]][[2]]

    # Filtering by discard.
    dlkeep <- .discardReads(cur.chr, left.pos, left.len, param$discard)
    drkeep <- .discardReads(cur.chr, right.pos, right.len, param$discard)
    dkeep <- dlkeep & drkeep

    # Computing fragment sizes (implicit truncation of overrun reads).
    all.sizes <- right.pos + right.len - left.pos
    fkeep <- all.sizes <= param$max.frag 

    # Reporting output.
    keep <- dkeep & fkeep
    output <- list(pos=left.pos[keep], size=all.sizes[keep])
    if (with.reads) {
        output$forward <- list(pos=left.pos[keep], qwidth=left.len[keep])
        output$reverse <- list(pos=right.pos[keep], qwidth=right.len[keep])
    }
    return(output)
}

# Aliases, for convenience.

.getSingleEnd <- .extractSE

.getPairedEnd <- .extractPE
