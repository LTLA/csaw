#' @importFrom Rsamtools index 
#' @importFrom BiocGenerics path
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

    bam.file <- .make_BamFile(bam.file)
    bam.index <- path.expand(index(bam.file))
    bam.file <- path.expand(path(bam.file))

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

    cur.discard <- .getDiscard(param, cur.chr)

    out <- .Call(cxx_extract_single_data, bam.file, bam.index, 
        cur.chr, start(where), end(where), 
        param$minq, param$dedup, param$forward, use.first,
        cur.discard$pos, cur.discard$id)

    names(out) <- c("forward", "reverse")
    names(out$forward) <- names(out$reverse) <- c("pos", "qwidth")
    return(out)
}

.getDiscard <- function(param, chr) {    
    cur.discard <- param$processed.discard[[chr]]
    if (is.null(cur.discard)) {
        cur.discard <- list(pos=integer(0), id=integer(0))
    }
    cur.discard
}

#' @importFrom Rsamtools index 
#' @importFrom BiocGenerics path
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

    bam.file <- .make_BamFile(bam.file)
    bam.index <- path.expand(index(bam.file))
    bam.file <- path.expand(path(bam.file))

    if (!identical(param$forward, NA)) { 
        stop("cannot specify read strand when 'pe=\"both\"'") 
    }

    cur.discard <- .getDiscard(param, cur.chr)

    out <- .Call(cxx_extract_pair_data, bam.file, bam.index, 
        cur.chr, start(where), end(where), 
        param$minq, param$dedup, 
        cur.discard$pos, cur.discard$id, 
        diagnostics)

    if (diagnostics) {
        names(out) <- c("forward", "reverse", "total", "single", "unoriented", "one.unmapped", "inter.chr")
        return(out)
    }

    left.pos <- out[[1]][[1]]
    left.len <- out[[1]][[2]]
    right.pos <- out[[2]][[1]]
    right.len <- out[[2]][[2]]

    # Computing fragment sizes.
    all.sizes <- right.pos + right.len - left.pos
    keep <- all.sizes <= param$max.frag 
    output <- list(pos=left.pos[keep], size=all.sizes[keep])
    if (with.reads) {
        output$forward <- list(pos=left.pos[keep], qwidth=left.len[keep])
        output$reverse <- list(pos=right.pos[keep], qwidth=right.len[keep])
    }
    return(output)
}
