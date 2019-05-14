#' @export
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
#' @importFrom utils head
getPESizes <- function(bam.file, param=readParam(pe="both")) 
# Reads a BAM file to determine the size of the PE fragments. 
# Returns a vector of sizes which can be plotted for diagnostics. 
# The length of the vector will also tell you how many read pairs were considered valid. 
# The number of interchromosomal pairs and unoriented reads are also reported.
# 
# written by Aaron Lun
# a long long time ago
{
    if (param$pe!="both") { stop("paired-end inputs required") }
    param <- reform(param, max.frag=Inf) # to get all fragment sizes.

    extracted.chrs <- .activeChrs(bam.file, param$restrict)
    norm.list <- vector("list", length(extracted.chrs))
    diagnostics <- list(inter.chr=0L, unoriented=0L, discarded=0L)

    for (i in seq_along(extracted.chrs)) {
        cur.chr <- names(extracted.chrs)[i]
        cur.gr <- GRanges(cur.chr, IRanges(1L, extracted.chrs[i]))
        out <- .extractPE(bam.file, cur.gr, param, diagnostics=TRUE)
        norm.list[[i]] <- pmin(out$size + out$pos, extracted.chrs[i] + 1L) - out$pos # capping insert by chromosome boundaries 

        so.far <- head(names(extracted.chrs), i-1L)
        inter.set <- out$diagnostics$inter.chr
        diagnostics$inter.chr <- diagnostics$inter.chr + sum(inter.set[!names(inter.set) %in% so.far])
        diagnostics$unoriented <- diagnostics$unoriented + out$diagnostics$unoriented
        diagnostics$discarded <- diagnostics$discarded + out$diagnostics$discarded
    }

    # Returning sizes and some diagnostic data.
    list(sizes=unlist(norm.list), diagnostics=unlist(diagnostics))
}
