#' @export
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
getPESizes <- function(bam.file, param=readParam(pe="both")) 
# Reads a BAM file to determine the size of the PE fragments. 
# Returns a vector of sizes which can be plotted for diagnostics. 
# The length of the vector will also tell you how many read pairs were considered valid. 
# The total number of reads, the number of singletons and the number of interchromosomal pairs are also reported.
# 
# written by Aaron Lun
# a long long time ago
{
    if (param$pe!="both") { stop("paired-end inputs required") }

    bam.file <- .make_BamFile(bam.file)
    extracted.chrs <- .activeChrs(list(bam.file), param$restrict)
    nchrs <- length(extracted.chrs)
    totals <- singles <- one.unmapped <- mapped <- unoriented <- 0L
    norm.list <- loose.names.1 <- loose.names.2 <- vector("list", nchrs)

    for (i in seq_len(nchrs)) { 
        cur.chr <- names(extracted.chrs)[i]
        output <- .extractPE(bam.file, GRanges(cur.chr, IRanges(1L, extracted.chrs[i])), param=param, diagnostics=TRUE)
        totals <- totals + output$total
        singles <- singles + output$single
        unoriented <- unoriented + output$unoriented
        one.unmapped <- one.unmapped + output$one.unmapped 

        # Valid read pairs get their sizes stored.
        all.sizes <- pmin(output$reverse[[1]] + output$reverse[[2]], extracted.chrs[i] + 1L) - output$forward[[1]] 
        norm.list[[i]] <- all.sizes

        # For inter-chromosomals; either mapped (and then store names), or implicitly unmapped.
        loose.names.1[[i]] <- output$inter.chr[[1]]
        loose.names.2[[i]] <- output$inter.chr[[2]]
    }

    # Checking whether a read is positively matched to a mapped counterpart on another chromosome.
    # If not, then it's just a read in an unmapped pair.
    loose.names.1 <- unlist(loose.names.1)
    loose.names.2 <- unlist(loose.names.2)
    inter.chr <- sum(loose.names.1 %in% loose.names.2)
    one.unmapped <- one.unmapped + length(loose.names.2) + length(loose.names.1) - inter.chr*2L

    bam.index <- path.expand(index(bam.file))
    bam.file <- path.expand(path(bam.file))
    out <- .Call(cxx_get_leftovers, bam.file, bam.index, names(extracted.chrs))
    totals <- totals + out

    norm.list <- unlist(norm.list)
    mapped <- singles + one.unmapped + 2L * (unoriented + inter.chr + length(norm.list))

    # Returning sizes and some diagnostic data.
    list(sizes=norm.list, 
        diagnostics=c(
            total.reads=totals, 
            mapped.reads=mapped, 
            single=singles, 
            mate.unmapped=one.unmapped, 
            unoriented=unoriented, 
            inter.chr=inter.chr
        )
    )
}
