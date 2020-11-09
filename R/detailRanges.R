#' @export
#' @importFrom IRanges promoters trim IRanges flank findOverlaps
#' @importFrom BiocGenerics strand start end
#' @importFrom GenomeInfoDb seqnames seqinfo
#' @importFrom GenomicRanges GRanges
#' @importFrom S4Vectors queryHits subjectHits
#' @importFrom methods is
detailRanges <- function(incoming, txdb, orgdb, dist=5000, promoter=c(3000, 1000), key.field="ENTREZID", name.field="SYMBOL", ignore.strand=TRUE)
# Gives three character vectors for each 'incoming'. 
# The first specifies which features are wholly or partially overlapped by the current range.
# The second specifies what features are within 'tol' of the 5' end of the incoming region.
# The last specifies what features are within 'tol' of the 3' end of the incoming region.
# Promoters and exons are included. 
# If there are multiple overlaps, the first and last exon matching is shown.
#
# written by Aaron Lun
# created 23 November 2013
{
    # Obtain exons.
    exon.ranges <- GenomicFeatures::exonsBy(txdb, by="gene")
    exon.ranges <- unlist(exon.ranges)
    exon.ranges$exon_id <- NULL
    exon.ranges$exon_name <- NULL

    # Obtain promoters.
    if (length(promoter)!=2L) {  
        stop("need upstream/downstream specification in 'promoter'") 
    }
    prom.ranges <- suppressWarnings(trim(promoters(txdb, upstream=promoter[1], downstream=promoter[2], column="gene_id")))
    expanded <- rep(seq_along(prom.ranges), lengths(prom.ranges$gene_id))
    prom.ranges <- prom.ranges[expanded]
    names(prom.ranges) <- unlist(prom.ranges$gene_id)
    prom.ranges$gene_id <- NULL

    # Obtain gene bodies.
    if (is(txdb, "TxDb")) { 
        gene.ranges <- GenomicFeatures::genes(txdb, single.strand.genes.only=FALSE)
        gene.ranges <- unlist(gene.ranges)
    } else {
        gene.ranges <- GenomicFeatures::genes(txdb)
    }

    # Assembling all ranges.
    all.ranges <- c(exon.ranges, prom.ranges, gene.ranges)
    range.type <- rep(c("E", "P", "G"), c(length(exon.ranges), length(prom.ranges), length(gene.ranges)))
    gene.names <- suppressMessages(AnnotationDbi::mapIds(orgdb, keys=names(all.ranges), column=name.field, keytype=key.field))  
    gene.names <- ifelse(is.na(gene.names), names(all.ranges), gene.names)

    all.ranges$symbol <- gene.names
    all.ranges$type <- range.type 
    all.ranges <- all.ranges[order(names(all.ranges))]
    
    # Returning the useful stuff, if no overlaps are requested.
    if (missing(incoming)) { 
        return(all.ranges)
    }

    # Computing overlaps, possibly with strandedness.
    full.lap <- findOverlaps(incoming, all.ranges, ignore.strand=ignore.strand)

    # Note that we don't do intronic or promoter overlaps to the flanks; we don't care that we're so-and-so base pairs away from the end of the promoter. 
    # We also rule out negative or zero distances i.e. those that would overlap the# region itself (zero distance means 1-based end and start are equal).
    # Flanking sets ignore.true so that we always get left/right flanks.
    flank.only <- all.ranges$type == "E"
    to.flank <- all.ranges[flank.only]

    left.flank <- suppressWarnings(trim(flank(incoming, dist, ignore.strand=TRUE)))
    left.lap <- findOverlaps(left.flank, to.flank, ignore.strand=ignore.strand)
    left.dist <- start(incoming)[queryHits(left.lap)] - end(to.flank)[subjectHits(left.lap)]
    left.nolap <- left.dist > 0L
    left.lap <- left.lap[left.nolap,]
    left.dist <- left.dist[left.nolap]

    right.flank <- suppressWarnings(trim(flank(incoming, dist, start=FALSE, ignore.strand=TRUE)))
    right.lap <- findOverlaps(right.flank, to.flank, ignore.strand=ignore.strand)
    right.dist <- start(to.flank)[subjectHits(right.lap)] - end(incoming)[queryHits(right.lap)]  
    right.nolap <- right.dist > 0L
    right.lap <- right.lap[right.nolap,]
    right.dist <- right.dist[right.nolap]
    
    # Collating the left-overs.
    re.name <- names(all.ranges)
    re.symbol <- all.ranges$symbol
    re.type <- all.ranges$type
    re.str <- as.character(strand(all.ranges))

    olap.str <- .Call(cxx_annotate_overlaps, length(incoming), queryHits(full.lap)-1L, subjectHits(full.lap)-1L, NULL, 
        re.name, re.symbol, re.type, re.str)
    left.str <- .Call(cxx_annotate_overlaps, length(incoming), queryHits(left.lap)-1L, which(flank.only)[subjectHits(left.lap)]-1L, left.dist, 
        re.name, re.symbol, re.type, re.str)
    right.str <- .Call(cxx_annotate_overlaps, length(incoming), queryHits(right.lap)-1L, which(flank.only)[subjectHits(right.lap)]-1L, right.dist, 
        re.name, re.symbol, re.type, re.str)
    return(list(overlap=olap.str, left=left.str, right=right.str))
}
