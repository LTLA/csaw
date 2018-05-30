#' @export
consolidateTests <- function(id.list, result.list, weight.list, FUN=combineTests, reindex="best", ...) 
# Consolidate results from multiple window widths.
# 
# written by Aaron Lun
# created 19 May 2018
{
    if (length(id.list)!=length(result.list)) {
        stop("lengths of 'id.list' and 'result.list' should be the same")
    }
    nids <- lengths(id.list)
    if (!all(nids==vapply(result.list, nrow, FUN.VALUE=0L))) {
        stop("length of each ID vector and nrows of each result table should match")
    }

    if (missing(weight.list)) {
        warning("'weight.list' should be specified or NULL")
        weight.list <- NULL
    } else if (!is.null(weight.list)) { 
        if (length(id.list)!=length(weight.list)) {
            stop("lengths of 'id.list' and 'weight.list' should be the same")
        }
        if (!all(nids==lengths(weight.list))) {
            stop("lengths of ID and weight vectors should match")
        }
    }

    output <- FUN(unlist(id.list), do.call(rbind, result.list), weight=unlist(weight.list), ...)
    output <- .reindex_entries(result.list, output, reindex=reindex)
    return(output)
}

#' @importFrom S4Vectors DataFrame
.reindex_entries <- function(result.list, output, reindex) 
# Reindexing to indicate the entry of 'result.list' of origin for 
# indices that refer to the 'rbinded' data.frame.
{
    if (any(reindex %in% colnames(output))) {
        per.hit <- c(0L, vapply(result.list, nrow, FUN.VALUE=0L))
        endpoint <- cumsum(per.hit)

        for (rdx in intersect(colnames(output), reindex)) {
            current <- output[[rdx]] 
            origin <- findInterval(current, endpoint, rightmost.closed=TRUE)
            output[[rdx]] <- DataFrame(origin=origin, row=current - per.hit[origin])
        }
    }
    return(output)
}

#' @export
#' @importFrom S4Vectors nRnode nLnode queryHits subjectHits Hits
consolidateOverlaps <- function(olap.list, result.list, weight.list=NULL, FUN=combineOverlaps, reindex="best", ...)
# Consolidate results from multiple overlap objects.
# 
# written by Aaron Lun
# created 19 May 2018
{
    if (length(olap.list)!=length(result.list)) {
        stop("lengths of 'olap.list' and 'result.list' should be the same")
    }
    nolaps <- lengths(olap.list)

    nLnodes <- vapply(olap.list, nLnode, FUN.VALUE=0L)
    nRnodes <- vapply(olap.list, nRnode, FUN.VALUE=0L)
    if (length(unique(nLnodes))>1L) { 
        stop("'olap.list' contains overlaps to different query sets")
    }
    if (!all(nRnodes==vapply(result.list, nrow, FUN.VALUE=0L))) {
        stop("'olap.list' subject set has different length to elements of 'result.list'")
    }

    if (missing(weight.list)) {
        warning("'weight.list' should be specified or NULL")
        weight.list <- NULL
    } else if (!is.null(weight.list)) { 
        if (length(olap.list)!=length(weight.list)) {
            stop("lengths of 'olap.list' and 'weight.list' should be the same")
        }
        if (!all(nolaps==lengths(weight.list))) {
            stop("lengths of ID and weight vectors should match")
        }
    }

    query.list <- lapply(olap.list, queryHits)
    subject.list <- vector("list", length(olap.list))
    last <- 0L
    for (x in seq_along(olap.list)) {
        subject.list[[x]] <- subjectHits(olap.list[[x]]) + last
        last <- last + nRnodes[x]
    }

    olap <- Hits(unlist(query.list), unlist(subject.list), nLnode=nLnodes[1], nRnode=last, sort.by.query=TRUE)
    output <- FUN(olap, do.call(rbind, result.list), o.weight=unlist(weight.list), ...)
    output <- .reindex_entries(result.list, output, reindex=reindex)
    return(output)
}
