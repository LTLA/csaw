.decideStrand <- function(param) 
# Decides what strand we should assign to the output GRanges in the
# SummarizedExperiment object, after counting.
{
    getfs <- param$forward
    if (length(getfs)==0L) { 
        stop("unspecified strandedness")
    } else if (is.na(getfs)) { 
        return("*") 
    } else if (getfs) { 
        return("+") 
    } else { 
        return("-") 
    }
}

#' @importFrom stats weighted.mean
#' @importFrom S4Vectors DataFrame
#' @importFrom BiocGenerics path
.formatColData <- function(bam.files, totals, ext.data, all.extras, param) 
# Formats the column data in the output SummarizedExperiment after counting.
{
    nbam <- length(bam.files)
    store.ext <- ext.data$ext
    store.rlen <- rep(NA_integer_, nbam)

    store.extras <- numeric(nbam)
    for (bf in seq_len(nbam)) {
        current.extras <- do.call(rbind, all.extras[[bf]])
        store.extras[bf] <- weighted.mean(current.extras[,1], current.extras[,2])
    }
    store.extras <- as.integer(round(store.extras))

    if (param$pe=="both") { 
        store.ext <- store.extras 
    } else { 
        store.rlen <- store.extras 
    }

    if (!is.character(bam.files)) {
        bam.files <- vapply(bam.files, path, "")
    }
    DataFrame(bam.files=bam.files, totals=totals, ext=store.ext, rlen=store.rlen)
}
