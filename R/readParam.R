# This defines the readParam class, which specifies the universal
# parameters for read loading in csaw. A class is used here so that
# there are continuous validity checks on the list values.

#' @export
#' @importClassesFrom BiocParallel BiocParallelParam
#' @importClassesFrom GenomicRanges GRanges
setClass("readParam", representation(
    pe="character", max.frag="integer",
    dedup="logical", minq="integer", forward="logical", 
	restrict="character", 
    discard="GRanges", 
    processed.discard="list"))

setValidity("readParam", function(object) {
    msg <- NULL

    if (length(object@pe)!=1L || ! object@pe %in%c("none", "both", "first", "second")) { 
		msg <- c(msg, "PE specification must be a character scalar of 'none', 'both', 'first' or 'second'") 
	}
   	if (length(object@max.frag)!=1L || object@max.frag <= 0L) {
		msg <- c(msg, "maximum fragment specifier must be a positive integer")
	} 

	if (length(object@dedup)!=1L || !is.logical(object@dedup)) { 
		msg <- c(msg, "duplicate removal specification must be a logical scalar")
	}
	if (length(object@minq)!=1L || !is.numeric(object@minq)) { 
		msg <- c(msg, "minimum mapping quality must be a numeric scalar")
	}

	if (length(object@forward)>1L || !is.logical(object@forward)) { 
		msg <- c(msg, "forward strand specification must be a logical scalar or NULL")
	} else if ((length(object@forward)==0L || !is.na(object@forward)) && object@pe == "both") {
		msg <- c(msg, "strand-specific extraction is not supported for paired-end data")
	}

    if (length(object@processed.discard)==0L && length(object@discard)!=0) {
        msg <- c(msg, "discard ranges has not been properly processed")
    }

    if (length(msg)) {
        return(msg)
    }
	return(TRUE)
})

#' @export
setMethod("$", signature("readParam"), function(x, name) { 
	slot(x, name)
})

#' @export
#' @importFrom BiocParallel bpnworkers
setMethod("show", signature("readParam"), function(object) {
	cat("    ", switch(object@pe,
 	   none="Extracting reads in single-end mode",
	   both="Extracting reads in paired-end mode",
	   first="Extracting the first read of each pair",
	   second="Extracting the second read of each pair"), "\n", sep="")

	if (object@pe=="both") {
		cat("        Maximum allowed distance between paired reads is", object@max.frag, "bp\n")
	}

	cat("    Duplicate removal is turned", ifelse(object@dedup, "on", "off"), "\n")
	if (is.na(object@minq)) { 
		cat("    No minimum threshold is set on the mapping score\n")
	} else {
		cat("    Minimum allowed mapping score is", object@minq, "\n")
	}

	if (length(object@forward)==0L) { 
		cat("    Reads are extracted from either strand separately\n")
	} else if (is.na(object@forward)) { 
		cat("    Reads are extracted from both strands\n")
	} else {
		cat("    Reads are extracted from the", ifelse(object@forward, "forward", "reverse"), "strand only\n")
	}

	rl <- length(object@restrict)
	if (rl) { 
		cat("    Read extraction is limited to", rl, ifelse(rl==1L, "sequence\n", "sequences\n"))
	} else {
		cat("    No restrictions are placed on read extraction\n")
	}

	dl <- length(object@discard)
	if (dl) { 
		cat("    Reads in", dl, ifelse(dl==1L, "region", "regions"), "will be discarded\n")
	} else {
		cat("    No regions are specified to discard reads\n")
	}
})

#' @export
#' @importFrom BiocParallel SerialParam
readParam <- function(pe="none", max.frag=500, dedup=FALSE, minq=NA, forward=NA, restrict=NULL, discard=GRanges())
# This creates a list of parameters, formally represented as a readParam
# object, specifying how reads should be extracted from the BAM files. The
# aim is to synchronize read loading throughout the package, such that
# you don't have to manually respecify them in each function.
#
# written by Aaron Lun
# created 1 September 2014
{
	max.frag <- as.integer(max.frag)
	dedup <- as.logical(dedup)
	forward <- as.logical(forward)
	minq <- as.integer(minq)
	restrict <- as.character(restrict) 
	new("readParam", pe=pe, max.frag=max.frag, 
		dedup=dedup, forward=forward, minq=minq, 
		restrict=restrict, 
        discard=discard, 
        processed.discard=.setupDiscard(discard))
}

#' @importFrom GenomeInfoDb seqnames
#' @importFrom BiocGenerics start end
.setupDiscard <- function(discard) {
    by.chr <- split(discard, seqnames(discard))
    output <- vector("list", length(by.chr))
    names(output) <- names(by.chr)

    for (x in names(output)) {
        cur.discard <- by.chr[[x]]
        cur.pos <- c(start(cur.discard), end(cur.discard)+1L) # 1-based positions.
        cur.ids <- rep(seq_along(cur.discard) - 1L, 2) # zero indexed elements.
        o <- order(cur.pos)
        output[[x]] <- list(pos=cur.pos[o], id=cur.ids[o])
    }

    output
}

#' @export
setGeneric("reform", function(x, ...) { standardGeneric("reform") })

#' @export
setMethod("reform", signature("readParam"), function(x, ...) {
	incoming <- list(...)
	sn <- slotNames(x)
	for (sx in names(incoming)) {
		val <- incoming[[sx]]
		sx <- match.arg(sx, sn)
		incoming[[sx]] <- switch(sx, 
			max.frag=as.integer(val),
			dedup=as.logical(val),
			forward=as.logical(val),
			minq=as.integer(val),
			restrict=as.character(val),
			val)
	}

    # Strictly internal.
    incoming$processed.discard <- NULL
    if (!is.null(incoming$discard)) {
        incoming$processed.discard <- .setupDiscard(incoming$discard)
    }

	do.call(initialize, c(x, incoming))
})

