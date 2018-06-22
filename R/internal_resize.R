.coerceFragments <- function(starts, ends, final, chrlen) 
# Coerces the fragments to the desired 'final.ext', and ensures
# that prior manipulations do not redefine fragment beyond chromosome 
# boundaries (e.g., due to read extension or rescaling).
#
# written by Aaron Lun
# created 13 February 2014
{
    if (!is.na(final)) { 
        remainders <- as.integer(ends - starts + 1L - final)
        if (any(remainders!=0L)) { 
            starts <- starts + as.integer(floor(remainders/2))
            ends <- ends - as.integer(ceiling(remainders/2))
        } 
    }

    chrlen <- as.integer(chrlen)
    starts <- pmin(pmax(1L, starts), chrlen)
    ends <- pmin(pmax(1L, ends), chrlen)
    return(list(start=starts, end=ends)) 
}

.collateExt <- function(nbam, ext)
# Collates the extension parameters into a set of ext and remainder values.
# The idea is to extend each read directionally to 'ext', and then extend in
# both directions by 'remainder' to reach the desired fragment length.
{
    if (is.list(ext)) {
        if (length(ext)!=2L) {
            stop("'ext' must be a list of length 2")
        } 
        final.ext <- ext[[2]]
        ext <- ext[[1]]
    } else {
        if (length(ext)!=1L) {
            stop("'ext' must be an integer scalar")
        }
        ext <- rep(ext, length.out=nbam)
        final.ext <- NA_integer_
    }

    ext <- as.integer(round(ext))
    if (length(ext)!=nbam) {
        stop("length of extension vector is not consistent with the number of libraries")
    } 
    if (any(!is.na(ext) & ext <= 0L)) { 
        stop("extension length must be NA or a positive integer") 
    }

    final.ext <- as.integer(round(final.ext))
    if (length(final.ext)!=1L || (!is.na(final.ext) && final.ext <= 0L)) { 
        stop("final extension length must be a positive integer or NA") 
    }

    list(ext=ext, final=final.ext)
}

.extendSEdir <- function(reads, ext, final, chrlen, forward=TRUE) {
    if (is.na(ext)) { 
        start <- reads$pos
        end <- reads$pos + reads$qwidth -1L
    } else {
        if (forward) {
            start <- reads$pos
            end <- reads$pos + ext - 1L
        } else {
            end <- reads$pos + reads$qwidth - 1L
            start <- end - ext + 1L
        }
    }
    out <- .coerceFragments(start, end, final=final, chrlen=chrlen)
    return(out)
}

.extendSE <- function(reads, ext, final, chrlen)
# This decides how long to extend reads. The addition of the remainder kicks
# out (or truncates) the fragments to reach the desired 'final.ext'. If 'ext'
# is NA, the read length is used instead.
#
# written by Aaron Lun
# created 12 December 2014
{
    fout <- .extendSEdir(reads$forward, ext, final, chrlen, forward=TRUE)
    rout <- .extendSEdir(reads$reverse, ext, final, chrlen, forward=FALSE)
    mapply(c, fout, rout, SIMPLIFY=FALSE)
}
