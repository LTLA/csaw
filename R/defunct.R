#' Defunct functions
#'
#' Functions that have passed on to the function afterlife.
#' Their successors are also listed.
#'
#' @param ... Ignored arguments.
#'
#' @details
#' \code{filterWindows} is succeeded by \code{\link{filterWindowsGlobal}} and related functions, which provide a more focused programmatic interface.
#'
#' \code{consolidateWindows} is succeeded by \code{\link{mergeWindowsList}} and \code{\link{findOverlapsList}}.
#'
#' \code{consolidateTests} and \code{consolidateOverlaps} are succeeded by \code{\link{mergeResultsList}} and \code{\link{overlapResultsList}}, respectively.
#'
#' @return All functions error out with a defunct message pointing towards its descendent (if available).
#'
#' @author Aaron Lun
#'
#' @examples
#' try(filterWindows())
#' @name defunct
NULL

#' @export
#' @rdname defunct
filterWindows <- function(...) {
    .Defunct(new="filterWindowsGlobal")
}

#' @export
#' @rdname defunct
consolidateWindows <- function(...) {
    .Defunct(new="mergeWindowsList")
}

#' @export
#' @rdname defunct
consolidateTests <- function(...) {
    .Defunct(new="mergeResultsList")
}

#' @export
#' @rdname defunct
consolidateOverlaps <- function(...) {
    .Defunct(new="mergeOverlapsList")
}
