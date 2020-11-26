#' @export
#' @importFrom utils browseURL
csawUsersGuide <- function()
# Show the csaw user's guide 
# 
# written by Aaron Lun.
{
    thing <- BiocStyle::Biocbook("csawBook")
    target <- sub(".*\\((.*)\\).*", "\\1", thing)
    browseURL(target)
}
