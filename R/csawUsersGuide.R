#' @export
#' @importFrom utils browseURL
csawUsersGuide <- function()
# Show the csaw user's guide 
# 
# written by Aaron Lun.
{
    browseURL(sprintf("http://bioconductor.org/%s/workflows/vignettes/csawUsersGuide/inst/doc/csawUsersGuide.pdf", BiocManager::version()))
}
