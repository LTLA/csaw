#' @export
#' @importFrom utils browseURL
csawUsersGuide <- function()
# Show the csaw user's guide 
# 
# written by Aaron Lun.
{
    browseURL(sprintf("https://bioconductor.org/packages/%s/csawUsersGuide/vignettes/csaw.pdf", BiocManager::version()))
}
