#' @export
#' @importFrom S4Vectors DataFrame
#' @importFrom stats p.adjust
getBestTest <- function(ids, tab, by.pval=TRUE, weights=NULL, pval.col=NULL, fc.col=NULL, fc.threshold=0.05, cpm.col=NULL)
# This uses Holms' method to provide strong control of the FWER within each
# cluster. The idea is that it returns the test with the lowest p-value in the
# cluster. You can then use one test as the representative of the entire cluster,
# which is more specific than Simes (where it's a vague statement of, the DB
# event is somewhere in there). 
#
# written by Aaron Lun
# created 17 April 2014
# last modified 8 January 2017
{
    if (by.pval) {
        minimalTests(ids, tab, weights=weights, pval.col=pval.col, fc.col=fc.col, 
            fc.threshold=fc.threshold, min.sig.n=0, min.sig.prop=0)
    } else {
        if (is.null(cpm.col)) { 
            cpm.col <- "logCPM" 
        }
        if (length(cpm.col)!=1L) {
            stop("absent or multiple logCPM columns are not supported")
        }
        abs<- tab[,cpm.col]

        .general_test_combiner(ids=ids, tab=tab, weights=weights, 
            pval.col=pval.col, fc.col=fc.col, fc.threshold=fc.threshold,
            FUN=function(...) {
                .Call(cxx_compute_cluster_maxed, ..., abs)
            }
        )
    }
}

# You can reduce the conservativeness of the Bonferroni method by using windows that 
# don't overlap. You can also just use combineTests with smaller clusters if you're getting
# lost within each cluster.
#    More complex approaches would attempt to estimate the correlation based on genome-wide
# data. I think the best approach would be to compute the 1st order statistic for p-values
# in each cluster of a given size across the dataset. You can then fit the observed statistics
# robustly (i.e., get rid of very low p-values) to a B(1, n) distribution to determine the
# effective number of independent tests 'n'. You can then apply that to the Bonferroni 
# correction for a cluster of that size.
