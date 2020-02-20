#' @export
empiricalFDR <- function(ids, tab, weights=NULL, pval.col=NULL, fc.col=NULL, fc.threshold=0.05, neg.down=TRUE) {
    fc.col <- .parseFCcol(fc.col, tab, multiple=FALSE)
    pval.col <- .getPValCol(pval.col, tab)
    com.out <- .get_two_one_sided_results(tab, pval.col=pval.col, fc.col=fc.col,
        weights=weights, fc.threshold=fc.threshold)

    if (neg.down) { 
        right <- "up"
        wrong <- "down"
    } else {
        right <- "down"
        wrong <- "up"
    }

    right.com <- com.out[[right]]
    wrong.com <- com.out[[wrong]]
    all.down <- sprintf("num.%s.%s", wrong, colnames(tab)[fc.col])
    right.com[,all.down] <- wrong.com[,all.down]

    # Computing the empirical FDR.
    pval.colname <- colnames(tab)[pval.col]
    right.comp <- right.com[,pval.colname]
    o <- order(right.comp)
    right.comp <- right.comp[o]
    empirical <- findInterval(right.comp, sort(wrong.com[,pval.colname]))/seq_along(right.comp)
    
    empirical <- pmin(1, empirical)
    empirical <- rev(cummin(rev(empirical)))
    empirical[o] <- empirical
    right.com$FDR <- empirical

    # Mopping up.
    right.com$direction <- right
    right.com
}

.make_one_sided <- function(tab, pval.col, fc.col) {
    cur.fc <- tab[,fc.col]
    going.up <- cur.fc > 0
    pval <- tab[,pval.col]
    
    # Calculating each set fresh, to avoid numeric 
    # imprecision from repeated "1-" operations
    up.p <- pval/2
    up.p[!going.up] <- 1 - up.p[!going.up]
    down.p <- pval/2
    down.p[going.up] <- 1 - down.p[going.up]

    list(up=up.p, down=down.p)
}
