# Tests the detailRanges() function.
# library(csaw); library(testthat); source("test-detail.R")

library(org.Mm.eg.db)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
X <- detailRanges(orgdb=org.Mm.eg.db, txdb=TxDb.Mmusculus.UCSC.mm10.knownGene, promoter=c(3000, 1000))

test_that("detailRanges works without any supplied regions", {
    expect_true(all(X$type %in% c("E", "P", "G")))
    all.ex <- X[X$type=="E"]
    all.prom <- X[X$type=="P"]
    all.gene <- X[X$type=="G"]

    # Names should be the same.
    ex.names <- sort(unique(names(all.ex)))
    prom.names <- sort(unique(names(all.prom)))
    gene.names <- sort(unique(names(all.gene)))

    expect_identical(ex.names, prom.names)
    expect_identical(gene.names, ex.names)
    expect_false(is.unsorted(names(X)))

    # Symbols should be unique within names.
    by.name <- split(X$symbol, names(X))
    by.name <- lapply(by.name, unique)
    expect_true(all(lengths(by.name)==1L))

    # Exons should be contained within genes.
    olap <- overlapsAny(all.ex, all.gene, type="within")
    expect_true(all(olap))

    # Promoters should match up to the gene start.
    # Not all of them, though, due to alternative TSS.
    # Using BiocGenerics:: as testthat is unhappy during CHECK.
    custom.prom <- suppressWarnings(trim(promoters(all.gene, upstream=3000, downstream=1000)))
    expect_true(all(!is.na(BiocGenerics::match(custom.prom, all.prom))))
})

set.seed(700001)
test_that("detailRanges works with a supplied region", {
    chromos <- seqlengths(X)
    chromos <- chromos[grepl("^chr[0-9XYM]+$", names(chromos))]
    all.win <- generateWindows(chromos*1.8, 1e2, 1000)
    Y <- detailRanges(all.win, orgdb=org.Mm.eg.db, txdb=TxDb.Mmusculus.UCSC.mm10.knownGene, dist=5000)

    all.ex <- X[X$type=="E"]
    all.prom <- X[X$type=="P"]
    all.gene <- X[X$type=="G"]

    # Checking direct overlaps.
    E.olap <- findOverlaps(all.win, all.ex)
    P.olap <- findOverlaps(all.win, all.prom)
    G.olap <- findOverlaps(all.win, all.gene)

    ref <- character(length(all.win))
    all.valid <- unique(c(queryHits(E.olap), queryHits(P.olap), queryHits(G.olap)))
    for (i in all.valid) { 
        curE <- subjectHits(E.olap)[queryHits(E.olap)==i]
        curP <- subjectHits(P.olap)[queryHits(P.olap)==i]
        curG <- subjectHits(G.olap)[queryHits(G.olap)==i]

        curE.id <- names(all.ex)[curE]
        curP.id <- names(all.prom)[curP]
        curG.id <- names(all.gene)[curG]
        universe <- sort(unique(c(curE.id, curP.id, curG.id)))

        has.ex <- universe %in% curE.id
        has.prom <- universe %in% curP.id
        has.gene <- universe %in% curG.id
        tag <- character(length(universe))
        tag[has.prom & !has.ex & !has.gene] <- "P"
        tag[has.prom & has.ex] <- "PE"
        tag[has.prom & !has.ex & has.gene] <- "PI"
        tag[!has.prom & has.ex] <- "E"
        tag[!has.prom & !has.ex & has.gene] <- "I"

        # Need to make sure we're getting the right strand, for multi-locus genes.
        sub.names <- c(curE.id, curP.id, curG.id)
        sub.symb <- c(all.ex$symbol[curE], all.prom$symbol[curP], all.gene$symbol[curG])
        sub.strand <- c(strand(all.ex)[curE], strand(all.prom)[curP], strand(all.gene)[curG])

        m <- match(universe, sub.names)
        cur.symb <- sub.symb[m]
        cur.str <- sub.strand[m]
        ref[i] <- paste(paste0(cur.symb, ":", cur.str, ":", tag), collapse=",")
    }

    expect_identical(ref, Y$overlap)

    # Checking left-flanking overlaps.
    left.flank <- suppressWarnings(trim(flank(all.win, 5000, ignore.strand = TRUE)))
    left.lap <- findOverlaps(left.flank, all.ex, ignore.strand=TRUE)
    left.dist <- start(all.win)[queryHits(left.lap)] - end(all.ex)[subjectHits(left.lap)]
    left.nolap <- left.dist > 0L
    left.lap <- left.lap[left.nolap, ]
    left.dist <- left.dist[left.nolap]

    ref <- character(length(all.win))
    for (i in unique(queryHits(left.lap))) {
        current <- queryHits(left.lap)==i
        curE <- subjectHits(left.lap)[current]
        by.gene <- by(left.dist[current], names(all.ex)[curE], FUN=min)
        mindist <- as.vector(by.gene)
        universe <- names(by.gene)

        # Need to make sure we're getting the right strand, for multi-locus genes.
        sub.names <- names(all.ex)[curE]
        sub.symb <- all.ex$symbol[curE]
        sub.strand <- strand(all.ex)[curE]
        m <- match(universe, sub.names)
        cur.symb <- sub.symb[m]
        cur.str <- sub.strand[m]
        ref[i] <- paste(paste0(cur.symb, ":", cur.str, ":", mindist), collapse=",")
    }

    expect_identical(ref, Y$left)

    # Checking right-flanking overlaps.
    right.flank <- suppressWarnings(trim(flank(all.win, 5000, start=FALSE, ignore.strand = TRUE)))
    right.lap <- findOverlaps(right.flank, all.ex, ignore.strand=TRUE)
    right.dist <- start(all.ex)[subjectHits(right.lap)] - end(all.win)[queryHits(right.lap)]
    right.nolap <- right.dist > 0L
    right.lap <- right.lap[right.nolap, ]
    right.dist <- right.dist[right.nolap]

    ref <- character(length(all.win))
    for (i in unique(queryHits(right.lap))) {
        current <- queryHits(right.lap)==i
        curE <- subjectHits(right.lap)[current]
        curE.id <- names(all.ex)[curE]

        by.gene <- by(right.dist[current], curE.id, FUN=min)
        mindist <- as.vector(by.gene)
        universe <- names(by.gene)

        # Need to make sure we're getting the right strand, for multi-locus genes.
        sub.symb <- all.ex$symbol[curE]
        sub.strand <- strand(all.ex)[curE]
        m <- match(universe, curE.id)
        cur.symb <- sub.symb[m]
        cur.str <- sub.strand[m]
        ref[i] <- paste(paste0(cur.symb, ":", cur.str, ":", mindist), collapse=",")
    }

    expect_identical(ref, Y$right)
})
