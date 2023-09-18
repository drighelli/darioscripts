#'
#' @param se
#' @param condCol
#' @param assay
#' @param method
#' @param isNormalized
#' @param cv.percentage
#' @param cpm.cutoff
#' @param seq.depth
#' @param lib.size
#' @param min.count
#' @param min.total.count
#' @param large.n
#' @param min.prop
#' @param verbose
#'
#' @return
#' @export
#' @import SummarizedExperiment
#' @importFrom NOISeq filtered.data
#' @importFrom methods is
#' @importFrom edgeR filterByExpr
#' @examples
filterLowCounts <- function(se,
    condCol=NULL, assay="counts",
    method=c("cpm", "wilc", "prop", "edger"),
    ##NOIseq
    isNormalized=FALSE, cv.percentage=100, cpm.cutoff=1,
    seq.depth=NULL,
    ##edgeR
    lib.size=NULL, min.count=10, min.total.count=15, large.n=10,
    min.prop = 0.7,
    verbose=TRUE)
{
    stopifnot(methods::is(se, "SummarizedExperiment"))
    method <- match.arg(method)
    stopifnot(is.logical(isNormalized))
    stopifnot(!is.null(condCol))
    stopifnot((condCol %in% colnames(colData(se))))

    conditions <- colData(se)[, condCol]
    cnts <- assays(se)[[assay]]
    if(method == "Proportion")
    {
        if(isNormalized) {
            ## calcolare seq.depth
            if(is.null(seq.depth))
                stop("Proportion test cannot be performed on normalized counts",
                     " without sequencing depth!\nYou need column totals ",
                     " before normalizing the data.")
        } else {
            seq.depth <- NULL
        }
    }

    switch( method,
            cpm = {
                methodn <- 1
                if(verbose) message("Using CPM test")
            },
            wilc = {
                methodn <- 2
                if(verbose) message("Using Wilcoxon test")
            },
            prop = {
                methodn <- 3
                if(verbose) message("Using Proportion test")
            },
            edger = {
                methodn <- 4
                if(verbose) message("Using edgeR")
            }
    )

    if(verbose) message("Features dimension before normalization: ",
                        dim(cnts)[1])
    if(methodn != 4)
    {
        invisible(capture.output(
        fcounts <- NOISeq::filtered.data(dataset=cnts,
            factor=conditions, norm=isNormalized, depth=seq.depth,
            method=methodn, cv.cutoff=cv.percentage, cpm=cpm.cutoff,
            p.adj="BH")
        ))
    } else {
        keepg <- edgeR::filterByExpr(y=cnts, design=colData(se),
            group=conditions, lib.size=lib.size, min.count=min.count,
            min.total.count=min.total.count, large.n=large.n, min.prop=min.prop)
        fcounts <- cnts[keepg,]

    }
    if(verbose) message("Features dimension after normalization: ", dim(fcounts)[1])
    rowData(se)[] #### use method name
    rowData(se)$kept <- FALSE
    rowData(se)$kept[which(rownames(se) %in% rownames(fcounts))] <- TRUE
    return(se)
}
