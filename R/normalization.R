#' SelectGenesFromHistBreaks
#'
#' @param data.counts
#' @param breakss
#' @param n.genes.per.break
#'
#' @return
#' @export
#'
#' @examples
SelectGenesFromHistBreaks <- function(data.counts, breakss=100,
                                    n.genes.per.break=1)
{
    data.counts.mean <- apply(data.counts, 1, mean)
    idxc <- which(data.counts.mean <= summary(data.counts.mean)["3rd Qu."])
    data.counts.mean.filtered <- data.counts.mean[idxc]
    h <- hist(data.counts.mean.filtered, breaks=breakss, plot=FALSE)
    data.counts.mean.sorted <- sort(data.counts.mean.filtered)
    gene.names <- c()
    j<-1
    for(i in 1:length(h$counts)) {
        names <- names(data.counts.mean.sorted[ j : (j+h$counts[i]-1) ])
        if(length(names) >= n.genes.per.break) {
            n.g.names <- names[c(1:n.genes.per.break)]
        } else {
            n.g.names <- names[c(1:length(names))]
        }
        gene.names <- c(gene.names, n.g.names)

        j <- j+h$counts[i]
    }
    return(gene.names)
}

#' EstimateNegativeControlGenesForRUVg
#' @description Extrapolates a list of negative genes for RUVg normalization.
#' It organizes the n.tail.genes genes in histogram sampling
#' n.genes.per.hist.break for each bin.
#' @param de.genes A dataframe of DE genes resulting from DESeq2 or NOISeq.
#' @param n.tail.genes the number of DE genes to take from the tail of the
#' dataframe (default is 2000).
#' @param counts.dataset the dataframe of counts of the genes
#' (requires same rownames as DE genes).
#' @param n.genes.per.hist.break number of genes to sample for each bin.
#' @param threshold the p.adjusted threshold for significative genes.
#'
#' @return the names of the negative genes
#' @export
#'
#' @examples
EstimateNegativeControlGenesForRUV <- function(de.genes, n.tail.genes=2000,
                                            counts.dataset,
                                            n.genes.per.hist.break=1,
                                            threshold=0.05)
{
    de.genes.tail <- TailDeGenesByPAdj(de.genes=de.genes, n.genes=n.tail.genes,
                                    threshold=threshold)
    counts.tail <- counts.dataset[which(rownames(counts.dataset) %in%
                                            de.genes.tail),]
    estimated.genes <- SelectGenesFromHistBreaks(counts.tail,
                                    n.genes.per.break=n.genes.per.hist.break)
    return(estimated.genes)
}


#' RUVgNormalizationFunction
#'
#' @param data.to.normalize
#' @param design.matrix
#' @param desMatColStr
#' @param estimated.gene.names
#' @param k
#' @param isLog
#' @param gene.names.col
#'
#' @return
#' @export
#' @importFrom EDASeq newSeqExpressionSet
#' @importFrom  RUVSeq RUVg
#' @examples
RUVgNormalizationFunction <- function(data.to.normalize,
                                      design.matrix,
                                      desMatColStr,
                                      estimated.gene.names,
                                      k=1, isLog=FALSE)
{

    if( length( which(colnames(design.matrix) == desMatColStr)) == 0 ) {
        stop("Please provide a design matrix with a named ", desMatColStr,
            " column")
    }
    # require("RUVSeq")

    genes <- rownames(design.matrix)

    data.to.normalize <- data.to.normalize[,
                which(colnames(data.to.normalize) %in% colnames(design.matrix))]
    if(!isLog) {
        dataset.set <- EDASeq::newSeqExpressionSet(as.matrix(data.to.normalize),
                                    phenoData=data.frame(
                                    design.matrix[[desMatColStr]],
                                    row.names=
                                    ))
    } else {
        dataset.set <- data.to.normalize
    }


    ruved.set <- RUVSeq::RUVg(x=dataset.set, cIdx=estimated.gene.names,
                              k=k, isLog=isLog)


    #return(as.data.frame(ruved.set@assayData$normalizedCounts))
    return(ruved.set)
}

RUVsNorm <- function(data, samples, nc, k=1, isLog=FALSE)
{
    stop("RUVs needs to be implemented")
}

#' NormalizeData
#'
#' @param data.to.normalize
#' @param norm.type
#' @param design.matrix
#' @param design.matrix.factors.column
#' @param estimated.genes
#' @param is.log
#'
#' @return
#' @export
#' @importFrom edgeR DGEList calcNormFactors estimateCommonDisp
#' estimateTagwiseDisp
#' @importFrom preprocessCore normalize.quantiles
#' @importFrom RUVSeq makeGroups
#' @examples
NormalizeData <- function(data.to.normalize,
    norm.type=c("fqua", "uqua", "tmm", "ruvg", "ruvs"),
    design.matrix=NULL, factors.column=NULL,
    estimated.genes=NULL, is.log=FALSE, ruv_k=1)
{
    ## @ norm.type can be uqua, tmm, fqua or ruvg

    x <- data.to.normalize
    norm.type<- match.arg(norm.type)

    if( all(is.null(design.matrix), is.null(estimated.genes),
        is.null(factors.column)) )
    {
        stop("Please select a design matrix and a list of negative",
             " control genes for RUVg/s normalization")
    }

    switch(norm.type,
        "fqua"=
        {
            x <- as.matrix(x)
            normalized.data <- as.data.frame(normalize.quantiles(x, copy=FALSE))
            colnames(normalized.data) <- colnames(x)
            rownames(normalized.data) <- rownames(x)
        },
        "ruvg"=
        {
            normalized.data <- RUVgNormalizationFunction(
                data.to.normalize=data.to.normalize,
                design.matrix=design.matrix,
                desMatColStr=design.matrix,
                estimated.gene.names=estimated.genes,
                k=ruv_k,
                isLog=is.log)
        },
        "ruvs"=
        {
            samples <- makeGroups(design.matrix[[factors.column]])
            normalized.data <- RUVsNorm(data=data.to.normalize,
                samples=samples, k=ruv_k, nc=estimated.genes, isLog=is.log)
        },
        "uqua"={
            normalized.data = .normedger(x, norm.type)
        },
        "tmm"=
        {
            normalized.data = .normedger(x, norm.type)
        }#,
        # {
        #     warning("No valid normalization selected, returning NULL")
        #     normalized.data=NULL
        # }
    )

    return(normalized.data)
}

.normedger <- function(x, norm.type)
{
    x <- edgeR::DGEList(counts=x)
    norm <- ifelse(norm.type=="tmm", "TMM", "upperquartile")
    x <- edgeR::calcNormFactors(x, method=norm)
    x <- edgeR::estimateCommonDisp(x, verbose=FALSE)
    x <- edgeR::estimateTagwiseDisp(x)
    normalized.data <- as.data.frame(x$pseudo.counts)
    return(normalized.data)
}
#' SortDeGenesByPAdj
#'
#' @param de.genes
#'
#' @return
#' @export
#'
#' @examples
SortDeGenesByPAdj <- function(de.genes)
{
    if("baseMean" %in% colnames(de.genes))
    { ## working on DESeq results
        ordered.de.genes <- de.genes[order(de.genes$padj),]
    } else if("theta" %in% colnames(de.genes)) { ## noiseqbio
        ordered.de.genes <- de.genes[order(de.genes$prob, decreasing=TRUE),]
    } else if("M" %in% colnames(de.genes)) { ## noiseq
        ordered.de.genes <- de.genes[order(de.genes$prob, decreasing=TRUE),]
    } else {
        stop("In SortDeGenesByPAdj not recognized de results table!")
    }

    # ordered.de.genes <- de.genes[order(de.genes$padj),]
    return(ordered.de.genes)
}

#' FilterOutGenesNAPAdj
#'
#' @param de.genes
#'
#' @return
#' @export
#'
#' @examples
FilterOutGenesNAPAdj <- function(de.genes)
{
    if("baseMean" %in% colnames(de.genes)) { ## working on DESeq results
        if(sum(is.na(de.genes$padj))>0)
        {
            de.genes.not.na <- de.genes[-which(is.na(de.genes$padj)),]
        }
        else de.genes.not.na <- de.genes
    } else {
        de.genes.not.na <- de.genes
    }

    # if(sum(is.na(de.genes$padj))>0) de.genes.not.na <- de.genes[-which(is.na(de.genes$padj)),]
    # else de.genes.not.na <- de.genes
    return(de.genes.not.na)
}



#' TailDeGenesByPAdj
#'
#' @param de.genes
#' @param n.genes
#' @param threshold
#'
#' @return
#' @export
#'
#' @examples
TailDeGenesByPAdj <- function(de.genes, n.genes=10, threshold=0.05)
{
    ord.de.genes <- SortDeGenesByPAdj(de.genes=de.genes)
    ord.de.genes.not.na <- FilterOutGenesNAPAdj(ord.de.genes)
    if("baseMean" %in% colnames(de.genes)) { ## working on DESeq results
        ord.not.de.genes <- ord.de.genes.not.na[ ord.de.genes.not.na$padj >=
                                                threshold,]
    } else if("prob" %in% colnames(de.genes)) { ## noiseqbio
        ord.not.de.genes <- ord.de.genes.not.na[ ord.de.genes.not.na$prob <=
                                                threshold,]
    }
    tail.genes <- rownames(tail(ord.not.de.genes, n=n.genes))
    return(tail.genes)
}


