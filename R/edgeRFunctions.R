
#' applyEdgeR
#'
#' @param counts
#' @param colData
#' @param factors a character string indicating where to find the contrasts into
#' the design.matrix argument
#' @param wnames a character indicating the colnames of the weights to
#' add to the model.matrix function
#' @param contrasts a character indicating the contrasts to test (i.e. "c1-c2"),
#' more than one can be passed as a list
#' @param edgeroffs offset for \link[edger]{glmQLFit}
#' @param useIntercept
#' @param p.threshold numeric cutoff for adjusted p-values (default is 1)
#' @param is.normalized a logical indicating if the counts are normalized
#' @param verbose logical for printing additional messages
#'
#' @return
#' @export
#'
#' @examples
applyEdgeR <- function(counts, colData, factors=NULL,
                       wnames=NULL, contrasts=NULL, edgeroffs=NULL,
                       useIntercept=FALSE, p.threshold=1, verbose=FALSE)
{
    stopifnot(is.character(factors))
    stopifnot(!is.null(contrasts))
    stopifnot(is.character(contrasts))
    stopifnot(sum(colnames(colData) %in% factors) > 0)
    fc <- factors
    factors <- colData[[factors]]

    if(verbose) message("Setting intercept to: ", useIntercept)

    modelformula <- ifelse(useIntercept, "~ 1 + factors", "~ 0 + factors")
    if(useIntercept) cn <- c("(Intercept)", cn[-1])

    if(!is.null(wnames))
    {
        stopifnot(is.character(wnames))
        stopifnot(sum(colnames(colData) %in% wnames)>0)
        weights <- as.matrix(colData[wnames])
        modelformula <- paste0(modelformula, " + weights")
    }

    design <- model.matrix(as.formula(modelformula))
    colnames(design) <- gsub("factors|weights", "", colnames(design))
    if(length(wnames) == 1 )
    {
        colnames(design) <- c(colnames(design)[1:length(colnames(design))-1],
                              wnames)
    }
    rownames(design) <- factors

    fit <- applyEdgeRQLFit(counts=counts, factors=factors, design=design,
                           edgeroffs=edgeroffs, verbose=verbose)

    resClist <- lapply(contrasts, function(c)
    {

        resC <- applyEdgeRContrast(contrast=c, design=design,
                                    fit=fit, p.threshold=p.threshold,
                                    verbose=verbose)

        cg <- gsub(pattern=" ", replacement="", x=c)
        cs <- strsplit(cg, split="-")[[1]]
        genes <- rownames(resC)
        # if(is.normalized)
        # {
            ctMeans <- computeMeans(counts=counts, design.matrix=colData,
                                    factors.column=fc, contrst=cs,
                                    genes=genes)
            resCMeans <- cbind(ctMeans, resC)
        # } else {
        #     resCMeans <- resC
        # }

        if(p.threshold != 1)
        {
            resCMeans <- resCMeans[(resCMeans$FDR < p.threshold),]
        } else {
            resCMeans <- resCMeans[(resCMeans$FDR <= p.threshold),]
        }


    })
    names(resClist) <- contrasts

    return(resClist)
}


#' attachMeans
#'
#' @param normalized.counts
#' @param design.matrix
#' @param factor.column
#' @param contrast.name
#' @param de.results
#'
#' @return
#' @export
#'
#' @examples
attachMeans <- function(normalized.counts, design.matrix, factor.column,
                        contrast.name, de.results)
{
    cg <- gsub(pattern=" ", replacement="", x=contrast.name)
    cs <- strsplit(cg, split="-")[[1]]
    genes <- rownames(de.results)

    ctMeans <- computeMeans(counts=normalized.counts,
                            design.matrix=design.matrix,
                            factors.column=factor.column,
                            contrst=cs, genes=genes)
    ctMeans <- ctMeans[order(rownames(ctMeans)),]
    de.results <- de.results[order(rownames(de.results)),]
    resCMeans <- cbind(ctMeans, de.results)

    return(resCMeans)

}

#' computeMeans
#'
#' @param counts
#' @param design.matrix
#' @param factors.column
#' @param contrst
#' @param genes
#'
#' @return
#' @export
#'
#' @examples
computeMeans <- function(counts, design.matrix, factors.column, contrst, genes)
{

    design.factors <- design.matrix[, factors.column, drop=FALSE]
    counts <- counts[match(genes, rownames(counts)),]
    contrMeans <- do.call( cbind,
                           lapply(contrst, function(cc)
                           {
                               idx <- which(design.factors[[factors.column]]
                                                                        %in% cc)
                               if(length(idx) == 0) stop("contrast ", cc,
                                                                " not found!")
                               samples <- rownames(design.factors)[idx]
                               subcounts <- counts[,samples]
                               submeans <- apply(X=subcounts, MARGIN=1, mean)
                               return(submeans)
                           })
    )
    colnames(contrMeans) <- contrst
    return(contrMeans)
}

#' applyEdgeRQLFit
#' @description creates a DGEList object, estimate dispersion and applies
#' glmQLFit
#'
#' @param counts matrix of counts of features on rows and samples on columns
#' @param factors grouping factors vector of the same length of the counts column
#' @param design a design matrix typically created with \link[stats]{model.matrix}
#' @param edgeroffs offset for \link[edger]{glmQLFit}
#' @param verbose
#'
#' @return a DGEGLM class object
#' @export
#' @importFrom edgeR DGEList calcNormFactors estimateDisp glmQLFit
#' @examples
#' TBD
applyEdgeRQLFit <- function(counts, factors, design,
        edgeroffs=NULL, verbose=FALSE)
{
    if(verbose) message("Fitting edgeR QL model")
    dgel <- edgeR::DGEList(counts=counts, group=factors)
    edisp <- edgeR::estimateDisp(y=dgel, design=design)
    # fito <- edgeR::glmQLFit(y=edisp, design=design, robust=TRUE)
    fit <- edgeR::glmQLFit(y=edisp$counts, design=design,
        dispersion=edisp$trended.dispersion, robust=TRUE, offset=edgeroffs)
    return(fit)
}

#' applyEdgeRGLMFit
#'
#' @param counts
#' @param factors
#' @param design
#' @param verbose
#' @param is.normalized
#'
#' @return
#' @export
#' @importFrom edgeR DGEList calcNormFactors estimateDisp glmFit
#' @examples
applyEdgeRGLMFit <- function(counts, factors, design,
                            is.normalized=FALSE, method="TMM", verbose=FALSE)
{
    if(verbose) message("Fitting edgeR GLM model")
    dgel <- edgeR::DGEList(counts=counts, group=factors)
    if(!is.normalized) dgel <- edgeR::calcNormFactors(dgel, method=method)
    edisp <- edgeR::estimateDisp(y=dgel, design=design)
    fit <- edgeR::glmFit(edisp, design, robust=TRUE)
    return(fit)
}

#' applyEdgeRContrast
#'
#' @param contrast
#' @param design
#' @param fit
#' @param p.threshold using 1 to be returned all the genes results
#' @param verbose
#'
#' @return
#' @export
#' @importFrom limma makeContrasts
#' @importFrom edgeR glmQLFTest topTags
#' @examples
applyEdgeRContrast <- function(contrast, design, fit, p.threshold=1,
                                verbose=FALSE)
{
    if(verbose) message("contrasting ", contrast)
    contr <- limma::makeContrasts(contrasts=contrast, levels=design)
    qlf <- edgeR::glmQLFTest(fit, contrast=contr)
    res <- edgeR::topTags(qlf, n=Inf, p.value=p.threshold)
    if(dim(res)[1]==0) message("pvalue threshold ", p.threshold,
        " too low, no DEGs found!\nTry with an higher value of p.threshold ...")
    all.gen.res <- res$table
    return(all.gen.res)
}

#' applyEdgeRLRT
#'
#' @param fit
#' @param interaction.matrix
#' @param interaction.term
#' @param p.threshold
#' @param verbose
#'
#' @return
#' @export
#' @importFrom edgeR glmLRT topTags
#' @examples
applyEdgeRLRT <- function(fit, interaction.matrix, interaction.term=NULL,
                          p.threshold=1, verbose=FALSE)
{
    if(is.null(interaction.term))
    {
        warning("interaction.term is missing, using last column of",
                " interaction design!")
        interaction.term <- colnames(interaction.matrix)[[
            ncol(interaction.matrix)]]
    }

    if(verbose) message("testing ", interaction.term)
    qlf <- edgeR::glmLRT(fit, coef=interaction.term)
    res <- edgeR::topTags(qlf, n=Inf, p.value=p.threshold)
    all.gen.res <- res$table
    return(all.gen.res)
}


#' Title
#'
#' @param colname
#' @param colnames
#'
#' @return
#' @export
#'
#' @examples
checkColName <- function(colname, colnames)
{
    if(missing(colname)) stop("Please provide a valid column!")
    if(is.character(colname))
    {
        if(sum(colnames %in% colname)==0)
        {
            stop("Please provide a valid colname!")
        }
    } else {
        if(colname > length(colnames))
        {
            stop("Please provide a valid colname index!")
        }
    }
    return(TRUE)
}


#' constructInteractionMatrix
#'
#' @param design.matrix a design matrix of the experiment where each row
#' describes a sample
#' @param genotype.col the column name/index describing the genotype in
#' the design matrix
#' @param genotype.ref an optional character indicating which factor name to
#' use as genotype reference
#' @param condition.col the column name/index describing the condition in
#' the design matrix
#' @param weight.col optional column name/index describing the weigths in
#' the design matrix (see also weigths parameter)
#' @param weights optional matrix with weights to use in the
#' interaction matrix (see also weigth.col parameter)
#' @param useIntercept
#'
#' @return a design matrix with interaction terms
#' @export
#'
#' @examples
constructInteractionMatrix <- function(design.matrix,
                                       genotype.col, genotype.ref=NULL,
                                       condition.col, weight.col, weights=NULL,
                                       useIntercept=TRUE)
{
    if(missing(design.matrix)) stop("Please provide a valid Design Matrix!")
    if(missing(genotype.col)) stop("Please provide a valid genotype column!")
    checkColName(genotype.col, colnames(design.matrix))
    checkColName(condition.col, colnames(design.matrix))

    if(!is.null(genotype.ref))
    {
        geno <- relevel(design.matrix[[genotype.col]], ref=genotype.ref)
    } else {
        geno <- design.matrix[[genotype.col]]
    }
    cond <- design.matrix[[condition.col]]

    weigths.flag <- TRUE
    if(is.null(weights))
    {
        if(!missing(weight.col))
        {
            checkColName(weight.col, colnames(design.matrix))
            weights <- design.matrix[[weight.col]]
        } else {
            weigths.flag <- FALSE
        }
    }

    if(weigths.flag)
    {
        if(useIntercept)
        {
            interactionFormula <- "~geno*cond+weights"
        } else {
            interactionFormula <- "~0+geno*cond+weights"
        }
    } else {
        if(useIntercept)
        {
            interactionFormula <- "~geno*cond"
        } else {
            interactionFormula <- "~0+geno*cond"
        }
    }

    interactionMatrix <- model.matrix(as.formula(interactionFormula))



    return(interactionMatrix)
}
