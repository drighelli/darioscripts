
#' @name PCAPlots
#'
#' @param counts a matrix or data.frame of genes x samples
#' @param colData a data.frame of samples metadata containing the column
#' for shapeBy and colorBy argument
#' @param shapeBy character indicating the column of the colMetadata to shape
#' samples by
#' @param colorBy character indicating the column of the colMetadata to color
#' samples by
#' @param scale logical indicating if to scale the columns of the counts see
#' \link[stats]{prcomp} scale. argument for additional information
#' @param xPCA,yPCA characters indicating which PCs to use for x and y axis in
#' the plot (defaults are PC1 and PC2), only for
#' @param pcas a list of principal components to plot (default is
#' c("PC1", "PC2", "PC3"))
#' @param outFolder character which if not null automatically saves the plot in
#' the indicated path
#' @param size integer indicating the size of the points (default is 2)
#' @param prefix character indicating a prefix to add to the title of the plot
#' @param ellipse logical indicating if to draw ellipses around the detected
#' groups (default is FALSE)
#' @param commleg logical indicating if to use a common legend for all the PCAs,
#' only for the matrix plot (default is TRUE)
#' @param legsize dimension of the legend (default is 14)
#' @param cowplot logical indicating if to use the cowplot theme (default is
#' TRUE)
#' @param plotly logical indicating if to use the plotly theme (default is
#' FALSE)
#' @param returnOnly logical which if TRUE only returns the plot without showing
#' it
#' @param title logical indicating if to include a title or not
#' @importFrom ggpubr ggarrange
#' @return a ggplot2 object

#'
#' @examples
#
# plotPCAsMatrix <- function(se, assayName="counts", shapeBy=NULL, colorBy=NULL,
#                             scale=FALSE,
#                             pcas=c("PC1", "PC2", "PC3"),
#                             prefix=NULL,
#                             ellipse=FALSE, outFolder=NULL,
#                             size=2, commleg=TRUE, cowplot=TRUE,
#                             plotly=FALSE,
#                             returnOnly=FALSE)
# {
#     stopifnot(is(se, "SummarizedExperiment"))
#     cnts <- assays(se)[grep(assayName, assays(se))]
#     plotPCAsMatrix(counts=cnts, colData=colData(se),
#                    shapeBy=shapeBy, colorBy=colorBy,
#                    pcas=pcas,
#                    plotly=plotly, returnOnly=TRUE, title=FALSE,
#                    size=size, scale=scale, prefix=prefix, cowplot=cowplot)
# }

#' @export
plotPCAsMatrix <- function(counts, colData,
                        shapeBy="condition", colorBy="genotype",
                        scale=FALSE,
                        pcas=c("PC1", "PC2", "PC3"),
                        prefix=NULL,
                        ellipse=FALSE, outFolder=NULL,
                        size=2, commleg=TRUE, cowplot=TRUE,
                        legsize=14,
                        plotly=FALSE,
                        returnOnly=FALSE)
{
    strings <- plotStrings(path=outFolder, prefix=prefix, plot="PCA")
    cbn <- t(combn(x=pcas, m=2))
    pca.list <- as.list(apply(cbn, MARGIN=1, function(coord)
    {
        ggp <- plotPCAs(counts=counts,
                        colData=colData, shapeBy=shapeBy, colorBy=colorBy,
                        xPCA=coord[1], yPCA=coord[2],
                        plotly=plotly, returnOnly=TRUE, title=FALSE,
                        size=size, scale=scale, prefix=prefix, cowplot=cowplot)
        ggp
    }))
    ggp <- ggpubr::ggarrange(plotlist=pca.list, common.legend=commleg)
    ggp <- ggpubr::annotate_figure(ggp, top=ggpubr::text_grob(strings$title,
                                color = "black", face = "bold", size = legsize))
    if(!returnOnly) plot(ggp)

    if(!is.null(outFolder))
    {
        if(!is.null(strings$file.name))
        {

            ggplot2::ggsave(filename=strings$file.name, plot=ggp, device="pdf",
                            path=strings$folder, dpi=600)
        }
    }
    return(ggp)
}

# returns TRUE if presence is in colnames(df)
.checkColnm <- function(df, presence)
{
    return(sum(colnames(df) %in% presence) != 0 )
}

#'@export
plotPCAs <- function(counts, colData, shapeBy="condition", colorBy="genotype",
    scale=FALSE, xPCA="PC1", yPCA="PC2", prefix=NULL, ellipse=FALSE,
    outFolder=NULL, size=2, cowplot=TRUE, plotly=FALSE, returnOnly=FALSE,
    title=TRUE)
    ## Add argument to plot colnames with points
{
    stopifnot(!is.null(colnames(counts)))
    if(!is.null(shapeBy)) {
        if(!.checkColnm(colData, shapeBy))
            stop(shapeBy, "column not present in colData")
    }
    if(!is.null(colorBy)) {
        if(!.checkColnm(colData, shapeBy))
            stop(colorBy, "column not present in colData")
    }
    # require("ggfortify")
    # require("plotly")
    strings <- plotStrings(path=outFolder, prefix=prefix, plot="PCA")
    subcounts <- counts[ , which(colnames(counts) %in% rownames(colData))]
    stopifnot(dim(subcounts)[2]!=0)
    PCA <- prcomp(t(subcounts), scale.=scale, center=TRUE)
    propvar <- summary(PCA)$importance[2,,drop=FALSE]
    xproppca <- paste0(as.character( propvar[,xPCA]*100 ),"%")
    yproppca <- paste0(as.character( propvar[,yPCA]*100 ),"%")
    sub.pca <- as.data.frame(PCA$x[,c(xPCA, yPCA)])
    sub.pca$samples <- rownames(sub.pca)
    sub.pca <- sub.pca[order(sub.pca$samples), ]
    colDatao <- colData[order(rownames(colData)), ,
                                    drop=FALSE]
    ncond <- length(unique(colData[[colorBy]]))
    nshape <- length(unique(colData[[shapeBy]]))
    mypalette <- colorRampPalette(RColorBrewer::brewer.pal(8, "Set2"))(ncond)
    if(dim(sub.pca)[1] == dim(colDatao)[1])
    {
        new.sub.pca <- cbind(sub.pca, colDatao)
    } else {
        subColData <- subset(colDatao,rownames(colDatao) %in% rownames(sub.pca))
        new.sub.pca <- cbind(sub.pca, subColData)
    }

    if(ellipse)
    {
        if(.checkColnm(new.sub.pca, colorBy)) ## cambiare per verificare se non sono NULL
        {
            setcolname <- colorBy
        } else {
            setcolname <- shapeBy
        }

        aesStrObj <- ggplot2::aes_string(x=xPCA, y=yPCA, color=setcolname,
                                         shape=shapeBy)
        aesStrObjEll <- ggplot2::aes_string(x=xPCA, y=yPCA, color=shapeBy)

        ggp <- ggplot2::ggplot(new.sub.pca) +
            ggplot2::geom_point(aesStrObj, size=size) +
            ggplot2::scale_shape_manual(values=1:nshape) +
            ggplot2::stat_ellipse(aesStrObjEll) +
            ggplot2::xlab(paste(xPCA, xproppca)) +
            ggplot2::ylab(paste(yPCA, yproppca))
        if(title) ggp <- ggp + ggplot2::ggtitle(strings$title)
        strings$file.name <- paste0(strings$file.name, "_ellipse")
    } else if(!ellipse) {
            aesStrObj <- ggplot2::aes_string(x=xPCA, y=yPCA, color=colorBy,
                                    shape=shapeBy)
            ggp <- ggplot2::ggplot(new.sub.pca) +
                ggplot2::geom_point(aesStrObj, size=size) +
                ggplot2::scale_shape_manual(values=1:nshape) +
                ggplot2::xlab(paste(xPCA, xproppca)) +
                ggplot2::ylab(paste(yPCA, yproppca)) +
                ggplot2::scale_fill_manual(values=mypalette)
            if(title) ggp <- ggp + ggplot2::ggtitle(strings$title)
    }
    if(cowplot) ggp <- ggp + cowplot::theme_cowplot()

    if(!is.null(outFolder))
    {
        if(!is.null(strings$file.name))
        {
            ggplot2::ggsave(filename=strings$file.name, plot=ggp, device="pdf",
                    path=strings$folder, dpi=600)
        }
    }

    if(!returnOnly)
    {
        if(plotly)
        {
            plotly::ggplotly(ggp)
        } else {
            plot(ggp)
        }
    }
    return(ggp)
}



#' PlotCountsAlongTimes
#'
#' @param normalized.counts
#' @param design.matrix
#' @param gene.name
#' @param gene.name.column.name
#' @param show.plot.flag
#' @param plotly.flag
#' @param save.plot
#' @param plot.folder
#' @param prefix.plot
#'
#' @return
#' @export
#'
#' @import plotly
#'
#' @examples
PlotCountsAlongTimes <- function(normalized.counts, design.matrix,
                                 gene.name, gene.name.column.name="gene.names",
                                 show.plot.flag=TRUE, plotly.flag=FALSE,
                                 save.plot=FALSE, plot.folder=NULL,
                                 prefix.plot="Counts Plot")
{

    strings <- GeneratePlotStrings(path=plot.folder,
                                   prefix=prefix.plot,
                                   plot.type="CountsTimePlot")

    sub.normalized.counts <- normalized.counts[,
                which(colnames(normalized.counts) %in% rownames(design.matrix))]
    sub.normalized.counts[,gene.name.column.name] <- normalized.counts[,gene.name.column.name]

    gene.name.norm.counts <- sub.normalized.counts[which(tolower(sub.normalized.counts[,gene.name.column.name]) %in% tolower(gene.name)), ]

    gene.name.norm.counts <- gene.name.norm.counts[, which(colnames(gene.name.norm.counts) %in% rownames(design.matrix))]

    if(dim(gene.name.norm.counts)[1] == 1) {
        processed.counts.df <- ProcessCountDataFrameForPlotCountsAcrossTimes(gene.name.norm.counts, design.matrix, gene.name)
    } else if(dim(gene.name.norm.counts)[1] == 0) {
        stop(gene.name, " gene Not Found!")
    } else if(dim(gene.name.norm.counts)[1] > 1) {
        stop(gene.name, " founded in more than one row!")
    }
    require("plotly")
    ggp <- ggplot(processed.counts.df, mapping=aes(x=Times, y=Counts, color=Conditions, group=Conditions)) + geom_point() + stat_smooth(se=FALSE, method="loess") + scale_y_log10() + ggtitle(paste( strings$title, gene.name, "gene", sep=" "))

    if(save.plot) {
        if(is.null(strings$plot.folder)) {
            stop("Please set a folder where to plot the boxplot!")
        }
        if(!is.null(strings$plot.file.name)){
            SaveGGplot(ggplot.to.save=ggp, plot.folder=strings$plot.folder, plot.file.name=paste(strings$plot.file.name, gene.name, sep="_"), plotly.flag=plotly.flag)
        }

    }

    if(show.plot.flag) {
        if(plotly.flag) {
            ggplotly(ggp)
        } else {
            plot(ggp)
        }
    }

}

