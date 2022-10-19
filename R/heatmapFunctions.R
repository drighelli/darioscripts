
#' plotPHeatmap
#'
#' @param expression.data
#' @param filename
#' @param conversion.map
#' @param width
#' @param height
#' @importFrom pheatmap pheatmap
#' @importFrom heatmaply heatmaply
#' @importFrom clusterExperiment seqPal3
#' @importFrom grDevices colorRampPalette
#' @return
#' @export
#'
#' @examples
plotPHeatmap <- function(expression.data, filename, conversion.map,
                    width=10, height=10)
{
    color.palette = clusterExperiment::seqPal3#c("black", "yellow")
    pal <- colorRampPalette(color.palette)(n = 1000)

    # heatmap_data_scaled <- t(scale(t(as.matrix(expression.data))))

    # pkm7 <- pheatmap(as.matrix(expression.data), color=pal, kmeans_k=5, scale="row", cluster_cols=FALSE)
    # expression.data <- expression.data[order(pkm7$kmeans$cluster),]
    #
    # pheatmap(as.matrix(expression.data), color=pal, scale="row", cluster_rows=FALSE, cluster_cols=FALSE, gaps_row=table(sort(pkm7$kmeans$cluster)))
    # ph1 <- pheatmap(as.matrix(expression.data), color=pal, scale="row",
    #                 cluster_cols=FALSE, filename=paste0(filename,"heatmap.pdf"))

    gene.ordered <- as.data.frame(rownames(expression.data))
    colnames(gene.ordered) <- "ID"
    if(!is.null(conversion.map))
    {
        idx <- match(gene.ordered$ID, conversion.map$ensembl_gene_id)
        gene.ordered$gene_name <- NA
        gene.ordered$gene_name <- conversion.map$external_gene_name[idx]
    }
    expression.data$gene_names <- gene.ordered$gene_name
    heatmaply::heatmaply(expression.data, colors=c("darkred","white","darkblue"), scale="row", file = paste0(filename,"heatmaply_plot.html"))
    # WriteDataFrameAsTsv(gene.ordered, paste0(filename, "_heatmap_hclust_genes"))
}


#' computeGeneMeansOverGroups
#'
#' @param counts
#' @param design
#' @param groupColumn
#'
#' @return
#' @export
#' @importFrom plyr ldply
#' @examples
computeGeneMeansOverGroups <- function(counts, design, groupColumn)
{
    groups <- unique(design[[groupColumn]])
    aa <- plyr::ldply(groups, function(x)
    {
        idx <- which(design[[groupColumn]] %in% x)
        apply(counts[,rownames(design)[idx]], 1, mean)
    })
    rownames(aa) <- groups
    aa <- t(aa)
}


#' savePheatmapPdf
#'
#' @param plot
#' @param filename
#' @param width
#' @param height
#'
#' @return
#' @export
#'
#' @importFrom grid grid.newpage grid.draw
#' @importFrom grDevices pdf dev.off
#'
#' @examples
savePheatmapPdf <- function(plot, filename, width=10, height=10)
{
    x <- plot
    stopifnot(!missing(x))
    stopifnot(!missing(filename))
    grDevices::pdf(filename, width=width, height=height)
    grid::grid.newpage()
    grid::grid.draw(x$gtable)
    grDevices::dev.off()
}


