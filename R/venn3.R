#' Venn3de
#' @description TBW
#'
#' @param x a list of elements
#' @param y a list of elements
#' @param z a list of elements
#' @param label1 the label for the x list
#' @param label2 the label for the y list
#' @param label3 the label for the z list
#' @param title an optional title for the Venn plot
#' @param intersection.flag a boolean for saving the lists of intersected
#' elements
#' @param intersection.exclusion.flag a boolean for saving the lists of unique
#' sets lists
#' @param plot.dir the directory to save the intersections
#' @param enrich.lists.flag If to perform the functional enrichment for the
#' produced lists
#' @param conversion.map a data.frame with gene ids conversions
#' (i.e. Ensembl_id, Gene_name)
#' @param plot.heatmap flag indicating if plotting a heatmap for the main
#' intersection
#' @param expression.data.list a list of dataframe/matrix with rownames with
#' expression values ordered for x,y,z
#' @param expression.colname the name of the column within the expression values
#'
#' @return NULL
#' @export
#'
#' @examples
#' TBD
Venn3de <- function(x, y, z, label1="x", label2="y", label3="z",
                    title="Venn Diagram",
                    fill_color=c("darkcyan", "darkcyan", "darkcyan"),
                    fill_alpha = 0.4,
                    stroke_size=0.5, stroke_alpha=0.5, stroke_color="black",
                    set_name_size=6,
                    intersection.flag=TRUE,
                    intersection.exclusion.flag=FALSE,
                    plot.dir=NULL,
                    enrich.lists.flag=FALSE,
                    conversion.map=NULL,
                    plot.heatmap=FALSE,
                    expression.data.list=NULL,
                    expression.colname="logFC")
{

    a=x
    b=y
    c=z
    plot.combined.string <- paste(label1, "_", label2, "_", label3, sep="")
    if(is.null(plot.dir))
    {
        stop("Please provide a directory where to save VENN results!")
    }

    # Lists <- list(a, b, c)  # put the word vectors into a list to supply lapply
    # Lists <- lapply(Lists, function(x) as.character(unlist(x)))
    # items <- sort(unique(unlist(Lists)))   #put in alphabetical order
    # MAT <- matrix(rep(0, length(items)*length(Lists)), ncol=3)  #make a matrix of 0's
    # names <- c(label1,label2,label3)
    # colnames(MAT) <- names
    # rownames(MAT) <- items
    # lapply(seq_along(Lists), function(i)
    # {   #fill the matrix
    #     MAT[items %in% Lists[[i]], i] <<- table(Lists[[i]])
    # })

    outputName <- UpdateFilename(filename="VennDiagram",
                                 label1, label2, label3, extension="pdf")
    out.path.name <- UpdateFolderPath(plot.dir, "Venn3")
    prefix=""

    file.path.name <- file.path(out.path.name, outputName)
    # limma::vennDiagram(MAT, circle.col=c("red", "green", "yellow"), main=title)
    l <- list(x, y, z)
    names(l) <- c(label1, label2, label3)
    ggp <- ggvenn::ggvenn(data=l,
                          fill_color = fill_color, fill_alpha = fill_alpha,
                          stroke_size=stroke_size, stroke_alpha=stroke_alpha,
                          stroke_color=stroke_color, set_name_size=set_name_size,
                          text_size=3)
    if(!is.null(title)) ggp <- ggp + ggtitle(title)
    ggplot2::ggsave(filename=file.path.name, width=297, height=210, units="mm")
    intersections <- list()
    if(intersection.flag)
    {
        if(is.null(plot.dir))
        {
            stop("Please provide a directory where to save VENN results!")
        }

        abc <- intersect(intersect(a, b), c)
        intersections <- list("XYZ"=abc)
        expression.data=NULL
        if(plot.heatmap && (!is.null(abc)))
        {
            expressions.list <- lapply(expression.data.list, function(ed)
            {
                ed <- ed[order(rownames(ed)),]
                ed[which(rownames(ed) %in% abc),expression.colname, drop=FALSE]
            })
            expression.data <- rlist::list.cbind(expressions.list)
            colnames(expression.data) <- c(label1, label2, label3)
        }


        SaveInteserctionsList(gene.list=abc,
                        conversion.map=conversion.map,
                        root.dir=out.path.name, prefix=prefix,
                        labels.list=c(label1, paste0("AND_", label2),
                                    paste0("AND_", label3)),
                        enrich.lists.flag=enrich.lists.flag,
                        heatmap.flag=plot.heatmap,
                        expression.data=expression.data)

        ab <- setdiff(intersect(a, b), abc)
        SaveInteserctionsList(gene.list=ab, conversion.map=conversion.map,
                            root.dir=out.path.name, prefix=prefix,
                            labels.list=c(label1, paste0("AND_",label2)),
                            enrich.lists.flag=enrich.lists.flag)

        bc <- setdiff(intersect(b, c), abc)
        SaveInteserctionsList(gene.list=bc, conversion.map=conversion.map,
                            root.dir=out.path.name, prefix=prefix,
                            labels.list=c(label2, paste0("AND_", label3)),
                            enrich.lists.flag=enrich.lists.flag)

        ac <- setdiff(intersect(a, c), abc)
        SaveInteserctionsList(gene.list=ac, conversion.map=conversion.map,
                            root.dir=out.path.name, prefix=prefix,
                            labels.list=c(label1, paste0("AND_", label3)),
                            enrich.lists.flag=enrich.lists.flag)
        intrs <- list("XY"=ab, "YZ"=bc, "XZ"=ac)
        intersections <- c(intersections, intrs)

    }

    if(intersection.exclusion.flag) {
        if(is.null(plot.dir)) {
            stop("Please provide a directory where to save VENN results!")
        }

        a.not.b <-  setdiff(a, b)

        b.not.a <-  setdiff(b, a)

        c.not.a <-  setdiff(c, a)

        a.not.bc <- setdiff(a.not.b, c)
        SaveInteserctionsList(gene.list=a.not.bc, conversion.map=conversion.map,
                            root.dir=out.path.name, prefix=prefix,
                            labels.list=c(label1, paste0("_NOT_", label2), paste0("_NOT_", label3)),
                            enrich.lists.flag=enrich.lists.flag)

        b.not.ac <- setdiff(b.not.a, c)
        SaveInteserctionsList(gene.list=b.not.ac, conversion.map=conversion.map,
                            root.dir=out.path.name, prefix=prefix,
                            labels.list=c(label2, paste0("_NOT_", label1), paste0("_NOT_", label3)),
                            enrich.lists.flag=enrich.lists.flag)

        c.not.ab <- setdiff(c.not.a, b)
        SaveInteserctionsList(gene.list=c.not.ab, conversion.map=conversion.map,
                            root.dir=out.path.name, prefix=prefix,
                            labels.list=c(label3, paste0("_NOT_", label1), paste0("_NOT_", label2)),
                            enrich.lists.flag=enrich.lists.flag)
        intrs <- list("XnotYZ"=a.not.bc, "YnotXZ"=b.not.ac, "ZnotXY"=c.not.ab)
        intersections <- c(intersections, intrs)
    }
    intersections <- list("ggp"=ggp, intersections)
    return(intersections)
}


#' SaveInteserctionsList
#'
#' @param gene.list
#' @param conversion.map
#' @param root.dir
#' @param prefix
#' @param labels.list
#' @param enrich.lists.flag
#' @param heatmap.flag
#' @param expression.data
#'
#' @return
#' @export
#'
#' @examples
SaveInteserctionsList <- function(gene.list, conversion.map,
                            root.dir, prefix, labels.list,
                            enrich.lists.flag=FALSE,
                            heatmap.flag=FALSE,
                            expression.data=NULL)
{
    if(length(gene.list) == 0 )
    {
        warning("skipping ", labels.list, " it's empty")
        return(NULL)
    }
    for(lbl in labels.list) {
        prefix <- UpdatePrefix(prefix, lbl)
    }

    out.dir <- UpdateFolderPath(root.dir, prefix)
    filename <- UpdateFilename(prefix, "genes")

    if(!is.null(conversion.map)) {
        gene.list.df <- CreateConvertedDataframe(gene.list, conversion.map)
    } else {
        gene.list.df <- as.data.frame(gene.list)
    }

    WriteDataFrameAsTsv(data.frame.to.save=gene.list.df,
                    file.name.path=file.path(out.dir, filename),
                    col.names=TRUE, row.names=FALSE)

    if(heatmap.flag)
    {
        # stopifnot(!is.null(expression.data))
        if(length(gene.list)>2)
            plotPHeatmap(expression.data, filename=file.path(out.dir, filename),
                    conversion.map=conversion.map)

    }

    if(enrich.lists.flag)
    {
        enrichSuitedSingleList(de.gene.list=gene.list,
                        functional.folder=file.path(out.dir, "functional"),
                        filename=filename)
    }

}

#' Create a Converted DataFrame from Ensembl IDs
#'
#' This function takes a list of Ensembl gene IDs and a conversion map (dataframe)
#' containing Ensembl IDs and corresponding gene symbols. It returns a dataframe
#' mapping Ensembl IDs to gene symbols.
#'
#' @param gene_list A character vector containing Ensembl gene IDs.
#' @param conversion_map A dataframe with two columns:
#'   - `"EnsemblID"`: Ensembl gene IDs.
#'   - `"GeneSymbol"`: Corresponding gene symbols.
#'
#' @return A dataframe with two columns:
#'   - `"EnsemblID"`: The input Ensembl IDs.
#'   - `"GeneSymbol"`: The corresponding gene symbols (NA if not found).
#'
#' @examples
#' gene_list <- c("ENSG000001", "ENSG000002", "ENSG000003")
#' conversion_map <- data.frame(
#'   EnsemblID = c("ENSG000001", "ENSG000002"),
#'   GeneSymbol = c("GeneA", "GeneB"),
#'   stringsAsFactors = FALSE
#' )
#'
#' df <- CreateConvertedDataframe(gene_list, conversion_map)
#' print(df)
#'
#' @export
CreateConvertedDataframe <- function(gene_list, conversion_map) {
    # Check if inputs are valid
    if (!is.vector(gene_list) || !is.character(gene_list)) {
        stop("gene_list must be a character vector containing Ensembl IDs.")
    }

    if (!is.data.frame(conversion_map) || !all(c("EnsemblID", "GeneSymbol") %in% colnames(conversion_map))) {
        stop("conversion_map must be a dataframe with columns 'EnsemblID' and 'GeneSymbol'.")
    }

    # Ensure that EnsemblID is unique in the conversion map
    if (any(duplicated(conversion_map$EnsemblID))) {
        stop("conversion_map contains duplicate Ensembl IDs.")
    }

    # Merge gene_list with conversion_map to get GeneSymbol
    df <- merge(
        data.frame(EnsemblID = gene_list, stringsAsFactors = FALSE),
        conversion_map,
        by = "EnsemblID",
        all.x = TRUE
    )

    # Ensure correct column order
    df <- df[, c("EnsemblID", "GeneSymbol")]

    return(df)
}
