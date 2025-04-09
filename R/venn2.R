#' Title
#'
#' @param x
#' @param y
#' @param label1
#' @param label2
#' @param title
#' @param plot.dir
#' @param conversion.map
#' @param intersection.flag
#' @param enrich.lists.flag
#' @param prefix
#' @param plot.heatmap
#'
#' @return
#' @export
#'
#' @examples
Venn2de <- function(x, y, label1, label2, title=NULL, plot.dir,
                    conversion.map=NULL, intersection.flag=FALSE,
                    enrich.lists.flag=FALSE,
                    fill_color = c("darkcyan", "darkcyan"), fill_alpha = 0.4,
                    stroke_size=0.5, stroke_alpha=0.5, stroke_color="black",
                    set_name_size=6,
                    prefix="", plot.heatmap=FALSE, expression.data=NULL)
{


    # out.path <- UpdateFolderPath(plot.dir, "venn2")
    outputName <- UpdateFilename(filename="VennDiagram",
                                 label1, label2, extension="pdf")
    out.path.name <- file.path(plot.dir, "venn2")
    file.path.name <- file.path(out.path.name, outputName)
    dir.create(out.path.name, recursive=TRUE)

    a15 <- x
    b15 <- y

    c15 <- intersect(a15, b15)     #common gene names

    ab <- setdiff(a15, b15)

    ba <- setdiff(b15, a15)


    Lists <- list(a15, b15)  #put the word vectors into a list to supply lapply
    Lists <- lapply(Lists, function(x) as.character(unlist(x)))
    items <- sort(unique(unlist(Lists)))   #put in alphabetical order

    SaveInteserctionsList(gene.list=c15,
                          conversion.map=conversion.map,
                          root.dir=out.path.name, prefix=prefix,
                          labels.list=c(label1, paste0("AND_", label2)),
                          enrich.lists.flag=enrich.lists.flag,
                          heatmap.flag=plot.heatmap,
                          expression.data=expression.data)
    if(intersection.flag)
    {
        if(is.null(plot.dir))
        {
            stop("Please provide a directory where to save VENN results!")
        }

        xy <- intersect(x, y)

        expression.data=NULL


        SaveInteserctionsList(gene.list=xy,
                              conversion.map=conversion.map,
                              root.dir=out.path.name, prefix=prefix,
                              labels.list=c(label1, paste0("AND_", label2)),
                              enrich.lists.flag=enrich.lists.flag,
                              heatmap.flag=plot.heatmap,
                              expression.data=expression.data)

        xx <- setdiff(x, xy)
        SaveInteserctionsList(gene.list=xx, conversion.map=conversion.map,
                              root.dir=out.path.name, prefix=prefix,
                              labels.list=label1,
                              enrich.lists.flag=enrich.lists.flag)

        yy <- setdiff(y, xy)
        SaveInteserctionsList(gene.list=yy, conversion.map=conversion.map,
                              root.dir=out.path.name, prefix=prefix,
                              labels.list=label2,
                              enrich.lists.flag=enrich.lists.flag)
    }

    # limma::vennDiagram(MAT, circle.col= c("red","green"), main=title)
    l <- list(x, y)
    names(l) <- c(label1, label2)
    ggp <- ggvenn::ggvenn(data=l,
                          fill_color = fill_color, fill_alpha = fill_alpha,
                          stroke_size=stroke_size, stroke_alpha=stroke_alpha,
                          stroke_color=stroke_color, set_name_size=set_name_size)
    if(!is.null(title)) ggp <- ggp + ggtitle(title)
    ggplot2::ggsave(filename=file.path.name, width=297, height=210, units="mm")
    return(list(p=ggp, int=c15, XnoY=ab, YnoX=ba))
}
