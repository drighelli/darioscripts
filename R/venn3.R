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
#' @param intersection.flag a boolean for saving the lists of intersected elements
#' @param intersection.exclusion.flag a boolean for saving the lists of unique sets lists
#' @param plot.dir the directory to save the intersections
#' @param enrich.lists.flag If to perform the functional enrichment for the produced lists
#' @param conversion.map a data.frame with gene ids conversions (i.e. Ensembl_id, Gene_name)
#'
#' @return NULL
#' @export
#'
#' @examples
#' TBD
Venn3de <- function(x, y, z, label1="x", label2="y", label3="z",
                    title="Venn Diagram",
                    intersection.flag=TRUE,
                    intersection.exclusion.flag=FALSE,
                    plot.dir=NULL,
                    enrich.lists.flag=FALSE,
                    conversion.map=NULL)
{

    a=x
    b=y
    c=z
    plot.combined.string <- paste(label1, "_", label2, "_", label3, sep="")
    # a=row.names(x)
    # b=row.names(y)
    # c=row.names(z)
    if(is.null(plot.dir))
    {
        stop("Please provide a directory where to save VENN results!")
    }

    Lists <- list(a, b, c)  # put the word vectors into a list to supply lapply
    Lists <- lapply(Lists, function(x) as.character(unlist(x)))
    items <- sort(unique(unlist(Lists)))   #put in alphabetical order
    MAT <- matrix(rep(0, length(items)*length(Lists)), ncol=3)  #make a matrix of 0's
    names <- c(label1,label2,label3)
    colnames(MAT) <- names
    rownames(MAT) <- items
    lapply(seq_along(Lists), function(i)
    {   #fill the matrix
        MAT[items %in% Lists[[i]], i] <<- table(Lists[[i]])
    })

    outputName <- UpdateFilename(filename="VennDiagram",
                                 label1, label2, label3, extension="pdf")
    out.path.name <- UpdateFolderPath(plot.dir, "Venn3")
    prefix="venn3"

    file.path.name <- file.path(out.path.name, outputName)
    # dev.new()
    # require("limma")
    # limma::vennDiagram(MAT, circle.col= c("red","green","yellow"), main=title)
    #
    # dev.print(device=pdf, file=file.path.name, width=10, height=10)
    # graphics.off()
    limma::vennDiagram(MAT, circle.col=c("red", "green", "yellow"), main=title)

    if(intersection.flag)
    {
        if(is.null(plot.dir))
        {
            stop("Please provide a directory where to save VENN results!")
        }

        abc <- intersect(intersect(a, b), c)
        SaveInteserctionsList(gene.list=abc, conversion.map=conversion.map, root.dir=out.path.name, prefix=prefix, labels.list=c(label1, paste0("AND_", label2), paste0("AND_", label3)), enrich.lists.flag=enrich.lists.flag)

        ab <- setdiff(intersect(a, b), abc)
        SaveInteserctionsList(gene.list=ab, conversion.map=conversion.map, root.dir=out.path.name, prefix=prefix, labels.list=c(label1, paste0("AND_",label2)), enrich.lists.flag=enrich.lists.flag)

        bc <- setdiff(intersect(b, c), abc)
        SaveInteserctionsList(gene.list=bc, conversion.map=conversion.map, root.dir=out.path.name, prefix=prefix, labels.list=c(label2, paste0("AND_", label3)), enrich.lists.flag=enrich.lists.flag)

        ac <- setdiff(intersect(a, c), abc)
        SaveInteserctionsList(gene.list=ac, conversion.map=conversion.map, root.dir=out.path.name, prefix=prefix, labels.list=c(label1, paste0("AND_", label3)), enrich.lists.flag=enrich.lists.flag)

        # union.of.intersections <- union(union(ab, bc), union(ac, abc))
        # SaveInteserctionsList(gene.list=union.of.intersections, conversion.map=conversion.map, root.dir=out.path.name, prefix=paste0(prefix, "_union_of_intersections"), labels.list=c(label1, label2, label3), enrich.lists.flag=enrich.lists.flag)

        # union.of.lugs <- union(bc, union(ac, abc))
        # SaveInteserctionsList(gene.list=union.of.lugs, conversion.map=conversion.map, root.dir=out.path.name, prefix=paste0(prefix, "_union_of_LUGs"), labels.list=c(label1, label2, label3), enrich.lists.flag=enrich.lists.flag)

    }

    if(intersection.exclusion.flag) {
        if(is.null(plot.dir)) {
            stop("Please provide a directory where to save VENN results!")
        }
        res.path <- file.path(plot.dir, "Venn3", paste0(plot.combined.string,"_gene_lists"))
        dir.create(res.path, recursive=TRUE, showWarnings=FALSE)

        a.not.b <-  setdiff(a, b)
        SaveInteserctionsList(gene.list=a.not.b, conversion.map=conversion.map, root.dir=out.path.name, prefix=prefix,
                              labels.list=c(label1, paste0("_NOT_",label2), enrich.lists.flag=enrich.lists.flag))

        b.not.a <-  setdiff(b, a)
        SaveInteserctionsList(gene.list=b.not.a, conversion.map=conversion.map, root.dir=out.path.name, prefix=prefix,
                              labels.list=c(label2, paste0("_NOT_",label1)), enrich.lists.flag=enrich.lists.flag)

        c.not.b <-  setdiff(c, b)
        SaveInteserctionsList(gene.list=c.not.b, conversion.map=conversion.map, root.dir=out.path.name, prefix=prefix,
                              labels.list=c(label3, paste0("_NOT_",label2)), enrich.lists.flag=enrich.lists.flag)

        b.not.c <-  setdiff(b, c)
        SaveInteserctionsList(gene.list=b.not.c, conversion.map=conversion.map, root.dir=out.path.name, prefix=prefix,
                              labels.list=c(label2, paste0("_NOT_", label3)), enrich.lists.flag=enrich.lists.flag)

        a.not.c <-  setdiff(a, c)
        SaveInteserctionsList(gene.list=a.not.c, conversion.map=conversion.map, root.dir=out.path.name, prefix=prefix,
                              labels.list=c(label1, paste0("_NOT_", label3)), enrich.lists.flag=enrich.lists.flag)

        c.not.a <-  setdiff(c, a)
        SaveInteserctionsList(gene.list=c.not.a, conversion.map=conversion.map, root.dir=out.path.name, prefix=prefix,
                              labels.list=c(label3, paste0("_NOT_", label1)), enrich.lists.flag=enrich.lists.flag)

        a.not.bc <- setdiff(a.not.b, c)
        SaveInteserctionsList(gene.list=a.not.bc, conversion.map=conversion.map, root.dir=out.path.name, prefix=prefix,
                              labels.list=c(label1, paste0("_NOT_", label2), paste0("_NOT_", label3)), enrich.lists.flag=enrich.lists.flag)

        b.not.ac <- setdiff(b.not.a, c)
        SaveInteserctionsList(gene.list=b.not.ac, conversion.map=conversion.map, root.dir=out.path.name, prefix=prefix,
                              labels.list=c(label2, paste0("_NOT_", label1), paste0("_NOT_", label3)), enrich.lists.flag=enrich.lists.flag)

        c.not.ab <- setdiff(c.not.a, b)
        SaveInteserctionsList(gene.list=c.not.ab, conversion.map=conversion.map, root.dir=out.path.name, prefix=prefix,
                              labels.list=c(label3, paste0("_NOT_", label1), paste0("_NOT_", label2)), enrich.lists.flag=enrich.lists.flag)

    }
}


SaveInteserctionsList <- function(gene.list, conversion.map,
                            root.dir, prefix, labels.list,
                            enrich.lists.flag=FALSE)
{

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

    if(enrich.lists.flag)
    {
        enrichSuitedSingleList(de.gene.list=gene.list,
                        functional.folder=file.path(out.dir, "functional"),
                        filename=filename)
    }

}
