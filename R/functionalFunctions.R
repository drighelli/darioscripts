enrichKEGGFunction <- function(gene.list.to.enrich, gene.type="ENSEMBL",
                               organism.db=c("org.Rn.eg.db", "org.Mm.eg.db"),
                               organism.code=c("rno", "mmu", "hsa"),
                               threshold=0.05, functional.folder, filename)
{
    # require(assign("organism.db", organism.db))
    require("clusterProfiler")

    kegg.folder <- UpdateFolderPath(functional.folder, "clusterProfiler", "KEGG")
    filename <- UpdateFilename(filename, "KEGG_PATHWAYS")
    sign.genes.entrez <- bitr(gene.list.to.enrich, fromType=gene.type, toType=c("ENTREZID", "SYMBOL"), OrgDb=organism.db)
    # total.genes.entrez <- bitr(total.gene.list, fromType=gene.type, toType=c("ENTREZID", "SYMBOL"), OrgDb=organism.db)
    kegg.results <- enrichKEGG(gene=sign.genes.entrez$ENTREZID,
                               organism=organism.code,
                               # universe=total.genes.entrez$ENTREZID,
                               pvalueCutoff=threshold)
    ksdf <- as.data.frame(kegg.results)
    ksdf$geneID <- gsub("/", "; ", ksdf$geneID)
    WriteDataFrameAsTsv(data.frame.to.save=ksdf, file.name.path=file.path(kegg.folder, filename))
    GenerateAndSaveNetwork(kegg.results, kegg.folder, filename)
    return(ksdf)
}


enrichGOFunction <- function(gene.list.to.enrich, gene.type="ENSEMBL", organism.db=c("org.Rn.eg.db", "org.Mm.eg.db"), organism.code=c("rno", "mmu"), threshold=0.05, functional.folder, filename, ontology =c("CC", "MF", "BP")) {
    # require(organism.db)
    require("clusterProfiler")

    go.folder <- UpdateFolderPath(functional.folder, "clusterProfiler", "GO")
    filename <- UpdateFilename(filename, paste0("clusterProfiler_GO_",ontology))
    sign.genes.entrez <- bitr(gene.list.to.enrich, fromType=gene.type, toType=c("ENTREZID", "SYMBOL"), OrgDb=organism.db)

    ego <- enrichGO(gene=sign.genes.entrez$ENTREZID,
                    keyType="ENTREZID",
                    # universe     =rownames(thoracic.de.uqua.notnorm2.03d$de.not.na),
                    OrgDb        =organism.db,
                    ont          =ontology,
                    pAdjustMethod="BH",
                    pvalueCutoff =threshold,
                    qvalueCutoff =threshold,
                    readable     =TRUE)

    ksdf <- as.data.frame(ego)
    ksdf$geneID <- gsub("/", "; ", ksdf$geneID)
    WriteDataFrameAsTsv(data.frame.to.save=ksdf, file.name.path=file.path(go.folder, filename))
    GenerateAndSaveHierarchicalGO(go.results=ego, go.folder=go.folder, filename=filename, ontology=ontology)
    return(ksdf)
}



GenerateAndSaveNetwork <- function(kegg.results, kegg.folder, filename) {
    require("clusterProfiler")

    plot.functional.folder <- UpdateFolderPath(kegg.folder, "KEGG_MAP")
    filename <- UpdateFilename(filename, "kegg_map.pdf")
    filename <- file.path(plot.functional.folder, filename)

    pdf(filename)
    enrichplot::emapplot(kegg.results)
    dev.off()
    message(filename, " Saved on disk!")
}



GenerateAndSaveHierarchicalGO <- function(go.results, go.folder, filename, ontology) {
    require("clusterProfiler")
    plot.functional.folder <- UpdateFolderPath(go.folder, paste0("GO_", ontology, "_Hierarchy"))
    filename <- UpdateFilename(filename, "_Hierarchy.pdf")
    filename <- file.path(plot.functional.folder, filename)

    pdf(filename)
    # tryCatch({
    plotGOgraph(go.results)
    # message("problem with plotgograph")
    # })
    dev.off()
    message(filename, " Saved on disk!")
}

########################################################################## GProfiler

# exclude.electr.flag <- FALSE, go.label <- "GO:BP", include.graph.flag <- TRUE, path.label <- "KEGG", organism.name="rnorvegicus", significant.flag=TRUE, min.isect.size=0, correction.method="fdr"

#' Title
#'
#' @param gene.list.to.enrich
#' @param organism.name
#' @param exclude.electr.flag
#' @param path.label
#' @param significant.flag
#' @param correction.method
#' @param functional.folder
#' @param filename
#'
#' @return
#' @export
#'
#' @examples
enrichPathwayGProfiler <- function(gene.list.to.enrich,
                                   organism.name=c("rnorvegicus", "mmusculus", "hsapiens"),
                                   exclude.electr.flag=FALSE,
                                   path.label=c("KEGG", "REAC"),
                                   significant.flag=FALSE,
                                   # min.isect.size=0,
                                   correction.method="fdr",
                                   functional.folder, filename) {
    # require(assign("organism.db", organism.db))
    require("gprofiler2")

    path.folder <- UpdateFolderPath(functional.folder, "gProfiler", path.label)
    filename <- UpdateFilename(filename, path.label, "_PATHWAYS")

    ksdf <- gost(query=gene.list.to.enrich,
                 organism=organism.name,
                 significant=significant.flag,
                 # min_isect_size=min.isect.size,
                 correction_method=correction.method,
                 sources=path.label,
                 exclude_iea=exclude.electr.flag, evcodes=TRUE)
    ksdf <- ksdf$result
    ksdf$intersection <- .elaborategProfilerListToPrint(ksdf$intersection)
    ksdf$evidence_codes <- .elaborategProfilerListToPrint(ksdf$evidence_codes)
    ksdf$parents <- .elaborategProfilerListToPrint(ksdf$parents)
    ksdf <- ksdf[order(ksdf$p_value), ]

    # print(head(ksdf))

    if(dim(ksdf)[1]!=0)
    {
        if(significant.flag)
        {
            WriteDataFrameAsTsv(data.frame.to.save=ksdf, file.name.path=file.path(path.folder, paste0(filename, "_significant")))
        } else {
            WriteDataFrameAsTsv(data.frame.to.save=ksdf[which(ksdf$significant==TRUE),], file.name.path=file.path(path.folder, paste0(filename, "_significant")))
            WriteDataFrameAsTsv(data.frame.to.save=ksdf, file.name.path=file.path(path.folder, paste0(filename, "_all")))
        }
    } else {
        message("no results produced for ", filename)
    }

    return(NULL)
    # return(ksdf)
}


#' Title
#'
#' @param gene.list.to.enrich
#' @param organism.name
#' @param exclude.electr.flag
#' @param ontology
#' @param significant.flag
#' @param min.isect.size
#' @param correction.method
#' @param functional.folder
#' @param filename
#'
#' @return
#' @export
#'
#' @examples
enrichGOGProfiler <- function(gene.list.to.enrich,
                              organism.name=c("rnorvegicus", "mmusculus", "hsapiens"),
                              exclude.electr.flag=FALSE,
                              ontology=c("BP", "MF", "CC"),
                              significant.flag=FALSE,
                              min.isect.size=0,
                              correction.method="fdr",
                              functional.folder,
                              filename)
{
    # require(organism.db)
    require("gprofiler2")

    go.folder <- UpdateFolderPath(functional.folder, "gProfiler", "GO")
    filename <- UpdateFilename(filename, paste0("gProfiler_GO_", ontology))

    gprof.ontology <- paste0("GO:", ontology)

    ego <- gost(query=gene.list.to.enrich,
                organism=organism.name,
                significant=significant.flag,
                # min_isect_size=min.isect.size,
                correction_method=correction.method,
                sources=gprof.ontology,
                exclude_iea=exclude.electr.flag, evcodes=TRUE)
    ego <- ego$result
    ego$intersection <- .elaborategProfilerListToPrint(ego$intersection)
    ego$evidence_codes <- .elaborategProfilerListToPrint(ego$evidence_codes)
    ego$parents <- .elaborategProfilerListToPrint(ego$parents)
    ego<-ego[order(ego$p_value),]

    if(dim(ego)[1]!=0) {
        if(significant.flag) {
            WriteDataFrameAsTsv(data.frame.to.save=ego, file.name.path=file.path(go.folder, paste0(filename, "_significant")))
        } else {
            WriteDataFrameAsTsv(data.frame.to.save=ego[which(ego$significant==TRUE),], file.name.path=file.path(go.folder, paste0(filename, "_significant")))
            WriteDataFrameAsTsv(data.frame.to.save=ego, file.name.path=file.path(go.folder, paste0(filename, "_all")))
        }
    } else {
        message("no results produced for ", filename)
    }

    # return(ego)
    return(NULL)
}

#' Title
#'
#' @param list
#'
#' @return
#' @keywords internal
#'
#' @examples
.elaborategProfilerListToPrint <- function(list)
{
    if(max(unlist(lapply(list, length)))>1)
    {
        list <- lapply(list, function(x)
        {
            x <- paste0(x, collapse=";")
            return(x)
        })
    } else {
        list <- gsub(",", ";", list)
    }

    list <- lapply(list, function(x)
    {
        if(identical(x, "character(0)")) { x <- "NA"}
        if(nchar(x)>32767) {x <- substr(x, 1, 32765)}
        return(x)
    })
    return(unlist(list))
}

PrepareDataForKeggMap <- function(counts.dataframe,
                                  design.matrix, factors.column,
                                  kegg.results, kegg.id,
                                  to.convert.flag=FALSE,
                                  genes.map=NULL)
{
    idx <- which(kegg.results@result$ID == "rno05010")
    genes.list <- strsplit(x=kegg.results[idx]$geneID, split="/")[[1]]
    if(to.convert.flag)
    {
        if(!is.null(genes.map))
        {

        } else {


        }
    }
    genes.ens <- sign.genes.entrez$ENSEMBL[which(sign.genes.entrez$ENTREZID %in% genes.list)]
    idx.gds <- which(rownames(counts.dataframe) %in% genes.ens)
    cds <- counts.dataframe[idx.gds, , drop=FALSE]

    factors.column <- "tc"
    factors <- unique(design.matrix[[factors.column]])

    lmeans<-lapply(factors, function(x){
        fact1 <- rownames(design.matrix)[design.matrix[[factors.column]] %in% x]
        return(apply(cds[,fact1], 1, mean))
    })

    mds <- rlist::list.cbind(lmeans)
    colnames(mds) <- factors
    lfc <- cbind(mds[,1]/mds[,2], mds[,3]/mds[,4], mds[,5]/mds[,6], mds[,7]/mds[,8])
    colnames(lfc) <- unique(design.matrix$Times)

}

PlotKeggMapTimeCoursePathview <- function(interested.genes.kegg.lfc,
                                          kegg.id,
                                          specie, gene.idtype,
                                          output.dir,
                                          prefix,
                                          plot.col.key.flag=TRUE,
                                          is.signal.flag=TRUE,
                                          low.color.list=c("gene"="red", cpd="red"),
                                          mid.color.list= c("gene"="grey", cpd="grey"),
                                          high.color.list=c("gene"="green", cpd="green"))
{

    require("pathview")
    pathview(gene.data=interested.genes.kegg.lfc, pathway.id=kegg.id, species=specie,
             gene.idtype=gene.idtype, kegg.dir=output.dir, out.suffix=prefix,
             plot.col.key=plot.col.key.flag,
             low=low.color.list,
             mid=mid.color.list,
             high=high.color.list,
             is.signal=is.signal.flag,
             # keys.align="y",
             kegg.native=TRUE, multi.state=TRUE, same.layer=TRUE)
}

PlotListKeggMapTimeCoursePathview <- function(interested.genes.kegg.lfc, kegg.id.list, specie, gene.idtype, output.dir, prefix, plot.col.key.flag=TRUE, is.signal.flag=TRUE, low.color.list=c("gene"="red", cpd="red"), mid.color.list= c("gene"="grey", cpd="grey"), high.color.list=c("gene"="green", cpd="green")) {
    for(kegg.id in kegg.id.list) {
        PlotKeggMapTimeCoursePathview(interested.genes.kegg.lfc, kegg.id, specie, gene.idtype, output.dir, prefix, plot.col.key.flag, is.signal.flag, low.color.list, mid.color.list, high.color.list)
    }
}


#' Title
#'
#' @param de.gene.list
#' @param functional.folder
#' @param filename
#' @param organism
#'
#' @return
#' @export
#'
#' @examples
enrichSuitedSingleList <- function(de.gene.list, functional.folder, filename,
                                   organism="hsapiens") {

    switch(organism,
           rnorvegicus={
               org.db="org.Rn.eg.db"
               org.code="rno"
           },
           mmusculus={
               org.db="org.Mm.eg.db"
               org.code="mmu"
           },
           hsapiens={
               org.db="org.Hs.eg.db"
               org.code="hsa"
           }
    )

    dir.create(functional.folder, showWarnings=FALSE, recursive=TRUE)
    tryCatch({
        # enrichKEGGFunction(gene.list.to.enrich=de.gene.list, functional.folder=functional.folder, filename=filename, organism.db=org.db, organism.code=org.code)
        for(path in c("KEGG", "REAC"))
        {
            enrichPathwayGProfiler(gene.list.to.enrich=de.gene.list,
                                   functional.folder=functional.folder,
                                   filename=filename, path.label=path,
                                   organism.name=organism)
        }
    }, error=function(e) {
        print(e)
    })

    for (go.label in c("BP", "CC", "MF"))
    {
        tryCatch({
            # enrichGOFunction(gene.list.to.enrich=de.gene.list,
            #                  functional.folder=functional.folder,
            #                  filename=filename, ontology=go.label,
            #                  organism.db=org.db, organism.code=org.code)
            enrichGOGProfiler(gene.list.to.enrich=de.gene.list,
                              functional.folder=functional.folder,
                              filename=filename, ontology=go.label,
                              organism.name=organism)
        }, error=function(e) {
            print(e)
        })
    }
    return(NULL)

}
#
# enrichSuited <- function(de.gene.list, functional.folder, filename,
#                          organism="rnorvegicus") {
#
#     switch(organism,
#            rnorvegicus={
#                org.db="org.Rn.eg.db"
#                org.code="rno"
#            },
#            mmusculus={
#                org.db="org.Mm.eg.db"
#                org.code="mmu"
#            }
#     )
#
#     for (df.label in c("SIGN", "UP", "DOWN")) {
#         switch(df.label,
#                SIGN={
#                    filename.suit <- UpdateFilename(filename, "SIGN")
#                },
#                UP={
#                    filename.suit <- UpdateFilename(filename, "UP")
#                },
#                DOWN={
#                    filename.suit <- UpdateFilename(filename, "DOWN")
#                }
#         )
#
#         tryCatch({
#             enrichKEGGFunction(gene.list.to.enrich=rownames(de.gene.list[[df.label]]), functional.folder=functional.folder, filename=filename.suit, organism.db=org.db, organism.code=org.code)
#             for(path in c("KEGG", "REAC")) {
#                 enrichPathwayGProfiler(gene.list.to.enrich=rownames(de.gene.list[[df.label]]), functional.folder=functional.folder, filename=filename.suit, path.label=path, organism.name=organism)
#             }
#
#         }, error=function(e) {
#             print(e)
#
#         }
#         )
#
#         for (go.label in c("BP", "CC", "MF")) {
#             tryCatch({
#                 enrichGOFunction(gene.list.to.enrich=rownames(de.gene.list[[df.label]]), functional.folder=functional.folder, filename=filename.suit, ontology=go.label, organism.db=org.db, organism.code=org.code)
#                 enrichGOGProfiler(gene.list.to.enrich=rownames(de.gene.list[[df.label]]), functional.folder=functional.folder, filename=filename.suit, ontology=go.label, organism.name=organism)
#             }, error=function(e) {
#                 print(e)
#             }
#             )
#         }
#
#     }
#
# }
#
# VennDE3ListsSuited <- function(de.df.list1, de.df.list2, de.df.list3, title, venn.plot.dir, intersection.flag=TRUE, intersection.exclusion.flag=TRUE){
#
#     for (df.label in c("SIGN", "UP", "DOWN")) {
#         switch(df.label,
#                SIGN={
#                    title.p=paste(title, "Significant DE Genes")
#                    label1.p="cervical"
#                    label2.p="thoracic"
#                    label3.p="cerv-thor"
#                },
#                UP={
#                    title.p=paste(title, "UP-regulated DE Genes")
#                    label1.p=paste("cervical", "UP", sep="_")
#                    label2.p=paste("thoracic", "UP", sep="_")
#                    label3.p=paste("cerv-thor", "UP", sep="_")
#                },
#                DOWN={
#                    title.p=paste(title, "DOWN-regulated DE Genes")
#                    label1.p=paste("cervical", "DOWN", sep="_")
#                    label2.p=paste("thoracic", "DOWN", sep="_")
#                    label3.p=paste("cerv-thor", "DOWN", sep="_")
#                }
#         )
#
#         Venn3de(x=rownames(de.df.list1[[df.label]]), y=rownames(de.df.list2[[df.label]]), z=rownames(de.df.list3[[df.label]]),
#                 label1=label1.p, label2=label2.p, label3=label3.p, title=title.p,
#                 intersection.flag=intersection.flag, intersection.exclusion.flag=intersection.exclusion.flag, plot.dir=venn.plot.dir)
#
#     }
# }
