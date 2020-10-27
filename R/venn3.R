Venn3de <- function(x, y, z, label1, label2, label3, title="Venn Diagram", intersection.flag=TRUE, intersection.exclusion.flag=FALSE, plot.dir=NULL, enrich.lists.flag=FALSE, conversion.map=NULL){
    
    a=x
    b=y
    c=z
    plot.combined.string <- paste(label1,"_",label2,"_",label3, sep="")
    # a = row.names(x)
    # b = row.names(y)
    # c = row.names(z)
    if(is.null(plot.dir)) {
      stop("Please provide a directory where to save VENN results!")
    }
    
    Lists <- list(a, b, c)  #put the word vectors into a list to supply lapply  
    Lists <- lapply(Lists, function(x) as.character(unlist(x)))
    items <- sort(unique(unlist(Lists)))   #put in alphabetical order
    MAT <- matrix(rep(0, length(items)*length(Lists)), ncol=3)  #make a matrix of 0's
    names <- c(label1,label2,label3)                
    colnames(MAT) <- names
    rownames(MAT) <- items
    lapply(seq_along(Lists), function(i) {   #fill the matrix
      MAT[items %in% Lists[[i]], i] <<- table(Lists[[i]])
    })
    
    outputName <- UpdateFilename(filename = "VennDiagram", label1, label2, label3, extension="pdf")
    out.path.name <- UpdateFolderPath(plot.dir, "Venn3")
    prefix="venn3"
    # outputName=paste(label1,"_",label2,"_",label3,"_VennDiagram.pdf",sep="")
    
    # file.path.name <- file.path(plot.dir, "Venn3")
    # dir.create(file.path.name, recursive = TRUE, showWarnings = FALSE)
    file.path.name <- file.path(out.path.name, outputName)
    dev.new()
    require("limma")
    limma::vennDiagram(MAT, circle.col= c("red","green","yellow"), main=title)
    
    dev.print(device = pdf, file=file.path.name, width=10, height=10)
    
    if(intersection.flag) {
      if(is.null(plot.dir)) {
        stop("Please provide a directory where to save VENN results!")
      }
      
      # res.path <- file.path(plot.dir, "Venn3", paste0(plot.combined.string,"_gene_lists"))
      # dir.create(res.path, recursive = TRUE, showWarnings = FALSE)
      ab <- intersect(a, b)
      
      # if(length(ab)==0) print("merd!!!!")
      SaveInteserctionsList(gene.list = ab, conversion.map = conversion.map, root.dir = out.path.name, prefix = prefix, labels.list = c(label1, paste0("AND_",label2)), enrich.lists.flag=enrich.lists.flag)
      # if(enrich.lists.flag) {
      #   enrichSuitedSingleList(de.gene.list = ab, functional.folder = file.path(res.path, "functional"), filename = filename)
      # }
      
      # if(!is.null(conversion.map)) {
      #   ab <- CreateConvertedDataframe(ab, ens.symb.biomart.map)
      # }
      # 
      # filename <- paste(label1,"_",label2,"_genes_in_intersection.txt",sep="")
      # filepathname <- file.path(res.path, filename)
      # write.table(ab, file = filepathname , quote=FALSE, sep="\t", row.names=FALSE)
      # 
      bc <- intersect(b, c)
      SaveInteserctionsList(gene.list = bc, conversion.map = conversion.map, root.dir = out.path.name, prefix = prefix, labels.list = c(label2, paste0("AND_", label3)), enrich.lists.flag=enrich.lists.flag)
      # bc<-AttachConvertedColumn(bc, ens.symb.biomart.map)
      # filename <- paste(label2,"_",label3,"_genes_in_intersection.txt",sep="")
      # filepathname <- file.path(res.path, filename)
      # write.table(bc, file = filepathname , quote=FALSE, sep="\t", row.names=FALSE)
      # if(enrich.lists.flag) {
      #   enrichSuitedSingleList(de.gene.list = bc, functional.folder = file.path(res.path, "functional"), filename = filename)
      # }
      
      ac <- intersect(a, c)
      SaveInteserctionsList(gene.list = ac, conversion.map = conversion.map, root.dir = out.path.name, prefix = prefix, labels.list = c(label1, paste0("AND_", label3)), enrich.lists.flag=enrich.lists.flag)
      # filename <- paste(label1,"_",label3,"_genes_in_intersection.txt",sep="")
      # filepathname <- file.path(res.path, filename)
      # # ac<-AttachConvertedColumn(ac, ens.symb.biomart.map)
      # write.table(ac, file = filepathname , quote=FALSE, sep="\t", row.names=FALSE)
      # if(enrich.lists.flag) {
      #   enrichSuitedSingleList(de.gene.list = ac, functional.folder = file.path(res.path, "functional"), filename = filename)
      # }
      
      abc <- intersect(ab, bc)
      SaveInteserctionsList(gene.list = abc, conversion.map = conversion.map, root.dir = out.path.name, prefix = prefix, labels.list = c(label1, paste0("AND_", label2), paste0("AND_", label3)), enrich.lists.flag=enrich.lists.flag)
      # abc<-AttachConvertedColumn(abc, ens.symb.biomart.map)
      # filename <- paste(label1, "_", label2, "_", label3, "_genes_in_intersection.txt",sep="")
      # filepathname <- file.path(res.path, filename)
      # write.table(abc, file = filepathname , quote=FALSE, sep="\t", row.names=FALSE)
      # if(enrich.lists.flag) {
      #   enrichSuitedSingleList(de.gene.list = abc, functional.folder = file.path(res.path, "functional"), filename = filename)
      # }
      
      union.of.intersections <- union(union(ab, bc), union(ac, abc))
      SaveInteserctionsList(gene.list = union.of.intersections, conversion.map = conversion.map, root.dir = out.path.name, prefix = paste0(prefix, "_union_of_intersections"), labels.list = c(label1, label2, label3), enrich.lists.flag=enrich.lists.flag)
      
      union.of.lugs <- union(bc, union(ac, abc))
      SaveInteserctionsList(gene.list = union.of.lugs, conversion.map = conversion.map, root.dir = out.path.name, prefix = paste0(prefix, "_union_of_LUGs"), labels.list = c(label1, label2, label3), enrich.lists.flag=enrich.lists.flag)
      # filename <- paste(label1, "_", label2, "_", label3, "_genes_in_union_of_intersections.txt",sep="")
      # filepathname <- file.path(res.path, filename)
      # write.table(union.of.intersections, file = filepathname , quote=FALSE, sep="\t", row.names=FALSE)
      # if(enrich.union.intersections) {
      #   enrichSuitedSingleList(de.gene.list = union.of.intersections, functional.folder = file.path(res.path, "functional"), filename = filename)
      # }
      
    }
    
    if(intersection.exclusion.flag) {
      if(is.null(plot.dir)) {
        stop("Please provide a directory where to save VENN results!")
      }
      res.path <- file.path(plot.dir, "Venn3", paste0(plot.combined.string,"_gene_lists"))
      dir.create(res.path, recursive = TRUE, showWarnings = FALSE)
      
      a.not.b <-  setdiff(a, b)
      SaveInteserctionsList(gene.list = a.not.b, conversion.map = conversion.map, root.dir = out.path.name, prefix = prefix, 
                            labels.list = c(label1, paste0("_NOT_",label2), enrich.lists.flag=enrich.lists.flag))
      # a.not.b<-AttachConvertedColumn(a.not.b, ens.symb.biomart.map)
      # filename <- paste(label1, "_NOT_", label2, "_genes.txt",sep="")
      # filepathname <- file.path(res.path, filename)
      # write.table(a.not.b, file = filepathname , quote=FALSE, sep="\t", row.names=FALSE)
      # if(enrich.lists.flag) {
      #   enrichSuitedSingleList(de.gene.list = a.not.b, functional.folder = file.path(res.path, "functional"), filename = filename)
      # }
      
      b.not.a <-  setdiff(b, a)
      SaveInteserctionsList(gene.list = b.not.a, conversion.map = conversion.map, root.dir = out.path.name, prefix = prefix, 
                            labels.list = c(label2, paste0("_NOT_",label1)), enrich.lists.flag=enrich.lists.flag)
      # 
      # filename <- paste(label2, "_NOT_", label1, "_genes.txt",sep="")
      # filepathname <- file.path(res.path, filename)
      # write.table(b.not.a, file = filepathname , quote=FALSE, sep="\t", row.names=FALSE)
      # if(enrich.lists.flag) {
      #   enrichSuitedSingleList(de.gene.list = b.not.a, functional.folder = file.path(res.path, "functional"), filename = filename)
      # }
      
      c.not.b <-  setdiff(c, b)
      SaveInteserctionsList(gene.list = c.not.b, conversion.map = conversion.map, root.dir = out.path.name, prefix = prefix, 
                            labels.list = c(label3, paste0("_NOT_",label2)), enrich.lists.flag=enrich.lists.flag)
      # filename <- paste(label3, "_NOT_", label2, "_genes.txt",sep="")
      # filepathname <- file.path(res.path, filename)
      # write.table(c.not.b, file = filepathname , quote=FALSE, sep="\t", row.names=FALSE)
      # if(enrich.lists.flag) {
      #   enrichSuitedSingleList(de.gene.list = c.not.b, functional.folder = file.path(res.path, "functional"), filename = filename)
      # }
      
      b.not.c <-  setdiff(b, c)
      SaveInteserctionsList(gene.list = b.not.c, conversion.map = conversion.map, root.dir = out.path.name, prefix = prefix, 
                            labels.list = c(label2, paste0("_NOT_", label3)), enrich.lists.flag=enrich.lists.flag)
      
      # filename <- paste(label2, "_NOT_", label3, "_genes.txt",sep="")
      # filepathname <- file.path(res.path, filename)
      # write.table(b.not.c, file = filepathname , quote=FALSE, sep="\t", row.names=FALSE)
      # if(enrich.lists.flag) {
      #   enrichSuitedSingleList(de.gene.list = b.not.c, functional.folder = file.path(res.path, "functional"), filename = filename)
      # }
      
    
      a.not.c <-  setdiff(a, c)
      SaveInteserctionsList(gene.list = a.not.c, conversion.map = conversion.map, root.dir = out.path.name, prefix = prefix, 
                            labels.list = c(label1, paste0("_NOT_", label3)), enrich.lists.flag=enrich.lists.flag)
      
      # filename <- paste(label1, "_NOT_", label3, "_genes.txt",sep="")
      # filepathname <- file.path(res.path, filename)
      # write.table(a.not.c, file = filepathname , quote=FALSE, sep="\t", row.names=FALSE)
      # if(enrich.lists.flag) {
      #   enrichSuitedSingleList(de.gene.list = a.not.c, functional.folder = file.path(res.path, "functional"), filename = filename)
      # }
      
      c.not.a <-  setdiff(c, a)
      SaveInteserctionsList(gene.list = c.not.a, conversion.map = conversion.map, root.dir = out.path.name, prefix = prefix, 
                            labels.list = c(label3, paste0("_NOT_", label1)), enrich.lists.flag=enrich.lists.flag)
      
      # filename <- paste(label3, "_NOT_", label1, "_genes.txt",sep="")
      # filepathname <- file.path(res.path, filename)
      # write.table(c.not.a, file = filepathname , quote=FALSE, sep="\t", row.names=FALSE)
      # if(enrich.lists.flag) {
      #   enrichSuitedSingleList(de.gene.list = a.not.c, functional.folder = file.path(res.path, "functional"), filename = filename)
      # }
      
      a.not.bc <- setdiff(a.not.b, c)
      SaveInteserctionsList(gene.list = a.not.bc, conversion.map = conversion.map, root.dir = out.path.name, prefix = prefix, 
                            labels.list = c(label1, paste0("_NOT_", label2), paste0("_NOT_", label3)), enrich.lists.flag=enrich.lists.flag)
      # 
      # filename <- paste(label1, "_NOT_", label2, "_NOT_", label3, "_genes.txt",sep="")
      # filepathname <- file.path(res.path, filename)
      # write.table(a.not.bc, file = filepathname , quote=FALSE, sep="\t", row.names=FALSE)
      # if(enrich.lists.flag) {
      #   enrichSuitedSingleList(de.gene.list = a.not.bc, functional.folder = file.path(res.path, "functional"), filename = filename)
      # }
      
      b.not.ac <- setdiff(b.not.a, c)
      SaveInteserctionsList(gene.list = b.not.ac, conversion.map = conversion.map, root.dir = out.path.name, prefix = prefix, 
                            labels.list = c(label2, paste0("_NOT_", label1), paste0("_NOT_", label3)), enrich.lists.flag=enrich.lists.flag)
      
      # filename <- paste(label2, "_NOT_", label1, "_NOT_", label3, "_genes.txt",sep="")
      # filepathname <- file.path(res.path, filename)
      # write.table(b.not.ac, file = filepathname , quote=FALSE, sep="\t", row.names=FALSE)
      # if(enrich.lists.flag) {
      #   enrichSuitedSingleList(de.gene.list = b.not.ac, functional.folder = file.path(res.path, "functional"), filename = filename)
      # }
      
      c.not.ab <- setdiff(c.not.a, b)
      SaveInteserctionsList(gene.list = c.not.ab, conversion.map = conversion.map, root.dir = out.path.name, prefix = prefix, 
                            labels.list = c(label3, paste0("_NOT_", label1), paste0("_NOT_", label2)), enrich.lists.flag=enrich.lists.flag)
      
      # 
      # filename <- paste(label3, "_NOT_", label1, "_NOT_", label2, "_genes.txt",sep="")
      # filepathname <- file.path(res.path, filename)
      # write.table(c.not.ab, file = filepathname , quote=FALSE, sep="\t", row.names=FALSE)
      # if(enrich.lists.flag) {
      #   enrichSuitedSingleList(de.gene.list = c.not.ab, functional.folder = file.path(res.path, "functional"), filename = filename)
      # }
    }
    
    # Lists <- list(a, b, c)  #put the word vectors into a list to supply lapply  
    # Lists <- lapply(Lists, function(x) as.character(unlist(x)))
    # items <- sort(unique(unlist(Lists)))   #put in alphabetical order
    # MAT <- matrix(rep(0, length(items)*length(Lists)), ncol=3)  #make a matrix of 0's
    # names <- c(label1,label2,label3)                
    # colnames(MAT) <- names
    # rownames(MAT) <- items
    # lapply(seq_along(Lists), function(i) {   #fill the matrix
    #   MAT[items %in% Lists[[i]], i] <<- table(Lists[[i]])
    # })
    # 
    # outputName=paste(label1,"_",label2,"_",label3,"_VennDiagram.pdf",sep="")
    # 
    # #b=paste(a,outputName,sep="")
    # 
    # file.path.name <- file.path(plot.dir, "Venn3")
    # dir.create(file.path.name, recursive = TRUE, showWarnings = FALSE)
    # file.path.name <- paste0(file.path.name, "/", outputName)
    # dev.new()
    # require("limma")
    # limma::vennDiagram(MAT, circle.col= c("red","green","yellow"), main=title)
    # 
    # dev.print(device = pdf, file=file.path.name, width=10, height=10)
    graphics.off()
  }


SaveInteserctionsList <- function(gene.list, conversion.map, root.dir, prefix, labels.list, enrich.lists.flag=FALSE) {
  
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
  
  WriteDataFrameAsTsv(data.frame.to.save = gene.list.df, file.name.path = file.path(out.dir, filename), col.names= TRUE, row.names = FALSE)

  if(enrich.lists.flag) {
    enrichSuitedSingleList(de.gene.list = gene.list, functional.folder = file.path(out.dir, "functional"), filename = filename)
  }

}
