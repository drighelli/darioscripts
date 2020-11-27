#' Title
#'
#' @param path
#' @param prefix
#' @param plot.type
#'
#' @return
#' @export
#'
#' @examples
GeneratePlotStrings <- function(path=NULL, prefix, plot.type) {
  title <- gsub(pattern = "_", replacement = " ", x = UpdatePrefix(prefix, plot.type))

  plot.folder <- gsub(pattern = " ", replacement = "_", x = file.path(path, plot.type))

  plot.file.name <- gsub(pattern = " ", replacement = "_", x = UpdatePrefix(prefix, plot.type))
  if(!is.null(path)) dir.create(plot.folder, showWarnings = FALSE, recursive = TRUE)

  return(list("title"= title, "plot.folder"=plot.folder, "plot.file.name"=plot.file.name))
}


#' Title
#'
#' @param prefix
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
UpdatePrefix <- function(prefix, ...) {
  # new.prefix <- paste(prefix, postix, sep=sep)
  dots <- list(...)
  if(length(dots) != 0) {
    for (str in dots) {
      # str <- gsub(pattern = ".", replacement = "_", str)
      prefix <- paste(prefix, str, sep = " " )
    }

  }else {
    stop("provide a string to append to ", new.prefix)
  }
  return(prefix)
}

#' Title
#'
#' @param path
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
UpdateFolderPath <- function(path, ...) {
  dots <- list(...)
  if(length(dots) != 0) {
    for (str in dots) {
        str <- gsub(pattern = " ", replacement = "_", str)
        path <- file.path(path, str)
      }
  } else {
    stop("provide a string to append to ", path)
  }
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  return(path)
}

#' Title
#'
#' @param filename
#' @param ...
#' @param extension
#'
#' @return
#' @export
#'
#' @examples
UpdateFilename <- function(filename, ..., extension=NULL) {
  dots <- list(...)
  filename <- gsub(pattern = " ", replacement = "_", x = filename)
  if(length(dots) != 0) {
    for (str in dots) filename <- paste(filename, str, sep = "_")
  } else {
    stop("provide a string to append to ", filename)
  }
  if(!is.null(extension)) {
    filename <- paste0(filename, ".", extension)
  }
  return(filename)
}
#
# ProduceOutputStrings <- function(prefix, out.dir, ..., filename=NULL) {
#   if(class(...)!="character") {
#     dots <- list(...)
#   } else {
#     dots <-
#   }
#
#
#   for(lbl in labels.list) {
#     prefix <- UpdatePrefix(prefix, lbl)
#   }
#
#   out.dir <- UpdateFolderPath(root.dir, prefix)
#   filename <- UpdateFilename(prefix, "genes")
#
#
#
# }

# ## global variables
# project.name <<- "cervical_thoracic_whole"
# results.folder <<- UpdateFolderPath(project.name, "results")
# counts.folder <<- UpdateFolderPath(results.folder, "counts")
# plots.folder <<- UpdateFolderPath(results.folder, "plots")
# de.plots.folder <<- UpdateFolderPath(plots.folder, "de_plots")
# boxplot.folder <<- UpdateFolderPath(plots.folder, "boxplots")
# pca.folder <<- UpdateFolderPath(plots.folder, "PCA")
# file.name <<- project.name



