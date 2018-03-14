

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





