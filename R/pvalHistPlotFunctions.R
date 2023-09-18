

#' generateGGHist
#'
#' @param processed.de.results
#' @param strings
#' @param nbins
#'
#' @return
#' @export
#' @import ggplot2
#' @examples
generateGGHist <- function (processed.de.results, strings, nbins=30,
    pval=c("pval", "padj"))
{
    pval <- match.arg(pval)
    hh <- ggplot(data=processed.de.results, aes_string(pval)) +
        geom_histogram(col="red", fill="white", bins=nbins) +
        labs(title=strings$title) +
        labs(x="PValues", y="Count")
    return(hh)
}




#' PlotHistPvalPlot
#'
#' @param de.results
#' @param design.matrix
#' @param show.plot.flag
#' @param plotly.flag
#' @param save.plot
#' @param plot.folder
#' @param prefix.plot
#' @param threshold
#' @param nbins
#' @param pvals
#' @param cowplot
#'
#' @return
#' @export
#' @importFrom plotly ggplotly
#' @importFrom cowplot theme_cowplot
#' @examples
PlotHistPvalPlot <- function(de.results, design.matrix,
                            show.plot.flag=TRUE, plotly.flag=FALSE,
                            save.plot=FALSE, plot.folder=NULL,
                            prefix.plot=NULL,
                            threshold=0.05, nbins=30, pvals=c("pval", "padj"),
                            cowplot=FALSE)
{
    pvals <- match.arg(pvals)
    strings <- GeneratePlotStrings(path=plot.folder,
                                    prefix=prefix.plot,
                                    plot.type="Pval Histogram")

    ## this function has to be generalized
    processed.de.results <- ProcessDEResultsForPlot(de.results,
                                                threshold=threshold,
                                                design.matrix=design.matrix)

    ggp <- generateGGHist(processed.de.results=processed.de.results,
                          strings=strings, nbins=nbins, pval=pvals)

    if(cowplot) ggp <- ggp+theme_cowplot()
    if(save.plot)
    {
        if(is.null(plot.folder)) {
            stop("Please set a folder where to plot the MA-Plot!")
        }
        if(!is.null(strings$plot.file.name)){
            # SaveGGplot(ggplot.to.save=ggp, plot.folder=strings$plot.folder,
            #    plot.file.name=strings$plot.file.name, plotly.flag=plotly.flag)
        }

    }

    if(show.plot.flag)
    {
        if(plotly.flag) {
            plotly::ggplotly(ggp)
        } else {
            plot(ggp)
        }
    } else {
        return(ggp)
    }
}
