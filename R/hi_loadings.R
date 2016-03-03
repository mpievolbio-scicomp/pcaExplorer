
# x <- rld_deplall
# rv <- rowVars(assay(x))
# select <- order(rv, decreasing = TRUE)[seq_len(min(nrow(x),length(rv)))]
# pcaobj <- prcomp(t(assay(x)[select, ]),scale = F, center = T)



#' Extract genes with highest loadings
#'
#' @param pcaobj A \code{prcomp} object
#' @param whichpc An integer number, corresponding to the principal component of
#' interest
#' @param topN Integer, number of genes with top and bottom loadings
#' @param exprTable A \code{matrix} object, e.g. the counts of a \code{\link{DESeqDataSet}}.
#' If not NULL, returns the counts matrix for the selected genes
#' @param annotation A \code{data.frame} object, with row.names as gene identifiers (e.g. ENSEMBL ids)
#' and a column, \code{gene_name}, containing e.g. HGNC-based gene symbols
#' @param title The title of the plot
#'
#' @return A base plot object, or a \code{matrix}, if \code{exprTable} is not null
#'
#'
#' @examples
#'
#'
#' @export
hi_loadings <- function(pcaobj, whichpc = 1, topN = 10, exprTable = NULL,
                        annotation = NULL, title="Top/bottom loadings - "){
  if(whichpc < 0)
    stop("Use a positive integer value for the principal component to select")
  if(whichpc > nrow(pcaobj$x))
    stop("You can not explore a principal component that is not in the data")

  geneloadings_sorted <- sort(pcaobj$rotation[,whichpc])
  geneloadings_extreme <- c(tail(geneloadings_sorted,topN),head(geneloadings_sorted,topN))

  if(!is.null(exprTable)) {
    tab <- exprTable[names(geneloadings_extreme),]
    if(!is.null(annotation))
      rownames(tab) <- annotation$gene_name[match(rownames(tab),rownames(annotation))]
    return(tab)
  }

  if(!is.null(annotation))
    names(geneloadings_extreme) <- annotation$gene_name[match(names(geneloadings_extreme),rownames(annotation))]

  barplot(geneloadings_extreme, las=2, col = c(rep("steelblue",topN),rep("coral",topN)),
          main=paste0(title, "PC", whichpc))
  # p <- ggplot(data.frame(loadings=geneloadings_extreme,geneID=names(geneloadings_extreme)),aes(geneID,weight=loadings)) + geom_bar()
}



