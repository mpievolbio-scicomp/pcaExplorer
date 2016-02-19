
# x <- rld_deplall
# rv <- rowVars(assay(x))
# select <- order(rv, decreasing = TRUE)[seq_len(min(nrow(x),length(rv)))]
# pcaobj <- prcomp(t(assay(x)[select, ]),scale = F, center = T)



#' Title
#'
#' @param pcaobj
#' @param whichpc
#' @param topN
#' @param exprTable
#' @param annotation
#' @param title
#'
#' @return
#' @export
#'
#' @examples
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



