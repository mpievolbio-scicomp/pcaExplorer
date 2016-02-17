#' Title
#'
#' @param x
#' @param intgroup
#' @param ntop
#' @param returnData
#' @param title
#' @param pcX
#' @param pcY
#' @param text_labels
#' @param point_size
#'
#' @return A value
#' @export
#'
#' @examples # An example
pcaplot <- function (x, intgroup = "condition", ntop = 500, returnData = FALSE,title=NULL,
                    pcX = 1, pcY = 2,text_labels=TRUE,point_size=3) # customized principal components
{
  library("DESeq2")
  library("genefilter")
  library("ggplot2")
  rv <- rowVars(assay(x))
  select <- order(rv, decreasing = TRUE)[seq_len(min(ntop,length(rv)))]
  pca <- prcomp(t(assay(x)[select, ]))

  percentVar <- pca$sdev^2/sum(pca$sdev^2)

  if (!all(intgroup %in% names(colData(x)))) {
    stop("the argument 'intgroup' should specify columns of colData(dds)")
  }
  intgroup.df <- as.data.frame(colData(x)[, intgroup, drop = FALSE])
  group <- factor(apply(intgroup.df, 1, paste, collapse = " : "))
  d <- data.frame(PC1 = pca$x[, pcX], PC2 = pca$x[, pcY], group = group,
                  intgroup.df, names = colnames(x))
  colnames(d)[1] <- paste0("PC",pcX)
  colnames(d)[2] <- paste0("PC",pcY)

  if (returnData) {
    attr(d, "percentVar") <- percentVar[1:2]
    return(d)
  }

  # clever way of positioning the labels
  d$hjust = ifelse((sign(d[,paste0("PC",pcX)])==1),0.9,0.1)# (1 + varname.adjust * sign(PC1))/2)

  g <- ggplot(data = d, aes_string(x = paste0("PC",pcX), y = paste0("PC",pcY), color = "group")) +
    geom_point(size = point_size) +
    xlab(paste0("PC",pcX,": ", round(percentVar[pcX] * 100,digits = 2), "% variance")) +
    ylab(paste0("PC",pcY,": ", round(percentVar[pcY] * 100,digits = 2), "% variance"))

  library("ggrepel")
  if(text_labels) g <- g + geom_label_repel(mapping = aes(label=names,fill=group),color="white", show.legend = T) +theme_bw()
  if(!is.null(title)) g <- g + ggtitle(title)
  g
}

library(topGO)

pcascree <- function(obj, type = c("pev", "cev"),pc_nr=NULL,title=NULL)
{
  type <- match.arg(type)
  d <- obj$sdev^2
  yvar <- switch(type, pev = d/sum(d), cev = cumsum(d)/sum(d))
  yvar.lab <- switch(type, pev = "proportion of explained variance",
                     cev = "cumulative proportion of explained variance")
  # df <- data.frame(PC = 1:length(d), yvar = yvar)

  if (!is.null(pc_nr)) {
    colsize <- pc_nr
    yvar <- yvar[1:pc_nr]
  } else {
    colsize <- length(d)
    yvar <- yvar[1:length(d)]
  }

  pc_df <- data.frame(PC_count = 1:colsize, var = yvar)

  if(type=="pev"){
    p <- ggplot(pc_df, aes(x = PC_count, y = var)) + geom_bar(stat = "identity")
    p <- p + scale_x_continuous(breaks = 1:length(d))
    p <- p + ylab(yvar.lab) + xlab("principal components")
    # p
  } else {
    p <- ggplot(pc_df, aes(x = PC_count, y = var)) + geom_point() + geom_path() + scale_x_continuous(breaks = 1:length(d))
    p <- p + ylab(yvar.lab) + xlab("principal components") + ylim(0,max(pc_df$var))
    # p
  }
  if(!is.null(title)) p <- p + ggtitle(title)
  p
}
