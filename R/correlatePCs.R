# p <- prcomp(t(assay(rld_airway)))
# cd <- as.data.frame(colData(rld_airway))



#' Title
#'
#' @param pcaobj
#' @param coldata
#' @param pcs
#'
#' @return
#' @export
#'
#' @examples
correlatePCs <- function(pcaobj,coldata,pcs=1:4){
  # split the analysis for continuous and categorial

  coldataTypes <- sapply(coldata,class)

  # extract the scores from the pc object
  x <- pcaobj$x

  # do it until 1:4 PCs

  res <- matrix(NA,nrow=length(pcs),ncol = ncol(coldata))

  colnames(res) <- colnames(cd)
  rownames(res) <- paste0("PC_",pcs)

  for (i in 1:ncol(res)){
    # for each covariate...
    for(j in pcs){
      if(coldataTypes[i] %in% c("factor","character"))
      {
        if(length(levels(coldata[,i])) > 1) {
          res[j,i] <- kruskal.test(x[,j],coldata[,i])$p.value
        }
      }else {
        res[j,i] <- cor.test(x[,j],coldata[,i],method="spearman")$p.value

      }
    }
  }
  res
}


#' Title
#'
#' @param pccorrs
#' @param pc
#' @param logp
#'
#' @return
#' @export
#'
#' @examples
plotPCcorrs <- function(pccorrs,pc=1,logp=TRUE) {
  selectedPC <- paste0("PC_",pc)
  pvals <- pccorrs[selectedPC,]

  if(logp) pvals <- -log10(pvals)

  barplot(pvals,las=2, col="steelblue",
          main=paste0("Correlations PC ",pc," vs covariates"),
          ylab=ifelse(logp,"-log10(pval)","pval"))
}



