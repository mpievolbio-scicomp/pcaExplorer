#' Title
#'
#' @param se
#' @param genelist
#' @param intgroup
#' @param plotZ
#'
#' @return
#' @export
#'
#' @examples
geneprofiler <- function(se, genelist = NULL, intgroup="condition", plotZ = FALSE){
  if(is.null(genelist))
    stop("Provide at least one gene to the genelist parameter")
  # check that at least one gene is found
  genelist <- unique(genelist)
  cat("you provided",length(genelist),"unique identifiers\n")
  inthedata <- genelist %in% rownames(se)
  if (sum(inthedata) == 0)
    stop("None of the provided genes were found in the experiment data")
  cat(sum(inthedata), "out of", length(genelist), "provided genes were found in the data")

  mydata <- as.data.frame(t(assay(se)[genelist,]))

  # resort the order of the rows according to the groups that are selected
  mygroups <- interaction(as.data.frame(colData(se)[intgroup]))

  mydata <- mydata[order(mygroups),]

  if(plotZ) {
    # remove 0 variance genes
    rv <- rowVars(t(mydata))
    mydata <- mydata[,rv >0]
    ## replace maybe with explicit code!
    mydata <- NMF:::scale_mat(mydata,"col")
  }


  mylabels <- colnames(se)[order(mygroups)]
  mycols <- scales::hue_pal()(length(levels(mygroups)))[sort(mygroups)]

  par(mar=c(7.1,4.1,2.1,2.1))
  plot(mydata[,1],type="l",xaxt="n",las=2,ylim=range(mydata),xlab="",ylab=ifelse(plotZ,"scaled expression value","expression value"))
  Map(function(x,y,z)
    axis(1,at=x,col.axis=y,labels=z,lwd=0,las=2),
    1:nrow(mydata),
    mycols,
    mylabels
  )
  axis(1,at=1:nrow(mydata),labels=FALSE)

  for (i in 2:(ncol(mydata)-1)){
      lines(mydata[,i],type="l",xaxt="n",las=2,col=i)
  }
  ## TODO: if desired, plot only the avg pro group -> maybe as boxplot?

}
