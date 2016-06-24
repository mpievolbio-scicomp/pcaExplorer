#' Title
#'
#' @param rld
#' @param plot_type
#'
#' @return
#' @export
#'
#' @examples
distro_rlog <- function(rld, plot_type="density") {
  allrld <- tidyr::gather(as.data.frame(assay(rld)))
  names(allrld) <- c("Sample","rlogExpression")

  if(plot_type=="boxplot"){
    p <- ggplot(allrld,aes(x=Sample,y=rlogExpression)) + geom_boxplot(aes(col=Sample,fill=Sample),alpha=0.5)
  }

  if(plot_type=="violin"){
    p <- ggplot(allrld,aes(x=Sample,y=rlogExpression)) + geom_violin(aes(col=Sample,fill=Sample),alpha=0.5)
  }

  if(plot_type=="density"){
    p <- ggplot(allrld,aes(x=rlogExpression)) + geom_density(aes(color=Sample),alpha=0.1)
  }
  return(p)
}
