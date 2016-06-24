#' Plot distribution of expression values
#'
#' @param rld A \code{\link{DESeqTransform}} object.
#' @param plot_type Character, choose one of \code{boxplot}, \code{violin} or
#' \code{density}. Defaults to \code{density}
#'
#' @return
#' @export
#'
#' @examples
distro_expr <- function(rld, plot_type="density") {
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
