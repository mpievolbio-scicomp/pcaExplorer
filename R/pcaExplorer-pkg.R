#' pcaExplorer: analyzing time-lapse microscopy imaging, from detection to tracking
#'
#' pcaExplorer provides functionality for interactive visualization of RNA-seq datasets
#' based on Principal Components Analysis. The methods provided allow for quick information
#' extraction and effective data exploration. A Shiny application encapsulates the whole analysis.
#'
#' pcaExplorer provides functionality for interactive visualization of RNA-seq datasets
#' based on Principal Components Analysis. The methods provided allow for quick information
#' extraction and effective data exploration. A Shiny application encapsulates the whole analysis.
#'
#' @import DESeq2
#' @import SummarizedExperiment
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
#' @importFrom S4Vectors DataFrame
#' @importFrom genefilter rowVars
#' @import d3heatmap
#' @importFrom scales brewer_pal hue_pal
#' @importFrom NMF aheatmap
#' @import plyr
#' @importFrom limma goana topGO
#' @importFrom AnnotationDbi select Term mapIds
#' @importMethodsFrom GOstats hyperGTest summary
#' @import GO.db
#' @import shiny
#' @import shinydashboard
#' @importFrom shinyBS bsTooltip
#' @import ggplot2
#' @import ggrepel
#' @importFrom DT datatable
#' @import methods
#'
#' @author
#' Federico Marini \email{marinif@@uni-mainz.de}, 2016
#'
#' Maintainer: Federico Marini \email{marinif@@uni-mainz.de}
#' @name pcaExplorer
#' @docType package
NULL
