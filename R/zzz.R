#' @importFrom shiny addResourcePath

.onLoad <- function(libname, pkgname) {
  # Create link to logo
  # shiny::addResourcePath("pcaExplorer", system.file("www", package="pcaExplorer"))
  
  shiny::addResourcePath("sbs", system.file("www", package = "shinyBS"))
}

