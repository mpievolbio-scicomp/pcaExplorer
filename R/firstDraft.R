## TODOS
# parameters instead of fixed values - MOSTLY DONE
# some other tweaks - DOING THEM; DONE MANY
# try one run from the scratch - DONE
# list whatever is required on top of the pca_SUPALIVE function - GUESS WE'RE DONE
### EMBED THIS INTO INTERACTIVE HTML REPORTS? e.g. a la maplots!!!!

library("DESeq2")
library("genefilter")
library("shiny")

# load("/Volumes/users$/marinif/linuxHome/032-ruf-macrophages/cm2.RData")
# load("goEnrichs.RData")
# load("/Volumes/users$/marinif/linuxHome/F07-schubert/seb_dds_rld.RData")
# load("/Volumes/users$/marinif/linuxHome/F07-schubert/goEnrichs_rld_deplall.RData")
# load("/Volumes/users$/marinif/Development/presentations/2015_12_07-08_bioconductor_developers_meeting/EXAMPLE_ruf.RData")
# load("/Volumes/users$/marinif/Development/presentations/2015_12_07-08_bioconductor_developers_meeting/EXAMPLE_kleinert.RData")
# load("/Volumes/users$/marinif/Development/presentations/2015_12_07-08_bioconductor_developers_meeting/goEnrichs_paired.RData")


# library(ggvis)
# library(shinythemes)


library("DT")
library("pheatmap")
library("d3heatmap")
library(shinydashboard)

# annotation <- data.frame(gene_id=rownames(cm2),gene_name=cm2$fromgtf,stringsAsFactors = FALSE,row.names = rownames(cm2))


footer<-function(){
  tags$div(
    class = "footer",
    style = "text-align:center",
    tags$div(class = "foot-inner",
             list(
               hr(),
               "pcaExplorer is a project developed by Federico Marini in the Bioinformatics division of the ",
               tags$a(href="http://www.unimedizin-mainz.de/imbei","IMBEI"),
               ". ",br(),
               "Development of the pcaExplorer package is on ",
               tags$a(href="https://github.com/federicomarini/pcaExplorer", "GitHub")
             )

    )
  )
}

#' Title
#'
#' @param obj
#' @param obj2
#' @param countmatrix
#' @param coldata
#' @param pca2go
#' @param annotation
#'
#' @import DESeq2
#' @import SummarizedExperiment
#' @import genefilter
#' @import pheatmap
#' @import d3heatmap
#' @import scales
#' @import NMF
#' @import plyr
#' @import grid
#' @import topGO
#' @import GO.db
#' @import shiny
#' @import shinydashboard
#' @import shinyBS
#' @import ggplot2
#' @import ggrepel
#' @import DT
#' @import methods
#'
#' @return
#' @export
#'
#' @examples
pcaExplorer <- function(obj=NULL,
                        obj2=NULL,
                        countmatrix=NULL,
                        coldata=NULL,
                        pca2go=NULL,
                        annotation=NULL){
  # stopifnot( is(obj, 'SummarizedExperiment') )
  ## plenty of tests on the selected object(s)
  if ( !require('shiny') ) {
    stop("pca_SUPALIVE requires 'shiny'. Please install it using
         install.packages('shiny')")
  }
  library("shinyURL")
  # here something like, if provided as a countmatrix, create dds and rld

#   if(!is.null(obj)){
#     poss_covars <- names(colData(obj))[] # exclude the size factor, not really required?
#   } else {
#     poss_covars <- c()
#   }
  ## TODO: colselection also passed as a palette?
  colSelection <- c("navyblue","steelblue","skyblue","darkred","coral3","darksalmon","green4","greenyellow","orange","gold")
  # alternative to evaluate: use other palettes, eg hue_pal
  ## place it in the server to be reactive on the ncol of the object!




#   cmcm <- counts(dds_deplall)[1:10,]
#   write.table(cmcm,file="minicm.txt",quote=F,sep="\t",row.names=TRUE,col.names=TRUE)
#   ddd <- colData(dds_deplall)
#   write.table(ddd,file="minicoldata.txt",quote=F,sep="\t",row.names=TRUE,col.names=TRUE)
  # write.table(counts(dds_deplall),file="fullcm.txt",quote=F,sep="\t",row.names=TRUE,col.names=TRUE)
  # write.table(annotation,file="anno_mouseEns.txt",quote=F,sep="\t",row.names=TRUE,col.names=TRUE)
  ## ------------------------------------------------------------------ ##
  ##                          Define UI                                 ##
  ## ------------------------------------------------------------------ ##

  # poss_covars <- names(colData(obj))
  # poss_covars <- NULL
  newuiui <-
    shinydashboard::dashboardPage(
      dashboardHeader(
        title = paste0("pcaExplorer - Interactive exploration of Principal Components",
                       "of Samples and Genes in RNA-seq data - version ",
                       packageVersion("pcaExplorer")),
        titleWidth = 900),

      dashboardSidebar(
        width = 250,
        menuItem("Data upload",icon = icon("upload"),
                 uiOutput("upload_count_matrix"),
                 shinyBS::bsTooltip(
                   "upload_count_matrix", paste0("Select file containing the count matrix"),
                   "right", options = list(container = "body")),
                 uiOutput("upload_metadata"),
                 shinyBS::bsTooltip(
                   "upload_metadata", paste0("Select file containing the samples metadata"),
                   "right", options = list(container = "body")),
                 uiOutput("upload_annotation"),
                 shinyBS::bsTooltip(
                   "upload_annotation", paste0("Select file containing the annotation data"),
                   "right", options = list(container = "body"))),
        menuItem("App settings",icon = icon("cogs"),
                 selectInput('pc_x', label = 'x-axis PC: ', choices = 1:8, selected = 1),
                 selectInput('pc_y', label = 'y-axis PC: ', choices = 1:8, selected = 2),
                 uiOutput("color_by"),
#                  selectInput('color_by', label = 'color by: ',
#                              choices = c(NULL, poss_covars()), selected = NULL,multiple = T),
                 numericInput('pca_nrgenes', label = 'Nr of (most variant) genes:', value = 300,min = 50,max = 20000),
                 numericInput('pca_point_alpha', label = 'alpha: ', value = 1,min = 0,max = 1,step = 0.01),
                 numericInput('pca_label_size', label = 'Labels size: ', value = 2,min = 1,max = 8),
                 numericInput('pca_point_size', label = 'Points size: ', value = 2,min = 1,max = 8),
                 numericInput('pca_varname_size', label = 'Varname size: ', value = 4,min = 1,max = 8),
                 numericInput('pca_scale_arrow', label = 'Scaling factor : ', value = 1,min = 0.01,max = 10)
                 # TODO custom color palette? see icobra's proposal

                 ),
        menuItem("Plot settings", icon = icon("paint-brush"),
                 # bstooltip: Use the widgets below to setup general parameters for exporting produced plots

                 # selectInput("col_palette","Color palette",choices = list("hue"=hue_pal()(10),"set1"=brewer_pal(palette = "Set1")(9),"rainbow"=rainbow(10)))  ## TODO: need to work on the palette selector: maybe use selectize? see other examples

                 numericInput("export_width",label = "Width of exported figures (cm)",value = 30,min = 2),
                 numericInput("export_height",label = "Height of exported figures (cm)",value = 30,min = 2),


                 # tooltips explanation to have less crowded ui and still good docu
                 shinyBS::bsTooltip(
                   "export_width",
                   paste0("Use the widgets below to setup general parameters for exporting produced plots"),
                   "right", options = list(container = "body"))


                 )
      ),

      dashboardBody(

        ## Define output size of error messages
        tags$head(
          tags$style(HTML("
                          .shiny-output-error-validation {
                          font-size: 15px;
                          color: forestgreen;
                          text-align: center;
                          }
                          "))
          ),

        tabBox(
          width=12,

          tabPanel(
            "About",
            includeMarkdown(system.file("extdata", "about.md",package = "pcaExplorer")),
            hr(),
            shiny::verbatimTextOutput("showuploaded1"),
            shiny::verbatimTextOutput("showuploaded2"),
            shiny::verbatimTextOutput("showuploaded3"),
            shiny::verbatimTextOutput("showuploaded4"),
            footer()
            ),

          tabPanel(
            "Instructions",
            includeMarkdown(system.file("extdata", "instructions.md",package = "pcaExplorer"))
            ),

          tabPanel(
            "Data Preview",
            h1("Here will go some head of the count dataset, the samples design/covariates and so"),

            h3("General information on the provided SummarizedExperiment/DESeqDataSet"),
            shiny::verbatimTextOutput("showdata"),
            h3("Available metadata"),
            DT::dataTableOutput("showcoldata"),
            h3("Number of million of reads per sample"),
            plotOutput("reads_barplot"),
            h3("Basic summary for the counts"),
            verbatimTextOutput("reads_summary"),





            footer()

            ),

          tabPanel(
            "Samples View",
            p(h3('principal component analysis'), "PCA projections of sample abundances onto any pair of components."),
            fluidRow(checkboxInput("sample_labels","Display sample labels",value = TRUE)),
            fluidRow(
              column(
                width = 6,
                plotOutput('samples_pca',brush = "pca_brush"),
                div(align = "right", style = "margin-right:15px; margin-bottom:10px",
                    downloadButton("download_samplesPca", "Download Plot"),
                    textInput("filename_samplesPca",label = "Save as...",value = "samplesPca.pdf"))
                ),
              column(
                width= 6,
                plotOutput("samples_scree"),
                div(align = "right", style = "margin-right:15px; margin-bottom:10px",
                    downloadButton("download_samplesScree", "Download Plot"),
                    textInput("filename_samplesScree",label = "Save as...",value = "samplesScree.pdf")),
                fluidRow(
                  column(
                    width = 6,
                    radioButtons("scree_type","Scree plot type:",choices=list("Proportion of explained variance"="pev","Cumulative proportion of explained variance"="cev"),"pev")
                    ),
                  column(
                    width = 6,
                    numericInput("scree_pcnr","Number of PCs to display",value=8,min=2)
                    )
                  )
              )
            ),
            hr(),
            fluidRow(
              column(
                width = 6,
                plotOutput("samples_pca_zoom")
              ),
              column(
                width = 6,
                numericInput("ntophiload", "Nr of genes to display (top & bottom)",value = 10, min = 1, max=40),
                plotOutput("geneshiload")
              )
            )
          ),

          tabPanel(
            "Genes View",
            p(h3('principal component analysis'), "PCA projections of sample abundances onto any pair of components."),

            shinyURL.ui(),

            fluidRow(checkboxInput("variable_labels","Display variable labels",value = TRUE)),
            fluidRow(
              column(
                width = 4,
                h4("Main Plot - interact!"),
                plotOutput('genes_biplot',brush = 'pcagenes_brush',click="pcagenes_click"),
                div(align = "right", style = "margin-right:15px; margin-bottom:10px",
                    downloadButton("download_genesPca", "Download Plot"),
                    textInput("filename_genesPca",label = "Save as...",value = "genesPca.pdf"))),
              column(
                width = 4,
                h4("Zoomed window"),
                plotOutput("genes_biplot_zoom",click="pcagenes_zoom_click"),
                div(align = "right", style = "margin-right:15px; margin-bottom:10px",
                    downloadButton("download_genesZoom", "Download Plot"),
                    textInput("filename_genesZoom",label = "Save as...",value = "genesPca_zoomed.pdf"))),
              column(
                width = 4,
                h4("Boxplot of selected gene"),
                plotOutput("genes_biplot_boxplot"))),

            fluidRow(
              column(
                width = 6,
                h4("Zoomed heatmap"),
                plotOutput("heatzoom"),
                div(align = "right", style = "margin-right:15px; margin-bottom:10px",
                    downloadButton("download_genesHeatmap","Download Plot"),
                    textInput("filename_genesHeatmap",label = "Save as...",value = "genesHeatmap.pdf"))),
            column(
              width = 6,
              h4("Zoomed interactive heatmap"),
              fluidRow(radioButtons("heatmap_colv","Cluster samples",choices = list("Yes"=TRUE,"No"=FALSE),selected = TRUE)),
              fluidRow(d3heatmapOutput("heatzoomd3")))),
            hr(),
            fluidRow(
              column(
                width = 6,
                h4("Points selected by brushing - clicking and dragging:"),
                DT::dataTableOutput("pca_brush_out"),
                downloadButton('downloadData_brush', 'Download brushed points'),
                textInput("brushedPoints_filename","File name...")),
              column(
                width = 6,
                h4("Points selected by clicking:"),
                DT::dataTableOutput("pca_click_out"),
                textInput("clickedPoints_filename","File name..."),
                downloadButton('downloadData_click', 'Download clicked (or nearby) points'))
              )
            ),

          tabPanel(
            "Gene finder",
            fluidRow(
              h1("GeneFinder"),
              textInput("genefinder",label = "type in the name of the gene to search",value = NULL),
              verbatimTextOutput("searchresult"),
              verbatimTextOutput("debuggene"),
              checkboxInput("ylimZero","Set y axis limit to 0",value=TRUE),
              plotOutput("genefinder_plot"))
            ),

          tabPanel(
            "PCA2GO",
            h1("Functions enriched in the genes with high loadings on the selected principal components"),
            verbatimTextOutput("enrichinfo"),
            checkboxInput("compact_pca2go","Display compact tables",value=TRUE),
            fluidRow(
              column(width = 3),
              column(
                width = 6,
                DT::dataTableOutput("dt_pcver_pos")),
              column(width = 3)
            ),

            fluidRow(
              column(4,
                     DT::dataTableOutput("dt_pchor_neg")),
              column(4,
                     plotOutput("pca2go")),
              column(4,
                     DT::dataTableOutput("dt_pchor_pos"))
              ),
            fluidRow(
              column(width = 3),
              column(
                width = 6,
                DT::dataTableOutput("dt_pcver_neg")),
              column(width = 3)
            )
          ),



          tabPanel("Multifactor exploration",
                   fluidRow(
                     column(
                       width = 6,
                       uiOutput("covar1")
                     ),
                     column(
                       width = 6,
                       uiOutput("covar2")
                     )
                   ),
                   fluidRow(
                     column(
                       width = 6,
                       uiOutput("c1levels")
                     ),
                     column(
                       width = 6,
                       uiOutput("c2levels")
                     )
                   ),
                   fluidRow(
                     column(
                       width = 6,
                       uiOutput("colnames1")
                     ),
                     column(
                       width = 6,
                       uiOutput("colnames2")
                     )
                   ),

#
#                    uiOutput("c1levels"),
#                    uiOutput("c2levels"),
#
#                    uiOutput("colnames1"),
#                    uiOutput("colnames2"),

                   actionButton("composemat","Compose the matrix",icon=icon("spinner")),

                   verbatimTextOutput("dd"),

                   fluidRow(
                     column(4,
                            selectInput('pc_x_ruf', label = 'x-axis PC: ', choices = 1:8,
                                        selected = 2)
                     ),
                     column(4,
                            selectInput('pc_y_ruf', label = 'y-axis PC: ', choices = 1:8,
                                        selected = 3)
                     )),

                   # fluidRow(verbatimTextOutput("rufdebug")),

                   fluidRow(
                     column(6,
                            plotOutput('pcaruf',brush = 'pcaruf_brush')),
                     column(6,
                            plotOutput("rufzoom"))
                   ),
                   fluidRow(downloadButton('downloadData_brush_ruf', 'Download brushed points'),
                            textInput("brushedPoints_filename_ruf","File name..."),
                            DT::dataTableOutput('pcaruf_out'))


          )

        )

      ),

      skin="blue"

    )

  ## ------------------------------------------------------------------ ##
  ##                          Define server                             ##
  ## ------------------------------------------------------------------ ##

  newserver <- shinyServer(function(input, output, session) {

    exportPlots <- reactiveValues(
      samplesPca=NULL,
      samplesZoom=NULL,
      samplesScree=NULL,
      genesPca=NULL,
      genesZoom=NULL,
      genesBoxplot=NULL,
      genesHeatmap=NULL,
      genefinder=NULL
    )

    values <- reactiveValues()
    values$mydds <- obj2
    values$myrlt <- obj
    values$mycountmatrix <- countmatrix
    values$mymetadata <- coldata
    values$mypca2go <- pca2go
    values$myannotation <- annotation

    user_settings <- reactiveValues(save_width = 45, save_height = 11)

    shinyURL.server()

    # compute only rlt if dds is provided but not cm&coldata
    if(!is.null(obj2) & (is.null(countmatrix) & is.null(coldata)) & is.null(obj))
      withProgress(message = "computing rlog transformed values...",
                   value = 0,
                   {
                     values$myrlt <- rlogTransformation(obj2)
                   })


    output$color_by <- renderUI({
      if(is.null(values$mydds))
        return(NULL)
      poss_covars <- names(colData(values$mydds))
      selectInput('color_by', label = 'color by: ',
                  choices = c(NULL, poss_covars), selected = NULL,multiple = TRUE)
    })



    ## Render the UI element to upload the count matrix
    output$upload_count_matrix <- renderUI({
      if (!is.null(obj2) | !is.null(countmatrix)) {
        NULL
      } else {
        return(fileInput(inputId = "uploadcmfile",
                         label = "Upload a count matrix file",
                         accept = c("text/csv", "text/comma-separated-values",
                                    "text/tab-separated-values", "text/plain",
                                    ".csv", ".tsv"), multiple = FALSE))
      }
    })

    readCountmatrix <- reactive({
      if (is.null(input$uploadcmfile))
        return(NULL)
      cm <- utils::read.delim(input$uploadcmfile$datapath, header = TRUE,
                              as.is = TRUE, sep = "\t", quote = "",
                              check.names = FALSE)

      return(cm)
    })




    output$upload_metadata <- renderUI({
      if (!is.null(obj2) | !is.null(coldata)) {
        NULL
      } else {
        return(fileInput(inputId = "uploadmetadatafile",
                         label = "Upload a sample metadata matrix file",
                         accept = c("text/csv", "text/comma-separated-values",
                                    "text/tab-separated-values", "text/plain",
                                    ".csv", ".tsv"), multiple = FALSE))
      }
    })

    readMetadata <- reactive({
      if (is.null(input$uploadmetadatafile))
        return(NULL)
      coldata <- utils::read.delim(input$uploadmetadatafile$datapath, header = TRUE,
                              as.is = TRUE, sep = "\t", quote = "",
                              check.names = FALSE)

      return(coldata)
    })


    output$upload_annotation <- renderUI({
      if (!is.null(annotation)) {
        NULL
      } else {
        return(fileInput(inputId = "uploadannotationfile",
                         label = "Upload an annotation file",
                         accept = c("text/csv", "text/comma-separated-values",
                                    "text/tab-separated-values", "text/plain",
                                    ".csv", ".tsv"), multiple = FALSE))
      }
    })

    readAnnotation <- reactive({
      if (is.null(input$uploadannotationfile))
        return(NULL)
      annodata <- utils::read.delim(input$uploadannotationfile$datapath, header = TRUE,
                                   as.is = TRUE, sep = "\t", quote = "",
                                   check.names = FALSE)

      return(annodata)
    })


    createDDS <- reactive({
      if(is.null(countmatrix) | is.null(coldata))
        return(NULL)

      dds <- DESeqDataSetFromMatrix(countData = countmatrix,
                                             colData = coldata,
                                             design=~1)

      return(dds)


    })

    createRLT <- reactive({
      if(is.null(countmatrix) | is.null(coldata))
        return(NULL)

      rlt <- rlogTransformation(values$mydds)

      return(rlt)


    })

    observeEvent(createDDS,
                 {
                   if(!is.null(values$mycountmatrix) & !is.null(values$mymetadata))
                     values$mydds <- createDDS()
                 })

    observeEvent(createRLT,
                 {
                   if(!is.null(values$mycountmatrix) & !is.null(values$mymetadata))
                     values$myrlt <- createRLT()
                 })




    # as in http://stackoverflow.com/questions/29716868/r-shiny-how-to-get-an-reactive-data-frame-updated-each-time-pressing-an-actionb
    observeEvent(input$uploadcmfile,
                 {
                   values$mycountmatrix <- readCountmatrix()
                   if(!is.null(values$mymetadata)){
                     withProgress(message="Computing the objects...",value = 0,{

                     values$mydds <- DESeqDataSetFromMatrix(countData = values$mycountmatrix,
                                                            colData = values$mymetadata,
                                                            design=~1)
                     values$myrlt <- rlogTransformation(values$mydds)})
                   }
                 })

    observeEvent(input$uploadmetadatafile,
                 {
                   values$mymetadata <- readMetadata()
                   if(!is.null(values$mycountmatrix)){
                     withProgress(message="Computing the objects...",value = 0,{

                     values$mydds <- DESeqDataSetFromMatrix(countData = values$mycountmatrix,
                                                            colData = values$mymetadata,
                                                            design=~1)
                     values$myrlt <- rlogTransformation(values$mydds)})
                   }
                 })

    observeEvent(input$uploadannotationfile,
                 {
                   values$myannotation <- readAnnotation()
                 })



    output$showuploaded1 <- renderPrint({
      head(values$mycountmatrix)
    })
    output$showuploaded2 <- renderPrint({
      values$mymetadata
    })
    output$showuploaded3 <- renderPrint({
      values$mydds
    })
    output$showuploaded4 <- renderPrint({
      values$myrlt
    })


    colSel <- reactive({
      hue_pal()(ncol(values$myrlt)/2) # or somewhat other way
    })




    output$showdata <- renderPrint({
      values$mydds
    })

    output$showcoldata <- DT::renderDataTable({
      datatable(as.data.frame(colData(values$mydds)))
    })

    output$reads_barplot <- renderPlot({
      rr <- colSums(counts(values$mydds))/1e6
      if(is.null(names(rr)))
        names(rr) <- paste0("sample_",1:length(rr))
      rrdf <- data.frame(Reads=rr,Sample=names(rr),stringsAsFactors = FALSE)
      if (!is.null(input$color_by)) {
        rrdf$Group <- colData(values$mydds)[input$color_by][[1]]
        p <- ggplot(rrdf,aes(Sample,weight=Reads)) + geom_bar(aes(fill=Group))
        p
      } else {
        p <- ggplot(rrdf,aes(Sample,weight=Reads)) + geom_bar()
        p
      }
    })

    output$reads_summary <- renderPrint({
      summary(colSums(counts(values$mydds))/1e6)
    })








    output$samples_pca <- renderPlot({
      res <- pcaplot(values$myrlt,intgroup = input$color_by,ntop = input$pca_nrgenes,
                     pcX = as.integer(input$pc_x),pcY = as.integer(input$pc_y),
                     text_labels = input$sample_labels,
                     point_size = input$pca_point_size, title="Samples PCA"
      )
      res <- res + theme_bw()
      exportPlots$samplesPca <- res
      res
    })

    output$samples_pca_zoom <- renderPlot({

      shiny::validate(
        need(!is.null(input$pca_brush),
             "Zoom in by brushing in the main plot panel above"
        )
      )
      # if(is.null(input$pca_brush))
        # return(ggplot() + annotate("text",label="zoom in by brushing",0,0) + theme_bw())

      res <- pcaplot(values$myrlt,intgroup = input$color_by,ntop = input$pca_nrgenes,
                     pcX = as.integer(input$pc_x),pcY = as.integer(input$pc_y),
                     text_labels = input$sample_labels,
                     point_size = input$pca_point_size, title="Samples PCA - zoom in"
)
      res <- res + xlim(input$pca_brush$xmin,input$pca_brush$xmax) + ylim(input$pca_brush$ymin,input$pca_brush$ymax)
      res <- res + theme_bw()
      exportPlots$samplesZoom <- res
      res
    })

    output$samples_scree <- renderPlot({
      rv <- rowVars(assay(values$myrlt))
      select <- order(rv, decreasing = TRUE)[seq_len(min(input$pca_nrgenes,length(rv)))]
      pca <- prcomp(t(assay(values$myrlt)[select, ]))

      res <- pcascree(pca,type = input$scree_type, pc_nr = input$scree_pcnr, title="Scree plot for the samples PCA")
      res <- res + theme_bw()
      exportPlots$samplesScree <- res
      res
    })


    output$geneshiload <- renderPlot({
      rv <- rowVars(assay(values$myrlt))
      select <- order(rv, decreasing = TRUE)[seq_len(min(input$pca_nrgenes,length(rv)))]
      pca <- prcomp(t(assay(values$myrlt)[select, ]))

      par(mfrow=c(2,1))
      hi_loadings(pca,whichpc = as.integer(input$pc_x),topN = input$ntophiload,annotation = values$myannotation)
      hi_loadings(pca,whichpc = as.integer(input$pc_y),topN = input$ntophiload,annotation = values$myannotation)

    })






    output$genes_biplot <- renderPlot({
      if(!is.null(input$color_by)) {
        expgroups <- as.data.frame(colData(values$myrlt)[,input$color_by])
        expgroups <- interaction(expgroups)
      } else {
        expgroups <- colnames(values$myrlt)
      }
      colGroups <- colSel()[factor(expgroups)]

      res <- genepca(values$myrlt,
                     ntop = input$pca_nrgenes,
                     choices = c(as.integer(input$pc_x),as.integer(input$pc_y)),
                     biplot = TRUE,
                     arrowColors = factor(colGroups),
                     groupNames = expgroups,
                     alpha=input$pca_point_alpha,coordEqual=FALSE,useRownamesAsLabels=FALSE,labels.size=input$pca_label_size,
                     point_size=input$pca_point_size,varname.size=input$pca_varname_size, scaleArrow = input$pca_scale_arrow,annotation=values$myannotation)
      exportPlots$genesPca <- res
      res
    })


    output$genes_biplot_zoom <- renderPlot({
      # if(is.null(input$pcagenes_brush)) return(ggplot() + annotate("text",label="zoom in by brushing",0,0) + theme_bw())

      shiny::validate(
        need(
          !is.null(input$pcagenes_brush),
          "Zoom in by brushing in the main panel - this will also allow displaying the gene names"
        )
      )

      if(!is.null(input$color_by)) {
        expgroups <- as.data.frame(colData(values$myrlt)[,input$color_by])
        expgroups <- interaction(expgroups)
      } else {
        expgroups <- colnames(values$myrlt)
      }
      colGroups <- colSel()[factor(expgroups)]

      res <- genepca(values$myrlt,
                     ntop = input$pca_nrgenes,
                     choices = c(as.integer(input$pc_x),as.integer(input$pc_y)),
                     biplot = TRUE,
                     arrowColors = factor(colGroups),
                     groupNames = expgroups,
                     alpha=input$pca_point_alpha,coordEqual=FALSE,
                     var.axes=input$variable_labels, # workaround for a ggplot2 bug/missing thing: here details: https://github.com/hadley/ggplot2/issues/905
                     labels.size=input$pca_label_size,varname.size=input$pca_varname_size,
                     scaleArrow = input$pca_scale_arrow,point_size=input$pca_point_size,annotation=values$myannotation)

      res <- res +
        xlim(input$pcagenes_brush$xmin,input$pcagenes_brush$xmax) +
        ylim(input$pcagenes_brush$ymin,input$pcagenes_brush$ymax)
      exportPlots$genesZoom <- res
      res
    })


    output$genes_biplot_boxplot <- renderPlot({
      # if(length(input$color_by)==0) return(ggplot() + annotate("text",label="select an experimental factor",0,0) + theme_bw())
      # if(is.null(input$pcagenes_zoom_click)) return(ggplot() + annotate("text",label="click to generate the boxplot\nfor the selected gene",0,0) + theme_bw())

      shiny::validate(
        need(
          length(input$color_by)>0,
          "Select an experimental factor"
        )
      )
      shiny::validate(
        need(
          !is.null(input$pcagenes_zoom_click),
          "Click to generate the boxplot for the selected gene"
        )
      )

      selectedGene <- curData_zoomClick()$ids
      selectedGeneSymbol <- values$myannotation$gene_name[match(selectedGene,rownames(values$myannotation))]
      # plotCounts(dds_cleaner,)

      shiny::validate(
        need(nrow(curData_zoomClick()) >0,message = "Click closer to a gene to get the boxplot")

        )

      genedata <- plotCounts(values$mydds,gene=selectedGene,intgroup = input$color_by,returnData = TRUE)

      onlyfactors <- genedata[,match(input$color_by,colnames(genedata))]
      genedata$plotby <- interaction(onlyfactors)

      res <- ggplot(genedata,aes(x=plotby,y=count,fill=plotby)) +
        geom_boxplot(outlier.shape = NA) + scale_y_log10(name="Normalized counts") +
        labs(title=paste0("Normalized counts for ",selectedGeneSymbol," - ",selectedGene)) +
        scale_x_discrete(name="") +
        geom_jitter(aes(x=plotby,y=count),position = position_jitter(width = 0.1)) +
        scale_fill_discrete(name="Experimental\nconditions")
      exportPlots$genesBoxplot <- res
      res
    })


    # for reading in the brushed/clicked points
    curData_brush <- reactive({
      df2 <- genepca(values$myrlt,
                     ntop = input$pca_nrgenes,
                     choices = c(as.integer(input$pc_x),as.integer(input$pc_y)),
                     biplot = TRUE,
                     # arrowColors = colGroups,
                     alpha=input$pca_point_alpha,
                     returnData=TRUE,annotation=values$myannotation)
      df2$geneName <- values$myannotation$gene_name[match(rownames(df2),rownames(values$myannotation))]
      res <- brushedPoints(df2, input$pcagenes_brush,xvar="xvar",yvar="yvar",)
      res
    })


    curData_click <- reactive({
      df2 <- genepca(values$myrlt,
                     ntop = input$pca_nrgenes,
                     choices = c(as.integer(input$pc_x),as.integer(input$pc_y)),
                     biplot = TRUE,
                     # arrowColors = colGroups,
                     alpha=input$pca_point_alpha,
                     returnData=TRUE,annotation=values$myannotation)
      df2$geneName <- values$myannotation$gene_name[match(rownames(df2),rownames(values$myannotation))]
      res <- nearPoints(df2, input$pcagenes_click,
                        threshold = 20, maxpoints = 3,
                        addDist = TRUE)
      # res <- brushedPoints(df2, input$pcagenes_brush,xvar="xvar",yvar="yvar",)
      res
    })


    curData_zoomClick <- reactive({
      df2 <- genepca(values$myrlt,
                     ntop = input$pca_nrgenes,
                     choices = c(as.integer(input$pc_x),as.integer(input$pc_y)),
                     biplot = TRUE,
                     # arrowColors = colGroups,
                     alpha=input$pca_point_alpha,
                     returnData=TRUE,annotation=values$myannotation)
      df2$geneName <- values$myannotation$gene_name[match(rownames(df2),rownames(values$myannotation))]
      res <- nearPoints(df2, input$pcagenes_zoom_click,
                        threshold = 20, maxpoints = 1,
                        addDist = TRUE)
      # res <- brushedPoints(df2, input$pcagenes_brush,xvar="xvar",yvar="yvar",)
      res
    })



    output$pca_brush_out <- DT::renderDataTable({
      datatable(curData_brush(),options = list(pageLength = 50))
    })

    output$pca_click_out <- DT::renderDataTable({
      datatable(curData_click(),options = list(pageLength = 50))
    })




    output$heatzoomd3 <- renderD3heatmap({
      shiny::validate(
        need(
          !is.null(input$pcagenes_brush),
          "Brush the main panel above to generate a heatmap"
        )
      )

      # if(is.null(input$pcagenes_brush)) return(NULL)

      brushedObject <- curData_brush()
      shiny::validate(
        need(
          nrow(brushedObject) > 1,
          "Brush to include at least two genes"
        )
      )

      selectedGenes <- brushedObject$ids
      toplot <- assay(values$myrlt)[selectedGenes,]
      rownames(toplot) <- values$myannotation$gene_name[match(rownames(toplot),rownames(values$myannotation))]

      mycolss <- c("#313695","#4575b4","#74add1","#abd9e9","#e0f3f8","#fee090","#fdae61","#f46d43","#d73027","#a50026") # to be consistent with red/blue usual coding

      d3heatmap(toplot,Colv = as.logical(input$heatmap_colv),colors = mycolss)
    })


    output$heatzoom <- renderPlot({
      # if(is.null(input$pcagenes_brush)) return(NULL)
      shiny::validate(
        need(
          !is.null(input$pcagenes_brush),
          "Brush the main panel above to generate a heatmap"
        )
      )

      brushedObject <- curData_brush()
      shiny::validate(
        need(
          nrow(brushedObject) > 1,
          "Brush to include at least two genes"
        )
      )
      selectedGenes <- brushedObject$ids
      toplot <- assay(values$myrlt)[selectedGenes,]
      rownames(toplot) <- values$myannotation$gene_name[match(rownames(toplot),rownames(values$myannotation))]
      # pheatmap(toplot,cluster_cols = as.logical(input$heatmap_colv))
      NMF::aheatmap(toplot,Colv = as.logical(input$heatmap_colv))
      ## aheatmap is actually consistent in displaying the clusters with most of other heatmap packages
      ## keep in mind: pheatmap does somehow a better job if scaling/centering
    })








    output$searchresult <- renderPrint({

      if(is.null(input$color_by)) return("Select a factor to plot your gene")
      if(input$genefinder=="")
        return("Type in the gene name/id you want to plot")

      foundGeneID <- input$genefinder %in% rownames(values$myrlt)
      foundGeneName <- input$genefinder %in% values$myannotation$gene_name
      if(!foundGeneID){
        foundGeneID <- toupper(input$genefinder) %in% toupper(rownames(values$myrlt))
        if(foundGeneID){
          return(paste0("Maybe you mis-spelled the name of your gene. Did you mean ",
                        unique(rownames(values$myannotation)[which(toupper(input$genefinder)==toupper(rownames(values$myannotation)))]),"?"))
        } else {
          foundGeneNAME <- input$genefinder %in% values$myannotation$gene_name
          if(!foundGeneNAME){
            foundGeneNAME <- toupper(input$genefinder) %in% toupper(values$myannotation$gene_name)
            if(foundGeneNAME){
              return(paste0("Maybe you mis-spelled the name of your gene. Did you mean ",
                            unique(values$myannotation$gene_name[which(toupper(input$genefinder)==toupper(values$myannotation$gene_name))]),"?"))
            } else {return("Could not find the gene you typed!")}
          } else {
            fgn <- values$myannotation$gene_name[which(values$myannotation$gene_name==input$genefinder)]
            if (length(fgn) > 1) return(paste0("Found more than one gene with the selected gene name. Select one of the following: ",paste(selectedGene,collapse=", ")))
            selectedGene <- rownames(values$myannotation)[which(values$myannotation$gene_name==input$genefinder)]

            fg <- rownames(values$myannotation)[match(fgn,values$myannotation$gene_name)]
            return(paste0("I found the gene! Plotting ", fg, " - ", values$myannotation$gene_name[match(fg,rownames(values$myannotation))],"..."))

          }}
      } else {
        fg <- rownames(values$myannotation)[match(input$genefinder,rownames(values$myrlt))]
        return(paste0("I found the gene! Plotting ", fg, " - ", values$myannotation$gene_name[match(fg,rownames(values$myannotation))],"..."))

      }
    })





    output$genefinder_plot <- renderPlot({
      anno_id <- rownames(values$myannotation)
      anno_gene <- values$myannotation$gene_name

      if(is.null(input$color_by) & input$genefinder!="")
        return(ggplot() + annotate("text",label="Select a factor to plot your gene",0,0) + theme_bw())
      if(is.null(input$color_by) & input$genefinder=="")
        return(ggplot() + annotate("text",label="Select a gene and a factor to plot gene",0,0) + theme_bw())
      if(input$genefinder=="")
        return(ggplot() + annotate("text",label="Type in a gene name/id",0,0) + theme_bw())
      if(!input$genefinder %in% anno_gene & !input$genefinder %in% anno_id)
        return(ggplot() + annotate("text",label="Gene not found...",0,0) + theme_bw())

      if (input$genefinder %in% anno_id) {
        selectedGene <- rownames(values$myrlt)[match(input$genefinder,rownames(values$myrlt))]
        selectedGeneSymbol <- values$myannotation$gene_name[match(selectedGene,rownames(values$myannotation))]
      }
      if (input$genefinder %in% anno_gene) {
        selectedGeneSymbol <- values$myannotation$gene_name[which(values$myannotation$gene_name==input$genefinder)]
        if (length(selectedGeneSymbol) > 1) return(ggplot() + annotate("text",label=paste0("Type in a gene name/id of the following:\n",paste(selectedGene,collapse=", ")),0,0) + theme_bw())
        selectedGene <- rownames(values$myannotation)[which(values$myannotation$gene_name==input$genefinder)]
      }
      genedata <- plotCounts(values$mydds,gene=selectedGene,intgroup = input$color_by,returnData = TRUE)

      onlyfactors <- genedata[,match(input$color_by,colnames(genedata))]
      genedata$plotby <- interaction(onlyfactors)

      p <- ggplot(genedata,aes(x=plotby,y=count,fill=plotby)) + geom_boxplot() + labs(title=paste0("Normalized counts for ",selectedGeneSymbol," - ",selectedGene)) +  scale_x_discrete(name="") + geom_jitter(aes(x=plotby,y=count),position = position_jitter(width = 0.1)) + scale_fill_discrete(name="Experimental\nconditions")

      if(input$ylimZero)
      {
        p <- p + scale_y_log10(name="Normalized counts - log10 scale",limits=c(1,max(genedata$count)))
      } else {
        p <- p + scale_y_log10(name="Normalized counts - log10 scale")
      }
      exportPlots$genefinder <- p

      p
    })






    output$pca2go <- renderPlot({
      shiny::validate(
        need(
          !is.null(pca2go),
          "Please provide a pca2go object to the app"
        )
      )
      # if(is.null(pca2go))
        # return(ggplot() + annotate("text",label="Provide a pca2go object to the app",0,0) + theme_bw())
      res <- pcaplot(values$myrlt,intgroup = input$color_by,
                     ntop = attr(pca2go,"n_genesforpca"),
                     pcX = as.integer(input$pc_x),pcY = as.integer(input$pc_y),text_labels = input$sample_labels,
                     point_size = input$pca_point_size, title=paste0("PCA on the samples - ",attr(pca2go,"n_genesforpca"), " genes used")

      )
      res
    })


    output$dt_pchor_pos <- DT::renderDataTable({
      if(is.null(pca2go)) return(datatable(NULL))
      goe <- pca2go[[paste0("PC",input$pc_x)]][["posLoad"]]
      if(input$compact_pca2go)
        return(datatable(goe[,c("GO.ID","Term","Significant","p.value_elim")],options = list(pageLength = 5)))
      datatable(goe)
    })

    output$dt_pchor_neg <- DT::renderDataTable({
      if(is.null(pca2go)) return(datatable(NULL))
      goe <- pca2go[[paste0("PC",input$pc_x)]][["negLoad"]]
      if(input$compact_pca2go)
        return(datatable(goe[,c("GO.ID","Term","Significant","p.value_elim")],options = list(pageLength = 5)))
      datatable(goe)
    })

    output$dt_pcver_pos <- DT::renderDataTable({
      if(is.null(pca2go)) return(datatable(NULL))
      goe <- pca2go[[paste0("PC",input$pc_y)]][["posLoad"]]
      if(input$compact_pca2go)
        return(datatable(goe[,c("GO.ID","Term","Significant","p.value_elim")],options = list(pageLength = 5)))
      datatable(goe)
    })

    output$dt_pcver_neg <- DT::renderDataTable({
      if(is.null(pca2go)) return(datatable(NULL))
      goe <- pca2go[[paste0("PC",input$pc_y)]][["negLoad"]]
      if(input$compact_pca2go)
        return(datatable(goe[,c("GO.ID","Term","Significant","p.value_elim")],options = list(pageLength = 5)))
      datatable(goe)
    })

    output$enrichinfo <- renderPrint({
      cat("enrich info:\n")
      # str(goEnrichs)
      class(input$pc_x)
      head(pca2go[[paste0("PC",input$pc_x)]][["posLoad"]])
      class(datatable(pca2go[[paste0("PC",input$pc_x)]][["posLoad"]]))
    })














    ## from here on, RUF APP


    output$covar1 <- renderUI({
      # if(is.null(values$myrlt))
        # return(NULL)
      poss_covars <- names(colData(values$mydds))
      selectInput('covar1', label = 'factor1: ',
                  choices = c(NULL, poss_covars), selected = NULL,multiple = F)
    })

    output$covar2 <- renderUI({
      # if(is.null(values$myrlt))
      # return(NULL)
      poss_covars <- names(colData(values$mydds))
      selectInput('covar2', label = 'factor2: ',
                  choices = c(NULL, poss_covars), selected = NULL,multiple = F)
    })


    output$c1levels <- renderUI({
      if(is.null(input$covar1))
        return(NULL)
      fac1lev <- levels(colData(values$myrlt)[[input$covar1]])
      selectInput('covar1levels', label = 'factor1 levels: ',
                  choices = c(NULL, fac1lev), selected = NULL,multiple = T) # actually 2
    })

    output$c2levels <- renderUI({
      if(is.null(input$covar2))
        return(NULL)
      fac2lev <- levels(colData(values$myrlt)[[input$covar2]])
      selectInput('covar2levels', label = 'factor2 levels: ',
                  choices = c(NULL, fac2lev), selected = NULL,multiple = T) # 2 or more are allowed!
    })

    output$colnames1 <- renderUI({
      if(is.null(values$myrlt))
        return(NULL)
      if(is.null(input$covar1))
        return(NULL)
      if(is.null(input$covar2))
        return(NULL)

      fac1 <- input$covar1
      fac2 <- input$covar2

      fac1_touse <- input$covar1levels
      fac2_touse <- input$covar2levels

      preselected_fac1 <- colnames(values$myrlt)[colData(values$myrlt)[[fac1]] %in% fac1_touse]
      preselected_fac2 <- colnames(values$myrlt)[colData(values$myrlt)[[fac2]] %in% fac2_touse]
      presel <- intersect(preselected_fac1,preselected_fac2)
      mysamples <- colData(values$myrlt)[presel,] # check that the repl are balanced

      presel1 <- colnames(values$myrlt)[(colData(values$myrlt)[[fac1]] %in% fac1_touse[1]) & colData(values$myrlt)[[fac2]] %in% fac2_touse]

      selectInput('picksamples1', label = 'combine these samples in the selected order: ',
                  choices = c(NULL, presel1), selected = NULL,multiple = TRUE)
    })


    output$colnames2 <- renderUI({
      if(is.null(values$myrlt))
        return(NULL)
      if(is.null(input$covar1))
        return(NULL)
      if(is.null(input$covar2))
        return(NULL)

      fac1 <- input$covar1
      fac2 <- input$covar2

      fac1_touse <- input$covar1levels
      fac2_touse <- input$covar2levels

      preselected_fac1 <- colnames(values$myrlt)[colData(values$myrlt)[[fac1]] %in% fac1_touse]
      preselected_fac2 <- colnames(values$myrlt)[colData(values$myrlt)[[fac2]] %in% fac2_touse]
      presel <- intersect(preselected_fac1,preselected_fac2)
      mysamples <- colData(values$myrlt)[presel,] # check that the repl are balanced

      presel2 <- colnames(values$myrlt)[(colData(values$myrlt)[[fac1]] %in% fac1_touse[2]) & colData(values$myrlt)[[fac2]] %in% fac2_touse]

      selectInput('picksamples2', label = 'combine these samples in the selected order: ',
                  choices = c(NULL, presel2), selected = NULL,multiple = TRUE)
    })





    # I WILL MODIFY HERE
#     ddsmf_clean
#     rld_global
#
#
#
#     kldds <- updateObject(dds_kl)
#     klrld <- updateObject(rld_kl)
#
#     ddsobj <- updateObject(ddsmf_clean)
#     rldobj <- updateObject(rld_global)


    composedMat <- eventReactive( input$composemat, {
      exprmat <- t(assay(values$myrlt))
      exprmat <- exprmat[,rowSums(counts(values$mydds) > 5)>2]

      withProgress(message = "Composing the matrix...",
                   value = 0,
                   {
                     pcmat <- cbind(exprmat[input$picksamples1,],
                                    exprmat[input$picksamples2,])
                   })
      pcmat
    })


    obj3 <- reactive({

      pcmat <- composedMat()

      library(scales)
      aval <- 0.3
      fac2pal <- alpha(c("green","red","blue","orange","violet"),aval) # 5 are enough

      # colData(values$myrlt)[input$covar2][rownames(pcmat),]
      max.type <- apply(pcmat[,1:(ncol(pcmat)/2)],2,which.max)

      fac2_col <- factor(colData(values$myrlt)[input$covar2][rownames(pcmat),],
                         levels=unique(as.character(colData(values$myrlt)[input$covar2][rownames(pcmat),])))
      tcol.justMax <- fac2pal[fac2_col][max.type]
      # tcol.justMax <- ifelse(max.type <= 4,"green",ifelse(max.type <= 8,"red",ifelse(max.type <= 12,"blue","orange")))

      max.type2 <- apply(pcmat[,((ncol(pcmat)/2)+1):ncol(pcmat)],2,which.max)
      # tcol2.justMax <- ifelse(max.type2 <= 4,alpha("green",aval),ifelse(max.type2 <= 8,alpha("red",aval),ifelse(max.type2 <= 12,alpha("blue",aval),alpha("orange",aval))))

      tcol2.justMax <- fac2pal[fac2_col][max.type2]



      # using the median across replicates
      celltypes <- gsub("_R.","",rownames(pcmat))

      # brutal way, with for loop...
      #         mediansByReps <- matrix(NA,4,ncol(pcmat))
      #         minipcmat <- pcmat[1:16,1:5]
      #         # mediansByReps <- matrix(NA,4,ncol(minipcmat))
      #
      #         starts <- c(1,5,9,13)
      #         for(j in 1:ncol(mediansByReps)){
      #           for(i in 1:4){
      #             tmp <- pcmat[(starts[i]:(starts[i]+3)),j]
      #             # cat(paste(tmp,collapse=","))
      #             res <- median(tmp)
      #             mediansByReps[i,j] <- res
      #           }
      #           # print(j)
      #         }
      #
      #         max.median.type <- apply(mediansByReps[,1:ncol(exprmat)],2,which.max)
      #         max.median.type2 <- apply(mediansByReps[,(ncol(exprmat)+1):ncol(mediansByReps)],2,which.max)
      #         tcol.median <- c(alpha("green",aval),alpha("red",aval),alpha("blue",aval),alpha("orange",aval))[max.median.type]
      #         tcol2.median <- c(alpha("green",aval),alpha("red",aval),alpha("blue",aval),alpha("orange",aval))[max.median.type2]
      #
      #         # optional step, in the end, to see if the median across replicates delivers different stuff
      #         tcol <- tcol.median
      #         tcol2 <- tcol2.median

      tcol <- tcol.justMax
      tcol2 <- tcol2.justMax
      # pcmat
      return(list(pcmat,tcol,tcol2))
    })


    output$pcaruf <- renderPlot({
    pcmat <- obj3()[[1]]
      tcol <- obj3()[[2]]
      tcol2 <- obj3()[[3]]
      pres <- prcomp(t(pcmat),scale=FALSE)

      plot.index <- c(as.integer(input$pc_x_ruf),as.integer(input$pc_y_ruf))
      offset <- ncol(pcmat)/2
      gene.no <- offset
      pcx <- pres$x
      # set.seed(11)
      # for (i in 1:ncol(pcx)) {
      #   pcx[,i] <- pcx[,i] + rnorm(nrow(pcx),sd=diff(range(pcx[,i]))/100)
      # }
      plot(pcx[(offset+1):ncol(pcmat),plot.index[1]][1:gene.no],pcx[(offset+1):ncol(pcmat),plot.index[2]][1:gene.no],xlim=range(pcx[,plot.index[1]]),ylim=range(pcx[,plot.index[2]]),pch=20,col=tcol,cex=0.3)#,type="n")
      #plot(0,type="n",xlim=range(pres$x[,plot.index]),ylim=range(pres$x[,plot.index]))
      lcol <- ifelse(tcol != tcol2,"black","grey")
      for (i in 1:gene.no) {
        lines(pcx[c(i,offset+i),plot.index[1]],pcx[c(i,offset+i),plot.index[2]],col=lcol[i])
      }
      points(pcx[1:offset,plot.index[1]][1:gene.no],pcx[1:offset,plot.index[2]][1:gene.no],pch=20,col=tcol,cex=0.3)
      points(pcx[(offset+1):ncol(pcmat),plot.index[1]][1:gene.no],pcx[(offset+1):ncol(pcmat),plot.index[2]][1:gene.no],pch=20,col=tcol2,cex=0.3)

#       ## ## ##
#       #         points(pcx[rownames(pcx) %in% rownames(cm2)[match(galon_macro_mouse,cm2$fromgtf)],plot.index[1]],
#       #                pcx[rownames(pcx) %in% rownames(cm2)[match(galon_macro_mouse,cm2$fromgtf)],plot.index[2]],
#       #                pch=20,col="darkviolet",cex=2)
#       #         points(pcx[rownames(pcx) %in% rownames(cm2)[match(galon_cd8_mouse,cm2$fromgtf)],plot.index[1]],
#       #                pcx[rownames(pcx) %in% rownames(cm2)[match(galon_cd8_mouse,cm2$fromgtf)],plot.index[2]],
#       #                pch=20,col="steelblue",cex=2)
#       #         # legend("topleft",fill = c("darkviolet","steelblue"),legend=c("galon Macro","galon CD8"))
#       mgenes_extended <- c("Tnf","Mmp12","Adam8","Mrc1","Cd36","Cd83","Itgam","Lyz1","Slamf8","Clec12a","Clec10a","Pdc","Cd274")
#       mgenes_extended_ENS <- rownames(cm2)[match(mgenes_extended,cm2$fromgtf)]
#       points(pcx[rownames(pcx) %in% mgenes_extended_ENS,plot.index[1]],
#              pcx[rownames(pcx) %in% mgenes_extended_ENS,plot.index[2]],
#              pch=20,col="salmon",cex=2)
      ## ## ##
    })

    output$dd <- renderPrint({
      input$picksamples1
    })


    output$rufzoom <- renderPlot({
      if(is.null(input$pcaruf_brush)) return(NULL)
      pcmat <- obj3()[[1]]
      tcol <- obj3()[[2]]
      tcol2 <- obj3()[[3]]
      pres <- prcomp(t(pcmat),scale=FALSE)

      plot.index <- c(as.integer(input$pc_x_ruf),as.integer(input$pc_y_ruf))
      offset <- ncol(pcmat)/2
      gene.no <- offset
      pcx <- pres$x
      # set.seed(11)
      # for (i in 1:ncol(pcx)) {
      #   pcx[,i] <- pcx[,i] + rnorm(nrow(pcx),sd=diff(range(pcx[,i]))/100)
      # }
      plot(pcx[(offset+1):ncol(pcmat),plot.index[1]][1:gene.no],
           pcx[(offset+1):ncol(pcmat),plot.index[2]][1:gene.no],
           xlim=c(input$pcaruf_brush$xmin,input$pcaruf_brush$xmax),
           ylim=c(input$pcaruf_brush$ymin,input$pcaruf_brush$ymax),
           pch=20,col=tcol,cex=0.3)#,type="n")
      #plot(0,type="n",xlim=range(pres$x[,plot.index]),ylim=range(pres$x[,plot.index]))
      lcol <- ifelse(tcol != tcol2,"black","grey")
      for (i in 1:gene.no) {
        lines(pcx[c(i,offset+i),plot.index[1]],pcx[c(i,offset+i),plot.index[2]],col=lcol[i])
      }
      points(pcx[1:offset,plot.index[1]][1:gene.no],pcx[1:offset,plot.index[2]][1:gene.no],pch=20,col=tcol,cex=0.3)
      points(pcx[(offset+1):ncol(pcmat),plot.index[1]][1:gene.no],pcx[(offset+1):ncol(pcmat),plot.index[2]][1:gene.no],pch=20,col=tcol2,cex=0.3)
#       #
#       #         ## ## ##
#       #         points(pcx[rownames(pcx) %in% rownames(cm2)[match(galon_macro_mouse,cm2$fromgtf)],plot.index[1]],
#       #                pcx[rownames(pcx) %in% rownames(cm2)[match(galon_macro_mouse,cm2$fromgtf)],plot.index[2]],
#       #                pch=20,col="darkviolet",cex=2)
#       #         points(pcx[rownames(pcx) %in% rownames(cm2)[match(galon_cd8_mouse,cm2$fromgtf)],plot.index[1]],
#       #                pcx[rownames(pcx) %in% rownames(cm2)[match(galon_cd8_mouse,cm2$fromgtf)],plot.index[2]],
#       #                pch=20,col="steelblue",cex=2)
#       #         # legend("topleft",fill = c("darkviolet","steelblue"),legend=c("galon Macro","galon CD8"))
#       mgenes_extended <- c("Tnf","Mmp12","Adam8","Mrc1","Cd36","Cd83","Itgam","Lyz1","Slamf8","Clec12a","Clec10a","Pdc","Cd274")
#       mgenes_extended_ENS <- rownames(cm2)[match(mgenes_extended,cm2$fromgtf)]
#       points(pcx[rownames(pcx) %in% mgenes_extended_ENS,plot.index[1]],
#              pcx[rownames(pcx) %in% mgenes_extended_ENS,plot.index[2]],
#              pch=20,col="salmon",cex=2)
      ## ## ##

    })



    #       output$plot_brushinfo <- renderPrint({
    #         cat("input$pcagenes_brush:\n")
    #         str(input$pcagenes_brush)
    #       })

    curData_brush_ruf <- reactive({
      pcmat <- obj3()[[1]]
      tcol <- obj3()[[2]]
      tcol2 <- obj3()[[3]]

      pres <- prcomp(t(pcmat),scale=FALSE)

      plot.index <- c(as.integer(input$pc_x_ruf),as.integer(input$pc_y_ruf))
      offset <- ncol(pcmat)/2
      gene.no <- offset
      pcx <- pres$x
      # set.seed(11)
      # for (i in 1:ncol(pcx)) {
      #   pcx[,i] <- pcx[,i] + rnorm(nrow(pcx),sd=diff(range(pcx[,i]))/100)
      # }

      firstPCselected <- c(
        pcx[1:offset,plot.index[1]][1:gene.no],
        pcx[(offset+1):ncol(pcmat),plot.index[1]][1:gene.no])

      # pcx[1:offset,plot.index[1]][1:gene.no],pcx[(offset+1):ncol(pcmat),plot.index[1]][1:gene.no],

      secondPCselected <- c(
        pcx[1:offset,plot.index[2]][1:gene.no],
        pcx[(offset+1):ncol(pcmat),plot.index[2]][1:gene.no]
      )

      # pcx[(offset+1):ncol(pcmat),plot.index[2]][1:gene.no]
      pcspcs <- data.frame(firstPC=firstPCselected,secondPC=secondPCselected,geneID=colnames(pcmat))
      rownames(pcspcs) <- c(paste0(colnames(pcmat)[1:gene.no],"_WT"),
                            paste0(colnames(pcmat)[(gene.no+1):(2*gene.no)],"_G37"))

      if(!is.null(values$myannotation))
        pcspcs$geneName <- values$myannotation$gene_name[match(pcspcs$geneID,rownames(values$myannotation))]


      res <- brushedPoints(pcspcs, input$pcaruf_brush,xvar="firstPC",yvar="secondPC",)
      res
    })



    output$pcaruf_out <- DT::renderDataTable({
      datatable(curData_brush_ruf())


    }) # IDEALLY HERE?,options = list(lengthMenu = c(25, 50, 100), pageLength = 100)
    #
    #
    #
    #
    output$downloadData_brush_ruf <- downloadHandler(
      filename = function() { paste(input$brushedPoints_filename_ruf, '.csv', sep='') },
      content = function(file) {
        if(length(input$pcaruf_out_rows_selected)){
          data <- curData_brush_ruf()[input$pcaruf_out_rows_selected,]
        } else {
          data <- curData_brush_ruf()
        }
        write.csv(data, file, quote=FALSE)
      }
    )



































    ## all download handlers
    output$downloadData_brush <- downloadHandler(
      filename = function() { paste(input$brushedPoints_filename, '.csv', sep='') },
      content = function(file) {
        if(length(input$pca_brush_out_rows_selected)){
          data <- curData_brush()[input$pca_brush_out_rows_selected,]
        } else {
          data <- curData_brush()
        }
        write.csv(data, file, quote=FALSE)
      }
    )

    output$downloadData_click <- downloadHandler(
      filename = function() { paste(input$clickedPoints_filename, '.csv', sep='') },
      content = function(file) {
        write.csv(curData_click(), file, quote=FALSE)
      }
    )

    output$download_samplesPca <- downloadHandler(
      filename = function() { input$filename_samplesPca },
      content = function(file) {
        ggsave(file, exportPlots$samplesPca, width = input$export_width, height = input$export_height, units = "cm")
      })

    output$download_samplesScree <- downloadHandler(
      filename = function() { input$filename_samplesScree },
      content = function(file) {
        ggsave(file, exportPlots$samplesScree, width = input$export_width, height = input$export_height, units = "cm")
      })

    output$download_genesPca <- downloadHandler(
      filename = function() { input$filename_genesPca },
      content = function(file) {
        ggsave(file, exportPlots$genesPca, width = input$export_width, height = input$export_height, units = "cm")
      })

    output$download_genesZoom <- downloadHandler(
      filename = function() { input$filename_genesZoom },
      content = function(file) {
        ggsave(file, exportPlots$genesZoom, width = input$export_width, height = input$export_height, units = "cm")
      })

    output$download_genesHeatmap <- downloadHandler(
      filename=function(){
        input$filename_genesHeatmap
      },
      content = function(file){
        pdf(file)
        brushedObject <- curData_brush()

        selectedGenes <- brushedObject$ids
        toplot <- assay(values$myrlt)[selectedGenes,]
        rownames(toplot) <- values$myannotation$gene_name[match(rownames(toplot),rownames(values$myannotation))]
        pheatmap(toplot,cluster_cols = as.logical(input$heatmap_colv))
        dev.off()
      }
    )













  })
  shinyApp(ui = newuiui, server = newserver)

}







