## TODOS
# parameters instead of fixed values
# some other tweaks
# try one run from the scratch
# list whatever is required on top of the pca_SUPALIVE function
### EMBED THIS INTO INTERACTIVE HTML REPORTS? e.g. a la maplots!!!!

library("DESeq2")
library("genefilter")
# load("/Volumes/users$/marinif/linuxHome/032-ruf-macrophages/cm2.RData")
# load("goEnrichs.RData")
library("shiny")
load("/Volumes/users$/marinif/linuxHome/F07-schubert/seb_dds_rld.RData")


# library(ggvis)
# library(shinythemes)


library("DT")
library("pheatmap")
library("d3heatmap")
library(shinydashboard)

annotation <- data.frame(gene_id=rownames(cm2),gene_name=cm2$fromgtf,stringsAsFactors = FALSE,row.names = rownames(cm2))

pcaExplorer <- function(obj,obj2,pca2go=NULL,annotation=NULL){
  # stopifnot( is(obj, 'SummarizedExperiment') )
  ## plenty of tests on the selected object(s)
  if ( !require('shiny') ) {
    stop("pca_SUPALIVE requires 'shiny'. Please install it using
         install.packages('shiny')")
  }

  # here something like, if provided as a countmatrix, create dds and rld

  poss_covars <- names(colData(obj))[] # exclude the size factor, not really required?
  ## TODO: colselection also passed as a palette?
  colSelection <- c("navyblue","steelblue","skyblue","darkred","coral3","darksalmon","green4","greenyellow","orange","gold")


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
        width = 350,
        menuItem("App settings",icon = icon("cogs"),
                 selectInput('pc_x', label = 'x-axis PC: ', choices = 1:8, selected = 1),
                 selectInput('pc_y', label = 'y-axis PC: ', choices = 1:8, selected = 2),
                 selectInput('color_by', label = 'color by: ',
                             choices = c(NULL, poss_covars), selected = NULL,multiple = T),
                 numericInput('pca_nrgenes', label = 'Nr of (most variant) genes:', value = 300,min = 50,max = 20000),




                 selectInput('pc_x_G', label = 'x-axis PC: ', choices = 1:8, selected = 1),
                 selectInput('pc_y_G', label = 'y-axis PC: ', choices = 1:8, selected = 2),
                 selectInput('color_by_G', label = 'color by: ',
                             choices = c(NULL, poss_covars), selected = NULL,multiple = T),
                 numericInput('pca_nrgenes_G', label = 'nr of genes: ', value = 300,min = 50,max = 20000),
                 numericInput('pca_point_alpha_G', label = 'alpha: ', value = 1,min = 0,max = 1,step = 0.01),
                 numericInput('pca_label_size_G', label = 'Labels size: ', value = 2,min = 1,max = 8),
                 numericInput('pca_point_size_G', label = 'Points size: ', value = 2,min = 1,max = 8),
                 numericInput('pca_varname_size_G', label = 'Varname size: ', value = 4,min = 1,max = 8),
                 numericInput('pca_scale_arrow_G', label = 'Scaling factor : ', value = 1,min = 0.01,max = 10),
                 checkboxInput("variable_labels","Display variable labels",value = TRUE),



                 selectInput('pc_x_go', label = 'x-axis PC: ', choices = 1:8, selected = 1),
                 selectInput('pc_y_go', label = 'y-axis PC: ', choices = 1:8, selected = 2)

                 ),
        menuItem("Plot settings", icon = icon("paint-brush"),
                 # bstooltip: Use the widgets below to setup general parameters for exporting produced plots
                 numericInput("export_width",label = "Width of exported figures (cm)",value = 30,min = 2),
                 numericInput("export_height",label = "Height of exported figures (cm)",value = 30,min = 2),


                 numericInput('pca_point_size', label = 'Point size:', value = 3,min = 1,max = 8),
                 checkboxInput("sample_labels","Display sample labels",value = TRUE),


                 # tooltips explanation to have less crowded ui and still good docu
                 shinyBS::bsTooltip(
                   "export_width",
                   paste0("Use the widgets below to setup general parameters for exporting produced plots"),
                   "right", options = list(container = "body"))

                 )
      ),

      dashboardBody(
        tabBox(
          width=12,

          tabPanel(
            "About",
            includeMarkdown(system.file("extdata", "about.md",package = "pcaExplorer"))
            ),

          tabPanel(
            "Instructions",
            includeMarkdown(system.file("extdata", "instructions.md",package = "pcaExplorer"))
            ),

          tabPanel(
            "Data Preview",
            h1("Here will go some head of the count dataset, the samples design/covariates and so")
            ),

          tabPanel(
            "Samples View",
            p(h3('principal component analysis'), "PCA projections of sample abundances onto any pair of components."),
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
            fluidRow(
              column(
                width = 6,
                plotOutput("samples_pca_zoom")
                )
              )
            ),

          tabPanel(
            "Genes View",
            p(h3('principal component analysis'), "PCA projections of sample abundances onto any pair of components."),

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

            fluidRow(
              column(
                width = 6,
                h4("Points selected by brushing - clicking and dragging:"),
                dataTableOutput("pca_brush_out"),
                downloadButton('downloadData_brush', 'Download brushed points'),
                textInput("brushedPoints_filename","File name...")),
              column(
                width = 6,
                h4("Points selected by clicking:"),
                dataTableOutput("pca_click_out"),
                textInput("clickedPoints_filename","File name..."),
                downloadButton('downloadData_click', 'Download clicked (or nearby) points'))
              )
            ),

          tabPanel(
            "Gene finder",
            fluidRow(
              h1("GeneFinder"),
              textInput("genefinder",label = "type in the name of the gene to search"),
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
                dataTableOutput("dt_pcver_pos")),
              column(width = 3)
            ),

            fluidRow(
              column(4,
                     dataTableOutput("dt_pchor_neg")),
              column(4,
                     plotOutput("pca2go")),
              column(4,
                     dataTableOutput("dt_pchor_pos"))
              ),
            fluidRow(
              column(width = 3),
              column(
                width = 6,
                dataTableOutput("dt_pcver_neg")),
              column(width = 3)
            )
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

    user_settings <- reactiveValues(save_width = 45, save_height = 11)


    output$samples_pca <- renderPlot({
      res <- pcaplot(obj,intgroup = input$color_by,ntop = input$pca_nrgenes,
                     pcX = as.integer(input$pc_x),pcY = as.integer(input$pc_y),
                     text_labels = input$sample_labels,
                     point_size = input$pca_point_size, title="Samples PCA"
      )
      res <- res + theme_bw()
      exportPlots$samplesPca <- res
      res
    })

    output$samples_pca_zoom <- renderPlot({
      if(is.null(input$pca_brush))
        return(ggplot() + annotate("text",label="zoom in by brushing",0,0) + theme_bw())

      res <- pcaplot(obj,intgroup = input$color_by,ntop = input$pca_nrgenes,
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
      rv <- rowVars(assay(obj))
      select <- order(rv, decreasing = TRUE)[seq_len(min(input$pca_nrgenes,length(rv)))]
      pca <- prcomp(t(assay(obj)[select, ]))

      res <- pcascree(pca,type = input$scree_type, pc_nr = input$scree_pcnr, title="Scree plot for the samples PCA")
      res <- res + theme_bw()
      exportPlots$samplesScree <- res
      res
    })






    output$genes_biplot <- renderPlot({
      if(!is.null(input$color_by)) {
        expgroups <- colData(obj)[input$color_by][[1]]
      } else {
        expgroups <- colnames(obj)
      }
      colGroups <- colSelection[factor(expgroups)]

      res <- genepca(obj,
                     ntop = input$pca_nrgenes_G,
                     choices = c(as.integer(input$pc_x_G),as.integer(input$pc_y_G)),
                     biplot = TRUE,
                     arrowColors = factor(colGroups),
                     alpha=input$pca_point_alpha_G,coordEqual=F,useRownamesAsLabels=FALSE,labels.size=input$pca_label_size_G,
                     point_size=input$pca_point_size_G,varname.size=input$pca_varname_size_G, scaleArrow = input$pca_scale_arrow_G,annotation=annotation)
      exportPlots$genesPca <- res
      res
    })


    output$genes_biplot_zoom <- renderPlot({
      if(is.null(input$pcagenes_brush)) return(ggplot() + annotate("text",label="zoom in by brushing",0,0) + theme_bw())
      if(!is.null(input$color_by)) {
        expgroups <- colData(obj)[input$color_by][[1]]
      } else {
        expgroups <- colnames(obj)
      }
      colGroups <- colSelection[factor(expgroups)]

      res <- genepca(obj,
                     ntop = input$pca_nrgenes_G,
                     choices = c(as.integer(input$pc_x_G),as.integer(input$pc_y_G)),
                     biplot = TRUE,
                     arrowColors = factor(colGroups),
                     alpha=input$pca_point_alpha_G,coordEqual=F,
                     var.axes=input$variable_labels, # workaround for a ggplot2 bug/missing thing: here details: https://github.com/hadley/ggplot2/issues/905
                     labels.size=input$pca_label_size_G,varname.size=input$pca_varname_size_G,
                     scaleArrow = input$pca_scale_arrow_G,point_size=input$pca_point_size_G,annotation=annotation)

      res <- res +
        xlim(input$pcagenes_brush$xmin,input$pcagenes_brush$xmax) +
        ylim(input$pcagenes_brush$ymin,input$pcagenes_brush$ymax)
      exportPlots$genesZoom <- res
      res
    })


    output$genes_biplot_boxplot <- renderPlot({
      if(length(input$color_by_G)==0) return(ggplot() + annotate("text",label="select an experimental factor",0,0) + theme_bw())
      if(is.null(input$pcagenes_zoom_click)) return(ggplot() + annotate("text",label="click to generate the boxplot\nfor the selected gene",0,0) + theme_bw())

      selectedGene <- curData_zoomClick()$ids
      selectedGeneSymbol <- annotation$gene_name[match(selectedGene,rownames(annotation))]
      # plotCounts(dds_cleaner,)
      genedata <- plotCounts(obj2,gene=selectedGene,intgroup = input$color_by_G,returnData = T)

      onlyfactors <- genedata[,match(input$color_by_G,colnames(genedata))]
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
      df2 <- genepca(obj,
                     ntop = input$pca_nrgenes_G,
                     choices = c(as.integer(input$pc_x_G),as.integer(input$pc_y_G)),
                     biplot = TRUE,
                     # arrowColors = colGroups,
                     alpha=input$pca_point_alpha_G,
                     returnData=T,annotation=annotation)
      df2$geneName <- annotation$gene_name[match(rownames(df2),rownames(annotation))]
      res <- brushedPoints(df2, input$pcagenes_brush,xvar="xvar",yvar="yvar",)
      res
    })


    curData_click <- reactive({
      df2 <- genepca(obj,
                     ntop = input$pca_nrgenes_G,
                     choices = c(as.integer(input$pc_x_G),as.integer(input$pc_y_G)),
                     biplot = TRUE,
                     # arrowColors = colGroups,
                     alpha=input$pca_point_alpha_G,
                     returnData=T,annotation=annotation)
      df2$geneName <- annotation$gene_name[match(rownames(df2),rownames(annotation))]
      res <- nearPoints(df2, input$pcagenes_click,
                        threshold = 20, maxpoints = 3,
                        addDist = TRUE)
      # res <- brushedPoints(df2, input$pcagenes_brush,xvar="xvar",yvar="yvar",)
      res
    })


    curData_zoomClick <- reactive({
      df2 <- genepca(obj,
                     ntop = input$pca_nrgenes_G,
                     choices = c(as.integer(input$pc_x_G),as.integer(input$pc_y_G)),
                     biplot = TRUE,
                     # arrowColors = colGroups,
                     alpha=input$pca_point_alpha_G,
                     returnData=T,annotation=annotation)
      df2$geneName <- annotation$gene_name[match(rownames(df2),rownames(annotation))]
      res <- nearPoints(df2, input$pcagenes_zoom_click,
                        threshold = 20, maxpoints = 1,
                        addDist = TRUE)
      # res <- brushedPoints(df2, input$pcagenes_brush,xvar="xvar",yvar="yvar",)
      res
    })



    output$pca_brush_out <- DT::renderDataTable({
      datatable(curData_brush(),options = list(pageLength = 50))
    })

    output$pca_click_out <- renderDataTable({
      datatable(curData_click(),options = list(pageLength = 50))
    })




    output$heatzoomd3 <- renderD3heatmap({
      if(is.null(input$pcagenes_brush)) return(NULL)

      brushedObject <- curData_brush()
      selectedGenes <- brushedObject$ids
      toplot <- assay(obj)[selectedGenes,]
      rownames(toplot) <- annotation$gene_name[match(rownames(toplot),rownames(annotation))]

      mycolss <- c("#313695","#4575b4","#74add1","#abd9e9","#e0f3f8","#fee090","#fdae61","#f46d43","#d73027","#a50026") # to be consistent with red/blue usual coding

      d3heatmap(toplot,Colv = as.logical(input$heatmap_colv),colors = mycolss)
    })


    output$heatzoom <- renderPlot({
      if(is.null(input$pcagenes_brush)) return(NULL)

      brushedObject <- curData_brush()
      selectedGenes <- brushedObject$ids
      toplot <- assay(obj)[selectedGenes,]
      rownames(toplot) <- annotation$gene_name[match(rownames(toplot),rownames(annotation))]
      # pheatmap(toplot,cluster_cols = as.logical(input$heatmap_colv))
      NMF::aheatmap(toplot,Colv = as.logical(input$heatmap_colv))
      ## aheatmap is actually consistent in displaying the clusters with most of other heatmap packages
      ## keep in mind: pheatmap does somehow a better job if scaling/centering
    })








    output$searchresult <- renderPrint({

      if(is.null(input$color_by_G)) return("Select a factor to plot your gene")
      if(is.null(input$genefinder))
        return("Type in the gene name/id you want to plot")

      foundGeneID <- input$genefinder %in% rownames(obj)
      foundGeneName <- input$genefinder %in% annotation$gene_name
      if(!foundGeneID){
        foundGeneID <- toupper(input$genefinder) %in% toupper(rownames(obj))
        if(foundGeneID){
          return(paste0("Maybe you mis-spelled the name of your gene. Did you mean ",
                        unique(rownames(annotation)[which(toupper(input$genefinder)==toupper(rownames(annotation)))]),"?"))
        } else {
          foundGeneNAME <- input$genefinder %in% annotation$gene_name
          if(!foundGeneNAME){
            foundGeneNAME <- toupper(input$genefinder) %in% toupper(annotation$gene_name)
            if(foundGeneNAME){
              return(paste0("Maybe you mis-spelled the name of your gene. Did you mean ",
                            unique(annotation$gene_name[which(toupper(input$genefinder)==toupper(annotation$gene_name))]),"?"))
            } else {return("Could not find the gene you typed!")}
          } else {
            fgn <- annotation$gene_name[which(annotation$gene_name==input$genefinder)]
            if (length(fgn) > 1) return(paste0("Found more than one gene with the selected gene name. Select one of the following: ",paste(selectedGene,collapse=", ")))
            selectedGene <- rownames(annotation)[which(annotation$gene_name==input$genefinder)]

            fg <- rownames(annotation)[match(fgn,annotation$gene_name)]
            return(paste0("I found the gene! Plotting ", fg, " - ", annotation$gene_name[match(fg,rownames(annotation))],"..."))

          }}
      } else {
        fg <- rownames(annotation)[match(input$genefinder,rownames(obj))]
        return(paste0("I found the gene! Plotting ", fg, " - ", annotation$gene_name[match(fg,rownames(annotation))],"..."))

      }
    })





    output$genefinder_plot <- renderPlot({
      anno_id <- rownames(annotation)
      anno_gene <- annotation$gene_name

      if(is.null(input$color_by_G))
        return(ggplot() + annotate("text",label="Select a factor to plot your gene",0,0) + theme_bw())
      if(is.null(input$genefinder))
        return(ggplot() + annotate("text",label="Type in a gene name/id",0,0) + theme_bw())
      if(!input$genefinder %in% anno_gene & !input$genefinder %in% anno_id)
        return(ggplot() + annotate("text",label="gene not found...",0,0) + theme_bw())

      if (input$genefinder %in% anno_id) {
        selectedGene <- rownames(obj)[match(input$genefinder,rownames(obj))]
        selectedGeneSymbol <- annotation$gene_name[match(selectedGene,rownames(annotation))]
      }
      if (input$genefinder %in% anno_gene) {
        selectedGeneSymbol <- annotation$gene_name[which(annotation$gene_name==input$genefinder)]
        if (length(selectedGeneSymbol) > 1) return(ggplot() + annotate("text",label=paste0("Type in a gene name/id of the following:\n",paste(selectedGene,collapse=", ")),0,0) + theme_bw())
        selectedGene <- rownames(annotation)[which(annotation$gene_name==input$genefinder)]
      }
      genedata <- plotCounts(obj2,gene=selectedGene,intgroup = input$color_by_G,returnData = T)

      onlyfactors <- genedata[,match(input$color_by_G,colnames(genedata))]
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
      if(is.null(pca2go))
        return(ggplot() + annotate("text",label="Provide a pca2go object to the app",0,0) + theme_bw())
      res <- pcaplot(obj,intgroup = input$color_by,
                     ntop = input$pca_nrgenes,pcX = as.integer(input$pc_x_go),pcY = as.integer(input$pc_y_go),text_labels = input$sample_labels,
                     point_size = input$pca_point_size, title="PCA on the samples"

      )
      res
    })


    output$dt_pchor_pos <- renderDataTable({
      if(is.null(pca2go)) return(datatable(NULL))
      goe <- pca2go[[paste0("PC",input$pc_x_go)]][["posLoad"]]
      if(input$compact_pca2go)
        return(datatable(goe[,c("GO.ID","Term","Significant","p.value_elim")]))
      datatable(goe)
    })

    output$dt_pchor_neg <- renderDataTable({
      if(is.null(pca2go)) return(datatable(NULL))
      goe <- pca2go[[paste0("PC",input$pc_x_go)]][["negLoad"]]
      if(input$compact_pca2go)
        return(datatable(goe[,c("GO.ID","Term","Significant","p.value_elim")]))
      datatable(goe)
    })

    output$dt_pcver_pos <- renderDataTable({
      if(is.null(pca2go)) return(datatable(NULL))
      goe <- pca2go[[paste0("PC",input$pc_y_go)]][["posLoad"]]
      if(input$compact_pca2go)
        return(datatable(goe[,c("GO.ID","Term","Significant","p.value_elim")]))
      datatable(goe)
    })

    output$dt_pcver_neg <- renderDataTable({
      if(is.null(pca2go)) return(datatable(NULL))
      goe <- pca2go[[paste0("PC",input$pc_y_go)]][["negLoad"]]
      if(input$compact_pca2go)
        return(datatable(goe[,c("GO.ID","Term","Significant","p.value_elim")]))
      datatable(goe)
    })

    output$enrichinfo <- renderPrint({
      cat("enrich info:\n")
      # str(goEnrichs)
      class(input$pc_x)
      head(pca2go[[paste0("PC",input$pc_x)]][["posLoad"]])
      class(datatable(pca2go[[paste0("PC",input$pc_x)]][["posLoad"]]))
    })








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
        toplot <- assay(obj)[selectedGenes,]
        rownames(toplot) <- annotation$gene_name[match(rownames(toplot),rownames(annotation))]
        pheatmap(toplot,cluster_cols = as.logical(input$heatmap_colv))
        dev.off()
      }
    )













  })
  shinyApp(ui = newuiui, server = newserver)

}







