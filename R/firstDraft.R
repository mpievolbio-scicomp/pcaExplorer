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


library(ggvis)
library(shinythemes)


library("DT")
library("pheatmap")
library("d3heatmap")
library(shinydashboard)

pca_SUPALIVE4 <- function(obj,obj2){
  # stopifnot( is(obj, 'SummarizedExperiment') )
  if ( !require('shiny') ) {
    stop("pca_SUPALIVE requires 'shiny'. Please install it using
         install.packages('shiny')")
  }

#   poss_covars <- dplyr::setdiff(
  #     colnames(obj$sample_to_covariates),
  #     'sample')
  #   samp_names <- obj$sample_to_covariates[['sample']]
  #   poss_models <- names(models(obj, verbose = FALSE))

  poss_covars <- names(colData(obj))[] # exclude the size factor, not really required
  # poss_models <- "full"
  # samp_names <- rownames(obj$rotation) # %>% rownames()

  # for the shinyUI
  # var_range <- function(id, label, variable) {
  # rng <- range(variable, na.rm = TRUE)
  # sliderInput(id, label, rng[1], rng[2], rng)
  # }

  newuiui <-
    shinydashboard::dashboardPage(
      dashboardHeader(
        title = paste0("pcaExplorer - Interactive exploration of Principal Components",
                       "of Samples and Genes in RNA-seq data - version ",
                       packageVersion("pcaExplorer")),
        titleWidth = 900),

      dashboardSidebar(
        width = 350,
        menuItem("Settings",icon = icon("cogs"),
                 selectInput('pc_x', label = 'x-axis PC: ', choices = 1:8,
                             selected = 1)
                 ),
        menuItem("Plot settings", icon = icon("paint-brush"),
                 numericInput("export_width",label = "Width of exported figures (cm)",value = 30,min = 2),
                 numericInput("export_height",label = "Height of exported figures (cm)",value = 30,min = 2)
                 )
                       ),

      dashboardBody(

      ),

      skin="blue"








    )

  newserver <- function(input, output) { }
  shinyApp(ui = newuiui, server = newserver)



  uiui <- (
    # var_range <- function(id, label, variable) {
    #   rng <- range(variable, na.rm = TRUE)
    #   sliderInput(id, label, rng[1], rng[2], rng)
    # }

    navbarPage(theme = shinytheme("journal"),"PCALIVE - beta version",

               tabPanel("Overview",
                        fluidRow(
                          div(h3('pcaLive'), align = 'center')
                        ),
                        fluidRow("A Shiny application for exploring expression data in different conditions and experimental factors.\nUse the widgets below to setup general parameters for exporting produced plots"),
                        fluidRow(
                          numericInput("export_width",label = "Width of exported figures (cm)",value = 30,min = 2),
                          numericInput("export_height",label = "Height of exported figures (cm)",value = 30,min = 2)
                        )
               ),
               navbarMenu("Analysis on...",
                          tabPanel('Samples',
                                   fluidRow(
                                     column(12,
                                            p(h3('principal component analysis'), "PCA projections of sample abundances onto any pair of components.")
                                     ),
                                     offset = 1),

                                   fluidRow(
                                     column(2,
                                            selectInput('pc_x', label = 'x-axis PC: ', choices = 1:8,
                                                        selected = 1)
                                     ),
                                     column(2,
                                            selectInput('pc_y', label = 'y-axis PC: ', choices = 1:8,
                                                        selected = 2)
                                     ),
                                     column(3,
                                            selectInput('color_by', label = 'color by: ',
                                                        choices = c(NULL, poss_covars), selected = NULL,multiple = T)
                                     ),
                                     #                                      column(2,
                                     #                                             numericInput('pca_point_alpha', label = 'Point alpha:', value = 1,min = 0,max = 1,step = 0.01)),
                                     column(2,
                                            numericInput('pca_point_size', label = 'Point size:', value = 3,min = 1,max = 8)),

                                     column(2,
                                            numericInput('pca_nrgenes', label = 'Nr of (most variant) genes:', value = 300,min = 50,max = 20000))
                                   ),
                                   fluidRow(
                                     column(3,
                                            checkboxInput("sample_labels","Display sample labels",value = TRUE))),
                                   fluidRow(
                                     column(6,
                                            plotOutput('pca_plt',brush = "pca_brush"),
                                            div(align = "right", style = "margin-right:15px; margin-bottom:10px",
                                                downloadButton("download_samplesPca", "Download Plot"),
                                                textInput("filename_samplesPca",label = "Save as...",value = "samplesPca.pdf"))),
                                     column(6,
                                            plotOutput("scree"),
                                            div(align = "right", style = "margin-right:15px; margin-bottom:10px",
                                                downloadButton("download_samplesScree", "Download Plot"),
                                                textInput("filename_samplesScree",label = "Save as...",value = "samplesScree.pdf")),
                                            fluidRow(
                                              column(6,
                                                     radioButtons("scree_type","Scree plot type:",choices=list("Proportion of explained variance"="pev","Cumulative proportion of explained variance"="cev"),"pev")),
                                              column(6,
                                                     numericInput("scree_pcnr","Number of PCs to display",value=8,min=2))
                                            )
                                     )
                                   ),
                                   fluidRow(
                                     column(6,
                                            plotOutput("pca_pltZoom"
                                            )
                                     ))
                                   # )obj
                          ),
                          tabPanel("Genes",
                                   fluidRow(
                                     column(12,
                                            p(h3('principal component analysis'), "PCA projections of gene abundances onto any pair of components.")
                                     ),
                                     offset = 1),
                                   fluidRow(
                                     column(2,
                                            selectInput('pc_x_G', label = 'x-axis PC: ', choices = 1:8,
                                                        selected = 1)
                                     ),
                                     column(2,
                                            selectInput('pc_y_G', label = 'y-axis PC: ', choices = 1:8,
                                                        selected = 2)
                                     ),
                                     column(1,
                                            selectInput('color_by_G', label = 'color by: ',
                                                        choices = c(NULL, poss_covars), selected = NULL,multiple = T)
                                     ),
                                     column(1,
                                            numericInput('pca_nrgenes_G', label = 'nr of genes: ', value = 300,min = 50,max = 20000)),
                                     column(1,
                                            numericInput('pca_point_alpha_G', label = 'alpha: ', value = 1,min = 0,max = 1,step = 0.01)),
                                     column(1,
                                            numericInput('pca_label_size_G', label = 'Labels size: ', value = 2,min = 1,max = 8)),
                                     column(1,
                                            numericInput('pca_point_size_G', label = 'Points size: ', value = 2,min = 1,max = 8)),
                                     column(1,
                                            numericInput('pca_varname_size_G', label = 'Varname size: ', value = 4,min = 1,max = 8)),
                                     column(1,
                                            numericInput('pca_scale_arrow_G', label = 'Scaling factor : ', value = 1,min = 0.01,max = 10)),
                                     column(1,
                                            checkboxInput("variable_labels","Display variable labels",value = TRUE))


                                   ),
                                   fluidRow(verbatimTextOutput('x4')),
                                   fluidRow(verbatimTextOutput('x5')),
                                   fluidRow(verbatimTextOutput('x6')),fluidRow(verbatimTextOutput('x8')),
                                   fluidRow(column(4,
                                                   h4("Main Plot - interact!"),
                                                   plotOutput('pca_plt_G',brush = 'pcagenes_brush',click="pcagenes_click"),
                                                   div(align = "right", style = "margin-right:15px; margin-bottom:10px",
                                                       downloadButton("download_genesPca", "Download Plot"),
                                                       textInput("filename_genesPca",label = "Save as...",value = "genesPca.pdf"))),
                                            #                           fluidRow(textOutput('debug'))),
                                            column(4,
                                                   h4("Zoomed window"),
                                                   plotOutput("testzoom",click="pcagenes_zoom_click"),
                                                   div(align = "right", style = "margin-right:15px; margin-bottom:10px",
                                                       downloadButton("download_genesZoom", "Download Plot"),
                                                       textInput("filename_genesZoom",label = "Save as...",value = "genesPca_zoomed.pdf"))),
                                            column(4,
                                                   h4("Boxplot of selected gene"),
                                                   plotOutput("genePlot"))

                                   ),
                                   fluidRow(h4("genefinder"),
                                            textInput("searchgene",label = "type in the name of the gene to search"),
                                            verbatimTextOutput("searchresult"),
                                            verbatimTextOutput("debuggene"),
                                            checkboxInput("ylimZero","Set y axis limit to 0",value=TRUE),
                                            plotOutput("searchgenePlot")),
                                   fluidRow(column(6,
                                                   h4("Zoomed heatmap"),
                                                   plotOutput("heatzoom"),
                                                   div(align = "right", style = "margin-right:15px; margin-bottom:10px",
                                                       downloadButton("download_genesHeatmap","Download Plot"),
                                                       textInput("filename_genesHeatmap",label = "Save as...",value = "genesHeatmap.pdf"))
                                   ),
                                   column(6,
                                          h4("Zoomed interactive heatmap"),
                                          fluidRow(radioButtons("heatmap_colv","Cluster samples",choices = list("Yes"=TRUE,"No"=FALSE),selected = TRUE)),
                                          fluidRow(d3heatmapOutput("heatzoomd3")))),
                                   # verbatimTextOutput("plot_brushinfo"),
                                   # dataTableOutput("restab"),
                                   fluidRow(column(6,
                                                   h4("Points selected by brushing - clicking and dragging:"),
                                                   downloadButton('downloadData_brush', 'Download brushed points'),
                                                   textInput("brushedPoints_filename","File name..."),
                                                   dataTableOutput("pca_brush_out")),
                                            column(6,
                                                   h4("Points selected by clicking:"),
                                                   downloadButton('downloadData_click', 'Download clicked (or nearby) points'),
                                                   textInput("clickedPoints_filename","File name..."),
                                                   dataTableOutput("pca_click_out"))
                                   )
                                   # fluidRow(plotOutput("testzoom"))
                          )

               ),
               tabPanel("interactive genes",headerPanel("Zooming demo"),
                        sidebarPanel(
                          # sliderInput("ngenes_sh","Number of genes to include", 30, nrow(assay(obj)), 500)
                          sliderInput("ngenes_sh","Number of genes to include", 30, 20000, 500)
                          #                  var_range("y_domain", "Y", mtcars$mpg)
                        ),
                        mainPanel(
                          ggvisOutput("p")
                        )),

               tabPanel("pca2go",
                        h2("Functions enriched in the genes with high loadings on the selected principal components"),
                        fluidRow(
                          column(2,
                                 selectInput('pc_x_go', label = 'x-axis PC: ', choices = 1:8,
                                             selected = 1)
                          ),
                          column(2,
                                 selectInput('pc_y_go', label = 'y-axis PC: ', choices = 1:8,
                                             selected = 2)
                          )),
                        verbatimTextOutput("enrichinfo"),
                        fluidRow(dataTableOutput("dt_pcver_pos")),
                        fluidRow(
                          column(4,
                                 dataTableOutput("dt_pchor_neg")),
                          column(4,
                                 plotOutput("pca2go")),
                          column(4,
                                 dataTableOutput("dt_pchor_pos")
                          )),
                        fluidRow(dataTableOutput("dt_pcver_neg"))),

               tabPanel("Multifactor exploration",


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
                                 dataTableOutput('pcaruf_out'))


               )

    )
  )





  serser <-

    shinyServer(function(input, output, session) {

      exportPlots <- reactiveValues(
        samplesPca=NULL,
        samplesZoom=NULL,
        samplesScree=NULL,
        genesPca=NULL,
        genesZoom=NULL,
        genesHeatmap=NULL
      )

      user_settings <- reactiveValues(save_width = 45, save_height = 11)

      reactive({
        rv <- rowVars(assay(obj))
        select <- order(rv, decreasing = TRUE)[seq_len(min(input$ngenes_sh,length(rv)))]
        pca <- prcomp((assay(obj)[select, ]))
        percentVar <- pca$sdev^2/sum(pca$sdev^2)
        nobs.factor <- sqrt(nrow(pca$x) - 1)
        devs <- pca$sdev
        pcast <- pca
        pcast$x <- sweep(pca$x, 2, 1/(devs * nobs.factor), FUN = "*") * nobs.factor
        d <- data.frame(PC1 = pcast$x[, 1], PC2 = pcast$x[, 2],
                        # group = group,
                        # intgroup.df,
                        names = rownames((assay(obj)[select, ])))
        # with the d object...
        d$name <- cm2$fromgtf[match(rownames(d),rownames(cm2))]
        library(ggvis)
        labs <- function(data){
          if(is.null(data)) return(NULL)
          data$name
        }
        # d %>% ggvis( ~PC1, ~ PC2,key:=~ name) %>% layer_points(size.hover := 200) %>% add_tooltip(labs, "hover")

        d %>% ggvis( ~PC1, ~ PC2,key:=~ name) %>% layer_points(size:=20,size.hover := 100,fill.hover:="steelblue") %>% add_tooltip(labs, "hover") }) %>%
        bind_shiny("p")





      output$debug <- reactive({      cat(paste(input$pca_x_G,input$pca_y_G))    })

      output$pca_plt <- renderPlot({
        res <- pcaplot(obj,intgroup = input$color_by,
                      ntop = input$pca_nrgenes,pcX = as.integer(input$pc_x),pcY = as.integer(input$pc_y),text_labels = input$sample_labels,
                      point_size = input$pca_point_size, title="PCA on the samples"

        )
        res <- res + theme_bw()
        exportPlots$samplesPca <- res
        res
      })

      output$pca_pltZoom <- renderPlot({
        if(is.null(input$pca_brush)) return(ggplot() + annotate("text",label="zoom in by brushing",0,0) + theme_bw())
        res <- pcaplot(obj,intgroup = input$color_by,
                      ntop = input$pca_nrgenes,pcX = as.integer(input$pc_x),pcY = as.integer(input$pc_y),text_labels = input$sample_labels,
                      point_size = input$pca_point_size, title="PCA on the samples"

        )
        res <- res + xlim(input$pca_brush$xmin,input$pca_brush$xmax) + ylim(input$pca_brush$ymin,input$pca_brush$ymax)
        res <- res + theme_bw()
        exportPlots$samplesZoom <- res
        res
      })


      output$scree <- renderPlot({
        rv <- rowVars(assay(obj))
        select <- order(rv, decreasing = TRUE)[seq_len(min(input$pca_nrgenes,length(rv)))]
        pca <- prcomp(t(assay(obj)[select, ]))

        res <- pcascree(pca,type = input$scree_type, pc_nr = input$scree_pcnr, title="Scree plot for the samples PCA")
        res <- res + theme_bw()
        exportPlots$samplesScree <- res
        res
      })

      output$pca_plt_G <- renderPlot({

        colSelection <- c("navyblue","steelblue","skyblue","darkred","coral3","darksalmon","green4","greenyellow","orange","gold")
        color_by <- ifelse(is.null(input$color_by), NULL,
                           as.character(input$color_by))
        expgroups <- gsub(pattern = "_R.*","",colnames(obj))
        colGroups <- colSelection[factor(expgroups,levels = unique(as.character(expgroups)))]


        # genepca(rld_global,ngenes,choices = c(1,2), biplot = TRUE,arrowColors = colGroups,alpha=0.3)
        res <- genepca(obj,
                       ntop = input$pca_nrgenes_G,
                       choices = c(as.integer(input$pc_x_G),as.integer(input$pc_y_G)),
                       biplot = TRUE,
                       arrowColors = factor(colGroups),
                       alpha=input$pca_point_alpha_G,coordEqual=F,useRownamesAsLabels=FALSE,labels.size=input$pca_label_size_G,
                       point_size=input$pca_point_size_G,varname.size=input$pca_varname_size_G, scaleArrow = input$pca_scale_arrow_G)
        exportPlots$genesPca <- res
        res

      })

      output$plot_brushinfo <- renderPrint({
        cat("input$pcagenes_brush:\n")
        str(input$pcagenes_brush)
      })

      curData_brush <- reactive({
        #         colSelection <- c("navyblue","steelblue","skyblue","darkred","coral3","darksalmon","green4","greenyellow","orange","gold")
        #         color_by <- ifelse(is.null(input$color_by), NULL,
        #                            as.character(input$color_by))
        #         expgroups <- gsub(pattern = "_R.*","",colnames(obj))
        #         colGroups <- colSelection[factor(expgroups,levels = unique(as.character(expgroups)))]
        #
        df2 <- genepca(obj,
                       ntop = input$pca_nrgenes_G,
                       choices = c(as.integer(input$pc_x_G),as.integer(input$pc_y_G)),
                       biplot = TRUE,
                       # arrowColors = colGroups,
                       alpha=input$pca_point_alpha_G,
                       returnData=T)
        df2$geneName <- cm2$fromgtf[match(rownames(df2),rownames(cm2))]


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
                       returnData=T)
        df2$geneName <- cm2$fromgtf[match(rownames(df2),rownames(cm2))]
        res <- nearPoints(df2, input$pcagenes_click,
                          threshold = 40, maxpoints = 10,
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
                       returnData=T)
        df2$geneName <- cm2$fromgtf[match(rownames(df2),rownames(cm2))]
        res <- nearPoints(df2, input$pcagenes_zoom_click,
                          threshold = 20, maxpoints = 1,
                          addDist = TRUE)
        # res <- brushedPoints(df2, input$pcagenes_brush,xvar="xvar",yvar="yvar",)
        res
      })



      output$pca_brush_out <- DT::renderDataTable({
        datatable(curData_brush())


      }) # IDEALLY HERE?,options = list(lengthMenu = c(25, 50, 100), pageLength = 100)

      output$pca_click_out <- renderDataTable({
        datatable(curData_click())
      })


      output$testzoom <- renderPlot({
        if(is.null(input$pcagenes_brush)) return(ggplot() + annotate("text",label="zoom in by brushing",0,0) + theme_bw())
        colSelection <- c("navyblue","steelblue","skyblue","darkred","coral3","darksalmon","green4","greenyellow","orange","gold")
        color_by <- ifelse(is.null(input$color_by), NULL,
                           as.character(input$color_by))
        expgroups <- gsub(pattern = "_R.*","",colnames(obj))
        colGroups <- colSelection[factor(expgroups,levels = unique(as.character(expgroups)))]


        # genepca(rld_global,ngenes,choices = c(1,2), biplot = TRUE,arrowColors = colGroups,alpha=0.3)
        res <- genepca(obj,
                       ntop = input$pca_nrgenes_G,
                       choices = c(as.integer(input$pc_x_G),as.integer(input$pc_y_G)),
                       biplot = TRUE,
                       arrowColors = factor(colGroups),
                       alpha=input$pca_point_alpha_G,coordEqual=F,
                       var.axes=input$variable_labels, # workaround for a ggplot2 bug/missing thing: here details: https://github.com/hadley/ggplot2/issues/905
                       labels.size=input$pca_label_size_G,varname.size=input$pca_varname_size_G, scaleArrow = input$pca_scale_arrow_G,point_size=input$pca_point_size_G) + xlim(input$pcagenes_brush$xmin,input$pcagenes_brush$xmax) + ylim(input$pcagenes_brush$ymin,input$pcagenes_brush$ymax)
        exportPlots$genesZoom <- res
        res
      })


      output$genePlot <- renderPlot({
        if(length(input$color_by_G)==0) return(ggplot() + annotate("text",label="select an experimental factor",0,0) + theme_bw())
        if(is.null(input$pcagenes_zoom_click)) return(ggplot() + annotate("text",label="click to generate the boxplot\nfor the selected gene",0,0) + theme_bw())

        selectedGene <- curData_zoomClick()$ids
        selectedGeneSymbol <- cm2$fromgtf[match(selectedGene,rownames(cm2))]
        # plotCounts(dds_cleaner,)
        genedata <- plotCounts(obj2,gene=selectedGene,intgroup = input$color_by_G,returnData = T)

        # onlyfactors <- genedata[match(input$color_by_G,colnames(genedata))]
        onlyfactors <- genedata[,match(input$color_by_G,colnames(genedata))]

        ## intgroup can be a vector of factors. then i need interactions of the two factors
        # plotCounts(ddsmf_global,gene="ENSMUSG00000026080",intgroup = c("tissue","condition"),returnData = T) -> dh
        # dh$f1f2 <- interaction(dh$tissue,dh$condition)
        # dh  %>% ggplot(aes(x=f1f2,y=count,fill=f1f2)) + geom_boxplot()


        ## TODO: make the intgroup/colr by also somehow multiple-selectable?

        # genedata$plotby <- lapply(1:ncol(onlyfactors),function(arg) onlyfactors[,arg]) %>% interaction()
        genedata$plotby <- interaction(onlyfactors)



        ggplot(genedata,aes(x=plotby,y=count,fill=plotby)) + geom_boxplot(outlier.shape = NA) + scale_y_log10(name="Normalized counts") + labs(title=paste0("Normalized counts for ",selectedGeneSymbol," - ",selectedGene)) +  scale_x_discrete(name="") + geom_jitter(aes(x=plotby,y=count),position = position_jitter(width = 0.1)) + scale_fill_discrete(name="Experimental\nconditions")
        # exportPlots$genesZoom <- res
        # res
      })


      output$heatzoomd3 <- renderD3heatmap({
        if(is.null(input$pcagenes_brush)) return(NULL)

        brushedObject <- curData_brush()

        selectedGenes <- brushedObject$ids
        toplot <- assay(obj)[selectedGenes,]
        rownames(toplot) <- cm2$fromgtf[match(rownames(toplot),rownames(cm2))]

        mycolss <- c("#313695","#4575b4","#74add1","#abd9e9","#e0f3f8","#fee090","#fdae61","#f46d43","#d73027","#a50026") # to be consistent with red/blue usual coding

        d3heatmap(toplot,Colv = as.logical(input$heatmap_colv),colors = mycolss)


      })

      output$heatzoom <- renderPlot({
        if(is.null(input$pcagenes_brush)) return(NULL)

        brushedObject <- curData_brush()

        selectedGenes <- brushedObject$ids
        toplot <- assay(obj)[selectedGenes,]
        rownames(toplot) <- cm2$fromgtf[match(rownames(toplot),rownames(cm2))]
        # pheatmap(toplot,cluster_cols = as.logical(input$heatmap_colv))
        NMF::aheatmap(toplot,Colv = as.logical(input$heatmap_colv))

        ## aheatmap is actually consistent in displaying the clusters with most of other heatmap packages
        ## keep in mind: pheatmap does somehow a better job if scaling/centering



      })



      output$plot_brushinfo <- renderPrint({
        cat("input$pcagenes_brush:\n")
        str(input$pcagenes_brush)
      })

      output$x4 = renderPrint({
        s = input$pca_brush_out_rows_selected
        if (length(s)) {
          cat('These rows were selected:\n\n')
          cat(s, sep = ', ')
        }
      })
      output$x5 = renderPrint({
        s = input$pca_brush_out_rows_selected
        curData_brush()[s,]
      })
      output$x6 = renderPrint({
        length(input$pca_brush_out_rows_selected)
      })

      output$x7 <- renderPrint({
        str(input$pca_brush)
      })

      output$x8 <- renderPrint({
        str(input$color_by_G)
      })



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
          rownames(toplot) <- cm2$fromgtf[match(rownames(toplot),rownames(cm2))]
          pheatmap(toplot,cluster_cols = as.logical(input$heatmap_colv))
          dev.off()
        }
      )




      output$pca2go <- renderPlot({
        res <- pcaplot(obj,intgroup = input$color_by,
                      ntop = input$pca_nrgenes,pcX = as.integer(input$pc_x_go),pcY = as.integer(input$pc_y_go),text_labels = input$sample_labels,
                      point_size = input$pca_point_size, title="PCA on the samples"

        )
        res
      })


      output$dt_pchor_pos <- renderDataTable({
        datatable(goEnrichs[[paste0("PC",input$pc_x_go)]][["posLoad"]])
      })

      output$dt_pchor_neg <- renderDataTable({
        datatable(goEnrichs[[paste0("PC",input$pc_x_go)]][["negLoad"]])

      })

      output$dt_pcver_pos <- renderDataTable({
        datatable(goEnrichs[[paste0("PC",input$pc_y_go)]][["posLoad"]])

      })

      output$dt_pcver_neg <- renderDataTable({
        datatable(goEnrichs[[paste0("PC",input$pc_y_go)]][["negLoad"]])

      })

      output$enrichinfo <- renderPrint({
        cat("enrich info:\n")
        # str(goEnrichs)
        class(input$pc_x)
        head(goEnrichs[[paste0("PC",input$pc_x)]][["posLoad"]])
        class(datatable(goEnrichs[[paste0("PC",input$pc_x)]][["posLoad"]]))
      })






      ## from here on, RUF APP


      obj3 <- reactive({
        # preliminary on the object to morph into obj3
        exprmat <- t(assay(obj))


        # removing non expressed genes in advance?
        # expressed <- counts(ddsmf_global)#...
        (rowSums(counts(obj2) > 5)>2) %>% sum
        (rowSums(counts(obj2) > 3)>1) %>% sum
        (rowSums(counts(obj2) > 3)>2) %>% sum
        (rowSums(counts(obj2) > 5)>1) %>% sum
        (rowSums(counts(obj2) > 5)>3) %>% sum
        sum(rowSums(counts(obj2))>0)

        exprmat <- exprmat[,rowSums(counts(obj2) > 5)>2]


        ## original
        pcmat <- cbind(exprmat[c("WT_macro_R1","WT_macro_R2","WT_macro_R3","WT_macro_R4",
                                 "WT_endo_R1","WT_endo_R2","WT_endo_R3","WT_endo_R4",
                                 "WT_CD11b_R1","WT_CD11b_R2","WT_CD11b_R3","WT_CD11b_R4",
                                 "WT_CD8_R1","WT_CD8_R2","WT_CD8_R3","WT_CD8_R4"),],
                       exprmat[c("G37I_macro_R1","G37I_macro_R2","G37I_macro_R3","G37I_macro_R4",
                                 "G37I_endo_R1","G37I_endo_R2","G37I_endo_R3","G37I_endo_R4",
                                 "G37I_CD11b_R1","G37I_CD11b_R2","G37I_CD11b_R5","G37I_CD11b_R4",
                                 "G37I_CD8_R1","G37I_CD8_R3","G37I_CD8_R3","G37I_CD8_R4"),])
        ## to check whether it is robust
        #         pcmat <- cbind(exprmat[c("WT_macro_R1","WT_macro_R2","WT_macro_R3","WT_macro_R4",
        #                                  "WT_endo_R1","WT_endo_R2","WT_endo_R3","WT_endo_R4",
        #                                  "WT_CD11b_R1","WT_CD11b_R2","WT_CD11b_R5","WT_CD11b_R4",
        #                                  "WT_CD8_R1","WT_CD8_R2","WT_CD8_R3","WT_CD8_R4"),],
        #                        exprmat[c("G37I_macro_R3","G37I_macro_R2","G37I_macro_R4","G37I_macro_R1",
        #                                  "G37I_endo_R4","G37I_endo_R3","G37I_endo_R1","G37I_endo_R2",
        #                                  "G37I_CD11b_R2","G37I_CD11b_R3","G37I_CD11b_R1","G37I_CD11b_R4",
        #                                  "G37I_CD8_R1","G37I_CD8_R4","G37I_CD8_R2","G37I_CD8_R3"),])


        library(scales)
        aval <- 0.3


        max.type <- apply(pcmat[,1:ncol(exprmat)],2,which.max)
        tcol.justMax <- ifelse(max.type <= 4,alpha("green",aval),ifelse(max.type <= 8,alpha("red",aval),ifelse(max.type <= 12,alpha("blue",aval),alpha("orange",aval))))
        max.type2 <- apply(pcmat[,(ncol(exprmat)+1):ncol(pcmat)],2,which.max)
        tcol2.justMax <- ifelse(max.type2 <= 4,alpha("green",aval),ifelse(max.type2 <= 8,alpha("red",aval),ifelse(max.type2 <= 12,alpha("blue",aval),alpha("orange",aval))))

        # using the median across replicates
        celltypes <- gsub("_R.","",rownames(pcmat))

        # attach this to the matrix to split accordingly
        # pcmat2 <- data.frame(pcmat,celltypes=celltypes,stringsAsFactors = F)


        # minipcmat <- pcmat[1:16,1:10]
        # pcmat2 <- data.frame(minipcmat,celltypes=celltypes,stringsAsFactors = F)
        # pcmat3 <- gather(as.data.frame(pcmat2),id,exprvalue,-celltypes)

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


        ## ## ##
        #         points(pcx[rownames(pcx) %in% rownames(cm2)[match(galon_macro_mouse,cm2$fromgtf)],plot.index[1]],
        #                pcx[rownames(pcx) %in% rownames(cm2)[match(galon_macro_mouse,cm2$fromgtf)],plot.index[2]],
        #                pch=20,col="darkviolet",cex=2)
        #         points(pcx[rownames(pcx) %in% rownames(cm2)[match(galon_cd8_mouse,cm2$fromgtf)],plot.index[1]],
        #                pcx[rownames(pcx) %in% rownames(cm2)[match(galon_cd8_mouse,cm2$fromgtf)],plot.index[2]],
        #                pch=20,col="steelblue",cex=2)
        #         # legend("topleft",fill = c("darkviolet","steelblue"),legend=c("galon Macro","galon CD8"))
        mgenes_extended <- c("Tnf","Mmp12","Adam8","Mrc1","Cd36","Cd83","Itgam","Lyz1","Slamf8","Clec12a","Clec10a","Pdc","Cd274")
        mgenes_extended_ENS <- rownames(cm2)[match(mgenes_extended,cm2$fromgtf)]
        points(pcx[rownames(pcx) %in% mgenes_extended_ENS,plot.index[1]],
               pcx[rownames(pcx) %in% mgenes_extended_ENS,plot.index[2]],
               pch=20,col="salmon",cex=2)
        ## ## ##
      })

      output$rufdebug <- renderPrint({
        str(obj3())
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
        #
        #         ## ## ##
        #         points(pcx[rownames(pcx) %in% rownames(cm2)[match(galon_macro_mouse,cm2$fromgtf)],plot.index[1]],
        #                pcx[rownames(pcx) %in% rownames(cm2)[match(galon_macro_mouse,cm2$fromgtf)],plot.index[2]],
        #                pch=20,col="darkviolet",cex=2)
        #         points(pcx[rownames(pcx) %in% rownames(cm2)[match(galon_cd8_mouse,cm2$fromgtf)],plot.index[1]],
        #                pcx[rownames(pcx) %in% rownames(cm2)[match(galon_cd8_mouse,cm2$fromgtf)],plot.index[2]],
        #                pch=20,col="steelblue",cex=2)
        #         # legend("topleft",fill = c("darkviolet","steelblue"),legend=c("galon Macro","galon CD8"))
        mgenes_extended <- c("Tnf","Mmp12","Adam8","Mrc1","Cd36","Cd83","Itgam","Lyz1","Slamf8","Clec12a","Clec10a","Pdc","Cd274")
        mgenes_extended_ENS <- rownames(cm2)[match(mgenes_extended,cm2$fromgtf)]
        points(pcx[rownames(pcx) %in% mgenes_extended_ENS,plot.index[1]],
               pcx[rownames(pcx) %in% mgenes_extended_ENS,plot.index[2]],
               pch=20,col="salmon",cex=2)
        ## ## ##

      })

      # output$debug <- reactive({      cat(paste(input$pca_x_G,input$pca_y_G))    })


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

        pcspcs$geneName <- cm2$fromgtf[match(pcspcs$geneID,rownames(cm2))]


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





      output$searchresult <- renderPrint({

        if(is.null(input$color_by_G)) return("Select a factor to plot your gene")

        if(is.null(input$searchgene))
          return("Type in the gene name/id you want to plot")


        anno <- cm2
        foundGeneID <- input$searchgene %in% rownames(obj)
        foundGeneName <- input$searchgene %in% anno$fromgtf
        if(!foundGeneID){
          foundGeneID <- toupper(input$searchgene) %in% toupper(rownames(obj))
          if(foundGeneID){
            return(paste0("Maybe you mis-spelled the name of your gene. Did you mean ",
                          unique(rownames(anno)[which(toupper(input$searchgene)==rownames(anno))]),"?")) # todo
          } else {
            foundGeneNAME <- input$searchgene %in% anno$fromgtf
            if(!foundGeneNAME){
              foundGeneNAME <- toupper(input$searchgene) %in% toupper(anno$fromgtf)
              if(foundGeneNAME){
                return(paste0("Maybe you mis-spelled the name of your gene. Did you mean ",
                              unique(anno$fromgtf[which(toupper(input$searchgene)==anno$fromgtf)]),"?"))
              } else {return("Could not find the gene you typed!")}
            } else {
              fgn <- anno$fromgtf[which(anno$fromgtf==input$searchgene)]
              if (length(fgn) > 1) return(paste0("Found more than one gene with the selected gene name. Select one of the following: ",paste(selectedGene,collapse=", ")))
              selectedGene <- rownames(obj)[which(cm2$fromgtf==input$searchgene)]

              fg <- rownames(obj)[match(fgn,anno$fromgtf)]
              return(paste0("I found the gene! Plotting ", fg, " - ", cm2$fromgtf[match(fg,rownames(cm2))],"..."))

            }}
        } else {
          fg <- rownames(obj)[match(input$searchgene,rownames(obj))]
          return(paste0("I found the gene! Plotting ", fg, " - ", cm2$fromgtf[match(fg,rownames(cm2))],"..."))

        }


        # foundGeneNAME <- input$searchgene %in% anno$fromgtf



      })

      ## todo: do it more elegantly with a double check on id + gene names?

      output$debuggene <- renderPrint({
        #         if(is.null(input$searchgene)) return(ggplot() + annotate("text",label="Type in a gene name/id",0,0) + theme_bw())
        #         if(as.character(input$searchgene) %in% rownames(obj)) return(ggplot() + annotate("text",label="gene not found...",0,0) + theme_bw())
        #
        #         selectedGene <- rownames(obj)[match(input$searchgene,rownames(obj))]
        #         selectedGeneSymbol <- cm2$fromgtf[match(selectedGene,rownames(cm2))]
        # #         # plotCounts(dds_cleaner,)
        # #         genedata <- plotCounts(obj2,gene=selectedGene,intgroup = input$color_by_G,returnData = T)
        #         as.character(selectedGene)
      })


      output$searchgenePlot <- renderPlot({
        annotation <- cm2

        anno_id <- rownames(cm2)
        anno_gene <- annotation$fromgtf

        if(is.null(input$color_by_G)) return(ggplot() + annotate("text",label="Select a factor to plot your gene",0,0) + theme_bw())

        if(is.null(input$searchgene)) return(ggplot() + annotate("text",label="Type in a gene name/id",0,0) + theme_bw())
        if(!input$searchgene %in% anno_gene & !input$searchgene %in% anno_id) return(ggplot() + annotate("text",label="gene not found...",0,0) + theme_bw())
        # rownames(obj)[match(as.character(input$searchgene),rownames(obj))]




        if (input$searchgene %in% anno_id) {
          selectedGene <- rownames(obj)[match(input$searchgene,rownames(obj))]
          selectedGeneSymbol <- cm2$fromgtf[match(selectedGene,rownames(cm2))]
        }

        if (input$searchgene %in% anno_gene) {
          selectedGeneSymbol <- cm2$fromgtf[which(cm2$fromgtf==input$searchgene)]
          if (length(selectedGeneSymbol) > 1) return(ggplot() + annotate("text",label=paste0("Type in a gene name/id of the following:\n",paste(selectedGene,collapse=", ")),0,0) + theme_bw())
          selectedGene <- rownames(obj)[which(cm2$fromgtf==input$searchgene)]
        }
        # plotCounts(dds_cleaner,)
        genedata <- plotCounts(obj2,gene=selectedGene,intgroup = input$color_by_G,returnData = T)

        # onlyfactors <- genedata[match(input$color_by_G,colnames(genedata))]
        onlyfactors <- genedata[,match(input$color_by_G,colnames(genedata))]

        ## intgroup can be a vector of factors. then i need interactions of the two factors
        # plotCounts(ddsmf_global,gene="ENSMUSG00000026080",intgroup = c("tissue","condition"),returnData = T) -> dh
        # dh$f1f2 <- interaction(dh$tissue,dh$condition)
        # dh  %>% ggplot(aes(x=f1f2,y=count,fill=f1f2)) + geom_boxplot()


        ## TODO: make the intgroup/colr by also somehow multiple-selectable?

        # genedata$plotby <- lapply(1:ncol(onlyfactors),function(arg) onlyfactors[,arg]) %>% interaction()
        genedata$plotby <- interaction(onlyfactors)



        p <- ggplot(genedata,aes(x=plotby,y=count,fill=plotby)) + geom_boxplot() + labs(title=paste0("Normalized counts for ",selectedGeneSymbol," - ",selectedGene)) +  scale_x_discrete(name="") + geom_jitter(aes(x=plotby,y=count),position = position_jitter(width = 0.1)) + scale_fill_discrete(name="Experimental\nconditions")

        if(input$ylimZero)
        {
          p <- p + scale_y_log10(name="Normalized counts - log10 scale",limits=c(1,max(genedata$count)))
        } else {
          p <- p + scale_y_log10(name="Normalized counts - log10 scale")
        }
        p

        # exportPlots$genesZoom <- res
        # res

      })







    })



  shinyApp(ui = uiui, server = serser)
  }



