## interpreting PCA.--..

# groups <- colData(dds_cleaner)$condition


pca2go <- function(se,
                   pca_ngenes = 10000,
                   annotation = NULL,
                   inputType = "geneSymbol",
                   organism = "Mm",
                   ensToGeneSymbol = FALSE,
                   loadings_ngenes = 500,
                   background_genes = NULL,
                   ... # further parameters to be passed to the topgo routine

                   ) {

  annopkg <- paste0("org.",organism,".eg.db")
  if (!require(annopkg,character.only=T)) {
    stop("The package",annopkg, "is not installed/available.")
  }
  exprsData <- assay(se)

  if(is.null(background_genes)) {
    BGids <- rownames(se)[rowSums(counts(se))>0]
  } else {
    BGids <- background_genes
  }

  exprsData <- assay(se)[order(rowVars(assay(se)),decreasing=TRUE),]
  exprsData <- exprsData[1:pca_ngenes,]

  ## convert ensembl to gene symbol ids
  if(ensToGeneSymbol & !(is.null(annotation))) {
    rownames(exprsData) <- annotation$gene_name[match(rownames(exprsData),rownames(annotation))]
    BGids <- annotation$gene_name[match(BGids,rownames(annotation))]
  }

  rv <- rowVars(exprsData)
  dropped <- sum(rv==0)
  if (dropped > 0)
    print(paste("Dropped", dropped, "genes with 0 variance"))

  exprsData <- exprsData[rv>0,]

  message("After subsetting/filtering for invariant genes, working on a ",nrow(exprsData),"x",ncol(exprsData)," expression matrix\n")

  p <- prcomp(t(exprsData), scale=TRUE, center=TRUE)
  pcaobj <- list(scores=p$x, loadings=p$rotation, pov=p$sdev^2/sum(p$sdev^2),
                 expressionData=NA)
  class(pcaobj) <- "pca" # to view it eventually with pcaGoPromoter
  # res

  # library("pcaGoPromoter")
  # pcaObj <- res


  print("Ranking genes by the loadings ...")
  # library("pcaGoPromoter")
  probesPC1pos <- rankedGeneLoadings(pcaobj, pc=1,decreasing=TRUE)[1:loadings_ngenes]
  probesPC1neg <- rankedGeneLoadings(pcaobj, pc=1,decreasing=FALSE)[1:loadings_ngenes]
  probesPC2pos <- rankedGeneLoadings(pcaobj, pc=2,decreasing=TRUE)[1:loadings_ngenes]
  probesPC2neg <- rankedGeneLoadings(pcaobj, pc=2,decreasing=FALSE)[1:loadings_ngenes]
  probesPC3pos <- rankedGeneLoadings(pcaobj, pc=3,decreasing=TRUE)[1:loadings_ngenes]
  probesPC3neg <- rankedGeneLoadings(pcaobj, pc=3,decreasing=FALSE)[1:loadings_ngenes]
  probesPC4pos <- rankedGeneLoadings(pcaobj, pc=4,decreasing=TRUE)[1:loadings_ngenes]
  probesPC4neg <- rankedGeneLoadings(pcaobj, pc=4,decreasing=FALSE)[1:loadings_ngenes]



  print("Ranking genes by the loadings ... done!")
  print("Extracting functional categories enriched in the gene subsets ...")
  topGOpc1pos <- topGOtable(probesPC1pos, BGids, annot = annFUN.org,mapping = annopkg)
  topGOpc1neg <- topGOtable(probesPC1neg, BGids, annot = annFUN.org,mapping = annopkg)
  topGOpc2pos <- topGOtable(probesPC2pos, BGids, annot = annFUN.org,mapping = annopkg)
  topGOpc2neg <- topGOtable(probesPC2neg, BGids, annot = annFUN.org,mapping = annopkg)
  topGOpc3pos <- topGOtable(probesPC3pos, BGids, annot = annFUN.org,mapping = annopkg)
  topGOpc3neg <- topGOtable(probesPC3neg, BGids, annot = annFUN.org,mapping = annopkg)
  topGOpc4pos <- topGOtable(probesPC4pos, BGids, annot = annFUN.org,mapping = annopkg)
  topGOpc4neg <- topGOtable(probesPC4neg, BGids, annot = annFUN.org,mapping = annopkg)

  goEnrichs <- list(PC1=list(posLoad=topGOpc1pos,negLoad=topGOpc1neg),
                    PC2=list(posLoad=topGOpc2pos,negLoad=topGOpc2neg),
                    PC3=list(posLoad=topGOpc3pos,negLoad=topGOpc3neg),
                    PC4=list(posLoad=topGOpc4pos,negLoad=topGOpc4neg)
  )
  print("Extracting functional categories enriched in the gene subsets ... done!")

  return(goEnrichs)
}





rankedGeneLoadings <- function (x, pc = 1, decreasing = TRUE)
{
  return(rownames(x$loadings)[order(x$loadings[, pc], decreasing = decreasing)])
}






#' Title
#'
#' @param DEgenes
#' @param BGgenes
#' @param ontology
#' @param maxP
#' @param desc
#' @param annot
#' @param mapping
#' @param geneID
#' @param topTablerows
#' @param fullNamesInRows
#' @param plotGraph
#' @param plotNodes
#' @param writeOutput
#' @param outputFile
#' @param addGeneToTerms
#'
#' @return A value
#' @export
#'
#' @examples # An example
topGOtable <- function(DEgenes,                  # Differentially expressed genes
                       BGgenes,                 # background genes, normally = rownames(cds) or filtering to genes
                       #  with at least 1 read - could also be ls(org.Mm.egGO)
                       ontology="BP",            # could use also "MF"
                       maxP = 0.001,             # use to subset the final table
                       desc="",                  # could be used for the output file name
                       annot = annFUN.org,       # parameters for creating topGO object
                       mapping = "org.Mm.eg.db",
                       geneID = "symbol" ,       # could also beID = "entrez")
                       topTablerows = 200,
                       fullNamesInRows = TRUE,
                       addGeneToTerms=TRUE,
                       plotGraph=FALSE, plotNodes= 10,
                       writeOutput=FALSE, outputFile="" #, outputToLatex=FALSE
) {
  # creating the vectors
  DEgenes_input <- factor(as.integer(BGgenes %in% DEgenes))
  names(DEgenes_input) <- BGgenes
  # instantiating the topGOdata object
  GOdata <- new("topGOdata",
                ontology = ontology,
                allGenes = DEgenes_input,
                nodeSize = 10,
                annot = annot,
                mapping = mapping,
                ID = geneID)
  # performing the test(s)
  resultFisher <- runTest(GOdata, algorithm = "elim", statistic = "fisher")
  resultClassic <- runTest(GOdata,algorithm="classic",statistic = "fisher")
  sTab <- GenTable(GOdata,
                   p.value_elim=resultFisher,
                   p.value_classic=resultClassic,
                   orderBy= "p.value_elim",
                   ranksOf= "p.value_classic",
                   topNodes=topTablerows)
  sTabSig <- subset(sTab, as.numeric(p.value_elim) < maxP)

  if(fullNamesInRows){
    library("GO.db")
    sTab$Term <- sapply(sTab$GO.ID ,function(go) { Term(GOTERM[[go]])})
  }

  if(addGeneToTerms) {
    # adapted from an elegant one liner found here: https://support.bioconductor.org/p/65856/
    # SignificantGenes <- genes with a score higher than 0.7 in ENTREZ ID!!
    SignificantGenes <- sigGenes(GOdata)
    # goresults <- GenTable(...your stuff...)
    sTab$genes <- sapply(sTab$GO.ID, function(x)
    {
      genes<-genesInTerm(GOdata, x)
      tmp <- genes[[1]][genes[[1]] %in% SignificantGenes]
    })
    # coerce the list to a comma separated vector
    sTab$genes <- unlist(lapply(sTab$genes,function(arg) paste(arg,collapse=",")))
  }

  # write all entries of the table
  if(writeOutput) write.table(sTab,file=outputFile,sep="\t",quote=F,col.names=T,row.names=F)
  if(plotGraph) showSigOfNodes(GOdata,topGO::score(resultFisher),firstSigNodes=plotNodes, useInfo="all")
  #   if(outputToLatex) sTabSig <- xtable(apply(sTabSig[1:15,], 2, as.character)) # take a smaller subset

  # and returns the significant ones # or all, like here
  return(sTab)
}

