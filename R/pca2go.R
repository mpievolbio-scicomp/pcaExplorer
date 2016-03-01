## interpreting PCA.--..

# groups <- colData(dds_cleaner)$condition
# bbgg <- rownames(dds_deplall)[rowSums(counts(dds_deplall))>0]
# seb_pca2go_unscaled <- pca2go(rld_deplall,annotation=annotation,ensToGeneSymbol = T,scale=F,background_genes=bbgg)
# seb_pca2go_scaled <- pca2go(rld_deplall,annotation=annotation,ensToGeneSymbol = T,scale=T,background_genes=bbgg)
#' Title
#'
#' @param se
#' @param pca_ngenes
#' @param annotation
#' @param inputType
#' @param organism
#' @param ensToGeneSymbol
#' @param loadings_ngenes
#' @param background_genes
#' @param scale
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
pca2go <- function(se,
                   pca_ngenes = 10000,
                   annotation = NULL,
                   inputType = "geneSymbol",
                   organism = "Mm",
                   ensToGeneSymbol = FALSE,
                   loadings_ngenes = 500,
                   background_genes = NULL,
                   scale = FALSE,
                   ... # further parameters to be passed to the topgo routine

                   ) {

  annopkg <- paste0("org.",organism,".eg.db")
  if (!require(annopkg,character.only=TRUE)) {
    stop("The package",annopkg, "is not installed/available. Try installing it with biocLite() ?")
  }
  exprsData <- assay(se)

  if(is.null(background_genes)) {
    BGids <- rownames(se)[rowSums(counts(se))>0] # TODO: at best some other way of doing it - maybe passing additionally the DDS? or if nothing provided, then use all rownames...
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

  p <- prcomp(t(exprsData), scale=scale, center=TRUE)
  # pcaobj <- list(scores=p$x, loadings=p$rotation, pov=p$sdev^2/sum(p$sdev^2),
                 # expressionData=NA)
  # class(pcaobj) <- "pca" # to view it eventually with pcaGoPromoter
  # res

  # library("pcaGoPromoter")
  # pcaObj <- res


  print("Ranking genes by the loadings ...")
  # library("pcaGoPromoter")
  probesPC1pos <- rankedGeneLoadings(p, pc=1,decreasing=TRUE)[1:loadings_ngenes]
  probesPC1neg <- rankedGeneLoadings(p, pc=1,decreasing=FALSE)[1:loadings_ngenes]
  probesPC2pos <- rankedGeneLoadings(p, pc=2,decreasing=TRUE)[1:loadings_ngenes]
  probesPC2neg <- rankedGeneLoadings(p, pc=2,decreasing=FALSE)[1:loadings_ngenes]
  probesPC3pos <- rankedGeneLoadings(p, pc=3,decreasing=TRUE)[1:loadings_ngenes]
  probesPC3neg <- rankedGeneLoadings(p, pc=3,decreasing=FALSE)[1:loadings_ngenes]
  probesPC4pos <- rankedGeneLoadings(p, pc=4,decreasing=TRUE)[1:loadings_ngenes]
  probesPC4neg <- rankedGeneLoadings(p, pc=4,decreasing=FALSE)[1:loadings_ngenes]



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

  attr(goEnrichs,"n_genesforpca") <- pca_ngenes

  return(goEnrichs)
}





rankedGeneLoadings <- function (x, pc = 1, decreasing = TRUE)
{
  # works directly on the prcomp object
  return(rownames(x$rotation)[order(x$rotation[, pc], decreasing = decreasing)])
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
    # library("GO.db")
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
  if(writeOutput) write.table(sTab,file=outputFile,sep="\t",quote=FALSE,col.names=TRUE,row.names=FALSE)
  if(plotGraph) showSigOfNodes(GOdata,topGO::score(resultFisher),firstSigNodes=plotNodes, useInfo="all")
  #   if(outputToLatex) sTabSig <- xtable(apply(sTabSig[1:15,], 2, as.character)) # take a smaller subset

  # and returns the significant ones # or all, like here
  return(sTab)
}




quickpca2go <- function(se,
                    pca_ngenes = 10000,
                    annotation = NULL,
                    inputType = "ensGene",
                    organism = "Mm",
                    loadings_ngenes = 500,
                    background_genes = NULL,
                    scale = FALSE,
                    ... # further parameters to be passed to the topgo routine

){
  annopkg <- paste0("org.",organism,".eg.db")
  if (!require(annopkg,character.only=TRUE)) {
    stop("The package",annopkg, "is not installed/available. Try installing it with biocLite() ?")
  }
  exprsData <- assay(se)

  if(is.null(background_genes)) {
    BGids <- rownames(se)[rowSums(counts(se))>0] # TODO: at best some other way of doing it - maybe passing additionally the DDS? or if nothing provided, then use all rownames...
  } else {
    BGids <- background_genes
  }

  exprsData <- assay(se)[order(rowVars(assay(se)),decreasing=TRUE),]
  exprsData <- exprsData[1:pca_ngenes,]



  rv <- rowVars(exprsData)
  dropped <- sum(rv==0)
  if (dropped > 0)
    print(paste("Dropped", dropped, "genes with 0 variance"))

  exprsData <- exprsData[rv>0,]

  message("After subsetting/filtering for invariant genes, working on a ",nrow(exprsData),"x",ncol(exprsData)," expression matrix\n")

  p <- prcomp(t(exprsData), scale=scale, center=TRUE)


  print("Ranking genes by the loadings ...")
  # library("pcaGoPromoter")
  probesPC1pos <- rankedGeneLoadings(p, pc=1,decreasing=TRUE)[1:loadings_ngenes]
  probesPC1neg <- rankedGeneLoadings(p, pc=1,decreasing=FALSE)[1:loadings_ngenes]
  probesPC2pos <- rankedGeneLoadings(p, pc=2,decreasing=TRUE)[1:loadings_ngenes]
  probesPC2neg <- rankedGeneLoadings(p, pc=2,decreasing=FALSE)[1:loadings_ngenes]
  probesPC3pos <- rankedGeneLoadings(p, pc=3,decreasing=TRUE)[1:loadings_ngenes]
  probesPC3neg <- rankedGeneLoadings(p, pc=3,decreasing=FALSE)[1:loadings_ngenes]
  probesPC4pos <- rankedGeneLoadings(p, pc=4,decreasing=TRUE)[1:loadings_ngenes]
  probesPC4neg <- rankedGeneLoadings(p, pc=4,decreasing=FALSE)[1:loadings_ngenes]

  ## convert ensembl to entrez ids

  probesPC1pos_ENTREZ <- AnnotationDbi::select(eval(parse(text=annopkg)), keys = probesPC1pos, columns=c("ENSEMBL","ENTREZID","GENENAME","SYMBOL"), keytype="ENSEMBL")
  probesPC1neg_ENTREZ <- AnnotationDbi::select(eval(parse(text=annopkg)), keys = probesPC1neg, columns=c("ENSEMBL","ENTREZID","GENENAME","SYMBOL"), keytype="ENSEMBL")
  probesPC2pos_ENTREZ <- AnnotationDbi::select(eval(parse(text=annopkg)), keys = probesPC2pos, columns=c("ENSEMBL","ENTREZID","GENENAME","SYMBOL"), keytype="ENSEMBL")
  probesPC2neg_ENTREZ <- AnnotationDbi::select(eval(parse(text=annopkg)), keys = probesPC2neg, columns=c("ENSEMBL","ENTREZID","GENENAME","SYMBOL"), keytype="ENSEMBL")
  probesPC3pos_ENTREZ <- AnnotationDbi::select(eval(parse(text=annopkg)), keys = probesPC3pos, columns=c("ENSEMBL","ENTREZID","GENENAME","SYMBOL"), keytype="ENSEMBL")
  probesPC3neg_ENTREZ <- AnnotationDbi::select(eval(parse(text=annopkg)), keys = probesPC3neg, columns=c("ENSEMBL","ENTREZID","GENENAME","SYMBOL"), keytype="ENSEMBL")
  probesPC4pos_ENTREZ <- AnnotationDbi::select(eval(parse(text=annopkg)), keys = probesPC4pos, columns=c("ENSEMBL","ENTREZID","GENENAME","SYMBOL"), keytype="ENSEMBL")
  probesPC4neg_ENTREZ <- AnnotationDbi::select(eval(parse(text=annopkg)), keys = probesPC4neg, columns=c("ENSEMBL","ENTREZID","GENENAME","SYMBOL"), keytype="ENSEMBL")
  bg_ENTREZ <- AnnotationDbi::select(eval(parse(text=annopkg)), keys = BGids, columns=c("ENSEMBL","ENTREZID","GENENAME","SYMBOL"), keytype="ENSEMBL")
  print("Ranking genes by the loadings ... done!")
  # print("Extracting functional categories enriched in the gene subsets ...")





#   goParams <- new("GOHyperGParams",geneIds = tt$ENTREZID,universeGeneIds = uu$ENTREZID,annotation ="org.Mm.eg" ,ontology = "BP",pvalueCutoff = 0.01,conditional = TRUE,testDirection = "over")
#   goResults <- hyperGTest(goParams)
#   summary(goResults)
#   goParams <- new("GOHyperGParams",geneIds = probesPC1pos_ENTREZ,universeGeneIds = bg_ENTREZ,annotation ="org.Mm.eg" ,ontology = "BP",pvalueCutoff = 0.01,conditional = TRUE,testDirection = "over")
#   goResults <- hyperGTest(goParams)
#   summary(goResults)

  print("Extracting functional categories enriched in the gene subsets ...")
  quickGOpc1pos <- GOenrich(probesPC1pos_ENTREZ$ENTREZID, bg_ENTREZ$ENTREZID, annotation = gsub(".db","",annopkg), ontology = "BP",pvalueCutoff = 0.01, conditional = TRUE)
  quickGOpc1neg <- GOenrich(probesPC1neg_ENTREZ$ENTREZID, bg_ENTREZ$ENTREZID, annotation = gsub(".db","",annopkg), ontology = "BP",pvalueCutoff = 0.01, conditional = TRUE)
  quickGOpc2pos <- GOenrich(probesPC2pos_ENTREZ$ENTREZID, bg_ENTREZ$ENTREZID, annotation = gsub(".db","",annopkg), ontology = "BP",pvalueCutoff = 0.01, conditional = TRUE)
  quickGOpc2neg <- GOenrich(probesPC2neg_ENTREZ$ENTREZID, bg_ENTREZ$ENTREZID, annotation = gsub(".db","",annopkg), ontology = "BP",pvalueCutoff = 0.01, conditional = TRUE)
  quickGOpc3pos <- GOenrich(probesPC3pos_ENTREZ$ENTREZID, bg_ENTREZ$ENTREZID, annotation = gsub(".db","",annopkg), ontology = "BP",pvalueCutoff = 0.01, conditional = TRUE)
  quickGOpc3neg <- GOenrich(probesPC3neg_ENTREZ$ENTREZID, bg_ENTREZ$ENTREZID, annotation = gsub(".db","",annopkg), ontology = "BP",pvalueCutoff = 0.01, conditional = TRUE)
  quickGOpc4pos <- GOenrich(probesPC4pos_ENTREZ$ENTREZID, bg_ENTREZ$ENTREZID, annotation = gsub(".db","",annopkg), ontology = "BP",pvalueCutoff = 0.01, conditional = TRUE)
  quickGOpc4neg <- GOenrich(probesPC4neg_ENTREZ$ENTREZID, bg_ENTREZ$ENTREZID, annotation = gsub(".db","",annopkg), ontology = "BP",pvalueCutoff = 0.01, conditional = TRUE)

  goEnrichs <- list(PC1=list(posLoad=quickGOpc1pos,negLoad=quickGOpc1neg),
                    PC2=list(posLoad=quickGOpc2pos,negLoad=quickGOpc2neg),
                    PC3=list(posLoad=quickGOpc3pos,negLoad=quickGOpc3neg),
                    PC4=list(posLoad=quickGOpc4pos,negLoad=quickGOpc4neg)
  )
  print("Extracting functional categories enriched in the gene subsets ... done!")

  attr(goEnrichs,"n_genesforpca") <- pca_ngenes

  return(goEnrichs)

}


GOenrich <- function(geneIds,universeGeneIds,annotation,ontology="BP",pvalueCutoff=0.01,conditional=TRUE,testDirection="over"){
  goParams <- new("GOHyperGParams",geneIds = geneIds,universeGeneIds = universeGeneIds,annotation =annotation ,ontology = ontology,pvalueCutoff = pvalueCutoff,conditional = conditional,testDirection = testDirection)
  goResults <- hyperGTest(goParams)
  summary(goResults)
}









