## interpreting PCA.--..

# groups <- colData(dds_cleaner)$condition
# bbgg <- rownames(dds_deplall)[rowSums(counts(dds_deplall))>0]
# seb_pca2go_unscaled <- pca2go(rld_deplall,annotation=annotation,ensToGeneSymbol = T,scale=F,background_genes=bbgg)
# seb_pca2go_scaled <- pca2go(rld_deplall,annotation=annotation,ensToGeneSymbol = T,scale=T,background_genes=bbgg)


#' Functional interpretation of the principal components
#'
#' Extracts the genes with the highest loadings for each principal component, and
#' performs functional enrichment analysis on them using routines and algorithms from
#' the \code{topGO} package
#'
#' @param se A \code{\link{DESeqTransform}} object, with data in \code{assay(se)},
#' produced for example by either \code{\link{rlog}} or
#' \code{\link{varianceStabilizingTransformation}}
#' @param pca_ngenes Number of genes to use for the PCA
#' @param annotation A \code{data.frame} object, with row.names as gene identifiers (e.g. ENSEMBL ids)
#' and a column, \code{gene_name}, containing e.g. HGNC-based gene symbols
#' @param inputType Input format type of the gene identifiers. Will be used by the routines of \code{topGO}
#' @param organism Character abbreviation for the species, using \code{org.XX.eg.db} for annotation
#' @param ensToGeneSymbol Logical, whether to expect ENSEMBL gene identifiers, to convert to gene symbols
#' with the \code{annotation} provided
#' @param loadings_ngenes Number of genes to extract the loadings (in each direction)
#' @param background_genes Which genes to consider as background.
#' @param scale Logical, defaults to FALSE, scale values for the PCA
#' @param ... Further parameters to be passed to the topGO routine
#'
#' @return A nested list object containing for each principal component the terms enriched
#' in each direction. This object is to be thought in combination with the displaying feature
#' of the main \code{\link{pcaExplorer}} function
#'
#' @examples
#'
#' @export
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
    if(is(se,"DESeqDataSet")) {
      BGids <- rownames(se)[rowSums(counts(se))>0]
    } else if(is(se,"DESeqTransform")){
      BGids <- rownames(se)[rowSums(assay(se))!=0]
    } else {
      BGids <- rownames(se)
    }
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
  topGOpc1pos <- topGOtable(probesPC1pos, BGids, annot = annFUN.org,mapping = annopkg,...)
  topGOpc1neg <- topGOtable(probesPC1neg, BGids, annot = annFUN.org,mapping = annopkg,...)
  topGOpc2pos <- topGOtable(probesPC2pos, BGids, annot = annFUN.org,mapping = annopkg,...)
  topGOpc2neg <- topGOtable(probesPC2neg, BGids, annot = annFUN.org,mapping = annopkg,...)
  topGOpc3pos <- topGOtable(probesPC3pos, BGids, annot = annFUN.org,mapping = annopkg,...)
  topGOpc3neg <- topGOtable(probesPC3neg, BGids, annot = annFUN.org,mapping = annopkg,...)
  topGOpc4pos <- topGOtable(probesPC4pos, BGids, annot = annFUN.org,mapping = annopkg,...)
  topGOpc4neg <- topGOtable(probesPC4neg, BGids, annot = annFUN.org,mapping = annopkg,...)

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






#' Extract functional terms enriched in the DE genes, based on topGO
#'
#' A wrapper for extracting functional GO terms enriched in the DE genes, based on
#' the algorithm and the implementation in the topGO package
#'
#' @param DEgenes A vector of (differentially expressed) genes
#' @param BGgenes A vector of background genes, e.g. all (expressed) genes in the assays
#' @param ontology Which Gene Ontology domain to analyze: \code{BP} (Biological Process), \code{MF} (Molecular Function), or \code{CC} (Cellular Component)
#' @param maxP Max p-value, used to subset the returned table
#' @param annot Which function to use for annotating genes to GO terms
#' @param mapping Which \code{org.XX.eg.db} to use for annotation - select according to the species
#' @param geneID Which format the genes are provided
#' @param topTablerows How many rows to report before any filtering
#' @param fullNamesInRows Logical, whether to display or not the full names for the GO terms
#' @param addGeneToTerms Logical, whether to add a column with all genes annotated to each GO term
#' @param plotGraph Logical, if TRUE additionally plots a graph on the identified GO terms
#' @param plotNodes Number of nodes to plot
#' @param writeOutput Logical, if TRUE additionally writes out the result to a file
#' @param outputFile Name of the file the result should be written into
#'
#' @import topGO
#'
#' @return A table containing the computed GO Terms and related enrichment scores
#'
#' @examples # An example
#'
#'
#' @export
topGOtable <- function(DEgenes,                  # Differentially expressed genes
                       BGgenes,                 # background genes, normally = rownames(cds) or filtering to genes
                       #  with at least 1 read - could also be ls(org.Mm.egGO)
                       ontology="BP",            # could use also "MF"
                       maxP = 0.001,             # use to subset the final table
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
    SignificantGenes <- sigGenes(GOdata)
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
    if(is(se,"DESeqDataSet")) {
      BGids <- rownames(se)[rowSums(counts(se))>0]
    } else if(is(se,"DESeqTransform")){
      BGids <- rownames(se)[rowSums(assay(se))!=0]
    } else {
      BGids <- rownames(se)
    }
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



#' Functional interpretation of the principal components, based on simple
#' overrepresentation analysis
#'
#' Extracts the genes with the highest loadings for each principal component, and
#' performs functional enrichment analysis on them using the simple and quick routine
#' provided by the \code{limma} package
#'
#'
#' @param se A \code{\link{DESeqTransform}} object, with data in \code{assay(se)},
#' produced for example by either \code{\link{rlog}} or
#' \code{\link{varianceStabilizingTransformation}}
#' @param pca_ngenes Number of genes to use for the PCA
#' @param inputType Input format type of the gene identifiers. Deafults to \code{ENSEMBL}, that then will
#' be converted to ENTREZ ids. Can assume values such as \code{ENTREZID},\code{GENENAME} or \code{SYMBOL},
#' like it is normally used with the \code{select} function of \code{AnnotationDbi}
#' @param organism Character abbreviation for the species, using \code{org.XX.eg.db} for annotation
#' @param loadings_ngenes Number of genes to extract the loadings (in each direction)
#' @param background_genes Which genes to consider as background.
#' @param scale Logical, defaults to FALSE, scale values for the PCA
#' @param ... Further parameters to be passed to the topGO routine
#'
#' @return A nested list object containing for each principal component the terms enriched
#' in each direction. This object is to be thought in combination with the displaying feature
#' of the main \code{\link{pcaExplorer}} function
#'
#' @examples
#'
#' @export
limmaquickpca2go <- function(se,
                        pca_ngenes = 10000,
                        inputType = "ENSEMBL",
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
    if(is(se,"DESeqDataSet")) {
      BGids <- rownames(se)[rowSums(counts(se))>0]
    } else if(is(se,"DESeqTransform")){
      BGids <- rownames(se)[rowSums(assay(se))!=0]
    } else {
      BGids <- rownames(se)
    }
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

  probesPC1pos_ENTREZ <- AnnotationDbi::select(eval(parse(text=annopkg)), keys = probesPC1pos, columns=c("ENSEMBL","ENTREZID","GENENAME","SYMBOL"), keytype=inputType)
  probesPC1neg_ENTREZ <- AnnotationDbi::select(eval(parse(text=annopkg)), keys = probesPC1neg, columns=c("ENSEMBL","ENTREZID","GENENAME","SYMBOL"), keytype=inputType)
  probesPC2pos_ENTREZ <- AnnotationDbi::select(eval(parse(text=annopkg)), keys = probesPC2pos, columns=c("ENSEMBL","ENTREZID","GENENAME","SYMBOL"), keytype=inputType)
  probesPC2neg_ENTREZ <- AnnotationDbi::select(eval(parse(text=annopkg)), keys = probesPC2neg, columns=c("ENSEMBL","ENTREZID","GENENAME","SYMBOL"), keytype=inputType)
  probesPC3pos_ENTREZ <- AnnotationDbi::select(eval(parse(text=annopkg)), keys = probesPC3pos, columns=c("ENSEMBL","ENTREZID","GENENAME","SYMBOL"), keytype=inputType)
  probesPC3neg_ENTREZ <- AnnotationDbi::select(eval(parse(text=annopkg)), keys = probesPC3neg, columns=c("ENSEMBL","ENTREZID","GENENAME","SYMBOL"), keytype=inputType)
  probesPC4pos_ENTREZ <- AnnotationDbi::select(eval(parse(text=annopkg)), keys = probesPC4pos, columns=c("ENSEMBL","ENTREZID","GENENAME","SYMBOL"), keytype=inputType)
  probesPC4neg_ENTREZ <- AnnotationDbi::select(eval(parse(text=annopkg)), keys = probesPC4neg, columns=c("ENSEMBL","ENTREZID","GENENAME","SYMBOL"), keytype=inputType)
  bg_ENTREZ <- AnnotationDbi::select(eval(parse(text=annopkg)), keys = BGids, columns=c("ENSEMBL","ENTREZID","GENENAME","SYMBOL"), keytype=inputType)
  print("Ranking genes by the loadings ... done!")
  # print("Extracting functional categories enriched in the gene subsets ...")





  #   goParams <- new("GOHyperGParams",geneIds = tt$ENTREZID,universeGeneIds = uu$ENTREZID,annotation ="org.Mm.eg" ,ontology = "BP",pvalueCutoff = 0.01,conditional = TRUE,testDirection = "over")
  #   goResults <- hyperGTest(goParams)
  #   summary(goResults)
  #   goParams <- new("GOHyperGParams",geneIds = probesPC1pos_ENTREZ,universeGeneIds = bg_ENTREZ,annotation ="org.Mm.eg" ,ontology = "BP",pvalueCutoff = 0.01,conditional = TRUE,testDirection = "over")
  #   goResults <- hyperGTest(goParams)
  #   summary(goResults)

  print("Extracting functional categories enriched in the gene subsets ...")
  quickGOpc1pos <- topGO(limma::goana(probesPC1pos_ENTREZ$ENTREZID, bg_ENTREZ$ENTREZID, species = organism),ontology="BP",number=200);message("1")
  quickGOpc1neg <- topGO(limma::goana(probesPC1neg_ENTREZ$ENTREZID, bg_ENTREZ$ENTREZID, species = organism),ontology="BP",number=200);message("2")
  quickGOpc2pos <- topGO(limma::goana(probesPC2pos_ENTREZ$ENTREZID, bg_ENTREZ$ENTREZID, species = organism),ontology="BP",number=200);message("3")
  quickGOpc2neg <- topGO(limma::goana(probesPC2neg_ENTREZ$ENTREZID, bg_ENTREZ$ENTREZID, species = organism),ontology="BP",number=200);message("4")
  quickGOpc3pos <- topGO(limma::goana(probesPC3pos_ENTREZ$ENTREZID, bg_ENTREZ$ENTREZID, species = organism),ontology="BP",number=200);message("5")
  quickGOpc3neg <- topGO(limma::goana(probesPC3neg_ENTREZ$ENTREZID, bg_ENTREZ$ENTREZID, species = organism),ontology="BP",number=200);message("6")
  quickGOpc4pos <- topGO(limma::goana(probesPC4pos_ENTREZ$ENTREZID, bg_ENTREZ$ENTREZID, species = organism),ontology="BP",number=200);message("7")
  quickGOpc4neg <- topGO(limma::goana(probesPC4neg_ENTREZ$ENTREZID, bg_ENTREZ$ENTREZID, species = organism),ontology="BP",number=200);message("8")

  quickGOpc1pos <- quickGOpc1pos[order(quickGOpc1pos$P.DE),]
  quickGOpc1neg <- quickGOpc1neg[order(quickGOpc1neg$P.DE),]
  quickGOpc2pos <- quickGOpc2pos[order(quickGOpc2pos$P.DE),]
  quickGOpc2neg <- quickGOpc2neg[order(quickGOpc2neg$P.DE),]
  quickGOpc3pos <- quickGOpc3pos[order(quickGOpc3pos$P.DE),]
  quickGOpc3neg <- quickGOpc3neg[order(quickGOpc3neg$P.DE),]
  quickGOpc4pos <- quickGOpc4pos[order(quickGOpc4pos$P.DE),]
  quickGOpc4neg <- quickGOpc4neg[order(quickGOpc4neg$P.DE),]

  goEnrichs <- list(PC1=list(posLoad=quickGOpc1pos,negLoad=quickGOpc1neg),
                    PC2=list(posLoad=quickGOpc2pos,negLoad=quickGOpc2neg),
                    PC3=list(posLoad=quickGOpc3pos,negLoad=quickGOpc3neg),
                    PC4=list(posLoad=quickGOpc4pos,negLoad=quickGOpc4neg)
  )
  print("Extracting functional categories enriched in the gene subsets ... done!")

  attr(goEnrichs,"n_genesforpca") <- pca_ngenes

  return(goEnrichs)

}







