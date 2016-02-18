
x <- rld_deplall
rv <- rowVars(assay(x))
select <- order(rv, decreasing = TRUE)[seq_len(min(nrow(x),length(rv)))]
pcaobj <- prcomp(t(assay(x)[select, ]),scale = F, center = T)



hi_loadings <- function(pcaobj, whichpc = 1, topN = 10, exprTable = NULL,
                        annotation = NULL, title="Top/bottom loadings - "){

  geneloadings_sorted <- sort(pcaobj$rotation[,whichpc])
  geneloadings_extreme <- c(tail(geneloadings_sorted,topN),head(geneloadings_sorted,topN))

  if(!is.null(exprTable)) {
    tab <- exprTable[names(geneloadings_extreme),]
    if(!is.null(annotation))
      rownames(tab) <- annotation$gene_name[match(rownames(tab),rownames(annotation))]
    return(tab)
  }

  if(!is.null(annotation))
    names(geneloadings_extreme) <- annotation$gene_name[match(names(geneloadings_extreme),rownames(annotation))]

  barplot(geneloadings_extreme, las=2,main=paste0(title, "PC", whichpc))
  # p <- ggplot(data.frame(loadings=geneloadings_extreme,geneID=names(geneloadings_extreme)),aes(geneID,weight=loadings)) + geom_bar()
}




#     par(mfrow=c(3,1))
#   for(i in 1:3){
#     load <- p$rotation[,i][order(p$rotation[,i])]
#     extreme <- c(tail(load), head(load))
#     extreme.ensg <- names(extreme)
#     ensembl = useMart("ensembl", dataset = "hsapiens_gene_ensembl") #select the ensembl database
#     extreme.symbols <- getBM(attributes=c("ensembl_gene_id", "hgnc_symbol"),
#                              filters = "ensembl_gene_id",
#                              values=extreme.ensg,
#                              mart=ensembl)
#     q <- extreme.symbols[,2]
#     names(q) <- extreme.symbols[,1]
#     fpkm <- cbind(q[extreme.ensg],fpkm.table[extreme.ensg,])
#     names(fpkm)[names(fpkm) == 'q[extreme.ensg]'] <- 'Gene Symbol'
#     barplot(extreme, names.arg=q[extreme.ensg],las=2,main=paste0(caption, ", PC", i))
#     print(fpkm)
#   }
# }