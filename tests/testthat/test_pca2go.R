library(pcaExplorer)

context("Checks on the functional enrichment of subset of genes/genes with hi loadings")

# testing this is quite lengthy... explore other possibilities?
library(airway)
library(DESeq2)
data(airway)
airway
dds_airway <- DESeqDataSet(airway, design= ~ cell + dex)

# Example, performing extraction of enriched functional categories in
# detected significantly expressed genes
dds_airway <- DESeq(dds_airway)
res_airway <- results(dds_airway)
library("AnnotationDbi")
library("org.Hs.eg.db")
res_airway$symbol <- mapIds(org.Hs.eg.db,
                            keys=row.names(res_airway),
                            column="SYMBOL",
                            keytype="ENSEMBL",
                            multiVals="first")
res_airway$entrez <- mapIds(org.Hs.eg.db,
                            keys=row.names(res_airway),
                            column="ENTREZID",
                            keytype="ENSEMBL",
                            multiVals="first")
expect_is(res_airway,"DESeqResults")
resOrdered <- as.data.frame(res_airway[order(res_airway$padj),])
de_df <- resOrdered[resOrdered$padj < .05 & !is.na(resOrdered$padj),]
de_symbols <- de_df$symbol
bg_ids <- rownames(dds_airway)[rowSums(counts(dds_airway)) > 0]
bg_symbols <- mapIds(org.Hs.eg.db,
                     keys=bg_ids,
                     column="SYMBOL",
                     keytype="ENSEMBL",
                     multiVals="first")
library(topGO)

expect_is(de_symbols,"character")
expect_is(bg_symbols,"character")

topgoDE_airway <- topGOtable(de_symbols, bg_symbols,
                             ontology = "BP",
                             mapping = "org.Hs.eg.db",
                             geneID = "symbol")

expect_is(topgoDE_airway,"data.frame")


rld_airway <- rlogTransformation(dds_airway)

ngenes_pca <- 10000

goquick_airway <- limmaquickpca2go(rld_airway,
                                   pca_ngenes = ngenes_pca,
                                   inputType = "ENSEMBL",
                                   organism = "Hs")

expect_error(limmaquickpca2go(rld_airway,
                              pca_ngenes = ngenes_pca,
                              inputType = "ENSEMBL",
                              organism = "foo")) # additionally throws a warning

expect_is(goquick_airway,"list")
expect_equal(length(goquick_airway),4) # ensure all pcs are there
sapply(goquick_airway,names)

expect_equal(attr(goquick_airway,"n_genesforpca"),ngenes_pca)


