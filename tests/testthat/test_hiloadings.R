library("pcaExplorer")
library("DESeq2")

context("Check that genes with hi loadings are extracted")

dds <- makeExampleDESeqDataSet(n=1000, m=8)
rlt <- rlogTransformation(dds)
pcaobj <- prcomp(t(SummarizedExperiment::assay(rlt)))
anno <- data.frame(gene_id=rownames(dds),gene_name=toupper(rownames(dds)),stringsAsFactors = FALSE,row.names = rownames(dds))

hi_loadings(pcaobj,1)
dev.off()

test_that("Hi loadings are computed", {
  expect_is(hi_loadings(pcaobj,1,exprTable = counts(dds)), "matrix")
  expect_true(all(rownames(hi_loadings(pcaobj,1,exprTable = counts(dds),annotation = NULL)) %in% rownames(dds)))
  expect_true(all(rownames(hi_loadings(pcaobj,1,exprTable = counts(dds),annotation = anno)) %in% anno$gene_name))
})