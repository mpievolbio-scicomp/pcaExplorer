library(pcaExplorer)

context("Check that correlations are computed and plotted")

library(DESeq2)
dds <- makeExampleDESeqDataSet_multifac(betaSD_condition = 3,betaSD_tissue = 1)
rlt <- rlogTransformation(dds)
pcaobj <- prcomp(t(assay(rlt)))

res <- correlatePCs(pcaobj,colData(dds))
expect_equal(dim(res),c(4,2))
expect_equal(colnames(res),colnames(colData(dds)))
plotPCcorrs(res)

plotPCcorrs(res,logp = FALSE,pc = 2)
