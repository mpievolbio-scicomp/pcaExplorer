library(pcaExplorer)

context("Checks on the pca on the genes")

library(DESeq2)
dds <- makeExampleDESeqDataSet_multifac(betaSD_condition = 3,betaSD_tissue = 1)
rlt <- rlogTransformation(dds)

groups <- colData(dds)$condition
cols <- scales::hue_pal()(2)[groups]
genespca(rlt,ntop=100,arrowColors=cols,groupNames=groups)
dat <- genespca(rlt,ntop=100,arrowColors=cols,groupNames=groups,returnData = TRUE)


genespca(rlt,ntop=100)

genespca(rlt,ntop=100,arrowColors = "green")

expect_error(genespca(rlt,ntop=100,arrowColors = c("green","red")))

groups_multi <- interaction(as.data.frame(colData(rlt)[,c("condition","tissue")]))
cols_multi <- scales::hue_pal()(length(levels(groups_multi)))[factor(groups_multi)]
genespca(rlt,ntop=100,arrowColors=cols_multi,groupNames=groups_multi)