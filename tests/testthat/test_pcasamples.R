library(pcaExplorer)

context("Checks on the pca on the samples")

library(DESeq2)
dds <- makeExampleDESeqDataSet_multifac(betaSD_condition = 3, betaSD_tissue = 1)
rlt <- rlogTransformation(dds)
pcaobj <- prcomp(t(assay(rlt)))

colData(dds)

pcaplot(rlt)
dat <- pcaplot(rlt, returnData = TRUE)

pcaplot(rlt, intgroup = c("condition", "tissue"))
dat <- pcaplot(rlt, intgroup = c("condition", "tissue"), returnData = TRUE)

expect_error(pcaplot(rlt, intgroup = "foo"))

pcascree(pcaobj)
pcascree(pcaobj, type = "cev")
expect_error(pcascree(pcaobj, type = "foo"))
