library("pcaExplorer")

context("Check that shiny app is generated")

dds <- makeExampleDESeqDataSet(n = 100, m = 8)
rlt <- rlogTransformation(dds)
cm <- counts(dds)
cd <- colData(dds)

test_that("Shiny app is generated", {
  expect_is(pcaExplorer(), "shiny.appobj")
  expect_is(pcaExplorer(dds, rlt), "shiny.appobj")
  expect_is(pcaExplorer(countmatrix = cm, coldata = cd), "shiny.appobj")
  expect_is(pcaExplorer(dds = dds), "shiny.appobj")
})
