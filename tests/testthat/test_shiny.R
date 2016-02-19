library("pcaExplorer")

context("Check that shiny app is generated")

test_that("Shiny app is generated", {
  expect_is(pcaExplorer(), "shiny.appobj")
  # expect_is(pcaExplorer(), "shiny.appobj")
  # expect_is(pcaExplorer(), "shiny.appobj")
  # expect_is(pcaExplorer(), "shiny.appobj")

  # expect_equal(2+2,4)
  # expect_is(COBRAapp(cobradata_example), "shiny.appobj")

  # expect_is(COBRAapp(autorun = TRUE), "shiny.appobj")
  # expect_is(COBRAapp(cobradata_example, autorun = TRUE), "shiny.appobj")
})