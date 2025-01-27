library(testthat)
library(metadeconfoundR)
library(ggraph)

test_that("Function `BuildHeatmap()` works correctly", {
  feature <- reduced_feature
  metaMat <- metaMatMetformin

  input <- read.table("2024_10_17_example_output.tsv", header = T, sep = "\t")
  input$feature <- as.factor(input$feature)
  input$metaVariable <- as.factor(input$metaVariable)
  expect_no_error(BuildConfounderMap(input))


  expect_warning(
    BuildConfounderMap(
      input,
      metaVariableNames = as.data.frame(cbind(colnames(metaMat), paste0(colnames(metaMat), "_nice"))),
      featureNames = cbind(colnames(feature), paste0(colnames(feature), "_nice"))
      ),
    'class(featureNames) was coerced to "data.frame"',
    fixed = TRUE
    )



  example_output_wide <- readRDS("2024_10_17_example_output_wide.rds")
  # plot from wide format data
  expect_no_error(
    BuildConfounderMap(
      example_output_wide
    )
  )

  # wrong number of colors
  expect_error(
    BuildConfounderMap(
      example_output_wide,
      d_col = c("red", "blue")
    ),
    "wrong number of colors in d_col!\nSupply colors for c(min, middle, max)!",
    fixed = TRUE
  )


  # no trusted
  expect_error(
    BuildHeatmap(
      example_output_wide,
      trusted = c()
    ),
    '"trusted" must contain at least one trusted status label',
    fixed = TRUE
  )

})

