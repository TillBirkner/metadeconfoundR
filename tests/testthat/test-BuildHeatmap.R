library(testthat)
library(metadeconfoundR)  # Your package name

test_that("Function `BuildHeatmap()` works correctly", {

  feature <- reduced_feature
  metaMat <- metaMatMetformin
  input <- readRDS("example_output.rds")

  # 2024 09 24 saveRDS(heatmapPlot, "tests/testthat/expected_BuildHeatmap_output.rds")
  expected_output <- readRDS("expected_BuildHeatmap_output.rds")
  #expected_output <- heatmapPlot


  # Call the function
  result <- BuildHeatmap(input)

  # Verify that the result matches the expected output
  expect_equal(result$data, expected_output$data)



  expect_warning(
    BuildHeatmap(
      input,
      cuneiform = T,
      keepMeta = colnames(metaMat),
      keepFeature = colnames(feature),
      metaVariableNames = as.data.frame(cbind(colnames(metaMat), paste0(colnames(metaMat), "_nice"))),
      featureNames = cbind(colnames(feature), paste0(colnames(feature), "_nice"))
      ),
    'class(featureNames) was coerced to "data.frame"',
    fixed = TRUE
    )

  # only do intermediate output
  result_intermed <- BuildHeatmap(input, intermedData = T)
  # plot from intermediate output
  expect_warning(
    BuildHeatmap(
      result_intermed
      ),
    "treating input as 'intermedData = T' Buildheatmap output!!"
    )

  example_output_wide <- readRDS("2024_10_17_example_output_wide.rds")
  # plot from wide format data
  expect_no_error(
    BuildHeatmap(
      example_output_wide
    )
  )

})
