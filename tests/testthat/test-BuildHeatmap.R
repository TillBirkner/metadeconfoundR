library(testthat)
library(metadeconfoundR)  # Your package name

test_that("Function `BuildHeatmap()` works correctly", {

  #feature <- reduced_feature
  #metaMat <- metaMatMetformin
  input <- readRDS("example_output.rds")

  # 2024 09 24 saveRDS(heatmapPlot, "tests/testthat/expected_BuildHeatmap_output.rds")
  expected_output <- readRDS("expected_BuildHeatmap_output.rds")
  #expected_output <- heatmapPlot


  # Call the function
  result <- BuildHeatmap(input)

  # Verify that the result matches the expected output
  expect_equal(result$data, expected_output$data)
})
