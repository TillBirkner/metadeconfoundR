library(testthat)
library(metadeconfoundR)

test_that("Function `BuildHeatmap()` works correctly", {
  feature <- reduced_feature
  metaMat <- metaMatMetformin

  input <- readRDS("example_output.rds")
  # 2024 09 24 saveRDS(heatmapPlot, "tests/testthat/expected_BuildHeatmap_output.rds")
  expected_output <- readRDS("expected_BuildHeatmap_output.rds")
  result <- BuildHeatmap(input)
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
      example_output_wide,
      d_range = "full"
    )
  )

  # wrong number of colors
  expect_error(
    BuildHeatmap(
      example_output_wide,
      d_col = c("red", "blue")
    ),
    "wrong number of colors in d_col!\nSupply colors for c(min, middle, max)!",
    fixed = TRUE
  )

  # wrong d_range
  expect_error(
    BuildHeatmap(
      example_output_wide,
      d_range = "partTime"
    ),
    'd_range must be either "fit" or "full"!',
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

test_that("Function `BuildHeatmap()` works correctly with partial Eff", {
  feature <- reduced_feature
  metaMat <- metaMatMetformin
  input <- readRDS("2024_10_10_example_output_partialRand.rds")

  expect_no_error(
    BuildHeatmap(
      input,
      plotPartial = "partial"
    )
  )

  expect_no_error(
    BuildHeatmap(
      input,
      plotPartial = "partialRel"
    )
  )

  expect_no_error(
    BuildHeatmap(
      input,
      plotPartial = "partialNorm"
    )
  )

})
