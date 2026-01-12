library(testthat)
library(metadeconfoundR)

test_that("Function `BuildHeatmap()` works correctly", {
  feature <- reduced_feature
  metaMat <- metaMatMetformin

  input <- read.table("2025_10_07_example_output.tsv", header = T, sep = "\t")
  input$feature <- as.factor(input$feature)
  input$metaVariable <- as.factor(input$metaVariable)
  # saveRDS(result, "tests/testthat/20241128_exampl_output_BuildHeatmap.rds")
  expected_output <- readRDS("20241128_exampl_output_BuildHeatmap.rds")
  result <- BuildHeatmap(input)
  expect_equal(result$data, expected_output$data)

  # "nice names" format checks
  expect_warning(
    BuildHeatmap(
      input,
      cuneiform = T,
      keepMeta = colnames(metaMat),
      keepFeature = colnames(input),
      metaVariableNames = as.data.frame(cbind(colnames(metaMat), paste0(colnames(metaMat), "_nice"))),
      featureNames = cbind(colnames(feature), paste0(colnames(feature), "_nice"))
    ),
    'class(featureNames) was coerced to "data.frame"',
    fixed = TRUE
  )
  expect_warning(
    BuildHeatmap(
      input,
      cuneiform = T,
      keepMeta = colnames(metaMat),
      keepFeature = colnames(input),
      metaVariableNames = cbind(colnames(metaMat), paste0(colnames(metaMat), "_nice")),
      featureNames = as.data.frame(cbind(colnames(feature), paste0(colnames(feature), "_nice")))
    ),
    'class(metaVariableNames) was coerced to "data.frame"',
    fixed = TRUE
  )
  expect_warning(
    BuildHeatmap(
      input,
      cuneiform = T,
      keepMeta = colnames(metaMat),
      keepFeature = colnames(input),
      metaVariableNames = as.data.frame(cbind(colnames(metaMat), paste0(colnames(metaMat), "_nice"))),
      featureNames = as.data.frame(cbind(colnames(feature), rep("niceNames", ncol(feature))))
    ),
    'non-unique human-readable feature names where made unique using base::make.unique',
    fixed = TRUE
  )
  expect_warning(
    BuildHeatmap(
      input,
      cuneiform = T,
      keepMeta = colnames(metaMat),
      keepFeature = colnames(input),
      metaVariableNames = as.data.frame(cbind(colnames(metaMat), rep("niceNames", ncol(metaMat)))),
      featureNames = as.data.frame(cbind(colnames(feature), paste0(colnames(feature), "_nice")))
    ),
    c('non-unique human-readable metaVariable names where made unique using base::make.unique'),
    fixed = TRUE
  )

  expect_error(
    BuildHeatmap(
      input,
      d_cutoff = 0.8
    ),
    'No associations pass current q_cutoff and/or d_cutoff filters!',
    fixed = TRUE
  )

  # test for weird settings, where only confounded signal remains in plot
  expect_error(
    BuildHeatmap(
      input[input$feature == "MS0006" & input$metaVariable != "Status", ],
      intermedData = F,
    ),
    "No unconfounded associations remain with the current cutoff values. Consider manually including categorical metaVariables into the Heatmap by listing them through the 'keepMeta' argument.",
    fixed = TRUE
  )
  expect_warning(
    BuildHeatmap(
      input[input$feature == "MS0006" & input$metaVariable != "Status", ],
      intermedData = T,
    ),
    "No unconfounded associations remain with the current cutoff values. ",
    fixed = TRUE
  )

  # test legend setup for plots without any confounded signal
  expect_no_error(
    BuildHeatmap(
      MetaDeconfound(feature,
                     metaMat[, c("Status"), drop = F],
                     logLevel = "ERROR"
                     )
      )
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

  example_output_wide <- readRDS("2025_10_07_example_output_wide.rds")
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
