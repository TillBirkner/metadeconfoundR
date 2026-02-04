library(testthat)
library(metadeconfoundR)  # Your package name

test_that("standard options", {
  # read in input and expected output
  feature <- reduced_feature
  metaMat <- metaMatMetformin
  # using metadeconfound output created 2024 09 24
  #saveRDS(ex_out_partial, "tests/testthat/2024_10_09_example_output_partial.rds")
  expected_output <- readRDS("2024_10_09_example_output_partial.rds")
  #expected_output$feature <- as.factor(expected_output$feature)
  #expected_output$metaVariable <- as.factor(expected_output$metaVariable)

  # Call the function
  resultA <- readRDS("example_output.rds")
  load("2024_10_09_example_output_partial.RData")
  result <- GetPartialEfSizes(featureMat = feature,
                              metaMat = metaMat,
                              metaDeconfOutput = resultA)

  result$feature <- as.character(result$feature)
  result$metaVariable <- as.character(result$metaVariable)
  expected_output$feature <- as.character(expected_output$feature)
  expected_output$metaVariable <- as.character(expected_output$metaVariable)


  # Verify that the result matches the expected output
  testthat::expect_equal(result, expected_output)
})

test_that('with randomVar and fixedVar', {
  # read in input and expected output
  feature <- reduced_feature
  metaMat <- metaMatMetformin
  #saveRDS(result, "tests/testthat/2026_02_04_example_output_partialFixRand.rds")
  expected_output <- readRDS("2026_02_04_example_output_partialFixRand.rds")
  #expected_output$feature <- as.factor(expected_output$feature)
  #expected_output$metaVariable <- as.factor(expected_output$metaVariable)

  # Call the function
  resultA <- readRDS("2026_02_04_example_output_fixRand.rds")

  result <- GetPartialEfSizes(featureMat = feature,
                              metaMat = metaMat,
                              metaDeconfOutput = resultA,
                              randomVar = "Dataset",
                              fixedVar = "continuous_dummy")

  result$feature <- as.character(result$feature)
  result$metaVariable <- as.character(result$metaVariable)
  expected_output$feature <- as.character(expected_output$feature)
  expected_output$metaVariable <- as.character(expected_output$metaVariable)


  # Verify that the result matches the expected output
  testthat::expect_equal(result, expected_output, tolerance = 0.0001)
})
