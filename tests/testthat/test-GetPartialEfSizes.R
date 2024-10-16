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

test_that('with randomVar = "Dataset"', {
  # read in input and expected output
  feature <- reduced_feature
  metaMat <- metaMatMetformin
  #saveRDS(ex_out_partialRand, "tests/testthat/2024_10_10_example_output_partialRand.rds")
  expected_output <- readRDS("2024_10_10_example_output_partialRand.rds")
  #expected_output$feature <- as.factor(expected_output$feature)
  #expected_output$metaVariable <- as.factor(expected_output$metaVariable)

  # Call the function
  resultA <- readRDS("2024_10_09_example_output_rand.rds")

  result <- GetPartialEfSizes(featureMat = feature,
                              metaMat = metaMat,
                              metaDeconfOutput = resultA,
                              randomVar = "Dataset")

  result$feature <- as.character(result$feature)
  result$metaVariable <- as.character(result$metaVariable)
  expected_output$feature <- as.character(expected_output$feature)
  expected_output$metaVariable <- as.character(expected_output$metaVariable)


  # Verify that the result matches the expected output
  testthat::expect_equal(result, expected_output, tolerance = testthat_tolerance())
})
