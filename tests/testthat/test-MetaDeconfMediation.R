library(testthat)
library(metadeconfoundR)

test_that("standard options", {
  feature <- reduced_feature[, 1:40]
  mediationMat <- reduced_feature[, 41:50]
  metaMat <- metaMatMetformin

  # using metadeconfound output created 2024 09 24
  # saveRDS(result, "tests/testthat/2024_11_18_example_output_mediation.rds")
  expected_output <- readRDS("2024_11_18_example_output_mediation.rds")

  result <- MetaDeconfound(featureMat = feature,
                           metaMat = metaMat,
                           mediationMat = mediationMat,
                           logLevel = "WARN",
                           returnLong = T
                           )

  result$feature <- as.character(result$feature)
  result$metaVariable <- as.character(result$metaVariable)
  expected_output$feature <- as.character(expected_output$feature)
  expected_output$metaVariable <- as.character(expected_output$metaVariable)
  expect_equal(result, expected_output)
})
