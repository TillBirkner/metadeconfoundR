library(testthat)
library(metadeconfoundR)

test_that("standard options", {
  feature <- reduced_feature[, 1:40]
  mediationMat <- reduced_feature[, 41:50]
  metaMat <- metaMatMetformin

  # using metadeconfound output created 2024 09 24
  # saveRDS(result, "tests/testthat/2025_10_07_example_output_mediation.rds")
  expected_output <- readRDS("2025_10_07_example_output_mediation.rds")

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

  # wrong adjustLevel
  expect_error(
    MetaDeconfound(featureMat = feature,
                   metaMat = metaMat,
                   adjustLevel = 3,
                   logLevel = "ERROR",
                   returnLong = T
    ),
    'adjustLevel == 3 not possible without supplying mediationMat.',
    fixed = TRUE
  )

  # wrong dimensions of mediationMat
  expect_error(
    MetaDeconfound(featureMat = feature,
                   metaMat = metaMat,
                   mediationMat = mediationMat[-750, ],
                   logLevel = "ERROR",
                   returnLong = T
    ),
    "mediationMat and metaMat don't have same number of rows.",
    fixed = TRUE
  )

  # wrong row order of mediationMat
  expect_error(
    MetaDeconfound(featureMat = feature,
                   metaMat = metaMat,
                   mediationMat = mediationMat[sample(c(1:nrow(mediationMat)), nrow(mediationMat)), ],
                   logLevel = "ERROR",
                   returnLong = T
    ),
    "Rownames of mediationMat and metaMat don't have same order.",
    fixed = TRUE
  )

})
