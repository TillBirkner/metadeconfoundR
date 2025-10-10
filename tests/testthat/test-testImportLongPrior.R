library(testthat)
library(metadeconfoundR)

test_that("standard options", {

  feature <- reduced_feature
  metaMat <- metaMatMetformin

  example_output <- read.table("2025_10_07_example_output.tsv", header = T, sep = "\t")

  minQValues <- ImportLongPrior(longPrior = example_output,
                                featureMat = feature,
                                metaMat = metaMat
                                )

  result <- MetaDeconfound(featureMat = feature,
                           metaMat = metaMat,
                           minQValues = minQValues,
                           returnLong = T
                           )
  result$feature <- as.character(result$feature)
  result$metaVariable <- as.character(result$metaVariable)
  expect_equal(result, example_output)

  # only names
  expect_no_error(
    ImportLongPrior(
      longPrior = example_output[, c("feature", "metaVariable")],
      featureMat = feature,
      metaMat = metaMat
    )
  )
  # no status labels
  expect_no_error(
    ImportLongPrior(
      longPrior = example_output[, c("feature", "metaVariable", "Qs")],
      featureMat = feature,
      metaMat = metaMat
    )
  )
  # no explicit Qs
  expect_no_error(
    ImportLongPrior(
      longPrior = example_output[, c("feature", "metaVariable", "status")],
      featureMat = feature,
      metaMat = metaMat
    )
  )
})
