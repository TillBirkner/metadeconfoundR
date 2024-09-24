library(testthat)
library(metadeconfoundR)  # Your package name

test_that("standard options", {
  # read in input and expected output
  feature <- reduced_feature
  metaMat <- metaMatMetformin
  # using metadeconfound output created 2024 09 24
  # write.table(example_output, "tests/testthat/example_output.tsv", sep = "\t", row.names = F)
  expected_output <- read.table("example_output.tsv", header = T, sep = "\t")
  #expected_output$feature <- as.factor(expected_output$feature)
  #expected_output$metaVariable <- as.factor(expected_output$metaVariable)

  # Call the function
  result <- MetaDeconfound(featureMat = feature,
                           metaMat = metaMat,
                           logLevel = "ERROR",
                           returnLong = T
                           )

  result$feature <- as.character(result$feature)
  result$metaVariable <- as.character(result$metaVariable)

  # Verify that the result matches the expected output
  expect_equal(result, expected_output)
})

test_that("random effects", {
  # read in input and expected output
  feature <- reduced_feature
  metaMat <- metaMatMetformin
  # using metadeconfound output created 2024 09 24
  # write.table(example_output_rand, "tests/testthat/example_outputrand.tsv", sep = "\t", row.names = F)
  expected_output <- read.table("example_output_rand.tsv", header = T, sep = "\t")
  #expected_output$feature <- as.factor(expected_output$feature)
  #expected_output$metaVariable <- as.factor(expected_output$metaVariable)

  # Call the function
  result <- MetaDeconfound(featureMat = feature,
                           metaMat = metaMat,
                           logLevel = "ERROR",
                           returnLong = T,
                           randomVar = c("Dataset")
  )

  result$feature <- as.character(result$feature)
  result$metaVariable <- as.character(result$metaVariable)

  # Verify that the result matches the expected output
  expect_equal(result, expected_output)
})
