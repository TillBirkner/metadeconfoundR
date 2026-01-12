library(testthat)
library(metadeconfoundR)
library(ggraph)

test_that("Function `BuildConfounderMap()` works correctly", {
  feature <- reduced_feature
  metaMat <- metaMatMetformin

  input <- read.table("2025_10_07_example_output.tsv", header = T, sep = "\t")
  input$feature <- as.factor(input$feature)
  input$metaVariable <- as.factor(input$metaVariable)
  expect_no_error(BuildConfounderMap(input))


  input_short <- input[input$feature %in% colnames(feature)[1:5], ]
  feat_short <- feature[, 1:5]

  expect_warning(
    BuildConfounderMap(
      input_short,
      metaVariableNames = as.data.frame(cbind(colnames(metaMat), paste0(colnames(metaMat), "_nice"))),
      featureNames = cbind(colnames(feat_short), paste0(colnames(feat_short), "_nice"))
      ),
    'class(featureNames) was coerced to "data.frame"',
    fixed = TRUE
    )

  expect_warning(
    BuildConfounderMap(
      input_short,
      metaVariableNames = cbind(colnames(metaMat), paste0(colnames(metaMat), "_nice")),
      featureNames = as.data.frame(cbind(colnames(feat_short), paste0(colnames(feat_short), "_nice")))
    ),
    'class(metaVariableNames) was coerced to "data.frame"',
    fixed = TRUE
  )

  expect_warning(
    BuildConfounderMap(
      input_short,
      metaVariableNames = as.data.frame(cbind(colnames(metaMat), paste0(colnames(metaMat), "_nice"))),
      featureNames = as.data.frame(cbind(colnames(feat_short), rep("niceNames", ncol(feat_short))))
    ),
    'non-unique human-readable feature names where made unique using base::make.unique',
    fixed = TRUE
  )

  expect_warning(
    BuildConfounderMap(
      input_short,
      metaVariableNames = as.data.frame(cbind(colnames(metaMat), rep("niceNames", ncol(metaMat)))),
      featureNames = as.data.frame(cbind(colnames(feat_short), paste0(colnames(feat_short), "_nice")))
    ),
    c('non-unique human-readable metaVariable names where made unique using base::make.unique'),
    fixed = TRUE
  )

  # wrong number of colors
  expect_error(
    BuildConfounderMap(
      input,
      d_col = c("red", "blue")
    ),
    "wrong number of colors in d_col!\nSupply colors for c(min, middle, max)!",
    fixed = TRUE
  )


  # no trusted
  expect_error(
    BuildConfounderMap(
      input,
      trusted = c()
    ),
    '"trusted" must contain at least one trusted status label',
    fixed = TRUE
  )


  featureColor <- sample(
    c("darkred", "red", "orange", "green", "cyan", "blue", "purple", "pink"),
    50,
    replace = T)

  # no named elements in featureColor vector
  expect_warning(
    BuildConfounderMap(
      input_short,
      featureColor = featureColor[1:5]
    ),
    'No names assigned to featureColor. Unique feature names are selected instead.',
    fixed = TRUE
  )

  # length of featureColor not correct
  expect_error(
    BuildConfounderMap(
      input,
      featureColor = featureColor[1:10]
    ),
    'featureColor should either have length == 1 or the same length as unique features available in metaDeconfOutput!',
    fixed = TRUE
  )


  # test importing wide-format data
  example_output_wide <- readRDS("2025_10_07_example_output_wide.rds")
  expect_no_error(
    BuildConfounderMap(
      example_output_wide
    )
  )

})

