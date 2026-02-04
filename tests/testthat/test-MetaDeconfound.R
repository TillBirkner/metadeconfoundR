library(testthat)
library(metadeconfoundR)

test_that("standard options", {
  feature <- reduced_feature
  metaMat <- metaMatMetformin
  # using metadeconfound output created 2025 10 07
  #write.table(result, "tests/testthat/2025_10_07_example_output.tsv", sep = "\t", row.names = F)
  expected_output <- read.table("2025_10_07_example_output.tsv", header = T, sep = "\t")

  result <- MetaDeconfound(featureMat = feature,
                           metaMat = metaMat,
                           logLevel = "INFO",
                           returnLong = T,
                           doRanks = c("altered_dummy") # test doRanks functionality
                           )

  result$feature <- as.character(result$feature)
  result$metaVariable <- as.character(result$metaVariable)
  expect_equal(result, expected_output)


  #single metadata input
  # saveRDS(resultOnlyStatus, "tests/testthat/2024_10_17_example_output_OnlyStatus.rds")
  expected_output_OnlyStatus <- readRDS("2024_10_17_example_output_OnlyStatus.rds")

  resultOnlyStatus <- MetaDeconfound(featureMat = feature,
                           metaMat = metaMat[, 1,drop = F],
                           logLevel = "INFO",
                           returnLong = T
  )
  resultOnlyStatus$feature <- as.character(resultOnlyStatus$feature)
  resultOnlyStatus$metaVariable <- as.character(resultOnlyStatus$metaVariable)
  expected_output_OnlyStatus$feature <- as.character(expected_output_OnlyStatus$feature)
  expected_output_OnlyStatus$metaVariable <- as.character(expected_output_OnlyStatus$metaVariable)
  expect_equal(resultOnlyStatus, expected_output_OnlyStatus)

  expect_no_error(
    resultOnlyStatus_partial <- GetPartialEfSizes(
      featureMat = feature,
      metaMat = metaMat[, 1, drop = F],
      metaDeconfOutput = resultOnlyStatus
    )
  )

  # saveRDS(result_wide, "tests/testthat/2025_10_07_example_output_wide.rds")
  expected_output_wide <- readRDS("2025_10_07_example_output_wide.rds")

  result_wide <- MetaDeconfound(featureMat = feature,
                           metaMat = metaMat,
                           logLevel = "ERROR",
                           returnLong = F
  )

  # result_wide$feature <- as.character(result_wide$feature)
  # result_wide$metaVariable <- as.character(result_wide$metaVariable)
  # expected_output_wide$feature <- as.character(expected_output_wide$feature)
  # expected_output_wide$metaVariable <- as.character(expected_output_wide$metaVariable)
  expect_equal(result_wide, expected_output_wide)


  result_naive <- MetaDeconfound(
    featureMat = feature,
    metaMat = metaMat,
    logLevel = "ERROR",
    returnLong = T,
    startStop = "naiveStop"
  )
  result_naive$feature <- as.character(result_naive$feature)
  result_naive$metaVariable <- as.character(result_naive$metaVariable)
  expect_equal(result_naive, expected_output[, 1:5])

  result_naive_wide <- MetaDeconfound(
    featureMat = feature,
    metaMat = metaMat,
    logLevel = "ERROR",
    returnLong = F,
    startStop = "naiveStop"
  )
  expect_equal(result_naive_wide, expected_output_wide[1:3])

  result_QValues <- MetaDeconfound(
    featureMat = feature,
    metaMat = metaMat,
    logLevel = "INFO",
    QValues = result_naive_wide$Qs,
    DValues = result_naive_wide$Ds,
    returnLong = T
  )
  result_QValues$feature <- as.character(result_QValues$feature)
  result_QValues$metaVariable <- as.character(result_QValues$metaVariable)
  #overwrite Ps, as they are simply a copy of Qs, and won't match expected_output like that
  result_QValues$Ps <- expected_output$Ps
  expect_equal(result_QValues, expected_output)

})


test_that("standard options parallel", {
  feature <- reduced_feature
  metaMat <- metaMatMetformin
  expected_output <- read.table("2025_10_07_example_output.tsv", header = T, sep = "\t")

  result <- MetaDeconfound(featureMat = feature,
                           metaMat = metaMat,
                           logLevel = "ERROR",
                           returnLong = T,
                           nnodes = 2
  )

  result$feature <- as.character(result$feature)
  result$metaVariable <- as.character(result$metaVariable)

  expect_equal(result, expected_output)


  # test adjust level 2 output
  expected_output <- readRDS("2025_10_07_example_output_adjustLevel2.rds")
  result <- MetaDeconfound(featureMat = feature,
                           metaMat = metaMat,
                           logLevel = "ERROR",
                           returnLong = T,
                           nnodes = 2,
                           adjustLevel = 2
  )
  #saveRDS(result, "tests/testthat/2025_10_07_example_output_adjustLevel2.rds")
  expect_equal(result, expected_output)

})

test_that("random and fixed effects", {
  feature <- reduced_feature
  metaMat <- metaMatMetformin
  # saveRDS(result, "tests/testthat/2024_10_16_example_output_rand.rds")
  expected_output <- readRDS("2024_10_16_example_output_rand.rds")

  temp_file <- tempfile()
  result <- MetaDeconfound(featureMat = feature,
                           metaMat = metaMat,
                           logLevel = "DEBUG",
                           returnLong = T,
                           doConfs = 1,
                           randomVar = c("Dataset"),
                           logfile = temp_file
  )
  unlink(temp_file)

  result$feature <- as.character(result$feature)
  result$metaVariable <- as.character(result$metaVariable)
  expected_output$feature <- as.character(expected_output$feature)
  expected_output$metaVariable <- as.character(expected_output$metaVariable)
  expect_equal(result, expected_output)


  # simply check, that collectMods is running without errors
  # saving all models to file would create >40 MB .rds file
  expect_no_error(MetaDeconfound(featureMat = feature[, 1:5],
                                 metaMat = metaMat,
                                 logLevel = "WARN",
                                 returnLong = T,
                                 collectMods = T,
                                 randomVar = c("Dataset")
  ))
  expect_no_error(MetaDeconfound(featureMat = feature[, 1:5],
                                 metaMat = metaMat,
                                 logLevel = "WARN",
                                 returnLong = F, # wide
                                 collectMods = T,
                                 randomVar = c("Dataset")
  ))

  # saveRDS(resultFix, "tests/testthat/2024_10_21_example_output_fix.rds")
  expected_output_fix <- readRDS("2024_10_21_example_output_fix.rds")
  resultFix <- MetaDeconfound(featureMat = feature,
                              metaMat = metaMat,
                              logLevel = "INFO",
                              returnLong = T,
                              fixedVar = c("continuous_dummy")
  )
  # 2024_10_18 TB: some small differences in optimazation calculations
  # within the model fitting steps of metadeconfoundR
  # (potentially caused by different BLAS/LAPACK implementations)
  # lead to differences in status label assignment in very rare cases.
  # the problematic case in this test data will be made equal manually
  # so that the tests run successfully on different OSs
  # print(lrtest(resultFix$collectedMods$MS0047$altered_dummy$Dataset$full,
  #              resultFix$collectedMods$MS0047$altered_dummy$Dataset$conf))
  # print(lrtest(resultFix$collectedMods$MS0047$altered_dummy$Dataset$full,
  #              resultFix$collectedMods$MS0047$altered_dummy$Dataset$cov))
  # print(lapply(resultFix$collectedMods$MS0047$altered_dummy$Dataset, summary))
  #
  # print(sessionInfo())

  #resultFix <- resultFix$stdOutput
  resultFix$feature <- as.character(resultFix$feature)
  resultFix$metaVariable <- as.character(resultFix$metaVariable)
  #resultFix$status[214] <- "OK_sd"
  #resultFix$status[204] <- "C: Dataset"
  expected_output_fix$feature <- as.character(expected_output_fix$feature)
  expected_output_fix$metaVariable <- as.character(expected_output_fix$metaVariable)
  #expect_equal(dim(resultFix), dim(expected_output_fix))
  #stop("Why are two status labels different between my machine and all unix github action runs????")
  expect_equal(resultFix, expected_output_fix, tolerance = 0.0001)

  #saveRDS(resultFixRand, "tests/testthat/2026_02_04_example_output_fixRand.rds")
  expected_output_fix_rand <- readRDS("2026_02_04_example_output_fixRand.rds")
  resultFixRand <- MetaDeconfound(featureMat = feature[, 1:15],
                                  metaMat = metaMat,
                                  logLevel = "INFO",
                                  returnLong = T,
                                  fixedVar = c("continuous_dummy"),
                                  randomVar = c("Dataset")
  )
  resultFixRand$feature <- as.character(resultFixRand$feature)
  resultFixRand$metaVariable <- as.character(resultFixRand$metaVariable)
  expected_output_fix_rand$feature <- as.character(expected_output_fix_rand$feature)
  expected_output_fix_rand$metaVariable <- as.character(expected_output_fix_rand$metaVariable)
  expect_equal(resultFixRand, expected_output_fix_rand)
})

test_that("logistic regression", {
  feature <- reduced_feature
  feature[feature > 0] <- 1
  metaMat <- metaMatMetformin
  # using metadeconfound output created 2024 09 26
  # write.table(example_output, "tests/testthat/2024_11_07_example_output_logistic.tsv", sep = "\t", row.names = F)
  # expected_output <- read.table("2024_11_07_example_output_logistic.tsv", header = T, sep = "\t")
  #
  # result <- MetaDeconfound(featureMat = feature,
  #                          metaMat = metaMat,
  #                          logLevel = "ERROR",
  #                          returnLong = T,
  #                          logistic = T
  # )
  #
  # result$feature <- as.character(result$feature)
  # result$metaVariable <- as.character(result$metaVariable)
  #
  # expect_equal(result, expected_output)

  expect_no_error(
    MetaDeconfound(featureMat = feature[, c("MS0035"), drop = F],
                   metaMat = metaMat,
                   logLevel = "WARN",
                   returnLong = T,
                   logistic = T,
                   collectMods = T
    )
  )

  log_file <- tempfile(fileext = ".txt")
  on.exit({
    if (file.exists(log_file)) unlink(log_file)
  }, add = TRUE)

  #saveRDS(result_loglog, "tests/testthat/2026_02_04_example_output_loglog.rds")
  expected_output_loglog <- readRDS("2026_02_04_example_output_loglog.rds")

  result_loglog <-  MetaDeconfound(featureMat = feature[, 5:10],
                 metaMat = metaMat,
                 logLevel = "WARN",
                 returnLong = T,
                 logistic = T,
                 logfile = log_file,
                 randomVar = "Dataset")
  # test output
  expect_equal(result_loglog, expected_output_loglog)

  # test logging output with correct number of separation warnings
  txt <- readLines(log_file, warn = FALSE)
  expect_equal(sum(
    grepl(
      "Separation for: MS0(035|037) and (Status|Dataset|Metformin)",
      txt
    )
  ), 4)

  # #combination of randomVars AND logistic mode
  # # saveRDS(resultlogRand, "tests/testthat/2026_01_07_example_output_log_rand.rds")
  # expected_output_logRand <- readRDS("2026_01_07_example_output_log_rand.rds")
  #
  # resultlogRand <- MetaDeconfound(featureMat = feature[, 1:15],
  #                          metaMat = metaMat,
  #                          logLevel = "INFO",
  #                          returnLong = T,
  #                          logistic = T,
  #                          randomVar = "Dataset"
  # )
  #
  # resultlogRand$feature <- as.character(resultlogRand$feature)
  # resultlogRand$metaVariable <- as.character(resultlogRand$metaVariable)
  # expected_output_logRand$feature <- as.character(expected_output_logRand$feature)
  # expected_output_logRand$metaVariable <- as.character(expected_output_logRand$metaVariable)
  #
  # expect_equal(resultlogRand, expected_output_logRand)

})

test_that("raw Counts mode", {
  feature <- round(reduced_feature)
  metaMat <- metaMatMetformin
  # saveRDS(result_rawCountsRand, "tests/testthat/2026_02_04_example_output_rawCountsRand.rds")
  expected_output_rawCountsRand <- readRDS("2026_02_04_example_output_rawCountsRand.rds")

  result_rawCountsRand <- MetaDeconfound(featureMat = feature[, 1:5],
                           metaMat = metaMat,
                           logLevel = "INFO",
                           returnLong = T,
                           rawCounts = T,
                           randomVar = "Dataset"
  )

  result_rawCountsRand$feature <- as.character(result_rawCountsRand$feature)
  result_rawCountsRand$metaVariable <- as.character(result_rawCountsRand$metaVariable)
  expected_output_rawCountsRand$feature <- as.character(expected_output_rawCountsRand$feature)
  expected_output_rawCountsRand$metaVariable <- as.character(expected_output_rawCountsRand$metaVariable)

  expect_equal(result_rawCountsRand, expected_output_rawCountsRand)
})


test_that("correct error handling", {
  feature <- reduced_feature
  metaMat <- metaMatMetformin

  # missing input: metaMat
  expect_error(
    MetaDeconfound(featureMat = feature,
                   logLevel = "ERROR",
                   returnLong = T
                   ),
    'Error - Necessary argument "metaMat" missing.',
    fixed = FALSE
  )

  # missing input: featureMat
  expect_error(
    MetaDeconfound(metaMat = metaMat,
                   logLevel = "ERROR",
                   returnLong = T
                   ),
    'Error - Necessary argument "featureMat" missing.',
    fixed = FALSE
  )

  # specific error
  expect_error(
    MetaDeconfound(featureMat = feature,
                   metaMat = metaMat[1:50, ],
                   logLevel = "ERROR",
                   returnLong = T
                   ),
    "featureMat and metaMat don't have same number of rows.",
    fixed = FALSE
  )

  # wrong order in metaMat/featureMat
  expect_error(
    MetaDeconfound(featureMat = feature,
                   metaMat = metaMat[nrow(metaMat):1, ],
                   logLevel = "ERROR",
                   returnLong = T
                   ),
    "Rownames of featureMat and metaMat don't have same order.
         (order(rownames(metaMat)) != order(rownames(featureMat)))",
    fixed = TRUE
  )

  # wrong metavar names in deconfT/deconfF
  expect_error(
    MetaDeconfound(featureMat = feature,
                   metaMat = metaMat,
                   logLevel = "INFO",
                   returnLong = T,
                   deconfT = c("notAMetaVariable")
                   ),
    "Elements of deconfT/deconfF are not present in colnames of metaMat.",
    fixed = FALSE
  )

  # metaVar in both deconfF AND deconfT
  expect_error(
    MetaDeconfound(featureMat = feature,
                   metaMat = metaMat,
                   logLevel = "ERROR",
                   returnLong = T,
                   deconfT = c("Dataset"),
                   deconfF = c("Dataset")
    ),
    "Some elements of deconfT and deconfF seem to be identical.",
    fixed = FALSE
  )

  # rawCounts AND logistic both set to true
  expect_error(
    MetaDeconfound(featureMat = feature,
                   metaMat = metaMat,
                   logLevel = "ERROR",
                   returnLong = T,
                   rawCounts = T,
                   logistic = T
    ),
    "rawCounts and logistic can not be both set to TRUE!",
    fixed = FALSE
  )

  # empty QValues argument
  expect_error(
    MetaDeconfound(featureMat = feature,
                   metaMat = metaMat,
                   logLevel = "ERROR",
                   returnLong = T,
                   QValues = NULL
    ),
    "QValues and/or DValues argument is supplied but seems to be empty (NULL).",
    fixed = TRUE
  )

  # old randomVar
  expect_error(
    MetaDeconfound(featureMat = feature,
                   metaMat = metaMat,
                   logLevel = "ERROR",
                   returnLong = T,
                   randomVar = list(a = "test", b = "expectingError")
    ),
    "randomVar does not need to be supplied as list anymore, please change to new syntax.",
    fixed = FALSE
  )

  # below robustCutoff
  expect_error(
    MetaDeconfound(featureMat = feature[1:11, ],
                   metaMat = metaMat[1:11, ],
                   logLevel = "ERROR",
                   returnLong = T
    ),
    'Not enough(robustCutoff = 5 ) samples in either case or controle group.',
    fixed = TRUE
  )

  # 2025 10 07
  # non-numeric columns in featureMat
  featureWithNonNum <- feature
  featureWithNonNum$nonNumCol <- rownames(featureWithNonNum)
  expect_error(
    MetaDeconfound(featureMat = featureWithNonNum,
                   metaMat = metaMat,
                   logLevel = "ERROR",
                   returnLong = T
    ),
    'Non-numeric columns detected in featureMat.',
    fixed = TRUE
  )

  # # specific error
  # expect_error(
  #   MetaDeconfound(featureMat = feature,
  #                  metaMat = metaMat,
  #                  logLevel = "ERROR",
  #                  returnLong = T
  #   ),
  #   "expected error message",
  #   fixed = FALSE
  # )
})
