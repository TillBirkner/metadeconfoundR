#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Build and Reload Package:  'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'

#' Naive Associations
#'
#' Naive Associations checks all feature <-> covariate combinations on significant associations with as few as possible assumptions about data structure.
#'
#' @param featureFile a tab delimited file or data frame with row(sample ID) and column(feature such as XXXXX) names listing features for all samples
#' @param metaFile a tab delimited file or data frame with row(sample ID) and column(meta data such as age,BMI and all possible confounders) names listing metadata for all samples. first column should be case status with case=1 and control=0.
#' @param nnodes number of nodes/cores to be used for parallel processing
#' @param cutoff minimamal number of sample size for each covariate in order to have sufficient power for association testing
#' @return XXXXdata frameXXX containing p-values and effect sizes for all simple feature/covariate associations
#' @importFrom foreach %dopar%
#' @export NaiveAssociation
#'

NaiveAssociation <- function(featureFile, metaFile, nnodes = 1, cutoff = 5, adjustMethod = "fdr") {

  featureMat <- featureFile # your input data here
  samples <- row.names (featureMat)
  features <- colnames (featureMat)
  noFeatures <- length (features)

  # read in matrix of metadata

  md <- metaFile # your input data here
  covariates <- colnames (md)# each covariate + the status category, specific to example
  noCovariates <- length (covariates)

  # define spaces to hold output of tests

  # Ps <- matrix (NA, noFeatures, noCovariates) # base P value
  # Ds <- matrix (NA, noFeatures, noCovariates) # effect size
  # Qs <- matrix (NA, noFeatures, noCovariates) # FDR
  # Ss <- matrix (NA, noFeatures, noCovariates) # status of post-hoc test
  #
  # colnames (Ps) <- covariates
  # colnames (Ds) <- covariates
  # colnames (Qs) <- covariates
  # colnames (Ss) <- covariates
  #
  # row.names (Ps) <- features
  # row.names (Ds) <- features
  # row.names (Qs) <- features
  # row.names (Ss) <- features

  # check sufficient power for all the covariates in metaFile

  print("NaiveAssociation  -  Computing Robustness")
  isRobust <- CheckSufficientPower(metaFile = md,
                                   nnodes = nnodes,
                                   cutoff = cutoff)


  # load parralel processing environment
  cl <- parallel::makeForkCluster(nnodes = nnodes, outfile = "")  # the parent process uses another core (so 4 cores will be used with this command)
  doParallel::registerDoParallel(cl)
  i <- 0


  write(paste
        ("counter",
          "aCovariate",
          "aFeature",
          "noCovariates" ,
          "isRobust[aCovariate, 4]",
          sep = "\t" ),
        file="progress.txt",
        append = FALSE)

  r = foreach::foreach(i= 1:5, .combine='rbind') %dopar% {

    somePs <- vector(length = noCovariates)
    someDs <- vector(length = noCovariates)

    #names(somePs) <- covariates
    #names(someDs) <- covariates



    for (j in 1:noCovariates) {

      aFeature <- as.character (features [i])
      aCovariate <- as.character (covariates [j])
      if (j>1) {
        write(paste
              (i,
                aCovariate,
                aFeature,
                noCovariates ,
                isRobust[aCovariate, 4],
                sep = "\t" ),
              file="progress.txt",
              append = TRUE)
      }

      aD <- NA_real_
      aP <- NA_real_

      if (!is.na(isRobust[aCovariate, 4]) && !isRobust[aCovariate, 4]) {
              somePs[j] <- aP
              someDs[j] <- aD
              next
      }

      subFeatures <- featureMat [,i]
      subMerge <- md
      subMerge$FeatureValue <- subFeatures # caveat - this works for this input, but in another case you may have to verify the samples have the same order...
      # results should be a data frame with FeatureValue and all predictors/covariates (columns) for each sample (rows)
      # now test very basic association without considering covariation, does aCovariate predict aFeature?

      if (length (unique (subMerge [, as.character (aCovariate)])) == 2 &&
          length (subMerge [subMerge [[as.character (aCovariate)]] == 0, "FeatureValue"]) > 1 &&
          length (subMerge [subMerge [[as.character (aCovariate)]] == 1, "FeatureValue"]) > 1) {

        # MWU test if binary
        aP <- stats::wilcox.test (subMerge [subMerge [[as.character (aCovariate)]] == 0, "FeatureValue"],
                                  subMerge [subMerge [[as.character (aCovariate)]] == 1, "FeatureValue"])$p.value
        aD <- orddom::orddom (subMerge [subMerge [[as.character (aCovariate)]] == 0, "FeatureValue"],
                              subMerge [subMerge [[as.character (aCovariate)]] == 1, "FeatureValue"]) [13]

      }

      else if (length (unique (subMerge [, as.character (aCovariate)])) > 2 &&
               length (subMerge [subMerge [[as.character (aCovariate)]] == 0, "FeatureValue"]) > 1 &&
               length (subMerge [subMerge [[as.character (aCovariate)]] == 1, "FeatureValue"]) > 1) {

        # spearman test if continuous
        aP <- stats::cor.test (subMerge [, as.character (aCovariate)],
                               subMerge [, "FeatureValue"],
                               method = "spearman")$p.value
        aD <- stats::cor.test (subMerge [, as.character (aCovariate)],
                               subMerge [, "FeatureValue"],
                               method = "spearman")$estimate
      }

      somePs[j] <- aP
      someDs[j] <- aD

    }

    return(c(as.numeric(somePs), as.numeric(someDs)))

  }

  parallel::stopCluster(cl)
  Ps <- r[, 1:(ncol(r)/2)]
  rownames(Ps) <- features[1:nrow(r)]
  colnames(Ps) <- covariates
  Ds <- r[, -(1:(ncol(r)/2))]
  rownames(Ds) <- features[1:nrow(r)]
  colnames(Ds) <- covariates

  #Qs <- matrix (NA, noFeatures, noCovariates)
  Qs <- matrix (NA, length(features[1:nrow(r)]), length(covariates))
  rownames(Qs) <- features[1:nrow(r)]
  colnames(Qs) <- covariates

  print(paste0
        ("NaiveAssociation  -  compute multiple testing p adjustment using ",
          adjustMethod))

  for (i in 1:ncol(Ps)) {
    Qs[, i] <- stats::p.adjust (Ps[, i], method = adjustMethod)
  }

  return(list(Ps=Ps, Ds=Ds, Qs=Qs))
}



