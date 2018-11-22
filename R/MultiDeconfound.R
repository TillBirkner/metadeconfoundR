#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Build and Reload Package:  'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'

#' Multi Deconfound
#'
#' Multi Deconfound checks all feature <-> covariate combinations for counfounding effects of covariates on feature <-> effect corellation
#'
#' @param featureMat a tab delimited file or data frame with row(sample ID) and column(feature such as XXXXX) names listing features for all samples
#' @param metaMat a tab delimited file or data frame with row(sample ID) and column(meta data such as age,BMI and all possible confounders) names listing metadata for all samples. first column should be case status with case=1 and control=0.
#' @param nnodes number of nodes/cores to be used for parallel processing
#' @param robustCutoff minimamal number of sample size for each covariate in order to have sufficient power for association testing
#' @return XXXXdata frameXXX containing information on reducability/confounding
#' @export MultiDeconfound
#'


MultiDeconfound <- function(featureMat,
                             metaMat,
                             nnodes = 1,
                             robustCutoff = 5,
                             adjustMethod = "fdr",
                             QCutoff = 0.1,
                             DCutoff = 0,
                             ...) {
  .MultiDeconfound(featureMat,
                   metaMat,
                   nnodes,
                   robustCutoff,
                   adjustMethod,
                   QCutoff,
                   DCutoff,
                   ...)
}

.MultiDeconfound <- function(featureMat,
                            metaMat,
                            nnodes = 1,
                            robustCutoff = 5,
                            adjustMethod = "fdr",
                            QCutoff = 0.1,
                            DCutoff = 0,
                            maintenance = FALSE,
                            verbosity = "silent") {

  #featureMat <- featureFile # your input data here
  samples <- row.names (featureMat)
  features <- colnames (featureMat)
  noFeatures <- length (features)

  # read in matrix of metadata

  #md <- metaFile # your input data here
  covariates <- colnames (metaMat)# each covariate + the status category, specific to example
  noCovariates <- length (covariates)

  isRobust <- CheckSufficientPower(metaFile = metaMat,
                                   nnodes = nnodes,
                                   cutoff = robustCutoff)

  if (verbosity == "debug") {
    print(isRobust[1:10, ])
  }


  naiveAssociation <- NaiveAssociation(featureFile = featureMat,
                                       metaFile = metaMat,
                                       isRobust = isRobust,
                                       adjustMethod = adjustMethod,
                                       nnodes = nnodes,
                                       maintenance = maintenance)
  if (verbosity == "debug") {
    print(naiveAssociation$Qs[1:5, 1:5])
    print(naiveAssociation$Ds[1:5, 1:5])
    print("now computing confounding status")
  }


  reducabilityStatus <- CheckReducibility(noFeatures,
                                          noCovariates,
                                          features,
                                          covariates,
                                          Qs = naiveAssociation$Qs,
                                          Ds = naiveAssociation$Ds,
                                          nnodes = nnodes,
                                          QCutoff = QCutoff,
                                          DCutoff = DCutoff,
                                          verbosity = verbosity)

  if (verbosity == "debug") {
    print("MultiDeconfound  --  All done!")
  }
  return(list(Ps = naiveAssociation$Ps, Qs = naiveAssociation$Qs, Ds = naiveAssociation$Ds, status=reducabilityStatus))
}
