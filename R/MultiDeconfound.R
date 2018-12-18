#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Build and Reload Package:  'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'

#' Multi Deconfound
#'
#' Multi Deconfound checks all feature <-> covariate combinations for
#' counfounding effects of covariates on feature <-> effect corellation
#'
#' @param featureMat a tab delimited file or data frame with row(sample ID) and
#' column(feature such as metabolite or gut microbial OTU )
#' names, listing features for all samples
#' @param metaMat a tab delimited file or data frame with row(sample ID) and
#' column(meta data such as age,BMI and all possible confounders)
#' names listing metadata for all samples. first column should be case status
#' with case=1 and control=0.
#' @param nnodes number of nodes/cores to be used for parallel processing
#' @param adjustMethod multiple testing p-value correction using one of the
#' methods of stats::p.adjust.methods
#' @param robustCutoff minimamal number of sample size for each covariate
#' in order to have sufficient power for association testing
#' @param QCutoff significance cutoff for q-value, DEFAULT = 0.1
#' @param DCutoff effect size cutoff
#' (either cliff's delta or spearman correlation test estimate), DEFAULT = 0
#' @param ... for additional arguments used internally (development/debugging)
#' @return list with elements Ds = effectsize,
#' Ps = uncorrected p-value for naive association, Qs = corrected p-value/fdr,
#' and status=  confounding/mediation status for all
#' feature <=> covariate combinations
#' @export MultiDeconfound
#'


MultiDeconfound <- function(featureMat,
                             metaMat,
                             nnodes = 1,
                             adjustMethod = "fdr",
                             robustCutoff = 5,
                             QCutoff = 0.1,
                             DCutoff = 0,
                             ...) {

  if (nrow(metaMat) != nrow(featureMat)) {
    stop("featureMat and metaMat don't have same number of rows.")
  }

  .MultiDeconfound(featureMat = featureMat,
                   metaMat = metaMat,
                   nnodes = nnodes,
                   adjustMethod = adjustMethod,
                   robustCutoff = robustCutoff,
                   QCutoff = QCutoff,
                   DCutoff = DCutoff,
                   ...
                   #maintenance = maintenance,
                   #verbosity = verbosity
                   )
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
  covariates <- colnames (metaMat)
    # each covariate + the status category, specific to example
  noCovariates <- length (covariates)

  isRobust <- CheckSufficientPower(metaMat = metaMat,
                                   nnodes = nnodes,
                                   robustCutoff = robustCutoff,
                                   maintenance = maintenance,
                                   verbosity = verbosity)

  if (verbosity == "debug") {
    print(isRobust[1:10, ])
  }


  naiveAssociation <- NaiveAssociation(featureMat = featureMat,
                                       metaMat = metaMat,
                                       isRobust = isRobust,
                                       adjustMethod = adjustMethod,
                                       nnodes = nnodes,
                                       maintenance = maintenance,
                                       verbosity = verbosity)
  if (verbosity == "debug") {
    print(naiveAssociation$Qs[1:5, 1:5])
    print(naiveAssociation$Ds[1:5, 1:5])
    print("now computing confounding status")
  }


  reducibilityStatus <- CheckReducibility(featureMat = featureMat,
                                          metaMat = metaMat,
                                          noFeatures = noFeatures,
                                          noCovariates = noCovariates,
                                          features = features,
                                          covariates = covariates,
                                          Qs = naiveAssociation$Qs,
                                          Ds = naiveAssociation$Ds,
                                          nnodes = nnodes,
                                          QCutoff = QCutoff,
                                          DCutoff = DCutoff,
                                          verbosity = verbosity)

  if (verbosity == "debug") {
    print("MultiDeconfound  --  All done!")
  }
  return(list(Ps = naiveAssociation$Ps,
              Qs = naiveAssociation$Qs,
              Ds = naiveAssociation$Ds,
              status=reducibilityStatus))
}
