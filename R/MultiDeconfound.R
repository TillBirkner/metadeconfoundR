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
#' @param featureMat a tab delimited file or data frame with row(sample ID)
#' and column(feature such as metabolite or microbial OTU )
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
#' @param PHS_cutoff PostHoc Significance cutoff
#' @param ... for additional arguments used internally (development/debugging)
#' @return list with elements Ds = effectsize,
#' Ps = uncorrected p-value for naive association,
#' Qs = multiple testing corrected p-value/fdr,
#' and status = confounding/mediation status for all
#' feature <=> covariate combinations
#' @examples
#'data(reduced_feature)
#'data(metaMatMetformin)
#'example_output <- MultiDeconfound(featureMat = reduced_feature,
#'                                   metaMat = metaMatMetformin)
#'
#' @export
#'


MultiDeconfound <- function(featureMat,
                             metaMat,
                             nnodes = 1,
                             adjustMethod = "fdr",
                             robustCutoff = 5,
                             QCutoff = 0.1,
                             DCutoff = 0,
                             PHS_cutoff = 0.05,
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
                   PHS_cutoff = PHS_cutoff,
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
                            PHS_cutoff = 0.05,
                            maintenance = FALSE,
                            verbosity = "silent") {


  samples <- row.names (featureMat)
  features <- colnames (featureMat)
  noFeatures <- length (features)
  covariates <- colnames (metaMat) # each covariate + the status category
  noCovariates <- length (covariates)

  if (nnodes > 1) {
    nnodes <- nnodes - 1
  }

  isRobust <- CheckSufficientPower(metaMat = metaMat,
                                   covariates = covariates,
                                   noCovariates = noCovariates,
                                   nnodes = nnodes,
                                   robustCutoff = robustCutoff,
                                   maintenance = maintenance,
                                   verbosity = verbosity)

  if (verbosity == "debug") {
    print("CheckSufficientPower -- head(isRobust):")
    print(utils::head(isRobust))
  }


  naiveAssociation <- NaiveAssociation(featureMat = featureMat,
                                       samples = samples,
                                       features = features,
                                       noFeatures = noFeatures,
                                       metaMat = metaMat,
                                       covariates = covariates,
                                       noCovariates = noCovariates,
                                       isRobust = isRobust,
                                       adjustMethod = adjustMethod,
                                       nnodes = nnodes,
                                       maintenance = maintenance,
                                       verbosity = verbosity)
  if (verbosity == "debug") {
    print(naiveAssociation$Ps[seq_len(5), seq_len(3)])
    print(naiveAssociation$Qs[seq_len(5), seq_len(3)])
    print(naiveAssociation$Ds[seq_len(5), seq_len(3)])
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
                                          PHS_cutoff = PHS_cutoff,
                                          verbosity = verbosity)

  if (verbosity == "debug") {
    print("MultiDeconfound  --  All done!")
  }
  return(list(Ps = naiveAssociation$Ps,
              Qs = naiveAssociation$Qs,
              Ds = naiveAssociation$Ds,
              status=reducibilityStatus))
}
