#' Multi Deconfound
#'
#' Multi Deconfound checks all feature <-> covariate combinations for
#' counfounding effects of covariates on feature <-> effect corellation
#'
#' @param featureMat a data frame with row(sample ID)
#' and column(feature such as metabolite or microbial OTU )
#' names, listing features for all samples
#' @param metaMat a data frame with row(sample ID) and
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
#' @param NA_imputation Missing data treatment. "remove": remove NA containing
#' data from analysis (default).
#' "others": Impute NAs  using methods from packages MICE or AMELIA
#' @param ... for additional arguments used internally (development/debugging)
#' @return list with elements Ds = effectsize,
#' Ps = uncorrected p-value for naive association,
#' Qs = multiple testing corrected p-value/fdr,
#' and status = confounding/mediation status for all
#' feature <=> covariate combinations with following categories:
#' (NS = not significant, SD = strictly deconfounded, LD = laxly deconfounded,
#' NC = no covariates, "covariate name" = confounded by this covariate)
#' @details for more details and explanations please see the vignette.
#' @examples
#'data(reduced_feature)
#'data(metaMatMetformin)
#'\donttest{
#'example_output <- MultiDeconfound(featureMat = reduced_feature,
#'                                   metaMat = metaMatMetformin)}
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
                             NA_imputation = "remove",
                             ...) {
  if (missing(metaMat)) {
    stop('Error - Necessary argument "metaMat" missing.')
  }
  if (missing(featureMat)) {
    stop('Error - Necessary argument "featureMat" missing.')
  }
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
                            NA_imputation = "remove",
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

  if (NA_imputation != "remove") {
    #impute using mice or somn like dat
    #featureMat <- mice(featuremat)
    #metaMat <- mice(metaMat)
  }

  isRobust <- CheckSufficientPower(metaMat = metaMat,
                                   covariates = covariates,
                                   noCovariates = noCovariates,
                                   nnodes = nnodes,
                                   robustCutoff = robustCutoff,
                                   NA_imputation = NA_imputation,
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
    print(naiveAssociation$Ps[seq_len(3), seq_len(noCovariates)])
    print(naiveAssociation$Qs[seq_len(3), seq_len(noCovariates)])
    print(naiveAssociation$Ds[seq_len(3), seq_len(noCovariates)])
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
