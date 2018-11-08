#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Build and Reload Package:  'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'


#' Check sufficiant power of all possible covariate groups
#'
#' @param featureFile a tab delimited file or data frame with row and colum names listing features for all samples
#' @param metaFile a tab delimited file or data frame with row and colum names listing metadata for all samples
#' @param nnodes number of nodes to be used for parallel computing
#' @return data frame containing p-values and effect sizes for all simple feature/covariate associations
#' @importFrom foreach %dopar%
#' @export NaiveAssociation
#'

NaiveAssociation <- function(featureFile, metaFile, nnodes=5) {

  featureMat <- featureFile # your input data here
  samples <- row.names (featureMat)
  features <- colnames (featureMat)
  noFeatures <- length (features)

  # read in matrix of metadata

  md <- metaFile # your input data here
  md$Case_status <- as.numeric (md$Case_status) - 1 # convert so that both drugs and other features has same scale, for simplicity, specific to example
  covariates <- colnames (md) [4:27] # each covariate + the status category, specific to example
  noCovariates <- length (covariates)

  # define spaces to hold output of tests

  Ps <- matrix (NA, noFeatures, noCovariates) # base P value
  Ds <- matrix (NA, noFeatures, noCovariates) # effect size
  Qs <- matrix (NA, noFeatures, noCovariates) # FDR
  Ss <- matrix (NA, noFeatures, noCovariates) # status of post-hoc test

  colnames (Ps) <- covariates
  colnames (Ds) <- covariates
  colnames (Qs) <- covariates
  colnames (Ss) <- covariates

  row.names (Ps) <- features
  row.names (Ds) <- features
  row.names (Qs) <- features
  row.names (Ss) <- features

  old <- Sys.time()
  # load parralel processing environment
  #cl<-makeCluster(3)
  cl <- parallel::makeForkCluster(nnodes=nnodes, outfile="")  # the parent process uses another core (so 4 cores will be used with this command)
  doParallel::registerDoParallel(cl)
  i <- 0
  r = foreach::foreach(i= 1:5, .combine = 'rbind') %dopar% {
    #print (i) # status marker to check it has not stalled or something
    somePs <-vector(length=noCovariates)
    someDs <-vector(length=noCovariates)
    print(i)
    for (j in 1:noCovariates) {

      aFeature <- as.character (features [i])
      aCovariate <- as.character (covariates [j])
      #rprintf(x = '%d Features processed\n',i)
      if (j==1) {
        write(paste(i, aCovariate, aFeature, noCovariates , sep = "\t" ), file="progress.txt", append = TRUE)
      }

      subFeatures <- featureMat [,i]
      subMerge <- md
      subMerge$FeatureValue <- subFeatures # caveat - this works for this input, but in another case you may have to verify the samples have the same order...
      # results should be a data frame with FeatureValue and all predictors/covariates (columns) for each sample (rows)
      # now test very basic association without considering covariation, does aCovariate predict aFeature?

      aD <- NA
      aP <- NA

      if (length (unique (subMerge [, as.character (aCovariate)])) == 2 && length (subMerge [subMerge [[as.character (aCovariate)]] == 0, "FeatureValue"]) > 1 && length (subMerge [subMerge [[as.character (aCovariate)]] == 1, "FeatureValue"]) > 1) {

        # MWU test if binary
        aP <- stats::wilcox.test (subMerge [subMerge [[as.character (aCovariate)]] == 0, "FeatureValue"], subMerge [subMerge [[as.character (aCovariate)]] == 1, "FeatureValue"])$p.value
        aD <- orddom::orddom (subMerge [subMerge [[as.character (aCovariate)]] == 0, "FeatureValue"], subMerge [subMerge [[as.character (aCovariate)]] == 1, "FeatureValue"]) [13]

      }

      else if (length (unique (subMerge [, as.character (aCovariate)])) > 2 && length (subMerge [subMerge [[as.character (aCovariate)]] == 0, "FeatureValue"]) > 1 && length (subMerge [subMerge [[as.character (aCovariate)]] == 1, "FeatureValue"]) > 1) {

        # spearman test if continuous
        aP <- stats::cor.test (subMerge [, as.character (aCovariate)], subMerge [, "FeatureValue"], method = "spearman")$p.value
        aD <- stats::cor.test (subMerge [, as.character (aCovariate)], subMerge [, "FeatureValue"], method = "spearman")$estimate
      }

      #		Ps [as.character (aFeature), as.character (aCovariate)] <- aP
      #		Ds [as.character (aFeature), as.character (aCovariate)] <- aD
      Ps[i,j] <- aP
      Ds[i,j] <- aD
      somePs[j] <- aP
      someDs[j] <- aD

    }

    return(list(somePs,someDs))

  }

  parallel::stopCluster(cl)
  return(r)
}



