#' @importFrom foreach %dopar%
# @export NaiveAssociation
#

NaiveAssociation <- function(featureMat,
                             metaMat,
                             isRobust,
                             nnodes,
                             maintenance,
                             adjustMethod,
                             verbosity) {


  if (maintenance == TRUE) {
    try(Qs <- utils::read.table(
      "~/Dropbox/UNI/Forslund_project/schizo/output_provided/FDR.r",
      header = TRUE,
      sep = "\t",
      row.names = 1))
    try(Ds <- utils::read.table(
      "~/Dropbox/UNI/Forslund_project/schizo/output_provided/D.r",
      header = TRUE,
      sep = "\t",
      row.names = 1))
    Ps <- "dummy"
    return(list(Ps=Ps, Ds=Ds, Qs=Qs))
  }
  featureMat <- featureMat # your input data here
  samples <- row.names (featureMat)
  features <- colnames (featureMat)
  noFeatures <- length (features)

  # read in matrix of metadata

  md <- metaMat
    # your input data here
  covariates <- colnames (md)
    # each covariate + the status category, specific to example
  noCovariates <- length (covariates)


  # load parralel processing environment
  cl <- parallel::makeForkCluster(nnodes = nnodes, outfile = "")
  # the parent process uses another core (nnodes - 1 might be good)
  doParallel::registerDoParallel(cl)
  i <- 0


  if (verbosity== "debug") {
    write(paste
        ("counter",
          "aCovariate",
          "aFeature",
          "noCovariates" ,
          "isRobust[aCovariate, 4]",
          sep = "\t" ),
        file="progress.txt",
        append = FALSE)
  }

  r = foreach::foreach(i= 1:noFeatures, .combine='rbind') %dopar% {

    somePs <- vector(length = noCovariates)
    someDs <- vector(length = noCovariates)

    #names(somePs) <- covariates
    #names(someDs) <- covariates



    for (j in 1:noCovariates) {

      aFeature <- as.character (features [i])
      aCovariate <- as.character (covariates [j])
      if (verbosity == "debug" && j>1) {
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
              if (verbosity == "debug") {
                write("skipped",
                      file="progress.txt",
                      append = TRUE)
              }

              next
      }

      subFeatures <- featureMat [,i]
      subMerge <- md
      subMerge$FeatureValue <- subFeatures # caveat *1*: see snippets.R

      con1 <- length (unique (subMerge [, aCovariate])) == 2
      # the distribution of the covariate is binary
      con2 <- length (unique (subMerge [, aCovariate])) > 2
      # the distribution of the covariate is continuous
      con3 <- length (subMerge[subMerge[[aCovariate]] == 0, "FeatureValue"]) > 1
      # feature has a measurement in more than one sample with
      #covaraite status 0
      con4 <- length (subMerge[subMerge[[aCovariate]] == 1, "FeatureValue"]) > 1
      # feature has a measurement in more than one sample with
      #covaraite status 1

      if (con1 && con3 && con4) {  # MWU test if binary

        aP <- stats::wilcox.test (subMerge [subMerge [[aCovariate]] == 0,
                                            "FeatureValue"],
                                  subMerge [subMerge [[aCovariate]] == 1,
                                            "FeatureValue"])$p.value
        aD <- orddom::orddom (subMerge [subMerge [[aCovariate]] == 0,
                                        "FeatureValue"],
                              subMerge [subMerge [[aCovariate]] == 1,
                                        "FeatureValue"]) [13]

      }

      else if (con2 && con3 && con4) {  # spearman test if continuous

        aP <- stats::cor.test (subMerge [, aCovariate],
                               subMerge [, "FeatureValue"],
                               method = "spearman")$p.value
        aD <- stats::cor.test (subMerge [, aCovariate],
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

  Qs <- matrix (NA, length(features[1:nrow(r)]), length(covariates))
  rownames(Qs) <- features[1:nrow(r)]
  colnames(Qs) <- covariates

  ##
  ##
  if (verbosity == "debug") {
    print(paste0
        ("NaiveAssociation  -  compute multiple testing p adjustment using ",
          adjustMethod))
  }
  ##
  ##

  for (i in 1:ncol(Ps)) {
    Qs[, i] <- stats::p.adjust (Ps[, i], method = adjustMethod)
  }

  return(list(Ps=Ps, Ds=Ds, Qs=Qs))
}
