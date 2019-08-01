#' @importFrom foreach %dopar%
# @importFrom stats na.exclude
#' @import stats
# @export NaiveAssociation
#

NaiveAssociation <- function(featureMat,
                             samples,
                             features,
                             noFeatures,
                             metaMat,
                             covariates,
                             noCovariates,
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

  # load parralel processing environment
  cl <- parallel::makeForkCluster(nnodes = nnodes, outfile = "")
  # the parent process uses another core (nnodes - 1 might be good)
  doParallel::registerDoParallel(cl)
  i <- 0

  ##
  ##
  if (verbosity== "debug") {
    cat("NaiveAssociation -- \n\tnoSamples: ",
        length(samples),
        "\n\tnoFeatures: ",
        noFeatures,
        "\n\tnoCovariates: ",
        noCovariates)
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
  ##
  ##
  r = foreach::foreach(i= seq_along(features), .combine='rbind') %dopar% {

    somePs <- vector(length = noCovariates)
    someDs <- vector(length = noCovariates)

    if (var(featureMat[, i], na.rm = TRUE) == 0 ) {
      somePs[seq_along(covariates)] <- 0
      someDs[seq_along(covariates)] <- 0
      return(c(as.numeric(somePs), as.numeric(someDs)))
    }



    # if (anyNA(featureMat [, i])) {  # if feature i contains NAs
    #                                   # exclude NA rows from feature and meta
    #   subFeatures <- na.exclude(featureMat [, i])
    #   subMerge <- metaMat[-(stats::na.action(subFeatures)),]
    #   subMerge$FeatureValue <- as.vector(subFeatures) # caveat *1*
    # } else {
    #   subFeatures <- featureMat [,i]
    #   subMerge <- metaMat
    #   subMerge$FeatureValue <- subFeatures # caveat *1*: see snippets.R
    # }

    #new approach
    subFeatures <- featureMat [,i]
    subMerge <- metaMat
    subMerge$FeatureValue <- subFeatures
    #subMerge <- as.data.frame(na.exclude(subMerge))


    for (j in seq_along(covariates)) {


      aFeature <- as.character (features [i])
      aCovariate <- as.character (covariates [j])

      ##
      ##
      if (verbosity == "debug") {
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
      }
      ##
      ##

      aD <- NA_real_
      aP <- NA_real_

      if (!is.na(isRobust[aCovariate, 2]) && !isRobust[aCovariate, 2]) {
              somePs[j] <- aP
              someDs[j] <- aD

              ##
              ##
              if (verbosity == "debug") {
                write("skipped",
                      file="progress.txt",
                      append = TRUE)
              }
              ##
              ##

              next
      }

      # subFeatures <- featureMat [,i]
      # subMerge <- md
      # subMerge$FeatureValue <- subFeatures # caveat *1*: see snippets.R

      #con1 <- length (unique (subMerge [, aCovariate])) == 2
      con1 <- length (unique (subMerge[[aCovariate]])) == 2
      # the distribution of the covariate is binary
      #con2 <- length (unique (subMerge [, aCovariate])) > 2
      con2 <- length (unique (subMerge[[aCovariate]])) > 2
      # the distribution of the covariate is continuous
      con3 <- length (
        na.exclude (
          subMerge[subMerge[[aCovariate]] == 0, "FeatureValue"])) > 1
      # feature has a measurement in more than one sample with
      #covaraite status 0
      con4 <- length (
        na.exclude (
          subMerge[subMerge[[aCovariate]] == 1, "FeatureValue"])) > 1
      # feature has a measurement in more than one sample with
      #covaraite status 1
      con5 <- is.numeric(subMerge[[aCovariate]])
      # covariate is true numeric (distinguishes between continuous
        # numeric data and "level" data, that is converted to numbers)

      if (con1 && con3 && con4) {  # MWU test if binary and

        aP <- stats::wilcox.test (
          na.exclude (
            subMerge [subMerge [[aCovariate]] == 0, "FeatureValue"]),
          na.exclude (
            subMerge [subMerge [[aCovariate]] == 1, "FeatureValue"]))$p.value
        aD <- orddom::orddom (
          as.vector (
            na.exclude (
              subMerge [subMerge [[aCovariate]] == 0, "FeatureValue"])),
          as.vector (
            na.exclude (
              subMerge [subMerge [[aCovariate]] == 1, "FeatureValue"]))) [13]
      }

      else if (con2 && con5) {  # spearman test if continuous and numerical

        aP <- stats::cor.test (subMerge [, aCovariate],
                               subMerge [, "FeatureValue"],
                               method = "spearman")$p.value
        aD <- stats::cor.test (subMerge [, aCovariate],
                               subMerge [, "FeatureValue"],
                               method = "spearman")$estimate
      }

      else if (con2 && !con5) {  # kruskal-wallis test if
                                  # not binary AND not numerical
        aP <- stats::kruskal.test (
          g = as.factor(subMerge [[aCovariate]]),
          x = subMerge [["FeatureValue"]])$p.value

        aD <- Inf
      }

      somePs[j] <- aP
      someDs[j] <- aD

    }

    return(c(as.numeric(somePs), as.numeric(someDs)))

  }
  parallel::stopCluster(cl)


  if (verbosity == "debug") {
    print(paste("NaiveAssociation -- ncol(parallelreturn):", ncol(r)))
  }

  Ps <- r[, seq_len(ncol(r)/2)]
  rownames(Ps) <- features[seq_len(nrow(r))]
  colnames(Ps) <- covariates

  Ds <- r[, -(seq_len(ncol(r)/2))]
  rownames(Ds) <- features[seq_len(nrow(r))]
  colnames(Ds) <- covariates

  Qs <- matrix (NA, length(features[seq_len(nrow(r))]), length(covariates))
  rownames(Qs) <- features[seq_len(nrow(r))]
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

  for (i in seq_len(ncol(Ps))) {
    Qs[, i] <- stats::p.adjust (Ps[, i], method = adjustMethod)
  }

  return(list(Ps=Ps, Ds=Ds, Qs=Qs))
}
