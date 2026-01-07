#' @importFrom foreach %dopar%
#' @importFrom  lmtest lrtest
#' @importFrom  stats glm
#' @importFrom  stats cor.test kruskal.test wilcox.test p.adjust update
#' @import logger
#' @import detectseparation
#' @importFrom reshape2 melt dcast

NaiveAssociation <- function(featureMat,
                             samples,
                             features,
                             noFeatures,
                             metaMat,
                             covariates,
                             noCovariates,
                             isRobust,
                             typeCategorical, # new SKF20200221
                             typeContinuous, # new SKF20200221
                             logistic, # new SKF20201017
                             nnodes,
                             rawCounts, # new TB20221129
                             adjustMethod,
                             adjustLevel,
                             mediationMat
                             ) {

  #new TB20240926
  `%toggleDoPar%` <- `%do%`
  if (nnodes != 1) {
    `%toggleDoPar%` <- `%dopar%`
  }

  #new TB20221125
  if (rawCounts == TRUE) {
    # compute totReadCount per sample and append to metaMat
    totReadCount <- as.data.frame(rowSums(featureMat, na.rm = T))
    metaMat <- merge(metaMat, totReadCount, by = 0, sort = FALSE)
    colnames(metaMat)[ncol(metaMat)] <- "totReadCount"
    row.names(metaMat) <- metaMat$Row.names
    metaMat$Row.names <- NULL
    featureMat <- featureMat/metaMat$totReadCount
  }

  isRobust <- isRobust[[1]]

  # load parralel processing environment
  if (.Platform$OS.type == "unix") {
    # unix
    cl <- parallel::makeForkCluster(nnodes = nnodes, outfile = "")
    doParallel::registerDoParallel(cl)
  } else {
    # windows
    cl <- parallel::makeCluster(nnodes, type = "PSOCK", outfile = "")
    doParallel::registerDoParallel(cl)
  }

  # compute steps for progress log.info #TB20240229
  # aim for no more than 100 steps
  progressSteps <- round(noFeatures/100)
  if (progressSteps == 0) {
    progressSteps <- 1
  }
  if (noFeatures/progressSteps > 100) {
    progressSteps <- progressSteps + 1
  }

  i <- 0
  logger::log_debug(namespace = "metadeconfoundR",
    paste(
      "counter",
      "aCovariate",
      "aFeature",
      "noCovariates" ,
      "isRobust[aCovariate]",
      sep = "\t"
    ))


  r = foreach::foreach(i= seq_along(features), .combine='rbind') %toggleDoPar% {

    somePs <- vector(length = noCovariates)
    someDs <- vector(length = noCovariates)

    if ((var(featureMat[, i], na.rm = TRUE) == 0 ||
         length(na.exclude(featureMat[, i])) < 2) ) {
      # if variance of a feature == 0 OR
        # less than 2 elements of a feature are non-NA:
        # set all associations for that feature to NA
      somePs[seq_along(covariates)] <- NA
      someDs[seq_along(covariates)] <- NA


      if ((i %% progressSteps) == 0) {#TB20240229
        progress <- paste0(round(x = ((i/length(features))*100),
                                 digits = 2), "%")
        logger::log_info(namespace = "metadeconfoundR", paste("NaiveAssociation -- processed",
                              progress,
                              "of features."))
      }

      return(c(as.numeric(somePs), as.numeric(someDs)))
    }

    subMerge <- metaMat
    subMerge$FeatureValue <- featureMat [,i]

    for (j in seq_along(covariates)) {

      aFeature <- as.character (features [i])
      aCovariate <- as.character (covariates [j])

      logger::log_debug(namespace = "metadeconfoundR",
        paste(
          i,
          aCovariate,
          aFeature,
          noCovariates ,
          isRobust[aCovariate],
          sep = "\t"))

      aD <- NA_real_
      aP <- NA_real_

      if (!is.na(isRobust[j]) && !isRobust[j]) {
        somePs[j] <- aP
        someDs[j] <- aD


        logger::log_debug(namespace = "metadeconfoundR", paste0("skipped ", aCovariate))

        next
      }

      variableType <- VarType (na.exclude(subMerge[[aCovariate]]), aCovariate, typeCategorical, typeContinuous)
      conVar <- TRUE
      varX <- na.exclude (cbind (subMerge [[aCovariate]], subMerge [["FeatureValue"]]))

      if (sum (! is.na (subMerge [["FeatureValue"]])) < 1 ||
          # TRUE if only NA-values
          nrow (varX) <= 2 ||
          # TRUE if if there are not at least three rows without NAs
          length (unique (varX [, 1])) < 2 ||
          # TRUE if there are not at least two different metadata values in non-NA subset
          length (unique (varX [, 2])) < 2) {
          # TRUE if there are not at least two different feature values in non-NA subset
        conVar <- FALSE
      }


      subSubMerge <- na.exclude(subMerge[, c(aCovariate, "FeatureValue")])

      if (logistic == TRUE && conVar) {


        formulaVar <- paste0 (
          "stats::glm (FeatureValue ~ ",
          aCovariate,
          ", data = subSubMerge, family = \"binomial\")",
          collapse = ""
        )
        lmVar <- eval (parse (text = as.character (formulaVar)))

        formulaNull <- paste0 ("stats::glm (FeatureValue ~ 1, data = subSubMerge, family = \"binomial\")", collapse = "")
        lmNull <- eval (parse (text = as.character (formulaNull)))

    aP <- NA
    glmSepTested <- update(lmVar, method="detect_separation")
    if (glmSepTested$outcome) {
      logger::log_warn(namespace = "metadeconfoundR",
        paste(
          "Separation for:", aFeature, "and",
          aCovariate))
    } else {
      aP <- lmtest::lrtest (lmNull, lmVar)$'Pr(>Chisq)' [2]
    }

    if (variableType == "categorical") {
      aD <- Inf
    } else if (variableType == "binary") {
      aD <- stats::cor.test (subMerge [, aCovariate],
                             subMerge [, "FeatureValue"],
                             )$estimate
    } else  if (variableType == "continuous") {
      aD <- CliffsDelta(
        as.vector (
          na.exclude (
            subMerge [subMerge [["FeatureValue"]] == 0, aCovariate])),
        as.vector (
          na.exclude (
            subMerge [subMerge [["FeatureValue"]] == 1, aCovariate])))
    }
    #aD <- lmVar$coef [2]

  }

  else if (variableType == "categorical" && conVar) {
    # KW test if false binary and 	# SKF20200221

    aP <- stats::kruskal.test (
      g = as.factor(subMerge [[aCovariate]]),
      x = subMerge [["FeatureValue"]])$p.value

    aD <- Inf
  }

  else if (variableType == "binary" && conVar) {
    # MWU test if binary and 	# SKF20200221

    aP <- suppressWarnings(stats::wilcox.test (
        subMerge [subMerge [[aCovariate]] == 0, "FeatureValue"],
        subMerge [subMerge [[aCovariate]] == 1, "FeatureValue"]))$p.value

    aD <- CliffsDelta(
      as.vector (
        na.exclude (
          subMerge [subMerge [[aCovariate]] == 0, "FeatureValue"])),
      as.vector (
        na.exclude (
          subMerge [subMerge [[aCovariate]] == 1, "FeatureValue"])))
  }

  else if (variableType == "continuous" && conVar) {
    # spearman test if continuous and numerical 	# SKF20200221

    corTestObj <- suppressWarnings(stats::cor.test (subMerge [, aCovariate],
                                                    subMerge [, "FeatureValue"],
                                                    method = "spearman"))
    aP <- corTestObj$p.value
    aD <- corTestObj$estimate
    # aP <- stats::cor.test (subMerge [, aCovariate],
    #                        subMerge [, "FeatureValue"],
    #                        method = "spearman")$p.value
    # aD <- stats::cor.test (subMerge [, aCovariate],
    #                        subMerge [, "FeatureValue"],
    #                        method = "spearman")$estimate
  }

#     else if (variableType == "categorical" && conVar) {  # now never happens, probably 	# SKF20200221
# #      else if (con2 && !con5) {  # kruskal-wallis test if
#                                   # not binary AND not numerical
#       aP <- stats::kruskal.test (
#         g = as.factor(subMerge [[aCovariate]]),
#         x = subMerge [["FeatureValue"]])$p.value
#
#       aD <- Inf
#     }

      somePs[j] <- aP
      someDs[j] <- aD

    } # for j


    if ((i %% progressSteps) == 0) {#TB20240229
      progress <- paste0(round(x = ((i/length(features))*100),digits = 2), "%")
      logger::log_info(namespace = "metadeconfoundR", paste("NaiveAssociation -- processed",
                            progress,
                            "of features."))
    }

    return(c(as.numeric(somePs), as.numeric(someDs)))

  } # for i (foreach)

  # close parallel processing environment
  parallel::stopCluster(cl)


  logger::log_info(namespace = "metadeconfoundR", paste("NaiveAssociation -- processed 100% of features."))
  logger::log_debug(namespace = "metadeconfoundR", paste("NaiveAssociation -- ncol(parallelreturn):", ncol(r)))

  # Ps <- r[, seq_len(ncol(r)/2), drop = F]
  # rownames(Ps) <- features[seq_len(nrow(r))]
  # colnames(Ps) <- covariates
  #
  # Ds <- r[, -(seq_len(ncol(r)/2)), drop = F]
  # rownames(Ds) <- features[seq_len(nrow(r))]
  # colnames(Ds) <- covariates
  #
  # Qs <- matrix (NA, length(features[seq_len(nrow(r))]), length(covariates))
  # rownames(Qs) <- features[seq_len(nrow(r))]
  # colnames(Qs) <- covariates

  if (is.null(ncol(r))) {
    r <- t(r)
  }

  Ps <- r[, seq_len(noCovariates), drop = F]
  rownames(Ps) <- features
  colnames(Ps) <- covariates

  Ds <- r[, -(seq_len(noCovariates)), drop = F]
  rownames(Ds) <- features
  colnames(Ds) <- covariates

  Qs <- matrix (NA, length(features), length(covariates))
  rownames(Qs) <- features
  colnames(Qs) <- covariates

  if (adjustLevel == 1) {
    for (i in seq_len(noCovariates)) {
      Qs[, i] <- stats::p.adjust (Ps[, i], method = adjustMethod)
    }
  }
  else if (adjustLevel == 2) {
    Ps_long <- reshape2::melt(Ps, varnames = c("feature", "metaVariable"), value.name = "Ps")
    Ps_long$Ps <- stats::p.adjust (Ps_long$Ps, method = adjustMethod)

    Qs <- reshape2::dcast(Ps_long, feature ~ metaVariable, value.var = "Ps")
    rownames(Qs) <- Qs$feature
    Qs$feature <- NULL
    Qs <- as.matrix(Qs)
  }
  else if (adjustLevel == 3) {
    Ps_long <- reshape2::melt(Ps, varnames = c("feature", "metaVariable"), value.name = "Ps")
    Ps_long$Qs <- Ps_long$Ps
    Ps_long$Qs[Ps_long$metaVariable %in% colnames(mediationMat)] <-
      stats::p.adjust (Ps_long$Qs[Ps_long$metaVariable %in% colnames(mediationMat)], method = adjustMethod)
    Ps_long$Ps <- Ps_long$Qs
    Ps_long$Qs <- NULL
    Qs <- reshape2::dcast(Ps_long, feature ~ metaVariable, value.var = "Ps")
    rownames(Qs) <- Qs$feature
    Qs$feature <- NULL
    Qs <- as.matrix(Qs)
  }

  return(list(Ps=Ps, Ds=Ds, Qs=Qs))
}
