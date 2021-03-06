#' @importFrom foreach %dopar%
# @importFrom stats na.exclude
#' @import stats
#' @importFrom  lmtest lrtest
#' @import futile.logger
# @export NaiveAssociation

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

  isRobust <- isRobust[[1]]

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
          "isRobust[aCovariate]",
          sep = "\t" ),
        file="progress.txt",
        append = FALSE)
  }
  ##
  ##
  r = foreach::foreach(i= seq_along(features), .combine='rbind') %dopar% {

    somePs <- vector(length = noCovariates)
    someDs <- vector(length = noCovariates)

    if ((var(featureMat[, i], na.rm = TRUE) == 0 ) ) {
      somePs[seq_along(covariates)] <- NA
      someDs[seq_along(covariates)] <- NA

      if ((i %% 10) == 0) {
        progress <- paste0(round(x = ((i/length(features))*100),
                                 digits = 2), "%")
        flog.info(msg = paste("NaiveAssociation -- processed",
                              progress,
                              "of features."),
                  name = "my.logger")
      }

      return(c(as.numeric(somePs), as.numeric(someDs)))
    }

    subMerge <- metaMat
    subMerge$FeatureValue <- featureMat [,i]

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
                  isRobust[aCovariate],
                  sep = "\t" ),
                file="progress.txt",
                append = TRUE)
        }
      }
      ##
      ##

      aD <- NA_real_
      aP <- NA_real_

      if (!is.na(isRobust[j]) && !isRobust[j]) {
              somePs[j] <- aP
              someDs[j] <- aD

              ##
              ##
              if (verbosity == "debug") {
                write(paste0("skipped ", aCovariate),
                      file="progress.txt",
                      append = TRUE)
              }
              ##
              ##

              next
      }

	getVariableType <- function (values, variable) {
		if (is.numeric (values) && all (values %in% c (0, 1)) && ! (! is.null (typeContinuous) && variable %in% typeContinuous) && ! (! is.null (typeCategorical) && variable %in% typeCategorical)) {

			return ("binary") # treat as binary

		} # this fulfils criteria of binary (0, 1) and not in the special cases

		else if ((is.numeric (values) || (! is.null (typeContinuous) && variable %in% typeContinuous)) && ! (! is.null (typeCategorical) && variable %in% typeCategorical)) {

			return ("continuous") # treat as continuous

		} # this fulfils criteria of not being restricted to 0, 1; still numeric, or guaranteed continuous (redundant?); and not categorical

		else {

			return ("categorical") # treat as categorical

		} # default as categorical

	}

	variableType <- getVariableType (na.exclude(subMerge[[aCovariate]]), aCovariate)

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
      formulaNull <- paste0 ("stats::glm (FeatureValue ~ 1, data = subSubMerge, family = \"binomial\")", collapse = "")
      formulaVar <- paste0 ("stats::glm (FeatureValue ~ ", aCovariate, ", data = subSubMerge, family = \"binomial\")", collapse = "")

      lmNull <- eval (parse (text = as.character (formulaNull)))
      lmVar <- eval (parse (text = as.character (formulaVar)))

      aP <- lmtest::lrtest (lmNull, lmVar)$'Pr(>Chisq)' [2]
      aD <- lmVar$coef [2]

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

      aP <- stats::wilcox.test (
          subMerge [subMerge [[aCovariate]] == 0, "FeatureValue"],
          subMerge [subMerge [[aCovariate]] == 1, "FeatureValue"])$p.value

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

      aP <- stats::cor.test (subMerge [, aCovariate],
                             subMerge [, "FeatureValue"],
                             method = "spearman")$p.value
      aD <- stats::cor.test (subMerge [, aCovariate],
                             subMerge [, "FeatureValue"],
                             method = "spearman")$estimate
    }

    else if (variableType == "categorical" && conVar) {  # now never happens, probably 	# SKF20200221
#      else if (con2 && !con5) {  # kruskal-wallis test if
                                  # not binary AND not numerical
      aP <- stats::kruskal.test (
        g = as.factor(subMerge [[aCovariate]]),
        x = subMerge [["FeatureValue"]])$p.value

      aD <- Inf
    }

      somePs[j] <- aP
      someDs[j] <- aD

    } # for j


    if ((i %% 10) == 0) {
      progress <- paste0(round(x = ((i/length(features))*100),digits = 2), "%")
      flog.info(msg = paste("NaiveAssociation -- processed",
                            progress,
                            "of features."),
                name = "my.logger")
    }

    return(c(as.numeric(somePs), as.numeric(someDs)))

  } # for i (foreach)

  parallel::stopCluster(cl)

  flog.info(msg = paste("NaiveAssociation -- processed 100% of features."),
            name = "my.logger")


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

  for (i in seq_len(ncol(Ps))) {
    Qs[, i] <- stats::p.adjust (Ps[, i], method = adjustMethod)
  }
  return(list(Ps=Ps, Ds=Ds, Qs=Qs))
}
