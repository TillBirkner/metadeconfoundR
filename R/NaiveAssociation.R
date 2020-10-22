#' @importFrom foreach %dopar%
# @importFrom stats na.exclude
#' @import stats
#' @importFrom  lmtest lrtest
#' @import futile.logger
# @export NaiveAssociation
#



# being able to specify a range of features (only columns 50-100) (low relevance)
# report all confounders, not only the first one # fixed quick and dirty, might need to be revisited when deciding on new data outuput formats

# whitelist/blacklist for metaVariables that should or schould not be rank transformed
  # DONE: within one for loop all variabels listed in doRanks will be rank transformed in the "subMerge" data.frame

# for the lm that is shared between forward and reverse call do optional: confint(lm(xxx)) and check if confint spans over 0
  # DONE (without categorical covariates!!)

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

    #new approach --> append i-th feature to the metaMat dataframe
    subMerge <- metaMat
    subMerge$FeatureValue <- featureMat [,i] # caveat *1*: see snippets.R
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

      # subFeatures <- featureMat [,i]
      # subMerge <- md
      # subMerge$FeatureValue <- subFeatures # caveat *1*: see snippets.R


	## BEGIN new section SKF20200221

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

	## END new section SKF20200221

	## BEGIN section commented out SKF20200221

#      con1 <- length (unique (na.exclude(subMerge[[aCovariate]]))) == 2
#        # the distribution of the covariate is binary
#      con2 <- length (unique (na.exclude(subMerge[[aCovariate]]))) > 2
#        # the distribution of the covariate is continuous (more than 2 levels)
#      con3 <- length (
#                na.exclude(
#                  subMerge[subMerge[[aCovariate]] == 0, "FeatureValue"])) > 1
#        # feature has a measurement in more than one sample with
#          #covaraite status == 0
#      con4 <- length (
#                na.exclude(
#                  subMerge[subMerge[[aCovariate]] == 1, "FeatureValue"])) > 1
#        # feature has a measurement in more than one sample with
#          #covaraite status == 1
#      con5 <- is.numeric(subMerge[[aCovariate]])
#        # covariate is true numeric (distinguishes between continuous
#          # numeric data and "level" data, that is converted to numbers)


	## END section commented out SKF20200221


      # if (verbosity == "debug" && i == 1) {
      #   write(paste(covariates[j], con1, con2, con3, con4, con5),
      #         file="data_type.txt",
      #         append = TRUE)
      # }

  conVar <- TRUE
  varX <- na.exclude (cbind (subMerge [[aCovariate]], subMerge [["FeatureValue"]]))

  if (sum (! is.na (subMerge [["FeatureValue"]])) < 1 || # TRUE if only NA-values
      nrow (varX) <= 2 || # TRUE if if there are not at least three rows without NAs
      length (unique (varX [, 1])) < 2 || # TRUE if there are not at least two different metadata values in non-NA subset
      length (unique (varX [, 2])) < 2) { # TRUE if there are not at least two different feature values in non-NA subset
        conVar <- FALSE
    }

      #varX <- as.data.frame(varX)
      #colnames(varX) <- c(aCovariate, "FeatureValue")
      subSubMerge <- na.exclude(subMerge[, c(aCovariate, "FeatureValue")])
      if (logistic == TRUE && conVar) {
        formulaNull <- paste0 ("stats::glm (FeatureValue ~ 1, data = subSubMerge)", collapse = "")
        formulaVar <- paste0 ("stats::glm (FeatureValue ~ ", aCovariate, ", data = subSubMerge)", collapse = "")

        lmNull <- eval (parse (text = as.character (formulaNull)))
        lmVar <- eval (parse (text = as.character (formulaVar)))

        aP <- lmtest::lrtest (lmNull, lmVar)$'Pr(>Chisq)' [2]
        aD <- Inf

      }

      else if (variableType == "categorical" && conVar) {  # KW test if false binary and 	# SKF20200221
#      if (con1 && !con5) {  # KW test if false binary and

        aP <- stats::kruskal.test (
          g = as.factor(subMerge [[aCovariate]]),
          x = subMerge [["FeatureValue"]])$p.value

        aD <- Inf

      }

      else if (variableType == "binary" && conVar) {  # MWU test if binary and 	# SKF20200221
#      else if (con1 && con3 && con4 && con5) {  # MWU test if binary and

        aP <- stats::wilcox.test (
            subMerge [subMerge [[aCovariate]] == 0, "FeatureValue"],
            subMerge [subMerge [[aCovariate]] == 1, "FeatureValue"])$p.value
        # aD <- orddom::orddom (
        #   as.vector (
        #     na.exclude (
        #       subMerge [subMerge [[aCovariate]] == 0, "FeatureValue"])),
        #   as.vector (
        #     na.exclude (
        #       subMerge [subMerge [[aCovariate]] == 1, "FeatureValue"]))) [13]
        aD <- CliffsDelta(
          as.vector (
            na.exclude (
              subMerge [subMerge [[aCovariate]] == 0, "FeatureValue"])),
          as.vector (
            na.exclude (
              subMerge [subMerge [[aCovariate]] == 1, "FeatureValue"])))
      }

      else if (variableType == "continuous" && conVar) {  # spearman test if continuous and numerical 	# SKF20200221
#      else if (con2 && con5) {  # spearman test if continuous and numerical

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

    }


    if ((i %% 10) == 0) {
      progress <- paste0(round(x = ((i/length(features))*100),digits = 2), "%")
      flog.info(msg = paste("NaiveAssociation -- processed",
                            progress,
                            "of features."),
                name = "my.logger")
    }

    return(c(as.numeric(somePs), as.numeric(someDs)))

  }
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
