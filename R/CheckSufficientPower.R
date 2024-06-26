#' @import futile.logger


CheckSufficientPower <- function(metaMat,
                                 covariates,
                                 noCovariates,
                                 nnodes,
                                 robustCutoff,
                                 robustCutoffRho, # new SKF20200221
                                 typeCategorical, # new SKF20200221
                                 typeContinuous, # new SKF20200221
                                 verbosity,
                                 RVnames,
                                 startStop,
                                 deconfF # new TB20220704
                                 )  {


  # global check for sufficient samples of each case and control groups
  noCondition <- nrow(metaMat[metaMat[,1] == 1, ])
    # number of samples with Status == 1
  noControl <- nrow(metaMat[metaMat[,1] == 0, ])
    # number of samples with Status == 0
  conditionMat <- metaMat[metaMat[,1] == 1, ]

  if (ncol(metaMat) == 1) {
    noCondition <- length(metaMat[metaMat[,1] == 1, ])
    noControl <- length(metaMat[metaMat[,1] == 0, ])
  }

  flog.debug(paste("CheckSufficientPower -- covariates: ",
                   paste(covariates, collapse = ", ")),
             name = "my.logger")
  flog.debug(paste("CheckSufficientPower -- noCovariates:",
                   noCovariates),
             name = "my.logger")


  if (noCondition < robustCutoff | noControl < robustCutoff) {
    flog.error(msg = paste("Not enough(robustCutoff =", robustCutoff, ") samples in either case or controle group."),
               name = "my.logger")
    stop(paste("Not enough(robustCutoff =", robustCutoff, ") samples in either case or controle group."))
  }

  robustSingle <- vector(length = length(covariates), mode = "logical")
  names(robustSingle) <- covariates

  robustCombination <- as.data.frame(matrix(nrow = length(covariates),
                                            ncol = length(covariates),
                                            data = FALSE))
  rownames(robustCombination) <- colnames(robustCombination) <- covariates

  ## ---- getVariableType
  getVariableType <- function (values, variable) {
    # variable == name of covariate
    if (is.numeric (values) &&
        all (values %in% c (0, 1)) &&
        ! (! is.null(typeContinuous) && variable %in% typeContinuous) &&
        ! (! is.null(typeCategorical) && variable %in% typeCategorical)) {
      return ("binary") # treat as binary
    }
    else if ((is.numeric (values) ||
              (! is.null (typeContinuous) && variable %in% typeContinuous)) &&
             ! (! is.null (typeCategorical) && variable %in% typeCategorical)) {
      return ("continuous") # treat as continuous
    }
    else {
      return ("categorical") # treat as categorical
    }
  }

  for (i in seq_along(covariates)) {

    if (!is.na(RVnames[[1]]) && covariates[i] %in% RVnames) {
      # all random effect variables are labeled as "not robust"
      next
    }

    variableType <- getVariableType (metaMat [, i], covariates [i])

    if (variableType == "binary" || variableType == "categorical") {
      # how many samples have a value != the majority value
        # if median(var) == 0: are there >robustCutoff samples with var == 1
      if (sum (as.numeric (as.factor(metaMat[, i])) !=
               median (as.numeric (as.factor(metaMat[, i])),
                       na.rm = TRUE),
               na.rm = TRUE) >= robustCutoff) {
        robustSingle[i] <- TRUE
      }

    }

    if (variableType == "continuous") {

      if (sum (metaMat[, i] != median (metaMat[, i], na.rm = TRUE),
               na.rm = TRUE) >= robustCutoffRho) {
        robustSingle[i] <- TRUE
      }

    }

    # skip combinations, if "naiveStop"
    if ((!is.na(startStop[[1]]) && startStop == "naiveStop")) {
      next
    }

    # new TB20220704
    # if a covariate is listed in deconfF, set its whole column to FALSE
    if (covariates[i] %in% deconfF) {
      robustCombination[i,] <- FALSE
      next
    }

    for (j in seq_along(covariates)) {
      if (i == j) {
        robustCombination[i,j] <- TRUE
        next
      }

      toCompare <- na.exclude(metaMat[, c(i,j)])

      variableType1 <- getVariableType (toCompare [, 1], covariates [i])
      variableType2 <- getVariableType (toCompare [, 2], covariates [j])

      robust1 <- FALSE

      if (variableType1 == "binary" || variableType1 == "categorical") {
    	  if (sum (as.numeric (as.factor(toCompare [, 1])) != median (as.numeric (as.factor(toCompare [, 1])), na.rm = TRUE), na.rm = TRUE) >= robustCutoff) {
    	    robust1 <- TRUE
    	  }
    	}

      if (variableType1 == "continuous") {

    	  if (sum (toCompare [, 1] != median (toCompare [, 1], na.rm = TRUE), na.rm = TRUE) >= robustCutoffRho) {
    	    robust1 <- TRUE # robustCutoffRho == 1 matches older behaviour of MetaCardis pipeline: at least two values for continuous parameter
    	  }
    	}

      robust2 <- FALSE

      if (variableType2 == "binary" || variableType2 == "categorical") {
    	  if (sum (as.numeric (as.factor(toCompare [, 2])) != median (as.numeric (as.factor(toCompare [, 2])), na.rm = TRUE), na.rm = TRUE) >= robustCutoff) {
    	    robust2 <- TRUE
    	  }
    	}

      if (variableType2 == "continuous") {
    	  if (sum (toCompare [, 2] != median (toCompare [, 2], na.rm = TRUE), na.rm = TRUE) >= robustCutoffRho) {
    	    robust2 <- TRUE # robustCutoffRho == 1 matches older behaviour of MetaCardis pipeline: at least two values for continuous parameter
    	  }
    	}

      if (robust1 && robust2) {
    		robustCombination[i,j] <- TRUE
            	robustCombination[j,i] <- TRUE
      }
    } # for j
  } # for i
  return(list(robustSingle, robustCombination))
}
