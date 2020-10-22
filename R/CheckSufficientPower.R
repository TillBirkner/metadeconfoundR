#' @importFrom foreach %dopar%
#' @import futile.logger
# @export CheckSufficientPower
#

## ---- test-a ----
CheckSufficientPower <- function(metaMat,
                                 covariates,
                                 noCovariates,
                                 nnodes,
                                 robustCutoff,
				 robustCutoffRho, # new SKF20200221
				 typeCategorical, # new SKF20200221
		 		 typeContinuous, # new SKF20200221
                                 NA_imputation,
                                 maintenance,
                                 verbosity) {

## ---- test-b ---
  #print(dim(metaMat))
  #print(metaMat[1:5, 1:5])
  #print(length(metaMat[1, ]) == length(covariates))
  #metaMat <- na.exclude(metaMat)
    # exclude NA entries from metaMat
  noCondition <- nrow(metaMat[metaMat[,1] == 1, ])
    # number of samples with Status == 1
  noControl <- nrow(metaMat[metaMat[,1] == 0, ])
    # number of samples with Status == 0
  conditionMat <- metaMat[metaMat[,1] == 1, ]

  ##
  ##
  if (verbosity == "debug") {
    print(paste("CheckSufficientPower -- covariates:", c(paste(covariates))))
    print(paste("CheckSufficientPower -- noCovariates:", noCovariates))
    #print(paste("CheckSufficientPower -- dim(conditionMat): ",
    #            nrow(conditionMat)))
  }
  ##
  ##

  if (noCondition < robustCutoff | noControl < robustCutoff) {
    flog.error(msg = paste("Not enough(robustCutoff =", robustCutoff, ") samples in either case or controle group."),
               name = "my.logger")
    stop(paste("Not enough(robustCutoff =", robustCutoff, ") samples in either case or controle group."))
  }

  robustSingle <- vector(length = length(covariates), mode = "logical")
  names(robustSingle) <- covariates

  robustCombination <- as.data.frame(matrix(nrow = length(covariates), ncol = length(covariates), data = FALSE))
  rownames(robustCombination) <- colnames(robustCombination) <- covariates

  #cl <- parallel::makeForkCluster(nnodes = nnodes, outfile="")
    # the parent process uses another core
      #(so 4 cores will be used with this command)
  #doParallel::registerDoParallel(cl)

  #i <- 0
  #parallelReturn <- foreach::foreach(
 #   i= seq_along(covariates),
 #   .combine = 'rbind') %dopar% {


  # dataType <- function(data) {
  #   con1 <- length (unique (na.exclude(data))) == 2
  #   # the distribution of the covariate is binary
  #   con2 <- length (unique (na.exclude(data))) > 2
  #   # the distribution of the covariate is continuous (more than 2 levels)
  #   con5 <- is.numeric(data)
  #   # covariate is true numeric (distinguishes between continuous
  #   # numeric data and "level" data, that is converted to numbers)
  #   if (con1 && con5) {
  #     return("binary")
  #   } else if (con2 && con5) {
  #     return("continuous")
  #   } else if ((con1 | con2) && !con5) {
  #     return("categorical")
  #   } else {
  #     return("empty")
  #     #print(paste0("could not compute data type for ", i))
  #     #print(data[1:5])
  #   }
  # }

  for (i in seq_along(covariates)) {


    #print(metaMat[1:5, i])
    #iDataType <- dataType(metaMat[, i])
    #print(iDataType)
    # if (iDataType == "continuous") {
    #   #if (length(which(metaMat[, i] > 0)) >= robustCutoff) {
    #   if (sum (as.numeric (metaMat[, i]) != median (as.numeric (metaMat[, i]))) >= robustCutoff)
    #     robustSingle[i] <- TRUE
    #   }
    # } else if (iDataType == "binary") {
    #   if ((length(which(metaMat[, i] == 0)) >= robustCutoff) &&
    #       (length(which(metaMat[, i] == 1)) >= robustCutoff)) {
    #     robustSingle[i] <- TRUE
    #   }
    # } else if (iDataType == "categorical") {
    #   #if (all(summary(as.factor(metaMat[, i])) >= robustCutoff)) {
    #   if (sort (summary(as.factor(metaMat[, i]))) [2] >= robustCutoff) {
    #     robustSingle[i] <- TRUE
    #   }
    # }

	## BEGIN new section SKF20200221
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
  ## ---- robustSingle

	variableType <- getVariableType (metaMat [, i], covariates [i])

	if (variableType == "binary" || variableType == "categorical") {

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

	## END new section SKF20200221

	## BEGIN section commented out SKF20200221


#	  if (sum (as.numeric (as.factor(metaMat[, i])) != median (as.numeric (as.factor(metaMat[, i])), na.rm = TRUE), na.rm = TRUE) >= robustCutoff) {
#	    robustSingle[i] <- TRUE
#	  }

	## END section commented out SKF20200221

    for (j in seq_along(covariates)) {
      if (i == j) {
        robustCombination[i,j] <- TRUE
        next
      }


      toCompare <- na.exclude(metaMat[, c(i,j)])

	## BEGIN new section SKF20200221

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

	## END new section SKF20200221

	## BEGIN section commented out SKF20200221

#      toCompare[, 1] <- as.factor(toCompare[, 1])
#      toCompare[, 2] <- as.factor(toCompare[, 2])


      # print(dim(toCompare))
      # print(toCompare[, 2])
      #
      # print(median (as.numeric (toCompare[, 1])))
      # print(sum (as.numeric (toCompare[, 1]) != median (as.numeric (toCompare[, 1]))))
      #
      # print(median (as.numeric (as.factor(toCompare[, 2]))))
      # print(sum (as.numeric (toCompare[, 2]) != median (as.numeric (toCompare[, 2]))))


#      if ((sum (as.numeric (toCompare[, 1]) != median (as.numeric (toCompare[, 1]))) >= robustCutoff) &&
#          (sum (as.numeric (toCompare[, 2]) != median (as.numeric (toCompare[, 2]))) >= robustCutoff)) {
#        robustCombination[i,j] <- TRUE
#        robustCombination[j,i] <- TRUE
#      }

	## END section commented out SKF20200221


      # iDataType <- dataType(toCompare[, 1])
      # jDataType <- dataType(toCompare[, 2])
      #
      # if ((iDataType == "binary") && (jDataType %in% c("binary", "categorical"))) {
      #   if (length(which(table(toCompare)[2, ] < robustCutoff)) == 0) {
      #     robustCombination[i,j] <- TRUE
      #     robustCombination[j,i] <- TRUE
      #   }
      # } else if ((iDataType == "categorical") && (jDataType == "categorical")) {
      #   if (length(which(table(toCompare) < robustCutoff)) == 0) {
      #     robustCombination[i,j] <- TRUE
      #     robustCombination[j,i] <- TRUE
      #   }
      # } else if ((iDataType == "binary") && (jDataType == "continuous")) {
      #   if (length(which(toCompare[toCompare[, 1] == 1, 2] > 0)) >= robustCutoff) {
      #     robustCombination[i,j] <- TRUE
      #     robustCombination[j,i] <- TRUE
      #   }
      # } else if (iDataType == jDataType) { # both are continuous
      #   robustCombination[i,j] <- TRUE
      # } else if (iDataType == "categorical" && jDataType == "continuous") {
      #   gatherer <- TRUE
      #   for (k in unique(toCompare[, 1])) {
      #     gatherer <- c(gatherer, (length(which(toCompare[toCompare[, 1] == k, 2] > 0)) >= robustCutoff))
      #   }
      #   if (all(gatherer)) {
      #     robustCombination[i,j] <- TRUE
      #     robustCombination[j,i] <- TRUE
      #   }
      # }
    } # for j

  } # for i

#
#     aCovariate <- as.character (covariates [i])
#
#
#     if (i == 1) {
#       robustCombination <- TRUE
#       return(data.frame
#              (row.names = aCovariate,
#                noControl,
#                robustCombination))
#     }
#
#     robustCombination <- FALSE
#
#     if ((length(unique (metaMat[[aCovariate]])) == 2) && # binary covariate
#          (length(which(table(metaMat[,c(1,i)])[2, ] < robustCutoff)) == 0)) {
#           # in second row of contingency table (status == 1), no entry
#             # must be smaller than robustCutoff
#             # --> there must be at least 10 samples each, of "diseased" that
#             # are positive or negative for the covariate
#       robustCombination <- TRUE
#     } else if (length(unique (metaMat[[aCovariate]])) > 2) {  # not binary
#       if (is.numeric(metaMat[[aCovariate]])) {  # continuous coavariate
#         if (length(which(conditionMat[[aCovariate]] > 0)) >= robustCutoff) {
#         # if at least ~robustCutoff~ samples of
#           # "deseased"-status are bigger than 0
#           robustCombination <- TRUE
#         }
#       }
#       else if(length(
#                 which(
#                   table(metaMat[,c(1, i)])[2, ] < robustCutoff)) == 0) {
#         # nominal variable
#         robustCombination <- TRUE
#       }
#
#     }
#
#
#     return(data.frame
#            (row.names = aCovariate,
#              noControl,
#              robustCombination))
#
#   }
#
#   parallel::stopCluster(cl)
#
#   if (sum(parallelReturn$robustCombination == FALSE) > 0) {
#
#     flog.info(
#       msg = paste(
#         c("following variables are filtered out due to insufficient data:",
#           rownames(parallelReturn)[parallelReturn$robustCombination == FALSE]),
#         collapse = " "),
#       name = "my.logger")
#
#   }

  #return(parallelReturn)
  return(list(robustSingle, robustCombination))
}
