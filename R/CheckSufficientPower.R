#' @importFrom foreach %dopar%
# @export CheckSufficientPower

## ---- test-a ----
CheckSufficientPower <- function(metaMat,
                                 nnodes,
                                 robustCutoff,
                                 maintenance,
                                 verbosity) {

## ---- test-b ---
  covariates <- colnames (metaMat)
  # each covariate + the status category, specific to example
  noCovariates <- length (covariates)

  conditionMat <- metaMat[metaMat[,1] == 1, ]  # extract all "positive" samples
  print(paste("CheckSufficientPower  -  dim(conditionMat): ",
              nrow(conditionMat)))
  noCondition <- nrow(conditionMat)
  noControl <- nrow(metaMat) - noCondition

  ##
  ##
  if (verbosity == "debug") {
    write(paste("covariate",
                "control",
                "condition_negative",
                "condition_positive",
                sep = "\t"),
          file = "testfile.txt",
          append = FALSE)
  }
  ##
  ##

  cl <- parallel::makeForkCluster(nnodes = nnodes, outfile="")
    # the parent process uses another core
      #(so 4 cores will be used with this command)
  doParallel::registerDoParallel(cl)

  i <- 0
  parallelReturn <- foreach::foreach(i= 2:noCovariates,
                                     .combine = 'rbind',
                                     .export = "conditionMat") %dopar% {

    aCovariate <- as.character (covariates [i])
    condition <- table(eval
                       (parse
                         (text = as.character
                           (paste0
                             ("conditionMat$",
                               aCovariate)))))
      # count number of individuals having the phenotype and not
        #taking/taking the covariate
    noConditionNegative <- condition[1]
    noConditionPositive <- condition[2]
    robustCombination <- FALSE

    if ((noControl >= robustCutoff) &&
        (noConditionNegative >= robustCutoff) &&
        (noConditionPositive >= robustCutoff)) {
      robustCombination <- TRUE
    }

    ##
    ##
    if (verbosity == "debug") {
      write(paste
            (aCovariate,
              noControl,
              noConditionNegative,
              noConditionPositive,
              robustCombination,
              sep = "\t"),
            file = "testfile.txt",
            append = TRUE)
    }
    ##
    ##

    return(data.frame
           (row.names = aCovariate,
             noControl,
             noConditionNegative,
             noConditionPositive,
             robustCombination))

  }

  parallel::stopCluster(cl)

  return(parallelReturn)
}
