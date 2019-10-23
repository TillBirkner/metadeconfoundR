#' @importFrom foreach %dopar%
# @export CheckSufficientPower
#

## ---- test-a ----
CheckSufficientPower <- function(metaMat,
                                 covariates,
                                 noCovariates,
                                 nnodes,
                                 robustCutoff,
                                 NA_imputation,
                                 maintenance,
                                 verbosity) {

## ---- test-b ---

  metaMat <- na.exclude(metaMat)
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

  if (noCondition < 10 | noControl < 10) {
    stop(paste("Number of samples with status (first column metaMat)
         0 or 1 below robustCutoff", noControl, noCondition))
  }

  cl <- parallel::makeForkCluster(nnodes = nnodes, outfile="")
    # the parent process uses another core
      #(so 4 cores will be used with this command)
  doParallel::registerDoParallel(cl)

  i <- 0
  parallelReturn <- foreach::foreach(i= seq_along(covariates),
                                     .combine = 'rbind') %dopar% {

    aCovariate <- as.character (covariates [i])

    if (i == 1) {
      robustCombination <- TRUE
      return(data.frame
             (row.names = aCovariate,
               noControl,
               robustCombination))
    }

    robustCombination <- FALSE

    if ((length(unique (metaMat[[aCovariate]])) == 2) && # binary covariate
         (length(which(table(metaMat[,c(1,i)])[2, ] < robustCutoff)) == 0)) {
          # in second row of contingency table (status == 1), no entry
            # must be smaller than robustCutoff
            # --> there must be at least 10 samples each, of "diseased" that
            # are positive or negative for the covariate
      robustCombination <- TRUE
    } else if (length(unique (metaMat[[aCovariate]])) > 2) {  # not binary
      if (is.numeric(metaMat[[aCovariate]])) {  # continuous coavariate
        if (length(which(conditionMat[[aCovariate]] > 0)) >= robustCutoff) {
        # if at least ~robustCutoff~ samples of
          # "deseased"-status are bigger than 0
          robustCombination <- TRUE
        }
      }
      else if(length(
                which(
                  table(metaMat[,c(1, i)])[2, ] < robustCutoff)) == 0) {
        # nominal variable
        robustCombination <- TRUE
      }

    }


    return(data.frame
           (row.names = aCovariate,
             noControl,
             robustCombination))

  }

  parallel::stopCluster(cl)

  return(parallelReturn)
}
