#' @importFrom foreach %dopar%
# @export CheckSufficientPower
#

## ---- test-a ----
CheckSufficientPower <- function(metaMat,
                                 covariates,
                                 noCovariates,
                                 nnodes,
                                 robustCutoff,
                                 maintenance,
                                 verbosity) {

## ---- test-b ---

  conditionMat <- metaMat[metaMat[,1] == 1, ] # extract all "positive" samples

  noCondition <- nrow(conditionMat)
  noControl <- nrow(metaMat) - noCondition

  ##
  ##
  if (verbosity == "debug") {
    print(paste("CheckSufficientPower -- covariates:", paste(covariates)))
    print(paste("CheckSufficientPower -- noCovariates:", noCovariates))
    print(paste("CheckSufficientPower -- dim(conditionMat): ",
                nrow(conditionMat)))
    # write(paste("covariate",
    #             "control",
    #             "condition_negative",
    #             "condition_positive",
    #             sep = "\t"),
    #       file = "testfile.txt",
    #       append = FALSE)
  }
  ##
  ##

  if (noCondition < 10 | noControl < 10) {
    stop("Number of samples with status (first column metaMat)
         0 or 1 below robustCutoff")
  }

  cl <- parallel::makeForkCluster(nnodes = nnodes, outfile="")
    # the parent process uses another core
      #(so 4 cores will be used with this command)
  doParallel::registerDoParallel(cl)

  i <- 0
  parallelReturn <- foreach::foreach(i= 2:noCovariates,
                                     .combine = 'rbind') %dopar% {

    aCovariate <- as.character (covariates [i])
    # condition <- table(eval
    #                    (parse
    #                      (text = as.character
    #                        (paste0
    #                          ("conditionMat$",
    #                            aCovariate)))))
      # count number of individuals having the phenotype and not
        #taking/taking the covariate
    # noConditionNegative <- condition[1]
    # noConditionPositive <- condition[2]
    robustCombination <- FALSE
#
#     if ((noConditionNegative >= robustCutoff) &&
#         (noConditionPositive >= robustCutoff)) {
#       robustCombination <- TRUE
#     }



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
        # if at least robustCutoff samples of
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


    ##
    ##
    # if (verbosity == "debug") {
    #   write(paste
    #         (aCovariate,
    #           noControl,
    #           noConditionNegative,
    #           noConditionPositive,
    #           robustCombination,
    #           sep = "\t"),
    #         file = "testfile.txt",
    #         append = TRUE)
    # }
    ##
    ##

    return(data.frame
           (row.names = aCovariate,
             noControl,
             #noConditionNegative,
             #noConditionPositive,
             robustCombination))

  }

  parallel::stopCluster(cl)

  return(parallelReturn)
}
