#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Build and Reload Package:  'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'


#' Check sufficiant power of all possible covariate groups
#'
#' @param metaFile a tab delimited file or data frame with row and colum names listing metadata for all samples
#' @param nnodes number of nodes to be used for parallel computing
#' @param cutoff minimamal number of sample size for each covariate in order to have sufficient power, default = 5
#' @return data frame containing sample sizes for all covariates and sufficientPower as yes(1)/no(0)
#' @importFrom foreach %dopar%
#' @export CheckSufficientPower
CheckSufficientPower <- function(metaFile, nnodes=1, cutoff=5) {

  metaMat <- ##### is this necessary? user should import file themselves
    if (is.character(metaFile)) {
      utils::read.table(file = metaFile, header = T, sep = "\t", row.names = 1)
    }
    else if (is.data.frame(metaFile)) {
      return(metaFile)
    } else {
      stop('Wrong metaFile argument! Needs to be either "path/to/file" or dataframe with row and column names.')
    }
  covariates <- colnames (metaMat) [4:27] # each covariate + the status category, specific to example
  noCovariates <- length (covariates)

  conditionMat <- metaMat[metaMat$Case_status == "Schizo", ]
  noCondition <- nrow(conditionMat)
  noControl <- nrow(metaMat) - noCondition

  write(paste("covariate",
              "control",
              "condition_negative",
              "condition_positive",
              sep = "\t"),
        file = "testfile.txt")
  cl <- parallel::makeForkCluster(nnodes = nnodes, outfile="")  # the parent process uses another core (so 4 cores will be used with this command)
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
                               aCovariate)))))  # count number of individuals having the phenotype and not taking/taking the covariate
    noConditionNegative <- condition[1]
    noConditionPositive <- condition[2]
    robustCombination <- 0
    if ((noControl >= cutoff) &&
        (noConditionNegative >= cutoff) &&
        (noConditionPositive >= cutoff)) {
      robustCombination <- 1
    }
    write(paste
          (aCovariate,
            noControl,
            noConditionNegative,
            noConditionPositive,
            robustCombination,
            sep = "\t"),
          file = "testfile.txt",
          append = TRUE)
    return(data.frame
           (row.names = aCovariate,
             noControl,
             noConditionNegative,
             noConditionPositive,
             robustCombination))

  }

  parallel::stopCluster(cl)
  #return(as.data.frame(parallelReturn,stringsAsFactors=FALSE))
  return(parallelReturn)
}
