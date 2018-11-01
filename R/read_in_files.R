# Hello, world!
#
# This is an example function named 'hello'
# which prints 'Hello, world!'.
#
# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Build and Reload Package:  'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'


#' @export checkSufficientPower
checkSufficientPower <- function(featureFile, metaFile, nnodes=1) {

  metaMat <- read.table(file = metaFile, header = T, sep = "\t", row.names = 1)
  covariates <- colnames (metaMat) [4:27] # each covariate + the status category, specific to example
  noCovariates <- length (covariates)

  conditionMat <- metaMat[metaMat$Case_status == "Schizo", ]
  noCondition <- nrow(conditionMat)
  noControl <- nrow(metaMat) - noCondition

  write(paste("covariate", "control", "condition_negative", "condition_positive", sep = "\t"), file = "testfile.txt")
  cl <- parallel::makeForkCluster(nnodes = nnodes, outfile="")  # the parent process uses another core (so 4 cores will be used with this command)
  doParallel::registerDoParallel(cl)

  parallelReturn <- foreach::foreach(i= 2:noCovariates, .combine = 'rbind') %dopar% {

    aCovariate <- as.character (covariates [i])
    condition <- table(eval(parse( text = as.character(paste0("schizoMat$", aCovariate)))))
    noConditionNegative <- condition[1]
    noConditionPositive <- condition[2]
    write(paste(aCovariate, noControl, noConditionNegative, noConditionPositive, sep = "\t"), file = "testfile.txt", append = TRUE)
    return(list(noquote(aCovariate), noControl, noConditionNegative, noConditionPositive))

  }

  parallel::stopCluster(cl)
  #return(as.data.frame(parallelReturn,stringsAsFactors=FALSE))
  return(parallelReturn)
}
