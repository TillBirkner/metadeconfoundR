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

CheckIndividualSignificance <- function(featureMat, metaMat, nnodes=1) {
  cl<-makeForkCluster(nnodes = nnodes, outfile="")  # the parent process uses another core (so 4 cores will be used with this command)
  registerDoParallel(cl)
  parallelReturn = foreach(i= 1:n, .combine = 'rbind') %dopar% {
    return(i)
  }
  stopCluster(cl)
}
