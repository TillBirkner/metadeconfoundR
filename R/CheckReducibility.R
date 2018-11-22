CheckReducibility <- function(noFeatures,
                              noCovariates,
                              features,
                              covariates,
                              Qs,
                              Ds,
                              nnodes,
                              QCutoff,
                              DCutoff,
                              maintenance,
                              verbosity) {

  # load parralel processing environment
  cl <- parallel::makeForkCluster(nnodes = nnodes, outfile = "")  # the parent process uses another core
  doParallel::registerDoParallel(cl)
  i <- 0

  r = foreach::foreach(i= 1:noFeatures, .combine='rbind') %dopar% {

    statusLine <- vector(length = noCovariates, mode = "character")

    # find all covariates which on their end has effect on the feature
    lCovariates <- covariates[which(Qs[i, ] < 0.1)]

    if (verbosity == "debug") {
      write(paste0(as.character(i), as.character(length(lCovariates)), lCovariates), file = "lCovariatesIsZero.txt", append = TRUE, sep = "\t")
    }

    if (length (lCovariates) == 0 ) {
      statusLine[1:noCovariates] <- "NS"
      write("returned whole NS line", file = "lCovariatesIsZero.txt", append = TRUE)
      return(statusLine)
    }

    if (length (lCovariates) == 0 ) {
      write("This should never be printed!!!", file = "lCovariatesIsZero.txt", append = TRUE)
    }

    for (j in 1:noCovariates) {

      status <- "NS"

      if (is.na (Qs [i, j]) |
          Qs [i, j] >= QCutoff |
          abs (Ds [i, j]) <= DCutoff) {

        statusLine[j] <- status
        next
      }

      aFeature <- as.character (features [i])
      aCovariate <- as.character (covariates [j])

      # only do post-hoc test for feature that is significant by itself
      # if ((! is.na (Qs [i, j])) &&  # if q exists
      #     Qs [i, j] < QCutoff &&  # q smaller then cutoff
      #     abs (Ds [i, j]) > DCutoff ) {  # effect size bigger than cutoff

        subFeatures <- featureMat [,i]
        subMerge <- md
        subMerge$FeatureValue <- subFeatures

        # find all covariates which on their end has effect on the feature

        # test for each of these covariates the forward and reverse formula, count if lax or strict status achieved
        if (length (lCovariates) > 0 &&
            paste0 (as.character (lCovariates), collapse = "") != as.character (aCovariate)) {

          status <- "STRICTLY DECONFOUNDED" # hardest to reach, means no covariate eliminates this signal

          for (anotherCovariate in lCovariates) {

            if (anotherCovariate == aCovariate) { next; }

            # if fail too far, set status to CONFOUNDED
            # otherwise set to LAX

            # significance of post-hoc test: does the tested covariate/predictor add anything beyond what this confounder/covariate does?
            aP_call_forward <- paste0 ("lmtest::lrtest (stats::lm (rank (FeatureValue) ~ ",
                                       aCovariate,
                                       " + ",
                                       anotherCovariate,
                                       ", data = subMerge), stats::lm (rank (FeatureValue) ~ ",
                                       anotherCovariate,
                                       ", data = subMerge))$'Pr(>Chisq)' [2]")

            # significance of post-hoc test: does the confounder/covariate add anything beyond the covariate/predictor?
            aP_call_reverse <- paste0 ("lmtest::lrtest (stats::lm (rank (FeatureValue) ~ ",
                                       aCovariate,
                                       " + ",
                                       anotherCovariate,
                                       ", data = subMerge), stats::lm (rank (FeatureValue) ~ ",
                                       aCovariate,
                                       ", data = subMerge))$'Pr(>Chisq)' [2]")

            aP_forward <- eval (parse (text = as.character (aP_call_forward)))
            aP_reverse <- eval (parse (text = as.character (aP_call_reverse)))

            if (! is.na (aP_forward) && aP_forward >= 0.05 && aP_reverse < 0.05) {
              status <- "CONFOUNDED"; break;
            } # another feature fully explains this

            if (! is.na (aP_forward) && aP_forward >= 0.05 && aP_reverse >= 0.05) {
              status <- "LAXLY DECONFOUNDED";
            } # cannot be ruled out another feature explains this

          }

        } else {
          status <- "NO COVARIATES" # trivially unconfounded because no other features plays a role
        }

      # }# status assignment for a sinlge feature <-> covariate combination
      statusLine[j] <- status
      #Ss [as.character (aFeature), as.character (aCovariate)] <- status


    }# inner nested for j loop
    statusLine
  }# foreach loop
  parallel::stopCluster(cl)
  rownames(r) <- features
  colnames(r) <- covariates
  r
}# function body
