#' @importFrom foreach %dopar%
#' @importFrom  lmtest lrtest

CheckReducibility <- function(featureMat,
                              metaMat,
                              noFeatures,
                              noCovariates,
                              features,
                              covariates,
                              Qs,
                              Ds,
                              nnodes,
                              QCutoff,
                              DCutoff,
                              PHS_cutoff,
                              maintenance,
                              verbosity) {

  # load parralel processing environment
  cl <- parallel::makeForkCluster(nnodes = nnodes, outfile = "")
  doParallel::registerDoParallel(cl)
  i <- 0

  r = foreach::foreach(i = seq_along(features), .combine='rbind') %dopar% {

    statusLine <- vector(length = noCovariates, mode = "character")

    # find all covariates which on their end have effect on the feature
    lCovariates <- covariates[which(Qs[i, ] < 0.1)]

    if (verbosity == "debug") {
      write(paste0(as.character(i),
                   as.character(length(lCovariates)),
                   lCovariates),
            file = "lCovariatesIsZero.txt",
            append = TRUE,
            sep = "\t")
    }

    if (length (lCovariates) == 0 ) {

      ##
      ##
      if(verbosity == "debug"){
        write("returned whole NS line",
                file = "lCovariatesIsZero.txt",
                sep = "\t",
                append = TRUE)
      }
      ##
      ##

      statusLine[seq_along(covariates)] <- "NS"
      return(statusLine)
    }

    for (j in seq_along(covariates)) {

      aFeature <- as.character (features [i])
      aCovariate <- as.character (covariates [j])

      status <- "NS"

      if (is.na (Qs [i, j]) |
          Qs [i, j] >= QCutoff |
          abs (Ds [i, j]) <= DCutoff) {
        statusLine[j] <- status

        if(verbosity == "debug"){
          write(paste(aFeature, aCovariate, status, sep = "\t"),
                file = "LRT_pValue.txt",
                append = TRUE)
        }

        next
      }



      # only do post-hoc test for feature that is significant by itself
      # if ((! is.na (Qs [i, j])) &&  # if q exists
      #     Qs [i, j] < QCutoff &&  # q smaller then cutoff
      #     abs (Ds [i, j]) > DCutoff ) {  # effect size bigger than cutoff

        #new approach --> append i-th feature to the metaMat dataframe
      subMerge <- metaMat
      subMerge$FeatureValue <- featureMat [,i]
      subMerge <- na.exclude(subMerge)
      #subMerge <- as.data.frame(na.exclude(subMerge))

        # find all covariates which on their end has effect on the feature

        # test for each of these covariates the forward and reverse formula,
          #count if lax or strict status achieved
        if (length (lCovariates) > 0 &&
            paste0 (lCovariates, collapse = "") != aCovariate) {

          #status <- "STRICTLY DECONFOUNDED"
          status <- "SD"
          # hardest to reach, means no covariate eliminates this signal

          for (anotherCovariate in lCovariates) {

            if (anotherCovariate == aCovariate) { next; }

            # if fail too far, set status to CONFOUNDED
            # otherwise set to LAX

            # significance of post-hoc test:
              #does the tested covariate/predictor add anything beyond what
              #this confounder/covariate does?
            aP_call_forward <-
              paste0 ("lmtest::lrtest (stats::lm (rank (FeatureValue) ~ ",
                      aCovariate,
                      " + ",
                      anotherCovariate,
                      ", data = subMerge), stats::lm (rank (FeatureValue) ~ ",
                      anotherCovariate,
                      ", data = subMerge))$'Pr(>Chisq)' [2]")

            # significance of post-hoc test:
              #does the confounder/covariate add anything beyond the
              #covariate/predictor?
            aP_call_reverse <-
              paste0 ("lmtest::lrtest (stats::lm (rank (FeatureValue) ~ ",
                      aCovariate,
                      " + ",
                      anotherCovariate,
                      ", data = subMerge), stats::lm (rank (FeatureValue) ~ ",
                      aCovariate,
                      ", data = subMerge))$'Pr(>Chisq)' [2]")

            if(verbosity == "debug"){
              write(paste("I have come this far", i, j,  sep = "\t"),
                    file = "LRT_pValue.txt",
                    append = TRUE)
            }

            tryCatch({
                        aP_forward <- eval (parse (text = as.character (aP_call_forward)))
                      },
                      error=function(cond){
                        if(verbosity == "debug"){
                          write(paste("there seems tpo be aproblem here", i, j, cond,  sep = "\t"),
                                file = "LRT_pValue.txt",
                                append = TRUE)
                        }
                        message("HEEEELP")
                        message(cond)
                        return(NA)
                      }
                     )

            tryCatch({
              aP_reverse <- eval (parse (text = as.character (aP_call_reverse)))
            },
            error=function(cond){
              if(verbosity == "debug"){
                write(paste("there seems tpo be aproblem here", i, j, cond,  sep = "\t"),
                      file = "LRT_pValue.txt",
                      append = TRUE)
              }
              message("HEEEELP")
              message(cond)
              return(NA)
            }
            )

            #aP_forward <- try(eval (parse (text = as.character (aP_call_forward))), outFile = "~/IsItHere.txt")
            #aP_reverse <- try(eval (parse (text = as.character (aP_call_reverse))), outFile = "~/IsItHere.txt")

            if (! is.na (aP_forward) &&
                aP_forward >= PHS_cutoff &&
                aP_reverse < PHS_cutoff) {
              #status <- "CONFOUNDED"
              status <- anotherCovariate
              #break
            } # another feature fully explains this

            if (! is.na (aP_forward) &&
                aP_forward >= PHS_cutoff &&
                aP_reverse >= PHS_cutoff) {
              #status <- "LAXLY DECONFOUNDED"
              status <- "LD"
            } # cannot be ruled out another feature explains this

            if(verbosity == "debug"){
              write(paste(aFeature,
                          aCovariate,
                          anotherCovariate,
                          aP_forward,
                          aP_reverse,
                          status,
                          sep = "\t"),
                    file = "LRT_pValue.txt",
                    append = TRUE)
            }

          }

        } else {
          #status <- "NO COVARIATES"
          status <- "NC"# trivially unconfounded because no
                                      #other features play a role
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
