#' @importFrom foreach %dopar%
#' @importFrom  lmtest lrtest
#' @importFrom lme4 glmer
#' @import futile.logger

CheckReducibility <- function(featureMat,
                              metaMat,
                              noFeatures,
                              noCovariates,
                              features,
                              covariates,
                              Qs,
                              Ds,
                              minQValues,
                              nnodes,
                              QCutoff,
                              DCutoff,
                              PHS_cutoff,
                              deconfT,
                              deconfF,
                              doConfs,
                              doRanks,
                              randomVar,
                              RVnames,
                              isRobust,
                              logistic, # new SKF20201017
                              maintenance,
                              verbosity) {

  # load parralel processing environment
  cl <- parallel::makeForkCluster(nnodes = nnodes, outfile = "")
  doParallel::registerDoParallel(cl)
  i <- 0
  isRobust <- isRobust[[2]]
  #isRobust[!isRobust] <- TRUE

  if(verbosity == "debug"){
    write(paste(length(features),
                length(covariates),
                sep = "\t"),
          file = "LRT_pValue.txt",
          append = TRUE)
  }
  r = foreach::foreach(i = seq_along(features), .combine='rbind') %dopar% {


    # if(verbosity == "debug"){
    #   write(paste("I have come this far_35", i,  sep = "\t"),
    #         file = "LRT_pValue.txt",
    #         append = TRUE)
    # }

    statusLine <- vector(length = noCovariates, mode = "character")

    # find all covariates which on their end have effect on the feature
    # add all those that shall always be tested, and remove those that shall never be tested
    lCovariates <- covariates[which(Qs[i, ] < 0.1)]
    lCovariates <- c(lCovariates, deconfT)
    lCovariates <- lCovariates[!(lCovariates %in% deconfF)]
    if (!is.na(RVnames)) {
      lCovariates <- lCovariates[!(lCovariates %in% RVnames)]
    }

    if(verbosity == "debug"){
      write(paste(length(features),
                  length(covariates),
                  length(lCovariates),
                  sep = "\t"),
            file = "LRT_pValue.txt",
            append = TRUE)
    }

    if (verbosity == "debug") {
      write(paste0(as.character(i),
                   as.character(length(lCovariates)),
                   lCovariates),
            file = "lCovariatesIsZero.txt",
            append = TRUE,
            sep = "\t")
    }

    # if(verbosity == "debug"){
    #   write(paste("I have come this far_58", i,  sep = "\t"),
    #         file = "LRT_pValue.txt",
    #         append = TRUE)
    # }


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


      progress <- paste0(round(x = ((i/length(features))*100),
                               digits = 2), "%")
      flog.info(msg = paste("Deconfounding -- processed",
                            progress,
                            "of features."),
                name = "my.logger")

      statusLine[seq_along(covariates)] <- "NS"

      return(statusLine)
    }

    for (j in seq_along(covariates)) {

      if(verbosity == "debug"){
        write(paste("106!",
                    i,
                    j,
                    sep = "\t"),
              file = "LRT_pValue.txt",
              append = TRUE)
      }


      aFeature <- as.character (features [i])
      aCovariate <- as.character (covariates [j])

      status <- "NS"

      if (is.na (Qs [i, j]) |
          Qs [i, j] >= QCutoff |
          is.na(Ds [i, j]) | # this little idiot was missing!
          abs (Ds [i, j]) <= DCutoff) {
        statusLine[j] <- status

        # if(verbosity == "debug"){
        #   write(paste("skipped:", i, j, status, sep = "\t"),
        #         file = "LRT_pValue.txt",
        #         append = TRUE)
        # }

        next
      }




      # only do post-hoc test for feature that is significant by itself
      # if ((! is.na (Qs [i, j])) &&  # if q exists
      #     Qs [i, j] < QCutoff &&  # q smaller then cutoff
      #     abs (Ds [i, j]) > DCutoff ) {  # effect size bigger than cutoff

        #new approach --> append i-th feature to the metaMat dataframe
      subMerge <- metaMat
      subMerge$FeatureValue <- featureMat [,i]
      #subMerge <- na.exclude(subMerge)

      subMerge <- subset (subMerge, ! is.na (FeatureValue))
      subMerge <- eval (parse (text = paste0 ("subset (subMerge, ! is.na (", aCovariate, "))")))

      # rank transfer the metavariables listed in doRanks
      if (!is.na(doRanks)) {
        for (toRank in doRanks) {
          subMerge[toRank] <- rank(subMerge[toRank])
        }
      }

      if (!is.na(randomVar)) {

        # standard behaviour for continuous features
        head <- "lme4::lmer (rank (FeatureValue) ~ "
        tail <- ", data = subMerge, REML = FALSE)"

        if (logistic == TRUE) { # alternative behaviour for binary features
          head <- "lme4::glmer (FeatureValue ~ "
          tail <- ", data = subMerge, family = \"binomial\")"
        }

        mixedmodel1 <- eval(
          parse(
            text = as.character(
              paste0 (
                head,
                aCovariate,
                randomVar,
                tail))))


        mixedmodel2 <- eval(
          parse(
            text = as.character(
              paste0 (
                head,
                substr(randomVar, 2, nchar(randomVar)), # remove the "+" from the string, as randomVar is first argument here
                tail))))


        aP_mixed <- lmtest::lrtest(mixedmodel1, mixedmodel2)$'Pr(>Chisq)' [2]
        if(verbosity == "debug"){
          write(paste("I have come this far_171", i, ", aP_mixed is", aP_mixed,  sep = "\t"),
                file = "LRT_pValue.txt",
                append = TRUE)
        }

        if (aP_mixed >= PHS_cutoff) {

          statusLine[j] <- status
          next
        }

      } # end randomVar

      # if(verbosity == "debug"){
      #   write(paste("I have come this far_162", i,  sep = "\t"),
      #         file = "LRT_pValue.txt",
      #         append = TRUE)
      # }

      #subMerge <- as.data.frame(na.exclude(subMerge))

        # find all covariates which on their end have effect on the feature

        # test for each of these covariates the forward and reverse formula,
          #count if lax or strict status achieved
      confounders <- NULL
        if (length (lCovariates) > 0 &&
            #paste0 (lCovariates, collapse = "") != aCovariate) {
            (paste0 (lCovariates, collapse = "") != aCovariate || length (covariates[which(minQValues[i, ] < 0.1)]) > 1)) {



          #status <- "STRICTLY DECONFOUNDED"
          status <- "SD"
          # hardest to reach, means no covariate eliminates this signal

          # create vector to collect all confounder names for this feature <-> covariate pair
          # put new dataframe of qs in here
          otherCovariates <- lCovariates
          if (!is.null(minQValues[[1]])) {
            #otherCovariates <- covariates[which(minQValues[i, ] < 0.1)]
            otherCovariates <- unique(c(covariates[which(minQValues[i, ] < 0.1)], lCovariates))
          }
          for (anotherCovariate in otherCovariates) {

            subsubMerge <- eval (parse (text = paste0 ("subset (subMerge, ! is.na (", anotherCovariate, "))")))



            if (anotherCovariate == aCovariate) { next; }

            if (!isRobust[aCovariate, anotherCovariate]) {
              #status <- "NS"
              #statusLine[j] <- status

              # if(verbosity == "debug"){
              #   write(paste("skipped:", i, j, anotherCovariate, status, sep = "\t"),
              #         file = "LRT_pValue.txt",
              #         append = TRUE)
              # }

              next
            }

            # if fail too far, set status to CONFOUNDED
            # otherwise set to LAX

            # significance of post-hoc test:
              #does the tested covariate/predictor add anything beyond what
              #this confounder/covariate does?

            # aP_call_forward <-
            #   paste0 ("lmtest::lrtest (stats::lm (rank (FeatureValue) ~ ",
            #           aCovariate,
            #           " + ",
            #           anotherCovariate,
            #           ", data = subMerge), stats::lm (rank (FeatureValue) ~ ",
            #           anotherCovariate,
            #           ", data = subMerge))$'Pr(>Chisq)' [2]")

            # significance of post-hoc test:
              #does the confounder/covariate add anything beyond the
              #covariate/predictor?

            # aP_call_reverse <-
            #   paste0 ("lmtest::lrtest (stats::lm (rank (FeatureValue) ~ ",
            #           aCovariate,
            #           " + ",
            #           anotherCovariate,
            #           ", data = subMerge), stats::lm (rank (FeatureValue) ~ ",
            #           aCovariate,
            #           ", data = subMerge))$'Pr(>Chisq)' [2]")

            # if(verbosity == "debug"){
            #   write(paste("I have come this far_198", i, j,  sep = "\t"),
            #         file = "LRT_pValue.txt",
            #         append = TRUE)
            # }


            modAlg <- "stats::lm (rank (FeatureValue) ~ "
            lastPart <- paste0(", data = subsubMerge)", collapse = "")

            if (logistic == TRUE) { # alternative behaviour for binary features
              modAlg <- "stats::glm (FeatureValue ~ "
              lastPart <- ", data = subsubMerge, family = \"binomial\")"
            }

            if (!is.na(randomVar)) { # switch to lmer and REML = FALSE when randomEffects are included
              modAlg <- "lme4::lmer (rank (FeatureValue) ~ "
              lastPart <- paste0(randomVar, ", data = subsubMerge, REML = FALSE)", collapse = "")

              if (logistic == TRUE) { # alternative behaviour for binary features
                modAlg <- "lme4::glmer (FeatureValue ~ "
                lastPart <- ", data = subsubMerge, family = \"binomial\")"
              }
            }
            # compute the three needed linear models
            lmBoth <- eval(
              parse(
                text = as.character(
                  paste0 (
                    modAlg,
                    aCovariate,
                    " + ",
                    anotherCovariate,
                    lastPart))))

            lmA <- eval(
              parse(
                text = as.character(
                  paste0 (
                    modAlg,
                    aCovariate,
                    lastPart))))

            lmAnother <- eval(
              parse(
                text = as.character(
                  paste0 (
                    modAlg,
                    anotherCovariate,
                    lastPart))))

            # make the two needed likelihood ratio tests to determine wether ...
              # one covariate has influence on feature ...
              # beyond that of the other covariate

            aP_forward <- lmtest::lrtest (lmBoth, lmAnother)$'Pr(>Chisq)' [2]
            aP_reverse <- lmtest::lrtest (lmBoth, lmA)$'Pr(>Chisq)' [2]

            # if(verbosity == "debug"){
            #   write(paste("I have come this far_301", i, j, anotherCovariate,  sep = "\t"),
            #         file = "LRT_pValue.txt",
            #         append = TRUE)
            # }

            # additonal control of confidence intervals for the covariates within the linear models
            conf_aCovariate <- TRUE
            conf_anotherCovariate <- TRUE
            if (doConfs > 0) { # doConfs = 1 --> just logging
              if (is.numeric(subsubMerge$aCovariate) && is.numeric(subsubMerge$anotherCovariate)) { # categorical variables are excluded for easier processing
                confints <- confint(lmBoth)
                conf_aCovariate <- !(sign(confints[1, 1]) == sign(confints[1, 2])) # signs are different, if confint is spanning 0
                conf_anotherCovariate  <- !(sign(confints[2, 1]) == sign(confints[2, 2]))

                if ((aP_forward < PHS_cutoff) && conf_aCovariate) { # if forward test is significant, but aCovariate confint spans 0
                  flog.warn(msg = paste('lrt: ',
                                        features[i],
                                        aCovariate,
                                        anotherCovariate,
                                        '-- forward linear model is < PHS_cutoff, but confidence intervall for',
                                        aCovariate, 'is spanning 0.'),
                            name = "my.logger")
                  if (doConfs > 1) { # if doConfs ==2 make lrt non-significant
                    aP_forward <- 1
                  }
                }

                if ((aP_reverse < PHS_cutoff) && conf_anotherCovariate) { # if reverse test is significant, but anotherCovariate confint spans 0
                  flog.warn(msg = paste('lrt: ',
                                        features[i],
                                        aCovariate,
                                        anotherCovariate,
                                        '-- reverse linear model is < PHS_cutoff, but confidence intervall for',
                                        anotherCovariate, 'is spanning 0.'),
                            name = "my.logger")
                  if (doConfs > 1) { # if doConfs ==2 make lrt non-significant
                    aP_reverse <- 1
                  }
                }

              }
            }


            # if(verbosity == "debug"){
            #   write(paste("I have come this far 330", i, j,  sep = "\t"),
            #         file = "LRT_pValue.txt",
            #         append = TRUE)
            # }

            # tryCatch({
            #             aP_forward <- eval (
            #               parse (
            #                 text = as.character (aP_call_forward)))
            #           },
            #           error=function(cond){
            #             if(verbosity == "debug"){
            #               write(paste("there seems tpo be aproblem here",
            #                           i,
            #                           j,
            #                           cond,
            #                           sep = "\t"),
            #                     file = "LRT_pValue.txt",
            #                     append = TRUE)
            #             }
            #             message("HEEEELP")
            #             message(cond)
            #             return(NA)
            #           }
            #          )
            #
            # tryCatch({
            #   aP_reverse <- eval (parse (text = as.character (aP_call_reverse)))
            # },
            # error=function(cond){
            #   if(verbosity == "debug"){
            #     write(paste("there seems tpo be aproblem here",
            #                 i,
            #                 j,
            #                 cond,
            #                 sep = "\t"),
            #           file = "LRT_pValue.txt",
            #           append = TRUE)
            #   }
            #   message("HEEEELP")
            #   message(cond)
            #   return(NA)
            # }
            # )

            #aP_forward <- try(eval (parse (text = as.character (aP_call_forward))), outFile = "~/IsItHere.txt")
            #aP_reverse <- try(eval (parse (text = as.character (aP_call_reverse))), outFile = "~/IsItHere.txt")

            if (! is.na (aP_forward) &&
                aP_forward >= PHS_cutoff &&
                aP_reverse < PHS_cutoff) {
              #status <- "CONFOUNDED"
              status <- anotherCovariate

              # all detected confounders are added to this vector
              confounders <- c(confounders, anotherCovariate)

              # if(verbosity == "debug"){
              #   write(paste("Eureka!",
              #               i,
              #               j,
              #               status,
              #               sep = "\t"),
              #         file = "LRT_pValue.txt",
              #         append = TRUE)
              # }

              #break
            } # another feature fully explains this

            if (! is.na (aP_forward) &&
                aP_forward >= PHS_cutoff &&
                aP_reverse >= PHS_cutoff) {
              #status <- "LAXLY DECONFOUNDED"
              status <- "LD"
            } # cannot be ruled out another feature explains this
#
#             if(verbosity == "debug"){
#               write(paste(aFeature,
#                           aCovariate,
#                           anotherCovariate,
#                           aP_forward,
#                           aP_reverse,
#                           status,
#                           sep = "\t"),
#                     file = "LRT_pValue.txt",
#                     append = TRUE)
#             }

            if(verbosity == "debug"){
              write(paste("436!",
                          i,
                          j,
                          anotherCovariate,
                          sep = "\t"),
                    file = "LRT_pValue.txt",
                    append = TRUE)
            }
          } #end "for (anotherCovariate in lCovariates) {"


          # if(verbosity == "debug"){
          #   write(paste("Karamba!",
          #               i,
          #               j,
          #               status,
          #               sep = "\t"),
          #         file = "LRT_pValue.txt",
          #         append = TRUE)
          # }

        } else {
          #status <- "NO COVARIATES"
          status <- "NC"# trivially unconfounded because no
                                      #other features play a role
          # if(verbosity == "debug"){
          #   write(paste("EASY!",
          #               i,
          #               j,
          #               status,
          #               sep = "\t"),
          #         file = "LRT_pValue.txt",
          #         append = TRUE)
          # }
        }

      # if confounders where detected, they will all be printed as a single string as status
      if (!is.null(confounders)) {
        status <- paste(confounders, collapse = ", ")
      }
      # }# status assignment for a sinlge feature <-> covariate combination


      statusLine[j] <- status
      #Ss [as.character (aFeature), as.character (aCovariate)] <- status
      # if(verbosity == "debug"){
      #   write(paste("468!",
      #               i,
      #               j,
      #               status,
      #               sep = "\t"),
      #         file = "LRT_pValue.txt",
      #         append = TRUE)
      # }


    }# inner nested for j loop


     if ((i %% 10) == 0) {
    progress <- paste0(round(x = ((i/length(features))*100),
                             digits = 2), "%")
    flog.info(msg = paste("Deconfounding -- processed",
                          progress,
                          "of features."),
              name = "my.logger")
     }

    if(verbosity == "debug"){
      write("one more line done!",
            file = "LRT_pValue.txt",
            append = TRUE)
    }

    statusLine



  }# foreach loop
  parallel::stopCluster(cl)

  if(verbosity == "debug"){
    write("Everything done in checkReducibility!",
          file = "LRT_pValue.txt",
          append = TRUE)
  }

  flog.info(msg = paste("Deconfounding -- processed 100% of features."),
            name = "my.logger")

  rownames(r) <- features
  colnames(r) <- covariates
  r
}# function body
