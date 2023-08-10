#' @import foreach
#' @import  lmtest
#' @import lme4
#' @import futile.logger
#' @import detectseparation


CheckReducibility_linear <- function(featureMat,
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
                              fixedVar, #new TB20230727
                              RVnames,
                              isRobust,
                              logistic, # new SKF20201017
                              rawCounts, # new TB20220202
                              maintenance,
                              verbosity,
                              nAGQ,
                              collectMods # new TB20220208
                              ) {


  if (collectMods) {
    collectedMods <- list()
  }
  # safe_lm <- function(lmText) {
  #   out <- tryCatch(
  #     {
  #       eval(parse(text = as.character(lmText)))
  #     }, error = function(cond) {
  #       print(cond)
  #       return(NA)
  #     }
  #   )
  #   return(out)
  # }
  #
  # safe_lrtest <- function(lm1, lm2) {
  #   out <- tryCatch(
  #     {
  #       lmtest::lrtest(lm1, lm2)$'Pr(>Chisq)' [2]
  #     }, error = function(cond) {
  #       print(cond)
  #       return(NA)
  #     }
  #   )
  #   return(out)
  # }

  # removed tryCatch from safe functions
  safe_lm <- function(lmText) {
    out <- eval(parse(text = as.character(lmText)))
    return(out)
  }

  safe_lrtest <- function(lm1, lm2) {
    out <- lmtest::lrtest(lm1, lm2)$'Pr(>Chisq)' [2]
    return(out)
  }

  #new TB20220202
  if (rawCounts == TRUE) {
    # compute totReadCount per sample and append to metaMat
    #print(head(featureMat))
    totReadCount <- as.data.frame(rowSums(featureMat, na.rm = T))
    #print(head(totReadCount))
    metaMat <- merge(metaMat, totReadCount, by = 0, sort = FALSE)
    #print(head(metaMat))
    colnames(metaMat)[ncol(metaMat)] <- "totReadCount"
    row.names(metaMat) <- metaMat$Row.names
    metaMat$Row.names <- NULL
    #print(head(metaMat))
  }

  if (!is.na(randomVar[[1]])) {
    randomVarLine <- paste0("+ (1|", randomVar, ")", collapse = ' ')
    if (!is.na(fixedVar[[1]])) {
      fixedVarLine <- paste0("+ ", fixedVar, collapse = ' ')
      randomVarLine <- paste0(randomVarLine, " ", fixedVarLine, collapse = ' ')
    }
  }

  if (is.na(randomVar[[1]]) & !is.na(fixedVar[[1]])) {
    randomVarLine <- paste0("+ ", fixedVar, collapse = ' ')
  }

  # load parralel processing environment
  if (.Platform$OS.type == "unix") {
    # unix
    cl <- parallel::makeForkCluster(nnodes = nnodes, outfile = "")
    doParallel::registerDoParallel(cl)
  } else {
    # windows
    cl <- snow::makeCluster(nnodes, type = "SOCK", outfile = "")
    doSNOW::registerDoSNOW(cl)
  }
  i <- 0
  isRobust <- isRobust[[2]]
  #isRobust[!isRobust] <- TRUE

  # if(verbosity == "debug"){
  #   write(paste(length(features),
  #               length(covariates),
  #               sep = "\t"),
  #         file = "LRT_pValue.txt",
  #         append = TRUE)
  # }

  r = foreach::foreach(i = seq_along(features), .combine='rbind') %do% {

    if (collectMods) {
      collectedMods[[features[i]]] <- list()
    }

    statusLine <- vector(length = noCovariates, mode = "character")
    # find all covariates which on their end have effect on the feature
    # add all those that shall always be tested, and remove those that shall never be tested
    lCovariates <- covariates[which(Qs[i, ] < 0.1)]
    lCovariates <- c(lCovariates, deconfT)
    lCovariates <- lCovariates[!(lCovariates %in% deconfF)]
    # remove names of random variables from this list
    if (!is.na(RVnames[[1]])) {
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



    if (length (lCovariates) == 0 ) {

      ##
      ##
      if (verbosity == "debug") {
        write(
          "returned whole NS line",
          file = "lCovariatesIsZero.txt",
          sep = "\t",
          append = TRUE
        )
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
      # if(verbosity == "debug"){
      #   write(paste("106!",
      #               i,
      #               j,
      #               sep = "\t"),
      #         file = "LRT_pValue.txt",
      #         append = TRUE)
      #
      # }

      aFeature <- as.character (features [i])
      aCovariate <- as.character (covariates [j])

      status <- "NS"
      if (collectMods) {
        collectedMods[[aFeature]][[aCovariate]] <- list()
      }

      if (is.na (Qs [i, j]) |
          Qs [i, j] >= QCutoff |
          is.na(Ds [i, j]) |
          abs (Ds [i, j]) <= DCutoff) {

        if (aCovariate %in% RVnames) {
          # set label to NA for all random vars
          statusLine[j] <- NA
          next
        }

        statusLine[j] <- status
        next
      }

      subMerge <- metaMat
      subMerge$FeatureValue <- featureMat [,i]

      #remove NAs only in feature and acovariate column
      subMerge <- subset (subMerge, ! is.na (FeatureValue))
      subMerge <- eval (parse (text = paste0 ("subset (subMerge, ! is.na (", aCovariate, "))")))

      # test for complete separation in model using only feature and covariate
      if (logistic == T & !is.na(randomVar[[1]])) {
        print("testin separation without pot conf")
        #test for separation in model without random part
        glmmodeltext <-          paste0 ("stats::glm (FeatureValue ~ ",
                                         aCovariate,
                                         ", data = subMerge, family = \"binomial\", method = \"detect_separation\")")
        isSeperated <- eval(parse(text = as.character(glmmodeltext)))
        if (isSeperated$outcome) {
          flog.warn(msg = paste("Separation for:",
                                aFeature,
                                "and",
                                aCovariate),
                    name = "my.logger")
          statusLine[j] <- "AD"
          next
        }
      }

      # rank transfer the metavariables listed in doRanks
      if (!is.na(doRanks[[1]])) {
        for (toRank in doRanks) {
          subMerge[toRank] <- rank(subMerge[toRank])
        }
      }

      if (!(is.na(randomVar[[1]]) & !(aCovariate %in% fixedVar)) |
          !(is.na(fixedVar[[1]]) & !(aCovariate %in% fixedVar))) {
        # if either random or fixed Var is supplied and aCovariate != fixedVar, do a naive test
          # whether association is reducible to these random/fixed Vars

        head <- "stats::lm (rank (FeatureValue) ~ "
        tail <- paste0(", data = subMerge)", collapse = "")

        if (!is.na(randomVar[[1]])) {
          head <- "lme4::lmer (rank (FeatureValue) ~ "
          tail <- ", data = subMerge, REML = FALSE)"
        }

        if (logistic == TRUE) { # alternative behaviour for binary features
          head <- "lme4::glmer (FeatureValue ~ "
          tail <- paste0(", data = subMerge, family = \"binomial\", nAGQ = ",
                         nAGQ,
                         ")",
                         collapse = "")
        }
        if (rawCounts == TRUE) { # alternative behavior for not rarefied abundances
          head <- "lme4::glmer (cbind(FeatureValue, totReadCount) ~ "
          tail <- paste0(", data = subMerge, family = \"binomial\", nAGQ = ",
                         nAGQ,
                         ")",
                         collapse = "")
        }

        if(verbosity == "debug"){
          write(paste("LRT_randOnly_for", aFeature, aCovariate,  sep = "\t"),
                file = "LRT_pValue.txt",
                append = TRUE)
        }

        mixedmodel1Text <- paste0 (head,
                                   aCovariate,
                                   randomVarLine,
                                   tail)

        mixedmodel1 <- safe_lm(mixedmodel1Text)

        mixedmodel2Text <- paste0 (head,
                                   substr(randomVarLine,
                                          2,
                                          nchar(randomVarLine)),
                                   # remove the "+" from the string, as randomVar is first argument here
                                   tail)

        mixedmodel2 <- safe_lm(mixedmodel2Text)


        aP_mixed <- safe_lrtest(mixedmodel1, mixedmodel2)
        if(verbosity == "debug"){
          write(paste("LRT_randOnly_for", aFeature, aCovariate, aP_mixed,  sep = "\t"),
                file = "LRT_pValue.txt",
                append = TRUE)
        }

        if (!is.na(aP_mixed) && aP_mixed >= PHS_cutoff) {

          statusLine[j] <- status
          next
        }

      } # end randomVar



      # find all covariates which on their end have effect on the feature

      # test for each of these covariates the forward and reverse formula,
      #count if lax or strict status achieved
      confounders <- NULL
      if (length (lCovariates) > 0 &&
          (paste0 (lCovariates, collapse = "") != aCovariate ||
           length (covariates[which(minQValues[i, ] < 0.1)]) > 1)) {

        #status <- "STRICTLY DECONFOUNDED"
        status <- "OK_sd"
        # hardest to reach, means no covariate eliminates this signal

        # create vector to collect all confounder names for this feature <-> covariate pair
        # put new dataframe of qs in here
        otherCovariates <- lCovariates
        if (!is.null(minQValues[[1]])) {
          otherCovariates <- unique(c(covariates[which(minQValues[i, ] < 0.1)], lCovariates))
        }
        for (anotherCovariate in otherCovariates) {

          if ((anotherCovariate == aCovariate) ||
              (!isRobust[aCovariate, anotherCovariate])) {
            next
          }

          if (collectMods) {
            collectedMods[[aFeature]][[aCovariate]][[anotherCovariate]] <- list()
          }

          # remove rows where anotherCovariate has NAs
          subsubMerge <- eval (parse (text = paste0 ("subset (subMerge, ! is.na (", anotherCovariate, "))")))

          if (logistic == T & !is.na(randomVar[[1]])) {
            print("testin separation")
            #test for separation in model without random part
            glmmodeltext <-          paste0 ("glm (FeatureValue ~ ",
                                             aCovariate,
                                             " + ",
                                             anotherCovariate,
                                             ", data = subMerge, family = \"binomial\")")
            glmmodel <- eval(parse(text = as.character(glmmodeltext)))
            isSeperated <- update(glmmodel, method="detect_separation")
            if (isSeperated$outcome) {
              flog.warn(msg = paste("Separation for:",
                                    aFeature,
                                    ",",
                                    aCovariate,
                                    "and",
                                    anotherCovariate),
                        name = "my.logger")
              separation <- T
              status <- NA
              if (collectMods) {
                collectedMods[[aFeature]][[aCovariate]][[anotherCovariate]][["full"]] <- glmmodel
              }
              next
            }
          }


          modAlg <- "stats::lm (rank (FeatureValue) ~ "
          lastPart <- paste0(", data = subsubMerge)", collapse = "")

          if (logistic == TRUE) { # alternative behavior for binary features
            modAlg <- "stats::glm (FeatureValue ~ "
            lastPart <- ", data = subsubMerge, family = \"binomial\")"
          }
          if (rawCounts == TRUE) { # alternative behavior for not rarefied abundances
            modAlg <- "stats::glm (cbind(FeatureValue, totReadCount) ~ "
            lastPart <- ", data = subsubMerge, family = \"binomial\")"
          }

          if (!is.na(fixedVar[[1]])) { # prefix the lastPart with fixedEffect names
            lastPart <- paste0(randomVarLine, lastPart, collapse = "")
          }

          if (!is.na(randomVar[[1]])) { # switch to lmer and REML = FALSE when randomEffects are included
            modAlg <- "lme4::lmer (rank (FeatureValue) ~ "
            lastPart <- paste0(randomVarLine,
                               ", data = subsubMerge, REML = FALSE)",
                               collapse = "")

            if (logistic == TRUE) { # alternative behavior for binary features
              modAlg <- "lme4::glmer (FeatureValue ~ "
              lastPart <- paste0(randomVarLine,
                                 ", data = subsubMerge, family = \"binomial\", nAGQ = ",
                                 nAGQ,
                                 ")",
                                 collapse = "")

            }
            if (rawCounts == TRUE) { # alternative behavior for not rarefied abundances
              modAlg <- "lme4::glmer (cbind(FeatureValue, totReadCount) ~ "
              lastPart <- paste0(randomVarLine,
                                 ", data = subsubMerge, family = \"binomial\", nAGQ = ",
                                 nAGQ,
                                 ")",
                                 collapse = "")
            }
          }

          if(verbosity == "debug"){
            write(paste("LRTsFwdRvsFor",
                        aFeature,
                        aCovariate,
                        anotherCovariate,
                        "\n",
                        sep = "\t"),
                  file = "LRT_pValue.txt",
                  append = TRUE)
          }

          # compute the three needed linear models
          lmBothText <- paste0 (modAlg,
                                aCovariate,
                                " + ",
                                anotherCovariate,
                                lastPart)
          lmBoth <- safe_lm(lmBothText)

          lmAText <- paste0 (modAlg,
                             aCovariate,
                             lastPart)
          lmA <- safe_lm(lmAText)

          lmAnotherText <- paste0 (modAlg,
                                   anotherCovariate,
                                   lastPart)
          lmAnother <- safe_lm(lmAnotherText)

          # collect the fitted model objects
          if (collectMods) {
            collectedMods[[aFeature]][[aCovariate]][[anotherCovariate]][["full"]] <- lmBoth
            collectedMods[[aFeature]][[aCovariate]][[anotherCovariate]][["cov"]] <- lmA
            collectedMods[[aFeature]][[aCovariate]][[anotherCovariate]][["conf"]] <- lmAnother
          }
          # make the two needed likelihood ratio tests to determine wether ...
          # one covariate has influence on feature ...
          # beyond that of the other covariate


          if (is(lmA, "logical") ||
              is(lmAnother, "logical") ||
              is(lmBoth, "logical") ) {
            # class is logical if lmX == NA

            if(verbosity == "debug"){
              write(paste("LRTsFwdRvsFor",
                          aFeature,
                          aCovariate,
                          anotherCovariate,
                          "found faulty models!\n",
                          lmAText,
                          "\n",
                          summary(lmA),
                          "\n",
                          lmAnotherText,
                          "\n",
                          summary(lmAnother),
                          "\n",
                          lmBothText,
                          "\n",
                          summary(lmBoth),
                          "\n",
                          sep = "\t"),
                    file = "LRT_pValue.txt",
                    append = TRUE)

            }

            #status <- paste0("AD_", anotherCovariate)
            status <- "AD"
            print("This should never be executed!")
            next
          }


          # get_aP_forward <- purrr::possibly(lmtest::lrtest(lmBoth, lmAnother)$'Pr(>Chisq)' [2],
          #                                 otherwise = NA)
          #
          # aP_forward <- get_aP_forward()
          aP_forward <- safe_lrtest(lmBoth, lmAnother)


          # get_aP_reverse <- purrr::possibly(lmtest::lrtest(lmBoth, lmA)$'Pr(>Chisq)' [2],
          #                                   otherwise = NA)
          # aP_reverse <- get_aP_reverse()
          aP_reverse <- safe_lrtest(lmBoth, lmA)


          # additonal control of confidence intervals for the covariates within the linear models
          conf_aCovariate <- TRUE
          conf_anotherCovariate <- TRUE
          if (doConfs >= 0 && !is(lmBoth, "logical")) { # doConfs = 1 --> just logging
            print("computing Confs")
            if (is.numeric(subsubMerge$aCovariate) && is.numeric(subsubMerge$anotherCovariate)) { # categorical variables are excluded for easier processing
              confints <- confint(lmBoth)
              conf_aCovariate <- !(sign(confints[1, 1]) == sign(confints[1, 2])) # signs are different, if confint is spanning 0
              conf_anotherCovariate  <- !(sign(confints[2, 1]) == sign(confints[2, 2]))

              if (! is.na (aP_forward) && (aP_forward < PHS_cutoff) && conf_aCovariate) { # if forward test is significant, but aCovariate confint spans 0
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
                status <- "OK_d"
              }

              if (! is.na (aP_forward) && (aP_reverse < PHS_cutoff) && conf_anotherCovariate) { # if reverse test is significant, but anotherCovariate confint spans 0
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

          if (! is.na (aP_forward) &&
              ! is.na (aP_reverse) &&
              aP_forward >= PHS_cutoff &&
              aP_reverse < PHS_cutoff) {
            #status <- "CONFOUNDED"
            status <- anotherCovariate

            # all detected confounders are added to this vector
            confounders <- c(confounders, anotherCovariate)

          } # another feature fully explains this

          if (! is.na (aP_forward) &&
              ! is.na (aP_reverse) &&
              aP_forward >= PHS_cutoff &&
              aP_reverse >= PHS_cutoff) {
            #status <- "LAXLY DECONFOUNDED"
            status <- "AD"
          } # cannot be ruled out another feature explains this

          if(verbosity == "debug"){
            write(paste("LRTsFwdRvsFor",
                        aFeature,
                        aCovariate,
                        anotherCovariate,
                        aP_forward,
                        aP_reverse,
                        sep = "\t"),
                  file = "LRT_pValue.txt",
                  append = TRUE)
          }
        } #end "for (anotherCovariate in lCovariates) {"


      } else {
        #status <- "NO COVARIATES"
        status <- "OK_nc"# trivially unconfounded because no
        #other features play a role
      }

      # if confounders where detected, they will all be printed as a single string as status
      if (!is.null(confounders)) {
        status <- paste0("C: ", paste(confounders, collapse = ", "))
      }
      # }# status assignment for a sinlge feature <-> covariate combination

      statusLine[j] <- status
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

  # close parallel processing environment
  if (.Platform$OS.type == "unix") {
    parallel::stopCluster(cl) # unix
  } else {
    snow::stopCluster(cl) # windows
  }

  if(verbosity == "debug"){
    write("Everything done in checkReducibility!",
          file = "LRT_pValue.txt",
          append = TRUE)
  }

  flog.info(msg = paste("Deconfounding -- processed 100% of features."),
            name = "my.logger")

  rownames(r) <- features
  colnames(r) <- covariates

  if (collectMods) {
    r <- list(r, collectedMods)
  }
  r
}# function body
