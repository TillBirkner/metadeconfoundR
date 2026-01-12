#' @import foreach
#' @import lmtest
#' @import lme4
#' @import logger
#' @importFrom methods is
#' @import detectseparation
#' @importFrom stats update confint
#' @importFrom stringr str_split

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
                              fixedVar, #new TB20230727
                              RVnames,
                              isRobust,
                              logistic, # new SKF20201017
                              rawCounts, # new TB20220202
                              maintenance,
                              nAGQ,
                              collectMods, # new TB20220208
                              noConfConfs # new TB20250827
                              ) {


  if (collectMods) {
    collectedMods <- list()
  }

  # removed tryCatch from safe functions
  safe_lm <- function(lmText) {
    #out <- eval(parse(text = as.character(lmText)))
    out <- tryCatch(eval(parse(text = as.character(lmText))),
                    error = function(cond){
      logger::log_warn(namespace = "metadeconfoundR", paste0("(g)lm(er) failed with error meessage: ", cond$message))
      NA
    })
    return(out)
  }

  safe_lrtest <- function(lm1, lm2) {
    #out <- lmtest::lrtest(lm1, lm2)$'Pr(>Chisq)' [2]

    out <- tryCatch(lmtest::lrtest(lm1, lm2)$'Pr(>Chisq)' [2],
                    error = function(cond){
      logger::log_warn(namespace = "metadeconfoundR", paste0("model lrt failed with error meessage: ", cond$message))
      NA
    })
    return(out)
  }

  #new TB20220202
  if (rawCounts == TRUE) {
    # compute totReadCount per sample and append to metaMat
    totReadCount <- as.data.frame(rowSums(featureMat, na.rm = T))
    metaMat <- merge(metaMat, totReadCount, by = 0, sort = FALSE)
    colnames(metaMat)[ncol(metaMat)] <- "totReadCount"
    row.names(metaMat) <- metaMat$Row.names
    metaMat$Row.names <- NULL
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


  # compute steps for progress log.info #TB20240229
  # aim for no more than 100 steps
  progressSteps <- round(noFeatures/100)
  if (progressSteps == 0) {
    progressSteps <- 1
  }
  if (noFeatures/progressSteps > 100) {
    progressSteps <- progressSteps + 1
  }


  #new TB20240926
  `%toggleDoPar%` <- `%do%`
  if (nnodes != 1) {
    `%toggleDoPar%` <- `%dopar%`
  }

  # load parallel processing environment
  if (.Platform$OS.type == "unix") {
    # unix
    cl <- parallel::makeForkCluster(nnodes = nnodes, outfile = "")
    doParallel::registerDoParallel(cl)
  } else {
    # windows
    cl <- parallel::makeCluster(nnodes, type = "PSOCK", outfile = "")
    doParallel::registerDoParallel(cl)
  }
  i <- 0
  isRobust <- isRobust[[2]]
  #isRobust[!isRobust] <- TRUE

  r = foreach::foreach(i = seq_along(features), .combine='rbind') %toggleDoPar% {

    if (collectMods) {
      collectedMods[[features[i]]] <- list()
    }

    statusLine <- vector(length = noCovariates, mode = "character")
    # find all covariates which on their end have effect on the feature
    # add all those that shall always be tested, and remove those that shall never be tested
    lCovariates <- covariates[which((Qs[i, ] < QCutoff) & (abs(Ds[i, ]) > DCutoff))]
    lCovariates <- c(lCovariates, deconfT)
    lCovariates <- lCovariates[!(lCovariates %in% deconfF)]
    # remove names of random variables from this list
    if (!is.na(RVnames[[1]])) {
      lCovariates <- lCovariates[!(lCovariates %in% RVnames)]
    }



    logger::log_debug(namespace = "metadeconfoundR", paste(
      as.character(i),
      length(features),
      length(covariates),
      length(lCovariates),
      paste(lCovariates, collapse = ", "),
      sep = "\t"
      ))

    if (length (lCovariates) == 0 ) {

      logger::log_debug(namespace = "metadeconfoundR", "returned whole NS line", name = "my.logger")

      if ((i %% progressSteps) == 0) {#TB20240229
        progress <- paste0(round(x = ((i/length(features))*100),
                                 digits = 2), "%")
        logger::log_info(namespace = "metadeconfoundR", paste("Deconfounding -- processed",
                              progress,
                              "of features."))
      }

      statusLine[seq_along(covariates)] <- "NS"

      return(statusLine)
    }

    aFeature <- as.character (features [i])
    for (j in seq_along(covariates)) {
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
          # randomVars will have NA Qs, so are always caught here
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
      subMerge <- subMerge[!is.na(subMerge$FeatureValue), ]
      subMerge <- eval (parse (text = paste0 ("subset (subMerge, ! is.na (", aCovariate, "))")))

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

        logger::log_debug(namespace = "metadeconfoundR", paste("LRT_randOnly_for", aFeature, aCovariate, sep = "\t"))

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

        logger::log_debug(namespace = "metadeconfoundR", paste("LRT_randOnly_for", aFeature, aCovariate, aP_mixed, sep = "\t"))

        if (collectMods) {
          collectedMods[[aFeature]][[aCovariate]][["randomFixedEffectsOnly"]] <- list()
          collectedMods[[aFeature]][[aCovariate]][["randomFixedEffectsOnly"]][["full"]] <- mixedmodel1
          collectedMods[[aFeature]][[aCovariate]][["randomFixedEffectsOnly"]][["small"]] <- mixedmodel2
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
           length (covariates[which(minQValues[i, ] < QCutoff)]) > 1)) {

        #status <- "STRICTLY DECONFOUNDED"
        status <- "OK_sd"
        # hardest to reach, means no covariate eliminates this signal

        # create vector to collect all confounder names for this feature <-> covariate pair
        # put new dataframe of qs in here
        otherCovariates <- lCovariates
        if (!is.null(minQValues[[1]])) {
          otherCovariates <- unique(c(covariates[which(minQValues[i, ] < QCutoff)], lCovariates))
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

          modAlg <- "stats::lm (rank (FeatureValue) ~ "
          lastPart <- paste0(", data = subsubMerge)", collapse = "")


          if (logistic == T || rawCounts == T) {
            modAlg <- "stats::glm (FeatureValue ~ "
            lastPart <- ", data = subsubMerge, family = \"binomial\")"

            if (rawCounts == T) {
              modAlg <- "stats::glm (cbind(FeatureValue, totReadCount) ~ "

            }

            logger::log_debug(namespace = "metadeconfoundR", "testing separation")
            glmmodeltext <- paste0 (modAlg,
                                    paste0(c(aCovariate, anotherCovariate),
                                           collapse = " + "),
                                    lastPart)
            glmmodel <- eval(parse(text = as.character(glmmodeltext)))
            if (update(glmmodel, method="detect_separation")$outcome) {
              logger::log_warn(namespace = "metadeconfoundR",
                paste(
                  "Separation for:", aFeature, ",",
                  aCovariate, "and", anotherCovariate))
              status <- NA
              if (collectMods) {
                collectedMods[[aFeature]][[aCovariate]][[anotherCovariate]][["full"]] <- glmmodel
              }
              next
            }

          }
          # else if (rawCounts == TRUE) { # alternative behavior for not rarefied abundances
          #   modAlg <- "stats::glm (cbind(FeatureValue, totReadCount) ~ "
          #   lastPart <- ", data = subsubMerge, family = \"binomial\")"
          # }

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
            } # end if (logistic == TRUE)
            else if (rawCounts == TRUE) { # alternative behavior for not rarefied abundances
              modAlg <- "lme4::glmer (cbind(FeatureValue, totReadCount) ~ "
              lastPart <- paste0(randomVarLine,
                                 ", data = subsubMerge, family = \"binomial\", nAGQ = ",
                                 nAGQ,
                                 ")",
                                 collapse = "")
            }
          } # end if (!is.na(randomVar[[1]]))

          logger::log_debug(namespace = "metadeconfoundR",
            paste("LRTsFwdRvsFor",
                  aFeature,
                  aCovariate,
                  anotherCovariate,
                  sep = "\t"
                  ))

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

            status <- "AD"
            logger::log_debug(namespace = "metadeconfoundR", "This should never be executed!")
            next
          }

          aP_forward <- safe_lrtest(lmBoth, lmAnother) # effect of adding aCovariate
          aP_reverse <- safe_lrtest(lmBoth, lmA) #effect of adding anotherCovariate

          if (is(lmBoth, "lm") && anyNA(lmBoth$coefficients)) {
            logger::log_warn(namespace = "metadeconfoundR", paste0("In full model containing: ",
                             paste0(names(lmBoth$coefficients),
                                    collapse = ", "),
                             ", NA coefficient(s) are present: ",
                             paste0(names(lmBoth$coefficients)[is.na(lmBoth$coefficients)],
                                    collapse = ", "),
                             ". Setting forward and reverse LRTs to non-significant."
                             ))
            aP_forward <- 1
            aP_reverse <- 1
          }

          # additonal control of confidence intervals for the covariates within the linear models
          conf_aCovariate <- TRUE
          conf_anotherCovariate <- TRUE
          if (doConfs > 0 && #          doConfs = 1 --> just logging
              !is(lmBoth, "logical") && #  the full model is working
              !is.na (aP_forward) && #     the lrts worked
              !is.na (aP_reverse) &&
              (aP_reverse < PHS_cutoff || # at least one of the lrts is significant
               aP_forward < PHS_cutoff) &&
              is.numeric(subsubMerge[, aCovariate]) && # the covariates are numeric
              is.numeric(subsubMerge[, anotherCovariate])
              ) {
            # if (is.numeric(subsubMerge[, aCovariate]) && is.numeric(subsubMerge[, anotherCovariate])) { # categorical variables are excluded for easier processing
            confints <- tryCatch({
              suppressMessages(confint(lmBoth, parm = c(aCovariate, anotherCovariate)))
            }, error = function(e) {
              NA
            })
            conf_aCovariate <- tryCatch({
              sign(confints[aCovariate, 1]) != sign(confints[aCovariate, 2])
              # signs are different, if confint is spanning 0
            }, error = function(e) {
              NA
            })

            conf_anotherCovariate  <- tryCatch({
              sign(confints[anotherCovariate, 1]) != sign(confints[anotherCovariate, 2])
            }, error = function(e) {
              NA
            })

            if (is.na(conf_aCovariate) || # if something went wrong with calculating confInts
                ((aP_forward < PHS_cutoff) && # if forward test is significant, but aCovariate confint spans 0
                conf_aCovariate)) {

              ending <- 'is spanning 0.'

              if (is.na(conf_aCovariate)) {
                ending <- 'is NA. Test for high collinearity of these metavariables
                (including fixed/random effects if applicable)!'
              }

              logger::log_warn(namespace = "metadeconfoundR",
                paste(
                  'lrt: ',
                  features[i],
                  aCovariate,
                  anotherCovariate,
                  '-- forward linear model is < PHS_cutoff, but confidence intervall for',
                  aCovariate,
                  ending
                ))
              if (doConfs > 1) {
                # if doConfs ==2 make lrt non-significant
                # aP_forward <- 1 # --> no unique effect of aCovariate
                # setting aP_forward <- 1 leads to labeling as confounded or AD,
                # so "OK_d" can never be kept
                status <- "OK_d"
              }

            }

            #if (! is.na (aP_forward) && (aP_reverse < PHS_cutoff) && conf_anotherCovariate) { # if reverse test is significant, but anotherCovariate confint spans 0
            if (is.na(conf_anotherCovariate) || ((aP_reverse < PHS_cutoff) &&
                conf_anotherCovariate)) {
              # if reverse test is significant, but anotherCovariate confint spans 0

              ending <- 'is spanning 0.'

              if (is.na(conf_anotherCovariate)) {
                ending <- 'is NA. Test for high collinearity of these metavariables
                (including fixed/random effects if applicable)!'
              }

              logger::log_warn(namespace = "metadeconfoundR",
                paste(
                  'lrt: ',
                  features[i],
                  aCovariate,
                  anotherCovariate,
                  '-- reverse linear model is < PHS_cutoff, but confidence intervall for',
                  anotherCovariate,
                  ending
                ))
              if (doConfs > 1) {
                aP_reverse <- 1 # --> no unique effect of anotherCovariate
                # aCovariate will not be labeled as being confounded by anotherCovariate anymore
              }
            }
          } # if doConfs

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

          logger::log_debug(namespace = "metadeconfoundR",
            paste(
              "LRTsFwdRvsFor",
              aFeature,
              aCovariate,
              anotherCovariate,
              aP_forward,
              aP_reverse,
              sep = "\t"
            ))
        } #end "for (anotherCovariate in lCovariates) {"


      } else {
        #status <- "NO COVARIATES"
        status <- "OK_nc"# trivially unconfounded because no
                        #   other features play a role
      }

      # if confounders where detected, they will all be printed as a single string as status
      if (!is.null(confounders)) {
        status <- paste0("C: ", paste(confounders, collapse = ", "))
      }
      # }# status assignment for a sinlge feature <-> covariate combination

      statusLine[j] <- status
    }# inner nested for j loop

    if ((i %% progressSteps) == 0) {#TB20240229
      progress <- paste0(round(x = ((i/length(features))*100),
                               digits = 2), "%")
      logger::log_info(namespace = "metadeconfoundR", paste("Deconfounding -- processed",
                            progress,
                            "of features."))
    }
    logger::log_debug(namespace = "metadeconfoundR", "one more line done!")
    statusLine
  }# foreach loop

  parallel::stopCluster(cl) # close parallel processing environment
  logger::log_debug(namespace = "metadeconfoundR", "Everything done in checkReducibility!")

  logger::log_info(namespace = "metadeconfoundR", paste("Deconfounding -- processed 100% of features."))

  if (is.null(ncol(r))) {
    r <- t(r)
  }
  rownames(r) <- features
  colnames(r) <- covariates

  # 2025 08 27 remove confounded confounders
  if (noConfConfs) {
    logger::log_info(namespace = "metadeconfoundR", paste("Removing 'confounded confounders' from status labels. Set noConfConfs = FALSE to keep them."))
    for (i in seq_along(rownames(r))) {
      if (any(grepl("^C:", r[i, ]))) {
        confounders <- sub("C: ", "", r[i, ])
        confounders <- stringr::str_split(confounders, ", ")
        names(confounders) <- colnames(r)
        for (j in seq_along(colnames(r))) {
          if (any(confounders[[j]] %in% colnames(r))) {
            confounded_confounder <- c()
            for (k in seq_along(confounders[[j]])) {
              ks_confounders <- confounders[[ confounders[[j]][k] ]]
              if (any(ks_confounders %in% confounders[[j]])) {
                confounded_confounder <- c(confounded_confounder, k)
              }
            }
            if (length(confounded_confounder) > 0) {
              confounders[[j]] <- confounders[[j]][-confounded_confounder]
              r[i, j] <- paste0("C: ", paste(confounders[[j]], collapse = ", "))
            }
          }
        }
      }
    }
  }

  if (collectMods) {
    r <- list(r, collectedMods)
  }
  r
}# function body
