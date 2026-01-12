#' GetPartialEfSizes
#'
#' GetPartialEfSizes takes \link[metadeconfoundR]{MetaDeconfound} output and
#' genarates partial effect sizes for all significant associations
#'
#' @param featureMat a data frame with row(sample ID)
#' and column(feature such as metabolite or microbial OTU )
#' names, listing features for all samples
#' @param metaMat a data frame with row(sample ID) and
#' column(meta data such as age,BMI and all possible confounders)
#' names listing metadata for all samples. first column should be case status
#' with case=1 and control=0. All binary variables need to be in 0/1 syntax!
#' @param metaDeconfOutput long format output of the MetaDeconfound output created for the
#' supplied featureMat and metaMat
#' @param doRanks optional vector of metavariable names, that should be rank
#' transformed when building linear models in the doconfounding step
#' @param randomVar optional vector of metavariable names to be treated as
#' random effect variables. These variables will not be tested for naive
#' associations and will not be included as potential confounders,
#' but will be added as random effects "+ (1|variable)" into any models being built.
#' Any associations reducible to the supplied random effect(s) will be labeled
#'  as "NS". Note: Ps, Qs, Ds are computed independently and thereby not changed
#'  through inclusion of random effects.
#' @param fixedVar optional vector of metavariable names to be treated as
#' fixed effect variables. These variabels will not be tested for naive
#' associations and will not be included as potential confounders,
#' but will be added as fixed effects "+ variable" into any models being built.
#' Any associations reducible to the supplied fixed effect(s) will be labeled
#' as "NS". Note: Ps, Qs, Ds are computed independently and thereby not changed
#' through inclusion of fixed effects.

#' @return long format data.frame similar to Metadeconfound() output
#' @details for more details and explanations please see the package vignette.
#' @examples
#'data(reduced_feature)
#'data(metaMatMetformin)
#'\donttest{
#'
#'example_output <- MetaDeconfound(featureMat = reduced_feature,
#'                                   metaMat = metaMatMetformin,
#'                                   logLevel = "ERROR",
#'                                   returnLong = TRUE)
#' #
#' ex_out_partial <- GetPartialEfSizes(featureMat = reduced_feature,
#'                                       metaMat = metaMatMetformin,
#'                                       metaDeconfOutput = example_output)
#'}
#'
#' @importFrom rlang .data
#' @import lme4
#' @export

GetPartialEfSizes <- function(featureMat,
                              metaMat,
                              metaDeconfOutput,
                              doRanks = NA,
                              randomVar = NA,
                              fixedVar = NA) {
  output <- metaDeconfOutput
  # add new columns to the metadeconfoundR output
  output$partial <- NA
  output$partialRel <- NA
  output$partialNorm <- NA
  output$metaVariable <- as.character(output$metaVariable)
  for (multiFeat in unique(output$feature)) {
    subTemp <- output[((output$feature == multiFeat) & (!is.na(output$status) & output$status != "NS")), ]
    if (nrow(subTemp) < 1) {
      # if no metavariable at all is associated, no partial effect size can be calculated
      next
    }
    subMerge <- metaMat
    subMerge$FeatureValue <- featureMat [, multiFeat]
    subTemp_metaVars <- as.character(subTemp$metaVariable)
    subSubMerge <- na.exclude(subMerge[, c(subTemp_metaVars, "FeatureValue")])
    fullModForm <- as.formula(paste(
      "rank (FeatureValue) ~",
      paste(subTemp_metaVars, collapse = " + ")
    ))
    if (!is.na(fixedVar[[1]])) {
      subSubMerge <- na.exclude(subMerge[, c(subTemp_metaVars, fixedVar, "FeatureValue")])
      fullModForm <- as.formula(paste(
        "rank (FeatureValue) ~",
        paste(c(subTemp_metaVars, fixedVar), collapse = " + ")
      ))
    }
    fullRsq <- summary(lm(fullModForm, data = subSubMerge))$r.squared
    if (!is.na(randomVar[[1]])) {
      subSubMerge <- na.exclude(subMerge[, na.exclude(c(subTemp_metaVars, fixedVar, randomVar, "FeatureValue"))])
      fullModForm <- as.formula(paste(
        "rank (FeatureValue) ~",
        paste(subTemp_metaVars, collapse = " + "),
        paste0("+ (1|", randomVar, ")")
      ))
      if (!is.na(fixedVar[[1]])) {
        fullModForm <- as.formula(paste(
          "rank (FeatureValue) ~",
          paste(c(subTemp_metaVars, fixedVar), collapse = " + "),
          paste0("+ (1|", randomVar, ")")
        ))
      }
      fullRsq <-  ConditionalR2(lme4::lmer(fullModForm, REML = F, data = subSubMerge))
    }

    if (nrow(subTemp) == 1) {# if only one metavariable is associated
      # we test either against a null model or a model only containing random+fixed effects
      reducedModForm <- as.formula("rank (FeatureValue) ~ 1")
      reducedRsq <- summary(lm(reducedModForm, data = subSubMerge))$r.squared
      if (!is.na(fixedVar[[1]])) {
        reducedModForm <- as.formula(paste(
          "rank (FeatureValue) ~",
          paste(fixedVar, collapse = " + ")
        ))
        reducedRsq <- summary(lm(reducedModForm, data = subSubMerge))$r.squared
        if (!is.na(randomVar[[1]])) {
          reducedModForm <- as.formula(paste(
            "rank (FeatureValue) ~",
            paste(fixedVar, collapse = " + "),
            paste0("+ (1|", randomVar, ")")
          ))
          reducedRsq <-  ConditionalR2(lme4::lmer(reducedModForm, REML = F, data = subSubMerge))
        }
      } else if (!is.na(randomVar[[1]])) {
        reducedModForm <- as.formula(paste(
          "rank (FeatureValue) ~", paste(paste0("(1|", randomVar, ")"), collapse = " + ")
        ))
        reducedRsq <-  ConditionalR2(lme4::lmer(reducedModForm, REML = F, data = subSubMerge))
      }
      partialRsq <- fullRsq - reducedRsq
      relPartialRsq <- partialRsq/fullRsq
      normPartialRsq <- partialRsq/(partialRsq + 1 - fullRsq)
      # 1- fullRSQ == residual R2
      # (partialRsq + 1 - fullRsq) is equal to (1 - reducedRsq)
      if (!is.na(subTemp$Ds) & (subTemp$Ds < 0)) {
        # negative effect size
        partialRsq <- partialRsq * (-1)
        relPartialRsq <- relPartialRsq * (-1)
        normPartialRsq <- normPartialRsq * (-1)
      }
      output[((output$feature == multiFeat) &
                (output$metaVariable == subTemp$metaVariable)), "partial"] <- partialRsq
      output[((output$feature == multiFeat) &
                (output$metaVariable == subTemp$metaVariable)), "partialRel"] <- relPartialRsq
      output[((output$feature == multiFeat) &
                (output$metaVariable == subTemp$metaVariable)), "partialNorm"] <- normPartialRsq
      output <- rbind(output, c(multiFeat, "maxRsq", 0.5, 0.5, Inf, "NS", fullRsq, 1, fullRsq))
      output$feature <- as.factor(output$feature)
      output$metaVariable <- as.factor(output$metaVariable)
      output$Ps <- as.numeric(output$Ps)
      output$Qs <- as.numeric(output$Qs)
      output$Ds <- as.numeric(output$Ds)
      output$partial <- as.numeric(output$partial)
      output$partialRel <- as.numeric(output$partialRel)
      output$partialNorm <- as.numeric(output$partialNorm)
      next
    }
    # iteratively remove one of the individually significant metavariables from the model
    for (metaVar in subTemp_metaVars) {
      if (metaVar %in% c(randomVar, fixedVar)) {
        next
      }
      reducedModForm <- as.formula(paste(
        "rank (FeatureValue) ~",
        paste(na.exclude(c(subTemp_metaVars[-which(subTemp_metaVars %in% c(metaVar, randomVar))], fixedVar)), collapse = " + ")
        # keeping randomVar in previous line might not be necessary, as randomVars should not be in subTemp anyways
      ))
      reducedRsq <- summary(lm(reducedModForm, data = subSubMerge))$r.squared
      if (!is.na(randomVar[[1]])) {
        reducedModForm <- as.formula(paste(
          "rank (FeatureValue) ~",
          paste(na.exclude(c(subTemp_metaVars[-which(subTemp_metaVars %in% c(metaVar, randomVar))], fixedVar)), collapse = " + "),
          paste0("+ (1|", randomVar, ")")
        ))
        reducedRsq <-  ConditionalR2(lme4::lmer(reducedModForm, REML = F, data = subSubMerge))
      }

      partialRsq <- fullRsq - reducedRsq
      relPartialRsq <- partialRsq/fullRsq
      normPartialRsq <- partialRsq/(partialRsq + 1 - fullRsq) # 1- fullRSQ == residual R2
      if (!is.na(subTemp[subTemp$metaVariable == metaVar, "Ds"]) &
          subTemp[subTemp$metaVariable == metaVar, "Ds"] < 0) {
        # negative effect size
        partialRsq <- partialRsq * (-1)
        relPartialRsq <- relPartialRsq * (-1)
        normPartialRsq <- normPartialRsq * (-1)
      }
      output[((output$feature == multiFeat) &
                (output$metaVariable == metaVar)), "partial"] <- partialRsq
      output[((output$feature == multiFeat) &
                (output$metaVariable == metaVar)), "partialRel"] <- relPartialRsq
      output[((output$feature == multiFeat) &
                (output$metaVariable == metaVar)), "partialNorm"] <- normPartialRsq
    }
    output <- rbind(output,      # add one row per feature, giving R2 of the full model
                    c(multiFeat,
                      "maxRsq",
                      0.5,
                      0.5,
                      Inf,
                      "NS",
                      fullRsq,
                      1,
                      fullRsq
                      )
                    )
    output$feature <- as.factor(output$feature)
    output$metaVariable <- as.factor(output$metaVariable)
    output$Ps <- as.numeric(output$Ps)
    output$Qs <- as.numeric(output$Qs)
    output$Ds <- as.numeric(output$Ds)
    output$partial <- as.numeric(output$partial)
    output$partialRel <- as.numeric(output$partialRel)
    output$partialNorm <- as.numeric(output$partialNorm)
  } # for multifeat
  return(output)
}
