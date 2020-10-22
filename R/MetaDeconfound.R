#' Metadeconfound
#'
#' Metadeconfound checks all feature <-> covariate combinations for
#' counfounding effects of covariates on feature <-> effect corellation
#'
#' @param featureMat a data frame with row(sample ID)
#' and column(feature such as metabolite or microbial OTU )
#' names, listing features for all samples
#' @param metaMat a data frame with row(sample ID) and
#' column(meta data such as age,BMI and all possible confounders)
#' names listing metadata for all samples. first column should be case status
#' with case=1 and control=0. All binary variables need to be in 0/1 syntax!
#' @param nnodes number of nodes/cores to be used for parallel processing
#' @param adjustMethod multiple testing p-value correction using one of the
#' methods of stats::p.adjust.methods
#' @param robustCutoff minimamal number of sample size for each covariate
#' in order to have sufficient power for association testing
#' @param QCutoff significance cutoff for q-value, DEFAULT = 0.1
#' @param DCutoff effect size cutoff
#' (either cliff's delta or spearman correlation test estimate), DEFAULT = 0
#' @param PHS_cutoff PostHoc Significance cutoff
#' @param NA_imputation Missing data treatment. "remove": remove NA containing
#' data from analysis (default).
#' "others": Impute NAs  using methods from packages MICE or AMELIA
#' (not yet implemented)
#' @param logfile name of optional logging file
#' @param intermediateOutput name base for optional intermediate files
#' (Ps, Qs and Ds from final ouput)
#' @param startStop vector of optional strings controlling which
#' parts of the pipeline should be executed.
#' ("naiveStop": only naive associations will be computed, no confounder analysis is done)
#' @param QValues optional data.frame containing pre-computed multiple-testing corrected p-values for naive associations
#' @param DValues optional data.frame containing pre-computed effect sizes for naive associations
#' @param minQValues pessimistic qvalues
#' @param deconfT vector of metavariable names *always* to be included as potenital confounder
#' @param deconfF vector of metavariable names *never* to be included as potenital confounder
#' @param doConfs optional parameter for additional computation of confidence
#' interval of linear models in the deconfounding step
#' (0 = no , 1 = logging, 2 = strict)
#' @param doRanks optional vector of metavariable names, that should be rank
#' transformed when building linear models in the doconfounding step
#' @param randomVar optional list with 2 elements. First: character string that should be appended to all
#' linear models used in the confounder analysis to account for random effect
#' variables. This string is part of the input to lme4::lmer() so syntax must
#' concur with that of lmer(). Second: vector of all variable names used in the string.
#' (e.g list("+ (1 | D) + (1 | E:D)" , c("D", "E")))
#' @param robustCutoffRho optional robustness cutoff for continuous variables
#' @param typeCategorical optional character vector of metavariable names to
#' always be treated as categorical
#' @param typeContinuous optional character vector of metavariable names to
#' always be treated as continuous
#' @param logistic optional logical parameter; DEFAULT = FALSE;
#' Set TRUE to treat supplied features as binary instead of continuous
#' @param ... for additional arguments used internally (development/debugging)
#' @return list with elements Ds = effectsize,
#' Ps = uncorrected p-value for naive association,
#' Qs = multiple testing corrected p-value/fdr,
#' and status = confounding/mediation status for all
#' feature <=> covariate combinations with following categories:
#' (NS = not significant, SD = strictly deconfounded, LD = laxly deconfounded,
#' NC = no covariates, "covariate name" = confounded by this covariate)
#' @details for more details and explanations please see the vignette.
#' @examples
#'data(reduced_feature)
#'data(metaMatMetformin)
#'\donttest{
#'example_output <- MetaDeconfound(featureMat = reduced_feature,
#'                                   metaMat = metaMatMetformin)}
#'
#' @import futile.logger
#' @importFrom utils write.table
#' @export
#'


MetaDeconfound <- function(featureMat,
                           metaMat,
                           nnodes = 1,
                           adjustMethod = "fdr",
                           robustCutoff = 5,
                           QCutoff = 0.1,
                           DCutoff = 0,
                           PHS_cutoff = 0.05,
                           NA_imputation = "remove",
                           logfile = NULL,
                           intermediateOutput = NULL,
                           startStop = NA,
                           QValues = NA,
                           DValues = NA,
                           minQValues= NULL,
                           deconfT = NULL,
                           deconfF = NULL,
                           doConfs = 0,
                           doRanks = NA,
                           randomVar = NA,
			   robustCutoffRho = NULL, # new SKF20200221
			   typeCategorical = NULL, # new SKF20200221
			   typeContinuous = NULL, # new SKF20200221
			   logistic = FALSE, # new SKF20201017
                           ...) {
  # create file for progress log


  flog.logger("my.logger", INFO, appender=appender.file(logfile))
  flog.appender(appender.tee(logfile), name = 'my.logger')
  #flog.threshold(INFO, name = 'my.logger')
  if (is.null(logfile)) {
    flog.threshold(FATAL, name = 'my.logger')
  }

  ###
  ###
  ### logging start and initial sanity checks on input paramters
  ###
  ###

  if (is.null(logfile) && (doConfs > 0)) {
    stop('Error - "doConfs" parameeter is set to > 0 but "logfile" is not specified. Can not log warnings for confidence intervalls spanning 0.')
  }

  flog.info(msg = '###',
            name = "my.logger")
  flog.info(msg = '###',
            name = "my.logger")
  flog.info(msg = 'Deconfounding run started',
            name = "my.logger")
  if ("naiveStop" %in% startStop) {
    flog.warn(msg = 'Detected "naiveStop" in "startStop" paramter. Will only return naive associations.',
              name = "my.logger")
  }

  if (missing(metaMat)) {
    flog.error(msg = 'Necessary argument "metaMat" missing.',
               name = "my.logger")
    stop('Error - Necessary argument "metaMat" missing.')
  }
  if (missing(featureMat)) {
    flog.error(msg = 'Necessary argument "featureMat" missing.',
               name = "my.logger")
    stop('Error - Necessary argument "featureMat" missing.')
  }
  if (nrow(metaMat) != nrow(featureMat)) {
    flog.error(msg = "featureMat and metaMat don't have same number of rows.",
               name = "my.logger")
    stop("featureMat and metaMat don't have same number of rows.")
  }
  if (!is.null(deconfT) | !is.null(deconfF)) {
    if ((sum(deconfT %in% colnames(metaMat)) < length(deconfT)) |
        (sum(deconfF %in% colnames(metaMat)) < length(deconfF))) {
      flog.error(msg = "Elements of deconfT/deconfF are not present in colnames of metaMat.",
                 name = "my.logger")
      flog.info(msg = "Check identical spelling of variable
                      names in deconfT/deconfF and metaMat",
                name = "my.logger")
      stop("Elements of deconfT/deconfF are not present in colnames of metaMat.")
    } else if (sum(deconfT %in% deconfF) > 0) {
      flog.error(msg = "Some elements of deconfT and deconfF seem to be identical.",
                 name = "my.logger")
      stop("Some elements of deconfT and deconfF seem to be identical.")
    }
  }

  if (is.null(QValues)) {
    flog.error(msg = "QValues argument is supplied but seems to be empty (NULL).",
               name = "my.logger")
    stop("QValues argument is supplied but seems to be empty (NULL).")
  }
  if (is.null(DValues)) {
    flog.error(msg = "DValues argument is supplied but seems to be empty (NULL).",
               name = "my.logger")
    stop("DValues argument is supplied but seems to be empty (NULL).")
  }



  .MetaDeconfound(featureMat = featureMat,
                   metaMat = metaMat,
                   nnodes = nnodes,
                   adjustMethod = adjustMethod,
                   robustCutoff = robustCutoff,
                   QCutoff = QCutoff,
                   DCutoff = DCutoff,
                   PHS_cutoff = PHS_cutoff,
                   logfile = logfile,
                   intermediateOutput = intermediateOutput,
                   startStop = startStop,
                  QValues = QValues,
                  DValues = DValues,
                  minQValues=minQValues,
                  deconfT = deconfT,
                  deconfF = deconfF,
                  doConfs = doConfs,
                  doRanks = doRanks,
                  randomVar = randomVar,
		  robustCutoffRho = robustCutoffRho, # new SKF20200221
		  typeCategorical = typeCategorical, # new SKF20200221
		  typeContinuous = typeContinuous, # new SKF20200221
		  logistic = logistic, # new SKF20201017
                   ...
                   #maintenance = maintenance,
                   #verbosity = verbosity
                   )
}

.MetaDeconfound <- function(featureMat,
                            metaMat,
                            nnodes = 1,
                            robustCutoff = 5,
                            adjustMethod = "fdr",
                            QCutoff = 0.1,
                            DCutoff = 0,
                            PHS_cutoff = 0.05,
                            NA_imputation = "remove",
                            logfile = NULL,
                            intermediateOutput = NULL,
                            startStop = NA,
                            QValues = NA,
                            DValues = NA,
                            minQValues=NULL,
                            deconfT = NULL,
                            deconfF = NULL,
                            doConfs = 0,
                            doRanks = NA,
                            randomVar = NA,
                            maintenance = FALSE,
			    robustCutoffRho = NULL, # new SKF20200221
			    typeCategorical = NULL, # new SKF20200221
 			    typeContinuous = NULL, # new SKF20200221
			    logistic = FALSE, # new SKF20201017
                            verbosity = "silent") {


  samples <- row.names (featureMat)
  features <- colnames (featureMat)
  noFeatures <- length (features)

  if (is.null (robustCutoffRho)) { robustCutoffRho <- robustCutoff } # new SKF20200221

  RVnames <- NA
  covariates <- colnames (metaMat) # each covariate + the status category
  if (!is.na(randomVar[[1]])) { # list input paramter is split for further use within pipeline
    RVnames <- randomVar[[2]]
    randomVar<- randomVar[[1]]
    #covariates <- covariates[!(covariates %in% RVnames)] # random effect variables are not treated as covariates

    flog.info(msg = paste0("The following parameters will be added to all linear models: '", randomVar, "'"),
              name = "my.logger")
    flog.warn(msg = paste0("the following random effect covariates will be excluded as potential donfounders: ", paste0(RVnames, collapse = ", ")),
              name = "my.logger")
  }
  noCovariates <- length (covariates)

  if (nnodes > 1) {
    nnodes <- nnodes - 1
  }

  if (NA_imputation != "remove") {
    #impute using mice or somn like dat
    #featureMat <- mice(featuremat)
    #metaMat <- mice(metaMat)
  }


  flog.info(msg = paste0("Checking robustness of data for covariates"),
            name = "my.logger")

  isRobust <- CheckSufficientPower(metaMat = metaMat,
                                   covariates = covariates,
                                   noCovariates = noCovariates,
                                   nnodes = nnodes,
                                   robustCutoff = robustCutoff,
                                   robustCutoffRho = robustCutoffRho, # new SKF20200221
                                   typeCategorical = typeCategorical, # new SKF20200221
                                   typeContinuous = typeContinuous, # new SKF20200221
                                   NA_imputation = NA_imputation,
                                   maintenance = maintenance,
                                   verbosity = verbosity)

  if (verbosity == "debug") {
    print("CheckSufficientPower -- head(isRobust):")
    print(dim(isRobust[[2]]))
  }

  if (anyNA(isRobust[[1]]) || anyNA(isRobust[[2]])) {
    flog.warn(msg = paste0(sum(is.na(isRobust[[1]])), " covariates where marked as too sparse and won't be considered in further analysis due to lack of sufficient data."),
              name = "my.logger")
    flog.warn(msg = paste0(sum(is.na(isRobust[[2]])), " covariate combinations where marked as too sparse and won't be considered as potential confunders of each other in the deconfounding step."),
              name = "my.logger")
  }



  if (is.na(QValues[[1]]) | is.na(DValues[[1]])) {
    flog.info(msg = paste0("Computation of naive associations started."),
              name = "my.logger")
    naiveAssociation <- NaiveAssociation(featureMat = featureMat,
                                         samples = samples,
                                         features = features,
                                         noFeatures = noFeatures,
                                         metaMat = metaMat,
                                         covariates = covariates,
                                         noCovariates = noCovariates,
                                         isRobust = isRobust,
					 typeCategorical = typeCategorical, # new SKF20200221
					 typeContinuous = typeContinuous, # new SKF20200221
					 logistic = logistic, # new SKF20201017
                                         adjustMethod = adjustMethod,
                                         nnodes = nnodes,
                                         maintenance = maintenance,
                                         verbosity = verbosity)
    if (verbosity == "debug") {
      print(naiveAssociation$Ps[seq_len(3), seq_len(noCovariates)])
      print(naiveAssociation$Qs[seq_len(3), seq_len(noCovariates)])
      print(naiveAssociation$Ds[seq_len(3), seq_len(noCovariates)])
      print("now computing confounding status")
    }

    if (!is.null(intermediateOutput)) {
      write.table(x = naiveAssociation$Ps,
                  file = paste0(intermediateOutput, "_Ps.csv"),
                  sep = "\t",
                  col.names=NA)
      write.table(x = naiveAssociation$Qs,
                  file = paste0(intermediateOutput, "_Qs.csv"),
                  sep = "\t",
                  col.names=NA)
      write.table(x = naiveAssociation$Ds,
                  file = paste0(intermediateOutput, "_Ds.csv"),
                  sep = "\t",
                  col.names=NA)

      flog.info(msg = paste('Naive Associations written to', intermediateOutput, "and returned as list within R."),
                name = "my.logger")

      flog.info(msg = 'Done!',
                name = "my.logger")
    }

    if ("naiveStop" %in% startStop) {
      flog.warn(msg = paste('Process stopped before computing confounding status because "startStop" parameter contained "naiveStop". '),
                name = "my.logger")

      return(list(Ps = naiveAssociation$Ps,
                  Qs = naiveAssociation$Qs,
                  Ds = naiveAssociation$Ds))
    }
  } else { # if precomputed Qs and Ds are supplied as arguments
    naiveAssociation <- list(Qs = QValues, Ds = DValues)
    flog.warn(msg = paste('Confonding status is computed based on Q-values and effect sizes supplied via "QValue" and "DValue" parameters. '),
              name = "my.logger")
  }



  reducibilityStatus <- CheckReducibility(featureMat = featureMat,
                                          metaMat = metaMat,
                                          noFeatures = noFeatures,
                                          noCovariates = noCovariates,
                                          features = features,
                                          covariates = covariates,
                                          Qs = naiveAssociation$Qs,
                                          Ds = naiveAssociation$Ds,
                                          minQValues= minQValues,
                                          nnodes = nnodes,
                                          QCutoff = QCutoff,
                                          DCutoff = DCutoff,
                                          PHS_cutoff = PHS_cutoff,
                                          deconfT = deconfT,
                                          deconfF = deconfF,
                                          doConfs = doConfs,
                                          doRanks = doRanks,
                                          randomVar = randomVar,
                                          RVnames = RVnames,
                                          isRobust = isRobust,
                                          logistic = logistic, # new SKF20201017
                                          verbosity = verbosity)

  if (verbosity == "debug") {
    print(utils::head(reducibilityStatus))
    print("MultiDeconfound  --  All done!")
  }

  flog.info(msg = "MetadecondoundR run completed successfully!",
            name = "my.logger")

  return(list(Ps = naiveAssociation$Ps,
              Qs = naiveAssociation$Qs,
              Ds = naiveAssociation$Ds,
              status=reducibilityStatus))
}
