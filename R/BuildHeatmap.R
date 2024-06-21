#' BuildHeatmap
#'
#' BuildHeatmap summarizes \link[metadeconfoundR]{MetaDeconfound} output in a heatmap or cuneiform plot
#'
#' @param metaDeconfOutput output of a metadeconfound run
#' @param q_cutoff optional FDR-value cutoff used to remove
#' low-significance entries from data
#' @param d_cutoff optional effect size cutoff used to remove
#' low effect size entries from data
#' @param d_range range of effect sizes shown; "full": (default) range from
#' -1 to +1;
#' "fit": range reduced according to maximum and minimum effect size
#' present in resulting plot
#' @param d_col set color range for effect size as c(minimum, middle, maximum),
#' default c("red", "white", "blue")
#' @param cuneiform optional logical parameter,
#' plot cuneiform instead of heatmap when cuneiform = TRUE
#' @param coloring optional, can be 0,1,2;
#' 0: color all tiles according to effectsize ;
#' 1: don't color not significant tiles
#' 2: like 1 but also don't color confounded signal tiles
#' @param showConfounded optional logical parameter;
#' set to FALSE to remove significance markers from confounded signals
#' @param intermedData only return intermediate data for plotting, default = FALSE
#' @param featureNames optional two-column-dataframe containing corresponding
#' "human-readable" names to the "machine-readable" feature names used as
#' row.names in metaDeconfOutput. These human readable
#' names will be displayed in the final plot. First column: machine-readable,
#' second column: human-readable.
#' @param metaVariableNames optional two-column-dataframe containing
#' corresponding  "human-readable" names to the "machine-readable" metadata
#' names used as column names in metaDeconfOutput. These human readable
#' names will be displayed in the final plot. First column: machine-readable,
#' second column: human-readable.
#' @param keepMeta character vector of metavariable names
#' (corresponding to names in metaDeconfOutput), that should be shown in
#' resulting plot, even when they have no associations
#' passing d_cutoff and q_cutoff
#' @param keepFeature character vector of metavariable names
#' (corresponding to names in metaDeconfOutput), that should be shown in
#' resulting plot, even when they have no associations
#' passing d_cutoff and q_cutoff
#' @param trusted character vector of confounding status labels to be treated
#' as trustworthy, not-confounded signal. default = c("OK_sd", "OK_nc", "OK_d", "AD")
#' @param tileBordCol tile border color of  heatmap tiles, default: "black"
#' @param reOrder reorder features and/or metadata? possible options: c("both", "feat", "meta", "none"), default: "both"
#' @return ggplot2 object
#' @details for more details and explanations please see the package vignette.
#' @examples
#'data(reduced_feature)
#'data(metaMatMetformin)
#'\donttest{
#'example_output <- MetaDeconfound(featureMat = reduced_feature,
#'                                   metaMat = metaMatMetformin,
#'                                   logLevel = "ERROR")
#'
#'plotObject <- BuildHeatmap(example_output)
#'
#'alternativePlot <- BuildHeatmap(example_output, coloring = 2, showConfounded = FALSE)
#'}
#'
#' @import ggplot2
#' @importFrom reshape2 melt
#' @importFrom methods is
#' @importFrom rlang .data
#' @export

BuildHeatmap <- function(metaDeconfOutput,
                         q_cutoff = 0.1,
                         d_cutoff = 0.01,
                         cuneiform = FALSE,
                         coloring = 0,
                         showConfounded = TRUE,
                         intermedData = FALSE,
                         featureNames = NULL,
                         metaVariableNames = NULL,
                         d_range = "fit",
                         d_col = c("blue", "white", "red"),
                         keepMeta = NULL,
                         keepFeature = NULL,
                         trusted = c("OK_sd", "OK_nc", "OK_d", "AD"),
                         tileBordCol = "black",
                         reOrder = "both"
                         ) {


  if (length(d_col) != 3) {
    stop("wrong number of colors in d_col!\nSupply colors for c(min, middle, max)!")
  }
  if (!(d_range %in% c("fit", "full"))) {
    stop('d_range must be either "fit" or "full"!')
  }

  if (length(trusted) == 0) {
    stop('"trusted" must contain at least one trusted status label')
  }

  fromIntermed <- FALSE

  allLables <- c ("OK_sd", "OK_nc", "OK_d", "AD", "NS")
  notTrusted <- allLables[!(allLables %in% trusted)]

  if (is(metaDeconfOutput, "list")) {
    # melt dataframes for fdr-valules, effectsizes, and confounding status
    effectSize <- reshape2::melt(data = metaDeconfOutput$Ds,
                       varnames = c("feature", "metaVariable"),
                       value.name = "Ds")

    # set NAs and Inf entries to 0
    effectSize$Ds[effectSize$Ds == Inf] <- 0
    effectSize$Ds[is.na(effectSize$Ds)] <- 0

    fdr <- reshape2::melt(data = metaDeconfOutput$Qs,
                varnames = c("feature", "metaVariable"),
                value.name = "Qs")
    # set NA Qvalues to 1
    fdr$Qs[is.na(fdr$Qs)] <- 1

    status <- reshape2::melt(data = metaDeconfOutput$status,
                   varnames = c("feature", "metaVariable"),
                   value.name = "status")
  } else if (ncol(metaDeconfOutput) == 9) {
    warning("treating input as 'intermedData = T' Buildheatmap output!!")
    fromIntermed <- TRUE
    effectSize <- metaDeconfOutput

  } else {
    effectSize <- metaDeconfOutput[, c("feature", "metaVariable", "Ds")]
    fdr <- metaDeconfOutput[, c("feature", "metaVariable", "Qs")]
    status <- metaDeconfOutput[, c("feature", "metaVariable", "status")]

    effectSize$Ds[effectSize$Ds == Inf] <- 0
    effectSize$Ds[is.na(effectSize$Ds)] <- 0
    fdr$Qs[is.na(fdr$Qs)] <- 1
  }

  if (!fromIntermed) {
  # identify all NS entries (both naive NS, and flagged "NS" by deconfounding step)
  insignificant <- unlist(lapply(strsplit(as.character(status$status), split = ", "),
                                 function(l) any(l %in% notTrusted)))
  # identify all non-confounded significant entries
  trueDeconf <- unlist(lapply(strsplit(as.character(status$status), split = ", "),
                              function(l) any(l %in% trusted)))


  # assign stars according to significance level
  effectSize$stars <- cut(fdr$Qs,
                          breaks = c(-Inf, 0.001, 0.01, 0.1, Inf),
                          label=c("***", "**", "*", ""))
  # remove NS entries
  effectSize$stars[insignificant] <- ""
  # pass on information on confounding status for plotting either * or \u00b0
  effectSize$status <- trueDeconf
  # change all confounded ( == !trueDeconf) fdr-value stars to "\u00b0" (degree symbol)
  # or set them to ""
  effectSize$stars <- as.character(effectSize$stars)
  for (m in seq_along(effectSize$stars)) {
    if (!effectSize$status[m] && length(effectSize$stars[m]) > 0) {
      if (showConfounded) {
        effectSize$stars[m] <- gsub("*",
                                    "\u00b0",
                                    effectSize$stars[m] ,
                                    fixed = TRUE)
      } else {
        effectSize$stars[m] <- ""
      }
    }
  }

  # reduce effectsizes for certain subsets so they are shown in white

  effectSize$insignificant <- insignificant
  effectSize$trueDeconf <- !trueDeconf

  if (coloring == 1) {
    effectSize$Ds[effectSize$insignificant] <- 0.000001
  }
  if (coloring == 2) {
    effectSize$Ds[effectSize$trueDeconf] <- 0.000001
  }


  #
  #
  #### fetching all rows/columns that have low Ds and Qs
  #### or don't contain any significant hits (i.e. all stars == "")
  #### and removing them
  #
  #
  remove_metavariables <- vector()
  for (i in unique(effectSize$metaVariable)) { # identify metavariables to be removed

    aMetaVariable <- fdr[fdr$metaVariable == i, ]
    aMetaVariableD <- effectSize[effectSize$metaVariable == i, ]
    if (sum(na.exclude(abs(aMetaVariable$Qs)) > q_cutoff) == length(na.exclude(aMetaVariable$Qs)) ||
        sum(na.exclude(abs(aMetaVariableD$Ds)) < d_cutoff) == length(na.exclude(aMetaVariableD$Ds)) ||
        all(aMetaVariableD$stars == "")
    ) {
      #remove <- c(remove, which(effectSize$metaVariable == i))
      remove_metavariables <- c(remove_metavariables, i)
    }
  }

  # remove metaVariabel names from the list of to-be-removed metaVariables
  if (!is.null(keepMeta)) {
    remove_metavariables <-
      remove_metavariables[!(remove_metavariables %in% keepMeta)]
  }
  # remove filtered out matavariables from the dataset
  effectSize <- effectSize[!(effectSize$metaVariable %in% remove_metavariables), ]

  remove <- vector()
  for (i in unique(effectSize$feature)) { # identify omics features to be removed
    aGenus <- fdr[fdr$feature == i, ]
    aGenusD <- effectSize[effectSize$feature == i, ]
    if (sum(na.exclude(abs(aGenus$Qs)) > q_cutoff) == length(na.exclude(aGenus$Qs)) |
        sum(na.exclude(abs(aGenusD$Ds)) < d_cutoff) == length(na.exclude(aGenusD$Ds)) |
        all(aGenusD$stars == "")
    ) {
      #remove <- c(remove, which(effectSize$feature == i))
      remove <- c(remove, i)
    }
  }

  remove <- as.vector(unique(remove))
  # check that there actually are features, that should be removed
  if (length(remove) > 0) {
    # remove feature names from the list of to-be-removed features
    if (!is.null(keepFeature)) {
      remove <- remove[!(remove %in% keepFeature)]
    }
    # check again for empty remove, because
    if (length(remove) > 0) {
      # remove features not passing q_cutoff and d_cutoff filtering
      effectSize <- effectSize[!(effectSize$feature %in% remove), ]
    }

  }
  effectSize <- droplevels(effectSize)


  if (length(unique(effectSize$metaVariable)) == 0) {
    stop("No associations pass current q_cutoff and/or d_cutoff filters!")
  }

  # cluster heatmap by reordering the factor levels for both dimensions of the heatmap
  # only if more than one feature or metavariable are present respectively
  eff_cast <- reshape2::dcast(effectSize,
                              effectSize[[1]] ~ metaVariable,
                              value.var = "Ds")
  rownames(eff_cast) <- eff_cast[[1]]
  eff_cast[[1]] <- NULL # move feature names to rownames

  if ((reOrder %in% c("both", "feat")) & (nrow(eff_cast) > 1)) {
    ord <-
      hclust(dist(eff_cast, method = "euclidean"), method = "ward.D")$order
    effectSize$feature <-
      factor(as.factor(effectSize$feature),
             levels = levels(as.factor(effectSize$feature)) [ord])
  }
  if ((reOrder %in% c("both", "meta")) & (ncol(eff_cast) > 1)) {
    eff_cast <- scale(t(eff_cast))
    ord2 <-
      hclust(dist(eff_cast, method = "euclidean"), method = "ward.D")$order
    effectSize$metaVariable <- droplevels(effectSize$metaVariable)
    effectSize$metaVariable <-
      factor(as.factor(effectSize$metaVariable),
             levels = levels(as.factor(effectSize$metaVariable)) [ord2])
  }


  # ord <- hclust(dist(eff_cast, method = "euclidean"), method = "ward.D")$order
  # eff_cast <- scale(t(eff_cast))
  # ord2 <- hclust(dist(eff_cast, method = "euclidean"), method = "ward.D")$order
  # effectSize$metaVariable <- droplevels(effectSize$metaVariable)
  #
  # effectSize$feature <- factor(as.factor(effectSize$feature),
  #                              levels = levels(as.factor(effectSize$feature)) [ord])
  # effectSize$metaVariable <- factor(as.factor(effectSize$metaVariable),
  #                                   levels = levels(as.factor(effectSize$metaVariable)) [ord2])


  effectSize$featureNames <- effectSize$feature
  effectSize$metaVariableNames <- effectSize$metaVariable

  if (!is.null(featureNames)) {
    if (!is(featureNames, "data.frame")) {
      warning('class(featureNames) was coerced to "data.frame"')
      featureNames <- as.data.frame(featureNames)
    }
    if (length(unique(featureNames[[2]])) != length(featureNames[[2]])) {
      featureNames[[2]] <- make.unique(featureNames[[2]])
      warning('non-unique human-readable feature names where made unique using base::make.unique')
    }
    map = stats::setNames(featureNames[[2]], featureNames[[1]])
    effectSize$featureNames <- map[as.vector(effectSize$feature)]


    effectSize$featureNames <- factor(as.factor(effectSize$featureNames),
                                      levels = map[levels(effectSize$feature)])
  }

  if (!is.null(metaVariableNames)) {
    if (!is(metaVariableNames, "data.frame")) {
      warning('class(metaVariableNames) was coerced to "data.frame"')
      metaVariableNames <- as.data.frame(metaVariableNames)
    }
    if (length(unique(metaVariableNames[[2]])) != length(metaVariableNames[[2]])) {
      metaVariableNames[[2]] <- make.unique(metaVariableNames[[2]])
      warning('non-unique human-readable metaVariable names where made unique using base::make.unique')
    }
    map = stats::setNames(metaVariableNames[[2]], metaVariableNames[[1]])
    effectSize$metaVariableNames <- map[as.vector(effectSize$metaVariable)]

      effectSize$metaVariableNames <- factor(as.factor(effectSize$metaVariableNames),
                                        levels = map[levels(effectSize$metaVariable)])
  }
}


  if (intermedData == TRUE) {
    if (!any(c("*", "**", "***") %in% unique(effectSize$stars))) {
      warning("No unconfounded associations remain with the current cutoff values. ")
    }
    return(effectSize)
  }

  if (!any(c("*", "**", "***") %in% unique(effectSize$stars))) {
    stop("No unconfounded associations remain with the current cutoff values. Consider manually including categorical metaVariables into the Heatmap by listing them through the 'keepMeta' argument.")
  }

  lowerLim <- min(effectSize$Ds)
  upperLim <- max(effectSize$Ds)

  if (d_range == "full") {
    lowerLim <- -1
    upperLim <- 1
  }


  signifCol <- c("gray45", "black")
  signifMeaning <- c("confounded", "deconfounded")
  legendShapes <- c(1,8)
  if (all(effectSize$status)) { # if no confounded signals in the output
    signifCol <- c("black")
    signifMeaning <- c("deconfounded")
    legendShapes <- c(8)
  }

  # include added name coluns into plots!!
  if (cuneiform) {

    # put together needed shapes and their meaning
    divShapes <- c()
    divShapesMeaning <- c()
    signs <- unique(sign(effectSize$Ds))
    if (-1 %in% signs) {
      divShapes <- c(divShapes, 25)
      divShapesMeaning <- c(divShapesMeaning, "negative association")
    }
    if (0 %in% signs) {
      divShapes <- c(divShapes, 23)
      divShapesMeaning <- c(divShapesMeaning, "no association/no data")
    }
    if (1 %in% signs) {
      divShapes <- c(divShapes, 24)
      divShapesMeaning <- c(divShapesMeaning, "positive association")
    }

    heatmapGGplot <- ggplot(effectSize, aes(x = .data$metaVariableNames, y = .data$featureNames)) +
      # do cuneiform plot with coloring based on effectsizes
      geom_point (aes (fill = .data$Ds,
                       shape = as.factor (sign (.data$Ds)),
                       color = .data$status)) +
      scale_shape_manual (name = "Direction",
                          values = divShapes,
                          labels = divShapesMeaning) +
      scale_fill_gradient2(low = d_col[1],
                           mid = d_col[2],
                           high = d_col[3],
                           midpoint = 0,
                           guide = guide_colorbar (raster = F),
                           limits = c(lowerLim,upperLim)) +
      # the shape lines color indicate confounding status
      scale_color_manual(name = "Confounding status",
                         values = signifCol,
                         labels = signifMeaning) +
      guides(#shape = FALSE,
             color = guide_legend(override.aes = list(shape  = 24))) +

      # make it pretty
      theme_classic() +
      theme(axis.text.x = element_text(size = 7,
                                       angle = 90,
                                       hjust = 1,
                                       vjust = 0.3),
            axis.text.y = element_text(size = 7,
                                       angle = 0,
                                       hjust = 1,
                                       vjust = 0.35),
            plot.title.position = "plot",
            plot.title = element_text(hjust = 0)) +
      labs(title="Summarizing cuneiform plot",
           #subtitle="FDR-values: < 0.001 = ***, < 0.01 = **, < 0.1 = * ",
           x = "Metadata variables",
           y = "Omics features")

  } else {
    heatmapGGplot <- ggplot(effectSize, aes(x = .data$metaVariableNames, y = .data$featureNames)) +
      # do the heatmap tile coloring based on effect sizes
      geom_tile(aes(fill = .data$Ds), color = tileBordCol) +
      scale_fill_gradient2(name = "effect size",
                           low = d_col[1],
                           mid = d_col[2],
                           high = d_col[3],
                           midpoint = 0,
                           guide = guide_colorbar (raster = F),
                           limits = c(lowerLim,upperLim)) +
      # add significance stars/circles for deconfounded/confounded associations
      geom_text(aes(label= .data$stars, colour = .data$status),
                size=2,
                key_glyph = "point") +
      scale_color_manual(name = "confounding status",
                         values = signifCol,
                         labels = signifMeaning) +
      guides(color = guide_legend(override.aes = list(shape = legendShapes) ) ) +

      # make it pretty
      theme_classic() +
      theme(axis.text.x = element_text(size = 7,
                                       angle = 90,
                                       hjust = 1,
                                       vjust = 0.3),
            axis.text.y = element_text(size = 7,
                                       angle = 0,
                                       hjust = 1,
                                       vjust = 0.35),
            plot.title.position = "plot",
            plot.title = element_text(hjust = 0),
            plot.subtitle=element_text(size=8)) +
      labs(title="Summarizing heatmap",
           subtitle="p.adjust-values: < 0.001 = ***, < 0.01 = **, < 0.1 = * ",
           x = "Metadata variables",
           y = "Omics features")

  }
  return(heatmapGGplot)
}

