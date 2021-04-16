#' BuildHeatmap
#'
#' BuildHeatmap summarizes Metadecofound output in a heatmap or cuneiform plot
#'
#' @param metaDeconfOutput output of a metadeconfound run
#' @param q_cutoff optional FDR-value cutoff used to remove
#' low-significance entries from data
#' @param d_cutoff optional effect size cutoff used to remove
#' low effect size entries from data
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
#' names used as col.names in metaDeconfOutput. These human readable
#' names will be displayed in the final plot. First column: machine-readable,
#' second column: human-readable.
#' @return ggplot2 object
#' @details for more details and explanations please see the package vignette.
#' @examples
#'data(reduced_feature)
#'data(metaMatMetformin)
#'\donttest{
#'example_output <- MetaDeconfound(featureMat = reduced_feature,
#'                                   metaMat = metaMatMetformin)
#'
#'plotObject <- BuildHeatmap(example_output)
#'print(plotObject)
#'
#'alternativePlot <- buildHeatmap(metadeconfoundR_output, coloring = 2, showConfounded = FALSE)
#'print(alternativePlot)
#'}
#'
#' @import ggplot2
#' @importFrom reshape2 melt
#' @export

BuildHeatmap <- function(metaDeconfOutput,
                         q_cutoff = 0.1,
                         d_cutoff = 0.01,
                         cuneiform = FALSE,
                         coloring = 0,
                         showConfounded = TRUE,
                         intermedData = FALSE,
                         featureNames = NULL,
                         metaVariableNames = NULL
                         ) {

  if (class(metaDeconfOutput) == "list") {
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
  } else {
    effectSize <- metaDeconfOutput[, c("feature", "metaVariable", "Ds")]
    fdr <- metaDeconfOutput[, c("feature", "metaVariable", "Qs")]
    status <- metaDeconfOutput[, c("feature", "metaVariable", "status")]

    effectSize$Ds[effectSize$Ds == Inf] <- 0
    effectSize$Ds[is.na(effectSize$Ds)] <- 0
    fdr$Qs[is.na(fdr$Qs)] <- 1
  }

  # identify all NS entries (both naive NS, and flagged "NS" by deconfounding step)
  insignificant <- unlist(lapply(strsplit(as.character(status$status), split = ", "),
                                 function(l) any(l %in% c("NS"))))
  # identify all non-confounded significant entries
  trueDeconf <- unlist(lapply(strsplit(as.character(status$status), split = ", "),
                              function(l) any(l %in% c("SD", "LD", "NC"))))


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
  #### fetching all rows/columns that have low Ds and Qs and removing them
  #
  #
  remove <- vector()
  for (i in unique(effectSize$feature)) { # identify omics features to be removed
    aGenus <- fdr[fdr$feature == i, ]
    aGenusD <- effectSize[effectSize$feature == i, ]
    if (sum(na.exclude(abs(aGenus$Qs)) > q_cutoff) == length(na.exclude(aGenus$Qs)) ||
        sum(na.exclude(abs(aGenusD$Ds)) < d_cutoff) == length(na.exclude(aGenusD$Ds))
    ) {
      remove <- c(remove, which(effectSize$feature == i))
    }
  }

  remove <- unique(remove)
  if (length(remove) > 0) {
    effectSize <- effectSize[-remove, ]
  }
  effectSize <- droplevels(effectSize)



  remove_metavariables <- vector()
  for (i in unique(effectSize$metaVariable)) { # identify metavariables to be removed

    aMetaVariable <- fdr[fdr$metaVariable == i, ]
    aMetaVariableD <- effectSize[effectSize$metaVariable == i, ]
    if (sum(na.exclude(abs(aMetaVariable$Qs)) > q_cutoff) == length(na.exclude(aMetaVariable$Qs)) ||
        sum(na.exclude(abs(aMetaVariableD$Ds)) < d_cutoff) == length(na.exclude(aMetaVariableD$Ds))
    ) {
      remove <- c(remove, which(effectSize$metaVariable == i))
      remove_metavariables <- c(remove_metavariables, i)
    }
  }

  # remove filtered out matavariables from the dataset
  effectSize <- effectSize[!(effectSize$metaVariable %in% remove_metavariables), ]

  # cluster heatmap by reordering the factor levels for both dimensions of the heatmap
  eff_cast <- reshape2::dcast(effectSize, effectSize[[1]]~metaVariable, value.var = "Ds")

  rownames(eff_cast) <- eff_cast[[1]]
  eff_cast[[1]] <- NULL # move feature names to rownames
  ord <- hclust(dist(eff_cast, method = "euclidean"), method = "ward.D")$order
  eff_cast <- scale(t(eff_cast))
  ord2 <- hclust(dist(eff_cast, method = "euclidean"), method = "ward.D")$order
  effectSize$metaVariable <- droplevels(effectSize$metaVariable)

  effectSize$feature <- factor(as.factor(effectSize$feature),
                               levels = levels(as.factor(effectSize$feature)) [ord])
  effectSize$metaVariable <- factor(as.factor(effectSize$metaVariable),
                                    levels = levels(as.factor(effectSize$metaVariable)) [ord2])


  effectSize$featureNames <- effectSize$feature
  effectSize$metaVariableNames <- effectSize$metaVariable

  if (!is.null(featureNames)) {
    if (length(unique(featureNames[[2]])) != ncol(featureNames)) {
      featureNames[[2]] <- make.unique(featureNames[[2]])
      warning('non-unique human-readable feature names where made unique using base::make.unique')
    }
    map = stats::setNames(featureNames[[2]], featureNames[[1]])
    effectSize$featureNames <- map[as.vector(effectSize$feature)]


    effectSize$featureNames <- factor(as.factor(effectSize$featureNames),
                                      levels = map[levels(effectSize$feature)])
  }

  if (!is.null(metaVariableNames)) {
    map = stats::setNames(metaVariableNames[[2]], metaVariableNames[[1]])
    effectSize$metaVariableNames <- map[as.vector(effectSize$metaVariable)]

      effectSize$metaVariableNames <- factor(as.factor(effectSize$metaVariableNames),
                                        levels = map[levels(effectSize$metaVariable)])
  }

  if (intermedData == TRUE) {
    return(effectSize)
  }

  lowerLim <- min(effectSize$Ds)
  upperLim <- max(effectSize$Ds)

  # include added name coluns into plots!!
  if (cuneiform) {
    heatmapGGplot <- ggplot(effectSize, aes(x = metaVariable, y = feature)) +
      # do cuneiform plot with coloring based on effectsizes
      geom_point (aes (fill = Ds,
                       shape = as.factor (sign (Ds)),
                       color = status)) +
      scale_shape_manual (name = "Direction",
                          values = c (25, 24),
                          labels = c("negative association",
                                     "positive association")) +
      scale_fill_gradient2(low = "red",
                           mid = "white",
                           high = "blue",
                           midpoint = 0,
                           guide = guide_colorbar (raster = F),
                           limits = c(lowerLim,upperLim)) +
      # the shape lines color indicate confounding status
      scale_color_manual(name = "Confounding status",
                         values = c("gray45", "black"),
                         labels = c("confounded", "deconfounded")) +
      guides(shape = FALSE,
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
                                       vjust = 0.35)) +
      labs(title="MetaDeconfoundR summarizing coneiform plot",
           subtitle="FDR-values: < 0.001 = ***, < 0.01 = **, < 0.1 = * ",
           x = "Metadata variables",
           y = "Omics features")

  } else {
    heatmapGGplot <- ggplot(effectSize, aes(x = metaVariableNames, y = featureNames)) +
      # do the heatmap tile coloring based on effect sizes
      geom_tile(aes(fill = Ds)) +
      scale_fill_gradient2(low = "red",
                           mid = "white",
                           high = "blue",
                           midpoint = 0,
                           guide = guide_colorbar (raster = F),
                           limits = c(lowerLim,upperLim)) +
      # add significance stars/circles for deconfounded/confounded associations
      geom_text(aes(label=stars, colour = status),
                size=2,
                key_glyph = "point") +
      scale_color_manual(name = "Confounding status",
                         values = c("gray45", "black"),
                         labels = c("confounded", "deconfounded"), ) +
      guides(color = guide_legend(override.aes = list(shape = c(1,8)) ) ) +

      # make it pretty
      theme_classic() +
      theme(axis.text.x = element_text(size = 7,
                                       angle = 90,
                                       hjust = 1,
                                       vjust = 0.3),
            axis.text.y = element_text(size = 7,
                                       angle = 0,
                                       hjust = 1,
                                       vjust = 0.35)) +
      labs(title="MetaDeconfoundR summarizing heatmap",
           subtitle="FDR-values: < 0.001 = ***, < 0.01 = **, < 0.1 = * ",
           x = "Metadata variables",
           y = "Omics features")

  }
  return(heatmapGGplot)
}

