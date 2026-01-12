#' BuildConfounderMap
#'
#' BuildConfounderMap summarizes confounder analysis of a \link[metadeconfoundR]{MetaDeconfound} output in a circle plot

#' @author Kilian Dahm
#' @param metaDeconfOutput output of a metadeconfound run
#' @param q_cutoff optional FDR-value cutoff used to remove
#' low-significance entries from data
#' @param featureColor optional vector of colors named after each unique feature in metaDeconfOutput
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
#' @param d_col set color range for effect size as c(minimum, middle, maximum),
#' default c("red", "white", "blue")
#' @param d_range range of effect size colors shown; "full": (default) range from
#' -1 to +1 (best for comparability between multiple plots);
#' "fit": range reduced according to maximum and minimum effect size
#' present in resulting plot (better color resolution for weaker effects)
#' @param trusted character vector of confounding status labels to be treated
#' as trustworthy, not-confounded signal. default = c("OK_sd", "OK_nc", "OK_d", "AD")
#' @return list of ggplot2 objects
#' @details for more details and explanations please see the package vignette.
#'
#' @examples
#'data(reduced_feature)
#'data(metaMatMetformin)
#'\donttest{
#'example_output <- MetaDeconfound(featureMat = reduced_feature,
#'                                   metaMat = metaMatMetformin,
#'                                   logLevel = "ERROR")
#'
#'plotObject <- BuildConfounderMap(example_output)
#'library(ggraph)
#'plotObject$MS0001
#'}
#'
#' @import ggplot2
#' @import ggraph
#' @importFrom ggraph guide_edge_colourbar scale_edge_colour_gradientn
#' @import igraph
#' @importFrom scales rescale
#' @importFrom stats setNames complete.cases
#' @importFrom dplyr filter pull mutate
#' @importFrom reshape2 melt
#' @import magrittr
#' @importFrom circlize colorRamp2
#' @export

BuildConfounderMap <- function(metaDeconfOutput,
                               q_cutoff = 0.1,
                               featureColor = c("black"),
                               featureNames = NULL,
                               metaVariableNames = NULL,
                               d_col = c("blue", "white", "red"),
                               d_range = "full",
                               trusted = c("OK_sd", "OK_nc", "OK_d", "AD")) {
  if (is(metaDeconfOutput, "list")) {
    long_out <- reshape2::melt(
      metaDeconfOutput$Ps,
      varnames = c("feature", "metaVariable"),
      value.name = "Ps"
    )
    long_out$Qs <- reshape2::melt(metaDeconfOutput$Qs)[, 3]
    long_out$Ds <- reshape2::melt(metaDeconfOutput$Ds)[, 3]
    long_out$status <- reshape2::melt(metaDeconfOutput$status)[, 3]

  } else {
    long_out <- metaDeconfOutput
  }

  if (length(d_col) != 3) {
    stop("wrong number of colors in d_col!\nSupply colors for c(min, middle, max)!")
  }

  if (length(trusted) == 0) {
    stop('"trusted" must contain at least one trusted status label')
  }

  allLables <- c ("OK_sd", "OK_nc", "OK_d", "AD", "NS")
  notTrusted <- allLables[!(allLables %in% trusted)]

  if (length(featureColor) == 1 |
      length(featureColor) == length(unique(long_out$feature))) {
    if (is.null(names(featureColor)) &
        length(featureColor) == length(unique(long_out$feature))) {
      warning("No names assigned to featureColor. Unique feature names are selected instead.")
      names(featureColor) <- unique(long_out$feature)
    }
  } else
    stop(
      "featureColor should either have length == 1 or the same length as unique features available in metaDeconfOutput!"
    )

  if (!is.null(featureNames)) {
    if (!is(featureNames, "data.frame")) {
      warning('class(featureNames) was coerced to "data.frame"')
      featureNames <- as.data.frame(featureNames)
    }
    if (length(unique(featureNames[[2]])) != length(featureNames[[2]])) {
      featureNames[[2]] <- make.unique(featureNames[[2]])
      warning(
        'non-unique human-readable feature names where made unique using base::make.unique'
      )
    }
    map = stats::setNames(featureNames[[2]], featureNames[[1]])
    long_out$feature <- map[as.vector(long_out$feature)]


    long_out$feature <- factor(as.factor(long_out$feature), levels = map)
  }

  if (!is.null(metaVariableNames)) {
    if (!is(metaVariableNames, "data.frame")) {
      warning('class(metaVariableNames) was coerced to "data.frame"')
      metaVariableNames <- as.data.frame(metaVariableNames)
    }
    if (length(unique(metaVariableNames[[2]])) != length(metaVariableNames[[2]])) {
      metaVariableNames[[2]] <- make.unique(metaVariableNames[[2]])
      warning(
        'non-unique human-readable metaVariable names where made unique using base::make.unique'
      )
    }
    map = stats::setNames(metaVariableNames[[2]], metaVariableNames[[1]])
    long_out$metaVariable <- map[match(as.vector(long_out$metaVariable), names(map))]

    long_out$metaVariable <- factor(as.factor(long_out$metaVariable), levels = map)
  }

  # create table of edges
  # edges <- data.frame(from = "origin",
  #                     to = unique(long_out[long_out$Qs < q_cutoff, ]$metaVariable) %>%
  #                       as.character()) %>%
  #   .[complete.cases(.), ]
  edges <- data.frame(
    from = "origin",
    to = as.character(unique(long_out[long_out$Qs < q_cutoff, ]$metaVariable))
  )
  edges <- edges[complete.cases(edges), ]

  uniqueNames <- unique(long_out$feature) #TB20250123
  names(uniqueNames) <- uniqueNames #TB20250123

  # plot_list <- lapply(unique(long_out$feature), function(x) {
  #   tmp <- long_out %>% dplyr::filter(Qs < q_cutoff & feature == x)
  plot_list <- lapply(uniqueNames, function(x) { #TB20250123
    tmp <- long_out %>% dplyr::filter(.data$Qs < q_cutoff & .data$feature == x & !(.data$status %in% notTrusted)) #TB20250123
    if (nrow(tmp) == 0) {#TB20250123
      return(NULL)#TB20250123
    }#TB20250123
    if (nrow(tmp %>% dplyr::filter(grepl("C:", .data$status))) >= 1) {
      tmp_edge <- lapply(tmp %>% dplyr::filter(grepl("C:", .data$status)) %>% dplyr::pull(.data$status), function(y) {
        unlist(strsplit(y, split = ": "))[2]
      }) %>% base::invisible()
      names(tmp_edge) <- tmp %>% dplyr::filter(grepl("C:", .data$status)) %>% dplyr::pull(.data$metaVariable)
      edge_df <- data.frame()

      lapply(names(tmp_edge), function(x) {
        tmp <- unlist(strsplit(tmp_edge[[x]], split = ", "))
        tmp_df <- data.frame(meta = rep(x, length(tmp)), confounder = tmp)
        if (nrow(edge_df) == 0)
          edge_df <<- tmp_df
        else
          edge_df <<- rbind(edge_df, tmp_df)
      }) %>% base::invisible()

      # create a data frame with connection between leaves (individuals)
      connect <- data.frame(from = edge_df$confounder , to = edge_df$meta)
      if (!is.null(metaVariableNames))
        connect$from <- map[match(connect$from, names(map))]
      connect$value <- 1
    }

    # create table of vertices
    vertices  <-  data.frame(name = unique(c(
      as.character(edges$from), as.character(edges$to)
    )), value = 1)
    vertices <- vertices %>% dplyr::mutate(
      group = edges$from[match(vertices$name, edges$to)],
      effectSize = tmp[match(vertices$name, tmp$metaVariable), ]$Ds,
      stroke = sapply(tmp[match(vertices$name, tmp$metaVariable), ]$status, function(x) {
        if (is.na(x))
          0.5
        else if (!x %in% notTrusted &
                 !grepl("C:", x))
          1.25
        else
          0.5
      }),
      id = NA
    )

    # orientate labels
    leaves <- which(is.na(match(vertices$name, edges$from)))
    nLeaves <- length(leaves)
    vertices$id[leaves] <- seq(1:nLeaves)
    vertices$angle <-  90 - 360 * vertices$id / nLeaves + 90 / nLeaves
    vertices$hjust <- ifelse(vertices$angle < -90, 1, 0)
    vertices$angle <- ifelse(vertices$angle < -90, vertices$angle + 180, vertices$angle)

    # transform into a graph object
    graph <- igraph::graph_from_data_frame(edges, vertices = vertices)

    if (nrow(tmp %>% dplyr::filter(grepl("C:", .data$status))) >= 1) {
      # prepare for get_con()
      from  <-  match(connect$from, vertices$name)
      to  <-  match(connect$to, vertices$name)

      # create color palettes
      if (length(featureColor) == 1) {
        color_graph_fun <- circlize::colorRamp2(c(0, 1), c(featureColor, "#F5F5F5"))
      } else {
        color_graph_fun <- circlize::colorRamp2(c(0, 1), c(featureColor[names(featureColor) == x], "#F5F5F5"))
      }
      color_graph <- color_graph_fun(seq(0, 1, length.out = 100))
      names(color_graph) <- seq(0, 1, length.out = 100)
    }


    color_node_fun <- circlize::colorRamp2(c(-1, 0, 1), d_col)
    color_node <- color_node_fun(seq(-1, 1, length.out = 100))
    names(color_node) <- seq(-1, 1, length.out = 100)

    # message(paste0("Plotting ConfounderMap for feature: ", x)) #TB20250123

    # generate plot
    graph_layout_original <- ggraph::create_layout(graph = graph,
                                           layout = "dendrogram",
                                           circular = TRUE)
    # order the nodes properly
    graph_layout <- graph_layout_original %>% dplyr::filter(.data$leaf == T) %>% dplyr::mutate(y_old = .data$y)
    # graph_layout$section <- sapply(graph_layout$x, function(x)
    #   ifelse(sign(x) == 1, 1, 2))

    graph_layout$y <- jitter(graph_layout$y, 0.01)
    graph_layout$x <- jitter(graph_layout$x, 0.01)
    graph_layout$y_old <- graph_layout$y

    graph_layout$section <- 1
    graph_layout$section[sign(graph_layout$x) != 1] <- 2

    graph_layout$y <- c(sapply(unique(graph_layout$section), function(x) {
      if (x == 1)
        graph_layout[graph_layout$section == x, ]$y[order(graph_layout[graph_layout$section == x, ]$y, decreasing = T)]
      else
        graph_layout[graph_layout$section == x, ]$y[order(graph_layout[graph_layout$section == x, ]$y, decreasing = F)]
    })) %>% unlist()
    graph_layout$x <- graph_layout[match(graph_layout$y, graph_layout$y_old), ]$x
    graph_layout$y_old <- NULL
    graph_layout$section <- NULL
    graph_layout_new <- rbind(graph_layout_original[1, ], graph_layout)
    attributes(graph_layout_new) <- attributes(graph_layout_original)

    lowerLim <- min(vertices$effectSize)
    upperLim <- max(vertices$effectSize)

    if (d_range == "full") {
      lowerLim <- -1
      upperLim <- 1
    }

    # plot the data
    if (nrow(tmp %>% dplyr::filter(grepl("C:", .data$status))) >= 1) {
      plot <- ggraph::ggraph(graph = graph_layout_new) +
        ggraph::geom_conn_bundle(
          data = ggraph::get_con(from = from, to = to),
          alpha = 0.9,
          #aes(colour = ..index..) ,
          aes(colour = ggplot2::after_stat(.data$index)) , #TB20250123
          width = 0.5,
          tension = 0.5
        ) +
        ggraph::scale_edge_colour_gradientn(
          name = "confounding direction",
          colors = color_graph,
          values = names(color_graph),
          breaks = c(0, 1),
          labels = c("confounding", "confounded")
        ) +
        ggraph::geom_node_point(
          aes(
            filter = .data$leaf,
            x = .data$x * 1.05,
            y = .data$y * 1.05,
            fill = .data$effectSize,
            stroke = .data$stroke
          ),
          shape = 21,
          color = "black",
          size = 4
        ) +
        scale_fill_gradient2(
          low = d_col[1],
          mid = d_col[2],
          high = d_col[3],
          n.breaks = 5,
          na.value = "white",
          name = "effect size",
          limits= c(lowerLim,upperLim)
        ) +
        ggraph::geom_node_text(
          aes(
            x = .data$x * 1.15,
            y = .data$y * 1.15,
            filter = .data$leaf,
            label = .data$name,
            angle = .data$angle,
            hjust = .data$hjust
          ),
          size = 3,
          alpha = 1
        ) +
        ggtitle(x) +
        theme_void() +
        theme(
          plot.margin = unit(c(0, 0, 0, 0), "cm"),
          plot.title = element_text(hjust = 0.5),
          legend.position = "bottom",
          legend.box = "vertical"
        ) +
        expand_limits(x = c(-2, 2), y = c(-2, 2)) +
        continuous_scale(
          aesthetics = "stroke",
          name = "confounding status",
          #name = "confounding status",
          palette = function(x) {
            scales::rescale(x, c(0.5, 1.25))
          },
          breaks = c(0.5, 1.25),
          labels = c("confounded", "deconfounded")
        ) +
        coord_fixed()
      plot
      #plot(plot)

      #ggsave(paste0("/data/Scripts/plot/Associations_", x, ".pdf")) #-->#TB20250123

    } else {
      plot <- ggraph(graph = graph_layout_new) +
        ggraph::geom_node_point(
          aes(
            filter = .data$leaf,
            x = .data$x * 1.05,
            y = .data$y * 1.05,
            fill = .data$effectSize,
            stroke = .data$stroke
          ),
          shape = 21,
          color = "black",
          size = 4
        ) +
        scale_fill_gradient2(
          low = d_col[1],
          mid = d_col[2],
          high = d_col[3],
          n.breaks = 5,
          na.value = "white",
          name = "effect size",
          limits= c(lowerLim,upperLim)
        ) +
        ggraph::geom_node_text(
          aes(
            x = .data$x * 1.15,
            y = .data$y * 1.15,
            filter = .data$leaf,
            label = .data$name,
            angle = .data$angle,
            hjust = .data$hjust
          ),
          size = 3,
          alpha = 1
        ) +
        ggtitle(x) +
        theme_void() +
        theme(
          plot.margin = unit(c(0, 0, 0, 0), "cm"),
          plot.title = element_text(hjust = 0.5),
          legend.position = "bottom",
          legend.box = "vertical"
        ) +
        expand_limits(x = c(-2, 2), y = c(-2, 2)) +
        continuous_scale(
          aesthetics = "stroke",
          name = "confounding status",
          #name = "confounding status",
          palette = function(x) {
            scales::rescale(x, c(0.5, 1.25))
          },
          breaks = c(0.5, 1.25),
          labels = c("confounded", "deconfounded")
        ) +
        coord_fixed()
      plot
      #plot(plot)
    }
  })
  #plot_list <- plot_list[sapply(plot_list, function(x) !is.null(x))]
  non_null_indices <- unlist(sapply(plot_list, function(x) !is.null(x) && length(x) > 0))
  plot_list <- plot_list[non_null_indices]
  return(plot_list)
}
