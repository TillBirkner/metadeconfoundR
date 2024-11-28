## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.align = "center"
)

## ---- echo=F, fig.cap = "Figure 1: metadeconfoundR pipeline overview.", out.width = "99%", dev = "png"----
knitr::include_graphics("Figures/metadeconfoundROverview_20240605.png")

## ---- echo=F, fig.show='hold', fig.cap = "Figure 2: detailed status labelling decision tree.", out.width = "99%", out.height = "99%", dev = "png"----
knitr::include_graphics("Figures/flowChartDecision_mixed_CI.png")

## ----example, eval=FALSE, echo=TRUE-------------------------------------------
#  # CRAN
#  install.packages("metadeconfoundR")
#  # github stable version
#  devtools::install_github("TillBirkner/metadeconfoundR")
#  # github developmental version
#  devtools::install_github("TillBirkner/metadeconfoundR@develop")
#  library(metadeconfoundR)

## ----hide, eval=TRUE, echo=FALSE, results="asis", results='hide', message=FALSE----
library(metadeconfoundR)
library(ggplot2)
library(gridExtra)
library(kableExtra)

## ----showTableF, eval=TRUE,  echo=FALSE---------------------------------------
kbl(reduced_feature[10:15, 1:5], caption = "Table 1: feature input format example")

## ----showTableM, eval=TRUE, echo=FALSE----------------------------------------
kbl(metaMatMetformin[10:15, 1:5], caption = "Table 2: metadata input format example")

## ----runExample, eval=TRUE, echo=TRUE, nobreak=TRUE---------------------------
data(reduced_feature)
data(metaMatMetformin)

# check correct ordering
all(rownames(metaMatMetformin) == rownames(reduced_feature))
all(order(rownames(metaMatMetformin)) == order(rownames(reduced_feature)))

example_output <- MetaDeconfound(
  featureMat = reduced_feature,
  metaMat = metaMatMetformin,
  returnLong = TRUE,
  logLevel = "ERROR"
)

## ----runExampleRand, eval=TRUE, echo=TRUE, nobreak=TRUE, message = FALSE------
RandDataset_output <- MetaDeconfound(
  featureMat = reduced_feature,
  metaMat = metaMatMetformin,
  randomVar = c("Dataset"),
  returnLong = TRUE,
  logLevel = "ERROR"
)

## ----showTableO, eval=TRUE, echo=FALSE----------------------------------------
kbl(example_output[1:5, 1:6], caption = "Table 3: example output of MetadDeconfound()")

## ----runMediation, eval=TRUE, echo=TRUE, nobreak=TRUE-------------------------
reduced_featureMedi <- reduced_feature[, 1:40]
mediationMat <- reduced_feature[, 41:50]

example_outputMedi <- MetaDeconfound(
  featureMat = reduced_featureMedi,
  metaMat = metaMatMetformin,
  mediationMat = mediationMat,
  returnLong = TRUE,
  logLevel = "ERROR"
)

## ----runMediPlot, echo=T, eval=FALSE------------------------------------------
#  BuildHeatmap(
#    example_outputMedi,
#    keepMeta = colnames(metaMatMetformin),
#    d_range = "full"
#  ) +
#    theme(strip.background = element_rect(fill = "red"))

## ----runMediPlotHidden, echo=F, fig.cap = "Figure 3: BuildHeatmap() output for mediation analysis data", fig.width = 4.5, fig.height=6----
BuildHeatmap(
  example_outputMedi,
  keepMeta = colnames(metaMatMetformin),
  d_range = "full"
) + 
  theme(strip.background = element_rect(fill = "red"))

## ----runHeatmap, echo=T, eval=FALSE-------------------------------------------
#  left <- BuildHeatmap(example_output)
#  right <- BuildHeatmap(RandDataset_output)
#  grid.arrange(left, right, ncol = 2)

## ----runHeatmapHidden, echo=F, fig.cap = "Figure 4: default output of the BuildHeatmap() function", fig.width = 5.5, fig.height=6----
left <- BuildHeatmap(example_output) + labs(title = "example_output")
right <- BuildHeatmap(RandDataset_output)  + labs(title = "RandDataset_output")
grid.arrange(left, right, ncol = 2)

## ----runCun, echo=T, eval=FALSE-----------------------------------------------
#  BuildHeatmap(
#    example_output,
#    cuneiform = TRUE,
#    keepMeta = colnames(metaMatMetformin),
#    d_range = "full"
#  )

## ----runCunHidden, echo=F, fig.cap = "Figure 5: alternative cuneiform output of the BuildHeatmap() function", fig.width = 3.5, fig.height=6----
BuildHeatmap(example_output, cuneiform = TRUE, keepMeta = colnames(metaMatMetformin), d_range = "full")

## ----runpostCustom, echo=T, fig.cap = "Figure 6: post-plotting ggplot2 alterations", fig.width = 1.9, fig.height=6----
BuildHeatmap(example_output) +
  theme(
    legend.position = "none",
    axis.text.y = element_text(face = "italic"),
    plot.title = element_blank(),
    plot.subtitle = element_blank()
  )

## ----produceLongOut, eval=FALSE-----------------------------------------------
#  print(example_output[101:105, ])

## ----headLongInput, echo=FALSE------------------------------------------------
knitr::kable(example_output[101:105, ])

## ----runImportLongPrior-------------------------------------------------------

minQValues <- ImportLongPrior(longPrior = example_output,
                                featureMat = reduced_feature,
                                metaMat = metaMatMetformin)

## ----fakeshowminQ, eval=FALSE-------------------------------------------------
#  print(minQValues[1:5, 1:5])

## ----showminQValues, echo=FALSE-----------------------------------------------
knitr::kable(minQValues[1:5, 1:5])

## ----runInformedMetadeconf, eval=FALSE----------------------------------------
#  
#  example_output2 <- MetaDeconfound(featureMat = reduced_feature,
#                                    metaMat = metaMatMetformin,
#                                    minQValues = minQValues)

## ----runPartial, eval=TRUE, echo=TRUE, nobreak=TRUE, message=FALSE------------
ex_out_partial <- GetPartialEfSizes(
  featureMat = reduced_feature,
  metaMat = metaMatMetformin,
  metaDeconfOutput = RandDataset_output,
  randomVar = c("Dataset")
)

## ----runPartialPlotting, eval=TRUE, echo=FALSE, message=FALSE, fig.width = 4.5, fig.height=5, fig.cap="Figure 7: BuildHeatmap() plotting of partial effect sizes calculated by GetPartialEfSizes()"----

partialHM_Rand <- BuildHeatmap(ex_out_partial,
                               plotPartial = "partial",
                               reOrder = "none",
                               keepMeta = colnames(metaMatMetformin),
                               d_range = "full",
                               d_cutoff = 0.0000001
                               ) +
  labs(title = "partial", subtitle = "") +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.6),
        axis.title.x = element_blank()
        )
partialRelHM_Rand <- BuildHeatmap(ex_out_partial,
                                  plotPartial = "partialRel",
                                  reOrder = "none",
                                  keepMeta = colnames(metaMatMetformin),
                                  d_range = "full",
                                  d_cutoff = 0.0000001
                                  ) +
  labs(title = "partialRel", subtitle = "") +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.6),
        axis.title.x = element_blank()
        )
partialNormHM_Rand <- BuildHeatmap(ex_out_partial,
                                   plotPartial = "partialNorm",
                                   reOrder = "none",
                                   keepMeta = colnames(metaMatMetformin),
                                   d_range = "full",
                                   d_cutoff = 0.0000001
                                   ) +
  labs(title = "partialNorm", subtitle = "") +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.6),
        axis.title.x = element_blank()
        )

grid.arrange(
  partialHM_Rand,
  partialRelHM_Rand,
  partialNormHM_Rand,
  nrow = 1,
  widths = c(1.05, 1, 1)
)

## ----sessionInfo, results="asis", echo=FALSE----------------------------------
pander::pander(sessionInfo())

