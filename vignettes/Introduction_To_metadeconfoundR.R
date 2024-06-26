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
#  # github
#  devtools::install_github("TillBirkner/metadeconfoundR")
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

## ----runExampleRand, eval=TRUE, echo=TRUE, nobreak=TRUE-----------------------
RandDataset_output <- MetaDeconfound(
  featureMat = reduced_feature,
  metaMat = metaMatMetformin,
  randomVar = c("Dataset"),
  returnLong = TRUE,
  logLevel = "ERROR"
)

## ----showTableO, eval=TRUE, echo=FALSE----------------------------------------
kbl(example_output[1:5, 1:6], caption = "Table 3: example output of MetadDeconfound()")

## ----runHeatmap, echo=T, eval=FALSE-------------------------------------------
#  left <- BuildHeatmap(example_output)
#  right <- BuildHeatmap(RandDataset_output)
#  grid.arrange(left, right, ncol = 2)

## ----runHeatmapHidden, echo=F, fig.cap = "Figure 3: default output of the BuildHeatmap() function", fig.width = 5.5, fig.height=6----
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

## ----runCunHidden, echo=F, fig.cap = "Figure 4: alternative cuneiform output of the BuildHeatmap() function", fig.width = 3.5, fig.height=6----
BuildHeatmap(example_output, cuneiform = TRUE, keepMeta = colnames(metaMatMetformin), d_range = "full")

## ----runpostCustom, echo=T, fig.cap = "Figure 5: post-plotting ggplot2 alterations", fig.width = 1.9, fig.height=6----
BuildHeatmap(example_output) +
  theme(
    legend.position = "none",
    axis.text.y = element_text(face = "italic"),
    plot.title = element_blank(),
    plot.subtitle = element_blank()
  )

## ----sessionInfo, results="asis", echo=FALSE----------------------------------
pander::pander(sessionInfo())

