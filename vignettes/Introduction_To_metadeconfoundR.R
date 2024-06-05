## ----style, eval=TRUE, echo=FALSE, results="asis"--------------------------
BiocStyle::latex()

## ----example, eval=FALSE, echo=TRUE, results="asis"------------------------
#  library(devtools)
#  install_github("TillBirkner/metadeconfoundR")
#  library(metadeconfoundR)

## ----hide, eval=TRUE, echo=FALSE, results="asis"---------------------------
library(metadeconfoundR)
library(ggplot2)
#library(cowplot)
library(gridExtra)
library(kableExtra)

## ----runExample, eval=TRUE, echo=TRUE, tidy=TRUE, tidy.opts = list(width.cutoff = 50), cache = TRUE----
		data(reduced_feature)
		data(metaMatMetformin)

		# check correct ordering
		all(rownames(metaMatMetformin) == rownames(reduced_feature))
		all(order(rownames(metaMatMetformin)) == order(rownames(reduced_feature)))

		example_output <- MetaDeconfound(featureMat = reduced_feature,
		metaMat = metaMatMetformin,
		returnLong = TRUE)


## ----runExampleRand, eval=TRUE, echo=TRUE, tidy=TRUE, tidy.opts = list(width.cutoff = 50), cache = TRUE----

		RandDataset_output <- MetaDeconfound(featureMat = reduced_feature,
		metaMat = metaMatMetformin,
		randomVar = c("Dataset"),
		returnLong = TRUE
		)

## ----runHeatmap, eval=FALSE, echo=TRUE, tidy=TRUE, tidy.opts = list(width.cutoff = 50)----
#  			left <- BuildHeatmap(example_output)
#  			right <- BuildHeatmap(RandDataset_output)
#  			#plot_grid(left, right)
#  			grid.arrange(left, right)

## ----runHeatmapHidden, eval=TRUE, echo=FALSE, tidy=TRUE, fig.width = 5.5, tidy.opts = list(width.cutoff = 50)----
left <- BuildHeatmap(example_output)
right <- BuildHeatmap(RandDataset_output)
#plot_grid(left, right)
grid.arrange(left, right)


## ----runCun, eval=FALSE, echo=TRUE, tidy=TRUE, tidy.opts = list(width.cutoff = 40)----
#  			BuildHeatmap(example_output,
#  			cuneiform = TRUE,
#  			keepMeta = colnames(example_output$status),
#  			d_range = "full")

## ----runCunHidden, eval=TRUE, echo=FALSE, tidy=TRUE, fig.width = 3.8, fig.height = 6.5, tidy.opts = list(width.cutoff = 50)----
BuildHeatmap(example_output, cuneiform = TRUE, keepMeta = colnames(example_output$status), d_range = "full")

## ----runpostCustom, eval=TRUE, echo=TRUE, tidy=TRUE, fig.width = 1.9, fig.height = 5, tidy.opts = list(width.cutoff = 40)----
BuildHeatmap(example_output) +
theme(legend.position = "none",
axis.text.y = element_text(face = "italic"),
plot.title = element_blank(),
plot.subtitle = element_blank())

## ----sessionInfo, results="asis"-------------------------------------------
toLatex(sessionInfo())

