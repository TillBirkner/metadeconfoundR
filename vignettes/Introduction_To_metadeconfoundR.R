## ----style, eval=TRUE, echo=FALSE, results="asis"--------------------------
BiocStyle::latex()

## ----chu_inf_run, cache = FALSE--------------------------------------------
library(metadeconfoundR)

data(reduced_feature)
data(metaMatMetformin)

example_output <- MetaDeconfound(featureMat = reduced_feature,
	metaMat = metaMatMetformin,
	nnodes = 1)

knitr::kable(example_output$Ds[1:3, 1:5], "latex", digits = 2, booktabs = TRUE)

## ----echo = FALSE----------------------------------------------------------
	knitr::kable(summary(example_output$status), "latex", digits = 2, booktabs = TRUE)

## ----echo = FALSE----------------------------------------------------------
	knitr::kable(t(summary(as.factor(metaMatMetformin$Dataset))), "latex", digits = 2, booktabs = TRUE)
	
	for (i in unique(metaMatMetformin$Dataset)) {
	binaryDummy <- rep(0, length(metaMatMetformin$Dataset))
	binaryDummy[metaMatMetformin$Dataset == i] <- 1
	metaMatMetformin[[i]] <- binaryDummy
	colnames(metaMatMetformin)[ncol(metaMatMetformin)] <- i
	}
	metaMatMetformin$Dataset <- NULL

## ----chu_split_run, cache = FALSE------------------------------------------

example_output2 <- MetaDeconfound(featureMat = reduced_feature,
	metaMat = metaMatMetformin, 
	nnodes = 1)


knitr::kable(example_output2$Ds[1:3, 1:7], "latex", digits = 2, booktabs = TRUE)

## ----sessionInfo, results="asis"-------------------------------------------
toLatex(sessionInfo())

