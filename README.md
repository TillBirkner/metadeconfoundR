# metadeconfoundR

<!-- badges: start -->
[![CRAN status](https://www.r-pkg.org/badges/version/metadeconfoundR)](https://CRAN.R-project.org/package=metadeconfoundR)
[![R-CMD-check](https://github.com/TillBirkner/metadeconfoundR/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/TillBirkner/metadeconfoundR/actions/workflows/R-CMD-check.yaml)
[![codecov](https://codecov.io/github/TillBirkner/metadeconfoundR/branch/develop/graph/badge.svg?token=GGLNXAB5B2)](https://app.codecov.io/github/TillBirkner/metadeconfoundR)
![CRAN/METACRAN License](https://img.shields.io/cran/l/metadeconfoundR)
<!-- badges: end -->

metadeconfoundR was developed to perform a confounder-aware biomarker search of cross-sectional multi-omics medical datasets. 
It first detects significant associations between individual supplied features and available metadata, using simple nonparametric tests like mann whitney u test. 
In a second step, potential confounding effects between different metadata variables are detected, using nested linear model comparison post-hoc tests.

metadeconfoundR is also able to incorporate prior knowledge about confounding effects into this second analysis step. Drug association knowledge gained and reported from analyses of the [MetaCardis](https://cordis.europa.eu/project/id/305312) cohort (Forslund et al., 2021) could, for example, now be used as additional input for future studies encompassing the same omics modalities and available metadata. Now, known confounders will be treated as such even if statistical power in the new dataset is not sufficient to detect them, thereby reducing the risk of drawing wrong conclusions based on undetected confounders. Details about this can be found in the latest release notes of metadeconfoundR.

## Instalation

```R
# in R
# from CRAN:
install.packages("metadeconfoundR")
# developmental version from GitHub:
library(devtools)
install_github("TillBirkner/metadeconfoundR@develop")
```
## Documentation

See [vignette](https://htmlpreview.github.io/?https://github.com/TillBirkner/metadeconfoundR/blob/main/vignettes/Introduction_To_metadeconfoundR.html) for example code and explanations.

## Bugs/issues

Please submit bugs and problems as an issue (A template with necessary contents of issues is supplied.)
