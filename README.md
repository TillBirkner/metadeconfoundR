# metadeconfoundR

metadeconfoundR was developed to perform a confounder-aware biomarker search of cross-sectional multi-omics medical datasets. 
It first detects significant associations between individual supplied features and available metadata, using simple nonparametric tests like mann whitney u test. 
In a second step, potential confounding effects between different metadata variables are detected, using nested linear model comparison post-hoc tests.

## Instalation

```
# in R
library(devtools)
install_github("TillBirkner/metadeconfoundR")
```
