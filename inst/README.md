[![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/MoEClust)](https://cran.r-project.org/package=MoEClust)
[![rstudio mirror downloads](https://cranlogs.r-pkg.org/badges/MoEClust?)](https://github.com/metacran/cranlogs.app)
[![rstudio mirror downloads](https://cranlogs.r-pkg.org/badges/grand-total/MoEClust?color=82b4e8)](https://github.com/metacran/cranlogs.app)

# MoEClust R Package
## Finite Gaussian Mixtures of Experts - 
## Parsimonious Model-Based Clustering 
## with Gating and Expert Network Covariates
### Written by Keefe Murphy

Fits _MoEClust_ models introduced by Murphy and Murphy (2017) <[arXiv:1711.05632](https://arxiv.org/abs/1711.05632)>, i.e. fits finite Gaussian mixture of experts models with gating and expert network covariates using parsimonious covariance parameterisations from __mclust__ via the EM algorithm. Also visualises mixture of experts models with parsimonious covariance structures using generalised pairs plots.

The most important function in the __MoEClust__ package is: `MoE_clust`, for fitting the model via EM with gating and/or expert network covariates, supplied via formula interfaces. Other functions also exist, e.g. `MoE_control`, `MoE_crit`, `MoE_dens`, `MoE_estep`, and `MoE_aitken`, which are all used within `MoE_clust` but are nonetheless made available for standalone use. `MoE_compare` is provided for conducting model selection between different results from `MoE_clust` using different covariate combinations &/or initialisation strategies, etc.

A dedicated plotting function exists for visualising the results using generalised pairs plots, for examining the gating network &/or log-likelihood, and/or graphing model selection criteria values. The generalised pairs plots (`MoE_gpairs`) visualise all pairwise relationships between clustered response variables and associated gating &/or expert network continuous &/or categorical variables, coloured according to the MAP classification, and also give the marginal distributions of each variable along the diagonal.

An `as.Mclust` method is provided to coerce the output of class `"MoEClust"` from `MoE_clust` to the `"Mclust"` class, to facilitate use of plotting and other functions for the `"Mclust"` class within the __mclust__ package. As per __mclust__, __MoEClust__ also facilitates modelling with an additional noise component (with or without the mixing proportion for the noise component depending on covariates).

The package also contains two data sets: `ais` and `CO2data`.

To install the development version of the package type:

```
# If required install devtools:  
# install.packages('devtools')  
devtools::install_github('Keefe-Murphy/MoEClust')
```

You can then explore the package with:

```
library(MoEClust)  
help(MoE_clust) # Help on the main modelling function
```
