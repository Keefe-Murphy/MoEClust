[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/MoEClust)](https://cran.r-project.org/package=MoEClust)
[![rstudio mirror downloads](http://cranlogs.r-pkg.org/badges/MoEClust?)](https://github.com/metacran/cranlogs.app)
[![rstudio mirror downloads](http://cranlogs.r-pkg.org/badges/grand-total/MoEClust?color=82b4e8)](https://github.com/metacran/cranlogs.app)

# MoEClust R Package
## Model-Based Clustering with Gating and Expert Covariates
### Written by Keefe Murphy

Fits finite Mixtures of Experts models with __mclust__ family covariance structures using the EM algorithm, ie. allows incorporation of covariates into the mixing proportions and/or Gaussian densities of finite mixture models under the various covariance parameterisations in the __mclust__ family. It also visualises mixture of experts models with parsimonious covariance parameterisations using generalised pairs plots.

The most important function in the __MoEClust__ package is: `MoE_clust`, for fitting the model via EM with gating and/or expert network covariates, supplied via formula interfaces. Other functions also exist, e.g. `MoE_control`, `MoE_crit`, `MoE_dens`, `MoE_estep`, and `MoE_aitken`, which are all used within `MoE_clust` but are nonetheless made available for standalone use. `MoE_compare` is provided for conducting model selection between different results from `MoE_clust` using different covariate combinations &/or initialisation strategies, etc.

A dedicated plotting function exists for visualising the results using generalised pairs plots, for examining the gating network, and/or graphing model selection criteria values.

An `as.Mclust` method is provided to coerce the output of class `"MoEClust"` from `MoE_clust` to the `"Mclust"` class, to facilitate use of plotting and other functions for the `"Mclust"` class within the __mclust__ package. As per __mclust__, __MoEClust__ also facilitates modelling with an additional noise component.

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
