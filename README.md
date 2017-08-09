[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/MoEClust)](https://cran.r-project.org/package=MoEClust)
[![rstudio mirror downloads](http://cranlogs.r-pkg.org/badges/MoEClust?)](https://github.com/metacran/cranlogs.app)
[![rstudio mirror downloads](http://cranlogs.r-pkg.org/badges/grand-total/MoEClust?color=82b4e8)](https://github.com/metacran/cranlogs.app)

# MoEClust R Package
## Model-Based Clustering with Gating and Expert Covariates
### Written by Keefe Murphy

Fits finite Mixtures of Experts models with mclust family covariance structures using the EM algorithm, ie. allows incorporation of covariates into the mixing proportions and/or Gaussian densities of finite mixture models under the various covariance parameterisations in the mclust family.

The most important function in the __MoEClust__ package is: `MoE_clust`, for fitting the model via EM with gating and/or expert network covariates, supplied via formula interfaces. Other functions also exist, e.g. `MoE_control`, `MoE_crit`, and `MoE_dens`, which are all used within `MoE_clust` but are nonetheless made available for standalone use.

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
