[![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/MoEClust)](https://cran.r-project.org/package=MoEClust)
[![rstudio mirror downloads](https://cranlogs.r-pkg.org/badges/MoEClust?)](https://github.com/r-hub/cranlogs.app)
[![rstudio mirror downloads](https://cranlogs.r-pkg.org/badges/grand-total/MoEClust?color=82b4e8)](https://github.com/r-hub/cranlogs.app)

# MoEClust R Package
## Gaussian Parsimonious Clustering Models
## with Gating and Expert Network Covariates
## and a Noise Component
### Written by Keefe Murphy

## Description

Fits _MoEClust_ models introduced by Murphy and Murphy (2020) <[doi:10.1007/s11634-019-00373-8](https://doi.org/10.1007/s11634-019-00373-8)>, i.e. fits finite Gaussian mixture of experts models with gating and/or expert network covariates supplied via formula interfaces using a range of parsimonious covariance parameterisations from the GPCM family via the EM/CEM algorithm. Visualisation of the results of such models using generalised pairs plots and the inclusion of an additional noise component is also facilitated.

The most important function in the __MoEClust__ package is: `MoE_clust`, for fitting the model via EM/CEM with gating and/or expert network covariates, supplied via formula interfaces. `MoE_compare` is provided for conducting model selection between different results from `MoE_clust` using different covariate combinations &/or initialisation strategies, etc. 

`MoE_stepwise` is provided for conducting a greedy forward stepwise search to identify the optimal model in terms of the number of components, GPCM covariance type, and the subsets of gating/expert network covariates.

`MoE_control` allows supplying additional arguments to `MoE_clust` and `MoE_stepwise` which govern, among other things, controls on the inclusion of an additional noise component and controls on the initialisation of the allocations for the EM/CEM algorithm.

A dedicated plotting function exists for visualising the results using generalised pairs plots, for examining the gating network, and/or log-likelihood, and/or clustering uncertainties, and/or graphing model selection criteria values. The generalised pairs plots (`MoE_gpairs`) visualise all pairwise relationships between clustered response variables and associated continuous, categorical, and/or ordinal covariates in the gating &/or expert networks, coloured according to the MAP classification, and also give the marginal distributions of each variable (incl. the covariates) along the diagonal.

An `as.Mclust` method is provided to coerce the output of class `"MoEClust"` from `MoE_clust` to the `"Mclust"` class, to facilitate use of plotting and other functions for the `"Mclust"` class within the __mclust__ package. As per __mclust__, __MoEClust__ also facilitates modelling with an additional noise component (with or without the mixing proportion for the noise component depending on covariates). Finally, a `predict` method is provided for predicting the fitted response and probability of cluster membership (and by extension the MAP classification) for new data, in the form of new covariates and new response data, or new covariates only.

Other functions also exist, e.g. `MoE_crit`, `MoE_dens`, `MoE_estep`, and `aitken`, which are all used within `MoE_clust` but are nonetheless made available for standalone use. The package also contains two data sets: `ais` and `CO2data`.

## Installation

You can install the latest stable official release of the `MoEClust` package from CRAN:

```
install.packages("MoEClust")
```

or the development version from GitHub:

```
# If required install devtools:  
# install.packages('devtools')  
devtools::install_github('Keefe-Murphy/MoEClust')
```

In either case, you can then explore the package with:

```
library(MoEClust)  
help(MoE_clust) # Help on the main modelling function
```

For a more thorough intro, the vignette document is available as follows:

```
vignette("MoEClust", package="MoEClust")
```

However, if the package is installed from GitHub, the vignette is not automatically created. It can be accessed when installing from GitHub with the code:

```
devtools::install_github('Keefe-Murphy/MoEClust', build_vignettes = TRUE)
```

Alternatively, the vignette is available on the package's CRAN page.

### References
Murphy, K. and Murphy, T. B. (2020). Gaussian parsimonious clustering models with covariates and a noise component. _Advances in Data Analysis and Classification_, 14(2): 293--325. <[doi:10.1007/s11634-019-00373-8](https://doi.org/10.1007/s11634-019-00373-8)>.
