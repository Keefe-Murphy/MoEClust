#' MoEClust: Gaussian Parsimonious Clustering Models with Covariates and a Noise Component
#'
#' Clustering via parsimonious Gaussian Mixtures of Experts using the \emph{MoEClust} models introduced by Murphy and Murphy (2020) <\doi{10.1007/s11634-019-00373-8}>. This package fits finite Gaussian mixture models with gating and/or expert network covariates using a range of parsimonious covariance parameterisations from the GPCM family via the EM/CEM algorithm. Visualisation of the results of such models using generalised pairs plots and the inclusion of an additional noise component is also facilitated.
#' @section Usage:
#' The most important function in the \pkg{MoEClust} package is: \code{\link{MoE_clust}}, for fitting the model via EM/CEM with gating and/or expert network covariates, supplied via formula interfaces.
#'
#' \code{\link{MoE_compare}} is provided for conducting model selection between different results from \code{\link{MoE_clust}} using different covariate combinations &/or initialisation strategies, etc.
#' 
#' \code{\link{MoE_stepwise}} is provided for conducting a greedy forward stepwise search to identify the optimal model in terms of the number of components, GPCM covariance type, and the subsets of gating/expert network covariates.
#' 
#' \code{\link{MoE_control}} allows supplying additional arguments to \code{\link{MoE_clust}} and \code{\link{MoE_stepwise}} which govern, among other things, controls on the inclusion of an additional noise component and controls on the initialisation of the allocations for the EM/CEM algorithm.
#'
#' A dedicated plotting function (\code{\link{plot.MoEClust}}) exists for visualising the results using generalised pairs plots, for examining the gating network, and/or log-likelihood, and/or clustering uncertainties, and/or similarity matrix, and/or graphing model selection criteria values. The generalised pairs plots (\code{\link{MoE_gpairs}}) visualise all pairwise relationships between clustered response variables and associated continuous, categorical, and/or ordinal covariates in the gating &/or expert networks, coloured according to the MAP classification, and also give the marginal distributions of each variable (incl. the covariates) along the diagonal.
#'
#' An \code{\link[=as.Mclust.MoEClust]{as.Mclust}} method is provided to coerce the output of class \code{"MoEClust"} from \code{\link{MoE_clust}} to the \code{"Mclust"} class, to facilitate use of plotting and other functions for the \code{"Mclust"} class within the \pkg{mclust} package. As per \pkg{mclust}, \pkg{MoEClust} also facilitates modelling with an additional noise component (with or without the mixing proportion for the noise component depending on covariates).
#'
#' Finally, a \code{\link[=predict.MoEClust]{predict}} method is provided for predicting the fitted response and probability of cluster membership (and by extension the MAP classification) for new data, in the form of new covariates and new response data, or new covariates only.
#'
#' Other functions also exist, e.g. \code{\link{MoE_crit}}, \code{\link{MoE_dens}}, \code{\link{MoE_estep}}, \code{\link{MoE_compare}}, and \code{\link{aitken}}, which are all used within \code{\link{MoE_clust}} but are nonetheless made available for standalone use. 
#' 
#' The package also contains two data sets: \code{ais} and \code{CO2data}.
#'
#' @section Details:
#' \describe{
#' \item{Type: }{Package}
#' \item{Package: }{MoEClust}
#' \item{Version: }{1.5.2}
#' \item{Date: }{2023-12-10 (this version), 2017-11-28 (original release)}
#' \item{Licence: }{GPL (>= 3)}
#' }
#'
#' @section See Also:
#' Further details and examples are given in the associated vignette document:\cr
#' \code{vignette("MoEClust", package = "MoEClust")}
#'
#' @author
#' Keefe Murphy [aut, cre], Thomas Brendan Murphy [ctb]
#'
#' \strong{Maintainer}: Keefe Murphy - <\email{keefe.murphy@@mu.ie}>
#' @references Murphy, K. and Murphy, T. B. (2020). Gaussian parsimonious clustering models with covariates and a noise component. \emph{Advances in Data Analysis and Classification}, 14(2): 293-325. <\doi{10.1007/s11634-019-00373-8}>.
#' @examples
#' \donttest{data(ais)
#' 
#' # Fit two sets of models
#' res1  <- MoE_clust(ais[,3:7], G=2, gating= ~ BMI, expert= ~ sex,
#'                    modelNames=c("VEE", "EVE", "VVE"), network.data=ais)
#' res2  <- MoE_clust(ais[,3:7], G=2, equalPro=TRUE, expert= ~ sex,
#'                    modelNames=c("VEE", "EVE", "VVE"), network.data=ais) 
#'         
#' # Compare the best model from each set of results
#' (comp <- MoE_compare(res1, res2, optimal.only=TRUE))
#' 
#' # Produce a plot for the optimal model                                                   
#' plot(comp$optimal, what="gpairs")
#' 
#' # Summarise its classification table, component parameters, and gating/expert networks
#' summary(comp$optimal, classification=TRUE, parameters=TRUE, networks=TRUE)
#'
#' data(CO2data)
#' CO2   <- CO2data$CO2
#' GNP   <- CO2data$GNP
#' 
#' # Fit a range of models 
#' m1    <- MoE_clust(CO2, G=1:3)
#' m2    <- MoE_clust(CO2, G=2:3, gating= ~ GNP)
#' m3    <- MoE_clust(CO2, G=1:3, expert= ~ GNP)
#' m4    <- MoE_clust(CO2, G=2:3, gating= ~ GNP, expert= ~ GNP)
#' m5    <- MoE_clust(CO2, G=2:3, equalPro=TRUE)
#' m6    <- MoE_clust(CO2, G=2:3, expert= ~ GNP, equalPro=TRUE)
#'
#' # Extract the model with highest BIC
#' (comp <- MoE_compare(m1, m2, m3, m4, m5, m6, criterion="bic"))
#'  
#' # See if a better model can be found using greedy forward stepwise selection
#' # Conduct a stepwise search on the same data
#' (mod1 <- MoE_stepwise(CO2, CO2data[,"GNP", drop=FALSE]))
#' 
#' # Conduct another stepwise search considering models with a noise component
#' (mod2 <- MoE_stepwise(CO2, CO2data[,"GNP", drop=FALSE], noise=TRUE))
#' 
#' # Compare all sets of results to choose the optimal model
#' (best <- MoE_compare(mod1, mod2, comp, pick=1)$optimal)}
#' @docType package
#' @keywords package
"_PACKAGE"

.onAttach <- function(lib, pkg) {
  path    <- file.path(lib, pkg, "DESCRIPTION")
  version <- read.dcf(path, "Version")
  name    <- read.dcf(path, "Package")
  if(interactive()) {
    packageStartupMessage(paste("\n___  ___      _____ _____ _           _   \n|  \\/  |     |  ___/ .__ \\ |         | |     Gaussian Parsimonious \n| .  . | ___ | |__ | |  \\/ |_   _ ___| |_   Clustering Models with\n| |\\/| |/ _ \\|  __|| |   | | | | / __| ._|\t  Covariates and a\n| |  | | (_) | |___| |__/\\ | |_| \\__ \\ |_\t   Noise Component\n\\_|  |_/\\___/\\____/ \\____/_|\\__,_|___/\\__|\t     version", version, "\n"))
  } else   {
    packageStartupMessage("\nPackage ", sQuote(name), " version ", version, ".\n")
  }
    packageStartupMessage(paste("Type", sQuote("?MoEClust"), "to see a brief guide to how to use this R package.\nType", sQuote(paste0("citation(", dQuote(name),")")) ,"for citing the package in publications.\nType", sQuote("MoE_news()"), "to see new features recent changes and bug fixes.\n"))
  if(interactive() &&
     name %in% utils::old.packages()[,1L]) {
    packageStartupMessage("\n   !!! A newer version of this package is available via CRAN !!!")
  }
}
