#' MoEClust: Gaussian Parsimonious Clustering Models with Covariates
#'
#' Clustering via parsimonious Gaussian Mixtures of Experts using the \emph{MoEClust} models introduced by Murphy and Murphy (2017) <\href{https://arxiv.org/abs/1711.05632}{arXiv:1711.05632}>. This package fits finite Gaussian mixture models with gating and/or expert network covariates using a range of parsimonious covariance parameterisations via the EM algorithm. Visualisation of the results of such models using generalised pairs plots is also facilitated.
#' @section Usage:
#' The most important function in the \pkg{MoEClust} package is: \code{\link{MoE_clust}}, for fitting the model via EM with gating and/or expert network covariates, supplied via formula interfaces.
#'
#' Other functions also exist, e.g. \code{\link{MoE_control}}, \code{\link{MoE_crit}}, \code{\link{MoE_dens}}, \code{\link{MoE_estep}}, \code{\link{MoE_compare}}, and \code{\link{aitken}}, which are all used within \code{\link{MoE_clust}} but are nonetheless made available for standalone use.
#'
#' \code{\link{MoE_compare}} is provided for conducting model selection between different results from \code{\link{MoE_clust}} using different covariate combinations &/or initialisation strategies, etc.
#'
#' A dedicated plotting function exists for visualising the results using generalised pairs plots, for examining the gating network, and/or log-likelihood, and/or clustering uncertaines, and/or graphing model selection criteria values. The generalised pairs plots (\code{\link{MoE_gpairs}}) visualise all pairwise relationships between clustered response variables and associated continuous, categorical, and/or ordinal covariates in the gating &/or expert networks, coloured according to the MAP classification, and also give the marginal distributions of each variable (incl. the covariates) along the diagonal.
#'
#' An \code{as.Mclust} method is provided to coerce the output of class \code{"MoEClust"} from \code{\link{MoE_clust}} to the \code{"Mclust"} class, to facilitate use of plotting and other functions for the \code{"Mclust"} class within the \pkg{mclust} package. As per \pkg{mclust}, \pkg{MoEClust} also facilitates modelling with an additional noise component (with or without the mixing proportion for the noise component depending on covariates).
#'
#' Finally, a \code{predict} method is provided for predicting the fitted response and probability of cluster membership (and by extension the MAP classification) for new data, in the form of new covariates and new response data, or new covariates only.
#'
#' The package also contains two data sets: \code{ais} and \code{CO2data}.
#'
#' @section Details:
#' \itemize{
#' \item{Type: }{Package}
#' \item{Package: }{MoEClust}
#' \item{Version: }{1.2.0}
#' \item{Date: }{2018-08-24 (this version), 2017-11-28 (original release)}
#' \item{Licence: }{GPL (>=2)}
#' }
#'
#' @section See Also:
#' Further details and examples are given in the associated vignette document:\cr
#' \code{vignette("MoEClust", package = "MoEClust")}
#'
#' @author
#' Keefe Murphy [aut, cre], Thomas Brendan Murphy [ctb]
#'
#' \strong{Maintainer}: Keefe Murphy - <\email{keefe.murphy@@ucd.ie}>
#' @references K. Murphy and T. B. Murphy (2017). Parsimonious Model-Based Clustering with Covariates. \emph{To appear}. <\href{https://arxiv.org/abs/1711.05632}{arXiv:1711.05632}>.
#' @examples
#' \dontrun{
#' data(ais)
#' res <- MoE_clust(ais[,3:7], G=2, gating=~BMI, expert=~sex,
#'                  modelNames=c("EVE", "VVE", "VEE"), network.data=ais)
#' plot(res, what="gpairs")
#'
#' data(CO2data)
#' GNP   <- CO2data[,1]
#' CO2   <- CO2data[,2]
#' m1    <- MoE_clust(CO2, G=1:2)
#' m2    <- MoE_clust(CO2, G=2, gating= ~ GNP)
#' m3    <- MoE_clust(CO2, G=1:2, expert= ~ GNP)
#' m4    <- MoE_clust(CO2, G=2, gating= ~ GNP, expert= ~ GNP)
#' MoE_compare(m1, m2, m3, m4)}
#' @docType package
#' @keywords package
"_PACKAGE"

.onAttach <- function(lib, pkg) {
  version <- read.dcf(file.path(lib, pkg, "DESCRIPTION"), "Version")
  if(interactive()) {
    packageStartupMessage(paste("\n___  ___      _____ _____ _           _   \n|  \\/  |     |  ___/  __ \\ |         | |     Gaussian Parsimonious\n| .  . | ___ | |__ | /  \\/ |_   _ ___| |_        Clustering Models\n| |\\/| |/ _ \\|  __|| |   | | | | / __| __|         with Covariates\n| |  | | (_) | |___| \\__/\\ | |_| \\__ \\ |_ \n\\_|  |_/\\___/\\____/ \\____/_|\\__,_|___/\\__|           version", version, "\n"))
  } else   {
    packageStartupMessage("\nPackage 'MoEClust' version ", version, ".")
  }
  packageStartupMessage(paste("Type '?MoEClust' to see a brief guide to how to use this R package.\nType", sQuote(paste0("citation(", dQuote("MoEClust"),")")) ,"for citing the package in publications.\nType 'MoE_news()' to see new features recent changes and bug fixes.\n"))
}
