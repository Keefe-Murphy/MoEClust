#' MoEClust: Finite Gaussian Mixtures of Experts: Model-Based Clustering with Covariates
#'
#' Fits finite Mixtures of Experts models with mclust family covariance structures using the EM algorithm, ie. allows incorporation of covariates into the mixing proportions and/or Gaussian densities of finite mixture models under the various covariance parameterisations in the mclust family.
#' @section Usage:
#' The most important function in the \pkg{MoEClust} package is: \code{\link{MoE_clust}}, for fitting the model via EM with gating and/or expert network covariates, supplied via formula interfaces.
#'
#' Other functions also exist, e.g. \code{\link{MoE_control}}, \code{\link{MoE_crit}}, \code{\link{MoE_dens}}, \code{\link{MoE_estep}}, \code{\link{MoE_compare}}, and \code{\link{MoE_aitken}}, which are all used within \code{\link{MoE_clust}} but are nonetheless made available for standalone use.
#'
#' @section Details:
#' Package: MoEClust
#'
#' Type: Package
#'
#' Version: 0.1.0
#'
#' Date: 2017-08-09
#'
#' Licence: GPL (>=2)
#'
#' @author
#' Keefe Murphy [aut, cre]
#'
#' Maintainer: Keefe Murphy - \href{keefe.murphy@ucd.ie}{<keefe.murphy@ucd.ie>}
#'
#' @docType package
#' @name MoEClust
NULL

.onAttach       <- function(...) if(interactive()) packageStartupMessage(paste("___  ___      _____ _____ _           _   \n|  \\/  |     |  ___/  __ \\ |         | |  \n| .  . | ___ | |__ | /  \\/ |_   _ ___| |_ \n| |\\/| |/ _ \\|  __|| |   | | | | / __| __|\n| |  | | (_) | |___| \\__/\\ | |_| \\__ \\ |_ \n\\_|  |_/\\___/\\____/ \\____/_|\\__,_|___/\\__|           version 0.1.0\nType '?MoECLUST' to see a brief guide to how to use this R package.\nType", sQuote(paste0("citation(", dQuote("MoEClust"),")")) ,"for citing the package in publications."))
