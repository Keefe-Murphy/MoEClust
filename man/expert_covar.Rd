% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Functions.R
\name{expert_covar}
\alias{expert_covar}
\title{Account for extra variability in covariance matrices with expert covariates}
\usage{
expert_covar(x,
             weighted = TRUE,
             ...)
}
\arguments{
\item{x}{An object of class \code{"MoEClust"} generated by \code{\link{MoE_clust}}, or an object of class \code{"MoECompare"} generated by \code{\link{MoE_compare}}. Models with a noise component are facilitated here too.}

\item{weighted}{A logical indicating whether the estimated cluster membership probabilities should be used to provide a weighted estimate of the variability due to the component means. Defaults to \code{TRUE}. The option \code{weighted=FALSE} is provided only so that previous behaviour under earlier versions of \pkg{MoEClust} can be recovered but is otherwise not recommended.}

\item{...}{Catches unused arguments.}
}
\value{
The \code{variance} component only from the \code{parameters} list from the output of a call to \code{\link{MoE_clust}}, modified accordingly.
}
\description{
In the presence of expert network covariates, this helper function modifies the component-specific covariance matrices of a \code{"MoEClust"} object, in order to account for the extra variability due to the component means, usually resulting in bigger shapes & sizes for the MVN ellipses in \code{\link{MoE_gpairs}} plots. The function also works for univariate response data.
}
\details{
This function is used internally by \code{\link{MoE_gpairs}}, \code{\link{plot.MoEClust}(x, what="gpairs")}, and \code{\link[=as.Mclust.MoEClust]{as.Mclust}}, for visualisation purposes.
}
\note{
The \code{modelName} of the resulting \code{variance} object may not correspond to the model name of the \code{"MoEClust"} object, in particular \code{scale}, \code{shape}, &/or \code{orientation} may no longer be constrained across clusters, and \code{cholsigma}, if it was in the input, will be discarded from the output. Usually, the \code{modelName} of the transformed \code{variance} object will be \code{"VVV"} for multivariate data and \code{"V"} for univariate data, but not always. Furthermore, the output will drop certain row and column names from the result.
}
\examples{
data(ais)
res   <- MoE_clust(ais[,3:7], G=2, gating= ~ 1, expert= ~ sex,
                   network.data=ais, modelNames="EEE", equalPro=TRUE)

# Extract the variance object
res$parameters$variance

# Modify the variance object
expert_covar(res)
}
\references{
Murphy, K. and Murphy, T. B. (2020). Gaussian parsimonious clustering models with covariates and a noise component. \emph{Advances in Data Analysis and Classification}, 14(2): 293-325. <\doi{10.1007/s11634-019-00373-8}>.
}
\seealso{
\code{\link{MoE_clust}}, \code{\link{MoE_gpairs}}, \code{\link{plot.MoEClust}}, \code{\link[=as.Mclust.MoEClust]{as.Mclust}}
}
\author{
Keefe Murphy - <\email{keefe.murphy@mu.ie}>
}
\keyword{utility}
