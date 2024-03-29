% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Plotting_Functions.R
\name{plot.MoEClust}
\alias{plot.MoEClust}
\title{Plot MoEClust Results}
\usage{
\method{plot}{MoEClust}(x,
     what = c("gpairs", "gating", "criterion", "loglik", "similarity", "uncertainty"),
     ...)
}
\arguments{
\item{x}{An object of class \code{"MoEClust"} generated by \code{\link{MoE_clust}}, or an object of class \code{"MoECompare"} generated by \code{\link{MoE_compare}}. Models with a noise component are facilitated here too.}

\item{what}{The type of graph requested:
\describe{
\item{\code{gpairs}}{A generalised pairs plot. To further customise this plot, arguments to \code{\link{MoE_gpairs}} can be supplied.}
\item{\code{gating}}{The gating network. To further customise this plot, arguments to \code{\link{MoE_plotGate}} and \code{\link[graphics]{matplot}} can be supplied.}
\item{\code{criterion}}{The model selection criteria. To further customise this plot, arguments to \code{\link{MoE_plotCrit}} and \code{\link[mclust]{plot.mclustBIC}} can be supplied.}
\item{\code{loglik}}{The log-likelihood vs. the iterations of the EM algorithm. To further customise this plot, arguments to \code{\link{MoE_plotLogLik}} and \code{\link{plot}} can be supplied.}
\item{\code{similarity}}{The similarity matrix constructed from \code{x$z} at convergence, in the form of a heatmap. To further customise this plot, arguments to \code{\link{MoE_Similarity}} can be supplied.}
\item{\code{uncertainty}}{The clustering uncertainty for every observation. To further customise this plot, arguments to \code{\link{MoE_Uncertainty}} can be supplied.}
}
By default, all of the above graphs are produced.}

\item{...}{Optional arguments to be passed to \code{\link{MoE_gpairs}}, \code{\link{MoE_plotGate}}, \code{\link{MoE_plotCrit}}, \code{\link{MoE_plotLogLik}}, \code{\link{MoE_Similarity}}, \code{\link{MoE_Uncertainty}}, \code{\link[graphics]{matplot}}, \code{\link[mclust]{plot.mclustBIC}} and \code{\link{plot}}. In particular, the argument \code{legendArgs} to \code{\link[mclust]{plot.mclustBIC}} can be passed to \code{\link{MoE_plotCrit}}.}
}
\value{
The visualisation according to \code{what} of the results of a fitted \code{MoEClust} model.
}
\description{
Plot results for fitted MoE_clust mixture models with gating &/or expert network covariates: generalised pairs plots, model selection criteria, the log-likelihood vs. the EM iterations, and the gating network are all currently visualisable.
}
\details{
For more flexibility in plotting, use \code{\link{MoE_gpairs}}, \code{\link{MoE_plotGate}}, \code{\link{MoE_plotCrit}}, \code{\link{MoE_plotLogLik}}, \code{\link{MoE_Similarity}}, and \code{\link{MoE_Uncertainty}} directly.
}
\note{
Caution is advised producing generalised pairs plots when the dimension of the data is large.

Other types of plots are available by first calling \code{\link{as.Mclust}} on the fitted object, and then calling \code{\link[mclust]{plot.Mclust}} on the results. These can be especially useful for univariate data.
}
\examples{
\donttest{data(ais)
res <- MoE_clust(ais[,3:7], gating= ~ BMI, expert= ~ sex,
                 G=2, modelNames="EVE", network.data=ais)

# Plot the gating network
plot(res, what="gating", x.axis=ais$BMI, xlab="BMI")

# Plot the log-likelihood
plot(res, what="loglik", col="blue")

# Plot the uncertainty profile
plot(res, what="uncertainty", type="profile")

# Produce a generalised pairs plot
plot(res, what="gpairs")

# Produce a heatmap of the similarity matrix
plot(res, what="similarity")

# Modify the gpairs plot by passing arguments to MoE_gpairs()
plot(res, what="gpairs", response.type="density", varwidth=TRUE,
     data.ind=c(5,3,4,1,2), jitter=FALSE, show.counts=FALSE)}
}
\references{
Murphy, K. and Murphy, T. B. (2020). Gaussian parsimonious clustering models with covariates and a noise component. \emph{Advances in Data Analysis and Classification}, 14(2): 293-325. <\doi{10.1007/s11634-019-00373-8}>.
}
\seealso{
\code{\link{MoE_clust}}, \code{\link{MoE_stepwise}}, \code{\link{MoE_gpairs}}, \code{\link{MoE_plotGate}}, \code{\link{MoE_plotCrit}}, \code{\link{MoE_plotLogLik}}, \code{\link{MoE_Similarity}}, \code{\link{MoE_Uncertainty}}, \code{\link{as.Mclust}}, \code{\link[mclust]{plot.Mclust}}
}
\author{
Keefe Murphy - <\email{keefe.murphy@mu.ie}>
}
\keyword{main}
\keyword{plotting}
