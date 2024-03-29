% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Functions.R
\name{MoE_control}
\alias{MoE_control}
\title{Set control values for use with MoEClust}
\usage{
MoE_control(init.z = c("hc", "quantile", "kmeans", "mclust", "random", "list"),
            noise.args = list(...),
            asMclust = FALSE,
            equalPro = FALSE,
            exp.init = list(...),
            algo = c("EM", "CEM", "cemEM"),
            criterion = c("bic", "icl", "aic"),
            stopping = c("aitken", "relative"),
            z.list = NULL, 
            nstarts = 1L,
            eps = .Machine$double.eps,
            tol = c(1e-05, sqrt(.Machine$double.eps), 1e-08),
            itmax = c(.Machine$integer.max, .Machine$integer.max, 1000L),
            hc.args = list(...),
            km.args = list(...),
            posidens = TRUE,
            init.crit = c("bic", "icl"),
            warn.it = 0L,
            MaxNWts = 1000L,
            verbose = interactive(),
            ...)
}
\arguments{
\item{init.z}{The method used to initialise the cluster labels for the \emph{non-noise} components. Defaults to \code{"hc"}, i.e. model-based agglomerative hierarchical clustering tree as per \code{\link[mclust]{hc}}, for multivariate data (see \code{hc.args}), or \code{"quantile"}-based clustering as per \code{\link{quant_clust}} for univariate data (unless there are expert network covariates incorporated via \code{exp.init$joint} &/or \code{exp.init$clustMD}, in which case the default is again \code{"hc"}). The \code{"quantile"} option is thus only available for univariate data when expert network covariates are not incorporated via \code{exp.init$joint} &/or \code{exp.init$clustMD}, or when expert network covariates are not supplied.

Other options include \code{"kmeans"} (see \code{km.args}), \code{"random"} initialisation (see \code{nstarts} below), a user-supplied \code{"list"}, and a full run of \code{\link[mclust]{Mclust}} (itself initialised via a model-based agglomerative hierarchical clustering tree, again see \code{hc.args}), although this last option \code{"mclust"} will be coerced to \code{"hc"} if there are no \code{gating} &/or \code{expert} covariates within \code{\link{MoE_clust}} (in order to better reproduce \code{\link[mclust]{Mclust}} output).

When \code{init.z="list"}, \code{exp.init$clustMD} is forced to \code{FALSE}; otherwise, when \code{isTRUE(exp.init$clustMD)} and the \code{\link[clustMD]{clustMD}} library is loaded, the \code{init.z} argument instead governs the method by which a call to \code{\link[clustMD]{clustMD}} is initialised. In this instance, \code{"quantile"} will instead default to \code{"hc"}, and the arguments to \code{hc.args} and \code{km.args} will be ignored (unless all \code{\link[clustMD]{clustMD}} model types fail for a given number of components).

When \code{init.z="mclust"} or \code{\link[clustMD]{clustMD}} is successfully invoked (via \code{exp.init$clustMD}), the argument \code{init.crit} (see below) specifies the model-selection criterion (\code{"bic"} or \code{"icl"}) by which the optimal \code{\link[mclust]{Mclust}} or \code{\link[clustMD]{clustMD}} model type to initialise with is determined, and \code{criterion} remains unaffected.

Finally, when the model includes expert network covariates and \code{isTRUE(exp.init$mahalanobis)}, the argument \code{exp.init$estart} (see below) can be used to modify the behaviour of \code{init.z="random"} when \code{nstarts > 1}, toggling between a full run of the EM algorithm for each random initialisation (i.e. \code{exp.init$estart=FALSE}, the default), or a single run of the EM algorithm starting from the best initial partition obtained among the random starts according to the iterative reallocation initialisation routine (i.e. \code{exp.init$estart=TRUE}).}

\item{noise.args}{A list supplying select named parameters to control inclusion of a noise component in the estimation of the mixture. If either or both of the arguments \code{tau0} &/or \code{noise.init} are supplied, a noise component is added to the the model in the estimation.
\describe{
\item{\code{tau0}}{Prior mixing proportion for the noise component. If supplied, a noise component will be added to the model in the estimation, with \code{tau0} giving the prior probability of belonging to the noise component for \emph{all} observations. Typically supplied as a scalar in the interval (0, 1), e.g. \code{0.1}. Can be supplied as a vector when gating covariates are present and \code{noise.args$noise.gate} is \code{TRUE}. This argument can be supplied instead of or in conjunction with the argument \code{noise.init} below.}
\item{\code{noise.init}}{A logical or numeric vector indicating an initial guess as to which observations are noise in the data. If numeric, the entries should correspond to row indices of the data. If supplied, a noise component will be added to the model in the estimation. This argument can be used in conjunction with \code{tau0} above, or can be replaced by that argument also.}
\item{\code{noise.gate}}{A logical indicating whether gating network covariates influence the mixing proportion for the noise component, if any. Defaults to \code{TRUE}, but leads to greater parsimony if \code{FALSE}. Only relevant in the presence of a noise component; only effects estimation in the presence of gating covariates.}
\item{\code{noise.meth}}{The method used to estimate the volume when a noise component is invoked. Defaults to \code{\link[mclust]{hypvol}}. For univariate data, this argument is ignored and the range of the data is used instead (unless \code{noise.vol} below is specified). The options \code{"convexhull"} and \code{"ellipsoidhull"} require loading the \code{geometry} and \code{cluster} libraries, respectively. This argument is only relevant if \code{noise.vol} below is not supplied.}
\item{\code{noise.vol}}{This argument can be used to override the argument \code{noise.meth} by specifying the (hyper)volume directly, i.e. specifying an improper uniform density. This will override the use of the range of the response data for univariate data if supplied. Note that the (hyper)volume, rather than its inverse, is supplied here. This can affect prediction and the location of the MVN ellipses for \code{\link{MoE_gpairs}} plots (see \code{\link{noise_vol}}).}
\item{\code{equalNoise}}{Logical which is \strong{only} invoked when \code{isTRUE(equalPro)} and gating covariates are not supplied. Under the default setting (\code{FALSE}), the mixing proportion for the noise component is estimated, and remaining mixing proportions are equal; when \code{TRUE} all components, including the noise component, have equal mixing proportions.}
\item{\code{discard.noise}}{A logical governing how the means are summarised in \code{parameters$mean} and by extension the location of the MVN ellipses in \code{\link{MoE_gpairs}} plots for models with \emph{both} expert network covariates and a noise component (otherwise this argument is irrelevant). 

The means for models with expert network covariates are summarised by the posterior mean of the fitted values. By default (\code{FALSE}), the mean of the noise component is accounted for in the posterior mean. Otherwise, or when the mean of the noise component is unavailable (due to having been manually supplied via \code{noise.args$noise.vol}), the \code{z} matrix is renormalised after discarding the column corresponding to the noise component prior to computation of the posterior mean. The renormalisation approach can be forced by specifying \code{noise.args$discard.noise=TRUE}, even when the mean of the noise component is available. For models with a noise component fitted with \code{algo="CEM"}, a small extra E-step is conducted for observations assigned to the non-noise components in this case.}
}
In particular, the argument \code{noise.meth} will be ignored for high-dimensional \code{n <= d} data, in which case the argument \code{noise.vol} \emph{must be} specified. Note that this forces \code{noise.args$discard.noise} to \code{TRUE}. See \code{\link{noise_vol}} for more details.

The arguments \code{tau0} and \code{noise.init} can be used separately, to provide alternative means to invoke a noise component. However, they can also be supplied together, in which case observations corresponding to \code{noise.init} have probability \code{tau0} (rather than 1) of belonging to the noise component.}

\item{asMclust}{The default values of \code{stopping} and \code{hc.args$hcUse} (see below) are such that results for models with \emph{no covariates in either network} are liable to differ from results for equivalent models obtained via \code{\link[mclust]{Mclust}}. \pkg{MoEClust} uses \code{stopping="aitken"} and \code{hcUse="VARS"} by default, while \pkg{mclust} always implicitly uses \code{stopping="relative"} and defaults to \code{hcUse="SVD"}.

\code{asMclust} is a logical variable (\code{FALSE}, by default) which functions as a simple convenience tool for overriding these two arguments (even if explicitly supplied!) such that they behave like the function \code{\link[mclust]{Mclust}}. Other \emph{user-specified} arguments which differ from \pkg{mclust} are not affected by \code{asMclust}, as their defaults already correspond to \pkg{mclust}. Results may still differ slightly as \pkg{MoEClust} calculates log-likelihood values with greater precision. Finally, note that \code{asMclust=TRUE} is invoked even for models with covariates which are not accommodated by \pkg{mclust}.}

\item{equalPro}{Logical variable indicating whether or not the mixing proportions are to be constrained to be equal in the model. Default: \code{equalPro = FALSE}. Only relevant when \code{gating} covariates are \emph{not} supplied within \code{\link{MoE_clust}}, otherwise ignored. In the presence of a noise component (see \code{noise.args}), only the mixing proportions for the non-noise components are constrained to be equal (by default, see \code{equalNoise}), after accounting for the noise component.}

\item{exp.init}{A list supplying select named parameters to control the initialisation routine in the presence of \emph{expert} network covariates (otherwise ignored):
\describe{
\item{\code{joint}}{A logical indicating whether the initial partition is obtained on the joint distribution of the response and expert network covariates (defaults to \code{TRUE}) or just the response variables (\code{FALSE}). By default, only continuous expert network covariates are considered (see \code{exp.init$clustMD} below). Only relevant when \code{init.z} is not \code{"random"} (unless \code{isTRUE(exp.init$clustMD)}, in which case \code{init.z} specifies the initialisation routine for a call to \code{\link[clustMD]{clustMD}}). This will render the \code{"quantile"} option to \code{init.z} for univariate data unusable if continuous expert network covariates are supplied &/or categorical/ordinal expert network covariates are supplied when \code{isTRUE(exp.init$clustMD)} and the \code{\link[clustMD]{clustMD}} library is loaded.}
\item{\code{mahalanobis}}{A logical indicating whether to iteratively reallocate observations during the initialisation phase to the component corresponding to the expert network regression to which it's closest to the fitted values of in terms of Mahalanobis distance (defaults to \code{TRUE}). This will ensure that each component can be well modelled by a single expert prior to running the EM/CEM algorithm.}
\item{\code{estart}}{A logical governing the behaviour of \code{init.z="random"} when \code{nstarts > 1} in the presence of expert network covariates. Only relevant when \code{isTRUE(exp.init$mahalanobis)}. Defaults to \code{FALSE}; i.e. all random starts are put through full runs of the EM algorithm. When \code{TRUE}, all random starts are put through the initial iterative reallocation routine prior to a full run of EM for only the single best random initial partition obtained. See the last set of \strong{Examples} below.}
\item{\code{\link[clustMD]{clustMD}}}{A logical indicating whether categorical/ordinal covariates should be incorporated when using the joint distribution of the response and expert network covariates for initialisation (defaults to \code{FALSE}). Only relevant when \code{isTRUE(exp.init$joint)}. Requires the use of the \code{\link[clustMD]{clustMD}} library. Note that initialising in this manner involves fitting all \code{\link[clustMD]{clustMD}} model types in parallel for all numbers of components considered, and may fail (especially) in the presence of nominal expert network covariates.

Unless \code{init.z="list"}, supplying this argument as \code{TRUE} when the \code{\link[clustMD]{clustMD}} library is loaded has the effect of superseding the \code{init.z} argument: this argument now governs instead how the call to \code{\link[clustMD]{clustMD}} is initialised (unless all \code{\link[clustMD]{clustMD}} model types fail for a given number of components, in which case \code{init.z} is invoked \emph{instead} to initialise for \code{G} values for which all \code{\link[clustMD]{clustMD}} model types failed). Similarly, the arguments \code{hc.args} and \code{km.args} will be ignored (again, unless all \code{\link[clustMD]{clustMD}} model types fail for a given number of components).}
\item{\code{max.init}}{The maximum number of iterations for the Mahalanobis distance-based reallocation procedure when \code{exp.init$mahalanobis} is \code{TRUE}. Defaults to \code{.Machine$integer.max}.}
\item{\code{identity}}{A logical indicating whether the identity matrix (corresponding to the use of the Euclidean distance) is used in place of the covariance matrix of the residuals (corresponding to the use of the Mahalanobis distance). Defaults to \code{FALSE} for multivariate response data but defaults to \code{TRUE} for univariate response data. Setting \code{identity=TRUE} with multivariate data may be advisable when the dimensions of the data are such that the covariance matrix cannot be inverted (otherwise, the pseudo-inverse is used under the \code{FALSE} default).}
\item{\code{drop.break}}{When \code{isTRUE(exp.init$mahalanobis)} observations will be completely in or out of a component during the initialisation phase. As such, it may occur that constant columns will be present when building a given component's expert regression (particularly for categorical covariates). It may also occur, due to this partitioning, that "unseen" data, when calculating the residuals, will have new factor levels. When \code{isTRUE(exp.init$drop.break)}, the Mahalanobis distance based initialisation phase will explicitly fail in either of these scenarios.

Otherwise, \code{\link{drop_constants}} and \code{\link{drop_levels}} will be invoked when \code{exp.init$drop.break} is \code{FALSE} (the default) to \emph{try} to remedy the situation. In any case, only a warning that the initialisation step failed will be printed, regardless of the value of \code{exp.init$drop.break}.}
}}

\item{algo}{Switch controlling whether models are fit using the \code{"EM"} (the default) or \code{"CEM"} algorithm. The option \code{"cemEM"} allows running the EM algorithm starting from convergence of the CEM algorithm.}

\item{criterion}{When either \code{G} or \code{modelNames} is a vector, \code{criterion} determines whether the \code{"bic"} (Bayesian Information Criterion), \code{"icl"} (Integrated Complete Likelihood), \code{"aic"} (Akaike Information Criterion) is used to determine the 'best' model when gathering output. Note that all criteria will be returned in any case.}

\item{stopping}{The criterion used to assess convergence of the EM/CEM algorithm. The default (\code{"aitken"}) uses Aitken's acceleration method via \code{\link{aitken}}, otherwise the \code{"relative"} change in log-likelihood is monitored (which may be less strict). The \code{"relative"} option corresponds to the stopping criterion used by \code{\link[mclust]{Mclust}}: see \code{asMclust} above. 

Both stopping rules are ultimately governed by \code{tol[1]}. When the \code{"aitken"} method is employed, the asymptotic estimate of the final converged maximised log-likelihood is also returned as \code{linf} for models with 2 or more components, though the largest element of the returned vector \code{loglik} still gives the log-likelihood value achieved by the parameters returned at convergence, under both \code{stopping} methods (see \code{\link{MoE_clust}}).}

\item{z.list}{A user supplied list of initial cluster allocation matrices, with number of rows given by the number of observations, and numbers of columns given by the range of component numbers being considered. In particular, \code{z.list} must only include columns corresponding to \emph{non-noise} components when using this method to initialise in the presence of a noise component. Only relevant if \code{init.z == "z.list"}. These matrices are allowed correspond to both soft or hard clusterings, and will be internally normalised so that the rows sum to 1. See \code{noise.init} and \code{tau0} above for details on initialisation in the presence of a noise component.}

\item{nstarts}{The number of random initialisations to use when \code{init.z="random"}. Defaults to \code{1}. When there are no expert covariates (or when \code{exp.init$mahalanobis=FALSE} or \code{exp.init$estart=FALSE}), the results will be based on the random start yielding the highest estimated log-likelihood after each initial partition is subjected to a full run of the EM algorithm. Note, in this case, that all \code{nstarts} random initialisations are affected by \code{exp.init$mahalanobis}, if invoked in the presence of expert network covariates, which may remove some of the randomness. 

Conversely, if \code{exp.init$mahalanobis=TRUE} and \code{exp.init$estart=TRUE}, all \code{nstarts} random starts are put through the initial iterative reallocation routine and only the single best initial partition uncovered is put through the full run of the EM algorithm. See \code{init.z} and \code{exp.init$estart} above for more details, though note that \code{exp.init$mahalanobis=TRUE} and \code{exp.init$estart=FALSE}, by default.}

\item{eps}{A scalar tolerance associated with deciding when to terminate computations due to computational singularity in covariances. Smaller values of \code{eps} allow computations to proceed nearer to singularity. The default is the relative machine precision \code{.Machine$double.eps}, which is approximately \emph{2e-16} on IEEE-compliant machines.}

\item{tol}{A vector of length three giving \emph{relative} convergence tolerances for 1) the log-likelihood of the EM/CEM algorithm, 2) parameter convergence in the inner loop for models with iterative M-step (\code{"VEI", "VEE", "EVE", "VVE", "VEV"}), and 3) optimisation in the multinomial logistic regression in the gating network, respectively. The default is \code{c(1e-05, sqrt(.Machine$double.eps), 1e-08)}. If only one number is supplied, it is used as the tolerance for all three cases given.}

\item{itmax}{A vector of length three giving integer limits on the number of iterations for 1) the EM/CEM algorithm, 2) the inner loop for models with iterative M-step (\code{"VEI", "VEE", "EVE", "VVE", "VEV"}), and 3) the multinomial logistic regression in the gating network, respectively.

The default is \code{c(.Machine$integer.max, .Machine$integer.max, 1000L)}, allowing termination to be completely governed by \code{tol[1]} & \code{tol[2]} for the inner and outer loops of the EM/CEM algorithm. If only one number is supplied, it is used as the iteration limit for the outer loop only and the other elements of \code{itmax} retain their usual defaults.

If, for any model with gating covariates, the multinomial logistic regression in the gating network fails to converge in \code{itmax[3]} iterations at any stage of the EM/CEM algorithm, an appropriate warning will be printed, prompting the user to modify this argument.}

\item{hc.args}{A list supplying select named parameters to control the initialisation of the cluster allocations when \code{init.z="hc"} (or when \code{init.z="mclust"}, which itself relies on \code{\link[mclust]{hc}}), unless \code{isTRUE(exp.init$clustMD)}, the \code{\link[clustMD]{clustMD}} library is loaded, and none of the \code{\link[clustMD]{clustMD}} model types fail (otherwise irrelevant):
\describe{
\item{\code{hcUse}}{A string specifying the type of input variables to be used. This defaults to \code{"VARS"} here, unlike \pkg{mclust} which defaults to \code{"SVD"}. Other allowable values are documented in \code{\link[mclust]{mclust.options}}. See \code{asMclust} above.}
\item{\code{hc.meth}}{A character string indicating the model to be used when hierarchical clustering (see \code{\link[mclust]{hc}}) is employed for initialisation (either when \code{init.z="hc"} or \code{init.z="mclust"}). Defaults to \code{"EII"} for high-dimensional data, or \code{"VVV"} otherwise.}
}}

\item{km.args}{A list supplying select named parameters to control the initialisation of the cluster allocations when \code{init.z="kmeans"}, unless \code{isTRUE(exp.init$clustMD)}, the \code{\link[clustMD]{clustMD}} library is loaded, and none of the \code{\link[clustMD]{clustMD}} model types fail (otherwise irrelevant):
\describe{
\item{\code{kstarts}}{The number of random initialisations to use. Defaults to 10.}
\item{\code{kiters}}{The maximum number of K-Means iterations allowed. Defaults to 10.}
}}

\item{posidens}{A logical governing whether to continue running the algorithm even in the presence of positive log-densities. Defaults to \code{TRUE}, but setting \code{posidens=FALSE} can help to safeguard against spurious solutions, which will be instantly terminated if positive log-densities are encountered. Note that versions of this package prior to and including version \code{1.3.1} always implicitly assumed \code{posidens=FALSE}.}

\item{init.crit}{The criterion to be used to determine the optimal model type to initialise with, when \code{init.z="mclust"} or when \code{isTRUE(exp.init$clustMD)} and the \code{\link[clustMD]{clustMD}} library is loaded (one of \code{"bic"} or \code{"icl"}). Defaults to \code{"icl"} when \code{criterion="icl"}, otherwise defaults to \code{"bic"}. The \code{criterion} argument remains unaffected.}

\item{warn.it}{A single number giving the iteration count at which a warning will be printed if the EM/CEM algorithm has failed to converge. Defaults to \code{0}, i.e. no warning (which is true for any \code{warn.it} value less than \code{3}), otherwise the message is printed regardless of the value of \code{verbose}. If non-zero, \code{warn.it} should be moderately large, but obviously less than \code{itmax[1]}. A warning will always be printed if one of more models fail to converge in \code{itmax[1]} iterations.}

\item{MaxNWts}{The maximum allowable number of weights in the call to \code{\link[nnet]{multinom}} for the multinomial logistic regression in the gating network. There is no intrinsic limit in the code, but increasing \code{MaxNWts} will probably allow fits that are very slow and time-consuming. It may be necessary to increase \code{MaxNWts} when categorical concomitant variables with many levels are included or the number of components is high.}

\item{verbose}{Logical indicating whether to print messages pertaining to progress to the screen during fitting. By default is \code{TRUE} if the session is interactive, and \code{FALSE} otherwise. If \code{FALSE}, warnings and error messages will still be printed to the screen, but everything else will be suppressed.}

\item{...}{Catches unused arguments.}
}
\value{
A named list in which the names are the names of the arguments and the values are the values supplied to the arguments.
}
\description{
Supplies a list of arguments (with defaults) for use with \code{\link{MoE_clust}}.
}
\details{
\code{\link{MoE_control}} is provided for assigning values and defaults within \code{\link{MoE_clust}} and \code{\link{MoE_stepwise}}.

While the \code{criterion} argument controls the choice of the optimal number of components and GPCM/\pkg{mclust} model type, \code{\link{MoE_compare}} is provided for choosing between fits with different combinations of covariates or different initialisation settings.
}
\note{
Note that successfully invoking \code{exp.init$clustMD} (though it defaults to \code{FALSE}) affects the role of the arguments \code{init.z}, \code{hc.args}, and \code{km.args}. Please read the documentation above carefully in this instance.

The initial allocation matrices before and after the invocation of the \code{exp.init} related arguments are both stored as attributes in the object returned by \code{\link{MoE_clust}} (named \code{"Z.init"} and \code{"Exp.init"}, respectively). If \code{init.z="random"} and \code{nstarts > 1}, the allocations corresponding to the best random start are stored (regardless of whether \code{exp.init$estart} is invoked or not). This can be useful for supplying \code{z.list} for future fits.
}
\examples{
ctrl1 <- MoE_control(criterion="icl", itmax=100, warn.it=15, init.z="random", nstarts=5)

data(CO2data)
GNP   <- CO2data$GNP
\donttest{res   <- MoE_clust(CO2data$CO2, G=2, expert = ~ GNP, control=ctrl1)

# Alternatively, specify control arguments directly
res2  <- MoE_clust(CO2data$CO2, G=2, expert = ~ GNP, stopping="relative")}

# Supplying ctrl1 without naming it as 'control' can throw an error
\dontrun{
res3  <- MoE_clust(CO2data$CO2, G=2, expert = ~ GNP, ctrl1)}

# Similarly, supplying control arguments via a mix of the ... construct
# and the named argument 'control' also throws an error
\dontrun{
res4  <- MoE_clust(CO2data$CO2, G=2, expert = ~ GNP, control=ctrl1, init.z="kmeans")}

\donttest{# Initialise via the mixed-type joint distribution of response & covariates
# Let the ICL criterion determine the optimal clustMD model type
# Constrain the mixing proportions to be equal
ctrl2 <- MoE_control(exp.init=list(clustMD=TRUE), init.crit="icl", equalPro=TRUE)
data(ais)
library(clustMD)
res4  <- MoE_clust(ais[,3:7], G=2, modelNames="EVE", expert= ~ sex,
                   network.data=ais, control=ctrl2)

# Include a noise component by specifying its prior mixing proportion
res5  <- MoE_clust(ais[,3:7], G=2, modelNames="EVE", expert= ~ sex,
                   network.data=ais, tau0=0.1)
                   
# Investigate the use of random starts
sex  <- ais$sex
# resA uses deterministic starting values (by default) for each G value
 system.time(resA <- MoE_clust(ais[,3:7], G=2, expert=~sex, equalPro=TRUE))
# resB passes each random start through the entire EM algorithm for each G value
 system.time(resB <- MoE_clust(ais[,3:7], G=2, expert=~sex, equalPro=TRUE,
                              init.z="random", nstarts=10))
# resC passes only the "best" random start through the EM algorithm for each G value
 system.time(resC <- MoE_clust(ais[,3:7], G=2, expert=~sex, equalPro=TRUE,
                               init.z="random", nstarts=10, estart=TRUE))
# Here, all three settings (listed here in order of speed) converge to the same model
 MoE_compare(resA, resC, resB)}
}
\seealso{
\code{\link{MoE_clust}}, \code{\link{MoE_stepwise}}, \code{\link{aitken}}, \code{\link[mclust]{Mclust}}, \code{\link[mclust]{hc}}, \code{\link[mclust]{mclust.options}}, \code{\link{quant_clust}}, \code{\link[clustMD]{clustMD}}, \code{\link{noise_vol}}, \code{\link[mclust]{hypvol}}, \code{\link[geometry]{convhulln}}, \code{\link[cluster]{ellipsoidhull}}, \code{\link{MoE_compare}}, \code{\link[nnet]{multinom}}
}
\author{
Keefe Murphy - <\email{keefe.murphy@mu.ie}>
}
\keyword{control}
