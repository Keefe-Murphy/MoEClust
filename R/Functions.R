#' MoEClust: Gaussian Parsimonious Clustering Models with Covariates and a Noise Component
#' 
#' Fits MoEClust models: Gaussian Mixture of Experts models with GPCM/\pkg{mclust}-family covariance structures. In other words, performs model-based clustering via the EM/CEM algorithm where covariates are allowed to enter neither, either, or both the mixing proportions (gating network) and/or component densities (expert network) of a Gaussian Parsimonious Clustering Model, with or without an additional noise component. Additional arguments are available via the function \code{\link{MoE_control}}, including the specification of a noise component, controls on the initialisation of the algorithm, and more. 
#' @param data A numeric vector, matrix, or data frame of observations. Categorical variables are not allowed. If a matrix or data frame, rows correspond to observations and columns correspond to variables.
#' @param G An integer vector specifying the numbers of mixture components (clusters) to fit. Defaults to \code{G=1:9}. Must be a strictly positive integer, unless a noise component is included in the estimation, in which case \code{G=0} is allowed and included by default. (see \code{\link{MoE_control}}).
#' @param modelNames A vector of character strings indicating the models to be fitted in the EM/CEM phase of clustering. With \code{n} observations and \code{d} variables, the defaults are:
#' \tabular{ll}{
#' for univariate data \tab \code{c("E", "V")}\cr
#' for multivariate data \eqn{n > d}{n > d} \tab \code{mclust.options("emModelNames")}\cr
#' for high-dimensional multivariate data \eqn{n \leq d}{n <= d} \tab \code{c("EII", "VII", "EEI", "EVI", "VEI", "VVI")}
#' }
#'
#' For single-component models these options reduce to:
#' \tabular{ll}{
#' for univariate data \tab \code{"E"}\cr
#' for multivariate data \eqn{n > d}{n > d} \tab \code{c("EII", "EEI", "EEE")}\cr
#' for high-dimensional multivariate data \eqn{n \leq d}{n <= d}  \tab \code{c("EII", "EEI")}
#' }
#' For zero-component models with a noise component only the \code{"E"} and \code{"EII"} models will be fit for univariate and multivariate data, respectively. The help file for \code{\link[mclust]{mclustModelNames}} further describes the available models (though the \code{"X"} in the single-component models will be coerced to \code{"E"} if supplied that way). For single-component models, other model names equivalent to those above can be supplied, but will be coerced to those above.
#' @param gating A \code{\link[stats]{formula}} for determining the model matrix for the multinomial logistic regression in the gating network when fixed covariates enter the mixing proportions. Defaults to \code{~1}, i.e. no covariates. This will be ignored where \code{G=1}. Continuous, categorical, and/or ordinal covariates are allowed. Logical covariates will be coerced to factors. Interactions, transformations, and higher order terms are permitted: the latter \strong{must} be specified explicitly using the \code{AsIs} operator (\code{\link{I}}). The specification of the LHS of the formula is ignored. Intercept terms are included by default.
#' @param expert A \code{\link[stats]{formula}} for determining the model matrix for the (multivariate) WLS in the expert network when fixed covariates are included in the component densities. Defaults to \code{~1}, i.e. no covariates. Continuous, categorical, and/or ordinal covariates are allowed. Logical covariates will be coerced to factors. Interactions, transformations, and higher order terms are permitted: the latter \strong{must} be specified explicitly using the \code{AsIs} operator (\code{\link{I}}). The specification of the LHS of the formula is ignored. Intercept terms are included by default.
#' @param control A list of control parameters for the EM/CEM and other aspects of the algorithm. The defaults are set by a call to \code{\link{MoE_control}}. In particular, arguments pertaining to the inclusion of an additional noise component are documented here.
#' @param network.data An optional data frame in which to look for the covariates in the \code{gating} &/or \code{expert} network formulas, if any. If not found in \code{network.data}, any supplied \code{gating} &/or \code{expert} covariates are taken from the environment from which \code{MoE_clust} is called. Try to ensure the names of variables in \code{network.data} do not match any of those in \code{data}.
#' @param ... An alternative means of passing control parameters directly via the named arguments of \code{\link{MoE_control}}. Do not pass the output from a call to \code{\link{MoE_control}} here! This argument is only relevant for the \code{\link{MoE_clust}} function and will be ignored for the associated \code{print} and \code{summary} functions.
#' @param x,object,digits Arguments required for the \code{print} and \code{summary} functions: \code{x} and \code{object} are objects of class \code{"MoEClust"} resulting from a call to \code{\link{MoE_clust}}, while \code{digits} gives the number of decimal places to round to for printing purposes (defaults to 3).

#' @importFrom matrixStats "colMeans2" "colSums2" "rowLogSumExps" "rowMaxs" "rowMins" "rowSums2"
#' @importFrom mclust "emControl" "hc" "hclass" "hcE" "hcEEE" "hcEII" "hcV" "hcVII" "hcVVV" "Mclust" "mclust.options" "mclustBIC" "mclustICL" "mclustModelNames" "mclustVariance" "mstep" "mstepE" "mstepEEE" "mstepEEI" "mstepEEV" "mstepEII" "mstepEVE" "mstepEVI" "mstepEVV" "mstepV" "mstepVEE" "mstepVEI" "mstepVEV" "mstepVII" "mstepVVE" "mstepVVI" "mstepVVV" "nVarParams" "unmap"
#' @importFrom mvnfast "dmvn"
#' @importFrom nnet "multinom"
#' @return A list (of class \code{"MoEClust"}) with the following named entries, mostly corresponding to the chosen optimal model (as determined by the \code{criterion} within \code{\link{MoE_control}}):
#' \item{\code{call}}{The matched call.}
#' \item{\code{data}}{The input data, as a \code{data.frame}.}
#' \item{\code{modelName}}{A character string denoting the GPCM/\pkg{mclust} model type at which the optimal \code{criterion} occurs.}
#' \item{\code{n}}{The number of observations in the \code{data}.}
#' \item{\code{d}}{The dimension of the \code{data}.}
#' \item{\code{G}}{The optimal number of mixture components.}
#' \item{\code{BIC}}{A matrix of \emph{all} BIC values with \code{length{G}} rows and \code{length(modelNames)} columns. May include missing entries: \code{NA} represents models which were not visited, \code{-Inf} represents models which were terminated due to error, for which a log-likelihood could not be estimated. Inherits the classes \code{"MoECriterion"} and \code{"mclustBIC"}, for which dedicated printing and plotting functions exist, respectively.}
#' \item{\code{ICL}}{A matrix of \emph{all} ICL values with \code{length{G}} rows and \code{length(modelNames)} columns. May include missing entries: \code{NA} represents models which were not visited, \code{-Inf} represents models which were terminated due to error, for which a log-likelihood could not be estimated. Inherits the classes \code{"MoECriterion"} and \code{"mclustICL"}, for which dedicated printing and plotting functions exist, respectively.}
#' \item{\code{AIC}}{A matrix of \emph{all} AIC values with \code{length{G}} rows and \code{length(modelNames)} columns. May include missing entries: \code{NA} represents models which were not visited, \code{-Inf} represents models which were terminated due to error, for which a log-likelihood could not be estimated. Inherits the classes \code{"MoECriterion"} and \code{"mclustAIC"}, for which dedicated printing and plotting functions exist, respectively.}
#' \item{\code{bic}}{The BIC value corresponding to the optimal model. May not necessarily be the optimal BIC.}
#' \item{\code{icl}}{The ICL value corresponding to the optimal model. May not necessarily be the optimal ICL.}
#' \item{\code{aic}}{The AIC value corresponding to the optimal model. May not necessarily be the optimal AIC.}
#' \item{\code{gating}}{An object of class \code{"MoE_gating"} and either \code{"multinom"} or \code{"glm"} (for single-component models) giving the \code{\link[nnet]{multinom}} regression coefficients of the \code{gating} network. If \code{gating} covariates were \emph{NOT} supplied (or the best model has just one component), this corresponds to a RHS of ~1, otherwise the supplied \code{gating} formula. As such, a fitted \code{gating} network is always returned even in the absence of supplied covariates. The number of parameters to penalise by for \code{\link{MoE_crit}} is given by \code{length(coef(gating))}, and the \code{gating} formula used is stored here as an attribute. If there is a noise component (and the option \code{noise.gate=TRUE} is invoked), its coefficients are those for the \emph{last} component. \strong{Users are cautioned against making inferences about statistical significance from summaries of the coefficients in the gating network}.}
#' \item{\code{expert}}{An object of class \code{"MoE_expert"} and \code{"lm"} giving the (multivariate) WLS regression coefficients of the \code{expert} network. If \code{expert} covariates were NOT supplied, this corresponds to a RHS of ~1, otherwise the supplied \code{expert} formula. As such, a fitted \code{expert} network is always returned even in the absence of supplied covariates. The number of parameters to penalise by for \code{\link{MoE_crit}} is given by \code{G * length(coef(expert[[1]]))}, and the \code{expert} formula used is stored here is an attribute. \strong{Users are cautioned against making inferences about statistical significance from summaries of the coefficients in the expert network}.}
#' \item{\code{LOGLIK}}{A matrix of \emph{all} maximal log-likelihood values with \code{length{G}} rows and \code{length(modelNames)} columns. May include missing entries: \code{NA} represents models which were not visited, \code{-Inf} represents models which were terminated due to error, for which a log-likelihood could not be estimated. Inherits the classes \code{"MoECriterion"} and \code{"mclustLoglik"}, for which dedicated printing and plotting functions exist, respectively.}
#' \item{\code{loglik}}{The vector of increasing log-likelihood values for every EM/CEM iteration under the optimal model. The last element of this vector is the maximum log-likelihood achieved by the parameters returned at convergence.}
#' \item{\code{linf}}{An asymptotic estimate of the final converged maximised log-likelihood. Returned when \code{stopping="aitken"} and \code{G > 1} (see \code{\link{MoE_control}} and \code{\link{aitken}}), otherwise the last element of \code{loglik} is returned instead.}
#' \item{\code{df}}{The number of estimated parameters in the optimal model (i.e. the number of 'used' degrees of freedom). Subtract this number from \code{n} to get the degrees of freedom. The number of parameters due to the gating network, expert network, and covariance matrices are also stored here as attributes of \code{df}.}
#' \item{\code{iters}}{The total number of EM/CEM iterations for the optimal model.}
#' \item{\code{hypvol}}{The hypervolume parameter for the noise component if required, otherwise set to \code{NA} (see \code{\link{MoE_control}}).}
#' \item{\code{parameters}}{A list with the following named components:
#' \describe{
#' \item{\code{pro}}{The mixing proportions: either a vector of length \code{G} or, if \code{gating} covariates were supplied, a matrix with an entry for each observation (rows) and component (columns).}
#' \item{\code{mean}}{The means of each component. If there is more than one component, this is a matrix whose \emph{k}-th column is the mean of the \emph{k}-th component of the mixture model.
#'
#' For models with expert network covariates, this is given by the posterior mean of the fitted values, otherwise the posterior mean of the response is reported. For models with expert network covariates, the \emph{observation-specific} means can be accessed by calling \code{predict} on each object in the list given by \code{expert}.}
#' \item{\code{variance}}{A list of variance parameters of each component of the model. The components of this list depend on the model type specification. See the help file for \code{\link[mclust]{mclustVariance}} for details. Also see \code{\link{expert_covar}} for an alternative approach to summarising the variance parameters in the presence of expert network covariates.}
#' \item{\code{Vinv}}{The inverse of the hypervolume parameter for the noise component if required, otherwise set to \code{NULL} (see \code{\link{MoE_control}}).}
#' }}
#' \item{\code{z}}{The final responsibility matrix whose \code{[i,k]}-th entry is the probability that observation \emph{i} belonds to the \emph{k}-th component. If there is a noise component, its values are found in the \emph{last} column.}
#' \item{\code{classification}}{The vector of cluster labels for the chosen model corresponding to \code{z}, i.e. \code{max.col(z)}. Observations belonging to the noise component will belong to component \code{0}.}
#' \item{\code{uncertainty}}{The uncertainty associated with the \code{classification}.}
#' \item{\code{net.covs}}{A data frame gathering the unique set of covariates used in the \code{gating} and \code{expert} networks, if any. Will contain zero columns in the absence of gating or expert network covariates. Supplied gating covariates will be exluded if the optimal model has only one component. May have fewer columns than covariates supplied via the \code{network.data} argument also, as only the included covariates are gathered here.}
#' \item{\code{resid.data}}{In the presence of expert network covariates, this is the augmented data actually used in the clustering at convergence, as a list of \code{G} matrices of WLS residuals of dimension \code{n * d}. Will contain zero columns in the absence of expert network covariates.}
#' \item{\code{DF}}{A matrix giving the numbers of estimated parameters (i.e. the number of 'used' degrees of freedom) for \emph{all} visited models, with \code{length{G}} rows and \code{length(modelNames)} columns. Subtract these numbers from \code{n} to get the degrees of freedom. May include missing entries: \code{NA} represents models which were not visited, \code{-Inf} represents models which were terminated due to error, for which parameters could not be estimated. Inherits the classes \code{"MoECriterion"} and \code{"mclustDF"}, for which dedicated printing and plotting functions exist, respectively.}
#' \item{\code{ITERS}}{A matrix giving the total number of EM/CEM iterations for \emph{all} visited models, with \code{length{G}} rows and \code{length(modelNames)} columns. May include missing entries: \code{NA} represents models which were not visited, \code{Inf} represents models which were terminated due to singularity/error and thus would never have converged. Inherits the classes \code{"MoECriterion"} and \code{"mclustITERS"}, for which dedicated printing and plotting functions exist, respectively.}
#' Dedicated \code{\link[=plot.MoEClust]{plot}}, \code{\link[=predict.MoEClust]{predict}}, \code{print} and \code{summary} functions exist for objects of class \code{"MoEClust"}. The results can be coerced to the \code{"Mclust"} class to access other functions from the \pkg{mclust} package via \code{\link{as.Mclust}}.
#' @details The function effectively allows 6 different types of Gaussian Mixture of Experts model (as well as the different models in the GPCM/\pkg{mclust} family, for each): i) the standard finite Gaussian mixture with no covariates, ii) fixed covariates only in the gating network, iii) fixed covariates only in the expert network, iv) the full Mixture of Experts model with fixed covariates entering both the mixing proportions and component densities. By constraining the mixing proportions to be equal (see \code{equalPro} in \code{\link{MoE_control}}) two extra special cases are facilitated when gating covariates are excluded. 
#' 
#' Note that having the same covariates in both networks is allowed. So too are interactions, transformations, and higher order terms (see \code{\link[stats]{formula}}): the latter \strong{must} be specified explicitly using the \code{AsIs} operator (\code{\link{I}}). Covariates can be continuous, categorical, logical, or ordinal, but the response must always be continuous.
#'
#' While model selection in terms of choosing the optimal number of components and the GPCM/\pkg{mclust} model type is performed within \code{\link{MoE_clust}}, using one of the \code{criterion} options within \code{\link{MoE_control}}, choosing between multiple fits with different combinations of covariates or different initialisation settings can be done by supplying objects of class \code{"MoEClust"} to \code{\link{MoE_compare}}.
#' @note Where \code{BIC}, \code{ICL}, \code{AIC}, \code{LOGLIK}, \code{DF} and \code{ITERS} contain \code{NA} entries, this corresponds to a model which was not run; for instance a VVV model is never run for single-component models as it is equivalent to EEE. As such, one can consider the value as not really missing, but equivalent to the EEE value. \code{BIC}, \code{ICL}, \code{AIC}, \code{LOGLIK}, \code{DF} and \code{ITERS} all inherit the classes \code{"MoECriterion"} and \code{"mclustBIC", "mclustICL", etc.}, for which dedicated printing and plotting functions exist, respectively.
#'
#' @seealso \code{\link{MoE_compare}}, \code{\link{plot.MoEClust}}, \code{\link{predict.MoEClust}}, \code{\link{MoE_control}}, \code{\link{as.Mclust}}, \code{\link{MoE_crit}}, \code{\link{MoE_estep}}, \code{\link{MoE_cstep}}, \code{\link{MoE_dens}}, \code{\link[mclust]{mclustModelNames}}, \code{\link[mclust]{mclustVariance}}, \code{\link{expert_covar}}, \code{\link{aitken}}, \code{\link{I}}
#' @export
#' @references K. Murphy and T. B. Murphy (2018). Gaussian Parsimonious Clustering Models with Covariates. \emph{To appear}. <\href{https://arxiv.org/abs/1711.05632}{arXiv:1711.05632}>.
#'
#' C. Fraley and A. E. Raftery (2002). Model-based clustering, discriminant analysis, and density estimation. \emph{Journal of the American Statistical Association}, 97:611-631.
#' @author Keefe Murphy - <\email{keefe.murphy@@ucd.ie}>
#' @keywords clustering main
#' @usage
#' MoE_clust(data,
#'           G = 1:9,
#'           modelNames = NULL,
#'           gating = ~1,
#'           expert = ~1,
#'           control = MoE_control(...),
#'           network.data = NULL,
#'           ...)
#' @examples
#' \dontrun{
#' data(ais)
#' hema  <- ais[,3:7]
#' sex   <- ais$sex
#' BMI   <- ais$BMI
#'
#' # Fit a standard finite mixture model
#' m1    <- MoE_clust(hema, G=2:3)
#'
#' # Allow covariates to enter the mixing proportions
#' m2    <- MoE_clust(hema, G=2:3, gating= ~ sex + BMI)
#'
#' # Allow covariates to enter the component densities
#' m3    <- MoE_clust(hema, G=2:3, expert= ~ sex)
#'
#' # Allow covariates to enter both the gating & expert network
#' m4    <- MoE_clust(hema, G=2:3, gating= ~ BMI, expert= ~ sex)
#' 
#' # Fit an equal mixing proportion model with an expert network covariate
#' m5    <- MoE_clust(hema, G=2:3, expert= ~ sex, equalPro=TRUE)
#' 
#' # Fit models with gating covariates & an additional noise component
#' m6    <- MoE_clust(hema, G=2:3, tau0=0.1, gating=~BMI, network.data=ais)
#'
#' # Extract the model with highest ICL
#' (comp <- MoE_compare(m1, m2, m3, m4, m5, m6, criterion="icl"))
#' (best <- comp$optimal)
#' (summ <- summary(best))
#'
#' # Examine the gating and expert networks in greater detail
#' # (but refrain from inferring statistical significance!)
#' summary(best$gating)
#' summary(best$expert)
#'
#' # Visualise the results, incl. the gating network and log-likelihood
#' plot(best, what="gpairs")
#' plot(best, what="gating")
#' plot(best, what="loglik")
#'
#' # Visualise the results using the 'lattice' library
#' require("lattice")
#' z     <- factor(best$classification, labels=paste0("Cluster", seq_len(best$G)))
#' splom(~ hema | sex, groups=z)
#' splom(~ hema | z, groups=sex)}
  MoE_clust       <- function(data, G = 1:9, modelNames = NULL, gating = ~1, expert = ~1, control = MoE_control(...), network.data = NULL, ...) {

  # Definitions and storage set-up
    call          <- match.call()
    multi         <- missing(modelNames)
    gate.x        <- !missing(gating)
    exp.x         <- !missing(expert)
    criterion     <- control$criterion
    stopaX        <- control$stopping == "aitken"
    init.z        <- control$init.z
    nstarts       <- switch(EXPR=init.z, random=control$nstarts, 1L)
    startseq      <- seq_len(nstarts)
    multstart     <- nstarts > 1
    algo          <- control$algo
    exp.init      <- control$exp.init
    do.joint      <- exp.init$joint   && init.z != "random"
    max.init      <- exp.init$max.init
    drop.exp      <- exp.init$drop.break
    tol           <- control$tol[1L]
    g.reltol      <- control$tol[3L]
    max.it        <- control$itmax[1L]
    max.it        <- ifelse(max.it == .Machine$integer.max, max.it, max.it + 2L)
    g.itmax       <- control$itmax[3L]
    MaxNWts       <- control$MaxNWts
    equalPro      <- control$equalPro
    noise.args    <- control$noise.args
    noise         <- noise.args$noise.init
    tau0          <- noise.args$tau0
    noise.gate    <- noise.args$noise.gate
    noise.vol     <- noise.args$noise.vol
    noise.meth    <- noise.args$noise.meth
    noise.meth    <- ifelse(is.null(noise.meth), ifelse(is.null(noise.vol), "hypvol", "manual"), noise.meth)
    hc.args       <- control$hc.args
    hcName        <- hc.args$hc.meth
    hcUse         <- hc.args$hcUse
    km.args       <- control$km.args
    kiters        <- km.args$kiters
    kstarts       <- km.args$kstarts
    init.crit     <- control$init.crit
    warnit        <- control$warn.it
    itwarn        <- warnit > 2L
    verbose       <- control$verbose
    miss.init     <- control$miss.init
    miss.hc       <- control$miss.hc
    ctrl          <- list(equalPro=control$equalPro, noise.gate=(is.null(noise.args$noise.gate) || isTRUE(noise.args$noise.gate)), equalNoise=(!is.null(noise.args$equalNoise) && isTRUE(noise.args$equalNoise)))
    control       <- control[names(control) %in% c("eps", "tol", "itmax", "equalPro")]
    control$itmax <- control$itmax[-3L]
    control$tol   <- control$tol[-3L]
    netmiss       <- missing(network.data)
    if(!multi     &&
       !all(is.character(modelNames)))            stop("'modelNames' must be a vector of character strings", call.=FALSE)
    if(gate.x     &&
       !inherits(gating, "formula"))              stop("'gating' must be a formula", call.=FALSE)
    if(exp.x      &&
       !inherits(expert, "formula"))              stop("'expert' must be a formula", call.=FALSE)
    if(missing(data))                             stop("'data' must be supplied!",   call.=FALSE)

    tmp.nam       <- as.character(substitute(data))
    data          <- as.data.frame(data)
    num.X         <- vapply(data, is.numeric, logical(1L))
    if(anyNA(data))    {
      if(isTRUE(verbose))                         message("Rows with missing values removed from data\n")
      comp.x      <- stats::complete.cases(data)
      data        <- data[comp.x,, drop=FALSE]
    } else comp.x <- TRUE
    if(sum(num.X) != ncol(data))    {
      if(isTRUE(verbose))                         message("Non-numeric columns removed from data\n")
      data        <- data[,num.X,  drop=FALSE]
    }
    X             <- as.matrix(data)
    n             <- nrow(X)
    if((d         <- ncol(X))  == 0L)             stop("'data' is empty!", call.=FALSE)
    Nseq          <- seq_len(n)
    comp.x        <- if(isTRUE(comp.x)) Nseq else comp.x
    noise.null    <- all(nnull <- is.null(noise), tnull <- is.null(tau0))
    equalNoise    <- !noise.null      && ctrl$equalNoise
    gate.noise    <- (!noise.null     && noise.gate)    || noise.null
    if(missing(G) && !noise.null) G   <- 0L:9L
    if(any(G      != floor(G))        &&
       any(G < ifelse(noise.null, 1, 0)))         stop(paste0("'G' must be ", ifelse(noise.null, "strictly positive", "strictly non-negative when modelling with a noise-component")), call.=FALSE)
    if(any(G      >= n))        {
      G           <- G[G <= n]
      if(length(G) > 1)         {                 warning("Removing G values >= the number of observations\n",  call.=FALSE, immediate.=TRUE)
      } else                                      stop("G values must be less than the number of observations", call.=FALSE)
    }

    mod.fam       <- mclust.options("emModelNames")
    range.G       <- sort(as.integer(unique(G)))
    anyg0         <- any(G == 0)
    allg0         <- all(G == 0)
    anyg1         <- any(G == 1)
    anyg0or1      <- anyg0  + anyg1
    if(anyg0 || !noise.null)    {
      if(!is.null(noise.vol))   {
        Vinv      <- ifelse(isTRUE(attr(noise.vol, "Inverse")), noise.vol, 1/noise.vol)
      } else if(n  > d)         {
        NoiseVol  <- noise_vol(X, noise.meth, reciprocal=TRUE)
        Vinv      <- NoiseVol$vol
        NoiseLoc  <- NoiseVol$loc
      } else                                      stop("'noise.args$noise.vol' must be specified directly for high-dimensional data when a noise component is included", call.=FALSE)
    }         
    len.G         <- length(range.G)
    Gall          <- ifelse(noise.null, all(G > 1), all(G[G > 0] > 1))
    Gany          <- ifelse(noise.null, any(G > 1), any(G[G > 0] > 1))
    if((uni <- d  == 1L))       {
      mfg         <- c("E", "V")
      mf1         <- mf0       <- "E"
      colnames(X) <- tmp.nam[length(tmp.nam)]
    } else        {
      mf0         <- "EII"
      if((low.dim <- n > d))    {
        mfg       <- mod.fam
        mf1       <- c("EII", "EEI", "EEE")
      } else      {
        mfg       <- mod.fam[seq_len(6L)]
        mf1       <- c("EII", "EEI")
      }
    }
    low.dim       <- !uni && low.dim
    x.names       <- colnames(X)
    if(!multi)    {
      mNs         <- toupper(modelNames)
      if(any(sX   <- grepl("X",     mNs)))      {
       mNs        <- gsub("X", "E", mNs)
       if(verbose &&
          all(is.element(mNs,             mfg)))  message(paste0("'modelNames' which contain 'X' coerced to ", paste(shQuote(mNs[sX]), collapse=" + "), "\n"))
      }
      if(Gany     && any(!is.element(mNs, mfg)))  stop(paste0("Invalid 'modelNames'", ifelse(uni, " for univariate data", ifelse(low.dim, "", " for high-dimensional data")), "!"), call.=FALSE)
      if(!Gall)   {
        if(any(sZ <- !is.element(mNs,     mf1))){
          mf1     <- tryCatch(unname(vapply(mNs,  function(x)  switch(EXPR=x, E=, V="E", EII=, VII="EII", EEI=, VEI=, EVI=, VVI="EEI", EEE=, EVE=, VEE=, VVE=, EEV=, VEV=, EVV=, VVV="EEE"), character(1L))),
                              error=function(e) { e$message <- paste0("Invalid 'modelNames' for single component models", ifelse(uni, " for univariate data", ifelse(low.dim, "", " for high-dimensional data")), "!")
                                                  stop(e, call.=FALSE) } )
          if(isTRUE(verbose))                     message(paste0("'modelNames'", ifelse(any(sX), " further", ""), " coerced from ", paste(shQuote(mNs[sZ]), collapse=" + "), " to ", paste(shQuote(mf1[sZ]), collapse=" + "), " where G=1\n"))
        }
      }
      mf1         <- mfg       <- mNs
    }
    mf1           <- unique(mf1)
    mfg           <- unique(mfg)
    all.mod       <- if(all(multi, !uni, Gany)) mclust.options("emModelNames") else unique(c(if(anyg0) mf0, if(any(G == 1)) mf1, if(any(G > 1)) mfg))
    multi         <- length(all.mod)    > 1L
    BICs          <- ICLs      <-
    AICs          <- DF.x      <- IT.x <- provideDimnames(matrix(NA, nrow=len.G, ncol=length(all.mod)), base=list(as.character(range.G), all.mod))
    LL.x          <- replicate(nstarts  + 1L, list(BICs))
    LL.x[[1L]][]  <- -Inf
    crit.tx       <- crit.gx   <- -sqrt(.Machine$double.xmax)
    if(!netmiss)   {
      netdat      <- network.data
      if(!is.data.frame(netdat))                  stop("'network.data' must be a data.frame if supplied",          call.=FALSE)
    }

  # Define the gating formula
    if(allg0 && gate.x)         { if(verbose)     message("Can't include gating network covariates in a noise-only model\n")
      gate.x      <- FALSE
    }
    gate.G        <- ifelse((range.G + !noise.null) > 1, gate.x, FALSE)
    if(gate.x)    {
      if(inherits(try(stats::terms(gating), silent=TRUE), "try-error")) {
        if(netmiss)                               stop("Can't use '.' in 'gating' formula without supplying 'network.data' argument", call.=FALSE)
        gating    <- setdiff(attr(stats::terms(gating, data=network.data), "term.labels"), x.names)
        gating    <- stats::reformulate(if(length(gating) == 0) "1" else gating, response="z")
      }
      gating      <- tryCatch(stats::update.formula(stats::as.formula(gating), zN ~ .),
                              error=function(e)   stop("Invalid 'gating' network formula supplied", call.=FALSE))
      environment(gating)      <- environment()
      if(gating[[3L]]   == 1)   { if(verbose)     message("Not including gating network covariates with only intercept on gating formula RHS\n")
        gate.x    <- FALSE
        gate.G    <- rep(gate.x, len.G)
      }
      Gn          <- G + !noise.null - !gate.noise
      if(gate.x   &&
         any(Gn   <= 1))        {
        if(all(Gn <= 1) && verbose)               message(paste0("Can't include gating network covariates ", ifelse(gate.noise, "in a single component mixture", "where G is less than 3 when 'noise.args$noise.gate' is FALSE\n")))
        gate.G[Gn <= 1]        <- FALSE
      }
      gate.names  <- stats::terms(gating)
      gate.names  <- labels(gate.names)[attr(gate.names, "order") <= 1]
    } else gating <- stats::as.formula(zN ~ 1)
    gate.noise    <- ifelse(gate.G, gate.noise, TRUE)
    if(equalPro   && gate.x)    { if(verbose)     message("Can't constrain mixing proportions to be equal when gating covariates are supplied\n")
      equalPro    <- FALSE
    }
    equal.tau     <- ifelse((range.G + !noise.null) == 1, TRUE, equalPro)   & !gate.G
    equal.noise   <- ifelse((range.G + !noise.null) == 1, TRUE, equalNoise) & equal.tau

  # Define the expert formula
    if(allg0 && exp.x)          { if(verbose)     message("Can't include expert network covariates in a noise-only model\n")
      exp.x       <- FALSE
    }
    if(exp.x)     {
      if(inherits(try(stats::terms(expert), silent=TRUE), "try-error"))   {
        if(netmiss)                               stop("Can't use '.' in 'expert' formula without supplying 'network.data' argument", call.=FALSE)
        expert    <- setdiff(attr(stats::terms(expert, data=network.data), "term.labels"), x.names)
        expert    <- stats::reformulate(if(length(expert) == 0) "1" else expert, response="X")
      }
      expert      <- tryCatch(stats::update.formula(stats::as.formula(expert), X ~ .),
                              error=function(e)   stop("Invalid 'expert' network formula supplied", call.=FALSE))
      environment(expert)      <- environment()
      if(expert[[3L]]   == 1)   { if(verbose)     message("Not including expert network covariates with only intercept on expert formula RHS\n")
        exp.x     <- FALSE
      }
      expx.names  <- stats::terms(expert)
      expx.names  <- labels(expx.names)[attr(expx.names, "order") <= 1]
    } else expert <- stats::as.formula(X ~ 1)
    
  # More managing of noise component
    if(!tnull)    {
      if(length(tau0)  > 1)     {
        if(all(gate.x, 
               noise.gate))     {
          if(length(tau0)      != n)              stop(paste0("'tau0' must be a scalar or a vector of length N=", n), call.=FALSE)
        } else                                    stop("'tau0' must be a scalar in the interval (0, 1)", call.=FALSE)
      }
    }
    if(!nnull)      {
      if(length(noise)         != n)              stop(paste0("'noise.args$noise.init' must be a vector of length N", n), call.=FALSE)
      if(!is.logical(noise))    {
        if(any(match(noise, Nseq,
               nomatch=0)      == 0))             stop("Numeric 'noise.args$noise.init' must correspond to row indices of data", call.=FALSE)
        noise      <- as.logical(match(Nseq, noise, nomatch=0))
      }
      if(!tnull)    {
        tau0       <- tau0 * noise
        noise      <- vector("logical", n)
        nnoise     <- 0L
        noisen     <- n
      } else  {
        nnoise     <- sum(as.numeric(noise))
        noisen     <- n - nnoise
        if(any(G    > noisen)) range.G <- range.G[range.G <= noisen]
      }
    } else   {
      noise       <- vector("logical", n)
      nnoise      <- 0L
      noisen      <- n
    }
    if(allg0)      {
      noise.null  <- FALSE
      noise       <- rep(TRUE, n)
      nnoise      <- n
      noisen      <- 0L
    }

  # Tell network formulas where to look for variables
    if(gate.x)     {
      gate.covs   <- eval(bquote(stats::model.frame(.(stats::update.formula(gating, NULL ~ .)), data=.(call$network.data), drop.unused.levels=TRUE)), envir=parent.frame(), enclos=environment())
      gate.names  <- colnames(gate.covs)
      gate.covs   <- cbind(gate.covs, eval(bquote(stats::model.frame(.(as.formula(paste("~", paste(eval(bquote(all.vars(.(gating))), envir=parent.frame())[-1L], collapse="+")))), data=.(call$network.data), drop.unused.levels=TRUE)), envir=parent.frame(), enclos=environment()))
      gate.covs   <- gate.covs[,unique(colnames(gate.covs)), drop=FALSE]
      netdat      <- gate.covs
    } 
    if(exp.x)      {
      expx.covs   <- eval(bquote(stats::model.frame(.(stats::update.formula(expert, NULL ~ .)), data=.(call$network.data), drop.unused.levels=TRUE)), envir=parent.frame(), enclos=environment())
      expx.names  <- colnames(expx.covs)
      expx.covs   <- cbind(expx.covs, eval(bquote(stats::model.frame(.(as.formula(paste("~", paste(eval(bquote(all.vars(.(expert))), envir=parent.frame())[-1L], collapse="+")))), data=.(call$network.data), drop.unused.levels=TRUE)), envir=parent.frame(), enclos=environment()))
      expx.covs   <- expx.covs[,unique(colnames(expx.covs)), drop=FALSE]
      netdat      <- if(netmiss || !gate.x) expx.covs else cbind(netdat, expx.covs)
      netdat      <- netdat[,unique(colnames(netdat)), drop=FALSE]
    }
    gate.names    <- if(gate.x) gate.names   else NA
    expx.names    <- if(exp.x)  expx.names   else NA
    netnames      <- unique(c(gate.names, expx.names))
    netnames      <- netnames[!is.na(netnames)]
    if(!netmiss)   {
      if(!all(netnames %in% colnames(netdat)))    stop("Supplied covariates not found in supplied 'network.data'", call.=FALSE)
      netdat      <- netdat[comp.x,netnames, drop=FALSE]
      gate.covs   <- if(gate.x)   gate.covs       else as.data.frame(matrix(0L, nrow=n, ncol=0L))
      expx.covs   <- if(exp.x)    expx.covs       else as.data.frame(matrix(0L, nrow=n, ncol=0L))
    } else {
      if(any(grepl("\\$", netnames)))             stop("Don't supply covariates to gating or expert networks using the $ operator: use the 'network.data' argument instead", call.=FALSE)
      gate.covs   <- if(gate.x)   stats::model.frame(gating[-2L], drop.unused.levels=TRUE)[comp.x,, drop=FALSE] else as.data.frame(matrix(0L, nrow=n, ncol=0L))
      expx.covs   <- if(exp.x)    stats::model.frame(expert[-2L], drop.unused.levels=TRUE)[comp.x,, drop=FALSE] else as.data.frame(matrix(0L, nrow=n, ncol=0L))
    }
    if(nrow(gate.covs) != n)                      stop("'gating' covariates must contain the same number of rows as 'data'", call.=FALSE)
    if(nrow(expx.covs) != n)                      stop("'expert' covariates must contain the same number of rows as 'data'", call.=FALSE)
    glogi         <- vapply(gate.covs, is.logical, logical(1L))
    elogi         <- vapply(expx.covs, is.logical, logical(1L))
    gate.covs[,glogi]          <- sapply(gate.covs[,glogi], as.factor)
    expx.covs[,elogi]          <- sapply(expx.covs[,elogi], as.factor)
    if(netmiss)    {
      netdat      <- cbind(gate.covs, expx.covs)
      netnames    <- unique(colnames(netdat))
      netdat      <- data.frame(if(ncol(netdat) > 0) netdat[,netnames, drop=FALSE]      else netdat, stringsAsFactors=TRUE)
      colnames(netdat)         <- netnames
    } else if(any(nlogi        <- unique(c(glogi, elogi))))      {
      netdat[,nlogi]           <- sapply(netdat[,nlogi],    as.factor)
    }
    attr(netdat, "Gating")     <- gate.names
    attr(netdat, "Expert")     <- expx.names
    attr(netdat, "Both")       <- if(length(intersect(gate.names, expx.names)) == 0) NA else intersect(gate.names, expx.names)
    if(!identical(gating,
       drop_constants(gate.covs, gating)))        stop("Constant columns exist in gating formula; remove offending gating covariate(s) and try again", call.=FALSE)
    if(!identical(expert,
       drop_constants(expx.covs, expert)))        stop("Constant columns exist in expert formula; remove offending expert covariate(s) and try again", call.=FALSE)
    jo.cts        <- names(which(!vapply(expx.covs, is.factor, logical(1L))))
    nct           <- length(jo.cts)
    g.range       <- range.G[range.G   > 1]
    do.joint      <- do.joint  && nct >= 1L
    XI            <- (if(exp.x && do.joint) cbind(X, expx.covs[,jo.cts, drop=FALSE]) else X)[!noise,, drop=FALSE]
    if(someG      <- !all(G    == 1)  && !allg0)   {
      init.var    <- ifelse(do.joint  && !allg0, d + nct, d)
      highd       <- init.var  >= n
      multv       <- init.var   > 1
      if(!multv)   {
        init.z    <- ifelse(miss.init, "quantile", init.z)
      } else if(init.z         == "quantile")     stop("Quantile-based initialisation of the allocations is only permitted for univariate data without expert network covariates", call.=FALSE)
      if(!any(gate.x, exp.x))   {
       if(verbose && isTRUE(exp.init$clustMD))    message("'exp.init$clustMD' not invoked - no covariates included!")
       exp.init$clustMD        <- FALSE
       if(init.z  == "mclust")  { if(verbose)     message("Initialisation method coerced from \"mclust\" to \"hc\" as there are no gating/expert network covariates\n")
         init.z   <- "hc"
       }
     }
    }
    if(clust.MD   <- exp.init$clustMD) {
      if(exp.init$joint)        {
        exp.fac   <- ncol(expx.covs)   - nct
        if((do.md <- exp.fac    > 0))  {
          if(!(has.md         <- suppressMessages(requireNamespace("clustMD",
                                 quietly=TRUE)))) warning("'exp.init$clustMD' not invoked - 'clustMD' library not loaded\n", call.=FALSE, immediate.=TRUE)
        } else if(isTRUE(verbose))                message("'exp.init$clustMD' not invoked - no categorical or ordinal expert network covariates\n")
      }   else if(isTRUE(verbose))                message("'exp.init$clustMD' not invoked - exp.init$joint not set to TRUE\n")
    }
    if((mdind     <- clust.MD  && do.md && has.md) && someG && exp.init$joint) {
      expx.facs   <- expx.covs[,setdiff(colnames(expx.covs), jo.cts), drop=FALSE]
      flevs       <- vapply(expx.facs, nlevels,           integer(1L))
      b.ind       <- which(flevs == 2L)
      o.ind       <- which(vapply(expx.facs, is.ordered, logical(1L)))
      n.ind       <- setdiff(seq_len(ncol(expx.facs)),  c(b.ind, o.ind))
      XY          <- cbind(X,  expx.covs[,jo.cts, drop=FALSE])
      XY          <- cbind(XY, vapply(expx.facs[,c(b.ind, o.ind, n.ind), drop=FALSE], as.numeric, numeric(n)))[!noise,, drop=FALSE]
      J           <- ncol(XY)
      CnsIndx     <- d + nct
      OrdIndx     <- J - length(n.ind)
      has.pkg     <- suppressMessages(requireNamespace("snow", quietly=TRUE))
      if(!has.pkg)                                stop("'snow' package note installed", call.=FALSE)
      mdx         <- utils::capture.output( {
        mds       <- clustMD::clustMDparallel(X=XY, G=g.range, CnsIndx=CnsIndx, OrdIndx=OrdIndx, Nnorms=25000, MaxIter=500, store.params=FALSE,
                                              model=c("EII", "VII", "EEI", "VEI", "EVI", "VVI", "BD"), autoStop=TRUE, stop.tol=1e-04, scale=FALSE,
                                              startCL=switch(EXPR=init.z, kmeans="kmeans", mclust="mclust", random="random", "hc_mclust")) })
      mdcrit      <- switch(EXPR=init.crit, bic=mds$BICarray, icl=mds$ICLarray)
      mderr       <- is.na(mdcrit)
      mdcrit      <- replace(mdcrit, mderr, -Inf)
      mdbest      <- rowMaxs(mdcrit, na.rm=TRUE)
    } else mderr  <- as.matrix(FALSE)
    if(is.element(init.z, c("hc", "mclust")) && someG) {
      if(miss.hc)  {
        hcName    <- ifelse(highd, "EII", "VVV")
      }
      if(multv    &&
         is.element(hcName, c("E",        "V")))  stop("'hc.args$hc.meth' can only be 'E' or 'V' for univariate data without expert network covariates", call.=FALSE)
      if(highd    && !is.element(hcName, c("EII",
                                "VII",  "EEE")))  warning("Consider a diagonal 'EII' or 'VII' model (or equal volume 'EEE' model) for 'hc.args$hc.meth' for initialising allocations for high-dimensional data\n", call.=FALSE)
      if(!multv   && !is.element(hcName, c("VVV",
                                "E",      "V")))  warning("Possibly invalid 'hc.args$hc.meth' for univariate data\n", call.=FALSE)
      Zhc         <- tryCatch(hc(XI, modelName=hcName, use=hcUse, minclus=min(g.range)), error=function(e) {
                     if(!mdind)                   stop(paste0("Hierarchical clustering initialisation failed",
                                                              ifelse(init.z == "hc", "", " (when initialising using 'mclust')")), call.=FALSE)
                                                  else try(stop(), silent=TRUE) })
      if(!(hcfail <- inherits(Zhc, "try-error")))       {
        if(init.z == "mclust")  {
          if(!mdind)            {
            mcarg <- list(data=XI, G=g.range, verbose=FALSE, control=emControl(equalPro=equalPro), initialization=list(hcPairs=Zhc))
            mcl   <- suppressWarnings(switch(EXPR=init.crit, icl=do.call(mclustICL, mcarg), bic=do.call(mclustBIC, mcarg)))
            mcerr <- apply(is.na(mcl), 1L, all)
            if(any(mcerr))                        stop(paste0("Mclust initialisation failed for the G=", paste(g.range[mcerr], collapse="/"), " model", ifelse(sum(mcerr) > 1, "s", "")), call.=FALSE)
            class(mcl)         <- "mclustBIC"
          }
          mcfail  <- rep(FALSE, len.G)
        } else if(init.z == "hc")             {
          hc1     <- any(range.G == 1)
          hcZ     <- hclass(Zhc, G=g.range)
        }
      }
    } else hcfail <- mcfail    <- FALSE

  # Loop over range of G values and initialise allocations
    G.last        <- range.G[len.G]
    for(g in range.G) {
      if(isTRUE(verbose))   {     cat(paste0("\n", g, " cluster model", ifelse(multi, "s", ""), " -\n"))
        last.G    <- g == G.last
      }
      x.dat       <- replicate(max(g, 1L), X, simplify=FALSE)
      h           <- which(range.G  == g)
      equal.pro   <- equal.tau[h]
      n0pro       <- equal.noise[h]
      gate.g      <- gate.G[h]
      exp.g       <- exp.x && g > 0
      init.exp    <- exp.g && exp.init$mahalanobis
      Gseq        <- seq_len(g)
      gN          <- max(g  + !noise.null, 1L)
      z           <- matrix(0L, n, gN)
      stopGx      <- g > 1 && stopaX
      noiseG      <- !noise.null || g == 0
      noise.gate  <- ifelse(noiseG, gate.noise[h], TRUE)
      Xinv        <- if(noiseG) Vinv

    # Initialise Expert Network & Allocations
      if(indmd    <- mdind &&   g > 1) {
        h2        <- h - anyg0or1
        mdh       <- mderr[h2,]
        if(all(mdh))              {              warning(paste0("\tInvoking 'exp.init$clustMD' failed for ALL clustMD model types where G=",  g, ",\n\t\tprobably due to the presence of nominal expert network covariates,\n\t\tor the 'init.z' method used to initialise the call to clustMD\n"), call.=FALSE, immediate.=TRUE)
        } else     {
          if(length(bmd        <- which(mdcrit[h2,] == mdbest[h2])) > 1)  {
            dfmd  <- .npars_clustMD(D=ifelse(length(n.ind) > 0, OrdIndx + sum(flevs - 1L), J), G=g, J=J, CnsIndx=CnsIndx, OrdIndx=OrdIndx, K=flevs)
            bmd   <- which(dfmd  == min(dfmd[bmd]))
          }
          z.tmp   <- unmap(mds$Output[[(bmd - 1L) * (len.G - anyg0or1)  + h2]]$cl, groups=Gseq)
          if(any(mdh))                           warning(paste0("\tInvoking 'exp.init$clustMD' failed for SOME clustMD model types where G=", g, ",\n\t\tprobably due to the presence of nominal expert network covariates,\n\t\tor the 'init.z' method used to initialise the call to clustMD\n"), call.=FALSE, immediate.=TRUE)
        }
      }
      if(!indmd   || all(mdh))    {
        if(hcfail)   next
        if(g > 1)  {
          if(all(indmd, init.z == "mclust"))   {
            mcarg <- list(data=XI, G=g, verbose=FALSE, control=emControl(equalPro=equal.pro), initialization=list(hcPairs=Zhc))
            mcl   <- suppressWarnings(switch(EXPR=init.crit,   icl=do.call(mclustICL, mcarg),   bic=do.call(mclustBIC, mcarg)))
            if(mcfail[h]       <- all(is.na(mcl))) next
            class(mcl)         <- "mclustBIC"
          }
          if(multstart)         {
            zg        <- replicate(nstarts, list(unmap(sample(x=Gseq, size=noisen, replace=TRUE), groups=Gseq)))
          } else       {
            z.tmp     <- unmap(switch(EXPR     = init.z,
                                      hc       = hcZ[,h - anyg0 - hc1],
                                      kmeans   = stats::kmeans(x=XI, centers=g, iter.max=kiters, nstart=kstarts)$cluster,
                                      mclust   = suppressWarnings(Mclust(XI, G=g, x=mcl))$classification,
                                      quantile = quant_clust(x=XI, G=g),
                                      random   = sample(x=Gseq, size=noisen, replace=TRUE)), groups=Gseq)
          }
        } else     {
          z.tmp   <- unmap(rep(1L, ifelse(noisen == 0 || g == 0, n, noisen)), groups=1L)
        }
      }
    for(i in if(g  > 1) startseq else 1L)      {
      if(isTRUE(multstart)     && g > 1)       {
        if(isTRUE(verbose))                      cat(paste0("\n\tRandom Start: #", i, "...\n"))
        z.tmp     <- zg[[i]]
      }
      zc          <- ncol(z.tmp)
      if(zc != g  && g > 0)     {
        z.tmp     <- cbind(z.tmp, matrix(0L, ncol=g - zc, nrow=n))
      }

    # Invoke Iterative Mahalanobis Step?
      if(exp.g)    {
        z.mat     <- z.alloc   <- matrix(0L, nrow=n * g,  ncol=g)
        muX       <- if(uni)      vector("numeric",   g)  else matrix(0L, nrow=d, ncol=g)
      } else {
        exp.pen   <- g * d
      }
      expold      <- init.exp
      if(init.exp) {
        tmp.z     <- matrix(NA, nrow=ifelse(noisen == 0, n, noisen), ncol=g)
        mahala    <- res.G     <- Efit <- list()
        xN        <- X[!noise,, drop=FALSE]
        expnoise  <- expx.covs[!noise,, drop=FALSE]
        expN      <- stats::update.formula(expert, xN ~ .)
        ix        <- 0L
        ne        <- ncol(expnoise)
        while(!identical(tmp.z, z.tmp) && ix <= max.init) {
          tmp.z   <- z.tmp
          ix      <- ix  + 1L
          for(k in Gseq) {
            sub   <- z.tmp[,k] == 1
            exp   <- tryCatch(stats::lm(expN, data=expnoise, subset=sub),              error=function(e) try(if(drop.exp) stop("DROP") else stats::lm(drop_constants(expnoise, expN, sub), data=expnoise, subset=sub), silent=TRUE))
            if(inherits(exp,        "try-error")) {
              init.exp         <- FALSE
              break
            } else Efit[[k]]   <- exp
            pred  <- tryCatch(suppressWarnings(stats::predict(exp, newdata=expnoise)), error=function(e) try(if(drop.exp) stop("DROP") else stats::predict(exp, newdata=drop_levels(exp, expnoise)), silent=TRUE))
            if(inherits(pred,       "try-error")) {
              init.exp         <- FALSE
            } else {
              pred             <- as.matrix(pred)
              if(sum(pna       <- !stats::complete.cases(pred)) >= 1) {
                if(drop.exp)    {
                  init.exp     <- FALSE
                } else   {
                  nexp         <- expnoise[,vapply(seq_len(ne), function(p, ep=expnoise[,p], fep=is.factor(ep)) {
                                 (!fep && !all(ep[sub] == ep[sub][1], na.rm=TRUE)) || (fep && (nlevels(droplevels(ep[sub])) == nlevels(ep))) }, logical(1L)), drop=FALSE]
                  px           <- try(stats::predict(stats::lm(drop_constants(nexp, expN, pna), data=nexp, subset=pna)), silent=TRUE)
                  if(!inherits(px,  "try-error")) {
                    pred[pna,] <- px
                  } else   {
                    init.exp   <- FALSE
                  }
                }
              }
              res              <- xN - pred
              res.G[[k]]       <- res
              mahala[[k]]      <- if(g > 1) MoE_mahala(exp, res, squared=TRUE)
            }
          }
          if(!init.exp) {
            break
          } else if(g   > 1)    {
            maha  <- do.call(cbind, mahala)
            if(anyNA(maha))     {
              init.exp         <- FALSE
              break
            } else   z.tmp     <- maha == rowMins(maha)
          }
        }
        if(ix     >= max.init)                    warning(paste0("\tMahalanobis initialisation step failed to converge in max.init=", max.init, " iterations for the ", g, " cluster models\n"), call.=FALSE, immediate.=TRUE)
        if(noiseG && init.exp  && tnull)  {
         nRG      <- replicate(g, matrix(NA, nrow=n, ncol=d), simplify=FALSE)
         nX       <- X[noise,, drop=FALSE]
         noisexp  <- expx.covs[noise,, drop=FALSE]
         for(k in Gseq) {
           nRG[[k]][!noise,]   <- res.G[[k]]
           exp    <- Efit[[k]]
           pred   <- tryCatch(stats::predict(exp, newdata=noisexp),       error=function(e) try(if(drop.exp) stop("DROP") else stats::predict(exp, newdata=drop_levels(exp, noisexp)), silent=TRUE))
           if(inherits(pred,        "try-error")) {
             init.exp          <- FALSE
           } else  {
             pred              <- as.matrix(pred)
             if(sum(pna        <- !stats::complete.cases(pred)) >= 1) {
               if(drop.exp)    {
                 init.exp      <- FALSE
               } else   {
                 sub           <- z.tmp[,k] == 1
                 nexp          <- expnoise[,vapply(seq_len(ne), function(p, ep=expnoise[,p], fep=is.factor(ep)) {
                                 (!fep && !all(ep[sub] == ep[sub][1], na.rm=TRUE)) || (fep && (nlevels(droplevels(ep[sub])) == nlevels(ep))) }, logical(1L)), drop=FALSE]
                 px            <- try(stats::predict(stats::lm(drop_constants(nexp, expN, which(noise)[pna]), data=nexp, subset=sub), newdata=noisexp[pna,colnames(nexp), drop=FALSE]), silent=TRUE)
                 if(!inherits(px,   "try-error")) {
                   pred[pna,]  <- px
                 } else {
                   init.exp    <- FALSE
                 }
               }
             }
             nRG[[k]][noise,]  <- nX - pred
           }
         }
         res.G    <- nRG
        }
        G.res     <- if(uni) as.matrix(do.call(base::c, res.G)) else do.call(rbind, res.G)
      }
      if(eNO      <- expold    != init.exp)       warning(paste0("\tExtra initialisation step with expert covariates failed where G=", g, ifelse(drop.exp, ": try setting 'drop.exp' to FALSE\n", ", even with 'drop_constants' and 'drop_levels' invoked:\n\t\tTry suppressing the initialisation step via 'exp.init$mahalanobis' or using other covariates\n")), call.=FALSE, immediate.=TRUE)

    # Account for Noise Component
      z.tmp       <- 0L + z.tmp
      if(noise.null) {
        z         <- z.init    <- z.tmp
      } else   {
        if(g   > 0)  {
          z[!noise, -gN]       <- z.tmp
          if(tnull)  {
            z[noise, gN]       <- 1L
          } else z             <- cbind(z[,-gN] * (1 - tau0), tau0)
        } else {
          z[]     <- 1L
        }
        z.init    <- z
      }
      if(all(noiseG, gate.g, !noise.gate, gN  > 1, algo != "EM")) {
        z[,-gN]   <- replace(z[,-gN], z[,-gN] > 0, 1L)
        z.init    <- z
      }
      if(noise.gate)      {
        zN        <- z
      } else   {
        zN        <- z[,-gN,     drop=FALSE]
        zN        <- zN[!noise,, drop=FALSE]
      }
      if(init.exp) {
        for(k in Gseq) z.alloc[(k - 1L) * n + Nseq,k] <- z.init[,k]
      }
      col.z       <- colSums2(z.init)
      emptyinit   <- FALSE
      if(any(col.z[Gseq]  < 1))   {               warning(paste0("\tFor the ", g, " component models, ", ifelse(gN > 1, "one or more", ""), " components were empty after initialisation\n"),          call.=FALSE, immediate.=TRUE)
        emptyinit <- TRUE
      } else if(any(col.z[Gseq]   < 2))           warning(paste0("\tFor the ", g, " component models, ", ifelse(gN > 1, "one or more", ""), " components were initialised with only 1 observation\n"), call.=FALSE, immediate.=TRUE)

    # Initialise gating network
      if(gate.g)  {
        if(noise.gate)    {
          g.init  <- multinom(gating, trace=FALSE, data=gate.covs, maxit=g.itmax, reltol=g.reltol, MaxNWts=MaxNWts)
          tau     <- g.init$fitted.values
        } else    {
          if(all(!noise) && algo != "EM" && !tnull) {
            zN[zN  > 0]  <- 1L
          }
          g.init  <- multinom(gating, trace=FALSE, data=gate.covs[!noise,, drop=FALSE], maxit=g.itmax, reltol=g.reltol, MaxNWts=MaxNWts)
          tau     <- .tau_noise(stats::predict(g.init, type="probs", newdata=gate.covs), z[,gN])
        }
        gate.pen  <- length(stats::coef(g.init)) + ifelse(noise.null, 0L, 1L) + ifelse(noise.gate, 0L, 1L)
       #g.init    <- glmnet::cv.glmnet(y=z, x=model.matrix(gating, data=gate.covs)[,-1L], family="multinomial", type.multinomial="grouped")
       #tau       <- stats::predict(g.init, type="response", newx=model.matrix(gating, data=gate.covs)[,-1L], s="lambda.1se")[,,1]
       #gate.pen  <- g.init$glmnet.fit$df[which(g.init$glmnet.fit$lambda == g.init$lambda.1se)] + ifelse(noise.null, 1, 2)
        ltau      <- log(tau)
      } else      {
        if(equal.pro && noiseG && !n0pro) {
          t0      <- mean(z.init[,gN])
          tau     <- c(rep((1 - t0)/g, g), t0)
        } else    {
          tau     <- if(equal.pro) rep(1/gN, gN) else col.z/n
        }
        ltau      <- .mat_byrow(log(tau), nrow=n, ncol=gN)
        gate.pen  <- ifelse(equal.pro && g > 1, noiseG - n0pro, gN - 1L) + ifelse(noiseG, 1L, 0L)
      }
      ltau.init   <- ltau
      failedM     <- NULL
      expinitG    <- init.exp

    # Loop over the mclust model type(s)
      modtypes    <- if(g > 1)    mfg     else if(g == 1) mf1 else mf0
      T.last      <- modtypes[length(modtypes)]
      for(modtype in modtypes)  {
        m0W       <- m0X       <- ERR  <- FALSE

      # Initialise parameters from allocations
        if(isTRUE(verbose))     { cat(paste0("\n\tModel: ", modtype, "\n"))
          last.T  <- modtype   == T.last
        }
        x.df      <- ifelse(g   > 0, nVarParams(modtype, d, g), 0L) + gate.pen
        if(g > 0  && expinitG)  {
         Mstep    <- try(mstep(modtype, G.res, z.alloc, control=control), silent=TRUE)
         init.exp <- ifelse(inherits(Mstep, "try-error"), FALSE, attr(Mstep, "returnCode") >= 0)
        }
        if(expold != init.exp  && !eNO) {
          failedM <- c(failedM, modtype)
        }
        if(g > 0  && !init.exp) {
          Mstep   <- try(mstep(modtype, X, if(noise.null) z.init else z.init[,-gN, drop=FALSE], control=control), silent=TRUE)
          ERR     <- inherits(Mstep, "try-error")           ||   attr(Mstep, "returnCode")  < 0
        }
        if(g > 0  && !ERR)      {
          mus     <- if(init.exp) muX     else Mstep$parameters$mean
          vari    <- Mstep$parameters$variance
        } else     {
          mus     <- matrix(NA, nrow=n, ncol=0L)
          vari    <- list(modelName=modtype, d=d, G=0L)
        }
        alG       <- ifelse(gN  > 1, algo, "EM")

        medens    <- try(MoE_dens(data=if(init.exp) res.G else x.dat, mus=mus, sigs=vari, log.tau=ltau.init, Vinv=Xinv), silent=TRUE)
        if((ERR   <- ERR || ((g > 0 && attr(Mstep, "returnCode") < 0) || (inherits(medens, "try-error")) || any(medens > 0)))) {
          ll      <- NA
          j       <- 1L
          if(isTRUE(verbose))     cat(paste0("\t\t# Iterations: ", ifelse(ERR, "stopped at ", ""), j, ifelse(last.G && last.T, "\n\n", "\n")))
          BICs[h,modtype]      <-
          ICLs[h,modtype]      <-
          AICs[h,modtype]      <-
          DF.x[h,modtype]      <- -Inf
          IT.x[h,modtype]      <- Inf
          LL.x[[i  + 1L]][h,modtype]   <- -Inf
          next
        } else     {
          Estep   <- switch(EXPR=alG, EM=MoE_estep(Dens=medens), MoE_cstep(Dens=medens))
          z       <- zN        <- Estep$z
          if((ERR <- any(is.nan(z))))             next
          ll      <- c(-Inf, ifelse(gN <= 1 && !exp.g, Estep$loglik, -sqrt(.Machine$double.xmax)))
          j       <- 2L
          stX     <- gN    > 1 || exp.g
        }
        nW        <- 1L
        alW       <- switch(EXPR=alG, cemEM="CEM", alG)

      # Run the EM/CEM algorithm
        while(stX)    {

        # Expert network
          if(exp.g)   {
           e.fit  <- e.res     <- list()
           for(k in Gseq)  {
            fitE  <- stats::lm(expert, weights=z[,k], data=expx.covs)
           #fitE  <- glmnet::cv.glmnet(y=X, x=model.matrix(expert, data=expx.covs)[,-1L], weights=z[,k], family="mgaussian")
            e.fit[[k]]        <- fitE
            e.res[[k]]        <- stats::residuals(fitE)
           #e.res[[k]]        <- X - stats::predict(fitE, weights=z[,k], type="response", newx=model.matrix(expert, data=expx.covs)[,-1L], s="lambda.1se")
            z.mat[(k - 1L) * n + Nseq,k]      <- z[,k]
           }
           res.x  <- if(uni) as.matrix(do.call(base::c, e.res))  else do.call(rbind, e.res)
          }

        # M-step
          Mstep   <- try(if(exp.g) mstep(modtype, res.x, z.mat, control=control) else mstep(modtype, X, if(noise.null) z else z[,-gN, drop=FALSE], control=control), silent=TRUE)
          ERR     <- (inherits(Mstep, "try-error") || attr(Mstep, "returnCode")  < 0)
          if(isTRUE(ERR))  {
            z.err <- if(exp.g) z.mat     else if(noise.null)  z  else z[,-gN, drop=FALSE]
            if(any(colSums2(z.err) == 0))         warning(paste0("\tThere were empty components: ", modtype, " (G=", g, ")\n"), call.=FALSE)
          } else   {
            mus   <- if(exp.g) muX       else Mstep$parameters$mean
            vari  <- Mstep$parameters$variance
          }

        # Gating Network
          if(gate.g && !ERR)    {
            if(noise.gate)      {
             fitG <- multinom(gating, trace=FALSE, data=gate.covs, maxit=g.itmax, reltol=g.reltol, MaxNWts=MaxNWts)
             tau  <- fitG$fitted.values
            } else {
             zN   <- .renorm_z(z[,-gN, drop=FALSE])
             zN[is.nan(zN)]    <- .Machine$double.eps
             fitG <- multinom(gating, trace=FALSE, data=gate.covs, maxit=g.itmax, reltol=g.reltol, MaxNWts=MaxNWts)
             tau  <- .tau_noise(fitG$fitted.values, z[,gN])
            }
           #fitG  <- glmnet::cv.glmnet(y=z, x=model.matrix(gating, data=gate.covs)[,-1L], family="multinomial", type.multinomial="grouped")
           #tau   <- stats::predict(fitG, type="response", newx=model.matrix(gating, data=gate.covs)[,-1L], s="lambda.1se")[,,1]
            ltau  <- log(tau)
          } else  {
            if(equal.pro && !noise.null && !n0pro)  {
              t0  <- mean(z[,gN])
              tau <- c(rep((1 - t0)/g, g), t0)
            } else   if(!equal.pro)      {
              tau <- if(noise.null && !ERR)   Mstep$parameters$pro else colMeans2(z)
            }
            tau   <- if(!exp.g     || !noise.null  || ERR)     tau else tau/sum(tau)
            ltau  <- if(equal.pro  && (noise.null  || n0pro)) ltau else .mat_byrow(log(tau), nrow=n, ncol=gN)
          }

        # E-step & record log-likelihood
          if(!ERR) {
           medens <- try(MoE_dens(data=if(exp.g) e.res else x.dat, mus=mus, sigs=vari, log.tau=ltau, Vinv=Xinv), silent=TRUE)
          }
          if((ERR <- ERR || (attr(Mstep, "returnCode") < 0  || (inherits(medens, "try-error")) || any(medens > 0)))) {
            ll    <- c(ll, NA)
            break
          } else   {
            Estep <- switch(EXPR=alW, EM=MoE_estep(Dens=medens), CEM=MoE_cstep(Dens=medens))
            z     <- zN        <- Estep$z
            ERR   <- any(is.nan(z))
            if(isTRUE(ERR))                       break
            ll    <- c(ll, Estep$loglik)
            j     <- j + 1L
            if(stopaX) {
             ait  <- aitken(ll[seq(j - 2L, j, 1L)])
             dX   <- ifelse(is.numeric(ait$a)  && ait$a < 0, 0L, abs(ait$linf - ll[j - 1L]))
             dX[is.nan(dX)]    <- Inf
            } else     {
             dX   <- abs(ll[j]  - ll[j - 1L])/(1 + abs(ll[j]))
            }
            stX   <- dX >= tol && j  < max.it  && gN    > 1
            if(itwarn && !m0X)  {
             m0W  <- ifelse(!m0X, warnit < j - 2L, m0X)
             if(m0W   && !m0X)  {                 tryCatch(warning("WARNIT", call.=FALSE), warning=function(w)
                                                  message(paste0("\t", algo, " algorithm for the ", modtype, " model has yet to converge in 'warn.it'=", warnit, " iterations\n")))
              m0X <- TRUE
             }
            }
            if(alG == "cemEM"  && !stX  && nW  == 1L)   {
              ll  <- c(ll[j - 1L], ll[j])
              alW <- "EM"
              j   <- nW        <- 2L
              stX <- TRUE
            }
          }
        } # while (j)

      # Store values corresponding to the maximum BIC/ICL/AIC so far
        j2        <- max(1L, j  - switch(EXPR=algo, cemEM=1L, 2L))
        if(isTRUE(verbose))       cat(paste0("\t\t# Iterations: ", ifelse(ERR, "stopped at ", ""), j2, ifelse(last.G && last.T, "\n\n", "\n")))
       #pen.exp   <- ifelse(exp.g, g * d * (fitE$glmnet.fit$df[which(fitE$glmnet.fit$lambda == fitE$lambda.1se)] + 1), exp.pen)
        pen.exp   <- ifelse(exp.g, g * length(stats::coef(fitE)), exp.pen)
        x.df      <- pen.exp + x.df
        max.ll    <- ll[j]
        choose    <- MoE_crit(modelName=modtype, loglik=max.ll, n=n, G=g, z=z, df=x.df)
        bics      <- choose["bic",]
        icls      <- choose["icl",]
        aics      <- choose["aic",]
        crit.t    <- switch(EXPR=criterion, bic=bics, icl=icls, aic=aics)
        crit.t    <- ifelse(is.na(crit.t) || ERR, -Inf, crit.t)
        if(crit.t  > crit.tx)   {
          crit.tx <- crit.t
          tau.x   <- tau
          z.x     <- z
          ll.x    <- ll
          gp.x    <- gate.pen
          ep.x    <- pen.exp
          sig.x   <- vari
          if(stopGx)   {
            linfx <- ait$linf
          }
          if(gate.g)   {
            gfit  <- fitG
          }
          if(exp.g)    {
            efit  <- e.fit
            eres  <- e.res
          } else   {
            mu.x  <- mus
          }
        }
        max.ll                 <- ifelse(ERR, -Inf, max.ll)
        if(!multstart || LL.x[[i]][h,modtype] <     max.ll) {
          BICs[h,modtype]      <- ifelse(ERR, -Inf, bics)
          ICLs[h,modtype]      <- ifelse(ERR, -Inf, icls)
          AICs[h,modtype]      <- ifelse(ERR, -Inf, aics)
          DF.x[h,modtype]      <- ifelse(ERR, -Inf, x.df)
          IT.x[h,modtype]      <- ifelse(ERR,  Inf, j2)
        }
        if(isTRUE(multstart))   {
          LL.x[[i + 1L]][h,modtype]  <- max(max.ll, LL.x[[i]][h,modtype])
        }
      } # for (modtype)
      if((faillen <- length(failedM)) > 1L)       warning(paste0("\tExtra initialisation step with expert covariates worked,\n\t\tbut expert networks themselves couldn't be properly initialised for the G=", g, " ", paste(shQuote(failedM), collapse=" + "), " model", ifelse(faillen == 1, "", "s"), ifelse(drop.exp, ": try setting 'drop.exp' to FALSE\n", ",\n\t\teven with 'drop_constants' and 'drop_levels' invoked:\n\t\tTry suppressing the initialisation step via 'exp.init$mahalanobis' or using other covariates\n")), call.=FALSE, immediate.=TRUE)

    # Pull out mclust model corresponding to highest BIC/ICL/AIC
      if(crit.tx   > crit.gx)   {
        crit.gx   <- crit.tx
        x.tau     <- tau.x
        x.z       <- z.x
        x.ll      <- ll.x
        x.gp      <- gp.x
        x.ep      <- ep.x
        x.sig     <- sig.x
        if(stopGx)     {
          x.linf  <- linfx
        }
        if(gate.g)     {
          x.fitG  <- gfit
        }
        if(exp.g)      {
          x.fitE  <- efit
          x.resE  <- eres
        } else     {
          x.mu    <- mu.x
        }
      }
    } # for (i)
    } # for (g)
    if(any(warnmd <- apply(mderr, 1L, all))) {
      mdwarn      <- paste0("\nInitialisation failed for ", ifelse(all(warnmd), "ALL", "SOME"), " G values due to invocation of 'exp.init$clustMD'")
      if(any(is.element(init.z,
         c("hc", "mclust"))    &&
         inherits(Zhc, "try-error"),
         init.z   == "mclust"  &&
         any(mcfail)))          {
        mdwarn    <- paste0(mdwarn, ":\nback-up option init.z=\"", init.z, "\" also failed")
        if(all(warnmd))         {                 stop(mdwarn, call.=FALSE)
        } else                                    warning(paste0(mdwarn, " for the G=", paste(range.G[if(any(mcfail)) mcfail else warnmd], collapse="/"), " model", ifelse(sum(if(any(mcfail)) mcfail else warnmd) > 1, "s\n", "\n")),       call.=FALSE, immediate.=FALSE)
      } else                                      warning(paste0(mdwarn, ":\ninitialisation defaulted to init.z=\"", init.z, "\" instead for the G=", paste(range.G[warnmd], collapse="/"), " model", ifelse(sum(warnmd) > 1, "s\n", "\n")), call.=FALSE, immediate.=TRUE)
    }
    if(all(is.infinite(BICs[!is.na(BICs)])))      stop("All models failed!", call.=FALSE)

  # Gather results + fit extra gating & expert networks
    CRITs         <- switch(EXPR=criterion, bic=BICs, icl=ICLs, aic=AICs)
    best.ind      <- which(CRITs == crit.gx, arr.ind=TRUE)
    if(nrow(best.ind) > 1)      {                 warning(paste0("Ties for the optimal model exist according to the '", criterion, "' criterion: choosing the most parsimonious model\n"), call.=FALSE, immediate.=TRUE)
      best.ind    <- which(DF.x  == min(DF.x[best.ind]), arr.ind=TRUE)
      best.ind    <- best.ind[which.min(best.ind[,1L]),]
    }
    best.G        <- best.ind[1L]
    G             <- G[best.G]
    if(len.G > 1  && verbose)   {
      if(G        == min(range.G))                message("Best model occurs at the min of the number of components considered\n")
      if(G        == max(range.G))                message("Best model occurs at the max of the number of components considered\n")
    }
    Gseq          <- seq_len(G)
    GN            <- max(G + !noise.null, 1L)
    zN       <- z <- x.z
    rownames(z)   <- as.character(Nseq)
    best.mod      <- colnames(CRITs)[best.ind[2L]]
    bic.fin       <- BICs[best.ind]
    icl.fin       <- ICLs[best.ind]
    aic.fin       <- AICs[best.ind]
    df.fin        <- DF.x[best.ind]
    uncert        <- if(GN > 1) 1     - rowMaxs(z) else vector("integer", n)
    exp.x         <- exp.x & G  != 0
    x.ll          <- x.ll[if(GN == 1 && !exp.x) 2L else if(GN == 1 && exp.x)    2L:3L else switch(EXPR=algo, EM=, CEM=-seq_len(2L), -1L)]
    x.ll          <- x.ll[!is.na(x.ll)]

    n0pro         <- equal.noise[best.G]
    equal.pro     <- equal.tau[best.G]
    noise.gate    <- ifelse(!noise.null, gate.noise[best.G], TRUE)
    if(!(bG       <- gate.G[best.G]))  {
      if(GN > 1)           {
        x.fitG    <- multinom(gating, trace=FALSE, data=gate.covs, maxit=g.itmax, reltol=g.reltol, MaxNWts=MaxNWts)
        if(equal.pro   && !noise.null &&   !n0pro) {
          t0      <- mean(z[,GN])
          x.tau   <- c(rep((1 - t0)/G, G), t0)
        } else     {
          x.tau   <- if(n0pro) rep(1/GN, GN) else if(equal.pro)           rep(1/G, G) else x.fitG$fitted.values[1L,]
        }
        x.tau     <- stats::setNames(x.tau, paste0("Cluster", if(noise.null)     Gseq else c(Gseq, 0L)))
      }   else     {
        x.fitG    <- suppressWarnings(stats::glm(z ~ 1, family=stats::binomial()))
      }
    }
    x.fitG$lab    <- if(noise.gate && !noise.null && GN > 1)              c(Gseq, 0L) else Gseq
    gnames        <- if(G >= 1) paste0("Cluster", Gseq)                               else "Cluster0"
    if(!exp.x)     {
      if(G > 0)    {
        x.fitE    <- list()
        for(k in Gseq)     {
          x.fitE[[k]]          <- stats::lm(expert, weights=z[,k], data=expx.covs)
        }
      } else       {
        x.fitE    <- NA
      }
      residX      <- stats::setNames(replicate(max(G, 1L), matrix(0L, nrow=n, ncol=0L)),   gnames)
    } else residX <- stats::setNames(x.resE, gnames)
    x.fitE        <- stats::setNames(x.fitE, gnames)
    colnames(z)   <- if(G != 0 && !noise.null)                  c(gnames, "Cluster0") else gnames
    exp.gate      <- c(exp.x, bG)
    net.msg       <- ifelse(any(exp.gate), paste0(" (incl. ", ifelse(all(exp.gate), "gating and expert", ifelse(exp.x, "expert", ifelse(bG, "gating", ""))), paste0(" network covariates", ifelse(bG, "", ifelse(equal.pro && G > 1, ", with equal mixing proportions", "")), 
                     ifelse(noise.null, ")", ", and a noise component)"))), ifelse(noise.null, ifelse(equal.pro && G > 1, " (and equal mixing proportions)", ""), ifelse(equal.pro && G > 1, " (with equal mixing proportions and a noise component)", " (and a noise component)")))
    attr(x.fitG, "Maxit")      <- g.itmax
    attr(x.fitG, "Reltol")     <- g.reltol
    attr(x.fitG, "EqualPro")   <- equal.pro
    attr(x.fitG, "EqualNoise") <- n0pro
    attr(x.fitG, "Formula")    <- Reduce(paste, deparse(gating[-2L]))
    attr(x.fitG, "NoiseGate")  <- noise.gate
    attr(x.fitG, "Noise")      <-
    attr(x.fitE, "Noise")      <- G == 0 || !noise.null
    attr(x.fitE, "Formula")    <- Reduce(paste, deparse(expert[-2L]))
    class(x.fitG) <- c("MoE_gating", class(x.fitG))
    class(x.fitE) <- c("MoE_expert", class(x.fitE))
    if(G > 0) {
      vari.fin    <- x.sig
      if(exp.x)    {
        if(noise.null || noise.meth == "manual") {
          z.norm  <- if(noise.null) z else .renorm_z(z[,-GN, drop=FALSE])
          fitdat  <- Reduce("+",  lapply(Gseq, function(g) z.norm[,g] * stats::predict(x.fitE[[g]])))
        } else     {
          fitdat  <- Reduce("+",  lapply(Gseq, function(g) z[,g]      * stats::predict(x.fitE[[g]])))
          fitdat  <- fitdat  +    z[,GN] * matrix(NoiseLoc, nrow=n, ncol=d, byrow=TRUE)
          z.norm  <- if(noise.null) z else z[,-GN, drop=FALSE]
        }
        mean.fin  <- sweep(crossprod(fitdat, z.norm), 2L, colSums2(z.norm), FUN="/", check.margin=FALSE)
      } else       {
        mean.fin  <- x.mu
      }
    } else    {
      mean.fin    <- vari.fin  <- NULL
    }

    if(any(l.warn <- x.ll      != cummax(x.ll)))           {
      if(which.max(l.warn)     != length(x.ll))   warning("Log-likelihoods are not strictly increasing\n", call.=FALSE)
    }
    if(any(IT.x[!is.na(IT.x)]  == max.it))        warning(paste0("One or more models failed to converge in the maximum number of allowed iterations (", max.it, ")\n"), call.=FALSE)
    if(isTRUE(verbose))           cat(paste0("\n\t\tBest Model", ifelse(length(CRITs) > 1, paste0(" (according to ", toupper(criterion), "): "), ": "), mclustModelNames(best.mod)$type, " (", best.mod, "), with ", ifelse(G == 0, "only a noise component",
                                      paste0(G, " component", ifelse(G > 1, "s", ""))), ifelse(any(exp.gate) || (!noise.null && G != 0) || equal.pro, paste0("\n\t\t\t   ", net.msg), ""), "\n\t\t",
                                      switch(EXPR=criterion, bic="BIC", icl="ICL", aic="AIC"), " = ", round(switch(EXPR=criterion, bic=bic.fin, icl=icl.fin, aic=aic.fin), 2), "\n"))
    if(gate.x     && (G + !noise.null - !noise.gate) <= 1) {
      if(attr(x.fitG, "Formula") != "None") tmpnet                    <- netdat
      netdat      <- expx.covs
      attr(netdat, "Gating")   <-
      attr(netdat, "Both")     <- NA
      attr(netdat, "Expert")   <- expx.names
      if(attr(x.fitG, "Formula") != "None") attr(netdat, "Discarded") <- tmpnet
    }
    LL.x          <- AICs/2 + DF.x
    class(BICs)   <- c("MoECriterion", "mclustBIC")
    class(ICLs)   <- c("MoECriterion", "mclustICL")
    class(AICs)   <- c("MoECriterion", "mclustAIC")
    class(DF.x)   <- c("MoECriterion", "mclustDF")
    class(IT.x)   <- c("MoECriterion", "mclustITER")
    class(LL.x)   <- c("MoECriterion", "mclustLoglik")
    attr(BICs, "G")            <-
    attr(ICLs, "G")            <-
    attr(AICs, "G")            <-
    attr(DF.x, "G")            <-
    attr(IT.x, "G")            <-
    attr(LL.x, "G")            <- rownames(BICs)
    attr(BICs, "modelNames")   <-
    attr(ICLs, "modelNames")   <-
    attr(AICs, "modelNames")   <-
    attr(DF.x, "modelNames")   <-
    attr(IT.x, "modelNames")   <-
    attr(LL.x, "modelNames")   <- colnames(BICs)
    if(G     == 0 || !noise.null) {
      hypvol      <- 1/Vinv
      attr(hypvol, "Hypvol0")  <- hypvol
      attr(hypvol, "Location") <- if(noise.meth != "manual") NoiseLoc
      attr(hypvol, "Meth")     <- noise.meth
      attr(hypvol, "Meth0")    <- noise.meth
      attr(hypvol, "g0only")   <- noise.null
      attr(BICs, "Vinv")       <-
      attr(ICLs, "Vinv")       <-
      attr(AICs, "Vinv")       <-
      attr(DF.x, "Vinv")       <-
      attr(IT.x, "Vinv")       <-
      attr(LL.x, "Vinv")       <- Vinv
    } else    {
      hypvol      <- NA
      if(any(range.G == 0))     {
       attr(hypvol, "Hypvol0") <- 1/Vinv
       attr(hypvol, "Meth0")   <- noise.meth
      }
      Vinv        <- NULL
    }
    attr(BICs, "algo")         <-
    attr(ICLs, "algo")         <-
    attr(AICs, "algo")         <-
    attr(DF.x, "algo")         <-
    attr(IT.x, "algo")         <-
    attr(LL.x, "algo")         <- algo
    attr(BICs, "control")      <-
    attr(ICLs, "control")      <-
    attr(AICs, "control")      <-
    attr(DF.x, "control")      <-
    attr(IT.x, "control")      <-
    attr(LL.x, "control")      <- control
    attr(BICs, "warn")         <-
    attr(ICLs, "warn")         <-
    attr(AICs, "warn")         <-
    attr(DF.x, "warn")         <-
    attr(IT.x, "warn")         <-
    attr(LL.x, "warn")         <- isTRUE(verbose)
    attr(BICs, "n")            <-
    attr(ICLs, "n")            <-
    attr(AICs, "n")            <-
    attr(DF.x, "n")            <-
    attr(IT.x, "n")            <-
    attr(LL.x, "n")            <- n
    attr(BICs, "d")            <-
    attr(ICLs, "d")            <-
    attr(AICs, "d")            <-
    attr(DF.x, "d")            <-
    attr(IT.x, "d")            <-
    attr(LL.x, "d")            <- d
    attr(BICs, "oneD")         <-
    attr(ICLs, "oneD")         <-
    attr(AICs, "oneD")         <-
    attr(DF.x, "oneD")         <-
    attr(IT.x, "oneD")         <-
    attr(LL.x, "oneD")         <- uni
    attr(BICs, "criterion")    <- "BIC"
    attr(ICLs, "criterion")    <- "ICL"
    attr(AICs, "criterion")    <- "AIC"
    attr(DF.x, "criterion")    <- "DF"
    attr(IT.x, "criterion")    <- "ITERS"
    attr(LL.x, "criterion")    <- "loglik"
    attr(BICs, "returnCodes")  <- provideDimnames(unname(ifelse(is.na(BICs) | is.infinite(BICs), -1L, 0L)), base=list(rownames(BICs), colnames(BICs)))
    attr(ICLs, "returnCodes")  <- provideDimnames(unname(ifelse(is.na(ICLs) | is.infinite(ICLs), -1L, 0L)), base=list(rownames(ICLs), colnames(ICLs)))
    attr(AICs, "returnCodes")  <- provideDimnames(unname(ifelse(is.na(AICs) | is.infinite(AICs), -1L, 0L)), base=list(rownames(AICs), colnames(AICs)))
    attr(DF.x, "returnCodes")  <- provideDimnames(unname(ifelse(is.na(DF.x) | is.infinite(DF.x), -1L, 0L)), base=list(rownames(DF.x), colnames(DF.x)))
    attr(IT.x, "returnCodes")  <- provideDimnames(unname(ifelse(is.na(IT.x) | is.infinite(IT.x), -1L, 0L)), base=list(rownames(IT.x), colnames(IT.x)))
    attr(LL.x, "returnCodes")  <- provideDimnames(unname(ifelse(is.na(LL.x) | is.infinite(LL.x), -1L, 0L)), base=list(rownames(LL.x), colnames(LL.x)))
    attr(BICs, "initialization")      <-
    attr(ICLs, "initialization")      <-
    attr(AICs, "initialization")      <-
    attr(DF.x, "initialization")      <-
    attr(IT.x, "initialization")      <-
    attr(LL.x, "initialization")      <- list(hcPairs = if((init.z == "hc" && someG && !hcfail) && (!mdind || isTRUE(warnmd[best.G]))) Zhc, subset = NULL, noise = if(!noise.null) noise)
    attr(df.fin, "Gate.Penalty")      <- x.gp
    attr(df.fin, "Expert.Penalty")    <- x.ep
    attr(df.fin, "nVar.Penalty")      <- nVarParams(best.mod, d, G)
    claX       <- max.col(z)
    claX[claX  == G + 1]       <- 0L
    results       <- list(call = call, data = as.data.frame(X), modelName = best.mod,
                          n = n, d = d, G = G, BIC = BICs, ICL = ICLs, AIC = AICs, bic = bic.fin,
                          icl = icl.fin, aic = aic.fin, gating = x.fitG, expert = x.fitE, LOGLIK = LL.x, loglik = x.ll,
                          linf = if(stopaX && G > 1) x.linf else x.ll[length(x.ll)], df = df.fin, iters = IT.x[best.ind],
                          hypvol = hypvol, parameters = list(pro = x.tau, mean = mean.fin, variance = vari.fin, Vinv = Vinv),
                          z = z, classification = stats::setNames(claX, Nseq), uncertainty = stats::setNames(uncert, Nseq),
                          net.covs = netdat, resid.data = residX, DF = DF.x, ITERS = IT.x)
    class(results)             <- "MoEClust"
    attr(results, "Algo")      <- algo
    attr(results, "Criterion") <- criterion
    attr(results, "Details")   <- paste0(best.mod, ": ", ifelse(G == 0, "only a noise component", paste0(G, " component", ifelse(G > 1, "s", ""), net.msg)))
    attr(results, "EqualNoise")<- ctrl$equalNoise
    attr(results, "EqualPro")  <- ctrl$equalPro
    attr(results, "Expert")    <- exp.x
    attr(results, "Gating")    <- bG
    attr(results, "Noise")     <- G == 0 || !noise.null
    attr(results, "NoiseGate") <- ctrl$noise.gate
      return(results)
  }

#' Density for MoEClust Mixture Models
#'
#' Computes densities (or log-densities) of observations in MoEClust mixture models.
#' @param data If there are no expert network covariates, \code{data} should be a numeric matrix or data frame, wherein rows correspond to observations (n) and columns correspond to variables (d). If there are expert network covariates, this should be a list of length G containing matrices/data.frames of (multivariate) WLS residuals for each component.
#' @param mus The mean for each of G components. If there is more than one component, this is a matrix whose k-th column is the mean of the k-th component of the mixture model. For the univariate models, this is a G-vector of means. In the presence of expert network covariates, all values should be equal to \code{0}.
#' @param sigs The \code{variance} component in the parameters list from the output to eg. \code{\link{MoE_clust}}. The components of this list depend on the specification of \code{modelName} (see \code{\link[mclust]{mclustVariance}} for details). The number of components \code{G}, the number of variables \code{d}, and the \code{modelName} are inferred from \code{sigs}.
#' @param log.tau If covariates enter the gating network, an n times G matrix of mixing proportions, otherwise a G-vector of mixing proportions for the components of the mixture. \strong{Must} be on the log-scale in both cases. The default of \code{0} effectively means densities (or log-densities) aren't scaled by the mixing proportions.
#' @param Vinv An estimate of the reciprocal hypervolume of the data region. See the function \code{\link{noise_vol}}. Used only if an initial guess as to which observations are noise is supplied. Mixing proportion(s) must be included for the noise component also.
#' @param logarithm A logical value indicating whether or not the logarithm of the component densities should be returned. This defaults to \code{TRUE}, otherwise component densities are returned, obtained from the component log-densities by exponentiation. The \strong{log}-densities can be passed to \code{\link{MoE_estep}} or \code{\link{MoE_cstep}}.
#'
#' @note This function is intended for joint use with \code{\link{MoE_estep}} or \code{\link{MoE_cstep}}, using the \strong{log}-densities. Note that models with a noise component are facilitated here too.
#' @importFrom mclust "mclustVariance"
#' @importFrom mvnfast "dmvn"
#' @return A numeric matrix whose \code{[i,k]}-th entry is the density or log-density of observation \emph{i} in component \emph{k}, scaled by the mixing proportions. These densities are unnormalised.
#' @keywords clustering
#' @export
#' @author Keefe Murphy - <\email{keefe.murphy@@ucd.ie}>
#'
#' @seealso \code{\link{MoE_estep}}, \code{\link{MoE_cstep}}, \code{\link{MoE_clust}}, \code{\link[mclust]{mclustVariance}}
#' @usage
#' MoE_dens(data,
#'          mus,
#'          sigs,
#'          log.tau = 0L,
#'          Vinv = NULL,
#'          logarithm = TRUE)
#' @examples
#' data(ais)
#' hema  <- ais[,3:7]
#' model <- MoE_clust(hema, G=3, gating= ~ BMI + sex, modelNames="EEE", network.data=ais)
#' Dens  <- MoE_dens(data=hema, mus=model$parameters$mean,
#'                   sigs=model$parameters$variance, log.tau=log(model$parameters$pro))
#'
#' # Construct the z matrix and compute the log-likelihood
#' Estep <- MoE_estep(Dens=Dens)
#' (ll   <- Estep$loglik)
#'
#' # Check that the z matrix & classification are the same as those from the model
#' identical(max.col(Estep$z), as.integer(unname(model$classification))) #TRUE
#' identical(Estep$z, model$z)                                           #TRUE
#'
#' # The same can be done for models with expert covariates &/or a noise component
#' # Note for models with expert covariates that the mean has to be supplied as 0
#' m2    <- MoE_clust(hema, G=2, expert= ~ sex, modelNames="EVE", network.data=ais, tau0=0.1)
#' Dens2 <- MoE_dens(data=m2$resid.data, sigs=m2$parameters$variance, mus=0, 
#'                   log.tau=log(m2$parameters$pro), Vinv=m2$parameters$Vinv)
  MoE_dens        <- function(data, mus, sigs, log.tau = 0L, Vinv = NULL, logarithm = TRUE) {
    ltau.miss     <- missing(log.tau)
    if(any(log.tau > 0))                          stop("'log.tau' cannot be greater than 0: mixing proportions must be supplied on the log scale", call.=FALSE)
    G             <- sigs$G
    Vnul          <- is.null(Vinv)
    Ldat          <- inherits(data, "list")
    if(!Ldat      || (Ldat &&
       length(data)        != max(G, 1L)))       {
      data        <- replicate(G, as.matrix(data), simplify=FALSE)
    } else data   <- lapply(data, as.matrix)
    dat1          <- data[[1L]]
    n             <- ifelse(is.matrix(dat1), nrow(dat1), length(dat1))
    bind2         <- if(n > 1) base::cbind else base::c
    if(G > 0) {
      modelName   <- sigs$modelName
      d           <- sigs$d
      mu.tmp      <- matrix(0L, nrow=d, ncol=G)
      mu.tmp[]    <- mus
      mus         <- mu.tmp
      Gseq        <- seq_len(G)
      sigmas      <- switch(EXPR=modelName, E=, V=, EII=, VII=sqrt(sigs$sigmasq), EEI=sigs$Sigma, EEE=sigs$cholSigma, VVV=sigs$cholsigma, sigs$sigma)
      sq_mat      <- if(d <= 50) sqrt else .sq_mat
      switch(EXPR=modelName, EEV=, EVE=, EVV=, VEE=, VEV=, VVE= {
        densi     <- vapply(Gseq, function(k) dmvn(data[[k]], mus[,k], sigmas[,,k],                 log=TRUE, isChol=FALSE), numeric(n))
      }, EVI=, VEI=, VVI =  {
        densi     <- vapply(Gseq, function(k) dmvn(data[[k]], mus[,k], sq_mat(sigmas[,,k]),         log=TRUE, isChol=TRUE),  numeric(n))
      }, EEE= {
         sigx     <- force_posiDiag(sigmas);
        densi     <- vapply(Gseq, function(k) dmvn(data[[k]], mus[,k], sigx,                        log=TRUE, isChol=TRUE),  numeric(n))
      }, EEI= {
        sigx      <- sq_mat(sigmas);
        densi     <- vapply(Gseq, function(k) dmvn(data[[k]], mus[,k], sigx,                        log=TRUE, isChol=TRUE),  numeric(n))
      }, EII= {
        sigx      <- .mat_sq(sigmas, d);
        densi     <- vapply(Gseq, function(k) dmvn(data[[k]], mus[,k], sigx,                        log=TRUE, isChol=TRUE),  numeric(n))
      }, VII= {
        densi     <- vapply(Gseq, function(k) dmvn(data[[k]], mus[,k], .mat_sq(sigmas[k], d),       log=TRUE, isChol=TRUE),  numeric(n))
      }, VVV= {
        densi     <- vapply(Gseq, function(k) dmvn(data[[k]], mus[,k], force_posiDiag(sigmas[,,k]), log=TRUE, isChol=TRUE),  numeric(n))
      }, E=   {
        densi     <- vapply(Gseq, function(k) stats::dnorm(data[[k]], mus[k], sigmas,    log=TRUE), numeric(n))
      }, V=   {
        densi     <- vapply(Gseq, function(k) stats::dnorm(data[[k]], mus[k], sigmas[k], log=TRUE), numeric(n))
      },                                          stop("Invalid 'modelName'", call.=FALSE))
      test        <- is.infinite(densi) & densi   > 0
      densi[test] <- 0L
    }
    densi         <- if(Vnul) densi else if(G     > 0) bind2(densi, log(Vinv)) else matrix(log(Vinv), nrow=n, ncol=G + !Vnul)
    if(ifelse(is.matrix(log.tau), ncol(log.tau), length(log.tau)) > G + !Vnul)                                  stop(paste0("Too many ", ifelse(is.matrix(log.tau), "columns", "entries"), " in 'log.tau'", ifelse(Vnul, ":\nPerhaps 'Vinv' needs to be supplied?", "")), call.=FALSE)
    densi         <- densi  + if(ltau.miss       || is.matrix(log.tau))                            log.tau else
                              if(length(log.tau) == G + !Vnul) .mat_byrow(log.tau, nrow=n, ncol=G + !Vnul) else stop(paste0("'log.tau' must be given for every component", ifelse(Vnul, "", ", incl. the noise component if 'Vinv' is supplied")), call.=FALSE)
      if(logarithm)  densi    else exp(densi)
  }

#' E-step for MoEClust Models
#'
#' Softmax function to compute the responsibility matrix z and the log-likelihood for MoEClust models, with the aid of \code{\link{MoE_dens}}.
#' @inheritParams MoE_dens
#' @param Dens (Optional) A numeric matrix whose \code{[i,k]}-th entry is the \strong{log}-density of observation \emph{i} in component \emph{k}, scaled by the mixing proportions, to which the softmax function is to be applied, typically obtained by \code{\link{MoE_dens}} but this is not necessary. If this is supplied, all other arguments are ignored, otherwise \code{\link{MoE_dens}} is called according to the other supplied arguments.
#'
#' @return A list containing two elements:
#' \item{\code{z}}{A matrix with \code{n} rows and \code{G} columns containing the probability of cluster membership for each of \code{n} observations and \code{G} clusters.}
#' \item{\code{loglik}}{The estimated log-likelihood, computed efficiently via \code{\link[matrixStats]{rowLogSumExps}}.}
#'
#' @importFrom matrixStats "rowLogSumExps"
#' @importFrom mclust "mclustVariance"
#' @export
#' @note This softmax function is intended for joint use with \code{\link{MoE_dens}}, using the \strong{log}-densities. Caution is advised using this function without explicitly naming the arguments. Models with a noise component are facilitated here too.
#'
#' The E-step can be replaced by a C-step, see \code{\link{MoE_cstep}} and the \code{algo} argument to \code{\link{MoE_control}}.
#' @author Keefe Murphy - <\email{keefe.murphy@@ucd.ie}>
#' @keywords clustering
#' @seealso \code{\link{MoE_dens}}, \code{\link{MoE_clust}}, \code{\link{MoE_cstep}}, \code{\link{MoE_control}}, \code{\link[mclust]{mclustVariance}}, \code{\link[matrixStats]{rowLogSumExps}}
#' @usage
#' MoE_estep(data,
#'           mus,
#'           sigs,
#'           log.tau = 0L,
#'           Vinv = NULL,
#'           Dens = NULL)
#' @examples
#' data(ais)
#' hema   <- ais[,3:7]
#' model  <- MoE_clust(hema, G=3, gating= ~ BMI + sex, modelNames="EEE", network.data=ais)
#' Dens   <- MoE_dens(data=hema, mus=model$parameters$mean,
#'                    sigs=model$parameters$variance, log.tau=log(model$parameters$pro))
#'
#' # Construct the z matrix and compute the log-likelihood
#' Estep  <- MoE_estep(Dens=Dens)
#' (ll    <- Estep$loglik)
#'
#' # Check that the z matrix & classification are the same as those from the model
#' identical(max.col(Estep$z), as.integer(unname(model$classification))) #TRUE
#' identical(Estep$z, model$z)                                           #TRUE
#'
#' # Call MoE_estep directly
#' Estep2 <- MoE_estep(data=hema, sigs=model$parameters$variance,
#'                     mus=model$parameters$mean, log.tau=log(model$parameters$pro))
#' identical(Estep2$loglik, ll)                                          #TRUE
#'
#' # The same can be done for models with expert covariates &/or a noise component
#' # Note for models with expert covariates that the mean has to be supplied as 0
#' m2     <- MoE_clust(hema, G=2, expert= ~ sex, modelNames="EVE", network.data=ais, tau0=0.1)
#' Estep3 <- MoE_estep(data=m2$resid.data, sigs=m2$parameters$variance, mus=0, 
#'                     log.tau=log(m2$parameters$pro), Vinv=m2$parameters$Vinv)
  MoE_estep       <- function(data, mus, sigs, log.tau = 0L, Vinv = NULL, Dens = NULL) {
    if(missing(Dens)) {
      Dens        <- do.call(MoE_dens, as.list(match.call())[-1L])
    } else if(!is.matrix(Dens) ||
              !is.numeric(Dens))                  stop("'Dens' must be a numeric matrix", call.=FALSE)
    norm          <- rowLogSumExps(Dens)
    z             <- provideDimnames(exp(Dens - norm), base=list("", paste0("Cluster", seq_len(ncol(Dens)))))
      return(list(z = z, loglik = sum(norm)))
  }

#' C-step for MoEClust Models
#'
#' Function to compute the assignment matrix z and the conditional log-likelihood for MoEClust models, with the aid of \code{\link{MoE_dens}}.
#' @inheritParams MoE_dens
#' @param Dens (Optional) A numeric matrix whose \code{[i,k]}-th entry is the \strong{log}-density of observation \emph{i} in component \emph{k}, scaled by the mixing proportions, to which the function is to be applied, typically obtained by \code{\link{MoE_dens}} but this is not necessary. If this is supplied, all other arguments are ignored, otherwise \code{\link{MoE_dens}} is called according to the other supplied arguments.
#'
#' @return A list containing two elements:
#' \item{\code{z}}{A matrix with \code{n} rows and \code{G} columns containing 1 where the observation belongs to the cluster indicated by the column number, and 0 otherwise.}
#' \item{\code{loglik}}{The estimated conditional log-likelihood.}
#'
#' @importFrom mclust "unmap"
#' @export
#' @note This function is intended for joint use with \code{\link{MoE_dens}}, using the \strong{log}-densities. Caution is advised using this function without explicitly naming the arguments. Models with a noise component are facilitated here too.
#'
#' The C-step can be replaced by an E-step, see \code{\link{MoE_estep}} and the \code{algo} argument to \code{\link{MoE_control}}.
#' @author Keefe Murphy - <\email{keefe.murphy@@ucd.ie}>
#' @keywords clustering
#' @seealso \code{\link{MoE_dens}}, \code{\link{MoE_clust}}, \code{\link{MoE_estep}}, \code{\link{MoE_control}}, \code{\link[mclust]{mclustVariance}}
#' @usage
#' MoE_cstep(data,
#'           mus,
#'           sigs,
#'           log.tau = 0L,
#'           Vinv = NULL,
#'           Dens = NULL)
#' @examples
#' # MoE_cstep can be invoked for fitting MoEClust models via the CEM algorithm
#' # via the 'algo' argument to MoE_control:
#' data(ais)
#' hema   <- ais[,3:7]
#' model  <- MoE_clust(hema, G=3, gating= ~ BMI + sex, modelNames="EEE", network.data=ais, algo="CEM")
#' Dens   <- MoE_dens(data=hema, mus=model$parameters$mean,
#'                    sigs=model$parameters$variance, log.tau=log(model$parameters$pro))
#'
#' # Construct the z matrix and compute the conditional log-likelihood
#' Cstep  <- MoE_cstep(Dens=Dens)
#' (ll    <- Cstep$loglik)
#'
#' # Check that the z matrix & classification are the same as those from the model
#' identical(max.col(Cstep$z), as.integer(unname(model$classification))) #TRUE
#' identical(Cstep$z, model$z)                                           #TRUE
#'
#' # Call MoE_cstep directly
#' Cstep2 <- MoE_cstep(data=hema, sigs=model$parameters$variance,
#'                     mus=model$parameters$mean, log.tau=log(model$parameters$pro))
#' identical(Cstep2$loglik, ll)                                          #TRUE
  MoE_cstep       <- function(data, mus, sigs, log.tau = 0L, Vinv = NULL, Dens = NULL) {
    if(missing(Dens)) {
      Dens        <- do.call(MoE_dens, as.list(match.call())[-1L])
    } else if(!is.matrix(Dens) ||
              !is.numeric(Dens))                  stop("'Dens' must be a numeric matrix", call.=FALSE)
    Gseq          <- seq_len(ncol(Dens))
    z             <- provideDimnames(unmap(max.col(Dens), groups=Gseq), base=list(as.character(seq_len(nrow(Dens))), paste0("Cluster", Gseq)))
      return(list(z = z, loglik = sum(z * Dens, na.rm=TRUE)))
  }

#' MoEClust BIC, ICL, and AIC Model-Selection Criteria
#'
#' Computes the BIC (Bayesian Information Criterion), ICL (Integrated Complete Likelihood), and AIC (Akaike Information Criterion) for parsimonious mixture of experts models given the log-likelihood, the dimension of the data, the number of mixture components in the model, the numbers of parameters in the gating and expert networks respectively, and, for the ICL, the numbers of observations in each component.
#' @param modelName A character string indicating the model. The help file for \code{\link[mclust]{mclustModelNames}} describes the available models.
#' @param loglik The log-likelihood for a data set with respect to the Gaussian mixture model specified in the \code{modelName} argument.
#' @param n,d,G The number of observations in the data, dimension of the data, and number of components in the Gaussian mixture model, respectively, used to compute \code{loglik}. \code{d} & \code{G} are not necessary if \code{df} is supplied.
#' @param gating.pen The number of parameters of the \emph{gating} network of the MoEClust model. Defaults to \code{G - 1}, which corresponds to no gating covariates. If covariates are included, this should be the number of regression coefficients in the fitted \code{gating} object. If there are no covariates and mixing proportions are further assumed to be present in equal proportion, \code{gating.pen} should be \code{0}. The number of parameters used in the estimation of the noise component, if any, should also be included. Not necessary if \code{df} is supplied.
#' @param expert.pen The number of parameters of the \emph{expert} network of the MoEClust model. Defaults to \code{G * d}, which corresponds to no expert covariates. If covariates are included, this should be the number of regression coefficients in the fitted \code{expert} object. Not necessary if \code{df} is supplied.
#' @param z The \code{n} times \code{G} responsibility matrix whose \code{[i,k]}-th entry is the probability that observation \emph{i} belonds to the \emph{k}-th component.. If supplied the ICL is also computed and returned, otherwise only the BIC and AIC.
#' @param df An alternative way to specify the number of estimated parameters (or 'used' degrees of freedom) exactly. If supplied, the arguments \code{d, G, gating.pen} and \code{expert.pen}, which are used to calculate the number of parameters, will be ignored. The number of parameters used in the estimation of the noise component, if any, should also be included.
#'
#' @details The function is vectorized with respect to the arguments \code{modelName} and \code{loglik}.
#'
#' If \code{model} is an object of class \code{"MoEClust"} with \code{G} components, the number of parameters for the \code{gating.pen} and \code{expert.pen} are \code{length(coef(model$gating))} and \code{G * length(coef(model$expert[[1]]))}, respectively.
#'
#' Models with a noise component are facilitated here too provided the extra number of parameters are accounted for by the user.
#' @importFrom matrixStats "rowMaxs"
#' @importFrom mclust "mclustModelNames" "nVarParams"
#' @return A simplified array containing the BIC, AIC, number of estimated parameters (\code{df}) and, if \code{z} is supplied, also the ICL, for each of the given input arguments.
#' @note In order to speed up repeated calls to the function inside \code{\link{MoE_clust}}, no checks take place.
#' @keywords clustering
#' @export
#' @author Keefe Murphy - <\email{keefe.murphy@@ucd.ie}>
#'
#' @references Biernacki, C., Celeux, G., Govaert, G. (2000). Assessing a mixture model for clustering with the integrated completed likelihood. \emph{IEEE Trans. Pattern Analysis and Machine Intelligence}, 22(7): 719-725.
#'
#' @seealso \code{\link{MoE_clust}}, \code{\link[mclust]{nVarParams}}, \code{\link[mclust]{mclustModelNames}}
#' @usage
#' MoE_crit(modelName,
#'          loglik,
#'          n,
#'          d,
#'          G,
#'          gating.pen = G - 1L,
#'          expert.pen = G * d,
#'          z = NULL,
#'          df = NULL)
#' @examples
#' MoE_crit(modelName=c("VVI", "VVE", "VVV"), n=120, d=8,
#'          G=3, loglik=c(-4036.99, -3987.12, -3992.45))
#'
#' data(CO2data)
#' GNP   <- CO2data$GNP
#' model <- MoE_clust(CO2data$CO2, G=1:2, expert= ~ GNP)
#' G     <- model$G
#' name  <- model$modelName
#' ll    <- max(model$loglik)
#' n     <- length(CO2data$CO2)
#' z     <- model$z
#'
#' # Compare BIC from MoE_crit to the BIC of the model
#' (bic2 <- MoE_crit(modelName=name, loglik=ll, n=n, d=1, G=G, z=z,
#'                   expert.pen=G * length(coef(model$expert[[1]])))["bic",])
#' identical(bic2, unname(model$bic)) #TRUE
#'
#' # Make the same comparison with the known number of estimated parameters
#' (bic3 <- MoE_crit(modelName=name, loglik=ll, n=n, G=G, df=model$df, z=z)["bic",])
#' identical(bic3, bic2)              #TRUE
  MoE_crit        <- Vectorize(function(modelName, loglik, n, d, G, gating.pen = G - 1L, expert.pen = G * d, z = NULL, df = NULL) {
    df            <- ifelse(!missing(df), df, nVarParams(modelName, d, G) + expert.pen + gating.pen)
    double.ll     <- 2 * loglik
    bic.x         <- double.ll  - df * log(n)
    aic.x         <- double.ll  - df * 2
      return(c(bic = bic.x, icl = if(!missing(z)) bic.x + 2L * sum(log(rowMaxs(z)), na.rm=TRUE), aic = aic.x, df = df))
  }, vectorize.args = c("modelName", "loglik"), SIMPLIFY="array")

#' Set control values for use with MoEClust
#'
#' Supplies a list of arguments (with defaults) for use with \code{\link{MoE_clust}}.
#' @param init.z The method used to initialise the cluster labels. Defaults to a model-based agglomerative hierarchical clustering tree as per "\code{\link[mclust]{hc}}" for multivariate data (see \code{hc.args}), or "\code{quantile}"-based clustering as per \code{\link{quant_clust}} for univariate data (unless there are expert network covariates incorporated via \code{exp.init$joint} &/or \code{exp.init$clustMD}, in which case the default is again "\code{\link[mclust]{hc}}"). The \code{"quantile"} option is thus only available for univariate data when expert network covariates are not incorporated via \code{exp.init$joint} &/or \code{exp.init$clustMD}, or when expert network covariates are not supplied.
#'
#' Other options include "\code{kmeans}" (see \code{km.args}), "\code{random}" initialisation, and a full run of \code{\link[mclust]{Mclust}} (itself initialised via a model-based agglomerative hierarchical clustering tree, again see \code{hc.args}), although this last option "\code{mclust}" will be coerced to "\code{hc}" if there are no \code{gating} &/or \code{expert} covariates within \code{\link{MoE_clust}} (in order to better reproduce \code{\link[mclust]{Mclust}} output).
#'
#' When \code{isTRUE(exp.init$clustMD)} and the \code{\link[clustMD]{clustMD}} library is loaded, the \code{init.z} argument instead governs the method by which a call to \code{\link[clustMD]{clustMD}} is initialised. In this instance, "\code{quantile}" will instead default to "\code{\link[mclust]{hc}}", and the arguments to \code{hc.args} and \code{km.args} will be ignored (unless all \code{\link[clustMD]{clustMD}} model types fail for a given number of components).
#'
#' When \code{init.z="mclust"} or \code{\link[clustMD]{clustMD}} is successfully invoked (via \code{exp.init$clustMD}), the argument \code{init.crit} (see below) specifies the model-selection criterion ("\code{bic}" or "\code{icl}") by which the optimal \code{\link[mclust]{Mclust}} or \code{\link[clustMD]{clustMD}} model type to initialise with is determined, and \code{criterion} remains unaffected.
#' @param noise.args A list supplying select named parameters to control inclusion of a noise component in the estimation of the mixture. If either or both of the arguments \code{tau0} &/or \code{noise.init} are supplied, a noise component is added to the the model in the estimation.
#' \describe{
#' \item{\code{tau0}}{Prior mixing proportion for the noise component. If supplied, a noise component will be added to the model in the estimation, with \code{tau0} giving the prior probability of belonging to the noise component for \emph{all} observations. Typically supplied as a scalar in the interval (0, 1), e.g. \code{0.1}. Can be supplied as a vector when gating covariates are present and \code{noise.args$noise.gate} is \code{TRUE}. This argument can be spplied instead of or in conjuction with the argument \code{noise.init} below.}
#' \item{\code{noise.init}}{A logical or numeric vector indicating an initial guess as to which observations are noise in the data. If numeric, the entries should correspond to row indices of the data. If supplied, a noise component will be added to the model in the estimation. This argument can be used in conjunction with \code{tau0} above, or can be replaced by that argument also.}
#' \item{\code{noise.gate}}{A logical indicating whether gating network covariates influence the mixing proportion for the noise component, if any. Defaults to \code{TRUE}, but leads to greater parsimony if \code{FALSE}. Only relevant in the presence of a noise component; only effects estimation in the presence of gating covariates.}
#' \item{\code{noise.meth}}{The method used to estimate the volume when a noise component is invoked. Defaults to \code{\link[mclust]{hypvol}}. For univariate data, this argument is ignored and the range of the data is used instead (unless \code{noise.vol} below is specified). The options "\code{convexhull}" and "\code{ellipsoidhull}" require loading the \code{geometry} and \code{cluster} libraries, respectively. This argument is only relevant if \code{noise.vol} below is not supplied.}
#' \item{\code{noise.vol}}{This argument can be used to override the argument \code{noise.meth} by specifying the (hyper)volume directly, i.e. specifying an improper uniform density. This will override the use of the range of the response data for univariate data if supplied. Note that the (hyper)volume, rather than its inverse, is supplied here. This can affect prediction and the location of the MVN ellipses for \code{\link{MoE_gpairs}} plots (see \code{\link{noise_vol}}).}
#' \item{\code{equalNoise}}{Logical which is only invoked when \code{isTRUE(equalPro)} and gating covariates are not supplied. Under the default setting (\code{FALSE}), the mixing proportion for the noise component is estimated, and remaining mixing proportions are equal; when \code{TRUE} all components, including the noise component, have equal mixing proportions.}
#' }
#' In particular, the argument \code{noise.meth} will be ignored for high-dimensional \code{n <= d} data, in which case the argument \code{noise.vol} \emph{must be} specified.
#' 
#' The arguments \code{tau0} and \code{noise.init} can be used separately, to provide alternative means to invoke a noise component. However, they can also be supplied together, in which case observations corresponding to \code{noise.init} have probability \code{tau0} (rather than 1) of belonging to the noise component.
#' @param equalPro Logical variable indicating whether or not the mixing proportions are to be constrained to be equal in the model. Default: \code{equalPro = FALSE}. Only relevant when \code{gating} covariates are \emph{not} supplied within \code{\link{MoE_clust}}, otherwise ignored. In the presence of a noise component (see \code{noise.args}), only the mixing proportions for the non-noise components are constrained to be equal, after accounting for the noise component.
#' @param exp.init A list supplying select named parameters to control the initialisation routine in the presence of \emph{expert} network covariates (otherwise ignored):
#' \describe{
#' \item{\code{joint}}{A logical indicating whether the initial partition is obtained on the joint distribution of the response and expert network covariates (defaults to \code{TRUE}) or just the response variables (\code{FALSE}). By default, only continuous expert network covariates are considered (see \code{exp.init$clustMD} below). Only relevant when \code{init.z} is not \code{"random"} (unless \code{isTRUE(exp.init$clustMD)}, in which case \code{init.z} specifies the initialisation routine for a call to \code{\link[clustMD]{clustMD}}). This will render the \code{"quantile"} option to \code{init.z} for univariate data unusable if continuous expert network covariates are supplied &/or categorical/ordinal expert network covariates are supplied when \code{isTRUE(exp.init$clustMD)} and the \code{\link[clustMD]{clustMD}} library is loaded.}
#' \item{\code{mahalanobis}}{A logical indicating whether to iteratively reallocate observations during the initialisation phase to the component corresponding to the expert network regression to which it's closest to the fitted values of in terms of Mahalanobis distance (defaults to \code{TRUE}). This will ensure that each component can be well modelled by a single expert prior to running the EM/CEM algorithm.}
#' \item{\code{\link[clustMD]{clustMD}}}{A logical indicating whether categorical/ordinal covariates should be incorporated when using the joint distribution of the response and expert network covariates for initialisation (defaults to \code{FALSE}). Only relevant when \code{isTRUE(exp.init$joint)}. Requires the use of the \code{\link[clustMD]{clustMD}} library. Note that initialising in this manner involves fiting all \code{\link[clustMD]{clustMD}} model types in parallel for all numbers of components considered, and may fail (especially) in the presence of nominal expert network covariates.
#'
#' Supplying this argument as \code{TRUE} when the \code{\link[clustMD]{clustMD}} library is loaded has the effect of superseding the \code{init.z} argument: this argument now governs instead how the call to \code{\link[clustMD]{clustMD}} is initialised (unless all \code{\link[clustMD]{clustMD}} model types fail for a given number of components, in which case \code{init.z} is invoked \emph{instead} to initialise for \code{G} values for which all \code{\link[clustMD]{clustMD}} model types failed). Similarly, the arguments \code{hc.args} and \code{km.args} will be ignored (again, unless all \code{\link[clustMD]{clustMD}} model types fail for a given number of components).}
#' \item{\code{max.init}}{The maximum number of iterations for the Mahalanobis distance-based reallocation procedure when \code{exp.init$mahalanobis} is \code{TRUE}. Defaults to \code{100}.}
#' \item{\code{drop.break}}{When \code{isTRUE(exp.init$mahalanobis)} observations will be completely in or out of a component during the initialisation phase. As such, it may occur that constant columns will be present when building a given componenet's expert regression (particularly for categorical covariates). It may also occur, due to this partitioning, that "unseen" data, when calculating the residuals, will have new factor levels. When \code{isTRUE(exp.init$drop.break)}, the Mahalanobis distance based initialisation phase will explicitly fail in either of these scenarios.
#'
#' Otherwise, \code{\link{drop_constants}} and \code{\link{drop_levels}} will be invoked when \code{exp.init$drop.break} is \code{FALSE} (the default) to \emph{try} to remedy the situation. In any case, only a warning that the initialisation step failed will be printed, regardless of the value of \code{exp.init$drop.break}.}
#' }
#' @param algo Switch controlling whether models are fit using the \code{"EM"} (the default) or \code{"CEM"} algorithm. The option \code{"cemEM"} allows running the EM algorithm starting from convergence of the CEM algorithm.
#' @param criterion When either \code{G} or \code{modelNames} is a vector, \code{criterion} determines whether the "\code{bic}" (Bayesian Information Criterion), "\code{icl}" (Integrated Complete Likelihood), "\code{aic}" (Akaike Information Criterion) is used to determine the 'best' model when gathering output. Note that all criteria will be returned in any case.
#' @param stopping The criterion used to assess convergence of the EM/CEM algorithm. The default (\code{"aitken"}) uses Aitken's acceleration method via \code{\link{aitken}}, otherwise the \code{"relative"} change in log-likelihood is monitored (which may be less strict). Both stopping rules are ultimately governed by \code{tol[1]}. When the \code{"aitken"} method is employed, the asymptotic estimate of the final converged maximised log-likelihood is also returned as \code{linf} for models with 2 or more components, though the largest element of the returned vector \code{loglik} still gives the log-likelihood value achieved by the parameters returned at convergence, under both \code{stopping} methods (see \code{\link{MoE_clust}}).
#' @param nstarts The number of random initialisations to use when \code{init.z="random"}. Defaults to 1. Results will be based on the random start yielding the highest estimated log-likelihood. Note that all \code{nstarts} random initialisations are affected by \code{exp.init$mahalanobis}, if invoked in the presence of expert network covariates, which may remove some of the randomness.
#' @param eps A scalar tolerance associated with deciding when to terminate computations due to computational singularity in covariances. Smaller values of \code{eps} allow computations to proceed nearer to singularity. The default is the relative machine precision \code{.Machine$double.eps}, which is approximately \emph{2e-16} on IEEE-compliant machines.
#' @param tol A vector of length three giving relative convergence tolerances for 1) the log-likelihood of the EM/CEM algorithm, 2) parameter convergence in the inner loop for models with iterative M-step (\code{"VEI", "EVE", "VEE", "VVE", "VEV"}), and 3) optimisation in the multinomial logistic regression in the gating network, respectively. The default is \code{c(1e-05, sqrt(.Machine$double.eps), 1e-08)}. If only one number is supplied, it is used as the tolerance for all three cases given.
#' @param itmax A vector of length three giving integer limits on the number of iterations for 1) the EM/CEM algorithm, 2) the inner loop for models with iterative M-step (\code{"VEI", "EVE", "VEE", "VVE", "VEV"}), and 3) the multinomial logistic regression in the gating network, respectively.
#'
#' The default is \code{c(.Machine$integer.max, .Machine$integer.max, 100)} allowing termination to be completely governed by \code{tol} for the inner and outer loops of the EM. If only one number is supplied, it is used as the iteration limit for the outer loop only.
#' @param hc.args A list supplying select named parameters to control the initialisation of the cluster allocations when \code{init.z="hc"} (or when \code{init.z="mclust"}, which itself relies on \code{\link[mclust]{hc}}), unless \code{isTRUE(exp.init$clustMD)}, the \code{\link[clustMD]{clustMD}} library is loaded, and none of the \code{\link[clustMD]{clustMD}} model types fail (otherwise irrelevant):
#' \describe{
#' \item{\code{hcUse}}{A string specifying the type of input variables to be used. Unlike \code{\link[mclust]{Mclust}}, this defaults to "\code{VARS}" here.}
#' \item{\code{hc.meth}}{A character string indicating the model to be used when hierarchical clustering (see \code{\link[mclust]{hc}}) is employed for initialisation (either when \code{init.z="hc"} or \code{init.z="mclust"}). Defaults to \code{"EII"} for high-dimensional data, or \code{"VVV"} otherwise.}
#' }
#' @param km.args A list supplying select named parameters to control the initialisation of the cluster allocations when \code{init.z="kmeans"}, unless \code{isTRUE(exp.init$clustMD)}, the \code{\link[clustMD]{clustMD}} library is loaded, and none of the \code{\link[clustMD]{clustMD}} model types fail (otherwise irrelevant):
#' \describe{
#' \item{\code{kstarts}}{The number of random initialisations to use. Defaults to 10.}
#' \item{\code{kiters}}{The maximum number of K-Means iterations allowed. Defaults to 10.}
#' }
#' @param init.crit The criterion to be used to determine the optimal model type to initialise with, when \code{init.z="mclust"} or when \code{isTRUE(exp.init$clustMD)} and the \code{\link[clustMD]{clustMD}} library is loaded (one of "\code{bic}" or "\code{icl}"). Defaults to "\code{icl}" when \code{criterion="icl"}, otherwise defaults to "\code{bic}". The \code{criterion} argument remains unaffected.
#' @param warn.it A single number giving the iteration count at which a warning will be printed if the EM/CEM algorithm has failed to converge. Defaults to \code{0}, i.e. no warning (which is true for any \code{warn.it} value less than \code{3}), otherwise the message is printed regardless of the value of \code{verbose}. If non-zero, \code{warn.it} should be moderately large, but obviously less than \code{itmax[1]}. A warning will always be printed if one of more models fail to converge in \code{itmax[1]} iterations.
#' @param MaxNWts The maximum allowable number of weights in the call to \code{\link[nnet]{multinom}} for the multinomial logistic regression in the gating network. There is no instrinsic limit in the code, but increasing \code{MaxNWts} will probably allow fits that are very slow and time-consuming. It may be necessary to increase \code{MaxNWts} when categorical concomitant variables with many levels are included or the number of components is high.
#' @param verbose Logical indicating whether to print messages pertaining to progress to the screen during fitting. By default is \code{TRUE} if the session is interactive, and \code{FALSE} otherwise. If \code{FALSE}, warnings and error messages will still be printed to the screen, but everything else will be suppressed.
#' @param ... Catches unused arguments.
#'
#' @details \code{\link{MoE_control}} is provided for assigning values and defaults within \code{\link{MoE_clust}}.
#'
#' While the \code{criterion} argument controls the choice of the optimal number of components and GPCM/\pkg{mclust} model type, \code{\link{MoE_compare}} is provided for choosing between fits with different combinations of covariates or different initialisation settings.
#' @importFrom mclust "hc" "hypvol" "mclust.options"
#' @importFrom nnet "multinom"
#' @note Note that successfully invoking \code{exp.init$clustMD} (though it defaults to \code{FALSE}) effects the role of the arguments \code{init.z}, \code{hc.args}, and \code{km.args}. Please read the documentation above carefully in this instance.
#' @return A named list in which the names are the names of the arguments and the values are the values supplied to the arguments.
#' @export
#' @keywords control
#' @author Keefe Murphy - <\email{keefe.murphy@@ucd.ie}>
#'
#' @seealso \code{\link{MoE_clust}}, \code{\link{aitken}}, \code{\link[mclust]{hc}}, \code{\link[mclust]{mclust.options}}, \code{\link{quant_clust}}, \code{\link[clustMD]{clustMD}}, \code{\link[mclust]{hypvol}}, \code{\link[geometry]{convhulln}}, \code{\link[cluster]{ellipsoidhull}}, \code{\link{MoE_compare}}, \code{\link[nnet]{multinom}}
#' @usage
#' MoE_control(init.z = c("hc", "quantile", "kmeans", "mclust", "random"),
#'             noise.args = list(...),
#'             equalPro = FALSE,
#'             exp.init = list(...),
#'             algo = c("EM", "CEM", "cemEM"),
#'             criterion = c("bic", "icl", "aic"),
#'             stopping = c("aitken", "relative"),
#'             nstarts = 1L,
#'             eps = .Machine$double.eps,
#'             tol = c(1e-05, sqrt(.Machine$double.eps), 1e-08),
#'             itmax = c(.Machine$integer.max, .Machine$integer.max, 100L),
#'             hc.args = list(...),
#'             km.args = list(...),
#'             init.crit = c("bic", "icl"),
#'             warn.it = 0L,
#'             MaxNWts = 1000L,
#'             verbose = interactive(),
#'             ...)
#' @examples
#' \dontrun{
#' ctrl1 <- MoE_control(criterion="icl", itmax=100, warn.it=15, init.z="random", nstarts=5)
#'
#' data(CO2data)
#' GNP   <- CO2data$GNP
#' res   <- MoE_clust(CO2data$CO2, G=2, expert = ~ GNP, control=ctrl1)
#'
#' # Alternatively, specify control arguments directly
#' res2  <- MoE_clust(CO2data$CO2, G=2, expert = ~ GNP, stopping="relative")
#'
#' # Supplying ctrl1 without naming it as control throws an error,
#' # when any of {modelNames, gating, expert} are not supplied
#' res3  <- MoE_clust(CO2data$CO2, G=2, expert = ~ GNP, ctrl1)
#'
#' # Initialise via the mixed-type joint distribution of response & covariates
#' # Let the ICL criterion determine the optimal clustMD model type
#' # Constrain the mixing proportions to be equal
#' ctrl2 <- MoE_control(exp.init=list(clustMD=TRUE), init.crit="icl", equalPro=TRUE)
#' data(ais)
#' library(clustMD)
#' res4  <- MoE_clust(ais[,3:7], G=2, modelNames="EVE", expert=~sex,
#'                    network.data=ais, control=ctrl2)
#'
#' # Include a noise component by specifying its prior mixing proportion
#' res5  <- MoE_clust(ais[,3:7], G=2, modelNames="EVE", expert=~sex,
#'                    network.data=ais, tau0=0.1)}
  MoE_control     <- function(init.z = c("hc", "quantile", "kmeans", "mclust", "random"), noise.args = list(...), equalPro = FALSE, exp.init = list(...), algo = c("EM", "CEM", "cemEM"), 
                              criterion = c("bic", "icl", "aic"), stopping = c("aitken", "relative"), nstarts = 1L, eps = .Machine$double.eps, tol = c(1e-05, sqrt(.Machine$double.eps), 1e-08), 
                              itmax = c(.Machine$integer.max, .Machine$integer.max, 100L), hc.args = list(...), km.args = list(...), init.crit = c("bic", "icl"), warn.it = 0L, MaxNWts = 1000L, verbose = interactive(), ...) {
    if(!missing(algo)       && length(algo) > 1 ||
       !is.character(algo))                       stop("'algo' must be a single character string",      call.=FALSE)
    algo          <- match.arg(algo)
    if(!missing(criterion)  && (length(criterion) > 1  ||
       !is.character(criterion)))                 stop("'criterion' must be a single character string", call.=FALSE)
    criterion     <- match.arg(criterion)
    if(!missing(stopping)   && (length(stopping)  > 1  ||
       !is.character(stopping)))                  stop("'stopping' must be a single character string",  call.=FALSE)
    stopping      <- match.arg(stopping)
    miss.init     <- missing(init.z)
    if(!miss.init && (length(init.z)        > 1 ||
       !is.character(init.z)))                    stop("'init.z' must be a single character string",    call.=FALSE)
    init.z        <- match.arg(init.z)
    if(init.z     == "random")      {
      if(length(nstarts)    != 1   ||
         !is.numeric(nstarts)      ||
         (nstarts            < 1   ||
          floor(nstarts)    != nstarts))          stop(paste0("'nstarts' must be a single integer >= 1 if when 'init.z'=", init.z), call.=FALSE)
    }

    if(is.null(exp.init$joint))     {
      exp.init$joint        <- TRUE
    } else if(length(exp.init$joint)        > 1 ||
              !is.logical(exp.init$joint))        stop("'exp.init$joint' must be a single logical indicator",           call.=FALSE)
    if(is.null(exp.init$clustMD))   {
      exp.init$clustMD      <- FALSE
    } else if(length(exp.init$clustMD)      > 1 ||
              !is.logical(exp.init$clustMD))      stop("'exp.init$clustMD' must be a single logical indicator",         call.=FALSE)
    if(is.null(exp.init$mahalanobis))       {
      exp.init$mahalanobis  <- TRUE
    } else if(length(exp.init$mahalanobis)  > 1 ||
              !is.logical(exp.init$mahalanobis))  stop("'exp.init$mahalanobis' must be a single logical indicator",     call.=FALSE)
    if(is.null(exp.init$max.init))          {
      exp.init$max.init     <- 100L
    } else if(isTRUE(exp.init$mahalanobis)      &&
             (length(exp.init$max.init)     > 1 ||
             ((!is.numeric(exp.init$max.init)   ||
              exp.init$max.init    <= 0))))       stop("'exp.init$max.init' must be a single strictly positive integer when 'exp.init$mahalanobis' is TRUE", call.=FALSE)
    if(is.null(exp.init$drop.break))        {
      exp.init$drop.break   <- FALSE
    } else if(isTRUE(exp.init$mahalanobis)      &&
              (length(exp.init$drop.break)  > 1 ||
              !is.logical(exp.init$drop.break)))  stop("'exp.init$break' must be a single logical indicator when 'exp.init$mahalanobis' is TRUE",            call.=FALSE)

    if(length(eps) > 2)                           stop("'eps' can be of length at most 2",   call.=FALSE)
    if(any(eps     < 0))                          stop("'eps' is negative",                  call.=FALSE)
    if(any(eps    >= 1))                          stop("'eps' is not less than 1",           call.=FALSE)
    if((len.tol   <- length(tol))   > 3)          stop("'tol' can be of length at most 3",   call.=FALSE)
    if(!is.numeric(tol))                          stop("'tol' must be numeric",              call.=FALSE)
    if(any(tol     < 0))                          stop("'tol' is negative",                  call.=FALSE)
    if(any(tol    >= 1))                          stop("'tol' is not less than 1",           call.=FALSE)
    if(len.tol    == 1)        tol <- rep(tol, 3L)
    if(!is.numeric(itmax)   ||
       any(itmax  != floor(itmax)))               stop("'itmax' must be of integer type",    call.=FALSE)
    if(any(itmax  <= 0))                          stop("'itmax' is not strictly positive",   call.=FALSE)
    if((len.itmax <- length(itmax)) > 3)          stop("'itmax' can be of length at most 3", call.=FALSE)
    if(len.itmax  == 1)      itmax <- c(itmax, .Machine$integer.max, 100L)
    inf           <- is.infinite(itmax)
    if(any(inf))        itmax[inf] <- .Machine$integer.max
    if(length(MaxNWts)  > 1 ||
       !is.numeric(MaxNWts) ||
       MaxNWts    <= 0)                           stop("'MaxNWts' must be a strictly positive scalar", call.=FALSE)
    if(length(equalPro) > 1 ||
       !is.logical(equalPro))                     stop("'equalPro' must be a single logical indicator", call.=FALSE)

    if(!is.null(noise.args$noise.init)          ||
       !is.null(noise.args$tau0))                {
      if(!is.null(noise.args$tau0)              &&
        (!is.numeric(noise.args$tau0)           ||
        any(noise.args$tau0 <= 0)  ||
        any(noise.args$tau0 >= 1)))               stop("'noise.args$tau0' must lie in the interval (0, 1)",  call.=FALSE)
      if(!is.null(noise.args$noise.vol))         {
        if(length(noise.args$noise.vol) > 1     ||
           noise.args$noise.vol    <= 0)          stop("Invalid 'noise.args$noise.vol'", call.=FALSE)
        noise.args$noise.meth      <- "manual"
      } else if(is.null(noise.args$noise.meth))  {
        noise.args$noise.meth      <- "hypvol"
      } else {
        if(length(noise.args$noise.meth)         > 1   ||
           !is.character(noise.args$noise.meth))  stop("'noise.args$noise.meth' must be a single character string",     call.=FALSE)
        if(!is.element(noise.args$noise.meth, c("hypvol",
                "ellipsoidhull", "convexhull")))  stop("'noise.args$noise.meth' must be one of 'hypvol', 'ellipsoidhull', or 'convexhull'",                  call.=FALSE)
      }
      if(is.null(noise.args$noise.gate))         {
        noise.args$noise.gate      <- TRUE
      } else if(length(noise.args$noise.gate)    > 1   ||
        !is.logical(noise.args$noise.gate))       stop("'noise.args$noise.gate' must be a single logical indicator",    call.=FALSE)
      if(is.null(noise.args$equalNoise))         {
        noise.args$equalNoise      <- FALSE
      } else if(length(noise.args$equalNoise)    > 1   ||
        !is.logical(noise.args$equalNoise))       stop("noise.args$equalNoise' must be a single logical indicator",     call.=FALSE)
      has.lib     <- switch(EXPR=noise.args$noise.meth, manual=, hypvol=TRUE, convexhull=suppressMessages(requireNamespace("geometry", quietly=TRUE)), ellipsoidhull=suppressMessages(requireNamespace("cluster", quietly=TRUE)))
      if(!has.lib)                                stop(paste0("Use of the ", noise.args$noise.meth, " option for 'noise.args$noise.meth' requires loading the ",
                                                              switch(EXPR=noise.args$noise.meth, hypvol="'mclust'", convexhull="'geometry'",
                                                              ellipsoidhull="'cluster'"), "library"),   call.=FALSE)
    }

    if(!(miss.hc  <- is.null(hc.args$hc.meth)))  {
      if(init.z   == "hc"   &&
         !is.element(hc.args$hc.meth, c("E", "V",
         mclust.options("hcModelNames"))))        stop("Invalid 'hc.args$hc.meth' selected for initialisation by agglomerative hierarchical clustering",     call.=FALSE)
    }
    if(is.null(hc.args$hcUse))      {
      hc.args$hcUse         <- "VARS"
    } else {
      hcUse                 <- hc.args$hcUse
      if(length(hcUse)  > 1 || (!is.character(hcUse)   ||
       !is.element(hcUse, c("VARS", "STD", "SPH",
                          "PCS", "PCR", "SVD")))) stop("Invalid 'hc.args$hcUse'",            call.=FALSE)
      hc.args$hcUse         <- hcUse
    }

    if(is.null(km.args$kiters))     {
      km.args$kiters        <- 10L
    } else {
      kiters      <- km.args$kiters
      if(length(kiters)     != 1   ||
       !is.numeric(kiters)  ||
       kiters      < 1      ||
       kiters     != floor(kiters))               stop("'km.args$kiters' must be a single strictly positive integer",   call.=FALSE)
      km.args$kiters        <- kiters
    }
    if(is.null(km.args$kstarts))    {
      km.args$kstarts       <- 10L
    } else {
      kstarts     <- km.args$kstarts
      if(length(kstarts)    != 1   ||
       !is.numeric(kstarts) ||
       kstarts     < 1      ||
       kstarts    != floor(kstarts))              stop("'km.argss$kstarts' must be a single strictly positive integer", call.=FALSE)
      km.args$kstarts       <- kstarts
    }

    if(!missing(init.crit)  && (length(init.crit) > 1   ||
       !is.character(init.crit)))                 stop("'init.crit' must be a single character string", call.=FALSE)
    init.crit     <- ifelse(missing(init.crit),   switch(EXPR=criterion, icl="icl", "bic"), init.crit)
    init.crit     <- match.arg(init.crit)
    if(length(warn.it)  > 1 ||
       !is.numeric(warn.it) ||
       warn.it     < 0      ||
       warn.it    != floor(warn.it))              stop("'warn.it' must be a single strictly non-negative integer",      call.=FALSE)
    if(length(verbose)  < 1 ||
       !is.logical(verbose))                      stop("'verbose' must be a single logical indicator",  call.=FALSE)
      list(algo = algo, criterion = criterion, stopping = stopping, init.z = init.z, nstarts = nstarts, exp.init = exp.init, eps = eps, tol = tol, itmax = itmax, MaxNWts = MaxNWts, 
           equalPro = equalPro, noise.args = noise.args, hc.args = hc.args, km.args = km.args, init.crit = init.crit, warn.it = warn.it, verbose = verbose, miss.init = miss.init, miss.hc = miss.hc)
  }

#' Aitken Acceleration
#'
#' Calculates the Aitken acceleration estimate of the final converged maximised log-likelihood under the EM/CEM framework.
#' @param loglik A vector of three consecutive log-likelihood values. These three values should be in ascending order, though this is not checked.
#'
#' @details The final converged maximised log-likelihood can be used to determine convergence of the EM/CEM algorithm within \code{\link{MoE_clust}}, i.e. by checking whether the absolute difference between the current log-likelihood estimate and the final converged maximised log-likelihood estimate is less than some tolerance.
#' @note Within \code{\link{MoE_clust}}, as specified by the \code{stopping} argument of \code{\link{MoE_control}}, \code{"aitken"} is the default method used to assess convergence. The other option monitors the \code{"relative"} change in log-likelihood against some tolerance. See \code{\link{MoE_control}}.
#'
#' @return A list with the following named components:
#' \item{\code{ll}}{The most current estimate for the log-likelihood.}
#' \item{\code{linf}}{The most current estimate of the final converged maxmised log-likelihood.}
#' \item{\code{a}}{The Aitken acceleration value where typically \code{0 <= a <= 1}. When \code{a < 0}, a numerical issue or bug has occured; when \code{a > 1}, the algorithm is accelerating and should not be stopped.}
#' When the \code{"aitken"} method is employed within \code{\link{MoE_clust}} (via \code{\link{MoE_control}}), \code{ll} at convergence gives the log-likelihood achieved by the estimated parameters, while \code{linf} at convergence estimates the log-likelihood that would be achieved after an infinite number of EM/CEM iterations.
#' @export
#' @keywords control
#' @author Keefe Murphy - <\email{keefe.murphy@@ucd.ie}>
#' @references Boehning, D., Dietz, E., Schaub, R., Schlattmann, P. and Lindsay, B. G. (1994). The distribution of the likelihood ratio for mixtures of densities from the one-parameter exponential family. \emph{Annals of the Institute of Statistical Mathematics}, 46(2): 373-388.
#'
#' @seealso \code{\link{MoE_control}}
#' @examples
#' (a1 <- aitken(-c(449.61534, 442.84221, 436.58999)))
#' a2  <- aitken(-c(442.84221, 436.58999, 436.58998))
#' abs(a2$linf - a1$linf) < 1e-05 #FALSE
#' a3  <- aitken(-c(436.58998, 436.58997, 436.58997))
#' abs(a3$linf - a2$linf) < 1e-05 #TRUE
#' (ll <- a3$linf)
#' (a  <- a3$a)
  aitken          <- function(loglik) {
    if(!is.numeric(loglik) ||
       length(loglik)      != 3)                  stop("'loglik' must be a numeric vector of length 3", call.=FALSE)
    l1            <- loglik[1L]
    l2            <- loglik[2L]
    l3            <- loglik[3L]
    if(any(is.infinite(loglik))) {
      linf        <- Inf
      a           <- NA
    } else {
      a           <- ifelse(l2 > l1, (l3 - l2) / (l2 - l1),    0L)
      denom       <- max(1L - a, .Machine$double.eps)
      linf        <- ifelse(a  < 1L,  l2 + (l3 - l2) / denom, Inf)
    }
      return(list(ll = l3, linf = linf, a = a))
  }

#' Choose the best MoEClust model
#'
#' Takes one or more sets of MoEClust models fitted by \code{\link{MoE_clust}} and ranks them according to the BIC, ICL, or AIC. It's possible to respect the internal ranking within each set of models, or to discard models within each set which were already deemed sub-optimal. This function can help with model selection via exhaustive or stepwise searches.
#' @param ... One or more objects of class \code{"MoEClust"} outputted by \code{\link{MoE_clust}}. All models must have been fit to the same data set. A single \emph{named} list of such objects can also be supplied. This argument is only relevant for the \code{\link{MoE_compare}} function and will be ignored for the associated \code{print} function.
#' @param criterion The criterion used to determine the ranking. Defaults to \code{"bic"}.
#' @param pick The (integer) number of models to be ranked and compared. Defaults to \code{10L}. Will be constrained by the number of models within the \code{"MoEClust"} objects supplied via \code{...} if \code{optimal.only} is \code{FALSE}, otherwise constrained simply by the number of \code{"MoEClust"} objects supplied. Setting \code{pick=Inf} is a valid way to select all models.
#' @param optimal.only Logical indicating whether to only rank models already deemed optimal within each \code{"MoEClust"} object (\code{TRUE}), or to allow models which were deemed suboptimal enter the final ranking (\code{FALSE}, the default). See \code{details}
#' @param x,index,digits Arguments required for the associated \code{print} function:
#' \describe{
#' \item{\code{x}}{An object of class \code{"MoECompare"} resulting from a call to \code{\link{MoE_compare}}.}
#' \item{\code{index}}{A logical or numeric vector giving the indices of the rows of the table of ranked models to print. This defaults to the full set of ranked models. It can be useful when the table of ranked models is large to examine a subset via this \code{index} argument, for display purposes.}
#' \item{\code{digits}}{The number of decimal places to round model selection criteria to (defaults to 3).}}
#' @note The \code{criterion} argument here need not comply with the criterion used for model selection within each \code{"MoEClust"} object, but be aware that a mismatch in terms of \code{criterion} \emph{may} require the optimal model to be re-fit in order to be extracted, thereby slowing down \code{\link{MoE_compare}}.
#' 
#' If random starts had been used via \code{init.z="random"} the \code{optimal} model may not necessarily correspond to the highest-ranking model in the presence of a criterion mismatch, due to the randomness of the initialisation. 
#'
#' A dedicated \code{print} function exists for objects of class \code{"MoECompare"}.
#' 
#' \code{\link{plot.MoEClust}} and \code{\link{as.Mclust}} can both also be called on objects of class \code{"MoECompare"}.
#'
#' @details The purpose of this function is to conduct model selection on \code{"MoEClust"} objects, fit to the same data set, with different combinations of gating/expert network covariates or different initialisation settings.
#'
#' Model selection will have already been performed in terms of choosing the optimal number of components and GPCM/\pkg{mclust} model type within each supplied set of results, but \code{\link{MoE_compare}} will respect the internal ranking of models when producing the final ranking if \code{optimal.only} is \code{FALSE}: otherwise only those models already deemed optimal within each \code{"MoEClust"} object will be ranked.
#'
#' As such if two sets of results are supplied when \code{optimal.only} is \code{FALSE}, the 1st, 2nd and 3rd best models could all belong to the first set of results, meaning a model deemed suboptimal according to one set of covariates could be superior to one deemed optimal under another set of covariates.
#' @return A list of class \code{"MoECompare"}, for which a dedicated print function exists, containing the following elements (each of length \code{pick}, and ranked according to \code{criterion}, where appropriate):
#' \item{\code{optimal}}{The single optimal model (an object of class \code{"MoEClust"}) among those supplied, according to the chosen \code{criterion}.}
#' \item{\code{pick}}{The final number of ranked models. May be different (i.e. less than) the supplied \code{pick} value.}
#' \item{\code{MoENames}}{The names of the supplied \code{"MoEClust"} objects.}
#' \item{\code{modelNames}}{The \code{\link[mclust]{mclustModelNames}}.}
#' \item{\code{G}}{The optimal numbers of components.}
#' \item{\code{df}}{The numbers of estimated parameters.}
#' \item{\code{iters}}{The numbers of EM/CEM iterations.}
#' \item{\code{bic}}{BIC values, ranked according to \code{criterion}.}
#' \item{\code{icl}}{TCL values, ranked according to \code{criterion}.}
#' \item{\code{aic}}{AIC values, ranked according to \code{criterion}.}
#' \item{\code{loglik}}{Maximal log-likelihood values, ranked according to \code{criterion}.}
#' \item{\code{gating}}{The gating formulas.}
#' \item{\code{expert}}{The expert formulas.}
#' \item{\code{algo}}{The algorithm used for fitting the model - either \code{"EM"}, \code{"CEM"}, \code{"cemEM"}.}
#' \item{\code{equalPro}}{Logical indicating whether mixing proportions were constrained to be equal across components.}
#' \item{\code{hypvol}}{Hypervolume parameters for the noise component if relevant, otherwise set to \code{NA} (see \code{\link{MoE_control}}).}
#' \item{\code{noise}}{Either a logical indicating the presence/absence of a noise component, or the type of noise component fitted (if any). Depends on the supplied value of \code{noise}. Only displayed if at least one of the compared models has a noise component.}
#' \item{\code{noise.gate}}{Logical indicating whether gating covariates were allowed to influence the noise component's mixing proportion. Only printed for models with a noise component, when at least one of the compared models has gating covariates, and even then only when \code{noise} is supplied as \code{TRUE}.}
#' \item{\code{equalNoise}}{Logical indicating whether the mixing proportion of the noise component for \code{equalPro} models is also equal (\code{TRUE}) or estimated (\code{FALSE}).}
#' @export
#' @keywords clustering main
#' @references K. Murphy and T. B. Murphy (2018). Gaussian Parsimonious Clustering Models with Covariates. \emph{To appear}. <\href{https://arxiv.org/abs/1711.05632}{arXiv:1711.05632}>.
#' @importFrom mclust "mclustModelNames"
#' @author Keefe Murphy - <\email{keefe.murphy@@ucd.ie}>
#'
#' @seealso \code{\link{MoE_clust}}, \code{\link[mclust]{mclustModelNames}}, \code{\link{plot.MoEClust}}, \code{\link{as.Mclust}}
#' @usage
#' MoE_compare(...,
#'             criterion = c("bic", "icl", "aic"),
#'             pick = 10L,
#'             optimal.only = FALSE)
#' @examples
#' data(CO2data)
#' CO2   <- CO2data$CO2
#' GNP   <- CO2data$GNP
#' 
#' m1    <- MoE_clust(CO2, G=1:3)
#' m2    <- MoE_clust(CO2, G=2:3, gating= ~ GNP)
#' m3    <- MoE_clust(CO2, G=1:3, expert= ~ GNP)
#' m4    <- MoE_clust(CO2, G=2:3, gating= ~ GNP, expert= ~ GNP)
#' m5    <- MoE_clust(CO2, G=2:3, equalPro=TRUE)
#' m6    <- MoE_clust(CO2, G=2:3, expert= ~ GNP, equalPro=TRUE)
#' m7    <- MoE_clust(CO2, G=2:3, expert= ~ GNP, tau0=0.1)
#'
#' # Rank only the optimal models and examine the best model
#' (comp <- MoE_compare(m1, m2, m3, m4, m5, m6, m7, pick=7, optimal.only=TRUE))
#' (best <- comp$optimal)
#' (summ <- summary(best))
#'
#' # Examine all models visited, including those already deemed suboptimal
#' # Only print models with expert covariates & more than one component
#' comp2 <- MoE_compare(m1, m2, m3, m4, m5, m6, m7, pick=Inf)
#' print(comp2, comp2$expert != "None" & comp2$G > 1)
  MoE_compare     <- function(..., criterion = c("bic", "icl", "aic"), pick = 10L, optimal.only = FALSE) {
    crit.miss     <- missing(criterion)
    if(!missing(criterion) && (length(criterion) > 1 ||
       !is.character(criterion)))                 stop("'criterion' must be a single character string", call.=FALSE)
    criterion     <- match.arg(criterion)
    num.miss      <- missing(pick)
    opt.miss      <- missing(optimal.only)
    if(length(pick)    != 1            ||
       !is.numeric(pick))                         stop("'pick' must be a single number", call.=FALSE)
    if(floor(pick)     != pick         ||
       pick        < 1)                           stop("'pick' must be a strictly positive integer", call.=FALSE)
    if(length(optimal.only)      > 1   ||
       !is.logical(optimal.only))                 stop("'optimal.only' must be a single logical indicator", call.=FALSE)
    call          <- match.call(expand.dots=TRUE)[-1L]
    call          <- if(crit.miss) call else call[-which(names(call) == "criterion")]
    call          <- if(num.miss)  call else call[-which(names(call) == "pick")]
    call          <- if(opt.miss)  call else call[-which(names(call) == "optimal.only")]
    len.call      <- length(as.list(call))
    if(len.call   == 1 && inherits(..., "list") && !inherits(..., "MoEClust")) {
      dots        <- as.list(...)
      mod.names   <- unique(names(dots))
      comparison  <- sapply(dots, inherits, "MoECompare", logical(1L))
      dots[comparison]          <- sapply(dots[comparison], "[", "optimal")
      MoEs        <- dots[mod.names]
      if(is.null(mod.names))                      stop("When supplying models as a list, every element of the list must be named", call.=FALSE)
    } else {
      lapply(call, function(x) if(inherits(x, "MoECompare")) x$optimal else x)
      dots        <- list(...)
      mod.names   <- vapply(call, deparse, character(1L))
      comparison  <- sapply(dots, inherits, "MoECompare", logical(1L))
      dots[comparison]          <- sapply(dots[comparison], "[", "optimal")
      MoEs        <- stats::setNames(dots, mod.names)
      mod.names   <- unique(mod.names)
      MoEs        <- MoEs[mod.names]
    }
    Mclass        <- vapply(MoEs, class,          character(1L))
    if(any(Mclass != "MoEClust"))                 stop("All models must be of class 'MoE_clust'!", call.=FALSE)
    if(length(unique(lapply(MoEs,
       "[[", "data")))          != 1)             stop("All models being compared must have been fit to the same data set!", call.=FALSE)
    title         <- "Comparison of Gaussian Parsimonious Clustering Models with Covariates"
    dat.name      <- deparse(MoEs[[1L]]$call$data)
    gate.x        <- lapply(MoEs, "[[", "gating")
    algo          <- sapply(MoEs,   attr, "Algo")
    equalNoise    <- sapply(MoEs,   attr, "EqualNoise")
    equalPro      <- sapply(MoEs,   attr, "EqualPro")
    noise.gate    <- sapply(MoEs,   attr, "NoiseGate")
    gating        <- lapply(gate.x, attr, "Formula")
    expert        <- lapply(lapply(MoEs, "[[", "expert"), attr, "Formula")
    hypvol        <- hypvol0    <- lapply(MoEs, "[[", "hypvol")
    noise.meth    <- sapply(hypvol, attr, "Meth")
    noise.g0only  <- sapply(hypvol, attr, "g0only")
    noise.null    <- vapply(noise.meth,   is.null,      logical(1L))
    noise.onlyg0  <- vapply(noise.g0only, isTRUE,       logical(1L))
    equalNoise[noise.null]      <-
    noise.gate[noise.null]      <- NA
    noise.meth[noise.null]      <- FALSE
    hypvol        <- unlist(hypvol)
    noise.meth    <- unlist(noise.meth)
    BICs          <- lapply(MoEs, "[[", "BIC")
    ICLs          <- lapply(MoEs, "[[", "ICL")
    AICs          <- lapply(MoEs, "[[", "AIC")
    LLxs          <- lapply(MoEs, "[[", "LOGLIK")
    DFxs          <- lapply(MoEs, "[[", "DF")
    ITxs          <- lapply(MoEs, "[[", "ITERS")
    choice        <- max(lengths(BICs))
    bics          <- lapply(BICs, function(x) .pick_MoECrit(x, choice)$crits)
    icls          <- lapply(ICLs, function(x) .pick_MoECrit(x, choice)$crits)
    aics          <- lapply(AICs, function(x) .pick_MoECrit(x, choice)$crits)
    llxs          <- lapply(LLxs, function(x) .pick_MoECrit(x, choice)$crits)
    dfxs          <- lapply(DFxs, function(x) .pick_MoECrit(x, choice)$crits)
    itxs          <- lapply(ITxs, function(x) .pick_MoECrit(x, choice)$crits)
    if(optimal.only) {
      opt.names   <- names(.crits_names(lapply(switch(EXPR=criterion, bic=bics, icl=icls, aic=aics), "[", 1L)))
    }
    bics          <- .crits_names(bics)
    icls          <- .crits_names(icls)
    aics          <- .crits_names(aics)
    llxs          <- .crits_names(llxs)
    dfxs          <- .crits_names(dfxs)
    itxs          <- .crits_names(itxs)
    if(optimal.only) {
      bics        <- bics[names(bics) %in% opt.names]
      icls        <- icls[names(icls) %in% opt.names]
      aics        <- aics[names(aics) %in% opt.names]
      llxs        <- llxs[names(llxs) %in% opt.names]
      dfxs        <- dfxs[names(dfxs) %in% opt.names]
      itxs        <- itxs[names(itxs) %in% opt.names]
    }
    crits         <- switch(EXPR=criterion, bic=bics, icl=icls, aic=aics)
    pick          <- min(pick, length(crits))
    max.crits     <- sort(crits, decreasing=TRUE)[seq_len(pick)]
    if(length(unique(max.crits)) < pick) {
      ties        <- max.crits  == max.crits[1L]
      if(any(ties[-1L]))         {                warning(paste0("Ties for the optimal model exist according to the '", criterion, "' criterion: choosing the most parsimonious model\n"), call.=FALSE, immediate.=TRUE)
        df.ties   <- dfxs[names(max.crits)][which(ties)]
        max.crits[ties]   <- max.crits[order(df.ties)]
        if(any((df.ties   == df.ties[1L])[-1L])) {
          max.crits[ties] <- max.crits[order(as.numeric(gsub(".*,", "", names(max.crits[ties]))))]
        }
      } else                                      warning(paste0("Ties exist according to the '", criterion, "' criterion\n"), call.=FALSE, immediate.=TRUE)
    }
    max.names     <- names(max.crits)
    crit.names    <- gsub("\\|.*", "",          max.names)
    G             <- as.numeric(gsub(".*,", "", max.names))
    gating        <- unname(unlist(gating[crit.names]))
    expert        <- unname(unlist(expert[crit.names]))
    modelNames    <- gsub(",.*", "", gsub(".*\\|", "", max.names))
    best.model    <- MoEs[[crit.names[1L]]]
    if(best.model$modelName     != modelNames[1L] || best.model$G != G[1L]) {
      bestGN      <- gating[1L] != "~1" && (best.model$G + !noise.null[crit.names[1L]] - isFALSE(noise.gate[crit.names[1L]])) <= 1
      best.model$net.covs       <- if(bestGN) attr(best.model$net.covs, "Discarded") else best.model$net.covs
      cat("Re-fitting optimal model due to mismatched 'criterion'...\n\n")
      old.call    <- best.model$call
      old.call    <- c(as.list(old.call)[1L], list(criterion=criterion), as.list(old.call)[-1L])
      old.call    <- as.call(old.call[!duplicated(names(old.call))])
      if(!is.null(old.call$init.z)      &&
         old.call$init.z == "random")             warning("Optimal model may differ slightly due to criterion mismatch and random starts used in the initialisation:\nPrinted output intended only as a guide", call.=FALSE, immediate.=TRUE)
      best.call   <- c(list(data=best.model$data, modelNames=modelNames[1L], G=G[1L], verbose=FALSE, network.data=best.model$net.covs), as.list(old.call[-1L]))
      best.mod    <- try(do.call(MoE_clust, best.call[!duplicated(names(best.call))]), silent=TRUE)
      if(!inherits(best.model, "try-error")) {
        best.model$call                <- old.call
        best.model$modelName           <- best.mod$modelName
        best.model$G                   <- best.mod$G
        best.model$bic                 <- best.mod$bic
        best.model$icl                 <- best.mod$icl
        best.model$aic                 <- best.mod$aic
        best.model$gating              <- best.mod$gating
        best.model$expert              <- best.mod$expert
        best.model$loglik              <- best.mod$loglik
        best.model$linf                <- best.mod$linf
        best.model$df                  <- best.mod$df
        best.model$iters               <- best.mod$iters
        best.model$hypvol              <- best.mod$hypvol
        best.model$parameters          <- best.mod$parameters
        best.model$z                   <- best.mod$z
        best.model$classification      <- best.mod$classification
        best.model$uncertainty         <- best.mod$uncertainty
        best.model$net.covs            <- best.mod$net.covs
        best.model$resid.data          <- best.mod$resid.data
        attributes(best.model)         <- attributes(best.mod)
      } else best.model                <- paste0("Failed to re-fit the optimal model: ", gsub("\"", "'", deparse(old.call, width.cutoff=500L), fixed=TRUE))
    }
    gating2       <- replace(gating, gating == "~1", "None")
    gating[gating == "~1" | G   <= 1]  <- "None"
    expert[expert == "~1"]             <- "None"
    if(any(G == 0,  !noise.null))       {
      noise.meth0 <- sapply(hypvol0, attr, "Meth0")
      hypvol0     <- sapply(hypvol0, attr, "Hypvol0")
    }
    hypvol        <- ifelse(G   == 0, hypvol0[crit.names], ifelse(noise.onlyg0[crit.names], NA, hypvol[crit.names]))
    noise.meth    <- ifelse(is.na(hypvol), "FALSE",        noise.meth0[crit.names])
    noise.gate    <- ifelse(is.na(hypvol), NA,             noise.gate[crit.names])
    Gtmp          <- G == 1  & noise.meth != "FALSE" & noise.gate
    gating[Gtmp]  <- gating2[Gtmp] 
    equalPro      <- replace(unname(equalPro[crit.names]), gating != "None" | G  <= 1, NA)
    equalNoise    <- ifelse(is.na(hypvol) | G <= 1, NA,    equalNoise[crit.names] & vapply(equalPro, isTRUE, logical(1L)))
    comp          <- list(title = title, data = dat.name, optimal = best.model, pick = pick, MoENames = crit.names, modelNames = modelNames, G = as.integer(G), df = as.integer(unname(dfxs[max.names])), iters = as.integer(unname(itxs[max.names])), 
                          bic = unname(bics[max.names]), icl = unname(icls[max.names]), aic = unname(aics[max.names]), loglik = unname(llxs[max.names]), gating = gating, expert = expert, algo = unname(algo[crit.names]), equalPro = equalPro, 
                          hypvol = unname(hypvol), noise = unname(noise.meth), noise.gate = unname(replace(noise.gate, gating == "None" | G <= 1, NA)), equalNoise = unname(replace(equalNoise, !equalPro | is.na(equalPro), NA)))
    class(comp)   <- c("MoECompare", "MoEClust")
    bic.tmp       <- sapply(BICs, as.vector)
    attr(comp, "Crit")   <- criterion
    attr(comp, "Opt")    <- optimal.only
    attr(comp, "NMods")  <- c(tried = sum(vapply(bic.tmp, function(x) length(x[!is.na(x)]),    numeric(1L))),
                              ran   = sum(vapply(bic.tmp, function(x) length(x[is.finite(x)]), numeric(1L))))
      comp
  }

#' Predictions for MoEClust models
#'
#' Predicts both cluster membership probability and fitted response values from a \code{MoEClust} model, using covariates and response data, or covariates only. The MAP classification is also reported in both cases.
#' @param object An object of class \code{"MoEClust"} generated by \code{\link{MoE_clust}}, or an object of class \code{"MoECompare"} generated by \code{\link{MoE_compare}}. Predictions for models with a noise component are facilitated here too.
#' @param newdata A list with two \emph{named} components, each of which must be a \code{data.frame} or \code{matrix} with named columns, giving the data for which predicitions are desired.
#' \describe{
#' \item{\code{new.x}}{The new covariates for the \code{gating} &/or \code{expert} networks. \strong{Must} be supplied when \code{newdata$new.y} is supplied.}
#' \item{\code{new.y}}{(Optional) response data. When supplied, cluster and response prediction is based on both \code{newdata$new.x} and \code{newdata$new.y}, otherwise only on the covariates in \code{newdata$new.x}.}
#' }
#' If supplied as a list with elements \code{new.x} and \code{new.y}, both \strong{must} have the same number of rows.
#'
#' Alternatively, a single \code{data.frame} or \code{matrix} can be supplied and an attempt will be made to extract & separate covariate and response columns (\emph{if any}) into \code{newdata$new.x} and \code{newdata$new.y} based on the variable names in \code{object$data} and \code{object$net.covs}.
#'
#' When \code{newdata} is not supplied in any way, the covariates and response variables used in the fitting of the model are used here.
#' @param resid A logical indicating whether to return the residuals also. Defaults to \code{FALSE}. Only allowed when response variables are supplied in some form. The function \code{residuals} is a wrapper to \code{predict} with the argument \code{resid} set to \code{TRUE}, with only the residuals returned.
#' @param ... Catches unused arguments.
#'
#' @return A list with the following named components, regardless of whether \code{newdata$new.x} and \code{newdata$new.y} were used, or \code{newdata$new.x} only.
#' \item{\code{y}}{Fitted values of the response variables.}
#' \item{\code{z}}{A matrix whose \code{[i,k]}-th entry is the probability that observation \emph{i} of the \code{newdata} belonds to the \emph{k}-th component. For models with a noise component, the final column gives the probability of belonging to the so-called \emph{Cluster0}.}
#' \item{\code{classification}}{The vector of predicted cluster labels for the \code{newdata}. \code{0} is returned for observations assigned to the noise component.}
#'
#' When \code{residuals} is called, only the residuals are returned; when \code{predict} is called with \code{resid=TRUE}, the list above will also contain the element \code{resids}, containing the residuals.
#'
#' @note Predictions can also be made for models with noise components, in which case \code{z} will include the probability of belonging to \code{"Cluster0"} & \code{classification} will include labels with the value \code{0} for observations classified as noise (if any).
#' @references K. Murphy and T. B. Murphy (2018). Gaussian Parsimonious Clustering Models with Covariates. \emph{To appear}. <\href{https://arxiv.org/abs/1711.05632}{arXiv:1711.05632}>.
#' @author Keefe Murphy - <\email{keefe.murphy@@ucd.ie}>
#' @seealso \code{\link{MoE_clust}}
#' @method predict MoEClust
#' @keywords clustering main
#' @importFrom matrixStats "rowSums2"
#' @export
#' @usage
#' \method{predict}{MoEClust}(object,
#'         newdata,
#'         resid = FALSE,
#'         ...)
#' @examples
#' data(ais)
#'
#' # Fit a MoEClust model and predict the same data
#' res     <- MoE_clust(ais[,3:7], G=2, gating=~BMI, expert=~sex,
#'                      modelNames="EVE", network.data=ais)
#' pred1   <- predict(res)
#' pred1$classification
#'
#' # Remove some rows of the data for prediction purposes
#' ind     <- sample(1:nrow(ais), 5)
#' dat     <- ais[-ind,]
#'
#' # Fit another MoEClust model to the retained data
#' res2    <- MoE_clust(dat[,3:7], G=3, gating=~BMI + sex,
#'                      modelNames="EEE", network.data=dat)
#'
#' # Predict held back data using the covariates & response variables
#' pred2   <- predict(res2, newdata=ais[ind,])
#' # pred2 <- predict(res2, newdata=list(new.y=ais[ind,3:7],
#' #                                     new.x=ais[ind,c("BMI", "sex")]))
#' pred2$y  # predicted response
#'
#' # Get the residuals
#' residuals(res2, newdata=ais[ind,])
#'
#' # Predict held back data using only the covariates
#' pred3   <- predict(res2, newdata=list(new.x=ais[ind,c("BMI", "sex")]))
#' # pred3 <- predict(res2, newdata=ais[ind,c("BMI", "sex")])
#' pred3$z  # predicted cluster membership probablities
predict.MoEClust  <- function(object, newdata = list(...), resid = FALSE, ...) {
  object          <- if(inherits(object, "MoECompare")) object$optimal else object
  if(length(resid) > 1   || !is.logical(resid))   stop("'resid' should be a single logical indicator", call.=FALSE)
  net             <- object$net.covs
  dat             <- object$data
  datnames        <- colnames(dat)
  hypvol          <- object$hypvol
  noise           <- !is.na(hypvol)
  yM              <- FALSE
  if(nmiss        <- ifelse(inherits(newdata,
                                     "list")    &&
                     length(newdata)  == 0,
                     missing(newdata), FALSE))   {
    newdata.x     <- net
    newdata.y     <- dat
    nL            <- FALSE
  } else if(inherits(newdata, "list") -> nL)     {
    if(is.null(newdata$new.y)  -> yM)  {
      if(length(newdata) != 1  ||
         names(newdata)  != "new.x")              stop("'newdata' must be a list with a single component named 'new.x' if it does not also contain the component named 'new.y'", call.=FALSE)
      newdata.x   <- newdata$new.x
      if(!is.matrix(newdata.x) &&
         !is.data.frame(newdata.x))               stop("'newdata$new.x' must be a 'matrix' or 'data.frame'", call.=FALSE)
      if(is.null(colnames(newdata.x)))            stop("'newdata$new.x' must have named columns", call.=FALSE)
    } else     {
      if(length(newdata) != 2  ||
         !is.element(names(newdata),
                     c("new.x", "new.y")))        stop("If 'newdata' is a list, it must be of length two, with components named 'new.x' and 'new.y'", call.=FALSE)
      if(any(vapply(newdata, function(x)
             !is.matrix(x)     &&
             !is.data.frame(x),
             logical(1L))))                       stop("Both 'new.x' and 'new.y' within 'newdata' must be of class 'matrix' or 'data.frame'", call.=FALSE)
      if(any(vapply(newdata, function(x)
             is.null(colnames(x)),
             logical(1L))))                       stop("Both 'new.x' and 'new.y' within 'newdata' must have named columns", call.=FALSE)
      newdata.x   <- newdata$new.x
      newdata.y   <- newdata$new.y
      if(nrow(newdata.x) != nrow(newdata.y))      stop("'new.x' and 'new.y' within 'newdata' must have the same number of rows", call.=FALSE)
    }
  } else       {
    if(!is.matrix(newdata)     &&
       !is.data.frame(newdata))                   stop("'newdata' must be either be a 'list', 'matrix', or 'data.frame'",   call.=FALSE)
    newdata.x     <- newdata.y <- newdata
  }
  Xexp            <- attr(object, "Expert")
  newdata.x       <- as.data.frame(newdata.x)
  x.names         <- names(newdata.x)
  netnames        <- names(net)
  newdata.x       <- newdata.x[,x.names %in% netnames,      drop=FALSE]
  if(!all(netnames     %in% x.names))             stop("Covariates missing in 'newdata'", call.=FALSE)
  gate            <- attr(net, "Gating")
  expx            <- attr(net, "Expert")
  cnew            <- stats::complete.cases(newdata.x)
  newdata.x       <- newdata.x[cnew,netnames,               drop=FALSE]
  if(!all(cnew))                                  message("Rows with missing values discarded from 'new.x'\n")
  newgate         <- newdata.x[,if(!any(is.na(gate))) gate, drop=FALSE]
  newexpx         <- newdata.x[,if(!any(is.na(expx))) expx, drop=FALSE]
  gatenames       <- colnames(newgate)
  expxnames       <- colnames(newexpx)
  if(any(vapply(seq_len(ncol(newgate)),
     function(p, gp=newgate[,p])  is.factor(gp) &&
     !identical(levels(gp),
     levels(net[,gatenames[p]])), logical(1L))))  warning("One of more categorical gating covariates in the unseen newdata has new factor levels\n", call.=FALSE, immediate.=TRUE)
  if(any(vapply(seq_len(ncol(newexpx)),
     function(p, ep=newexpx[,p])  is.factor(ep) &&
     !identical(levels(ep),
     levels(net[,expxnames[p]])), logical(1L))))  warning("One of more categorical expert covariates in the unseen newdata has new factor levels\n", call.=FALSE, immediate.=TRUE)
  rownames(newdata.x)          <- NULL
  nr              <- nrow(newdata.x)
  nrseq           <- seq_len(nr)
  G               <- object$G
  GN              <- G + noise
  Gseq            <- seq_len(G)
  params          <- object$parameters
  if(G == 0)   {
    retval        <- list(ystar=matrix(attr(hypvol, "Location"), nrow=nr, ncol=object$d, byrow=TRUE),
                          classification=rep(0L, nr),
                          zstar=provideDimnames(matrix(1L, nrow=nr, ncol=1L), base=list(as.character(nrseq), "Cluster0")))
      return(if(isTRUE(resid)) c(retval, list(resids=newdata.y - retval$ystar)) else retval)
  }
  pred.exp        <- lapply(object$expert, stats::predict, newdata=newexpx)
  if(nmiss)    {
    zstar         <- object$z
  } else       {
    gating        <- object$gating
    if(attr(gating, "EqualPro"))  {
      new.tau     <- matrix(params$pro, nrow=nrow(newgate), ncol=GN)
    } else if(noise      && attr(gating, "NoiseGate"))  {
      new.tau     <- matrix(stats::predict(gating, type=ifelse(GN > 1, "probs", "response"), newdata=newgate), ncol=GN)
    } else     {
      new.tau     <- matrix(stats::predict(gating, type=ifelse(G  > 1, "probs", "response"), newdata=newgate), ncol=G)
      new.tau     <- if(noise) .tau_noise(new.tau, if(attr(object, "Gating")) params$pro[1L,GN] else params$pro[GN]) else new.tau
    }
    if(!yM)    {
      newdata.y   <- as.data.frame(newdata.y)
      y.names     <- names(newdata.y)
      newdata.y   <- newdata.y[cnew,y.names %in% datnames,  drop=FALSE]
      if(!all(datnames %in% y.names))  {
        if(nL) {                                  stop("Response variables missing in 'newdata'",      call.=FALSE)
        } else {                                  warning("Response variables missing in 'newdata'\n", call.=FALSE, immediate.=TRUE)
          yM      <- TRUE
        }
      } else   {
        newdata.y <- newdata.y[,datnames,                   drop=FALSE]
        rownames(newdata.y)    <- NULL
      }
    }
    if(resid  && !(resid <- !yM))                 warning("'resid' can only be TRUE when response variables are supplied\n", call.=FALSE, immediate.=TRUE)
    attr(net, "Gating")  <-
    attr(net, "Expert")  <-
    attr(net, "Both")    <-
    rownames(net)        <-
    rownames(dat)        <- NULL
    if(identical(newdata.x, net)      &&
       !yM    && 
       identical(newdata.y, dat))      {
      zstar       <- object$z
      nmiss       <- TRUE
    } else if(yM)         {
      zstar       <- new.tau
    } else     {
      Xexp        <- attr(object, "Expert")
      newdata.y   <- if(Xexp) lapply(pred.exp, "-", newdata.y)    else newdata.y
      mus         <- if(Xexp) 0L                                  else params$mean
      zstar       <- MoE_estep(Dens=MoE_dens(data=newdata.y, mus=mus, sigs=params$variance, log.tau=log(new.tau), Vinv=params$Vinv))$z
    }
  }
  if(!noise || attr(hypvol, "Meth") == "manual") {
    zstar2        <- if(noise) .renorm_z(zstar[,-GN, drop=FALSE]) else zstar
    ystar         <- as.matrix(Reduce("+", lapply(Gseq, function(g) zstar2[,g] * pred.exp[[g]])))
  } else     {
    ystar         <- as.matrix(Reduce("+", lapply(Gseq, function(g) zstar[,g]  * pred.exp[[g]])))
    ystar         <- ystar + zstar[,GN] *  matrix(attr(hypvol, "Location"), nrow=nr, ncol=object$d, byrow=TRUE)
  }
  claX            <- stats::setNames(max.col(zstar), nrseq)
  claX[claX   == G + 1]  <- 0L
  rownames(ystar) <- nrseq
  colnames(ystar) <- datnames
  gnames          <- paste0("Cluster", Gseq)
  colnames(zstar) <- if(noise) c(gnames, "Cluster0") else gnames
  retval          <- list(y=ystar, classification=claX, z=zstar)
    return(if(isTRUE(resid)) c(retval, list(resids=newdata.y - ystar)) else retval)
}

#' @rdname predict.MoEClust
#' @method residuals MoEClust
#' @keywords clustering utility
#' @importFrom matrixStats "rowSums2"
#' @usage
#' \method{residuals}{MoEClust}(object,
#'           newdata,
#'           ...)
#' @export
  residuals.MoEClust     <- function(object, newdata = list(...), ...) {
    resids <- do.call(predict.MoEClust, c(as.list(match.call())[-1L], list(resid=TRUE)))$resids
      if(!is.null(resids))  return(resids)
  }

#' Convert MoEClust objects to the Mclust class
#'
#' Converts an object of class \code{"MoEClust"} generated by \code{\link{MoE_clust}} and converts it to an object of class \code{"Mclust"} as generated by fitting \code{\link[mclust]{Mclust}}, to facilitate use of plotting and other functions for the \code{"Mclust"} class within the \pkg{mclust} package. Some caution is advised when converting models with gating &/or expert covariates (see Note below).
#' @param x An object of class \code{"MoEClust"} generated by \code{\link{MoE_clust}} or an object of class \code{"MoECompare"} generated by \code{\link{MoE_compare}}. Models with a noise component are facilitated here too.
#' @param expert.covar Logical (defaults to \code{TRUE}) governing whether the extra variabilitity in the component means is added to the MVN ellipses corresponding to the component covariance matrices in the presence of expert network covariates. See the function \code{\link{expert_covar}}.
#' @param signif Significance level for outlier removal. Must be a single number in the interval [0, 1). Corresponds to the percentage of data to be considered extreme and therefore removed (half of \code{signif} at each endpoint, on a column-wise basis). The default, \code{0}, corresponds to no outlier removal. \strong{Only} invoke this argument as an aid to visualisation via \code{\link[mclust]{plot.Mclust}}.
#' @param ... Further arguments to be passed to other methods.
#'
#' @return An object of class \code{"Mclust"}. See \code{methods(class="Mclust")} for a (non-exhaustive) list of functions which can be applied to this class.
#' @details Of course, the user is always encouraged to use the dedicated \code{\link[=plot.MoEClust]{plot}} function for objects of the \code{"MoEClust"} class instead, but calling \code{plot} after converting via \code{\link{as.Mclust}} can be particularly useful for univariate mixtures.
#'
#' In the presence of expert network covariates, the component-specific covariance matrices are (by default, via the argument \code{expert.covar}) modified for plotting purposes via the function \code{\link{expert_covar}}, in order to account for the extra variability of the means, usually resulting in bigger shapes & sizes for the MVN ellipses.
#'
#' The \code{signif} argument is intended only to aid visualisation via \code{\link[mclust]{plot.Mclust}}, as plots therein can be sensitive to outliers, particularly with regard to axis limits.
#' @note Of the functions which can be applied to the result of the conversion, \code{\link[mclust]{logLik.Mclust}} shouldn't be trusted in the presence of either expert network covariates, or (for models with more than 1 component) gating network covariates.
#'
#' Mixing proportions are averaged over observations in components in the presence of gating network covariates during the coercion.
#'
#' Plots may be quite misleading in the presence of gating &/or expert covariates when the \code{what} argument is \code{"density"} within \code{\link[mclust]{plot.Mclust}}; users are \strong{strongly} encouraged to use \code{\link{MoE_gpairs}} with \code{response.type="density"} instead.
#'
#' The functions \code{\link[mclust]{clustCombi}} and \code{\link[mclust]{clustCombiOptim}} can be safely used (provided \code{as.Mclust(x)} is supplied as the \code{object} argument to \code{\link[mclust]{clustCombi}}), as they only rely on \code{x$z} and \code{x$G} only. See the examples below.
#' @importFrom mclust "as.densityMclust.Mclust" "clustCombi" "clustCombiOptim" "logLik.Mclust" "icl" "plot.Mclust" "plot.mclustBIC" "plot.mclustICL" "predict.Mclust" "print.Mclust" "summary.Mclust"
#' @export
#' @seealso \code{\link[mclust]{Mclust}}, \code{\link[mclust]{plot.Mclust}}, \code{\link{MoE_clust}}, \code{\link{plot.MoEClust}}, \code{\link{expert_covar}}
#' @author Keefe Murphy - <\email{keefe.murphy@@ucd.ie}>
#' @references C. Fraley and A. E. Raftery (2002). Model-based clustering, discriminant analysis, and density estimation. \emph{Journal of the American Statistical Association}, 97:611-631.
#' @keywords utility
#' @usage
#' as.Mclust(x,
#'           expert.covar = TRUE,
#'           signif = 0L,
#'           ...)
#' @examples
#' \dontrun{
#' library(mclust)
#' 
#' # Fit a gating network mixture of experts model to the ais data
#' data(ais)
#' mod   <- MoE_clust(ais[,3:7], G=1:9, gating= ~ BMI + sex, network.data=ais)
#'
#' # Convert to the "Mclust" class and examine the classification
#' mod2  <- as.Mclust(mod)
#' plot(mod2, what="classification")
#'
#' # Examine the uncertainty
#' plot(mod2, what="uncertainty")
#'
#' # Return the optimal number of clusters according to entropy
#' combi <- mclust::clustCombi(object=mod2)
#' optim <- mclust::clustCombiOptim(combi)
#' table(mod2$classification, ais$sex)
#' table(optim$cluster.combi, ais$sex)
#'
#' # While we could have just used plot.MoEClust above,
#' # plot.Mclust is especially useful for univariate data
#' data(CO2data)
#' res <- MoE_clust(CO2data$CO2, G=3, equalPro=TRUE, expert = ~ GNP, network.data=CO2data)
#' plot(as.Mclust(res))}
  as.Mclust       <- function(x, expert.covar = TRUE, signif = 0L, ...) {
      UseMethod("as.Mclust")
  }

#' @method as.Mclust MoEClust
#' @importFrom matrixStats "colMeans2"
#' @importFrom mclust "sigma2decomp"
#' @export
  as.Mclust.MoEClust      <- function(x, expert.covar = TRUE, signif = 0L, ...) {
    x             <- if(inherits(x, "MoECompare")) x$optimal else x
    if(length(signif) > 1 || !is.numeric(signif)     ||
       signif < 0 || signif   >= 1)               stop("'signif' must be a single number in the interval [0, 1)", call.=FALSE)
    uni           <- x$d  == 1
    gating        <- attr(x, "Gating")
    expert        <- attr(x, "Expert")
    x$data        <- as.matrix(x$data)
    x$loglik      <- x$loglik[length(x$loglik)]
    x$BIC         <- replace(x$BIC, !is.finite(x$BIC), NA)
    class(x$BIC)  <- "mclustBIC"
    x$uncertainty         <- if(uni)                 unname(x$uncertainty) else x$uncertainty
    x$classification      <- if(uni)              unname(x$classification) else x$classification
    x$parameters$pro      <- if(gating)        colMeans2(x$parameters$pro) else x$parameters$pro
    x$parameters$variance <- if(isTRUE(expert.covar)    &&
                                expert)  suppressWarnings(expert_covar(x)) else x$parameters$variance
    x$data        <- if(signif > 0)   apply(x$data, 2L, .trim_out, signif) else x$data
    colnames(x$z) <- NULL
    x             <- x[-which(is.element(names(x), c("ICL", "icl", "AIC", "aic", "gating", "expert", "LOGLIK", "linf", "iters", "net.covs", "resid.data", "DF", "ITERS")))]
    name.x        <- names(x)
    attributes(x) <- NULL
    names(x)      <- name.x
    class(x)      <- "Mclust"
      x
  }

#' Account for extra variability in covariance matrices with expert covariates
#'
#' In the presence of expert network covariates, this helper function modifies the component-specific covariance matrices of a \code{"MoEClust"} object, in order to account for the extra variability of the means, usually resulting in bigger shapes & sizes for the MVN ellipses. The function also works for univariate response data.
#' @param x An object of class \code{"MoEClust"} generated by \code{\link{MoE_clust}}, or an object of class \code{"MoECompare"} generated by \code{\link{MoE_compare}}. Models with a noise component are facilitated here too.
#'
#' @details This function is used internally by \code{\link{plot.MoEClust}} and \code{\link{as.Mclust}}, for visualisation purposes.
#' @note The \code{modelName} of the resulting \code{variance} object may not correspond to the model name of the \code{"MoEClust"} object, in particular scale, shape, &/or orientation may no longer be constrained across clusters. Usually, the \code{modelName} of the transformed \code{variance} object will be "\code{VVV}".
#' @return The \code{variance} component only from the \code{parameters} list from the output of a call to \code{\link{MoE_clust}}, modified accordingly.
#' @seealso \code{\link{MoE_clust}}, \code{\link{MoE_gpairs}}, \code{\link{plot.MoEClust}}, \code{\link{as.Mclust}}
#' @references K. Murphy and T. B. Murphy (2018). Gaussian Parsimonious Clustering Models with Covariates. \emph{To appear}. <\href{https://arxiv.org/abs/1711.05632}{arXiv:1711.05632}>.
#' @author Keefe Murphy - <\email{keefe.murphy@@ucd.ie}>
#' @keywords utility
#' @export
#'
#' @examples
#' data(ais)
#' res   <- MoE_clust(ais[,3:7], G=2, gating= ~ 1, expert= ~ sex,
#'                    network.data=ais, modelNames="EEE", equalPro=TRUE)
#'
#' # Extract the variance object
#' res$parameters$variance
#'
#' # Modify the variance object
#' expert_covar(res)
  expert_covar    <- function(x) {
      UseMethod("expert_covar")
  }

#' @method expert_covar MoEClust
#' @importFrom mclust "sigma2decomp"
#' @export
  expert_covar.MoEClust   <- function(x) {
    x             <- if(inherits(x, "MoECompare")) x$optimal else x
    x.sig         <- x$parameters$variance
    d             <- x$d
    G             <- x$G
    if(attr(x, "Expert"))  {
      pred.var    <- unlist(lapply(x$expert, function(expert) stats::cov(as.matrix(stats::predict(expert)))))
     #pred.var    <- unlist(lapply(seq_len(G), function(g) stats::cov.wt(stats::predict(experts[[g]]), x$z[,g])$cov))
      if(d  == 1)  {
        x.sig$sigmasq     <- unname(x.sig$sigmasq + sqrt(pred.var))
        if(x$modelName    == "V" || length(unique(x.sig$sigmasq)) > 1) {
          x.sig$scale     <- x.sig$sigmasq
        }
      } else {
        x.sig     <- suppressWarnings(sigma2decomp(x.sig$sigma + array(pred.var, dim=c(d, d, G))))
      }
    } else                                        message("No expert covariates: returning the variance object without modification\n")
      return(x.sig)
  }

#' Force diagonal elements of a triangular matrix to be positive
#'
#' This function ensures that the triangular matrix in a QR (or other) decomposition has positive values along its diagonal.
#' @param x A matrix, which must be either upper-triangular or lower-triangular.
#'
#' @return An upper or lower triangular matrix with positive diagonal entries such that the matrix is still a valid decomposition of the matrix the input \code{x} is a decomposition of.
#' @export
#' @author Keefe Murphy - <\email{keefe.murphy@@ucd.ie}>
#' @keywords utility
#'
#' @examples
#' data(ais)
#' res <- MoE_clust(ais[,3:7], G=3, modelNames="EEE")
#' sig <- res$parameters$variance
#' a   <- force_posiDiag(sig$cholSigma)
#' b   <- chol(sig$Sigma)
#' round(sum(a - b), 10) == 0          #TRUE
#' sum(crossprod(a) != sig$Sigma) == 0 #TRUE
#' sum(crossprod(b) != sig$Sigma) == 0 #TRUE
  force_posiDiag  <- function(x) {
    if(!is.matrix(x) ||
       any(x[row(x)   > col(x)] != 0) &&
       any(x[col(x)   > row(x)] != 0))            stop("'x' must be an upper or lower triangular matrix")
      diag(sign(diag(x))) %*% x
  }

#' Quantile-Based Clustering for Univariate Data
#'
#' Returns a quantile-based clustering for univariate data.
#' @param x A vector of numeric data.
#' @param G The desired number of clusters.
#'
#' @return The vector of cluster labels.
#' @export
#' @keywords clustering
#' @usage
#' quant_clust(x,
#'             G)
#' @examples
#' data(CO2data)
#' quant_clust(CO2data$CO2, G=2)
  quant_clust     <- function(x, G) {
    if((is.data.frame(x) || is.matrix(x)) &&
       (ncol(x)    > 1   || !all(is.numeric(x)))) stop("'x' must be univariate", call.=FALSE)
    x             <- as.vector(x)
    eps           <- stats::sd(x) * sqrt(.Machine$double.eps)
    q             <- NA
    n             <- G
    while(length(q) < (G + 1L))   {
      n           <- n   + 1L
      q           <- unique(stats::quantile(x, seq(from=0L, to=1L, length=n)))
    }
    if(length(q)   > (G + 1L))    {
      q           <- q[-order(diff(q))[seq_len(length(q) - G - 1L)]]
    }
    q[1L]         <- min(x) - eps
    q[length(q)]  <- max(x) + eps
    cl            <- vector("integer", length(x))
    for(i in seq_len(G)) {
      cl[x >= q[i] & x < q[i + 1L]] <- i
    }
      cl
  }

#' Drop constant variables from a formula
#'
#' Drops constant variables from the RHS of a formula taking the data set (\code{dat}), the formula (\code{formula}), and an optional subset vector (\code{sub}) as arguments.
#' @param dat A \code{data.frame} where rows correspond to observations and columns correspond to variables. Ideally column names should be present.
#' @param formula An object of class \code{"\link[stats]{formula}"}: a symbolic description of the model to be fitted. Variables in the \code{formula} not present in the columns of \code{dat} will automatically be discarded. The \code{formula} may include interactions, transformations, or higher order terms: the latter \strong{must} be specified explicitly using the \code{AsIs} operator (\code{\link{I}}).
#' @param sub An optional vector specifying a subset of observations to be used in the fitting process.
#'
#' @return The updated formula with constant variables removed.
#' @note Formulas with and without intercepts are accommodated.
#' @author Keefe Murphy - <\email{keefe.murphy@@ucd.ie}>
#' @keywords utility
#' @seealso \code{\link{drop_levels}}, \code{\link{I}}
#' @export
#' @usage
#' drop_constants(dat,
#'                formula,
#'                sub = NULL)
#' @examples
#' data(ais)
#' hema  <- as.matrix(ais[,3:7])
#' sex   <- ais$sex
#' BMI   <- ais$BMI
#'
#' # Set up a no-intercept regression formula with constant column 'sex'
#' form1 <- as.formula(hema ~ sex + BMI + I(BMI^2) - 1)
#' sub   <- ais$sex == "male"
#'
#' # Try fitting a linear model
#' mod1  <- try(lm(form1, data=ais, subset=sub), silent=TRUE)
#' mod1  <- try(lm(form1, data=ais, subset=sub), silent=TRUE)
#' inherits(mod1, "try-error") # TRUE
#'
#' # Remove redundant variables from formula & try again
#' form2 <- drop_constants(ais, form1, sub)
#' mod2  <- try(lm(form2, data=ais, subset=sub), silent=TRUE)
#' inherits(mod2, "try-error") # FALSE
  drop_constants  <- function(dat, formula, sub = NULL) {
    if(!is.data.frame(dat))                       stop("'dat' must be a data.frame", call.=FALSE)
    Nseq          <- seq_len(nrow(dat))
    sub           <- if(missing(sub)) Nseq else sub
    numsubs       <- all(is.numeric(sub))
    if(!any(numsubs, all(is.logical(sub)) &&
      length(sub) == nrow(dat)))                  stop("'sub' must be a numeric vector, or logical vector with length equal to the number of rows in 'dat'", call.=FALSE)
    if(numsubs    &&
       any(match(sub, Nseq, nomatch = 0)  == 0))  stop("Numeric 'sub' must correspond to row indices of data", call.=FALSE)
    if(!inherits(formula, "formula"))             stop("'formula' must actually be a formula!", call.=FALSE)
    intercept     <- attr(stats::terms(formula), "intercept")
    dat           <- dat[sub,colnames(dat) %in% attr(stats::terms(stats::update.formula(formula, NULL ~ .)), "term.labels"), drop=FALSE]
    ind           <- names(which(!apply(dat, 2L, function(x) all(x == x[1L], na.rm=TRUE))))
    fterms        <- attr(stats::terms(formula), "term.labels")
    ind           <- unique(c(ind[ind %in% fterms], fterms[grepl(":", fterms) | grepl("I\\(", fterms)]))
    response      <- all.vars(stats::update.formula(formula, . ~ NULL))
    form          <- if(length(ind) > 0) stats::reformulate(ind, response=response) else stats::as.formula(paste0(response, " ~ 1"))
    form          <- if(intercept  == 0) stats::update.formula(form, ~ . -1)        else form
    environment(form) <- environment(formula)
      form
  }

#' Drop unused factor levels to predict from unseen data
#'
#' Drops unseen factor levels in \code{newdata} for which predictions are required from a \code{\link[stats]{lm}} model \code{fit}.
#' @param fit A fitted \code{\link[stats]{lm}} model.
#' @param newdata A \code{data.frame} containing variables with which to predict.
#'
#' @return A \code{data.frame} like \code{newdata} with unseen factor levels replaced by \code{NA}.
#' @note This function is untested for models other than \code{\link[stats]{lm}}.
#' @author Keefe Murphy - <\email{keefe.murphy@@ucd.ie}>
#' @keywords utility
#' @seealso \code{\link{drop_constants}}
#' @export
#' @usage
#' drop_levels(fit,
#'             newdata)
#' @examples
#' data(ais)
#' hema  <- as.matrix(ais[,3:7])
#' BMI   <- ais$BMI
#' sport <- ais$sport
#' sub   <- ais$sport != "Row"
#'
#' # Fit a linear model
#' mod   <- lm(hema ~ BMI + sport, data=ais, subset=sub)
#'
#' # Make predictions
#' pred1 <- try(predict(mod, newdata=ais), silent=TRUE)
#' inherits(pred1, "try-error") #TRUE
#'
#' # Remove unused levels and try again
#' pred2 <- try(predict(mod, newdata=drop_levels(mod, ais)), silent=TRUE)
#' inherits(pred2, "try-error") #FALSE
#' anyNA(pred2)                 #TRUE
  drop_levels     <- function(fit, newdata) {
    if(!is.data.frame(newdata))                   stop("'newdata' must be a data.frame", call.=FALSE)
    dat.fac       <- vapply(newdata, is.factor, logical(1L))
    if(!any(dat.fac))                             return(newdata)
    factors       <- rep(names(fit$xlevels), vapply(fit$xlevels, length, integer(1L)))
    factorLevels  <- unname(unlist(fit$xlevels))
    modelFactors  <- cbind.data.frame(factors, factorLevels)
    predictors    <- names(newdata[names(newdata) %in% factors])
    for(i in seq_along(predictors))  {
      ind         <- newdata[,predictors[i]]      %in% modelFactors[modelFactors$factors == predictors[i],]$factorLevels
      if(any(!ind)) {
        newdata[!ind,predictors[i]] <- NA
        newdata[,predictors[i]]     <- factor(newdata[,predictors[i]], levels=modelFactors[modelFactors$factors == predictors[i],]$factorLevels)
      }
    }
      newdata
  }

#' Mahalanobis Distance Outlier Detection for Multivariate Response
#'
#' Computes the Mahalanobis distance between the fitted values and residuals of linear regression models with multivariate or univariate responses.
#' @param fit A fitted \code{\link[stats]{lm}} model, inheriting either the \code{"mlm"} or \code{"lm"} class.
#' @param resids The residuals. Can be residuals for observations included in the model, or residuals arising from predictions on unseen data.
#' @param squared A logical. By default (\code{FALSE}), the generalized interpoint distance is computed. Set this flag to \code{TRUE} for the squared value.
#'
#' @return A vector giving the Mahalanobis distance (or squared Mahalanobis distance) between fitted values and residuals for each observation.
#' @author Keefe Murphy - <\email{keefe.murphy@@ucd.ie}>
#' @importFrom matrixStats "rowSums2"
#' @keywords utility
#' @export
#' @usage
#' MoE_mahala(fit,
#'            resids,
#'            squared = FALSE)
#' @examples
#' data(ais)
#' hema <- as.matrix(ais[,3:7])
#' mod  <- lm(hema ~ sex + BMI, data=ais)
#' res  <- hema - predict(mod)
#' MoE_mahala(mod, res)
  MoE_mahala      <- function(fit, resids, squared = FALSE)    {
    if(length(squared) > 1   ||
       !is.logical(squared))                      stop("'squared' must be a single logical indicator", call.=FALSE)
    if(inherits(fit, "mlm"))  {
      covar       <- stats::estVar(fit)
      covar[!stats::complete.cases((covar))]   <- .Machine$double.eps
      if(diff(dim(resids))   >= 0)     {
        covsvd    <- svd(covar)
        posi      <- covsvd$d > max(sqrt(.Machine$double.eps) * covsvd$d[1L], 0L)
        icov      <- if(all(posi)) covsvd$v   %*% (t(covsvd$u)/covsvd$d) else if(any(posi))
        covsvd$v[,posi, drop=FALSE]   %*% (t(covsvd$u[,posi, drop=FALSE])/covsvd$d[posi]) else array(0L, dim(covar)[2L:1L])
      } else icov <- chol2inv(.chol(covar))
      res         <- rowSums2(resids  %*% icov * resids)
    } else if(inherits(fit, "lm"))     {
      res         <- (resids  * resids)/summary(fit)$sigma
    } else                                        stop("'fit' must be of class 'mlm' or 'lm'", call.=FALSE)
      return(if(isTRUE(squared)) res else sqrt(res))
  }

#' Approximate Hypervolume Estimate
#'
#' Computes simple appproximations to the hypervolume of univariate and multivariate data sets. Also returns the location of the centre of mass.
#' @param data A numeric vector, matrix, or data frame of observations. Categorical variables are not allowed, and covariates should not be included. If a matrix or data frame, rows correspond to observations and columns correspond to variables. There must be more observations than variables.
#' @param method The method used to estimate the hypervolume. The default method uses the function \code{\link[mclust]{hypvol}}. The "\code{convexhull}" and "\code{ellipsoidhull}" options require loading the \code{geometry} and \code{cluster} libraries, respectively. This argument is only relevant for multivariate data; for univariate data, the range of the data is used.
#' @param reciprocal A logical variable indicating whether or not the reciprocal hypervolume is desired rather than the hypervolume itself. The default is to return the hypervolume.
#'
#' @importFrom matrixStats "colMeans2" "colRanges" "rowDiffs" "rowMeans2"
#' @importFrom mclust "hypvol"
#' @note This function is called when adding a noise component to \code{MoEClust} models via the function \code{MoE_control}, specifically it's argument \code{noise.meth}. The function internally only uses the response variables, and not the covariates. However, one can bypass the invocation of this function by specificying its \code{noise.vol} argument directly. This is explicitly necessary for models for high-dimensional data which include a noise component for which this function cannot estimate a (hyper)volume.
#' 
#' Note that supplying the volume manually to \code{\link{MoE_clust}} can affect prediction and by extension the location of the MVN ellipses in \code{\link{MoE_gpairs}} plots. The location cannot be estimated when the volume is supplied manually; in this case, prediction is made on the basis of renormalising the \code{z} matrix after discarding the column corresponding to the noise component.
#' @return A list with the following two elements:
#' \describe{
#' \item{\code{vol}}{A hypervolume estimate (or its inverse). 
#' 
#' This can be used as the hypervolume parameter for the noise component when observations are designated as noise in \code{\link{MoE_clust}}.}
#' \item{\code{loc}}{A vector of length \code{ncol(data)} giving the location of the centre of mass.
#' 
#' This can help in predicting the fitted values of models fitted with noise components via \code{\link{MoE_clust}}.}}
#' @author Keefe Murphy - <\email{keefe.murphy@@ucd.ie}>
#' @seealso \code{\link[mclust]{hypvol}}, \code{\link[geometry]{convhulln}}, \code{\link[cluster]{ellipsoidhull}}
#' @keywords control
#' @export
#' @usage
#' noise_vol(data,
#'           method = c("hypvol", "convexhull", "ellipsoidhull"),
#'           reciprocal = FALSE)
#' @examples
#' data(ais)
#' noise_vol(ais[,3:7], reciprocal=TRUE)
#' 
#' noise_vol(ais[,3:7], reciprocal=FALSE, method="convexhull")
  noise_vol       <- function(data, method = c("hypvol", "convexhull", "ellipsoidhull"), reciprocal = FALSE) {
    data          <- as.matrix(data)
    method        <- match.arg(method)
    has.lib       <- switch(EXPR=method, hypvol=TRUE, convexhull=suppressMessages(requireNamespace("geometry", quietly=TRUE)), ellipsoidhull=suppressMessages(requireNamespace("cluster", quietly=TRUE)))
    if(!has.lib)                                  stop(paste0("Use of the ", method, " 'method' option requires loading the ", switch(EXPR=method, hypvol="'mclust'", convexhull="'geometry'", ellipsoidhull="'cluster'"), "library"), call.=FALSE)
    if(length(reciprocal) != 1 ||
       !is.logical(reciprocal))                   stop("'reciprocal' must be a single logical indicator", call.=FALSE)
    if(ncol(data) == 1)    {
      vol         <- ifelse(reciprocal, 1/abs(diff(range(data))), abs(diff(range(data))))
      loc         <- vol/2
    } else         {
      switch(EXPR=method, 
             hypvol=       {
        bdvlog    <- .SLDC(data)
        PCA       <- stats::princomp(data)
        scores    <- PCA$scores
        pcvlog    <- .SLDC(scores)
        vlog      <- min(bdvlog, pcvlog)
        if(reciprocal)     {
          minlog  <- log(.Machine$double.xmin)
          if(-vlog < minlog) {                    warning("hypervolume smaller than smallest machine representable positive number", call.=FALSE, immediate.=TRUE)
            vol   <- 0L
          } else vol      <- exp(-vlog)
        }   else   {
          maxlog  <- log(.Machine$double.xmax)
          if(vlog  > maxlog) {                    warning("hypervolume greater than largest machine representable number", call.=FALSE, immediate.=TRUE)
            vol   <- Inf
          } else vol      <- exp(vlog)
        } 
        loc       <- if(bdvlog <= pcvlog) rowMeans2(colRanges(data)) else colMeans2(data) + tcrossprod(rowMeans2(colRanges(scores)), PCA$loadings)
      }, ellipsoidhull=    {
        hull      <- cluster::ellipsoidhull(data)
        vol       <- ifelse(reciprocal, 1/.vol_ellipsoid(hull), .vol_ellipsoid(hull))
        loc       <- hull$loc
      }, convexhull=       {
        hull      <- geometry::convhulln(data, options=c("Pp", "FA", "Fx"), output.options=TRUE)
        vol       <- ifelse(reciprocal, 1/hull$vol, hull$vol)
        loc       <- colMeans2(hull$p)
      })
    }
    attr(vol, "Inverse")  <- reciprocal
      return(list(vol=vol, loc=loc))
  }

#' Show the NEWS file
#'
#' Show the \code{NEWS} file of the \code{MoEClust} package.
#' @return The \code{MoEClust} \code{NEWS} file, provided the session is interactive.
#' @export
#' @keywords utility
#'
#' @usage MoE_news()
#' @examples
#' MoE_news()
  MoE_news   <- function() {
    newsfile <- file.path(system.file(package  = "MoEClust"), "NEWS.md")
       if(interactive()) file.show(newsfile) else message("The session is not interactive\n")
  }

# Hidden/Print/Summary Functions
  .chol           <- function(x) tryCatch(chol(x), error=function(e) {
    d             <- nrow(x)
    eigs          <- eigen(x, symmetric = TRUE)
    eval          <- eigs$values
    evec          <- eigs$vectors
     return(chol(x + evec %*% tcrossprod(diag(pmax.int(0L, 2 * max(abs(eval)) * d * .Machine$double.eps - eval), d), evec)))
    }
  )

  .crits_names    <- function(x) {
      unlist(lapply(seq_along(x), function(i) stats::setNames(x[[i]], paste0(names(x[i]), "|", names(x[[i]])))))
  }

  .mat_byrow      <- function(x, nrow, ncol) {
      matrix(x, nrow=nrow, ncol=ncol, byrow=any(dim(as.matrix(x)) == 1))
  }

  .mat_sq         <- function(x, d) diag(x, d)

  .npars_clustMD  <- function(D, G, J, CnsIndx, OrdIndx, K) {
    G1            <- G  - 1L
    GD            <- G  * D
    OC            <- OrdIndx   - CnsIndx
    if(Jind       <- J  > OrdIndx)  {
      O1          <- OrdIndx   - 1L
      Dord        <- D  - OrdIndx
      G1dO        <- G1 * Dord
      Gord        <- G  * OrdIndx
      GO1dO       <- G  + Gord + G1dO
      npars       <- c(EII=GO1dO,
                       VII=2L  * G1 + GO1dO,
                       EEI=O1  + GO1dO,
                       VEI=O1  + 2L * G1  + GO1dO,
                       EVI=O1  * G  - G1  + G1dO + GO1dO,
                       VVI=G1  + G  * O1  + G1dO + GO1dO)
    } else         {
      D1          <- D  - 1L
      GGD         <- GD + G
      npars       <- c(EII=GGD,
                       VII=GGD + G1,
                       EEI=D   + G1 + GD,
                       VEI=GGD + D1 + G1,
                       EVI=G   * D1 + GGD,
                       VVI=G1  + 2L * GD)
    }
      c(npars, BD=G * CnsIndx  * (CnsIndx + 1L)/2L + G1 + GD +
        ifelse(J > CnsIndx, G  * OC * (OC + 1L)/2L, 0L) +
        ifelse(Jind,  G * sum(K/2L  * (K  - 1L)),   0L))
  }

  .pick_MoECrit   <- function(x, pick = 3L) {
    if(!inherits(x, "MoECriterion"))              stop("'x' must be an object of class 'MoECriterion'", call.=FALSE)
    x             <- replace(x, !is.finite(x), NA)
    pick          <- min(pick,        length(x[!is.na(x)]))
    decrease      <- !is.element(attr(x, "criterion"), c("DF", "ITERS"))
    x.sx          <- sort(x,          decreasing=decrease)[pick]
    x.crit        <- if(decrease)     x   >= x.sx else x <= x.sx
    x.ind         <- which(x.crit,    arr.ind=TRUE)
    x.val         <- sort(x[x.ind],   decreasing=decrease)
    ind.x         <- order(x[x.ind],  decreasing=decrease)
    x.ind         <- x.ind[ind.x,,    drop=FALSE]
    x.ind[,1L]    <- gsub(".*= ", "", rownames(x)[x.ind[,1L]])
    x.ind[,2L]    <- colnames(x)[as.numeric(x.ind[,2L])]
      return(list(crits = stats::setNames(x.val[seq_len(pick)], vapply(seq_len(pick), function(p, b=x.ind[p,]) paste0(b[2L], ",", b[1L]), character(1L))), pick = pick))
  }
  
  #' @importFrom matrixStats "rowSums2"
  .renorm_z       <- function(z) z/rowSums2(z)
  
  #' @importFrom matrixStats "colRanges" "rowDiffs"
  .SLDC           <- function(x) sum(log(abs(rowDiffs(colRanges(x)))))

  .sq_mat         <- function(x) diag(sqrt(diag(x)))

  .tau_noise      <- function(tau, z0) {
    t0            <- mean(z0)
      cbind(tau * (1 - t0), unname(t0))
  }
  
  .trim_out       <- function(x, signif = 0.01, na.rm = TRUE, replace = TRUE, ...) {
    qnt           <- stats::quantile(x, probs=c(signif, 2 - signif)/2, na.rm=na.rm, ...)
    H             <- 1.5     * stats::IQR(x, na.rm=na.rm)
    y             <- x
    li.qnt        <- qnt[1L] - H
    ui.qnt        <- qnt[2L] + H
    y[x < li.qnt] <- ifelse(replace, li.qnt, NA)
    y[x > ui.qnt] <- ifelse(replace, ui.qnt, NA)
      y
  }
  
  .vol_ellipsoid  <- function(x)  {
   exp(ifelse((p2 <- length(x$loc)/2) > 1,
              p2   * log(x$d2) + determinant(x$cov)$modulus/2 + p2 * log(base::pi) - lgamma(p2) - log(p2),
              log(base::pi)    + determinant(x$cov)$modulus/2 + log(x$d2))[1L])
  }

#' @method print MoEClust
#' @importFrom mclust "mclustModelNames"
#' @rdname MoE_clust
#' @usage
#' \method{print}{MoEClust}(x,
#'       digits = 3L,
#'       ...)
#' @export
  print.MoEClust  <- function(x, digits = 3L, ...) {
    cat("Call:\t");  print(x$call)
    if(length(digits)  > 1 || !is.numeric(digits) ||
       digits     <= 0)                           stop("Invalid 'digits'", call.=FALSE)
    name          <- x$modelName
    G             <- x$G
    gating        <- attr(x$gating, "Formula")
    expert        <- attr(x$expert, "Formula")
    gate.x        <- !attr(x, "Gating")
    exp.x         <- !attr(x, "Expert")
    net.x         <- !c(gate.x, exp.x)
    crit          <- round(unname(c(x$bic, x$icl, x$aic)), digits)
    hypvol        <- x$hypvol
    noise         <- !is.na(hypvol)
    equalP        <- G <= 1  || attr(x$gating, "EqualPro")
    equalN        <- noise   && attr(x$gating, "EqualNoise") && equalP
    cat(paste0("\nBest Model", ifelse(length(x$BIC) > 1, paste0(" (according to ", toupper(attr(x, "Criterion")), "): "), ": "), mclustModelNames(name)$type, " (", name, "), with ",
               ifelse(G == 0, "only a noise component", paste0(G, " component", ifelse(G   > 1, "s",   ""))),
               ifelse(G == 0 || !noise,   "\n", " (and a noise component)\n"),
               ifelse(!equalP |
                      G <= 1,  "",   paste0("Equal Mixing Proportions", ifelse(equalN | G <= 1 | !noise, "\n", " (with estimated noise component mixing proportion)\n"))),
               ifelse(!noise,  "",   paste0("Hypervolume of Noise Component: ", round(hypvol, digits), "\n")),
               "BIC = ", crit[1L], " | ICL = ", crit[2L], " | AIC = ",  crit[3L],
               ifelse(any(net.x),    paste0("\nIncluding",  ifelse(all(net.x), " gating and expert", ifelse(!gate.x, " gating", ifelse(!exp.x, " expert", ""))), " network covariates:\n"), "\nNo covariates\n"),
               ifelse(gate.x,  "",   paste0("\tGating: ",  gating,  ifelse(exp.x,  "", "\n"))),
               ifelse(exp.x,   "",   paste0("\tExpert: ",  expert, ""))))
      invisible()
  }

#' @method summary MoEClust
#' @rdname MoE_clust
#' @usage
#' \method{summary}{MoEClust}(object,
#'         ...)
#' @export
  summary.MoEClust        <- function(object, ...) {
    G             <- object$G
    attr(G, "range")      <- eval(object$call$G)
    params        <- object$parameters
    hypvol        <- object$hypvol
    equalPro      <- G <= 1 || attr(object$gating, "EqualPro")
    equalN        <- !is.na(hypvol) && attr(object$gating, "EqualNoise") && equalPro
    summ          <- list(title = "Gaussian Parsimonious Clustering Model with Covariates", data = deparse(object$call$data), n = object$n, d = object$d, G = G, modelName = object$modelName, algo=attr(object, "Algo"),
                          loglik = object$loglik[length(object$loglik)], df = object$df, iters = object$iters, gating = object$gating, expert = object$expert, bic=unname(object$bic), icl = unname(object$icl), aic = unname(object$aic),
                          pro = params$pro, mean = params$mean, variance = params$variance$sigma, Vinv = params$Vinv, hypvol = hypvol, z = object$z, equalPro = equalPro, equalNoise = equalN, classification = object$classification)
    class(summ)   <- "summary_MoEClust"
      summ
 }

#' @method print summary_MoEClust
#' @importFrom mclust "mclustModelNames"
#' @export
  print.summary_MoEClust  <- function(x, digits = 3L, ...) {
    if(length(digits)  > 1 || !is.numeric(digits) ||
       digits     <= 0)                           stop("Invalid 'digits'", call.=FALSE)
    tmp           <- data.frame(log.likelihood = round(x$loglik, digits), n = x$n, d = x$d, df = x$df, iters = x$iters,
                                BIC = round(x$bic, digits), ICL = round(x$icl, digits), AIC = round(x$aic, digits))
    tmp           <- if(is.na(x$hypvol))   tmp    else cbind(tmp, HypVol = x$hypvol)
    tmp           <- cbind(tmp, Algo = x$algo)
    rownames(tmp) <- NULL
    name          <- x$modelName
    G             <- x$G
    range.G       <- attr(G, "range")
    if(G          == min(range.G))                message("Best model occurs at the min of the number of components considered\n")
    if(G          == max(range.G))                message("Best model occurs at the max of the number of components considered\n")
    noise         <- !is.na(x$hypvol)
    gating        <- attr(x$gating, "Formula")
    expert        <- attr(x$expert, "Formula")
    gate.x        <- gating == "~1"
    exp.x         <- expert == "~1"
    equalP        <- x$equalPro && gate.x
    equalN        <- noise  && x$equalNoise && equalP
    zs            <- table(x$classification)
    cat(paste0("------------------------------------------------------\n", x$title, "\nData: ",
               x$data,"\n", "------------------------------------------------------\n\n",
               "MoEClust ",   name, " (", mclustModelNames(name)$type, "), with ",
               ifelse(G == 0, "only a noise component", paste0(G, " component", ifelse(G > 1, "s", ""))),
               ifelse(G == 0  || is.na(x$hypvol), "\n", " (and a noise component)\n"),
               paste0("\nGating Network Covariates:  ", ifelse(gate.x, "None", gating)),
               paste0("\nExpert Network Covariates:  ", ifelse(exp.x,  "None", expert)),
               ifelse(G  > 1, paste0("\nEqual Mixing Proportions:   ",  equalP), ""),
               paste0("\nNoise Component:            ", noise, "\n\n"),
               ifelse(G  > 1  && noise, paste0("\nNoise Proportion Estimated: ", !equalN), "")))
    print(tmp, row.names = FALSE)
    cat("\nClustering table:")
    print(zs,  row.names = FALSE)
      invisible()
  }

#' @method print MoECompare
#' @rdname MoE_compare
#' @usage
#' \method{print}{MoECompare}(x,
#'       index = seq_len(x$pick),
#'       digits = 3L,
#'       ...)
#' @export
  print.MoECompare       <- function(x, index=seq_len(x$pick), digits = 3L, ...) {
    index                <- if(is.logical(index)) which(index) else index
    if(length(index) < 1 || (!is.numeric(index) &&
       (any(index    < 1  | index > x$pick))))    stop("Invalid 'index'",  call.=FALSE)
    if(length(digits)     > 1    ||
       !is.numeric(digits)       ||
       digits            <= 0)                    stop("Invalid 'digits'", call.=FALSE)
    n.all                <- all(is.na(x$hypvol))
    x$hypvol             <- NULL
    x$noise              <- if(n.all)                   NULL else x$noise
    noise.gate           <- if(n.all)                   NULL else replace(x$noise.gate, is.na(x$noise.gate), "")
    x$noise.gate         <- NULL
    x$noise.gate         <- if(all(x$gating == "None")) NULL else noise.gate
    equalPro             <- if(all(is.na(x$equalPro)))  NULL else replace(x$equalPro,   is.na(x$equalPro),   "")
    x$equalPro           <- NULL
    x$equalPro           <- equalPro
    na.equalNoise        <- is.na(x$equalNoise)
    equalNoise           <- replace(x$equalNoise, na.equalNoise,    "")
    x$equalNoise         <- NULL
    x$equalNoise         <- if(all(na.equalNoise))      NULL else equalNoise
    x$bic                <- round(x$bic,    digits)
    x$icl                <- round(x$icl,    digits)
    x$aic                <- round(x$aic,    digits)
    x$loglik             <- round(x$loglik, digits)
    cat(paste0("------------------------------------------------------------------------------\n", 
               x$title, "\nData: ", x$data, "\nRanking Criterion: ", toupper(attr(x, "Crit")), "\nOptimal Only: ", attr(x, "Opt"), 
              "\n------------------------------------------------------------------------------\n\n"))
    comp.res             <- data.frame(do.call(cbind, x[-seq_len(4L)]))[index,, drop=FALSE]
    comp.res             <- comp.res[,!vapply(comp.res, function(x) all(x == ""), logical(1L)), drop=FALSE]
    comp.res             <- cbind(rank = rownames(comp.res), comp.res)
    rownames(comp.res)   <- NULL
    print(comp.res, row.names = FALSE)
      invisible()
  }

#' @method print MoECriterion
#' @export
  print.MoECriterion     <- function(x, pick = 3L, ...) {
    if(length(pick)      != 1    ||
       !is.numeric(pick))                         stop("'pick' must be a single number", call.=FALSE)
    if(floor(pick)       != pick ||
       pick        < 1)                           stop("'pick' be a strictly positive integer", call.=FALSE)
    algo          <- attr(x, "algo")
    crit          <- attr(x, "criterion")
    choice        <- .pick_MoECrit(x, pick)
    pick          <- choice$pick
    dim1          <- attr(x, "dim")
    dim2          <- attr(x, "dimnames")
    attributes(x) <- NULL
    attr(x, "dim")       <- dim1
    attr(x, "dimnames")  <- dim2
    cat(switch(EXPR= crit, BIC="Bayesian Information Criterion (BIC):\n",       ICL="Integrated Completed Likelihood (ICL):\n",
                           AIC="Akaike Information Criterion (AIC):\n",          DF="Number of Estimated Parameters (Residual DF):\n",
                         ITERS=paste0("Number of ", algo, " Iterations:\n"), loglik="Maximal Log-Likelihood:\n"))
    print(unclass(x))
    cat(paste0("\nTop ", ifelse(pick > 1, paste0(pick, " models"), "model"), " based on the ", crit, " criterion:\n"))
    print(choice$crits)
      invisible()
  }

#' @method print MoE_gating
#' @export
  print.MoE_gating       <- function(x, ...) {
    noise         <- attr(x, "Noise")
    equalpro      <- attr(x, "EqualPro")
    formula       <- attr(x, "Formula")
    equalNoise    <- noise && equalpro
    gateNoise     <- noise && !equalpro && formula != "~1"
    class(x)      <- class(x)[class(x) != "MoE_gating"]
    print(x, ...)
    cat(paste("Formula:",  formula,  "\n"))
    cat(paste("Noise:",    noise,    "\n"))
    if(gateNoise)    cat(paste("Noise Component Gating:", attr(x, "NoiseGate"), "\n"))
    cat(paste("EqualPro:", equalpro, ifelse(equalNoise, "\n", "")))
    if(equalNoise)   cat(paste("Noise Proportion Estimated:", attr(x, "EqualNoise")))
    message("\n\n\nUsers are cautioned against making inferences about statistical significance from summaries of the coefficients in the gating network\n")
      invisible(x)
  }

#' @method print MoE_expert
#' @export
  print.MoE_expert       <- function(x, ...) {
    if(is.na(x)   || (names(x) == "Cluster0"))    stop("No expert network exists for models with only a noise component", call.=FALSE)
    formula       <- attr(x, "Formula")
    attributes(x)[-1L]   <- NULL
    class(x)      <- "listof"
    print(x, ...)
    cat(paste("Formula:", formula))
    message("\n\n\nUsers are cautioned against making inferences about statistical significance from summaries of the coefficients in the expert network\n")
      invisible(x)
  }

#' @method summary MoE_gating
#' @export
  summary.MoE_gating     <- function(object, ...) {
    noise         <- attr(object, "Noise")
    noise.gate    <- attr(object, "NoiseGate")
    equalpro      <- attr(object, "EqualPro")
    equalnoise    <- attr(object, "EqualNoise")
    formula       <- attr(object, "Formula")
    class(object) <- "multinom"
    summ          <- summary(object, ...)
    class(summ)   <- "summary_MoEgate"
    attr(summ, "Noise")  <- noise
    attr(summ, "NoiseG") <- noise.gate
    attr(summ, "EqualP") <- equalpro
    attr(summ, "EqualN") <- equalnoise
    attr(summ, "Form")   <- formula
      summ
  }

#' @method summary MoE_expert
#' @export
  summary.MoE_expert     <- function(object, clusters = seq_along(object), ...) {
    if(length(object)    == 1 &&
      (is.na(object)     || 
      (names(object)     == "Cluster0")))        {
      summ               <- stats::setNames(NA, "Cluster0")
    } else                {
      if(any(!is.numeric(clusters), any(clusters < 1),
             any(clusters > length(object))))     stop("Invalid 'clusters'", call.=FALSE)
      summ        <- lapply(object[clusters], summary, ...)
    }
    class(summ)   <- "summary_MoEexp"
    attr(summ, "Form")   <- attr(object, "Formula") 
      summ
  }

#' @method print summary_MoEgate
#' @export
  print.summary_MoEgate  <- function(x, ...) {
    formula       <- attr(x, "Form")
    noise         <- attr(x, "Noise")
    equalpro      <- attr(x, "EqualP")
    equalNoise    <- noise && equalpro
    gateNoise     <- noise && !equalpro && formula != "~1"
    class(x)      <- "summary.multinom"
    print(x, ...)
    cat(paste("Formula:",  formula,  "\n"))
    cat(paste("Noise:",    noise,    "\n"))
    if(gateNoise)    cat(paste("Noise Component Gating:", attr(x, "NoiseG"), "\n"))
    cat(paste("EqualPro:", equalpro, ifelse(equalNoise, "\n", "")))
    if(equalNoise)   cat(paste("Noise Proportion Estimated:", attr(x, "EqualN")))
    message("\n\n\nUsers are cautioned against making inferences about statistical significance from summaries of the coefficients in the gating network\n")
      invisible(x)
  }

#' @method print summary_MoEexp
#' @export
  print.summary_MoEexp   <- function(x, ...) {
    if(is.na(x)   || (names(x) == "Cluster0"))    stop("No expert network exists for models with only a noise component", call.=FALSE)
    class(x)      <- "listof"
    print(x, ...)
    cat(paste("Formula:", attr(x, "Form"), "\n"))
    message("\n\n\nUsers are cautioned against making inferences about statistical significance from summaries of the coefficients in the expert network\n")
      invisible(x)
  }
