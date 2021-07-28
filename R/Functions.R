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
#' For zero-component models with a noise component only the \code{"E"} and \code{"EII"} models will be fitted for univariate and multivariate data, respectively, although this is clearly for naming consistency only. The help file for \code{\link[mclust]{mclustModelNames}} further describes the available models (though the \code{"X"} in the single-component models will be coerced to \code{"E"} if supplied that way). For single-component models, other model names equivalent to those above can be supplied, but will be coerced to those above.
#' @param gating A \code{\link[stats]{formula}} for determining the model matrix for the multinomial logistic regression in the gating network when fixed covariates enter the mixing proportions. Defaults to \code{~1}, i.e. no covariates. This will be ignored where \code{G=1}. Continuous, categorical, and/or ordinal covariates are allowed. Logical covariates will be coerced to factors. Interactions, transformations, and higher order terms are permitted: the latter \strong{must} be specified explicitly using the \code{AsIs} operator (\code{\link{I}}). The specification of the LHS of the formula is ignored. Intercept terms are included by default.
#' @param expert A \code{\link[stats]{formula}} for determining the model matrix for the (multivariate) WLS in the expert network when fixed covariates are included in the component densities. Defaults to \code{~1}, i.e. no covariates. Continuous, categorical, and/or ordinal covariates are allowed. Logical covariates will be coerced to factors. Interactions, transformations, and higher order terms are permitted: the latter \strong{must} be specified explicitly using the \code{AsIs} operator (\code{\link{I}}). The specification of the LHS of the formula is ignored. Intercept terms are included by default.
#' @param control A list of control parameters for the EM/CEM and other aspects of the algorithm. The defaults are set by a call to \code{\link{MoE_control}}. In particular, arguments pertaining to the inclusion of an additional noise component are documented here.
#' @param network.data An optional data frame (or a matrix with named columns) in which to look for the covariates in the \code{gating} &/or \code{expert} network formulas, if any. If not found in \code{network.data}, any supplied \code{gating} &/or \code{expert} covariates are taken from the environment from which \code{MoE_clust} is called. Try to ensure the names of variables in \code{network.data} do not match any of those in \code{data}.
#' @param ... An alternative means of passing control parameters directly via the named arguments of \code{\link{MoE_control}}. Do not pass the output from a call to \code{\link{MoE_control}} here! This argument is only relevant for the \code{\link{MoE_clust}} function and will be ignored for the associated \code{print} and \code{summary} functions.
#' @param x,object,digits,classification,parameters,networks Arguments required for the \code{print} and \code{summary} functions: \code{x} and \code{object} are objects of class \code{"MoEClust"} resulting from a call to \code{\link{MoE_clust}}, while \code{digits} gives the number of decimal places to round to for printing purposes (defaults to 3). \code{classification}, \code{parameters}, and \code{networks} are logicals which govern whether a table of the MAP classification of observations, the mixture component parameters, and the gating/expert network coefficients are printed, respectively.

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
#' \item{\code{G}}{The optimal number of mixture components according to \code{criterion}.}
#' \item{\code{BIC}}{A matrix of \emph{all} BIC values with \code{length{G}} rows and \code{length(modelNames)} columns. May include missing entries: \code{NA} represents models which were not visited, \code{-Inf} represents models which were terminated due to error, for which a log-likelihood could not be estimated. Inherits the classes \code{"MoECriterion"} and \code{"mclustBIC"}, for which dedicated printing and plotting functions exist, respectively.}
#' \item{\code{ICL}}{A matrix of \emph{all} ICL values with \code{length{G}} rows and \code{length(modelNames)} columns. May include missing entries: \code{NA} represents models which were not visited, \code{-Inf} represents models which were terminated due to error, for which a log-likelihood could not be estimated. Inherits the classes \code{"MoECriterion"} and \code{"mclustICL"}, for which dedicated printing and plotting functions exist, respectively.}
#' \item{\code{AIC}}{A matrix of \emph{all} AIC values with \code{length{G}} rows and \code{length(modelNames)} columns. May include missing entries: \code{NA} represents models which were not visited, \code{-Inf} represents models which were terminated due to error, for which a log-likelihood could not be estimated. Inherits the classes \code{"MoECriterion"} and \code{"mclustAIC"}, for which dedicated printing and plotting functions exist, respectively.}
#' \item{\code{bic}}{The BIC value corresponding to the optimal model. May not necessarily be the optimal BIC.}
#' \item{\code{icl}}{The ICL value corresponding to the optimal model. May not necessarily be the optimal ICL.}
#' \item{\code{aic}}{The AIC value corresponding to the optimal model. May not necessarily be the optimal AIC.}
#' \item{\code{gating}}{An object of class \code{"MoE_gating"} (for which dedicated \code{print}, \code{summary}, and \code{\link[=predict.MoE_gating]{predict}} methods exist) and either \code{"multinom"} or \code{"glm"} (only for single-component models or noise-only models) giving the \code{\link[nnet]{multinom}} regression coefficients of the \code{gating} network. If \code{gating} covariates were \emph{NOT} supplied (or the best model has just one component), this corresponds to a RHS of \code{~1}, otherwise the supplied \code{gating} formula. As such, a fitted \code{gating} network is always returned even in the absence of supplied covariates or clusters. The number of parameters to penalise by for \code{\link{MoE_crit}} is given by \code{length(coef(gating))}, and the \code{gating} formula used is stored here as an attribute. If there is a noise component (and the option \code{noise.gate=TRUE} is invoked), its coefficients are those for the \emph{last} component. \strong{Users are cautioned against making inferences about statistical significance from summaries of the coefficients in the gating network}.}
#' \item{\code{expert}}{An object of class \code{"MoE_expert"} (for which dedicated \code{print}, \code{summary}, and \code{\link[=predict.MoE_expert]{predict}} methods exist) and \code{"lm"} giving the (multivariate) WLS regression coefficients of the \code{expert} network. If \code{expert} covariates were NOT supplied, this corresponds to a RHS of \code{~1}, otherwise the supplied \code{expert} formula. As such, a fitted \code{expert} network is always returned even in the absence of supplied covariates. The number of parameters to penalise by for \code{\link{MoE_crit}} is given by \code{G * length(coef(expert[[1]]))}, and the \code{expert} formula used is stored here is an attribute. \strong{Users are cautioned against making inferences about statistical significance from summaries of the coefficients in the expert network}.}
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
#' For models with expert network covariates, this is given by the posterior mean of the fitted values, otherwise the posterior mean of the response is reported. For models with expert network covariates, the \emph{observation-specific} component means can be accessed by calling \code{\link[=predict.MoE_expert]{predict}} on the \code{expert} object above.}
#' \item{\code{variance}}{A list of variance parameters of each component of the model. The components of this list depend on the model type specification. See the help file for \code{\link[mclust]{mclustVariance}} for details. Also see \code{\link{expert_covar}} for an alternative approach to summarising the variance parameters in the presence of expert network covariates.}
#' \item{\code{Vinv}}{The inverse of the hypervolume parameter for the noise component if required, otherwise set to \code{NULL} (see \code{\link{MoE_control}}).}
#' }}
#' \item{\code{z}}{The final responsibility matrix whose \code{[i,k]}-th entry is the probability that observation \emph{i} belonds to the \emph{k}-th component. If there is a noise component, its values are found in the \emph{last} column.}
#' \item{\code{classification}}{The vector of cluster labels for the chosen model corresponding to \code{z}, i.e. \code{max.col(z)}. Observations belonging to the noise component, if any, will belong to component \code{0}.}
#' \item{\code{uncertainty}}{The uncertainty associated with the \code{classification}.}
#' \item{\code{net.covs}}{A data frame gathering the unique set of covariates used in the \code{gating} and \code{expert} networks, if any. Will contain zero columns in the absence of gating or expert network covariates. Supplied gating covariates will be excluded if the optimal model has only one component. May have fewer columns than covariates supplied via the \code{network.data} argument also, as only the included covariates are gathered here.}
#' \item{\code{resid.data}}{In the presence of expert network covariates, this is the augmented data actually used in the clustering at convergence, as a list of \code{G} matrices of WLS residuals of dimension \code{n * d}. Will contain zero columns in the absence of expert network covariates.}
#' \item{\code{DF}}{A matrix giving the numbers of estimated parameters (i.e. the number of 'used' degrees of freedom) for \emph{all} visited models, with \code{length{G}} rows and \code{length(modelNames)} columns. Subtract these numbers from \code{n} to get the degrees of freedom. May include missing entries: \code{NA} represents models which were not visited, \code{-Inf} represents models which were terminated due to error, for which parameters could not be estimated. Inherits the classes \code{"MoECriterion"} and \code{"mclustDF"}, for which dedicated printing and plotting functions exist, respectively.}
#' \item{\code{ITERS}}{A matrix giving the total number of EM/CEM iterations for \emph{all} visited models, with \code{length{G}} rows and \code{length(modelNames)} columns. May include missing entries: \code{NA} represents models which were not visited, \code{Inf} represents models which were terminated due to singularity/error and thus would never have converged. Inherits the classes \code{"MoECriterion"} and \code{"mclustITERS"}, for which dedicated printing and plotting functions exist, respectively.}
#' Dedicated \code{\link[=plot.MoEClust]{plot}}, \code{\link[=predict.MoEClust]{predict}}, \code{print}, and \code{summary} functions exist for objects of class \code{"MoEClust"}. The results can be coerced to the \code{"Mclust"} class to access other functions from the \pkg{mclust} package via \code{\link[=as.Mclust.MoEClust]{as.Mclust}}.
#' @details The function effectively allows 6 different types of Gaussian Mixture of Experts model (as well as the different models in the GPCM/\pkg{mclust} family, for each): i) the standard finite Gaussian mixture with no covariates, ii) fixed covariates only in the gating network, iii) fixed covariates only in the expert network, iv) the full Mixture of Experts model with fixed covariates entering both the mixing proportions and component densities. By constraining the mixing proportions to be equal (see \code{equalPro} in \code{\link{MoE_control}}) two extra special cases are facilitated when gating covariates are excluded. 
#' 
#' Note that having the same covariates in both networks is allowed. So too are interactions, transformations, and higher order terms (see \code{\link[stats]{formula}}): the latter \strong{must} be specified explicitly using the \code{AsIs} operator (\code{\link{I}}). Covariates can be continuous, categorical, logical, or ordinal, but the response must always be continuous.
#'
#' While model selection in terms of choosing the optimal number of components and the GPCM/\pkg{mclust} model type is performed within \code{\link{MoE_clust}}, using one of the \code{criterion} options within \code{\link{MoE_control}}, choosing between multiple fits with different combinations of covariates or different initialisation settings can be done by supplying objects of class \code{"MoEClust"} to \code{\link{MoE_compare}}.
#' @note Where \code{BIC}, \code{ICL}, \code{AIC}, \code{LOGLIK}, \code{DF} and \code{ITERS} contain \code{NA} entries, this corresponds to a model which was not run; for instance a VVV model is never run for single-component models as it is equivalent to EEE. As such, one can consider the value as not really missing, but equivalent to the EEE value. \code{BIC}, \code{ICL}, \code{AIC}, \code{LOGLIK}, \code{DF} and \code{ITERS} all inherit the classes \code{"MoECriterion"} and \code{"mclustBIC", "mclustICL", etc.}, for which dedicated printing and plotting functions exist, respectively.
#'
#' @seealso See \code{\link{MoE_stepwise}} for identifying the optimal model and its covariates via greedy forward stepwise selection.\cr
#' 
#' \code{\link{MoE_control}}, \code{\link{MoE_compare}}, \code{\link{plot.MoEClust}}, \code{\link{predict.MoEClust}}, \code{\link{predict.MoE_gating}}, \code{\link{predict.MoE_expert}}, \code{\link[=as.Mclust.MoEClust]{as.Mclust}}, \code{\link{MoE_crit}}, \code{\link{MoE_estep}}, \code{\link{MoE_cstep}}, \code{\link{MoE_dens}}, \code{\link[mclust]{mclustModelNames}}, \code{\link[mclust]{mclustVariance}}, \code{\link{expert_covar}}, \code{\link{aitken}}, \code{\link{I}}
#' @export
#' @references Murphy, K. and Murphy, T. B. (2020). Gaussian parsimonious clustering models with covariates and a noise component. \emph{Advances in Data Analysis and Classification}, 14(2): 293-325. <\doi{10.1007/s11634-019-00373-8}>.
#'
#' Fraley, C. and Raftery, A. E. (2002). Model-based clustering, discriminant analysis, and density estimation. \emph{Journal of the American Statistical Association}, 97(458): 611-631.
#' @author Keefe Murphy - <\email{keefe.murphy@@mu.ie}>
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
#' \donttest{data(ais)
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
#' m5    <- MoE_clust(hema, G=2:3, expert= ~ sex + BMI, equalPro=TRUE)
#' 
#' # Fit models with gating covariates & an additional noise component
#' m6    <- MoE_clust(hema, G=2:3, tau0=0.1, gating= ~ BMI, network.data=ais)
#'
#' # Extract the model with highest BIC
#' (comp <- MoE_compare(m1, m2, m3, m4, m5, m6, criterion="bic"))
#'  
#' # See if a better model can be found using greedy forward stepwise selection
#' (step <- MoE_stepwise(ais[,3:7], ais))
#' (comp <- MoE_compare(comp, step, optimal.only=TRUE))
#' (best <- comp$optimal)
#' (summ <- summary(best, classification=TRUE, parameters=TRUE, networks=TRUE))
#'
#' # Examine the expert network in greater detail
#' # (but refrain from inferring statistical significance!)
#' summary(best$expert)
#'
#' # Visualise the results, incl. the gating network and log-likelihood
#' plot(best, what="gpairs")
#' plot(best, what="gating") # equal mixing proportions!
#' plot(best, what="loglik")
#'
#' # Visualise the results using the 'lattice' library
#' require("lattice")
#' z     <- factor(best$classification, labels=paste0("Cluster", seq_len(best$G)))
#' splom(~ hema | sex, groups=z)
#' splom(~ hema | z, groups=sex)}
  MoE_clust       <- function(data, G = 1:9, modelNames = NULL, gating = ~1, expert = ~1, control = MoE_control(...), network.data = NULL, ...) {

  # Definitions and storage set-up
    call          <- call2     <- match.call()
    multi         <- is.null(modelNames)
    gate.x        <- !missing(gating)
    exp.x         <- !missing(expert)
    dots          <- list(...)
    if(!missing(control)       &&
      length(dots[names(dots) %in%
                  names(control)])     > 0)       stop("Arguments cannot be supplied via the '...' construct when the named argument 'control' is supplied", call.=FALSE) 
    criterion     <- control$criterion
    stopaX        <- control$stopping == "aitken"
    init.z        <- control$init.z
    z.list        <- control$z.list
    nstarts       <- switch(EXPR=init.z, random=control$nstarts, 1L)
    startseq      <- seq_len(nstarts)
    multstart     <- nstarts > 1
    algo          <- control$algo
    exp.init      <- control$exp.init
    do.joint      <- exp.init$joint   && init.z != "random"
    max.init      <- exp.init$max.init
    Identity      <- exp.init$identity
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
    equalNoise    <- noise.args$equalNoise
    discard.noise <- noise.args$discard.noise
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
    miss.list     <- control$miss.list
    miss.hc       <- control$miss.hc
    posdens       <- control$posidens
    ctrl          <- list(equalPro=control$equalPro, 
                          asMclust=control$asMclust,
                          noise.gate=(is.null(noise.gate)        || isTRUE(noise.gate)), 
                          equalNoise=(!is.null(equalNoise)       && isTRUE(equalNoise)),
                          discard.noise=(!is.null(discard.noise) && isTRUE(discard.noise)))
    control       <- control[names(control) %in% c("eps", "tol", "itmax", "equalPro")]
    control$itmax <- control$itmax[-3L]
    control$tol   <- control$tol[-3L]
    if(!miss.list) {
     if(!inherits(z.list, "list")  ||
        !all(sapply(z.list, inherits, "matrix"))) stop("'z.list' must be a list of matrices if supplied", call.=FALSE)
     if(miss.init &&
        init.z    != "list")        { 
       init.z     <- "list"
       if(isTRUE(verbose))                        message("'init.z' set to 'list' as 'z.list' was supplied\n")
      }
    }
    netmiss       <- is.null(network.data)
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
    anyg0         <- any(G == 0)
    allg0         <- all(G == 0)
    anyg1         <- any(G == 1)
    anyg0or1      <- anyg0  + anyg1
    comp.x        <- if(isTRUE(comp.x)) Nseq else comp.x
    noise.null    <- all(nnull <- is.null(noise), tnull   <- is.null(tau0), !anyg0)
    if(!is.null(noise.vol) && noise.null)         stop("Initial noise volume supplied without initial guess of noise allocations or initial guess of noise component proportion 'tau0'", call.=FALSE)
    equalNoise    <- !noise.null      && ctrl$equalNoise
    gate.noise    <- (!noise.null     && ctrl$noise.gate) || noise.null
    if(missing(G) && !noise.null) G   <- 0L:9L
    if(any(G      != floor(G))        ||
       any(G       < as.integer(noise.null)))     stop(paste0("'G' must be ", ifelse(noise.null, "strictly positive", "strictly non-negative when modelling with a noise-component")),   call.=FALSE)
    if(any(G      >= n))        {
      G           <- G[G <= n]
      if(length(G) > 1)         {                 warning("Removing G values >= the number of observations\n",  call.=FALSE, immediate.=TRUE)
      } else                                      stop("G values must be less than the number of observations", call.=FALSE)
    }

    mod.fam       <- mclust.options("emModelNames")
    range.G       <- sort(as.integer(unique(G)))
    if(anyg0 || !noise.null)    {
      if(!is.null(noise.vol))   {
        if(inherits(noise.vol, "NoiseVol"))   {
          NoiseLoc             <- noise.vol$loc
          noise.vol            <- noise.vol$vol
        }
        Vinv      <- ifelse(isTRUE(attr(noise.vol, "Inverse")), noise.vol, 1/noise.vol)
      } else if(n  > d)         {
        NoiseVol  <- noise_vol(X, noise.meth, reciprocal=TRUE)
        Vinv      <- NoiseVol$vol
        NoiseLoc  <- NoiseVol$loc
      } else                                      stop("'noise.args$noise.vol' must be specified directly for high-dimensional data when a noise component is included", call.=FALSE)
      Linv        <- log(Vinv)
      attr(Linv, "LogV")       <- TRUE
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
    sq_maha       <- !uni
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
      mfg         <- mNs
      if(!Gall)   {
        if(any(sZ <- !is.element(mNs,     mf1))){
          mf1     <- tryCatch(unname(vapply(mNs,  function(x)  switch(EXPR=x, E=, V="E", EII=, VII="EII", EEI=, VEI=, EVI=, VVI="EEI", EEE=, EVE=, VEE=, VVE=, EEV=, VEV=, EVV=, VVV="EEE"), character(1L))),
                              error=function(e) { e$message <- paste0("Invalid 'modelNames' for single component models", ifelse(uni, " for univariate data", ifelse(low.dim, "", " for high-dimensional data")), "!")
                                                  stop(e, call.=FALSE) } )
          if(isTRUE(verbose))                     message(paste0("'modelNames'", ifelse(any(sX), " further", ""), " coerced from ", paste(shQuote(mNs[sZ]), collapse=" + "), " to ", paste(shQuote(mf1[sZ]), collapse=" + "), " where G=1\n"))
        }
      } else mf1  <- mfg
    }
    mf1           <- unique(mf1)
    mfg           <- unique(mfg)
    all.mod       <- if(all(multi, !uni, Gany)) mclust.options("emModelNames") else unique(c(if(anyg0) mf0, if(any(G == 1)) mf1, if(any(G > 1)) mfg))
    multi         <- length(all.mod)    > 1L
    if(!miss.list) {
      if(length(z.list)   != len.G)               stop(paste0("'z.list' must be a list of length ", len.G),              call.=FALSE)  
      if(!all(pmax(G, 1L) ==
              vapply(z.list, ncol, numeric(1L)))) stop("Each element of 'z.list' must have 'G' columns",                 call.=FALSE)
      if(!all(n           == 
              vapply(z.list, nrow, numeric(1L)))) stop(paste0("Each element of 'z.list' must have N=", n, " rows"),      call.=FALSE) 
      exp.init$clustMD         <- FALSE
    }
    if(all(miss.list, init.z   == "list"))        stop(paste0("'z.list' must be supplied if 'init.z' is set to 'list'"), call.=FALSE)
    BICs          <- ICLs      <- AICs <- 
    DF.x          <- IT.x      <- PD.x <- provideDimnames(matrix(NA, nrow=len.G, ncol=length(all.mod)), base=list(as.character(range.G), all.mod))
    LL.x          <- replicate(nstarts  + 1L, list(BICs))
    LL.x[[1L]][]  <- -Inf
    crit.tx       <- crit.gx   <- -sqrt(.Machine$double.xmax)
    if(!netmiss)   {
      netdat      <- network.data
      if((!is.matrix(netdat)   &&
       !is.data.frame(netdat)) ||
         (ncol(netdat)   > 0   &&
       is.null(colnames(netdat))))                stop("'network.data' must be a data.frame or a matrix with named columns if supplied", call.=FALSE)
      netdat      <- call$network.data <- as.data.frame(network.data)
    }

  # Define the gating formula
    if(allg0 && gate.x)         { if(verbose)     message("Can't include gating network covariates in a noise-only model\n")
      gate.x      <- FALSE
    }
    gate.G        <- (range.G   + !noise.null) > 1 & gate.x
    if(gate.x)     {
      gvars       <- setdiff(all.vars(gating), c(".", x.names))
      if(!netmiss && length(gvars)   > 0 &&
         !all(gvars %in% names(netdat)))          stop("One or more variables in 'gating' formula not found in 'network.data'", call.=FALSE)
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
      if(gating[[3L]]   == "1 - 1")               stop("'gating' formula must include an intercept when it doesn't include covariates", call.=FALSE)
      Gn          <- G + !noise.null - !gate.noise
      if(gate.x   &&
         any(Gn   <= 1))        {
        if(all(Gn <= 1) && verbose)               message(paste0("Can't include gating network covariates ", ifelse(gate.noise, "in a single component mixture", "where G is less than 2 when 'noise.args$noise.gate' is FALSE\n")))
        gate.G[Gn <= 1]        <- FALSE
      }
      gate.names  <- stats::terms(gating)
      gate.names  <- labels(gate.names)[attr(gate.names, "order") <= 1]
    } else gating <- stats::as.formula(zN ~ 1)
    gate.noise    <- !gate.G    | gate.noise
    if(equalPro   && gate.x)    { if(verbose)     message("Can't constrain mixing proportions to be equal when gating covariates are supplied\n")
      equalPro    <- FALSE
    }
    equal.tau     <- (((range.G + !noise.null) == 1) | isTRUE(equalPro))   & !gate.G
    equal.noise   <- (((range.G + !noise.null) == 1) | isTRUE(equalNoise)) & equal.tau

  # Define the expert formula
    if(allg0 && exp.x)          { if(verbose)     message("Can't include expert network covariates in a noise-only model\n")
      exp.x       <- FALSE
    }
    if(exp.x)      {
      evars       <- setdiff(all.vars(expert), c(".", x.names))
      if(!netmiss && length(evars)   > 0 &&
         !all(evars %in% names(netdat)))          stop("One or more variables in 'expert' formula not found in 'network.data'", call.=FALSE)
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
      if(expert[[3L]]   == "1 - 1")               stop("'expert' formula must include an intercept when it doesn't include covariates", call.=FALSE)
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
    if(any(gate.x,
           exp.x) && isTRUE(ctrl$asMclust))       message("'asMclust=TRUE' invoked despite the inclusion of covariates\n")
    if(gate.x)     {
      gate.covs   <- eval(bquote(stats::model.frame(.(stats::update.formula(gating, NULL ~ .)), data=.(call$network.data), drop.unused.levels=TRUE)), envir=parent.frame(), enclos=environment())
      gate.names  <- colnames(gate.covs)
      if(any(gate.names %in% colnames(X)))        warning("Gating covariates found in response data!\n", call.=FALSE, immediate.=TRUE)
      gate.covs   <- cbind(gate.covs, eval(bquote(stats::model.frame(.(stats::as.formula(paste("~", paste(eval(bquote(all.vars(.(gating))), envir=parent.frame())[-1L], collapse="+")))), data=.(call$network.data), drop.unused.levels=TRUE)), envir=parent.frame(), enclos=environment()))
      gate.covs   <- gate.covs[,unique(colnames(gate.covs)), drop=FALSE]
      gate.char   <- vapply(gate.covs, is.character, logical(1L))
      gate.covs[gate.char]     <- lapply(gate.covs[gate.char], factor)
      netdat      <- gate.covs
    } 
    if(exp.x)      {
      expx.covs   <- eval(bquote(stats::model.frame(.(stats::update.formula(expert, NULL ~ .)), data=.(call$network.data), drop.unused.levels=TRUE)), envir=parent.frame(), enclos=environment())
      expx.names  <- colnames(expx.covs)
      if(any(expx.names %in% colnames(X)))        warning("Expert covariates found in response data!\n", call.=FALSE, immediate.=TRUE)
      expx.covs   <- cbind(expx.covs, eval(bquote(stats::model.frame(.(stats::as.formula(paste("~", paste(eval(bquote(all.vars(.(expert))), envir=parent.frame())[-1L], collapse="+")))), data=.(call$network.data), drop.unused.levels=TRUE)), envir=parent.frame(), enclos=environment()))
      expx.covs   <- expx.covs[,unique(colnames(expx.covs)), drop=FALSE]
      expx.char   <- vapply(expx.covs, is.character, logical(1L))
      expx.covs[expx.char]     <- lapply(expx.covs[expx.char], factor)
      netdat      <- if(netmiss || !gate.x) expx.covs else cbind(netdat, expx.covs)
      netdat      <- netdat[,unique(colnames(netdat)), drop=FALSE]
    }
    if((gate.x    && any(gate.char))   ||
       (exp.x     && any(expx.char)))             message("Character covariates coerced to factors\n")
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
      if(any(grepl("\\$", netnames)))             stop("Don't supply covariates to gating or expert networks using the '$' operator: use the 'network.data' argument instead", call.=FALSE)
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
    if(exp.crit   <- exp.x    && exp.init$mahalanobis) {
      crit.exp    <- 
      iter.exp    <- stats::setNames(rep(NA, len.G), paste0("G=", range.G))
    }
    if(clust.MD   <- exp.init$clustMD) {
      if((do.md   <- exp.init$joint))  {
        exp.fac   <- ncol(expx.covs)   - nct
        if((do.md <- exp.fac   > 0))   {
          if(!(has.md         <- suppressMessages(requireNamespace("clustMD", quietly=TRUE)) &&
             .version_above("clustMD", "1.2.1"))) warning("'exp.init$clustMD' not invoked - 'clustMD' library not loaded\n", call.=FALSE, immediate.=TRUE)
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
      if(!has.pkg)                                stop("'snow' package must be installed to use 'exp.init$clustMD=TRUE'", call.=FALSE)
      mdx         <- utils::capture.output( {
        mds       <- clustMD::clustMDparallel(X=XY, G=g.range, CnsIndx=CnsIndx, OrdIndx=OrdIndx, Nnorms=25000, MaxIter=500, store.params=FALSE,
                                              models=c("EII", "VII", "EEI", "VEI", "EVI", "VVI", "BD"), autoStop=TRUE, stop.tol=1e-04, scale=FALSE,
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
      Zhc         <- tryCatch(hc(data=XI, modelName=hcName, use=hcUse, minclus=min(g.range)), error=function(e) {
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
          hcZ     <- hclass(hcPairs=Zhc, G=g.range)
        }
      }
    } else hcfail <- mcfail    <- FALSE
    if(isTRUE(verbose))           message("\n################\n")

  # Loop over range of G values and initialise allocations
    G.last        <- range.G[len.G]
    MLRcon        <- TRUE
    for(g in range.G) {
      if(isTRUE(verbose))   {     message(paste0("\n", g, " cluster model", ifelse(multi, "s", ""), " -\n"))
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
      noise.gate  <- !noiseG     || gate.noise[h]
      Xinv        <- if(noiseG) Linv

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
          z.tmp   <- unmap(classification=mds$Output[[(bmd - 1L) * (len.G - anyg0or1)  + h2]]$cl, groups=Gseq)
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
            zg        <- replicate(nstarts, list(unmap(classification=sample(x=Gseq, size=noisen, replace=TRUE), groups=Gseq)))
          } else       {
            switch(EXPR=init.z, list={
            z.tmp     <- .renorm_z(z.list[[h]])
            },         {
            z.tmp     <- unmap(switch(EXPR     = init.z,
                                      hc       = hcZ[,h - anyg0 - hc1],
                                      kmeans   = stats::kmeans(x=XI, centers=g, iter.max=kiters, nstart=kstarts)$cluster,
                                      mclust   = suppressWarnings(Mclust(data=XI, G=g, x=mcl))$classification,
                                      quantile = quant_clust(x=XI, G=g),
                                      random   = sample(x=Gseq, size=noisen, replace=TRUE)), groups=Gseq)
            })
          }
        } else     {
          z.tmp   <- unmap(rep(1L, ifelse(noisen == 0 || g == 0, n, noisen)), groups=1L)
        }
      }
    for(i in if(g  > 1) startseq else 1L)      {
      if(isTRUE(multstart)     && g > 1)       {
        if(isTRUE(verbose))                      message(paste0("\n\tRandom Start: #", i, "...\n"))
        z.tmp     <- zg[[i]]
      }
      zc          <- ncol(z.tmp)
      if(zc != g  && g > 0)     {
        z.tmp     <- cbind(z.tmp, matrix(0L, ncol=g - zc, nrow=n))
      }
      z1start     <- z.tmp

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
        oldcrit   <- Inf
        newcrit   <- .Machine$double.xmax
        while(!identical(tmp.z, z.tmp) &&
              newcrit   <= oldcrit     && ix <= max.init) {
          old.z   <- tmp.z
          tmp.z   <- z.tmp
          oldcrit <- newcrit
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
                                 (!fep && !all(ep[sub] == ep[sub][1L], na.rm=TRUE)) || (fep && (nlevels(droplevels(ep[sub])) == nlevels(ep))) }, logical(1L)), drop=FALSE]
                  px           <- try(stats::lm(drop_constants(nexp, expN, pna), data=nexp, subset=pna)$fitted.values, silent=TRUE)
                  if(!inherits(px,  "try-error")) {
                    pred[pna,] <- px
                  } else {
                    init.exp   <- FALSE
                  }
                }
              }
              res              <- 
              res.G[[k]]       <- xN   - pred
              mahala[[k]]      <- MoE_mahala(exp, res, squared=sq_maha, identity=Identity)
            }
          }
          if(!init.exp)  {
            break
          } else   {
            maha  <- do.call(cbind, mahala)
            if(anyNA(maha))     {
              init.exp         <- FALSE
              break
            } else     {
              mahamin <- rowMins(maha)
              newcrit <- pmin(sum(mahamin), oldcrit)
              z.tmp   <- maha  == mahamin
              if(identical(z.tmp, old.z))         break
            }
          }
          if(g    == 1)  break
        }
        if(exp.crit) newcrit   -> crit.exp[which(g == range.G)]
        if(exp.crit) ix        -> iter.exp[which(g == range.G)]
        if(ix     >= max.init)                    warning(paste0("\tMahalanobis initialisation step failed to converge in max.init=", max.init, " iterations for the ", g, " cluster models\n"), call.=FALSE, immediate.=TRUE)
        if(noiseG && init.exp  && tnull)  {
         nRG      <- replicate(g, matrix(NA, nrow=n, ncol=d), simplify=FALSE)
         nX       <- X[noise,, drop=FALSE]
         noisexp  <- expx.covs[noise,, drop=FALSE]
         for(k in Gseq)  {
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
               } else    {
                 sub           <- z.tmp[,k] == 1
                 nexp          <- expnoise[,vapply(seq_len(ne), function(p, ep=expnoise[,p], fep=is.factor(ep)) {
                                 (!fep && !all(ep[sub] == ep[sub][1L], na.rm=TRUE)) || (fep && (nlevels(droplevels(ep[sub])) == nlevels(ep))) }, logical(1L)), drop=FALSE]
                 px            <- try(stats::predict(stats::lm(drop_constants(nexp, expN, which(noise)[pna]), data=nexp, subset=sub), newdata=noisexp[pna,colnames(nexp), drop=FALSE]), silent=TRUE)
                 if(!inherits(px,   "try-error")) {
                   pred[pna,]  <- px
                 } else  {
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
      z2start     <- z.tmp
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
        gate.pen  <- length(stats::coef(g.init)) + as.integer(!noise.null) + as.integer(!noise.gate)
       #g.init    <- glmnet::cv.glmnet(y=z, x=model.matrix(gating, data=gate.covs)[,-1L, drop=FALSE], family="multinomial", type.multinomial="ungrouped")
       #tau       <- stats::predict(g.init, type="response", newx=model.matrix(gating, data=gate.covs)[,-1L], s="lambda.1se")[,,1]
       #gate.pen  <- g.init$glmnet.fit$df[which(g.init$glmnet.fit$lambda == g.init$lambda.1se)] + as.integer(!noise.null) + 1L
        ltau      <- log(tau)
        MLRcon    <- MLRcon    && g.init$convergence == 0
      } else      {
        if(equal.pro && noiseG && !n0pro)  {
          t0      <- mean(z.init[,gN])
          tau     <- c(rep((1 - t0)/g, g), t0)
        } else    {
          tau     <- if(equal.pro) rep(1/gN, gN) else col.z/n
        }
        ltau      <- .mat_byrow(log(tau), nrow=n, ncol=gN)
        gate.pen  <- ifelse(equal.pro && g > 1, noiseG - n0pro, gN - 1L) + as.integer(noiseG)
      }
      ltau.init   <- ltau
      failedM     <- NULL
      expinitG    <- init.exp

    # Loop over the mclust model type(s)
      modtypes    <- if(g > 1)    mfg     else if(g == 1) mf1 else mf0
      T.last      <- modtypes[length(modtypes)]
      for(modtype in modtypes)  {
        m0W       <- m0X       <- ERR  <- FALSE
        denswarn  <- TRUE

      # Initialise parameters from allocations
        if(isTRUE(verbose))     { message(paste0("\n\tModel: ", modtype, "\n"))
          last.T  <- modtype   == T.last
        }
        x.df      <- ifelse(g   > 0, nVarParams(modelName=modtype, d=d, G=g), 0L) + gate.pen
        if(g > 0  && expinitG)  {
         Mstep    <- try(mstep(data=G.res, modelName=modtype, z=z.alloc, control=control), silent=TRUE)
         init.exp <- !inherits(Mstep, "try-error") && attr(Mstep, "returnCode")  >= 0
        }
        if(expold != init.exp  && !eNO) {
          failedM <- c(failedM, modtype)
        }
        if(g > 0  && !init.exp) {
          Mstep   <- try(mstep(data=X, modelName=modtype, z=if(noise.null) z.init else z.init[,-gN, drop=FALSE], control=control), silent=TRUE)
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
        ERR       <- 
        ERR2      <- ERR || ((g > 0    && attr(Mstep, "returnCode") < 0) || inherits(medens, "try-error"))
        if(denswarn      &&
           !ERR          &&
           !inherits(medens, "try-error")    &&
           any(medens     > 0)) {
          if(isTRUE(posdens))   {
            if(isTRUE(verbose))                   message("\t\tPositive log-densities occured: consider setting 'posidens' to FALSE\n")
          } else                {
            if(isTRUE(verbose))                   message("\t\tPositive log-densities occured: consider setting 'posidens' to TRUE\n")
            ERR   <- TRUE
          }
          denswarn             <- FALSE
        }
        if(isTRUE(ERR)) {
          ll      <- NA
          j       <- 1L
          if(isTRUE(verbose))     message(paste0("\t\t# Iterations: ", ifelse(ERR, "stopped at ", ""), j, ifelse(last.G && last.T, "\n\n", "\n")))
          BICs[h,modtype]      <-
          ICLs[h,modtype]      <-
          AICs[h,modtype]      <-
          DF.x[h,modtype]      <- 
          LL.x[[i  + 1L]][h,modtype]   <- -Inf
          PD.x[h,modtype]      <- !denswarn  
          IT.x[h,modtype]      <- Inf
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
           #fitE  <- glmnet::cv.glmnet(y=X, x=model.matrix(expert, data=expx.covs)[,-1L, drop=FALSE], weights=z[,k], family="mgaussian")
            e.fit[[k]]        <- fitE
            e.res[[k]]        <- stats::residuals(fitE)
           #e.res[[k]]        <- X - stats::predict(fitE, weights=z[,k], type="response", newx=model.matrix(expert, data=expx.covs)[,-1L], s="lambda.1se")
            z.mat[(k - 1L) * n + Nseq,k]      <- z[,k]
           }
           res.x  <- if(uni) as.matrix(do.call(base::c, e.res))  else do.call(rbind, e.res)
          }

        # M-step
          Mstep   <- try(if(exp.g) mstep(data=res.x, modelName=modtype, z=z.mat, control=control) else mstep(data=X, modelName=modtype, z=if(noise.null) z else z[,-gN, drop=FALSE], control=control), silent=TRUE)
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
           #fitG  <- glmnet::cv.glmnet(y=z, x=model.matrix(gating, data=gate.covs)[,-1L, drop=FALSE], family="multinomial", type.multinomial="ungrouped")
           #tau   <- stats::predict(fitG, type="response", newx=model.matrix(gating, data=gate.covs)[,-1L], s="lambda.1se")[,,1]
            ltau  <- log(tau)
            MLRcon             <- MLRcon && fitG$convergence == 0
          } else  {
            if(equal.pro && !noise.null  && !n0pro) {
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
          ERR     <- 
          ERR2    <- ERR || (attr(Mstep, "returnCode") < 0  || inherits(medens, "try-error"))
          if(denswarn    &&
             !ERR        &&
             !inherits(medens, "try-error")    &&
             any(medens   > 0)) {
            if(isTRUE(posdens)) {
              if(isTRUE(verbose))                 message("\t\tPositive log-densities occured: consider setting 'posidens' to FALSE\n")
            } else              {
              if(isTRUE(verbose))                 message("\t\tPositive log-densities occured: consider setting 'posidens' to TRUE\n")
              ERR   <- TRUE
            }
            denswarn           <- FALSE
          }
          if(isTRUE(ERR)) {
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
             dX   <- ifelse(is.numeric(ait$a) && ait$a < 0, 0L, abs(ait$linf - ll[j - 1L]))
             dX[is.nan(dX)]    <- Inf
            } else     {
             dX   <- abs(ll[j]  - ll[j - 1L])/(1 + abs(ll[j]))
            }
            stX   <- dX >= tol && j  < max.it && gN    > 1
            if(itwarn && !m0X)  {
             m0W  <- ifelse(!m0X, warnit < j - 2L, m0X)
             if(m0W   && !m0X)  {                 tryCatch(warning("WARNIT", call.=FALSE), warning=function(w)
                                                  message(paste0("\t", algo, " algorithm for the ", modtype, " model has yet to converge in 'warn.it'=", warnit, " iterations\n")))
              m0X <- TRUE
             }
            }
            if(alG == "cemEM"  && !stX  && nW == 1L)   {
              ll  <- c(ll[j - 1L], ll[j])
              alW <- "EM"
              j   <- nW        <- 2L
              stX <- TRUE
            }
          }
        } # while (j)

      # Store values corresponding to the maximum BIC/ICL/AIC so far
        j2        <- max(1L, j  - switch(EXPR=algo, cemEM=1L, 2L))
        if(isTRUE(verbose))       message(paste0("\t\t# Iterations: ", ifelse(ERR, "stopped at ", ""), j2, ifelse(last.G && last.T, "\n\n", "\n")))
       #pen.exp   <- ifelse(exp.g, g * d * (fitE$glmnet.fit$df[which(fitE$glmnet.fit$lambda == fitE$lambda.1se)] + 1), exp.pen)
        pen.exp   <- ifelse(exp.g, g * length(stats::coef(fitE)), exp.pen)
        x.df      <- pen.exp + x.df
        max.ll    <- ll[j]
        choose    <- MoE_crit(modelName=modtype, loglik=max.ll, n=n, d=d, G=g, z=z, df=x.df)
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
          PD.x[h,modtype]      <- ifelse(ERR2,  NA, !denswarn)
          x1start <- z1start
          x2start <- z2start
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
      }   else                                    warning(paste0(mdwarn, ":\ninitialisation defaulted to init.z=\"", init.z, "\" instead for the G=", paste(range.G[warnmd], collapse="/"), " model", ifelse(sum(warnmd) > 1, "s\n", "\n")), call.=FALSE, immediate.=TRUE)
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
    G             <- range.G[best.G]
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
    noise.gate    <- noise.null || gate.noise[best.G]
    if(!(bG       <- gate.G[best.G]))  {
      if(GN > 1)           {
       x.fitG     <- multinom(gating, trace=FALSE, data=gate.covs, maxit=g.itmax, reltol=g.reltol, MaxNWts=MaxNWts)
       if(equal.pro    && !noise.null &&   !n0pro)    {
         t0       <- mean(z[,GN])
         x.tau    <- c(rep((1 - t0)/G, G), t0)
       } else      {
         x.tau    <- if(n0pro) rep(1/GN, GN) else if(equal.pro)           rep(1/G, G) else x.fitG$fitted.values[1L,]
       }
       x.tau      <- stats::setNames(x.tau, paste0("Cluster", if(noise.null)     Gseq else c(Gseq, 0L)))
       MLRcon     <- MLRcon    && x.fitG$convergence == 0
      }  else      {
       x.fitG     <- suppressWarnings(stats::glm(z ~ 1, family=stats::binomial(), maxit=g.itmax))
       MLRcon     <- MLRcon    && isTRUE(x.fitG$converged)
      }
      attr(x.fitG, "Formula")  <- "~1"
      if(equal.pro)        {
       if(!noise.null          && 
          !equal.noise)    {
        x.fitG$wts[-(GN * 2L)] <- 0L  
       } else x.fitG$wts[]     <- 0L  
       x.fitG$fitted.values    <- matrix(x.tau, nrow=n, ncol=GN, byrow=TRUE)
       x.fitG$residuals        <- zN - x.fitG$fitted.values
      }
    } 
    if(isFALSE(MLRcon))                           warning(paste0("\tFor one or more models, in one or more ", algo, " iterations, the multinomial logistic regression\n\t\tin the gating network failed to converge in ", g.itmax, " iterations:\n\t\tmodify the 3rd element of the 'itmax' argument to MoE_control()\n"), call.=FALSE, immediate.=TRUE)
    if(is.matrix(x.tau))   {
      colnames(x.tau)          <- paste0("Cluster", if(!noise.null)       c(Gseq, 0L) else Gseq)  
    } else x.tau  <- stats::setNames(x.tau, paste0("Cluster", if(!noise.null) c(Gseq, 0L) else Gseq))
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
    attr(x.fitG, "Data")       <- z
    attr(x.fitG, "Maxit")      <- g.itmax
    attr(x.fitG, "Reltol")     <- g.reltol
    attr(x.fitG, "EqualPro")   <- equal.pro
    attr(x.fitG, "EqualNoise") <- n0pro
    attr(x.fitG, "Formula")    <- ifelse(is.null(attr(x.fitG, "Formula")), Reduce(paste, deparse(gating[-2L])), attr(x.fitG, "Formula"))
    attr(x.fitG, "G")          <- G
    attr(x.fitG, "NoGate")     <- if(attr(x.fitG, "Formula") == "~1") x.tau
    attr(x.fitG, "NoiseGate")  <- noise.gate
    attr(x.fitG, "NoisePro")   <- if(!noise.null && !noise.gate) ifelse(is.matrix(x.tau), x.tau[1L,GN], x.tau[GN])
    attr(x.fitG, "Noise")      <-
    attr(x.fitE, "Noise")      <- G == 0 || !noise.null
    attr(x.fitE, "Data")       <- as.data.frame(X)
    attr(x.fitE, "Formula")    <- Reduce(paste, deparse(expert[-2L]))
    attr(x.fitE, "Criterion")  <- if(exp.crit) crit.exp
    attr(x.fitE, "d")          <- d
    attr(x.fitE, "Iterations") <- if(exp.crit) iter.exp
    attr(x.fitE, "NoExp")      <- if(attr(x.fitE, "Formula") == "~1") x.mu
    class(x.fitG) <- c("MoE_gating", class(x.fitG))
    class(x.fitE) <- c("MoE_expert", class(x.fitE))
    if(G > 0) {
      vari.fin    <- x.sig
      if(exp.x)    {
        if(noise.null || noise.meth == "manual" || isTRUE(ctrl$discard.noise))    {
          z.norm  <- if(noise.null) z else .renorm_z(z[,-GN, drop=FALSE])
          if(algo == "CEM"     &&
             any(nan           <- apply(z.norm, 1L, function(x) all(is.nan(x))))) {
            z.norm[nan,]       <- .renorm_z(MoE_estep(data=lapply(x.resE, "[", nan, TRUE), mus=mus, sigs=vari.fin, 
                                            log.tau=log(if(is.matrix(x.tau)) x.tau[nan,, drop=FALSE] else x.tau), Vinv=Xinv)$z[,-GN, drop=FALSE])
          }
          if(any(nan           <- apply(z.norm, 1L, function(x) all(is.nan(x))))) {
            z.norm[nan,]       <- 1/G
          }
          fitdat  <- Reduce("+",  lapply(Gseq, function(g) z.norm[,g] * x.fitE[[g]]$fitted.values))
        } else     {
          fitdat  <- Reduce("+",  lapply(Gseq, function(g) z[,g]      * x.fitE[[g]]$fitted.values))
          fitdat  <- fitdat  +    z[,GN] * matrix(NoiseLoc, nrow=n, ncol=d, byrow=TRUE)
          z.norm  <- if(noise.null) z else z[,-GN, drop=FALSE]
        }
        mean.fin  <- crossprod(fitdat, z.norm) 
        mean.fin  <- if(GN  == 1) mean.fin/n else sweep(mean.fin, 2L, colSums2(z.norm), FUN="/", check.margin=FALSE)
      } else       {
        mean.fin  <- x.mu
      }
      if(is.matrix(mean.fin))   {
        colnames(mean.fin)     <- gnames
      } else         mean.fin  <- stats::setNames(mean.fin, gnames)      
    } else    {
      mean.fin    <- vari.fin  <- NULL
    }

    if(any(l.warn <- x.ll      != cummax(x.ll)))           {
      if(which.max(l.warn)     != length(x.ll))   warning("Log-likelihoods are not strictly increasing\n", call.=FALSE)
    }
    denswarn      <- PD.x[best.ind]
    if(isTRUE(posdens)         &&
       isTRUE(denswarn))                          warning("Optimal model contains positive log-densities; consider setting 'posidens' to FALSE\n", call.=FALSE, immediate.=TRUE)
    if(any(IT.x[!is.na(IT.x)]  == max.it))        warning(paste0("One or more models failed to converge in the maximum number of allowed iterations (", max.it, ")\n"), call.=FALSE)
    if(isTRUE(verbose))           message(paste0("\n\t\tBest Model", ifelse(length(CRITs) > 1, paste0(" (according to ", toupper(criterion), "): "), ": "), ifelse(G == 0, "single noise component", paste0(mclustModelNames(best.mod)$type, " (", best.mod, "), with ",
                                          G, " component", ifelse(G > 1, "s", ""))), ifelse(any(exp.gate) || (!noise.null && G != 0) || (equal.pro  && G != 0), paste0("\n\t\t\t  ", net.msg), ""), "\n\t\t",
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
    attr(LL.x, "posidens")     <- PD.x
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
    attr(df.fin, "nVar.Penalty")      <- nVarParams(modelName=best.mod, d=d, G=G)
    claX       <- max.col(z)
    claX[claX  == G + 1L]      <- 0L
    results       <- list(call = call2, data = as.data.frame(X), modelName = ifelse(G == 0, "", best.mod),
                          n = n, d = d, G = G, BIC = BICs, ICL = ICLs, AIC = AICs, bic = bic.fin,
                          icl = icl.fin, aic = aic.fin, gating = x.fitG, expert = x.fitE, LOGLIK = LL.x, loglik = x.ll,
                          linf = if(stopaX && G > 1) x.linf else x.ll[length(x.ll)], df = df.fin, iters = IT.x[best.ind],
                          hypvol = hypvol, parameters = list(pro = x.tau, mean = mean.fin, variance = vari.fin, Vinv = Vinv),
                          z = z, classification = stats::setNames(claX, Nseq), uncertainty = stats::setNames(uncert, Nseq),
                          net.covs = netdat, resid.data = residX, DF = DF.x, ITERS = IT.x)
    class(results)             <- "MoEClust"
    attr(results, "Algo")      <- algo
    attr(results, "Criterion") <- criterion
    attr(results, "Details")   <- paste0(best.mod, ": ", ifelse(G == 0, "single noise component", paste0(G, " component", ifelse(G > 1, "s", ""), net.msg)))
    attr(results, "Discard.N") <- ctrl$discard.noise
    attr(results, "EqualNoise")<- ctrl$equalNoise
    attr(results, "EqualPro")  <- ctrl$equalPro
    attr(results, "Expert")    <- exp.x
    attr(results, "Exp.init")  <- x2start
    attr(results, "Gating")    <- bG
    attr(results, "Init.Crit") <- if(exp.crit) crit.exp
    attr(results, "Init.Iter") <- if(exp.crit) iter.exp
    attr(results, "Noise")     <- G == 0 || !noise.null
    attr(results, "NoiseGate") <- ctrl$noise.gate
    attr(results, "Posdens")   <- denswarn
    attr(results, "Z.init")    <- x1start
      return(results)
  }

#' Density for MoEClust Mixture Models
#'
#' Computes densities (or log-densities) of observations in MoEClust mixture models.
#' @param data If there are no expert network covariates, \code{data} should be a numeric matrix or data frame, wherein rows correspond to observations (n) and columns correspond to variables (d). If there are expert network covariates, this should be a list of length G containing matrices/data.frames of (multivariate) WLS residuals for each component.
#' @param mus The mean for each of G components. If there is more than one component, this is a matrix whose k-th column is the mean of the k-th component of the mixture model. For the univariate models, this is a G-vector of means. In the presence of expert network covariates, all values should be equal to \code{0}.
#' @param sigs The \code{variance} component in the parameters list from the output to e.g. \code{\link{MoE_clust}}. The components of this list depend on the specification of \code{modelName} (see \code{\link[mclust]{mclustVariance}} for details). The number of components \code{G}, the number of variables \code{d}, and the \code{modelName} are inferred from \code{sigs}.
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
#' @author Keefe Murphy - <\email{keefe.murphy@@mu.ie}>
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
    Vinv          <- if(is.null(Vinv) || isTRUE(attr(Vinv, "LogV"))) Vinv else log(Vinv)
    Ldat          <- inherits(data, "list")
    if(!Ldat      || (Ldat &&
       length(data)        != max(G, 1L)))       {
      data        <- replicate(G, as.matrix(data), simplify=FALSE)
    } else data   <- lapply(data, as.matrix)
    dat1          <- data[[1L]]
    n             <- NROW(dat1)
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
        sigx      <- diag(sigmas, d);
        densi     <- vapply(Gseq, function(k) dmvn(data[[k]], mus[,k], sigx,                        log=TRUE, isChol=TRUE),  numeric(n))
      }, VII= {
        densi     <- vapply(Gseq, function(k) dmvn(data[[k]], mus[,k], diag(sigmas[k], d),          log=TRUE, isChol=TRUE),  numeric(n))
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
    densi         <- if(Vnul) densi else if(G     > 0) bind2(densi, Vinv) else matrix(Vinv, nrow=n, ncol=G + !Vnul)
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
#' @author Keefe Murphy - <\email{keefe.murphy@@mu.ie}>
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
      Dens        <- do.call(MoE_dens, list(data=data, mus=mus, sigs=sigs, log.tau=log.tau, Vinv=Vinv))
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
#' @author Keefe Murphy - <\email{keefe.murphy@@mu.ie}>
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
      Dens        <- do.call(MoE_dens, list(data=data, mus=mus, sigs=sigs, log.tau=log.tau, Vinv=Vinv))
    } else if(!is.matrix(Dens) ||
              !is.numeric(Dens))                  stop("'Dens' must be a numeric matrix", call.=FALSE)
    Gseq          <- seq_len(ncol(Dens))
    z             <- provideDimnames(unmap(max.col(Dens), groups=Gseq), base=list(as.character(seq_len(nrow(Dens))), paste0("Cluster", Gseq)))
      return(list(z = z, loglik = sum(z * Dens, na.rm=TRUE)))
  }

#' MoEClust BIC, ICL, and AIC Model-Selection Criteria
#'
#' Computes the BIC (Bayesian Information Criterion), ICL (Integrated Complete Likelihood), and AIC (Akaike Information Criterion) for parsimonious mixture of experts models given the log-likelihood, the dimension of the data, the number of mixture components in the model, the numbers of parameters in the gating and expert networks respectively, and, for the ICL, the numbers of observations in each component.
#' @param modelName A character string indicating the model. The help file for \code{\link[mclust]{mclustModelNames}} describes the available models. Not necessary if \code{df} is supplied.
#' @param loglik The log-likelihood for a data set with respect to the Gaussian mixture model specified in the \code{modelName} argument.
#' @param n,d,G The number of observations in the data, dimension of the data, and number of components in the Gaussian mixture model, respectively, used to compute \code{loglik}. \code{d} & \code{G} are not necessary if \code{df} is supplied.
#' @param gating.pen The number of parameters of the \emph{gating} network of the MoEClust model. Defaults to \code{G - 1}, which corresponds to no gating covariates. If covariates are included, this should be the number of regression coefficients in the fitted \code{gating} object. If there are no covariates and mixing proportions are further assumed to be present in equal proportion, \code{gating.pen} should be \code{0}. The number of parameters used in the estimation of the noise component, if any, should also be included. Not necessary if \code{df} is supplied.
#' @param expert.pen The number of parameters of the \emph{expert} network of the MoEClust model. Defaults to \code{G * d}, which corresponds to no expert covariates. If covariates are included, this should be the number of regression coefficients in the fitted \code{expert} object. Not necessary if \code{df} is supplied.
#' @param z The \code{n} times \code{G} responsibility matrix whose \code{[i,k]}-th entry is the probability that observation \emph{i} belongs to the \emph{k}-th component.. If supplied the ICL is also computed and returned, otherwise only the BIC and AIC.
#' @param df An alternative way to specify the number of estimated parameters (or 'used' degrees of freedom) exactly. If supplied, the arguments \code{modelName}, \code{d}, \code{G}, \code{gating.pen}, and \code{expert.pen}, which are used to calculate the number of parameters, will be ignored. The number of parameters used in the estimation of the noise component, if any, should also be included.
#'
#' @details The function is vectorised with respect to the arguments \code{modelName} and \code{loglik}.
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
#' @author Keefe Murphy - <\email{keefe.murphy@@mu.ie}>
#'
#' @references Biernacki, C., Celeux, G. and Govaert, G. (2000). Assessing a mixture model for clustering with the integrated completed likelihood. \emph{IEEE Trans. Pattern Analysis and Machine Intelligence}, 22(7): 719-725.
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
#' identical(unname(bic2), model$bic) #TRUE
#'
#' # Make the same comparison with the known number of estimated parameters
#' (bic3 <- MoE_crit(loglik=ll, n=n, df=model$df, z=z)["bic",])
#' identical(bic3, bic2)              #TRUE
  MoE_crit        <- Vectorize(function(modelName, loglik, n, d, G, gating.pen = G - 1L, expert.pen = G * d, z = NULL, df = NULL) {
    df            <- ifelse(!missing(df), df, nVarParams(modelName=modelName, d=d, G=G) + expert.pen + gating.pen)
    double.ll     <- 2 * loglik
    bic.x         <- double.ll  - df * log(n)
    aic.x         <- double.ll  - df * 2
      return(c(bic = bic.x, icl = if(!missing(z)) bic.x + 2L * sum(log(rowMaxs(z)), na.rm=TRUE), aic = aic.x, df = df))
  }, vectorize.args = c("modelName", "loglik"), SIMPLIFY="array")

#' Set control values for use with MoEClust
#'
#' Supplies a list of arguments (with defaults) for use with \code{\link{MoE_clust}}.
#' @param init.z The method used to initialise the cluster labels. Defaults to \code{"hc"}, i.e. model-based agglomerative hierarchical clustering tree as per \code{\link[mclust]{hc}}, for multivariate data (see \code{hc.args}), or \code{"quantile"}-based clustering as per \code{\link{quant_clust}} for univariate data (unless there are expert network covariates incorporated via \code{exp.init$joint} &/or \code{exp.init$clustMD}, in which case the default is again \code{"hc"}). The \code{"quantile"} option is thus only available for univariate data when expert network covariates are not incorporated via \code{exp.init$joint} &/or \code{exp.init$clustMD}, or when expert network covariates are not supplied.
#'
#' Other options include \code{"kmeans"} (see \code{km.args}), \code{"random"} initialisation, a user-supplied \code{"list"}, and a full run of \code{\link[mclust]{Mclust}} (itself initialised via a model-based agglomerative hierarchical clustering tree, again see \code{hc.args}), although this last option \code{"mclust"} will be coerced to \code{"hc"} if there are no \code{gating} &/or \code{expert} covariates within \code{\link{MoE_clust}} (in order to better reproduce \code{\link[mclust]{Mclust}} output).
#'
#' When \code{init.z="list"}, \code{exp.init$clustMD} is forced to \code{FALSE}; otherwise, when \code{isTRUE(exp.init$clustMD)} and the \code{\link[clustMD]{clustMD}} library is loaded, the \code{init.z} argument instead governs the method by which a call to \code{\link[clustMD]{clustMD}} is initialised. In this instance, \code{"quantile"} will instead default to \code{"hc"}, and the arguments to \code{hc.args} and \code{km.args} will be ignored (unless all \code{\link[clustMD]{clustMD}} model types fail for a given number of components).
#'
#' When \code{init.z="mclust"} or \code{\link[clustMD]{clustMD}} is successfully invoked (via \code{exp.init$clustMD}), the argument \code{init.crit} (see below) specifies the model-selection criterion (\code{"bic"} or \code{"icl"}) by which the optimal \code{\link[mclust]{Mclust}} or \code{\link[clustMD]{clustMD}} model type to initialise with is determined, and \code{criterion} remains unaffected.
#' @param noise.args A list supplying select named parameters to control inclusion of a noise component in the estimation of the mixture. If either or both of the arguments \code{tau0} &/or \code{noise.init} are supplied, a noise component is added to the the model in the estimation.
#' \describe{
#' \item{\code{tau0}}{Prior mixing proportion for the noise component. If supplied, a noise component will be added to the model in the estimation, with \code{tau0} giving the prior probability of belonging to the noise component for \emph{all} observations. Typically supplied as a scalar in the interval (0, 1), e.g. \code{0.1}. Can be supplied as a vector when gating covariates are present and \code{noise.args$noise.gate} is \code{TRUE}. This argument can be supplied instead of or in conjunction with the argument \code{noise.init} below.}
#' \item{\code{noise.init}}{A logical or numeric vector indicating an initial guess as to which observations are noise in the data. If numeric, the entries should correspond to row indices of the data. If supplied, a noise component will be added to the model in the estimation. This argument can be used in conjunction with \code{tau0} above, or can be replaced by that argument also.}
#' \item{\code{noise.gate}}{A logical indicating whether gating network covariates influence the mixing proportion for the noise component, if any. Defaults to \code{TRUE}, but leads to greater parsimony if \code{FALSE}. Only relevant in the presence of a noise component; only effects estimation in the presence of gating covariates.}
#' \item{\code{noise.meth}}{The method used to estimate the volume when a noise component is invoked. Defaults to \code{\link[mclust]{hypvol}}. For univariate data, this argument is ignored and the range of the data is used instead (unless \code{noise.vol} below is specified). The options \code{"convexhull"} and \code{"ellipsoidhull"} require loading the \code{geometry} and \code{cluster} libraries, respectively. This argument is only relevant if \code{noise.vol} below is not supplied.}
#' \item{\code{noise.vol}}{This argument can be used to override the argument \code{noise.meth} by specifying the (hyper)volume directly, i.e. specifying an improper uniform density. This will override the use of the range of the response data for univariate data if supplied. Note that the (hyper)volume, rather than its inverse, is supplied here. This can affect prediction and the location of the MVN ellipses for \code{\link{MoE_gpairs}} plots (see \code{\link{noise_vol}}).}
#' \item{\code{equalNoise}}{Logical which is \strong{only} invoked when \code{isTRUE(equalPro)} and gating covariates are not supplied. Under the default setting (\code{FALSE}), the mixing proportion for the noise component is estimated, and remaining mixing proportions are equal; when \code{TRUE} all components, including the noise component, have equal mixing proportions.}
#' \item{\code{discard.noise}}{A logical governing how the means are summarised in \code{parameters$mean} and by extension the location of the MVN ellipses in \code{\link{MoE_gpairs}} plots for models with \emph{both} expert network covariates and a noise component (otherwise this argument is irrelevant). 
#' 
#' The means for models with expert network covariates are summarised by the posterior mean of the fitted values. By default (\code{FALSE}), the mean of the noise component is accounted for in the posterior mean. Otherwise, or when the mean of the noise component is unavailable (due to having been manually supplied via \code{noise.args$noise.vol}), the \code{z} matrix is renormalised after discarding the column corresponding to the noise component prior to computation of the posterior mean. The renormalisation approach can be forced by specifying \code{noise.args$discard.noise=TRUE}, even when the mean of the noise component is available. For models with a noise component fitted with \code{algo="CEM"}, a small extra E-step is conducted for observations assigned to the non-noise components in this case.}
#' }
#' In particular, the argument \code{noise.meth} will be ignored for high-dimensional \code{n <= d} data, in which case the argument \code{noise.vol} \emph{must be} specified. Note that this forces \code{noise.args$discard.noise} to \code{TRUE}. See \code{\link{noise_vol}} for more details.
#' 
#' The arguments \code{tau0} and \code{noise.init} can be used separately, to provide alternative means to invoke a noise component. However, they can also be supplied together, in which case observations corresponding to \code{noise.init} have probability \code{tau0} (rather than 1) of belonging to the noise component.
#' @param asMclust The default values of \code{stopping} and \code{hc.args$hcUse} (see below) are such that results for models with \emph{no covariates in either network} are liable to differ from results for equivalent models obtained via \code{\link[mclust]{Mclust}}. \pkg{MoEClust} uses \code{stopping="aitken"} and \code{hcUse="VARS"} by default, while \pkg{mclust} always implicitly uses \code{stopping="relative"} and defaults to \code{hcUse="SVD"}.
#' 
#' \code{asMclust} is a logical variable (\code{FALSE}, by default) which functions as a simple convenience tool for overriding these two arguments (even if explicitly supplied!) such that they behave like the function \code{\link[mclust]{Mclust}}. Other \emph{user-specified} arguments which differ from \pkg{mclust} are not affected by \code{asMclust}, as their defaults already correspond to \pkg{mclust}. Results may still differ slightly as \pkg{MoEClust} calculates log-likelihood values with greater precision. Finally, note that \code{asMclust=TRUE} is invoked even for models with covariates which are not accommodated by \pkg{mclust}.
#' @param equalPro Logical variable indicating whether or not the mixing proportions are to be constrained to be equal in the model. Default: \code{equalPro = FALSE}. Only relevant when \code{gating} covariates are \emph{not} supplied within \code{\link{MoE_clust}}, otherwise ignored. In the presence of a noise component (see \code{noise.args}), only the mixing proportions for the non-noise components are constrained to be equal (by default, see \code{equalNoise}), after accounting for the noise component.
#' @param exp.init A list supplying select named parameters to control the initialisation routine in the presence of \emph{expert} network covariates (otherwise ignored):
#' \describe{
#' \item{\code{joint}}{A logical indicating whether the initial partition is obtained on the joint distribution of the response and expert network covariates (defaults to \code{TRUE}) or just the response variables (\code{FALSE}). By default, only continuous expert network covariates are considered (see \code{exp.init$clustMD} below). Only relevant when \code{init.z} is not \code{"random"} (unless \code{isTRUE(exp.init$clustMD)}, in which case \code{init.z} specifies the initialisation routine for a call to \code{\link[clustMD]{clustMD}}). This will render the \code{"quantile"} option to \code{init.z} for univariate data unusable if continuous expert network covariates are supplied &/or categorical/ordinal expert network covariates are supplied when \code{isTRUE(exp.init$clustMD)} and the \code{\link[clustMD]{clustMD}} library is loaded.}
#' \item{\code{mahalanobis}}{A logical indicating whether to iteratively reallocate observations during the initialisation phase to the component corresponding to the expert network regression to which it's closest to the fitted values of in terms of Mahalanobis distance (defaults to \code{TRUE}). This will ensure that each component can be well modelled by a single expert prior to running the EM/CEM algorithm.}
#' \item{\code{\link[clustMD]{clustMD}}}{A logical indicating whether categorical/ordinal covariates should be incorporated when using the joint distribution of the response and expert network covariates for initialisation (defaults to \code{FALSE}). Only relevant when \code{isTRUE(exp.init$joint)}. Requires the use of the \code{\link[clustMD]{clustMD}} library. Note that initialising in this manner involves fitting all \code{\link[clustMD]{clustMD}} model types in parallel for all numbers of components considered, and may fail (especially) in the presence of nominal expert network covariates.
#'
#' Unless \code{init.z="list"}, supplying this argument as \code{TRUE} when the \code{\link[clustMD]{clustMD}} library is loaded has the effect of superseding the \code{init.z} argument: this argument now governs instead how the call to \code{\link[clustMD]{clustMD}} is initialised (unless all \code{\link[clustMD]{clustMD}} model types fail for a given number of components, in which case \code{init.z} is invoked \emph{instead} to initialise for \code{G} values for which all \code{\link[clustMD]{clustMD}} model types failed). Similarly, the arguments \code{hc.args} and \code{km.args} will be ignored (again, unless all \code{\link[clustMD]{clustMD}} model types fail for a given number of components).}
#' \item{\code{max.init}}{The maximum number of iterations for the Mahalanobis distance-based reallocation procedure when \code{exp.init$mahalanobis} is \code{TRUE}. Defaults to \code{.Machine$integer.max}.}
#' \item{\code{identity}}{A logical indicating whether the identity matrix (corresponding to the use of the Euclidean distance) is used in place of the covariance matrix of the residuals (corresponding to the use of the Mahalanobis distance). Defaults to \code{FALSE}; only relevant for multivariate response data.}
#' \item{\code{drop.break}}{When \code{isTRUE(exp.init$mahalanobis)} observations will be completely in or out of a component during the initialisation phase. As such, it may occur that constant columns will be present when building a given component's expert regression (particularly for categorical covariates). It may also occur, due to this partitioning, that "unseen" data, when calculating the residuals, will have new factor levels. When \code{isTRUE(exp.init$drop.break)}, the Mahalanobis distance based initialisation phase will explicitly fail in either of these scenarios.
#'
#' Otherwise, \code{\link{drop_constants}} and \code{\link{drop_levels}} will be invoked when \code{exp.init$drop.break} is \code{FALSE} (the default) to \emph{try} to remedy the situation. In any case, only a warning that the initialisation step failed will be printed, regardless of the value of \code{exp.init$drop.break}.}
#' }
#' @param algo Switch controlling whether models are fit using the \code{"EM"} (the default) or \code{"CEM"} algorithm. The option \code{"cemEM"} allows running the EM algorithm starting from convergence of the CEM algorithm.
#' @param criterion When either \code{G} or \code{modelNames} is a vector, \code{criterion} determines whether the \code{"bic"} (Bayesian Information Criterion), \code{"icl"} (Integrated Complete Likelihood), \code{"aic"} (Akaike Information Criterion) is used to determine the 'best' model when gathering output. Note that all criteria will be returned in any case.
#' @param stopping The criterion used to assess convergence of the EM/CEM algorithm. The default (\code{"aitken"}) uses Aitken's acceleration method via \code{\link{aitken}}, otherwise the \code{"relative"} change in log-likelihood is monitored (which may be less strict). The \code{"relative"} option corresponds to the stopping criterion used by \code{\link[mclust]{Mclust}}: see \code{asMclust} above. 
#' 
#' Both stopping rules are ultimately governed by \code{tol[1]}. When the \code{"aitken"} method is employed, the asymptotic estimate of the final converged maximised log-likelihood is also returned as \code{linf} for models with 2 or more components, though the largest element of the returned vector \code{loglik} still gives the log-likelihood value achieved by the parameters returned at convergence, under both \code{stopping} methods (see \code{\link{MoE_clust}}).
#' @param z.list A user supplied list of initial cluster allocation matrices, with number of rows given by the number of observations, and numbers of columns given by the range of component numbers being considered. Only relevant if \code{init.z == "z.list"}. These matrices are allowed correspond to both soft or hard clusterings, and will be internally normalised so that the rows sum to 1.
#' @param nstarts The number of random initialisations to use when \code{init.z="random"}. Defaults to \code{1}. Results will be based on the random start yielding the highest estimated log-likelihood. Note that all \code{nstarts} random initialisations are affected by \code{exp.init$mahalanobis}, if invoked in the presence of expert network covariates, which may remove some of the randomness.
#' @param eps A scalar tolerance associated with deciding when to terminate computations due to computational singularity in covariances. Smaller values of \code{eps} allow computations to proceed nearer to singularity. The default is the relative machine precision \code{.Machine$double.eps}, which is approximately \emph{2e-16} on IEEE-compliant machines.
#' @param tol A vector of length three giving \emph{relative} convergence tolerances for 1) the log-likelihood of the EM/CEM algorithm, 2) parameter convergence in the inner loop for models with iterative M-step (\code{"VEI", "VEE", "EVE", "VVE", "VEV"}), and 3) optimisation in the multinomial logistic regression in the gating network, respectively. The default is \code{c(1e-05, sqrt(.Machine$double.eps), 1e-08)}. If only one number is supplied, it is used as the tolerance for all three cases given.
#' @param itmax A vector of length three giving integer limits on the number of iterations for 1) the EM/CEM algorithm, 2) the inner loop for models with iterative M-step (\code{"VEI", "VEE", "EVE", "VVE", "VEV"}), and 3) the multinomial logistic regression in the gating network, respectively.
#'
#' The default is \code{c(.Machine$integer.max, .Machine$integer.max, 1000L)}, allowing termination to be completely governed by \code{tol[1]} & \code{tol[2]} for the inner and outer loops of the EM/CEM algorithm. If only one number is supplied, it is used as the iteration limit for the outer loop only and the other elements of \code{itmax} retain their usual defaults.
#' 
#' If, for any model with gating covariates, the multinomial logistic regression in the gating network fails to converge in \code{itmax[3]} iterations at any stage of the EM/CEM algorithm, an appropriate warning will be printed, prompting the user to modify this argument.
#' @param hc.args A list supplying select named parameters to control the initialisation of the cluster allocations when \code{init.z="hc"} (or when \code{init.z="mclust"}, which itself relies on \code{\link[mclust]{hc}}), unless \code{isTRUE(exp.init$clustMD)}, the \code{\link[clustMD]{clustMD}} library is loaded, and none of the \code{\link[clustMD]{clustMD}} model types fail (otherwise irrelevant):
#' \describe{
#' \item{\code{hcUse}}{A string specifying the type of input variables to be used. This defaults to \code{"VARS"} here, unlike \pkg{mclust} which defaults to \code{"SVD"}. Other allowable values are documented in \code{\link[mclust]{mclust.options}}. See \code{asMclust} above.}
#' \item{\code{hc.meth}}{A character string indicating the model to be used when hierarchical clustering (see \code{\link[mclust]{hc}}) is employed for initialisation (either when \code{init.z="hc"} or \code{init.z="mclust"}). Defaults to \code{"EII"} for high-dimensional data, or \code{"VVV"} otherwise.}
#' }
#' @param km.args A list supplying select named parameters to control the initialisation of the cluster allocations when \code{init.z="kmeans"}, unless \code{isTRUE(exp.init$clustMD)}, the \code{\link[clustMD]{clustMD}} library is loaded, and none of the \code{\link[clustMD]{clustMD}} model types fail (otherwise irrelevant):
#' \describe{
#' \item{\code{kstarts}}{The number of random initialisations to use. Defaults to 10.}
#' \item{\code{kiters}}{The maximum number of K-Means iterations allowed. Defaults to 10.}
#' }
#' @param init.crit The criterion to be used to determine the optimal model type to initialise with, when \code{init.z="mclust"} or when \code{isTRUE(exp.init$clustMD)} and the \code{\link[clustMD]{clustMD}} library is loaded (one of \code{"bic"} or \code{"icl"}). Defaults to \code{"icl"} when \code{criterion="icl"}, otherwise defaults to \code{"bic"}. The \code{criterion} argument remains unaffected.
#' @param posidens A logical governing whether to continue running the algorithm even in the presence of positive log-densities. Defaults to \code{TRUE}, but setting \code{posidens=FALSE} can help to safeguard against spurious solutions, which will be instantly terminated if positive log-densities are encountered. Note that versions of this package prior to and including version 1.3.1 always implicitly assumed \code{posidens=FALSE}.
#' @param warn.it A single number giving the iteration count at which a warning will be printed if the EM/CEM algorithm has failed to converge. Defaults to \code{0}, i.e. no warning (which is true for any \code{warn.it} value less than \code{3}), otherwise the message is printed regardless of the value of \code{verbose}. If non-zero, \code{warn.it} should be moderately large, but obviously less than \code{itmax[1]}. A warning will always be printed if one of more models fail to converge in \code{itmax[1]} iterations.
#' @param MaxNWts The maximum allowable number of weights in the call to \code{\link[nnet]{multinom}} for the multinomial logistic regression in the gating network. There is no intrinsic limit in the code, but increasing \code{MaxNWts} will probably allow fits that are very slow and time-consuming. It may be necessary to increase \code{MaxNWts} when categorical concomitant variables with many levels are included or the number of components is high.
#' @param verbose Logical indicating whether to print messages pertaining to progress to the screen during fitting. By default is \code{TRUE} if the session is interactive, and \code{FALSE} otherwise. If \code{FALSE}, warnings and error messages will still be printed to the screen, but everything else will be suppressed.
#' @param ... Catches unused arguments.
#'
#' @details \code{\link{MoE_control}} is provided for assigning values and defaults within \code{\link{MoE_clust}} and \code{\link{MoE_stepwise}}.
#'
#' While the \code{criterion} argument controls the choice of the optimal number of components and GPCM/\pkg{mclust} model type, \code{\link{MoE_compare}} is provided for choosing between fits with different combinations of covariates or different initialisation settings.
#' @importFrom mclust "hc" "hypvol" "Mclust" "mclust.options"
#' @importFrom nnet "multinom"
#' @note Note that successfully invoking \code{exp.init$clustMD} (though it defaults to \code{FALSE}) affects the role of the arguments \code{init.z}, \code{hc.args}, and \code{km.args}. Please read the documentation above carefully in this instance.
#' 
#' The initial allocation matrices before and after the invocation of the \code{exp.init} related arguments are both stored as attributes in the object returned by \code{\link{MoE_clust}} (named \code{"Z.init"} and \code{"Exp.init"}, respectively). If \code{init.z="random"} and \code{nstarts > 1}, the allocations corresponding to the best random start are stored. This can be useful for supplying \code{z.list} for future fits.
#' @return A named list in which the names are the names of the arguments and the values are the values supplied to the arguments.
#' @export
#' @keywords control
#' @author Keefe Murphy - <\email{keefe.murphy@@mu.ie}>
#'
#' @seealso \code{\link{MoE_clust}}, \code{\link{MoE_stepwise}}, \code{\link{aitken}}, \code{\link[mclust]{Mclust}}, \code{\link[mclust]{hc}}, \code{\link[mclust]{mclust.options}}, \code{\link{quant_clust}}, \code{\link[clustMD]{clustMD}}, \code{\link{noise_vol}}, \code{\link[mclust]{hypvol}}, \code{\link[geometry]{convhulln}}, \code{\link[cluster]{ellipsoidhull}}, \code{\link{MoE_compare}}, \code{\link[nnet]{multinom}}
#' @usage
#' MoE_control(init.z = c("hc", "quantile", "kmeans", "mclust", "random", "list"),
#'             noise.args = list(...),
#'             asMclust = FALSE,
#'             equalPro = FALSE,
#'             exp.init = list(...),
#'             algo = c("EM", "CEM", "cemEM"),
#'             criterion = c("bic", "icl", "aic"),
#'             stopping = c("aitken", "relative"),
#'             z.list = NULL, 
#'             nstarts = 1L,
#'             eps = .Machine$double.eps,
#'             tol = c(1e-05, sqrt(.Machine$double.eps), 1e-08),
#'             itmax = c(.Machine$integer.max, .Machine$integer.max, 1000L),
#'             hc.args = list(...),
#'             km.args = list(...),
#'             posidens = TRUE,
#'             init.crit = c("bic", "icl"),
#'             warn.it = 0L,
#'             MaxNWts = 1000L,
#'             verbose = interactive(),
#'             ...)
#' @examples
#' ctrl1 <- MoE_control(criterion="icl", itmax=100, warn.it=15, init.z="random", nstarts=5)
#'
#' data(CO2data)
#' GNP   <- CO2data$GNP
#' \donttest{res   <- MoE_clust(CO2data$CO2, G=2, expert = ~ GNP, control=ctrl1)
#'
#' # Alternatively, specify control arguments directly
#' res2  <- MoE_clust(CO2data$CO2, G=2, expert = ~ GNP, stopping="relative")}
#'
#' # Supplying ctrl1 without naming it as 'control' can throw an error
#' \dontrun{
#' res3  <- MoE_clust(CO2data$CO2, G=2, expert = ~ GNP, ctrl1)}
#' 
#' # Similarly, supplying control arguments via a mix of the ... construct
#' # and the named argument 'control' also throws an error
#' \dontrun{
#' res4  <- MoE_clust(CO2data$CO2, G=2, expert = ~ GNP, control=ctrl1, init.z="kmeans")}
#'
#' \donttest{# Initialise via the mixed-type joint distribution of response & covariates
#' # Let the ICL criterion determine the optimal clustMD model type
#' # Constrain the mixing proportions to be equal
#' ctrl2 <- MoE_control(exp.init=list(clustMD=TRUE), init.crit="icl", equalPro=TRUE)
#' data(ais)
#' library(clustMD)
#' res4  <- MoE_clust(ais[,3:7], G=2, modelNames="EVE", expert= ~ sex,
#'                    network.data=ais, control=ctrl2)
#'
#' # Include a noise component by specifying its prior mixing proportion
#' res5  <- MoE_clust(ais[,3:7], G=2, modelNames="EVE", expert= ~ sex,
#'                    network.data=ais, tau0=0.1)}
  MoE_control     <- function(init.z = c("hc", "quantile", "kmeans", "mclust", "random", "list"), noise.args = list(...), asMclust = FALSE, equalPro = FALSE, exp.init = list(...), algo = c("EM", "CEM", "cemEM"), 
                              criterion = c("bic", "icl", "aic"), stopping = c("aitken", "relative"), z.list = NULL, nstarts = 1L, eps = .Machine$double.eps, tol = c(1e-05, sqrt(.Machine$double.eps), 1e-08),
                              itmax = c(.Machine$integer.max, .Machine$integer.max, 1000L), hc.args = list(...), km.args = list(...), posidens = TRUE, init.crit = c("bic", "icl"), warn.it = 0L, MaxNWts = 1000L, verbose = interactive(), ...) {
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
    miss.list     <- missing(z.list)
    if(!miss.init && (length(init.z)        > 1 ||
       !is.character(init.z)))                    stop("'init.z' must be a single character string",    call.=FALSE)
    init.z        <- match.arg(init.z)
    if(init.z     == "random")      {
      if(length(nstarts)    != 1   ||
         !is.numeric(nstarts)      ||
         (nstarts            < 1   ||
          floor(nstarts)    != nstarts))          stop(paste0("'nstarts' must be a single integer >= 1 when 'init.z'=", init.z), call.=FALSE)
    } else if(!missing(nstarts))    {
      if(length(nstarts)    == 1   &&
         is.numeric(nstarts)       &&
         nstarts             > 1)   {             warning("'nstarts' is > 1 but 'init.z' is not set to \"random\"\n",      call.=FALSE, immediate.=TRUE)
      } else                                      warning("'nstarts' is supplied but 'init.z' is not set to \"random\"\n", call.=FALSE, immediate.=TRUE)
    }
    
    exp.init      <- if(!is.null(exp.init)      && inherits(exp.init,   "list")) exp.init
    hc.args       <- if(!is.null(hc.args)       && inherits(hc.args,    "list")) hc.args
    km.args       <- if(!is.null(km.args)       && inherits(km.args,    "list")) km.args
    noise.args    <- if(!is.null(noise.args)    && inherits(noise.args, "list")) noise.args
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
    if(is.null(exp.init$identity))          {
      exp.init$identity     <- FALSE
    } else if(length(exp.init$identity)     > 1 ||
              !is.logical(exp.init$identity))     stop("'exp.init$identity' must be a single logical indicator",        call.=FALSE)
    if(is.null(exp.init$max.init))          {
      exp.init$max.init     <- .Machine$integer.max
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
    if(len.itmax  == 1)      itmax <- c(itmax, .Machine$integer.max, 1000L)
    if(len.itmax  == 2)                           stop("'itmax' must be of length 1 or 3",   call.=FALSE)
    inf           <- is.infinite(itmax)
    if(any(inf))        itmax[inf] <- .Machine$integer.max
    if(length(MaxNWts)  > 1 ||
       !is.numeric(MaxNWts) ||
       MaxNWts    <= 0)                           stop("'MaxNWts' must be a strictly positive scalar",  call.=FALSE)
    if(length(asMclust) > 1 ||
       !is.logical(asMclust))                     stop("'asMclust' must be a single logical indicator", call.=FALSE)
    if(length(equalPro) > 1 ||
       !is.logical(equalPro))                     stop("'equalPro' must be a single logical indicator", call.=FALSE)
    if(length(posidens) > 1 ||
       !is.logical(posidens))                     stop("'posidens' must be a single logical indicator", call.=FALSE)

    if(isTRUE(asMclust))     {
      if(stopping != "relative")                  message("'stopping' forced to \"relative\" as a result of 'asMclust=TRUE'\n")
      stopping    <- "relative"
      if(!is.null(hc.args$hcUse)   &&
         hc.args$hcUse      != "SVD")             message("'hcUse' forced to \"SVD\" as a result of 'asMclust=TRUE'\n")
      hc.args$hcUse                <- "SVD"
    }
    if(!is.null(noise.args$noise.init)          ||
       !is.null(noise.args$tau0))                {
      if(!is.null(noise.args$tau0)              &&
        (!is.numeric(noise.args$tau0)           ||
        any(noise.args$tau0 <= 0)  ||
        any(noise.args$tau0 >= 1)))               stop("'noise.args$tau0' must lie in the interval (0, 1)",  call.=FALSE)
      if(!is.null(noise.args$noise.vol))         {
        NoiseVol  <- inherits(noise.args$noise.vol, "NoiseVol")
        if(!NoiseVol        &&
          (length(noise.args$noise.vol) > 1     ||
           noise.args$noise.vol    <= 0))         stop("Invalid 'noise.args$noise.vol'", call.=FALSE)
        noise.args$noise.meth      <- ifelse(NoiseVol, attr(noise.args$noise.vol, "Method"), "manual")
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
      if(noise.args$equalNoise     && !equalPro &&
         isTRUE(verbose))                         message("'equalNoise' forced to FALSE as 'equalPro' is FALSE\n")
      noise.args$equalNoise        <- equalPro  && noise.args$equalNoise
      if(is.null(noise.args$discard.noise))      {
        noise.args$discard.noise   <- FALSE
      } else if(length(noise.args$discard.noise) > 1   ||
        !is.logical(noise.args$discard.noise))    stop("noise.args$discard.noise' must be a single logical indicator",  call.=FALSE)
      has.lib     <- switch(EXPR=noise.args$noise.meth, 
                            manual=, hypvol=TRUE, 
                            convexhull=          {
                              suppressMessages(requireNamespace("geometry", quietly=TRUE)) && .version_above("geometry", "0.4.0")
                            }, 
                            ellipsoidhull=       {
                              suppressMessages(requireNamespace("cluster",  quietly=TRUE)) && .version_above("cluster",  "1.4.0") })
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
    }

    if(is.null(km.args$kiters))     {
      km.args$kiters        <- 10L
    } else {
      kiters      <- km.args$kiters
      if(length(kiters)     != 1   ||
       !is.numeric(kiters)  ||
       kiters      < 1      ||
       kiters     != floor(kiters))               stop("'km.args$kiters' must be a single strictly positive integer",   call.=FALSE)
    }
    if(is.null(km.args$kstarts))    {
      km.args$kstarts       <- 10L
    } else {
      kstarts     <- km.args$kstarts
      if(length(kstarts)    != 1   ||
       !is.numeric(kstarts) ||
       kstarts     < 1      ||
       kstarts    != floor(kstarts))              stop("'km.argss$kstarts' must be a single strictly positive integer", call.=FALSE)
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
      list(algo = algo, criterion = criterion, stopping = stopping, init.z = init.z, nstarts = nstarts, exp.init = exp.init, eps = eps, tol = tol, itmax = itmax, MaxNWts = MaxNWts, equalPro = equalPro, noise.args = noise.args, 
           hc.args = hc.args, km.args = km.args, init.crit = init.crit, warn.it = warn.it, verbose = verbose, z.list = z.list, miss.init = miss.init, miss.list = miss.list, miss.hc = miss.hc, asMclust = asMclust, posidens = posidens)
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
#' \item{\code{linf}}{The most current estimate of the final converged maximised log-likelihood.}
#' \item{\code{a}}{The Aitken acceleration value where typically \code{0 <= a <= 1}. When \code{a < 0}, a numerical issue or bug has occurred; when \code{a > 1}, the algorithm is accelerating and should not be stopped.}
#' When the \code{"aitken"} method is employed within \code{\link{MoE_clust}} (via \code{\link{MoE_control}}), \code{ll} at convergence gives the log-likelihood achieved by the estimated parameters, while \code{linf} at convergence estimates the log-likelihood that would be achieved after an infinite number of EM/CEM iterations.
#' @export
#' @keywords control
#' @author Keefe Murphy - <\email{keefe.murphy@@mu.ie}>
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
#' Takes one or more sets of MoEClust models fitted by \code{\link{MoE_clust}} (or \code{\link{MoE_stepwise}}) and ranks them according to the BIC, ICL, or AIC. It's possible to respect the internal ranking within each set of models, or to discard models within each set which were already deemed sub-optimal. This function can help with model selection via exhaustive or stepwise searches.
#' @param ... One or more objects of class \code{"MoEClust"} outputted by \code{\link{MoE_clust}}. All models must have been fit to the same data set. A single \emph{named} list of such objects can also be supplied. Additionally, objects of class \code{"MoECompare"} outputted by this very function can also be supplied here.
#' 
#' This argument is only relevant for the \code{\link{MoE_compare}} function and will be ignored for the associated \code{print} function.
#' @param criterion The criterion used to determine the ranking. Defaults to \code{"bic"}.
#' @param pick The (integer) number of models to be ranked and compared. Defaults to \code{10L}. Will be constrained by the number of models within the \code{"MoEClust"} objects supplied via \code{...} if \code{optimal.only} is \code{FALSE}, otherwise constrained simply by the number of \code{"MoEClust"} objects supplied. Setting \code{pick=Inf} is a valid way to select all models.
#' @param optimal.only Logical indicating whether to only rank models already deemed optimal within each \code{"MoEClust"} object (\code{TRUE}), or to allow models which were deemed suboptimal enter the final ranking (\code{FALSE}, the default). See \code{details}.
#' @param x,index,posidens,rerank,digits,details,maxi Arguments required for the associated \code{print} function:
#' \describe{
#' \item{\code{x}}{An object of class \code{"MoECompare"} resulting from a call to \code{\link{MoE_compare}}.}
#' \item{\code{index}}{A logical or numeric vector giving the indices of the rows of the table of ranked models to print. This defaults to the full set of ranked models. It can be useful when the table of ranked models is large to examine a subset via this \code{index} argument, for display purposes. See \code{rerank}.}
#' \item{\code{posidens}}{A logical indicating whether models which have been flagged for having positive log-densities should be included in the comparison (defaults to \code{TRUE}). Such models may correspond to spurious solutions and can be discarded by specifying \code{posidens=FALSE}. Only relevant if any of the \code{"MoEClust"} objects being compared were themselves run with \code{posidens=TRUE}.}
#' \item{\code{rerank}}{A logical indicating whether the ranks should be recomputed when subsetting using \code{index}. Defaults to \code{FALSE}. Only relevant when \code{details=TRUE}.}
#' \item{\code{digits}}{The number of decimal places to round model selection criteria to (defaults to 3).}
#' \item{\code{details}}{Logical indicating whether some additional details should be printed, defaults to \code{TRUE}. Exists to facilitate \code{\link{MoE_stepwise}} printing.}
#' \item{\code{maxi}}{A number specifying the maximum number of rows/models to print. Defaults to \code{length(index)}.}}
#' @note The \code{criterion} argument here need not comply with the criterion used for model selection within each \code{"MoEClust"} object, but be aware that a mismatch in terms of \code{criterion} \emph{may} require the optimal model to be re-fit in order to be extracted, thereby slowing down \code{\link{MoE_compare}}.
#' 
#' If random starts had been used via \code{init.z="random"} the \code{optimal} model may not necessarily correspond to the highest-ranking model in the presence of a criterion mismatch, due to the randomness of the initialisation. 
#'
#' A dedicated \code{print} function exists for objects of class \code{"MoECompare"}.
#' 
#' \code{\link{plot.MoEClust}} and \code{\link[=as.Mclust.MoEClust]{as.Mclust}} can both also be called on objects of class \code{"MoECompare"}.
#'
#' @details The purpose of this function is to conduct model selection on \code{"MoEClust"} objects, fit to the same data set, with different combinations of gating/expert network covariates or different initialisation settings.
#'
#' Model selection will have already been performed in terms of choosing the optimal number of components and GPCM/\pkg{mclust} model type within each supplied set of results, but \code{\link{MoE_compare}} will respect the internal ranking of models when producing the final ranking if \code{optimal.only} is \code{FALSE}: otherwise only those models already deemed optimal within each \code{"MoEClust"} object will be ranked.
#'
#' As such if two sets of results are supplied when \code{optimal.only} is \code{FALSE}, the 1st, 2nd and 3rd best models could all belong to the first set of results, meaning a model deemed suboptimal according to one set of covariates could be superior to one deemed optimal under another set of covariates.
#' @return A list of class \code{"MoECompare"}, for which a dedicated print function exists, containing the following elements (each of length \code{pick}, and ranked according to \code{criterion}, where appropriate):
#' \item{\code{data}}{The name of the data set to which the models were fitted.}
#' \item{\code{optimal}}{The single optimal model (an object of class \code{"MoEClust"}) among those supplied, according to the chosen \code{criterion}.}
#' \item{\code{pick}}{The final number of ranked models. May be different (i.e. less than) the supplied \code{pick} value.}
#' \item{\code{MoENames}}{The names of the supplied \code{"MoEClust"} objects.}
#' \item{\code{modelNames}}{The \code{\link[mclust]{mclustModelNames}}.}
#' \item{\code{G}}{The optimal numbers of components.}
#' \item{\code{df}}{The numbers of estimated parameters.}
#' \item{\code{iters}}{The numbers of EM/CEM iterations.}
#' \item{\code{bic}}{BIC values, ranked according to \code{criterion}.}
#' \item{\code{icl}}{ICL values, ranked according to \code{criterion}.}
#' \item{\code{aic}}{AIC values, ranked according to \code{criterion}.}
#' \item{\code{loglik}}{Maximal log-likelihood values, ranked according to \code{criterion}.}
#' \item{\code{gating}}{The gating formulas.}
#' \item{\code{expert}}{The expert formulas.}
#' \item{\code{algo}}{The algorithm used for fitting the model - either \code{"EM"}, \code{"CEM"}, \code{"cemEM"}.}
#' \item{\code{equalPro}}{Logical indicating whether mixing proportions were constrained to be equal across components.}
#' \item{\code{hypvol}}{Hypervolume parameters for the noise component if relevant, otherwise set to \code{NA} (see \code{\link{MoE_control}}).}
#' \item{\code{noise}}{The type of noise component fitted (if any). Only displayed if at least one of the compared models has a noise component.}
#' \item{\code{noise.gate}}{Logical indicating whether gating covariates were allowed to influence the noise component's mixing proportion. Only printed for models with a noise component, when at least one of the compared models has gating covariates.}
#' \item{\code{equalNoise}}{Logical indicating whether the mixing proportion of the noise component for \code{equalPro} models is also equal (\code{TRUE}) or estimated (\code{FALSE}).}
#' @export
#' @keywords clustering main
#' @references Murphy, K. and Murphy, T. B. (2020). Gaussian parsimonious clustering models with covariates and a noise component. \emph{Advances in Data Analysis and Classification}, 14(2): 293-325. <\doi{10.1007/s11634-019-00373-8}>.
#' @importFrom mclust "mclustModelNames"
#' @author Keefe Murphy - <\email{keefe.murphy@@mu.ie}>
#'
#' @seealso See \code{\link{MoE_stepwise}} for identifying the optimal model and its covariates via greedy forward stepwise selection.\cr
#' 
#' \code{\link{MoE_clust}}, \code{\link[mclust]{mclustModelNames}}, \code{\link{plot.MoEClust}}, \code{\link[=as.Mclust.MoEClust]{as.Mclust}}
#' @usage
#' MoE_compare(...,
#'             criterion = c("bic", "icl", "aic"),
#'             pick = 10L,
#'             optimal.only = FALSE)
#' @examples
#' \donttest{data(CO2data)
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
#' m7    <- MoE_clust(CO2, G=2:3, expert= ~ GNP, tau0=0.1)
#'
#' # Rank only the optimal models and examine the best model
#' (comp <- MoE_compare(m1, m2, m3, m4, m5, m6, m7, optimal.only=TRUE))
#' (best <- comp$optimal)
#' (summ <- summary(best, classification=TRUE, parameters=TRUE, networks=TRUE))
#'
#' # Examine all models visited, including those already deemed suboptimal
#' # Only print models with expert covariates & more than one component
#' comp2 <- MoE_compare(m1, m2, m3, m4, m5, m6, m7, pick=Inf)
#' print(comp2, index=comp2$expert != "None" & comp2$G > 1)
#' 
#' # Conduct a stepwise search on the same data
#' (mod1 <- MoE_stepwise(CO2, GNP))
#' 
#' # Conduct another stepwise search considering models with a noise component
#' (mod2 <- MoE_stepwise(CO2, GNP, noise=TRUE))
#' 
#' # Compare both sets of results to choose the optimal model
#' (best <- MoE_compare(mod1, mod2, optimal.only=TRUE)$optimal)}
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
      dat.name    <- if(any(comparison)) dots[[1L]]$data
      dots[comparison]          <- sapply(dots[comparison], "[", "optimal")
      MoEs        <- dots[mod.names]
      if(is.null(mod.names))                      stop("When supplying models as a list, every element of the list must be named", call.=FALSE)
    } else {
      dots        <- list(...)
      mod.names   <- vapply(call, deparse, character(1L))
      comparison  <- sapply(dots, inherits, "MoECompare", logical(1L))
      dat.name    <- if(any(comparison)) dots[[1L]]$data
      dots[comparison]          <- sapply(dots[comparison], "[", "optimal")
      MoEs        <- stats::setNames(dots, mod.names)
      mod.names   <- unique(mod.names)
      MoEs        <- MoEs[mod.names]
    }
    Mclass        <- vapply(MoEs, class,          character(1L))
    if(any(Mclass != "MoEClust"))                 stop("All models must be of class 'MoEClust'!", call.=FALSE)
    data          <- lapply(MoEs,  "[[", "data")
    data          <- lapply(data, unname)
    if(length(data) > 1   && !.unique_list(data)) stop("All models being compared must have been fit to the same data set!", call.=FALSE)
    dat.name      <- if(is.null(dat.name)) deparse(MoEs[[1L]]$call$data) else dat.name
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
    PDxs          <- lapply(LLxs, attr, "posidens")
    choice        <- max(lengths(BICs))
    bics          <- lapply(BICs, function(x) .pick_MoECrit(x, choice)$crits)
    icls          <- lapply(ICLs, function(x) .pick_MoECrit(x, choice)$crits)
    aics          <- lapply(AICs, function(x) .pick_MoECrit(x, choice)$crits)
    llxs          <- lapply(LLxs, function(x) .pick_MoECrit(x, choice)$crits)
    dfxs          <- lapply(DFxs, function(x) .pick_MoECrit(x, choice)$crits)
    itxs          <- lapply(ITxs, function(x) .pick_MoECrit(x, choice)$crits)
    pdxs          <- lapply(PDxs, .pick_posidens)
    if(optimal.only) {
      opt.names   <- names(.crits_names(lapply(switch(EXPR=criterion, bic=bics, icl=icls, aic=aics), "[", 1L)))
    }
    bics          <- .crits_names(bics)
    icls          <- .crits_names(icls)
    aics          <- .crits_names(aics)
    llxs          <- .crits_names(llxs)
    dfxs          <- .crits_names(dfxs)
    itxs          <- .crits_names(itxs)
    pdxs          <- .crits_names(pdxs)
    if(optimal.only) {
      bics        <- bics[names(bics) %in% opt.names]
      icls        <- icls[names(icls) %in% opt.names]
      aics        <- aics[names(aics) %in% opt.names]
      llxs        <- llxs[names(llxs) %in% opt.names]
      dfxs        <- dfxs[names(dfxs) %in% opt.names]
      itxs        <- itxs[names(itxs) %in% opt.names]
      pdxs        <- pdxs[names(pdxs) %in% opt.names]
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
    modelNames    <- replace(gsub(",.*", "", gsub(".*\\|", "", max.names)), G == 0, "")
    best.model    <- MoEs[[crit.names[1L]]]
    if(best.model$modelName     != modelNames[1L] || best.model$G != G[1L]) {
      bestGN      <- gating[1L] != "~1" && (best.model$G + !noise.null[crit.names[1L]] - isFALSE(noise.gate[crit.names[1L]])) <= 1
      best.model$net.covs       <- if(bestGN) attr(best.model$net.covs, "Discarded") else best.model$net.covs
      message("Re-fitting optimal model due to mismatched 'criterion'...\n\n")
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
        attr(best.model, "Criterion")  <- criterion
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
    comp          <- list(data = dat.name, optimal = best.model, pick = pick, MoENames = crit.names, modelNames = modelNames, G = as.integer(G), 
                          df = as.integer(unname(dfxs[max.names])), iters = as.integer(unname(itxs[max.names])), bic = unname(bics[max.names]), 
                          icl = unname(icls[max.names]), aic = unname(aics[max.names]), loglik = unname(llxs[max.names]), posidens = as.logical(unname(pdxs[max.names])), 
                          gating = gating, expert = expert, algo = unname(algo[crit.names]), equalPro = equalPro, hypvol = unname(hypvol), noise = unname(noise.meth), 
                          noise.gate = unname(replace(noise.gate, gating == "None" | G <= 1L - !is.na(hypvol), NA)), equalNoise = unname(replace(equalNoise, !equalPro | is.na(equalPro), NA)))
    if(any(comp$posidens))                        warning("Potentially spurious solutions with positive log-densities are included in the comparison\n", call.=FALSE)
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
#' Predicts both cluster membership probabilities and fitted response values from a \code{MoEClust} model, using covariates and response data, or covariates only. The predicted MAP classification, mixing proportions, and component means are all also reported in both cases, as well as the predictions of the expert network corresponding to the most probable component.
#' @param object An object of class \code{"MoEClust"} generated by \code{\link{MoE_clust}}, or an object of class \code{"MoECompare"} generated by \code{\link{MoE_compare}}. Predictions for models with a noise component are facilitated here too (see \code{discard.noise}).
#' @param newdata A list with two \emph{named} components, each of which must be a \code{data.frame} or \code{matrix} with named columns, giving the data for which predictions are desired.
#' \describe{
#' \item{\code{new.x}}{The new covariates for the \code{gating} &/or \code{expert} networks. \strong{Must} be supplied when \code{newdata$new.y} is supplied.}
#' \item{\code{new.y}}{(Optional) response data (see \code{use.y} below). When supplied, cluster and response prediction is based on both \code{newdata$new.x} and \code{newdata$new.y}, otherwise only on the covariates in \code{newdata$new.x}.}
#' }
#' If supplied as a list with elements \code{new.x} and \code{new.y}, both \strong{must} have the same number of rows.
#'
#' Alternatively, a single \code{data.frame} or \code{matrix} can be supplied and an attempt will be made to extract & separate covariate and response columns (\emph{if any}) into \code{newdata$new.x} and \code{newdata$new.y} based on the variable names in \code{object$data} and \code{object$net.covs}.
#'
#' When \code{newdata} is not supplied in any way, the covariates and response variables used in the fitting of the model are used here. It is possible to not supply \code{new.y} and to supply an empty \code{data.frame} or \code{matrix} for \code{new.x} (or to equivalently supply an empty \code{data.frame} or \code{matrix} for \code{newdata} itself) for models with no covariates of any kind, which effectively predicts the weighted mean of the component means.
#' @param resid A logical indicating whether to return the residuals also. Defaults to \code{FALSE}. Only allowed when response variables are supplied in some form. The function \code{residuals} is a wrapper to \code{predict} with the argument \code{resid} set to \code{TRUE}, with only the residuals returned.
#' @param discard.noise A logical governing how predictions of the responses are made for models with a noise component (otherwise this argument is irrelevant). By default (\code{FALSE}), the mean of the noise component is accounted for. Otherwise, or when the mean of the noise component is unavailable (due to having been manually supplied through \code{\link{MoE_control}} via \code{noise.args$noise.vol}), prediction of the responses is performed using a \code{z} matrix which is renormalised after discarding the column corresponding to the noise component. The renormalisation approach can be forced by specifying \code{TRUE}, even when the mean of the noise component is available. For models with a noise component fitted with \code{algo="CEM"}, a small extra E-step is conducted for observations assigned to the non-noise components in this case.
#' @param MAPresids A logical indicating whether residuals are computed against \code{y} (\code{TRUE}, the default) or \code{MAPy} when \code{FALSE}. Not relevant for models with equal mixing proportions when only \code{new.x} is available. See \strong{Value} below for more details.
#' @param use.y A logical indicating whether the response variables (if any are supplied either via \code{new.y} or via \code{newdata} itself) are actually used in the prediction. Defaults to \code{TRUE}, but useful when \code{FALSE} for computing residuals as though only the covariates in \code{new.x} were supplied. For out-of-sample prediction, typically \code{new.y} would not be supplied anyway and so the \code{use.y=TRUE} default becomes irrelevant.
#' @param ... Catches unused arguments (and allows the \code{predict} arguments \code{discard.noise} &/or \code{use.y} to be passed through \code{fitted} or the \code{discard.noise}, \code{MAPresids}, and/or \code{use.y} arguments to be passed through \code{residuals}).
#'
#' @return A list with the following named components, regardless of whether \code{newdata$new.x} and \code{newdata$new.y} were used, or \code{newdata$new.x} only.
#' \item{\code{y}}{Aggregated fitted values of the response variables.}
#' \item{\code{z}}{A matrix whose \code{[i,k]}-th entry is the probability that observation \emph{i} of the \code{newdata} belongs to the \emph{k}-th component. For models with a noise component, the final column gives the probability of belonging to the so-called \emph{Cluster0}.}
#' \item{\code{classification}}{The vector of predicted cluster labels for the \code{newdata}. \code{0} is returned for observations assigned to the noise component.}
#' \item{\code{pro}}{The predicted mixing proportions for the \code{newdata}, i.e. predicted values of the gating network. \code{object$parameters$pro} is returned for models without gating network covariates. See \code{\link{predict.MoE_gating}}.}
#' \item{\code{mean}}{The predicted component means for the \code{newdata}, i.e. predicted values of the expert network. Given as a 3-dimensional array with dimensions given by the number of new observations, the number of variables, and the number of clusters. The first dimension is of length \code{1} when there are no expert network covariates, in which case the entries correspond to \code{object$parameters$mean}. See \code{\link{predict.MoE_expert}}.}
#' \item{\code{MAPy}}{Fitted values of the single expert network to which each observation is most probably assigned. Not returned for models with equal mixing proportions when only \code{new.x} is available. Likely to only be of use for models with gating and expert covariates when only \code{new.x} is supplied. Note that \code{MAPy} and \code{y} will coincide for models fitted via the CEM algorithm (see \code{\link{MoE_control}} and its argument \code{algo}).}
#'
#' When \code{residuals} is called, only the residuals (governed by \code{MAPresids}) are returned; when \code{predict} is called with \code{resid=TRUE}, the list above will also contain the element \code{resids}, containing the residuals.
#' 
#' The returned values of \code{pro} and \code{mean} are always the same, regardless of whether \code{newdata$new.x} and \code{newdata$new.y} were used, or \code{newdata$new.x} only.
#' 
#' Finally, \code{fitted} is simply a wrapper to \code{predict.MoEClust(object)$y} without any \code{newdata}, and with the \code{resid} and \code{MAPresids} arguments also ignored.
#'
#' @details Predictions can also be made for models with a noise component, in which case \code{z} will include the probability of belonging to \code{"Cluster0"} & \code{classification} will include labels with the value \code{0} for observations classified as noise (if any). The argument \code{discard.noise} governs how the responses are predicted in the presence of a noise component (see \code{\link{noise_vol}} for more details).
#' 
#' Note that the argument \code{discard.noise} is invoked for any models with a noise component, while the similar \code{\link{MoE_control}} argument \code{noise.args$discard.noise} is only invoked for models with both a noise component and expert network covariates.
#' 
#' Please be aware that a model considered optimal from a clustering point of view may not necessarily be optimal from a prediction point of view. In particular, full MoE models with covariates in both networks (for which both the cluster membership probabilities and component means are observation-specific) are recommended for out-of-sample prediction when only new covariates are observed (see \code{new.x} and \code{new.y} above, as well as \code{use.y}).
#' @note Note that a dedicated \code{\link[=predict.MoE_gating]{predict}} function is also provided for objects of class \code{"MoE_gating"} (typically \code{object$gating}, where \code{object} is of class \code{"MoEClust"}). This function is effectively a shortcut to \code{predict(object, ...)$pro}, which (unlike the \code{predict} method for \code{\link[nnet]{multinom}} on which it is based) accounts for the various ways of treating gating covariates and noise components, although its \code{type} argument defaults to \code{"probs"} rather than \code{"class"}. Notably, its \code{keep.noise} argument behaves differently from the \code{discard.noise} argument here; here, the noise component is \strong{only} discarded in the computation of the predicted responses. See \code{\link{predict.MoE_gating}} for further details.
#' 
#' Similarly, a dedicated \code{\link[=predict.MoE_expert]{predict}} function is also provided for objects of class \code{"MoE_expert"} (typically \code{object$expert}, where \code{object} is of class \code{"MoE_expert"}). This function is effectively a wrapper to \code{predict(object, ...)$mean}, albeit it returns a list (by default) rather than a 3-dimensional array and also \emph{always} preserves the dimensions of \code{newdata}, even for models without expert network covariates. See \code{\link{predict.MoE_expert}} for further details. 
#' @references Murphy, K. and Murphy, T. B. (2020). Gaussian parsimonious clustering models with covariates and a noise component. \emph{Advances in Data Analysis and Classification}, 14(2): 293-325. <\doi{10.1007/s11634-019-00373-8}>.
#' @author Keefe Murphy - <\email{keefe.murphy@@mu.ie}>
#' @seealso \code{\link{MoE_clust}}, \code{\link{MoE_control}}, \code{\link{noise_vol}}, \code{\link{predict.MoE_gating}}, \code{\link{predict.MoE_expert}}
#' @method predict MoEClust
#' @keywords prediction main
#' @importFrom matrixStats "rowSums2"
#' @export
#' @usage
#' \method{predict}{MoEClust}(object,
#'         newdata,
#'         resid = FALSE,
#'         discard.noise = FALSE,
#'         MAPresids = FALSE,
#'         use.y = TRUE,
#'         ...)
#' @examples
#' data(ais)
#' # Fit a MoEClust model and predict the same data
#' res     <- MoE_clust(ais[,3:7], G=2, gating= ~ BMI, expert= ~ sex,
#'                      modelNames="EVE", network.data=ais)
#' pred1   <- predict(res)
#' 
#' # Get only the fitted responses
#' fits    <- fitted(res)
#' all.equal(pred1$y, fits) #TRUE
#'
#' # Remove some rows of the data for prediction purposes
#' ind     <- sample(1:nrow(ais), 5)
#' dat     <- ais[-ind,]
#'
#' # Fit another MoEClust model to the retained data
#' res2    <- MoE_clust(dat[,3:7], G=3, gating= ~ BMI + sex,
#'                      modelNames="EEE", network.data=dat)
#'
#' # Predict held back data using the covariates & response variables
#' (pred2  <- predict(res2, newdata=ais[ind,]))
#' # pred2 <- predict(res2, newdata=list(new.y=ais[ind,3:7],
#' #                                     new.x=ais[ind,c("BMI", "sex")]))
#'
#' # Get the residuals
#' residuals(res2, newdata=ais[ind,])
#'
#' # Predict held back data using only the covariates
#' (pred3  <- predict(res2, newdata=ais[ind,], use.y=FALSE))
#' # pred3 <- predict(res2, newdata=list(new.x=ais[ind,c("BMI", "sex")]))
#' # pred3 <- predict(res2, newdata=ais[ind,c("BMI", "sex")])
predict.MoEClust  <- function(object, newdata = list(...), resid = FALSE, discard.noise = FALSE, 
                              MAPresids = FALSE, use.y = TRUE, ...) {
  object          <- if(inherits(object, "MoECompare")) object$optimal else object
  dcard           <- discard.noise
  if(length(resid) > 1   || !is.logical(resid))   stop("'resid' should be a single logical indicator",         call.=FALSE)
  if(length(dcard) > 1   || !is.logical(dcard))   stop("'discard.noise' should be a single logical indicator", call.=FALSE)
  if(length(MAPresids)    > 1  || 
     !is.logical(MAPresids))                      stop("'MAPresids' should be a single logical indicator",     call.=FALSE)
  if(length(use.y) > 1   || !is.logical(use.y))   stop("'use.y' should be a single logical indicator",         call.=FALSE)
  algo            <- attr(object, "Algo")
  net             <- object$net.covs
  dat             <- object$data
  datnames        <- colnames(dat)
  nc              <- length(datnames)
  hypvol          <- object$hypvol
  noise           <- !is.na(hypvol)
  yM              <- FALSE
  dot             <- list(...)
  if(any(names(dot) == "MAPWARN"))     {
    MAPWARN       <- dot$MAPWARN
    dot           <- dot[names(dot)         != "MAPWARN"]
    newdata       <- newdata[names(newdata) != "MAPWARN"]
    newdata       <- if(length(newdata)     != 0) newdata
  } else MAPWARN  <- TRUE
  if(any(!(names(dot) %in% c("new.x", "new.y")))) stop("Invalid arguments passed through '...' construct", call.=FALSE)
  if(nmiss        <- ifelse(inherits(newdata,
                                     "list")    &&
                     length(newdata)  == 0,
                     missing(newdata), is.null(newdata))) {
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
      if(ncol(newdata.x)  > 0  &&
         is.null(colnames(newdata.x)))            stop("'newdata$new.x' must have named columns", call.=FALSE)
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
  } else           {
    if(!is.matrix(newdata)     &&
       !is.data.frame(newdata))                   stop("'newdata' must be either be a 'list', 'matrix', or 'data.frame'",        call.=FALSE)
    newdata.x     <- newdata.y <- newdata
  }
  Xexp            <- attr(object, "Expert")
  Xgat            <- attr(object, "Gating")
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
     levels(net[,gatenames[p]])), logical(1L))))  warning("One or more categorical gating covariates in the unseen newdata has new factor levels\n", call.=FALSE, immediate.=TRUE)
  if(any(vapply(seq_len(ncol(newexpx)),
     function(p, ep=newexpx[,p])  is.factor(ep) &&
     !identical(levels(ep),
     levels(net[,expxnames[p]])), logical(1L))))  warning("One or more categorical expert covariates in the unseen newdata has new factor levels\n", call.=FALSE, immediate.=TRUE)
  rownames(newdata.x)          <- NULL
  nr              <- nrow(newdata.x)
  nrseq           <- seq_len(nr)
  G               <- object$G
  GN              <- G + noise
  Gseq            <- seq_len(G)
  params          <- object$parameters
  if(!yM)          {
    newdata.y     <- as.data.frame(newdata.y)
    y.names       <- names(newdata.y)
    newdata.y     <- newdata.y[cnew,y.names %in% datnames,  drop=FALSE]
    if(!all(datnames %in% y.names))  {
      if(nL)       {                              stop("Response variables missing in 'newdata'",      call.=FALSE)
      } else       {                              warning("Response variables missing in 'newdata'\n", call.=FALSE, immediate.=TRUE)
        yM        <- TRUE
      }
    } else         {
      newdata.y   <- newdata.y[,datnames,                   drop=FALSE]
      rownames(newdata.y)      <- NULL
    }
  }
  noise.loc       <- attr(hypvol, "Location")
  noise.loc       <- if(is.null(noise.loc)) rep(NA, nc) else noise.loc
  if(G == 0)   {
    if(all(is.na(noise.loc)))                     warning("Can't predict the response; mean of noise component unavailable\n", call.=FALSE, immediate.=TRUE)
    retval        <- list(ystar=matrix(noise.loc, nrow=nr, ncol=object$d, byrow=TRUE),
                          classification=rep(0L, nr),
                          pro=provideDimnames(matrix(1L,   nrow=1L, ncol=1L), base=list("pro", "Cluster0")),
                          zstar=provideDimnames(matrix(1L, nrow=nr, ncol=1L), base=list(as.character(nrseq), "Cluster0")))
    retval        <- c(retval, list(MAPy=retval$ystar))
    retval        <- if(isTRUE(resid)) c(retval, list(resids=newdata.y - retval$ystar)) else retval
    class(retval) <- "listof"
      return(retval)
  }
  pred.exp        <- predict.MoE_expert(object$expert, newdata=newexpx, simplify=FALSE, droplevels=FALSE)
  new.exp         <- 
  mus             <- NULL 
  if(isTRUE(nmiss))       {
    new.tau       <- if(isTRUE(Xgat))  params$pro else matrix(params$pro, nrow=nrow(dat), ncol=length(params$pro), byrow=TRUE)
    zstar         <- if(isTRUE(use.y))   object$z else new.tau
  } else       {
    if(isFALSE(nmiss))            {
      new.tau     <- predict.MoE_gating(object$gating, newdata=newgate, type="probs", keep.noise=TRUE, droplevels=FALSE)
    }
    if(resid  && !(resid <- !yM))                 warning("'resid' can only be TRUE when response variables are supplied\n", call.=FALSE, immediate.=TRUE)
    attr(net, "Gating")  <-
    attr(net, "Expert")  <-
    attr(net, "Both")    <-
    rownames(net)        <-
    rownames(dat)        <- NULL
    if(identical(newdata.x, net)      &&
       (!yM   && use.y)  && 
       identical(newdata.y, dat))      {
      zstar       <- object$z
      nmiss       <- TRUE
    } else if(yM  || !use.y)           {
      zstar       <- new.tau
    } else     {
      new.exp     <- if(Xexp) lapply(pred.exp, "-", newdata.y)    else newdata.y
      mus         <- if(Xexp) 0L                                  else params$mean
      zstar       <- MoE_estep(Dens=MoE_dens(data=new.exp, mus=mus, sigs=params$variance, log.tau=log(new.tau), Vinv=params$Vinv))$z
    }
  }
  if(ZD           <- !noise || attr(hypvol, "Meth") == "manual" || isTRUE(dcard))  {
    zstar2        <- if(noise) .renorm_z(zstar[,-GN, drop=FALSE]) else zstar
    if(algo == "CEM"     &&
       any(nan           <- apply(zstar2, 1L, function(x) all(is.nan(x)))))    {
      new.exp     <- if(!is.null(new.exp)) new.exp else if(Xexp) lapply(pred.exp, "-", newdata.y)  else newdata.y
      mus         <- if(!is.null(mus))     mus     else if(Xexp) 0L                                else params$mean
      zstar2[nan,]       <- .renorm_z(MoE_estep(data=lapply(new.exp, "[", nan, TRUE), mus=mus, sigs=params$variance, 
                                      log.tau=log(if(is.matrix(new.tau)) new.tau[nan,, drop=FALSE] else new.tau), Vinv=params$Vinv)$z[,-GN, drop=FALSE])
    }
    if(any(nan           <- apply(zstar2, 1L, function(x) all(is.nan(x)))))    {
      zstar2[nan,]       <- 1/G
    }
    ystar         <- as.matrix(Reduce("+", lapply(Gseq, function(g) zstar2[,g] * pred.exp[[g]])))
  } else     {
    ystar         <- as.matrix(Reduce("+", lapply(Gseq, function(g) zstar[,g]  * pred.exp[[g]])))
    ystar         <- ystar + zstar[,GN] *  matrix(noise.loc, nrow=nr, ncol=object$d, byrow=TRUE)
  }
  claX            <- stats::setNames(max.col(zstar), nrseq)
  claX[claX   == G + 1]  <- 0L
  dimnames(ystar) <- list(nrseq, datnames)
  gnames          <- paste0("Cluster", Gseq)
  gnames0         <- if(noise) c(gnames, "Cluster0") else gnames
  dimnames(zstar) <- list(nrseq, gnames0)
  mu              <- if(Xexp) array(do.call(cbind, pred.exp), dim=c(nr, nc, G))      else array(params$mean, dim=c(1L, nc, G))
  dimnames(mu)    <- list(if(Xexp) nrseq             else "mean", datnames, gnames)
  if(Xgat)         {
    dimnames(new.tau)    <- list(nrseq, gnames0)
  } else           {
    new.tau       <- provideDimnames(if(is.matrix(new.tau)) new.tau[1L,, drop=FALSE] else t(new.tau), base=list("pro", gnames0))
  }
  retval          <- list(y=ystar, classification=claX, z=zstar, pro=new.tau, mean=mu)
  if((isTRUE(yM)  || isFALSE(use.y))  && 
     attr(object, "EqualPro"))         {
    if(isTRUE(MAPWARN))                {           
      warnMAP     <- "'MAPy' not available due to randomness induced by the equal mixing proportions and the absent response variables\n"
      if(isTRUE(MAPresids)) {                      
        MAPresids <- FALSE
        warnMAP   <- paste0(warnMAP, "\nHence, 'MAPresids' is forced to FALSE\n")
      }
                                                  warning(warnMAP, call.=FALSE, immediate.=TRUE)
    }
  } else           {
    MAPy          <- ystar
    cD            <- if(ZD && noise) max.col(zstar2) else claX
    MAPy[]        <- t(vapply(seq_len(nr), function(i, z=cD[i]) if(z == 0) noise.loc else mu[ifelse(Xexp, i, 1L),,z], numeric(nc)))
    retval        <- c(retval, list(MAPy = MAPy))
  }
  retval          <- if(isTRUE(resid)) c(retval, list(resids=provideDimnames(as.matrix(newdata.y - if(isTRUE(MAPresids)) MAPy else ystar), 
                                         base=list(as.character(nrseq), datnames)))) else retval
  class(retval)   <- "listof"
    return(retval)
}

#' @rdname predict.MoEClust
#' @method fitted MoEClust
#' @keywords prediction utility
#' @importFrom matrixStats "rowSums2"
#' @usage 
#' \method{fitted}{MoEClust}(object,
#'        ...)
#' @export
  fitted.MoEClust        <- function(object, ...) {
    args   <- c(list(object=object, newdata=NULL), as.list(match.call())[-1L])
    fits   <- do.call(predict.MoEClust, args[unique(names(args))])$y
      if(!is.null(fits))    return(fits)
  }

#' @rdname predict.MoEClust
#' @method residuals MoEClust
#' @keywords prediction utility
#' @importFrom matrixStats "rowSums2"
#' @usage
#' \method{residuals}{MoEClust}(object,
#'           newdata,
#'           ...)
#' @export
  residuals.MoEClust     <- function(object, newdata = list(...), ...) {
    dots   <- list(...)
    MAPW   <- any(names(dots) == "MAPresids") && isTRUE(dots$MAPresids)
    args   <- c(list(object=object, resid=TRUE, MAPWARN=MAPW),  as.list(match.call())[-1L])
    resids <- do.call(predict.MoEClust, args[unique(names(args))])$resids
      if(!is.null(resids))  return(resids)
  }
  
#' Predictions from MoEClust expert networks
#' 
#' Predictions (point estimates) of observation-specific component means from each (non-noise) component's expert network linear regression.
#' @param object An object of class \code{"MoE_expert"} (typically \code{x$expert}, where \code{x} is of class \code{"MoEClust"}).
#' @param newdata A matrix or data frame of test examples. If omitted, the fitted values are used.
#' @param simplify Logical indicating whether to simplify the output (in the form of a list) to a 3-dimensional array with dimensions given by the number of new observations, the number of variables, and the number of clusters. The first dimension of such an array is of length \code{1} when there are no expert network covariates, in which case the entries correspond to \code{object$parameters$mean}. Defaults to \code{FALSE}.
#' @param droplevels A logical indicating whether unseen factor levels in categorical variables within \code{newdata} should be dropped (with \code{NA} predicted in their place). Defaults to \code{FALSE}. See \code{\link{drop_levels}}.
#' @param ... Catches unused arguments or allows the \code{simplify} argument to be passed through \code{fitted} and \code{residuals}.
#' 
#' @return For \code{simplify=FALSE}, either a list of vectors or predictions (for univariate data) or a list of matrices of predictions (for multivariate data). These lists are of the same length as number of non-noise components in the fitted model. When \code{simplify=TRUE}, a 3-dimensional array of predictions is returned, with respective dimensions given by the number of observations, variables, and non-noise components.
#' @details This function is effectively just a shortcut to \code{lapply(x$expert, predict.lm, newdata=...)}. It can also be thought of as a wrapper to \code{\link{predict.MoEClust}(x, ...)$mean}, although it returns a list (by default) rather than a 3-dimensional array and also \emph{always} preserves the dimensions of \code{newdata}, even for models without expert network covariates.
#' @author Keefe Murphy - <\email{keefe.murphy@@mu.ie}>
#' @seealso \code{\link{predict.MoEClust}}, \code{\link[stats]{lm}}, \code{\link{predict.MoE_gating}}, \code{\link{drop_levels}}
#' @method predict MoE_expert
#' @keywords prediction utility
#' @export
#' @usage 
#' \method{predict}{MoE_expert}(object,
#'         newdata = NULL,
#'         simplify = FALSE,
#'         droplevels = FALSE,
#'         ...)
#' @examples 
#' data(CO2data)
#' res <- MoE_clust(CO2data$CO2, G=3, equalPro=TRUE, expert= ~ GNP, network.data=CO2data)
#' predict(res$expert)
#' 
#' # Try with newdata and simplify=TRUE
#' predict(res$expert, newdata=CO2data[1:5,"GNP", drop=FALSE], simplify=TRUE)
  predict.MoE_expert     <- function(object, newdata = NULL, simplify=FALSE, droplevels = FALSE, ...) {
    if(!is.null(newdata) &&
       all(!is.matrix(newdata), 
           !is.data.frame(newdata)))               stop("'newdata' must be a matrix or data frame if supplied", call.=FALSE)
    if(length(simplify)   > 1    ||
       !is.logical(simplify))                      stop("'simplify' must be a single logical indicator",        call.=FALSE)
    if(length(droplevels) > 1    ||
       !is.logical(droplevels))                    stop("'droplevels' must be a single logical indicator",      call.=FALSE)
    fits   <- lapply(object, "[[", "fitted.values")
    nr     <- ifelse(is.null(newdata), nrow(object[[1]]$model), NROW(newdata))
    G      <- length(object)
    d      <- attr(object, "d")
    gnames <- paste0("Cluster", seq_along(fits))
    dnames <- if(d == 1) names(fits[[1L]]) else dimnames(fits[[1L]])
    dnames <- if(is.null(newdata))  dnames else as.character(seq_len(NROW(newdata)))
    exp    <- attr(object, "Formula")
    if(exp == "~1")       {
      mus  <- attr(object, "NoExp")
      if(d == 1) {
        fits    <- lapply(seq_len(G), function(g) stats::setNames(rep(mus[g], nr), dnames))
      } else     {
        fits    <- lapply(seq_len(G), function(g) provideDimnames(matrix(mus[,g], nrow=nr, ncol=d, byrow=TRUE), base=dnames))
      }
      fits      <- stats::setNames(fits, gnames)
    } else if(!is.null(newdata))  {
      newdata   <- if(isTRUE(droplevels)) drop_levels(object[[1L]], newdata) else newdata
      fits      <- lapply(object, stats::predict, newdata=newdata, type="response")
    }
      if(isTRUE(simplify)) provideDimnames(array(unlist(fits), dim=c(nr, d, G)), base=if(d == 1) list(dnames, if(exp != "~1") all.vars(stats::as.formula(exp)) else "", gnames) else c(dnames, list(gnames))) else fits
  }
  
#' @rdname predict.MoE_expert
#' @method fitted MoE_expert
#' @keywords prediction utility
#' @usage 
#' \method{fitted}{MoE_expert}(object,
#'        ...)
#' @export
  fitted.MoE_expert      <- function(object, ...) {
    args   <- c(list(object=object, newdata=NULL), as.list(match.call())[-1L])
    fits   <- do.call(predict.MoE_expert, args[unique(names(args))])
      if(!is.null(fits))    return(fits)
  }
  
#' @rdname predict.MoE_expert
#' @method residuals MoE_expert
#' @keywords prediction utility
#' @usage 
#' \method{residuals}{MoE_expert}(object,
#'           ...)
#' @export
  residuals.MoE_expert   <- function(object, ...)  {
    args   <- c(list(object=object, simplify=FALSE), as.list(match.call())[-1L])
    simple <- any(names(list(...)) == "simplify") && isTRUE(args$simplify)
    fits   <- do.call(fitted.MoE_expert, args[unique(names(args))])
    dat    <- attr(object, "Data")
    dat    <- if(ncol(dat) == 1) dat[[1L]] else dat
    fits   <- lapply(fits, function(x) dat - x)
      if(isTRUE(simple)) provideDimnames(array(unlist(fits), dim=c(dim(fits[[1L]]), length(fits))), 
                                         base=c(dimnames(fits[[1L]]), list(names(fits)))) else fits
  }

#' Predictions from MoEClust gating networks
#' 
#' Predicts mixing proportions from MoEClust gating networks. Effectively akin to predicting from a multinomial logistic regression via \code{\link[nnet]{multinom}}, although here the noise component (if any) is properly accounted for. So too are models with no gating covariates at all, or models with the equal mixing proportion constraint. Prior probabilities are returned by default.
#' @param object An object of class \code{"MoE_gating"} (typically \code{x$gating}, where \code{x} is of class \code{"MoEClust"}).
#' @param newdata A matrix or data frame of test examples. If omitted, the fitted values are used.
#' @param type The type of output desired. The default (\code{"probs"}) returns prior probabilities, while \code{"class"} returns labels indicating the most likely group \emph{a priori}. Note that observations classified assigned the noise component (if any) are given a label of \code{0}.
#' @param keep.noise A logical indicating whether the output should acknowledge the noise component (if any). Defaults to \code{TRUE}; when \code{FALSE}, this column is discarded and the matrix of probabilities is renormalised accordingly.
#' @param droplevels A logical indicating whether unseen factor levels in categorical variables within \code{newdata} should be dropped (with \code{NA} predicted in their place). Defaults to \code{FALSE}. See \code{\link{drop_levels}}.
#' @param ... Catches unused arguments or allows the \code{type} and \code{keep.noise} arguments to be passed through \code{fitted} and the \code{keep.noise} argument to be passed through \code{residuals}.
#' 
#' @return The return value depends on whether \code{newdata} is supplied or not and whether the model includes gating covariates to begin with. When \code{newdata} is not supplied, the fitted values are returned (as a matrix if the model contained gating covariates, otherwise as a vector as per \code{x$parameters$pro}). If \code{newdata} is supplied, the output is always a matrix with the same number of rows as the \code{newdata}.
#' @details This function is effectively a shortcut to \code{\link{predict.MoEClust}(x, ...)$pro}, which (unlike the \code{predict} method for \code{\link[nnet]{multinom}} on which \code{predict.MoE_gating} is based) accounts for the various ways of treating gating covariates, equal mixing proportion constraints, and noise components, although its \code{type} argument defaults to \code{"probs"} rather than \code{"class"}.
#' @note Note that the \code{keep.noise} argument does \strong{not} correspond in any way to the \code{discard.noise} argument to \code{\link{predict.MoEClust}}; there, the noise component is respected in the computation of the mixing proportions and only discarded (if at all) in the prediction of the responses.
#' 
#' @references Murphy, K. and Murphy, T. B. (2020). Gaussian parsimonious clustering models with covariates and a noise component. \emph{Advances in Data Analysis and Classification}, 14(2): 293-325. <\doi{10.1007/s11634-019-00373-8}>.
#' @author Keefe Murphy - <\email{keefe.murphy@@mu.ie}>
#' @seealso \code{\link{predict.MoEClust}}, \code{\link[nnet]{multinom}}, \code{\link{predict.MoE_expert}}, \code{\link{drop_levels}}
#' @method predict MoE_gating
#' @keywords prediction utility
#' @importFrom nnet "multinom"
#' @export
#' @usage 
#' \method{predict}{MoE_gating}(object,
#'         newdata = NULL,
#'         type = c("probs", "class"),
#'         keep.noise = TRUE,
#'         droplevels = FALSE,
#'         ...)
#' @examples 
#' data(ais)
#' mod    <- MoE_clust(ais[,3:7], G=2, modelNames="EEE", gating= ~ SSF + Ht,
#'                  expert= ~ sex, network.data=ais, tau0=0.1, noise.gate=FALSE)
#' (preds <- predict(mod$gating, newdata=ais[1:5,]))
#' 
#' all.equal(preds, predict(mod, newdata=ais[1:5,])$pro) #TRUE
#' 
#' # Note that the predictions are not the same as the multinom predict method
#' # in this instance, owing to the invocation of noise.gate=FALSE above
#' mod2   <- mod 
#' class(mod2$gating) <- c("multinom", "nnet")
#' predict(mod2$gating, newdata=ais[1:5,], type="probs")
#' 
#' # We can make this function behave in the same way by invoking keep.noise=FALSE
#' predict(mod$gating, keep.noise=FALSE, newdata=ais[1:5,])
#' 
#' # ... although keep.noise=FALSE in predict.MoE_gating does not
#' # yield the same output as discard.noise=TRUE in predict.MoEClust
#' predict(mod, discard.noise=TRUE, newdata=ais[1:5,])$pro
  predict.MoE_gating     <- function(object, newdata = NULL, type = c("probs", "class"), keep.noise = TRUE, droplevels = FALSE, ...) {
    if(!is.null(newdata) &&
       all(!is.matrix(newdata), 
           !is.data.frame(newdata)))               stop("'newdata' must be a matrix or data frame if supplied", call.=FALSE)
    if(!missing(type)    && 
       length(type)  > 1 || !is.character(type))   stop("'type' must be a single character string",             call.=FALSE)
    if(length(keep.noise) > 1   ||
       !is.logical(keep.noise))                    stop("'keep.noise' must be a single logical indicator",      call.=FALSE)
    if(length(droplevels) > 1   ||
       !is.logical(droplevels))                    stop("'droplevels' must be a single logical indicator",      call.=FALSE)
    class(object)        <- class(object)[-1L]
    fits   <- object$fitted.values
    noise  <- attr(object, "Noise")
    G      <- attr(object, "G")
    GN     <- G + noise
    gnames <- paste0("Cluster",      if(noise) replace(seq_len(GN), GN, "0")             else seq_len(GN))
    gat    <- attr(object, "Formula")
    if(gat == "~1")       {
      fits <- matrix(attr(object, "NoGate"), nrow=ifelse(is.null(newdata), 1L, NROW(newdata)), ncol=GN, byrow=TRUE)
    } else if(!is.null(newdata)) { 
      fits <- stats::predict(object, if(isTRUE(droplevels)) drop_levels(object, newdata) else newdata, type="probs")
      fits <- if(is.matrix(fits)) fits else matrix(fits, nrow=1L, ncol=length(fits), byrow=TRUE)
    }
    if(all(noise, !attr(object, "NoiseGate"))) {
      fits <- .tau_noise(fits, attr(object, "NoisePro"))
    }
    colnames(fits)       <- NULL
    if(isFALSE(keep.noise)      && noise)      {
      if(G == 0)                                   stop("Nothing to return as the model has only a noise component: use keep.noise=TRUE", call.=FALSE)
      fits <- .renorm_z(fits[,-GN, drop=FALSE])
    }
      switch(EXPR=match.arg(type), 
             probs=provideDimnames(fits, base=list(ifelse(nrow(fits) == 1, "pro", ""), gnames)), {
        if(attr(object, "EqualPro") && GN > 1) {
          if(!all(noise, keep.noise, !attr(object, "EqualNoise"),
                  max(fits)     == fits[GN]))      message("class predicted at random due to the equal mixing proportion constraint\n")  
        }
        CL <- max.col(fits)
          if(all(noise, isTRUE(keep.noise))) replace(CL, CL == GN, 0L) else CL
      })
  }

#' @rdname predict.MoE_gating
#' @method fitted MoE_gating
#' @keywords prediction utility
#' @importFrom nnet "multinom"
#' @usage 
#' \method{fitted}{MoE_gating}(object,
#'        ...)
#' @export
  fitted.MoE_gating      <- function(object, ...) {
    args   <- c(list(object=object, newdata=NULL), as.list(match.call())[-1L])
    fits   <- do.call(predict.MoE_gating, args[unique(names(args))])
      if(!is.null(fits))    return(fits)
  }
  
#' @rdname predict.MoE_gating
#' @method residuals MoE_gating
#' @keywords prediction utility
#' @importFrom nnet "multinom"
#' @usage 
#' \method{residuals}{MoE_gating}(object,
#'           ...)
#' @export
  residuals.MoE_gating   <- function(object, ...) {
    dat.z  <- attr(object, "Data")
    args   <- c(list(object=object, type="probs", newdata=dat.z), as.list(match.call())[-1L])
    keep   <- !any(names(list(...)) == "keep.noise") || isTRUE(args$keep.noise)
    fits   <- do.call(fitted.MoE_gating, args[unique(names(args))])
    dat.z  <- if(isTRUE(keep) || !attr(object, "Noise")) dat.z else .renorm_z(dat.z[,-ncol(dat.z), drop=FALSE])
    dat.z[is.nan(dat.z)]      <- 0L
    fits   <- if(nrow(fits)   == nrow(dat.z))             fits else matrix(fits, nrow=nrow(dat.z), ncol=ncol(dat.z), byrow=TRUE)
      tryCatch(dat.z - fits, error=function(e) 1L - fits)
  }
  
#' Stepwise model/variable selection for MoEClust models
#'
#' Conducts a greedy forward stepwise search to identify the optimal \code{MoEClust} model according to some \code{criterion}. Components and/or \code{gating} covariates and/or \code{expert} covariates are added to new \code{\link{MoE_clust}} fits at each step, while each step is evaluated for all valid \code{modelNames}.
#' @param data A numeric vector, matrix, or data frame of observations. Categorical variables are not allowed. If a matrix or data frame, rows correspond to observations and columns correspond to variables.
#' @param network.data An optional matrix or data frame in which to look for the covariates specified in the \code{gating} &/or \code{expert} networks, if any. Must include column names. Columns in \code{network.data} corresponding to columns in \code{data} will be automatically removed. While a single covariate can be supplied as a vector (provided the '\code{$}' operator or '\code{[]}' subset operator are not used), it is safer to supply a named 1-column matrix or data frame in this instance.
#' @param gating A vector giving the names of columns in \code{network.data} used to define the scope of the gating network. By default, the initial model will contain no covariates (unless \code{initialModel} is supplied with gating covariates), thereafter all variables in \code{gating} (save for those in \code{initialModel}, if any) will be considered for inclusion where appropriate.
#' 
#' If \code{gating} is not supplied (or set to \code{NULL}), \emph{all} variables in \code{network.data} will be considered for the gating network. \code{gating} can also be supplied as \code{NA}, in which case \emph{no} gating network covariates will ever be considered (save for those in \code{initialModel}, if any). Supplying \code{gating} and \code{expert} can be used to ensure different subsets of covariates enter different parts of the model.
#' @param expert A vector giving the names of columns in \code{network.data} used to define the scope of the expert network. By default, the initial model will contain no covariates (unless \code{initialModel} is supplied with expert covariates), thereafter all variables in \code{expert} (save for those in \code{initialModel}, if any) will be considered for inclusion where appropriate.
#' 
#' If \code{expert} is not supplied (or set to \code{NULL}), \emph{all} variables in \code{network.data} will be considered for the expert network. \code{expert} can also be supplied as \code{NA}, in which case \emph{no} expert network covariates will ever be considered (save for those in \code{initialModel}, if any). Supplying \code{expert} and \code{gating} can be used to ensure different subsets of covariates enter different parts of the model.
#' @param modelNames A character string of valid model names, to be used to restrict the size of the search space, if desired. By default, \emph{all} valid model types are explored. Rather than considering the changing of the model type as an additional step, every step is evaluated over all entries in \code{modelNames}. See \code{\link{MoE_clust}} for more details. 
#' 
#' Note that if \code{initialModel} is supplied (see below), \code{modelNames} will be augmented with \code{initialModel$modelName} if needs be. 
#' @param fullMoE A logical which, when \code{TRUE}, ensures that only models where the same covariates enter both parts of the model (the gating and expert networks) are considered. This restricts the search space to exclude models where covariates differ across networks. Thus, the search is likely to be faster, at the expense of potentially missing out on optimal models. Defaults to \code{FALSE}. 
#' 
#' Furthermore, when \code{TRUE}, the set of candidate covariates is automatically taken to be the \strong{union} of the \emph{named} covariates in \code{gating} and \code{expert}, for convenience. In other words, \code{gating=NA} will only work if \code{expert=NA} also, and both should be set to \code{NULL} in order to consider all potential covariates. 
#' 
#' In addition, caution is advised using this argument in conjunction with \code{initialModel}, which must satisfy the constraint that the same set of covariates be used in both parts of the model, for initial models where gating covariates are allowable. Finally, note that this argument does not preclude a model with only expert covariates included if the number of components is such that the inclusion of gating covariates is infeasible.
#' @param noise A logical indicating whether to assume all models contain an additional noise component (\code{TRUE}) or not (\code{FALSE}, the default). If \code{initialModel} or \code{initialG} is not specified, the search starts from a \code{G=0} noise-only model when \code{noise} is \code{TRUE}, otherwise the search starts from a \code{G=1} model with no covariates when \code{noise} is \code{FALSE}. See \code{\link{MoE_control}} for more details. Note, however, that if the model specified in \code{initialModel} contains a noise component, the value of the \code{noise} argument will be overridden to \code{TRUE}; similarly, if the \code{initialModel} model does not contain a noise component, \code{noise} will be overridden to \code{FALSE}.
#' @param initialModel An object of class \code{"MoEClust"} generated by \code{\link{MoE_clust}} or an object of class \code{"MoECompare"} generated by \code{\link{MoE_compare}}. This gives the initial model to use at the first step of the selection algorithm, to which components and/or covariates etc. can be added. Especially useful if the model is expected to have more than one component \emph{a priori} (see \code{initialG} below as an alternative). The \code{initialModel} model must have been fitted to the same data in \code{data}. 
#' 
#' If \code{initialModel} is not specified, the search starts from a \code{G=0} noise-only model when \code{noise} is \code{TRUE}, otherwise the search starts from a \code{G=1} model with no covariates when \code{noise} is \code{FALSE}. If \code{initialModel} \emph{is} supplied and it contains a noise component, only models with a noise component will be considered thereafter (i.e. the \code{noise} argument can be overridden by the \code{initialModel} argument). If \code{initialModel} contains gating &/or expert covariates, these covariates will be included in all subsequent searches, with covariates in \code{expert} and \code{gating} still considered as candidates for additional inclusion, as normal.
#' 
#' However, while \code{initialModel} \emph{can} include covariates not specified in \code{gating} &/or \code{expert}, the \code{initialModel$modelName} \strong{should} be included in the specified \code{modelNames}; if it is not, \code{modelNames} will be forcibly augmented with \code{initialModel$modelName} (as stated above). Furthermore, it is assumed that \code{initialModel} is already optimal with respect to the model type. If it is not, the algorithm may be liable to converge to a sub-optimal model, and so a warning will be printed if the function suspects that this \emph{might} be the case.
#' @param initialG A single (positive) integer giving the number of mixture components (clusters) to initialise the stepwise search algorithm with. This is a simpler alternative to the \code{initialModel} argument, to be used when the only prior knowledge relates to the number of components, and not other features of the model (e.g. the covariates which should be included). Consequently, \code{initialG} is only relevant when \code{initialModel} is not supplied. When neither \code{initialG} nor \code{initialModel} is specified, the search starts from a \code{G=0} noise-only model when \code{noise} is \code{TRUE}, otherwise the search starts from a \code{G=1} model with no covariates when \code{noise} is \code{FALSE}. See \code{stepG} below for fixing the number of components at this \code{initialG} value.
#' @param stepG A logical indicating whether the algorithm should consider incrementing the number of components at each step. Defaults to \code{TRUE}; use \code{FALSE} when searching only over configurations with the same number of components is of interest. Setting \code{stepG} to \code{FALSE} is possible with or without specifying \code{initialModel} or \code{initialG}, but is primarily intended for use when one of these arguments is supplied, otherwise the algorithm will be stuck forever with only one component.
#' @param criterion The model selection criterion used to determine the optimal action at each step. Defaults to \code{"bic"}.
#' @param equalPro A character string indicating whether models with equal mixing proportions should be considered. \code{"both"} means models with both equal and unequal mixing proportions will be considered, \code{"yes"} means only models with equal mixing proportions will be considered, and \code{"no"} means only models with unequal mixing proportions will be considered. Notably, no setting for \code{equalPro} is enough to rule out models with \code{gating} covariates from consideration.
#' 
#' The default (\code{"all"}) is equivalent to \code{"both"} with the addition that all possible mixing proportion constraints will be tried for the \code{initialModel} (if any, provided it doesn't contain gating covariate(s)) or \code{initialG} \emph{before} adding a component or additional covariates; otherwise, this \code{equalPro} argument only governs whether mixing proportion constraints are considered as components are added.
#' 
#' Considering \code{"all"} (or \code{"both"}) equal and unequal mixing proportion models increases the search space and the computational burden, but this argument becomes irrelevant after a model, if any, with gating network covariate(s) is considered optimal for a given step. The \code{"all"} default is \strong{strongly} recommended so that viable candidate models are not missed out on, particularly when \code{initialModel} or \code{initialG} are given. However, this does not guarantee that an optimal model will not be skipped; if \code{equalPro} is restricted via \code{"yes"} or \code{"no"}, a suboptimal model at one step may ultimately lead to a better final model, in some edge cases. See \code{\link{MoE_control}} for more details.
#' @param noise.gate A character string indicating whether models where the gating network for the noise component depends on covariates are considered. \code{"yes"} means only models where this is the case will be considered, \code{"no"} means only models for which the noise component's mixing proportion is constant will be considered and \code{"both"} means both of these scenarios will be considered.
#' 
#' The default (\code{"all"}) is equivalent to \code{"both"} with the addition that all possible gating network noise settings will be tried for the \code{initialModel} (if any, provided it contains gating covariates and a noise component) \emph{before} adding a component or additional covariates; otherwise, this \code{noise.gate} argument only governs the inclusion/exclusion of this constraint as components or covariates are added.
#' 
#' Considering \code{"all"} (or \code{"both"}) settings increases the search space and the computational burden, but this argument is only relevant when \code{noise=TRUE} and \code{gating} covariates are being considered. The \code{"all"} default is \strong{strongly} recommended so that viable candidate models are not missed out on, particularly when \code{initialModel} or \code{initialG} are given. However, this does not guarantee that an optimal model will not be skipped; if \code{noise.gate} is restricted via \code{"yes"} or \code{"no"}, a suboptimal model at one step may ultimately lead to a better final model, in some edge cases. See \code{\link{MoE_control}} for more details.
#' @param verbose Logical indicating whether to print messages pertaining to progress to the screen during fitting. By default is \code{TRUE} if the session is interactive, and \code{FALSE} otherwise. If \code{FALSE}, warnings and error messages will still be printed to the screen, but everything else will be suppressed.
#' @param ... Additional arguments to \code{\link{MoE_control}}. Note that these arguments will be supplied to \emph{all} candidate models for every step.
#'
#' @return An object of class \code{"MoECompare"} containing information on all visited models and the optimal model (accessible via \code{x$optimal}).
#' @details The arguments \code{modelNames}, \code{equalPro}, and \code{noise.gate} are provided for computational convenience. They can be used to reduce the number of models under consideration at each stage. 
#' 
#' The same is true of the arguments \code{gating} and \code{expert}, which can each separately (or jointly, if \code{fullMoE} is \code{TRUE}) be made to consider all variables in \code{network.data}, or a subset, or none at all. 
#' 
#' Finally, \code{initialModel} or \code{initialG} can be used to kick-start the search algorithm by incorporating prior information in a more direct way; in the latter case, only in the form of the number of components; in the former case, a full model with a given number of components, certain included gating and expert network covariates, and a certain model type can give the model an even more informed head start. In either case, the \code{stepG} argument can be used to fix the number of components and only search over different configurations of covariates.
#' 
#' Without any prior information, it is best to accept the defaults at the expense of a longer run-time.
#' @note It is advised to run this function once with \code{noise=FALSE} and once with \code{noise=TRUE} and then choose the optimal model across both sets of results.
#' 
#' At present, only additions (of components and covariates) are considered. In future updates, it will be possible to allow both additions and removals.
#' 
#' The function will attempt to remove duplicate variables found in both \code{data} and \code{network.data}; in particular, they will be removed from \code{network.data}. Users are however advised to careful specify \code{data} and \code{network.data} such that there are no duplicates, especially if the desired variable(s) should belong to \code{network.data}.
#' 
#' Finally, if the user intends to search for the best model according to the \code{"icl"} \code{criterion}, then specifying either \code{initialModel} or \code{initialG} is advisable. This is because the algorithm otherwise starts with a single component and thus there is no entropy term, meaning the stepwise search can quickly and easily get stuck at \code{G=1}. See the examples below.
#' @export
#' @author Keefe Murphy - <\email{keefe.murphy@@mu.ie}>
#' @references Murphy, K. and Murphy, T. B. (2020). Gaussian parsimonious clustering models with covariates and a noise component. \emph{Advances in Data Analysis and Classification}, 14(2): 293-325. <\doi{10.1007/s11634-019-00373-8}>.
#' @seealso \code{\link{MoE_clust}}, \code{\link{MoE_compare}}, \code{\link{MoE_control}}
#' @keywords clustering main
#' @usage
#' MoE_stepwise(data,
#'              network.data = NULL,
#'              gating = NULL,
#'              expert = NULL,
#'              modelNames = NULL,
#'              fullMoE = FALSE,
#'              noise = FALSE,
#'              initialModel = NULL,
#'              initialG = NULL,
#'              stepG = TRUE,
#'              criterion = c("bic", "icl", "aic"),
#'              equalPro = c("all", "both", "yes", "no"),
#'              noise.gate = c("all", "both", "yes", "no"),
#'              verbose = interactive(),
#'              ...)
#' @examples
#' \donttest{# data(CO2data)
#' # Search over all models where the single covariate can enter either network
#' # (mod1  <- MoE_stepwise(CO2data$CO2, CO2data[,"GNP", drop=FALSE]))
#' #
#' # data(ais)
#' # Only look for EVE & EEE models with at most one expert network covariate
#' # Do not consider any gating covariates and only consider models with equal mixing proportions
#' # (mod2  <- MoE_stepwise(ais[,3:7], ais, gating=NA, expert="sex",
#' #                        equalPro="yes", modelNames=c("EVE", "EEE")))
#' #
#' # Look for models with noise & only those where the noise component's mixing proportion is constant
#' # Speed up the search with an initialModel, fix G, and restrict the covariates & model type
#' # init   <- MoE_clust(ais[,3:7], G=2, modelNames="EEE", 
#' #                     expert= ~ sex, network.data=ais, tau0=0.1)
#' # (mod3  <- MoE_stepwise(ais[,3:7], ais, noise=TRUE, expert="sex",
#' #                        gating=c("SSF", "Ht"), noise.gate="no", 
#' #                        initialModel=init, stepG=FALSE, modelNames="EEE"))
#' #
#' # Compare both sets of results (with & without a noise component) for the ais data
#' # (comp1 <- MoE_compare(mod2, mod3, optimal.only=TRUE))
#' # comp1$optimal
#' #
#' # Target a model for the AIS data which is optimal in terms of ICL, without any restrictions
#' # mod4   <- MoE_stepwise(ais[,3:7], ais, criterion="icl")
#' # 
#' # This gets stuck at a G=1 model, so specify an initial G value as a head start
#' # mod5   <- MoE_stepwise(ais[,3:7], ais, criterion="icl", initialG=2)
#' #
#' # Check that specifying an initial G value enables a better model to be found
#' # (comp2 <- MoE_compare(mod4, mod5, optimal.only=TRUE, criterion="icl"))
#' 
#' # Finally, restrict the search to full MoE models only
#' # Notice that the candidate covariates are the union of gating and expert
#' # Notice also that the algorithm initially traverses models with only
#' #   expert covariates when the inclusion of gating covariates is infeasible
#' # mod6   <- MoE_stepwise(ais[,3:7], ais, fullMoE=TRUE, gating="BMI", expert="Bfat")}
  MoE_stepwise    <- function(data, network.data = NULL, gating = NULL, expert = NULL, modelNames = NULL, 
                              fullMoE = FALSE, noise = FALSE, initialModel = NULL, initialG = NULL, stepG = TRUE, 
                              criterion = c("bic", "icl", "aic"), equalPro = c("all", "both", "yes", "no"), 
                              noise.gate = c("all", "both", "yes", "no"), verbose = interactive(), ...) {
    
    call          <- match.call(expand.dots=TRUE)
    compare       <- list()
    gate.x        <- is.null(gating)
    exps.x        <- is.null(expert)
    if(!is.null(network.data))         {
      if((!is.matrix(network.data)    &&
        !is.data.frame(network.data)) &&
        NCOL(network.data) == 1)       {  
        tmpN      <- deparse(substitute(network.data))
        if(any(grepl("\"",  tmpN),
               grepl("\\[", tmpN)))    {          stop("Invalid 'network.data': must be a matrix or data frame with named columns",      call.=FALSE)
        }     else {
          network.data <- provideDimnames(as.matrix(network.data), base=list("", tmpN))  
        }
      }
      if(is.null(colnames(network.data)))         stop("Invalid 'network.data': must be a matrix or data frame with named columns",      call.=FALSE)
      if(NROW(data)    != nrow(network.data))     stop("Invalid 'network.data': must contain the same number of observations as 'data'", call.=FALSE)
    }
    
    tmp.nam       <- as.character(substitute(data))
    data          <- as.data.frame(data)
    if(any(grepl("\\$", colnames(network.data)))) stop("'network.data' column names cannot contain the '$' operator", call.=FALSE)
    if(!is.null(gating)               &&
      ((length(gating) != 1           || 
       !all(is.na(gating)))           && 
       !is.character(gating)))                    stop("Invalid 'gating': must be NA or a vector of variable names",  call.=FALSE)
    if(!is.null(expert)               &&
      ((length(expert) != 1           ||
       !all(is.na(expert)))           && 
       !is.character(expert)))                    stop("Invalid 'expert': must be NA or a vector of variable names",  call.=FALSE)
    if(length(fullMoE)  > 1           ||
       !is.logical(fullMoE))                      stop("'fullMoE' must be a single logical indicator",     call.=FALSE)
    if(length(noise)    > 1           ||
       !is.logical(noise))                        stop("'noise' must be a single logical indicator",       call.=FALSE)
    both.x        <- all(gate.x, exps.x)
    b.nets        <- if(!both.x       && 
                        all(!gate.x   && 
                            is.na(gating), 
                            !exps.x   &&
                            is.na(expert)))    NA else unique(c(gating[!is.na(gating)], 
                                                                expert[!is.na(expert)]))
    b.nets        <- if(length(b.nets) > 0) b.nets
    if(isTRUE(fullMoE)) {
      gating      <- 
      expert      <- b.nets
    }
    if((has.init  <- 
        !is.null(initialModel))  &&
        !is.null(initialG))                       stop("Only one of 'initialModel' and 'initialG' can be supplied",            call.=FALSE)
    if(has.init   &&
      (!inherits(initialModel, "MoEClust")  &&
       !inherits(initialModel, "MoECompare")))    stop("'initialModel' must be an object of class 'MoEClust' or 'MoECompare'", call.=FALSE)
    initialModel  <- if(inherits(initialModel, "MoECompare")) initialModel$optimal else initialModel
    if(!is.null(initialG)        &&
      (length(initialG)          != 1 ||
       !is.numeric(initialG)     ||
       initialG   <= 0           ||
       floor(initialG) != initialG))              stop("'initialG' must be a single positive integer",     call.=FALSE)
    if(length(stepG)    > 1           ||
       !is.logical(stepG))                        stop("'stepG' must be a single logical indicator",       call.=FALSE)
    if(!missing(criterion)            && 
       length(criterion)          > 1 ||
       !is.character(criterion))                  stop("'criterion' must be a single character string",    call.=FALSE)
    if(!missing(equalPro)             && 
       length(equalPro)           > 1 ||
       !is.character(equalPro))                   stop("'equalPro' must be a single character string",     call.=FALSE)
    if(!missing(noise.gate)           && 
       length(noise.gate)         > 1 ||
       !is.character(noise.gate))                 stop("'noise.gate' must be a single character string",   call.=FALSE)
    if(length(verbose)  > 1           ||
       !is.logical(verbose))                      stop("'verbose' must be a single logical indicator",     call.=FALSE)
    num.X         <- vapply(data, is.numeric, logical(1L))
    if(anyNA(data))     {
      if(isTRUE(verbose))                         message("Rows with missing values removed from data\n")
      data        <- data[stats::complete.cases(data),, drop=FALSE]
    }
    if(sum(num.X) != ncol(data))       {
      if(isTRUE(verbose))                         message("Non-numeric columns removed from data\n")
      data        <- data[,num.X,                       drop=FALSE]
    }
    if(ncol(data) == 1) {
      colnames(data)   <- tmp.nam[length(tmp.nam)]
    }
    if(isTRUE(has.init))          {
     if(!identical(initialModel$data, data))      stop("The 'initialModel' model must have been fitted to the same data!",                       call.=FALSE)
      if(!is.null(initialModel$call$modelName))   warning("Are you sure that initialModel$modelName is optimal for the given model?\n",          call.=FALSE, immediate.=TRUE)
      if(!is.null(modelNames)    &&
         !(initialModel$modelName         %in% 
           modelNames))           {               warning("Set of candidate 'modelNames' expanded with initialModel$modelName\n",                call.=FALSE, immediate.=TRUE)
       modelNames <- unique(c(initialModel$modelName, modelNames))
      }
      init.exp    <- stats::as.formula(attr(initialModel$expert, "Formula"))
      init.gate   <- stats::as.formula(attr(initialModel$gating, "Formula"))
      if(fullMoE  && 
         initialModel$G > (1L + !attr(initialModel, "NoiseGate")) &&
         init.exp != init.gate)                   stop("The 'initialModel' must have the same covariates in both networks if 'fullMoE' is TRUE", call.=FALSE)
      o.exps      <- all.vars(init.exp)
      o.gate      <- all.vars(init.gate)
      expert      <- if(is.null(expert)     || all(is.na(expert))) expert else unique(c(o.exps, expert)) 
      gating      <- if(is.null(gating)     || all(is.na(gating))) gating else unique(c(o.gate, gating))
      o.nets      <- unique(c(o.exps, o.gate))
      b.nets      <- unique(c(gating, expert))
      b.nets      <- if(!is.null(b.nets)    && 
                        anyNA(b.nets))             b.nets[!is.na(b.nets)] else b.nets
      network.data               <- cbind(initialModel$net.covs[,o.nets, drop=FALSE], 
                                          network.data[,setdiff(colnames(network.data), o.nets), drop=FALSE])
    }
    if(!is.null(network.data))    {
     dup.ind      <- if(any(is.matrix(data), is.data.frame(data)) && !any(grepl("\\$", colnames(data)))) (colnames(network.data) %in% colnames(data)) else vapply(seq_len(ncol(network.data)), function(j) isTRUE(all.equal(network.data[,j], unname(unlist(data)))), logical(1L))
     if(any(dup.ind))                             warning("Removing covariates found in response data\n", call.=FALSE, immediate.=TRUE)  
     network.data <- network.data[,!dup.ind, drop=FALSE]  
    }
    network.data  <- as.data.frame(network.data)
    network.char  <- vapply(network.data, is.character, logical(1L))
    network.data[network.char]   <- lapply(network.data[network.char], factor)
    if(any(network.char))                         message("Character covariates coerced to factors\n")
    netnames      <- colnames(network.data)
    if(suppressWarnings(!all(is.na(
       as.numeric(netnames)))))                   stop("Invalid variable names in 'network.data'",    call.=FALSE)
    if(has.init   && !is.null(initialG))          warning("'initialG' overruled by 'initialModel'\n", call.=FALSE, immediate.=TRUE)
    G             <- ifelse(has.init, initialModel$G, ifelse(is.null(initialG), 0L + isFALSE(noise), initialG))
    noise         <- ifelse(has.init, attr(initialModel, "Noise")    || G == 0, noise)
    na.gate       <- !gate.x     && all(is.na(gating))
    na.expx       <- !exps.x     && all(is.na(expert))
    na.both       <- !both.x     && !is.null(b.nets) && all(is.na(b.nets))
    if(!gate.x    && !na.gate         &&
      (!is.character(gating)          ||
       !all(gating   %in% netnames)))             stop("Invalid 'gating': must be a vector of variable names in 'network.data'", call.=FALSE)
    if(!exps.x    && !na.expx         &&
      (!is.character(expert)          ||
       !all(expert   %in% netnames)))             stop("Invalid 'expert': must be a vector of variable names in 'network.data'", call.=FALSE)
    netnames      <- if(gate.x   || exps.x) colnames(network.data) else unique(c(gating[!na.gate], expert[!na.expx]))
    network.data  <- network.data[,netnames, drop=FALSE]
    criterion     <- match.arg(criterion)
    equalPro      <- match.arg(equalPro)
    allPro        <- equalPro    == "all"
    unequalPro    <- is.element(equalPro,   c("all", "both", "no"))
    equalPro      <- is.element(equalPro,   c("all", "both", "yes"))
    noise.gate    <- match.arg(noise.gate)
    allNgate      <- noise.gate  == "all"
    noNgate       <- is.element(noise.gate, c("all", "both", "no"))  && isTRUE(noise)
    doNgate       <- is.element(noise.gate, c("all", "both", "yes")) || isFALSE(noise)
    gcov          <- if((gate.x  && 
                       !fullMoE) || is.null(gating)) paste0("~", netnames) else if(!na.gate) paste0("~", gating)
    ecov          <- if((exps.x  && 
                       !fullMoE) || is.null(expert)) paste0("~", netnames) else if(!na.expx) paste0("~", expert)
    bcov          <- if(both.x   || is.null(b.nets)) paste0("~", netnames) else if(!na.both) paste0("~", b.nets)
    na.gate       <- ifelse(fullMoE, !is.null(gating) && all(is.na(gating)), na.gate)
    na.expx       <- ifelse(fullMoE, !is.null(expert) && all(is.na(expert)), na.expx)
    if(isTRUE(has.init))          {
      if(length(o.exps)  == 0)    {
        init.exp  <- ~1
        o.exps    <- NULL
      }       else {
        ecov      <- setdiff(ecov, paste0("~", o.exps))
        any.e     <- is.element(ecov, paste0("~", as.character(init.exp)[-1L]))
        ecov      <- paste0("~", paste(paste(o.exps, collapse=" + "), setdiff(all.vars(init.exp),  o.exps), sep=" + "), 
                            gsub(paste(o.exps, collapse=" \\+ "), '', gsub('\\~', '', ecov[!any.e])))
        ecov      <- gsub(" $", "", gsub("\\+$", "", gsub(" $", "", ecov)))
      }
      if(length(o.gate)  == 0)    {
        init.gate <- ~1
        o.gate    <- NULL
      }       else {
        gcov      <- setdiff(gcov, paste0("~", o.gate))
        any.g     <- is.element(gcov, paste0("~", as.character(init.gate)[-1L]))
        gcov      <- paste0("~", paste(paste(o.gate, collapse=" + "), setdiff(all.vars(init.gate), o.gate), sep=" + "), 
                            gsub(paste(o.gate, collapse=" \\+ "), '', gsub('\\~', '', gcov[!any.g])))
        gcov      <- gsub(" $", "", gsub("\\+$", "", gsub(" $", "", gcov)))
      }
      if(isTRUE(fullMoE))         {
       if(length(o.nets) == 0)    {
        init.exp  <-
        init.gate <- ~1
        o.exps    <- 
        o.gate    <- 
        o.nets    <- NULL 
       }      else {
         bcov     <- setdiff(bcov, paste0("~", o.nets))
         any.b    <- is.element(bcov, paste0("~", as.character(init.exp)[-1L]))
         bcov     <- paste0("~", paste(paste(o.nets, collapse=" + "), setdiff(all.vars(init.exp), o.nets), sep=" + "), 
                            gsub(paste(o.nets, collapse=" \\+ "), '', gsub('\\~', '', bcov[!any.b])))
         bcov     <- gsub(" $", "", gsub("\\+$", "", gsub(" $", "", bcov)))
       }
      }
    }         else {
      init.gate   <- 
      init.exp    <- ~1
      o.exps      <-
      o.gate      <- 
      o.nets      <- NULL
    }
    dots          <- list(...)
    if(any(names(dots) == "control")  &&
       length(dots[names(dots)       %in%
              names(MoE_control())])   > 0)       stop("Arguments cannot be supplied via the '...' construct when the named argument 'control' is supplied", call.=FALSE) 
    dots          <- dots[names(dots) != "z.list"]
    args          <- c(list(data=data, G=G, gating=init.gate, expert=init.exp, network.data=network.data, modelNames=modelNames, equalPro=equalPro == "yes", 
                       noise.gate=noise.gate != "no", verbose=FALSE, tau0=if(any(names(dots) == "tau0")) dots$tau0 else 0.1), criterion=criterion)
    args          <- if(length(dots)  >= 1) c(args, dots) else args
    args          <- args[unique(names(args))]
    args$tau0     <- if(isTRUE(noise)) args$tau0
    if(isTRUE(verbose))                           message("\n\tStep ", 1, "...\n")
    res           <- if(isTRUE(has.init))    initialModel else suppressWarnings(do.call(MoE_clust, args))
    args$z.list   <- INIT        <- list(if(isTRUE(has.init) && attr(res, "Expert")) attr(res, "Exp.init")         else attr(res, "Z.init"))
    args$z.init   <- dots$init.z
    args$init.z   <- "list"
    compare[[j    <- 1L]]        <- res
    if(isTRUE(verbose))            suppressWarnings(print(suppressWarnings(MoE_compare(res, optimal.only=TRUE, pick=1L, criterion=criterion)), details=FALSE))
    crit          <- switch(EXPR=criterion, bic=res$bic, icl=res$icl, aic=res$aic)
    crit.old      <- -Inf
    if(isTRUE(has.init)          ||
       !is.null(initialG))        {
      addC        <- stepG       && (isTRUE(unequalPro) || (has.init   && attr(initialModel, "Gating"))  ||  G == 0)
      addE        <- G  > 0      && length(ecov) > 0    && !na.expx    && (is.null(o.exps) || !identical(expert,                     if(fullMoE) o.nets else o.exps))
      addG        <- !na.gate    && length(gcov) > 0    && G > !noise  && isTRUE(doNgate)  && (is.null(o.gate) || !identical(gating, if(fullMoE) o.nets else o.gate))
      addGN       <- !na.gate    && length(gcov) > 0    && G > 1       && isTRUE(noNgate)  && (is.null(o.gate) || !identical(gating, if(fullMoE) o.nets else o.gate))
      addNG       <- G  > 1      && isTRUE(noNgate)     && stepG       && (has.init        &&  attr(initialModel, "Gating"))
      addQ        <- G  > 0      && isTRUE(equalPro)    && stepG       && (!has.init       || !attr(initialModel, "Gating"))
      addB        <- fullMoE     && length(bcov) > 0    && G > 0       && !na.both         && (is.null(o.nets) || !identical(b.nets, o.nets))
      if(isTRUE(allPro)          &&
        (isFALSE(has.init)       || 
         !attr(initialModel, 
               "Gating"))        &&
         G         > 1)           {
        if(isTRUE(has.init)      &&
           attr(initialModel, 
                 "EqualPro"))     {
          args$equalPro          <- FALSE
          temp    <- do.call(MoE_clust, args)
          comp    <- suppressWarnings(MoE_compare(res, temp, optimal.only=TRUE, pick=1L, criterion=criterion))
          if(comp$MoENames == "temp")      {
            j     <- j  + 1L
            if(isTRUE(verbose))            {      message("\n\tStep ", j, "...\n")
              suppressWarnings(print(comp, details=FALSE))
            }
          }
          compare[[j]]           <- res   <- comp$optimal
          crit    <- switch(EXPR=criterion, bic=res$bic, icl=res$icl, aic=res$aic)
        } else if((isFALSE(has.init)      ||
                  !attr(initialModel, 
                        "EqualPro")))      {
          args$equalPro          <- TRUE
          temp    <- do.call(MoE_clust, args)
          args$equalPro          <- FALSE
          comp    <- suppressWarnings(MoE_compare(res, temp, optimal.only=TRUE, pick=1L, criterion=criterion))
          if(comp$MoENames == "temp")      {
            j     <- j  + 1L  
            if(isTRUE(verbose))            {      message("\n\tStep ", j, "...\n")
              suppressWarnings(print(comp, details=FALSE))
            }
          }
          compare[[j]]           <- res   <- comp$optimal
          crit    <- switch(EXPR=criterion, bic=res$bic, icl=res$icl, aic=res$aic)
        }
      }
      if(isTRUE(allNgate)        &&
         isTRUE(has.init)        &&
         attr(initialModel,   "Noise")    &&
         attr(initialModel,  "Gating"))    {
        if(!attr(initialModel, 
                 "NoiseGate"))    {
          args$noise.gate        <- TRUE
          temp    <- do.call(MoE_clust, args)
          comp    <- suppressWarnings(MoE_compare(res, temp, optimal.only=TRUE, pick=1L, criterion=criterion))
          if(comp$MoENames == "temp")      {
            j     <- j  + 1L
            if(isTRUE(verbose))            {      message("\n\tStep ", j, "...\n")
              suppressWarnings(print(comp, details=FALSE))
            }
          }
          compare[[j]]           <- res   <- comp$optimal
          crit    <- switch(EXPR=criterion, bic=res$bic, icl=res$icl, aic=res$aic)
        } else if(G     > 1)      { 
          args$noise.gate        <- FALSE
          temp    <- do.call(MoE_clust, args)
          args$noise.gate        <- TRUE
          comp    <- suppressWarnings(MoE_compare(res, temp, optimal.only=TRUE, pick=1L, criterion=criterion))
          if(comp$MoENames == "temp")      {
            j     <- j  + 1L  
            if(isTRUE(verbose))            {      message("\n\tStep ", j, "...\n")
              suppressWarnings(print(comp, details=FALSE))
            }
          }
          compare[[j]]           <- res   <- comp$optimal
          crit    <- switch(EXPR=criterion, bic=res$bic, icl=res$icl, aic=res$aic)
        }
      }
    }         else {
      addC        <- stepG
      addE        <- isFALSE(noise)   && !na.expx && length(ecov)     > 0
      addQ        <- isTRUE(equalPro) && stepG    && (isFALSE(noise) || (!is.null(args$equalNoise) && isTRUE(args$equalNoise)))
      addG        <- 
      addGN       <- 
      addNG       <- FALSE
      addB        <- fullMoE          && addE
    }
    
    while(crit     > crit.old)    {
      crit.old    <- crit
      resold      <- res
      j           <- j  + 1L
      if(isTRUE(verbose))                         message("\n\tStep ", j, "...\n")
      if(!fullMoE) {
        if(addE)   {
          args$init.z            <- args$z.init
          args$z.list            <- NULL
          exps    <- try(lapply(seq_along(ecov), function(x) { args$expert <- stats::as.formula(ecov[x]); suppressWarnings(do.call(MoE_clust, args)) }), silent=TRUE)
          if(inherits(exps,    "try-error")) {
            exps  <- NULL
          }   else   exps        <- suppressWarnings(MoE_compare(stats::setNames(exps,  seq_along(exps)),  optimal.only=TRUE, pick=1L, criterion=criterion)$optimal)
          args$init.z            <- "list"
          args$z.list            <- INIT
        }     else   exps        <- NULL
        if(addG)   {
          gates   <- try(lapply(seq_along(gcov), function(x) { args$gating <- stats::as.formula(gcov[x]); suppressWarnings(do.call(MoE_clust, args)) }), silent=TRUE)
          if(inherits(gates,   "try-error")) {
            gates <- NULL
          }   else   gates       <- suppressWarnings(MoE_compare(stats::setNames(gates, seq_along(gates)), optimal.only=TRUE, pick=1L, criterion=criterion)$optimal)
        }     else   gates       <- NULL
        if(addGN)  {
          args$noise.gate        <- FALSE
          Ngate   <- try(lapply(seq_along(gcov), function(x) { args$gating <- stats::as.formula(gcov[x]); suppressWarnings(do.call(MoE_clust, args)) }), silent=TRUE)
          if(inherits(Ngate,   "try-error")) {
            Ngate <- NULL
          }   else   Ngate       <- suppressWarnings(MoE_compare(stats::setNames(Ngate, seq_along(Ngate)), optimal.only=TRUE, pick=1L, criterion=criterion)$optimal)
          args$noise.gate        <- TRUE
        }     else   Ngate       <- NULL
      }       else if(addB)       {
        args$init.z              <- args$z.init
        args$z.list              <- NULL
        if(addE   &&
           !any(addG, addGN))     {
          exps    <- try(lapply(seq_along(ecov), function(x) { args$expert <- stats::as.formula(ecov[x]); suppressWarnings(do.call(MoE_clust, args)) }), silent=TRUE)
          if(inherits(exps,    "try-error")) {
            exps  <- NULL
          }   else   exps        <- suppressWarnings(MoE_compare(stats::setNames(exps,  seq_along(exps)),  optimal.only=TRUE, pick=1L, criterion=criterion)$optimal)
          gates   <-
          Ngate   <- NULL
        }     else {
         exps     <- NULL
         if(addG)  {
           gates  <- try(lapply(seq_along(bcov), function(x) { args$expert <- 
                                                               args$gating <- stats::as.formula(bcov[x]); suppressWarnings(do.call(MoE_clust, args)) }), silent=TRUE)
           if(inherits(gates,  "try-error")) {
            gates <- NULL
           }  else   gates       <- suppressWarnings(MoE_compare(stats::setNames(gates, seq_along(gates)), optimal.only=TRUE, pick=1L, criterion=criterion)$optimal)
         }    else   gates       <- NULL
         if(addGN) {
           args$noise.gate       <- FALSE
           Ngate  <- try(lapply(seq_along(bcov), function(x) { args$expert <- 
                                                               args$gating <- stats::as.formula(bcov[x]); suppressWarnings(do.call(MoE_clust, args)) }), silent=TRUE)
           if(inherits(Ngate,  "try-error")) {
            Ngate <- NULL
           }  else   Ngate       <- suppressWarnings(MoE_compare(stats::setNames(Ngate, seq_along(Ngate)), optimal.only=TRUE, pick=1L, criterion=criterion)$optimal)
           args$noise.gate       <- TRUE
         }    else   Ngate       <- NULL
        }
        args$init.z              <- "list"
        args$z.list              <- INIT
      }       else {
        exps      <- 
        gates     <-
        Ngate     <- NULL
      }
      if(fullMoE)  {
        if(addC)   {
         args$init.z             <- args$z.init
         args$z.list             <- NULL
         args$G   <- G  + 1L
         if(G     == 0           ||
            args$G > 
            args$noise.gate)      {
           args$gating           <- args$expert
           clplus <- try(suppressWarnings(do.call(MoE_clust, args)), silent=TRUE)
           if(inherits(clplus, "try-error")) {
             clplus              <- NULL    
           }  else   clplus      <- suppressWarnings(MoE_compare(clplus, optimal.only=TRUE, pick=1L, criterion=criterion)$optimal)
         }    else   clplus      <- NULL
         if(addNG) {
           args$noise.gate       <- FALSE
           Nplus  <- try(suppressWarnings(do.call(MoE_clust, args)), silent=TRUE)
           if(inherits(Nplus,  "try-error")) {
            Nplus <- NULL    
           }  else   Nplus       <- suppressWarnings(MoE_compare(Nplus,  optimal.only=TRUE, pick=1L, criterion=criterion)$optimal)
           args$noise.gate       <- TRUE
         }    else   Nplus       <- NULL
         args$G   <- G
         args$init.z             <- "list"
         args$z.list             <- INIT
        }     else {
          clplus  <- 
          Nplus   <- NULL
        }
      }       else {
        if(addC)   { 
          args$init.z            <- args$z.init
          args$z.list            <- NULL
          args$G  <- G  + 1L
          clplus  <- try(suppressWarnings(do.call(MoE_clust, args)), silent=TRUE)
          if(inherits(clplus,  "try-error")) {
           clplus <- NULL    
          }   else   clplus      <- suppressWarnings(MoE_compare(clplus, optimal.only=TRUE, pick=1L, criterion=criterion)$optimal)
          if(addNG)  {
            args$noise.gate      <- FALSE
            Nplus <- try(suppressWarnings(do.call(MoE_clust, args)), silent=TRUE)
            if(inherits(Nplus, "try-error")) {
              Nplus              <- NULL    
            } else   Nplus       <- suppressWarnings(MoE_compare(Nplus,  optimal.only=TRUE, pick=1L, criterion=criterion)$optimal)
            args$noise.gate      <- TRUE
          }   else   Nplus       <- NULL
          args$G  <- G
          args$init.z            <- "list"
          args$z.list            <- INIT
        }     else {
          clplus  <- 
          Nplus   <- NULL
        }
      }
      if(addQ)     {
        args$init.z              <- args$z.init
        args$z.list              <- NULL
        args$equalPro            <- TRUE
        args$G    <- G  + 1L
        eqplus    <- try(suppressWarnings(do.call(MoE_clust, args)), silent=TRUE)
        if(inherits(eqplus,    "try-error")) {
          eqplus  <- NULL
        }     else   eqplus      <- suppressWarnings(MoE_compare(eqplus, optimal.only=TRUE, pick=1L, criterion=criterion)$optimal)
        args$G    <- G
        if(addE   && G != 1      &&
           !fullMoE)              {
          expQ    <- try(lapply(seq_along(ecov), function(x) { args$expert <- stats::as.formula(ecov[x]); suppressWarnings(do.call(MoE_clust, args)) }), silent=TRUE)
          if(inherits(expQ,    "try-error")) {
            expQ  <- NULL
          }   else   expQ        <- suppressWarnings(MoE_compare(stats::setNames(expQ,  seq_along(expQ)),  optimal.only=TRUE, pick=1L, criterion=criterion)$optimal)
        }     else   expQ        <- NULL
        args$equalPro            <- FALSE
        args$init.z              <- "list"
        args$z.list              <- INIT
      }       else {
        eqplus    <-
        expQ      <- NULL
      }
      models      <- c(list(exps), list(gates), list(clplus), list(eqplus), list(expQ), list(Ngate), list(Nplus), list(resold))
      models      <- models[!vapply(models, is.null, logical(1L))]
      names(models)              <- seq_along(models)
      comp        <- suppressWarnings(MoE_compare(models, optimal.only=TRUE, pick=1L, criterion=criterion))
      if(isTRUE(verbose))           suppressWarnings(print(comp, details=FALSE))
      compare[[j]]     <- res    <- comp$optimal
      crit        <- switch(EXPR=criterion, bic=res$bic, icl=res$icl, aic=res$aic)
      if(res$G    != args$G)      {
        args$G    <- G <- res$G 
        args$init.z    <- "list" 
        args$z.list    <- INIT   <- list(attr(res, "Exp.init"))
      }
      args$expert      <- stats::as.formula(attr(res$expert, "Formula"))
      args$gating      <- stats::as.formula(attr(res$gating, "Formula"))
      if(length(ecov)   > 1      &&
        (any(any.e     <- is.element(ecov, 
        paste0("~", as.character(args$expert)[-1L]))))) {
        if(is.null(o.exps))       {
          ecov    <- paste0("~", paste(all.vars(args$expert), gsub('\\~', '', ecov[!any.e]), sep=" + "))  
        }     else {
          ecov    <- paste0("~", paste(paste(o.exps, collapse=" + "), setdiff(all.vars(args$expert), o.exps), sep=" + "), 
                            gsub(paste(o.exps, collapse=" \\+ "), '', gsub('\\~', '', ecov[!any.e])))
        }
        ecov      <- gsub(" $", "", gsub("\\+$", "", gsub(" $", "", ecov)))
        o.exps    <- all.vars(args$expert)
      }
      if(length(gcov)   > 1      &&
        (any(any.g     <- is.element(gcov, 
        paste0("~", as.character(args$gating)[-1L]))))) {
        if(is.null(o.gate))       {
          gcov    <- paste0("~", paste(all.vars(args$gating), gsub('\\~', '', gcov[!any.g]), sep=" + "))
        }     else { 
          gcov    <- paste0("~", paste(paste(o.gate, collapse=" + "), setdiff(all.vars(args$gating), o.gate), sep=" + "), 
                            gsub(paste(o.gate, collapse=" \\+ "), '', gsub('\\~', '', gcov[!any.g])))
        }
        gcov      <- gsub(" $", "", gsub("\\+$", "", gsub(" $", "", gcov)))
        o.gate    <- all.vars(args$gating)
      }
      if(isTRUE(fullMoE)         &&
         length(bcov)   > 1      &&
        (any(any.b     <- is.element(bcov, 
        paste0("~", as.character(args$gating)[-1L]))))) {
        if(is.null(o.nets))       {
          bcov    <- paste0("~", paste(all.vars(args$gating), gsub('\\~', '', bcov[!any.b]), sep=" + "))
        }     else { 
          bcov    <- paste0("~", paste(paste(o.nets, collapse=" + "), setdiff(all.vars(args$gating), o.nets), sep=" + "), 
                            gsub(paste(o.nets, collapse=" \\+ "), '', gsub('\\~', '', bcov[!any.b])))
        }
        bcov      <- gsub(" $", "", gsub("\\+$", "", gsub(" $", "", bcov)))
        o.nets    <- all.vars(args$gating)
      }
      addC        <- stepG       && (isTRUE(unequalPro)       || attr(res, "Gating"))
      addE        <- G  > 0      && length(ecov) > 0
      addQ        <- G  > 0      && isTRUE(equalPro) && stepG && !attr(res, "Gating")
      addG        <- G  > !noise && length(gcov) > 0 && isTRUE(doNgate)
      addGN       <- G  > 1      && length(gcov) > 0 && isTRUE(noNgate)
      addNG       <- G  > 1      && isTRUE(noNgate)  && stepG &&  attr(res, "Gating")
      addQ        <- addQ        && !all(fullMoE, attr(res, "Expert"))
      addB        <- G  > 0      && length(bcov) > 0
      compare[[j]]$call          <- NULL
    }
    compare       <- compare[-j]
    for(i in seq_along(compare))  {
      compare[[i]]$call          <- call
      names(compare[[i]]$data)   <- colnames(data)
    }
    names(compare)               <- paste0("Step_", seq_along(compare))
    cat("\n")
    res           <- suppressWarnings(MoE_compare(compare, optimal.only=TRUE, criterion=criterion))
    if(any(res$posidens))                         warning("Potentially spurious solutions with positive log-densities were chosen at one or more steps\n", call.=FALSE)
      return(res)
  }
  
#' Convert MoEClust objects to the Mclust class
#'
#' Converts an object of class \code{"MoEClust"} generated by \code{\link{MoE_clust}} and converts it to an object of class \code{"Mclust"} as generated by fitting \code{\link[mclust]{Mclust}}, to facilitate use of plotting and other functions for the \code{"Mclust"} class within the \pkg{mclust} package. Some caution is advised when converting models with gating &/or expert covariates (see Note below).
#' @param x An object of class \code{"MoEClust"} generated by \code{\link{MoE_clust}} or an object of class \code{"MoECompare"} generated by \code{\link{MoE_compare}}. Models with a noise component are facilitated here too.
#' @param expert.covar Logical (defaults to \code{TRUE}) governing whether the extra variability in the component means is added to the MVN ellipses corresponding to the component covariance matrices in the presence of expert network covariates. See the function \code{\link{expert_covar}}.
#' @param signif Significance level for outlier removal. Must be a single number in the interval [0, 1). Corresponds to the percentage of data to be considered extreme and therefore removed (half of \code{signif} at each endpoint, on a column-wise basis). The default, \code{0}, corresponds to no outlier removal. \strong{Only} invoke this argument as an aid to visualisation via \code{\link[mclust]{plot.Mclust}}.
#' @param ... Further arguments to be passed to other methods.
#'
#' @return An object of class \code{"Mclust"}. See \code{methods(class="Mclust")} for a (non-exhaustive) list of functions which can be applied to this class.
#' @details Of course, the user is always encouraged to use the dedicated \code{\link[=plot.MoEClust]{plot}} function for objects of the \code{"MoEClust"} class instead, but calling \code{plot} after converting via \code{\link[=as.Mclust.MoEClust]{as.Mclust}} can be particularly useful for univariate mixtures.
#'
#' In the presence of expert network covariates, the component-specific covariance matrices are (by default, via the argument \code{expert.covar}) modified for plotting purposes via the function \code{\link{expert_covar}}, in order to account for the extra variability of the means, usually resulting in bigger shapes & sizes for the MVN ellipses.
#'
#' The \code{signif} argument is intended only to aid visualisation via \code{\link[mclust]{plot.Mclust}}, as plots therein can be sensitive to outliers, particularly with regard to axis limits.
#' @note Mixing proportions are averaged over observations in components in the presence of gating network covariates during the coercion.
#'
#' Plots may be quite misleading in the presence of gating &/or (especially) expert network covariates when the \code{what} argument is \code{"density"} within \code{\link[mclust]{plot.Mclust}}; users are \strong{strongly} encouraged to use \code{\link{MoE_gpairs}} with \code{response.type="density"} instead.
#' 
#' Predictions (via \code{\link[mclust]{predict.Mclust}}) will also be misleading in the presence of covariates of any kind when \code{newdata} is supplied; thus, users are \strong{strongly} encouraged to use \code{\link{predict.MoEClust}} instead. 
#'
#' The functions \code{\link[mclust]{clustCombi}} and \code{\link[mclust]{clustCombiOptim}} can be safely used (provided \code{as.Mclust(x)} is supplied as the \code{object} argument to \code{\link[mclust]{clustCombi}}), as they only rely on \code{x$z} and \code{x$G} only. See the examples below.
#' 
#' Users may expect MoEClust models with no covariates of any kind to be identical to models fitted via \pkg{mclust}, but this is not necessarily true: see the \code{\link{MoE_control}} argument \code{asMclust}.
#' @importFrom matrixStats "colMeans2"
#' @importFrom mclust "as.Mclust" "clustCombi" "clustCombiOptim" "logLik.Mclust" "icl" "plot.Mclust" "plot.mclustBIC" "plot.mclustICL" "predict.Mclust" "print.Mclust" "sigma2decomp" "summary.Mclust"
#' @method as.Mclust MoEClust
#' @seealso \code{\link[mclust]{Mclust}}, \code{\link[mclust]{plot.Mclust}}, \code{\link{MoE_clust}}, \code{\link{plot.MoEClust}}, \code{\link{expert_covar}}, \code{\link{MoE_control}}
#' @author Keefe Murphy - <\email{keefe.murphy@@mu.ie}>
#' @references Fraley, C. and Raftery, A. E. (2002). Model-based clustering, discriminant analysis, and density estimation. \emph{Journal of the American Statistical Association}, 97(458): 611-631.
#' 
#' Scrucca L., Fop M., Murphy T. B. and Raftery A. E. (2016). mclust 5: clustering, classification and density estimation using Gaussian finite mixture models. \emph{The R Journal}, 8(1): 289-317.
#' @keywords utility
#' @usage
#' \method{as.Mclust}{MoEClust}(x,
#'          expert.covar = TRUE,
#'          signif = 0L,
#'          ...)
#' @export
#' @name as.Mclust
#' @examples
#' \donttest{# library(mclust)
#' 
#' # Fit a gating network mixture of experts model to the ais data
#' # data(ais)
#' # mod   <- MoE_clust(ais[,3:7], G=1:9, gating= ~ BMI + sex, network.data=ais)
#'
#' # Convert to the "Mclust" class and examine the classification
#' # mod2  <- as.Mclust(mod)
#' # plot(mod2, what="classification")
#'
#' # Examine the uncertainty
#' # plot(mod2, what="uncertainty")
#'
#' # Return the optimal number of clusters according to entropy
#' # combi <- mclust::clustCombi(object=mod2)
#' # optim <- mclust::clustCombiOptim(object=combi)
#' # table(mod2$classification, ais$sex)
#' # table(optim$cluster.combi, ais$sex)
#'
#' # While we could have just used plot.MoEClust above,
#' # plot.Mclust is especially useful for univariate data
#' # data(CO2data)
#' # res <- MoE_clust(CO2data$CO2, G=3, equalPro=TRUE, expert = ~ GNP, network.data=CO2data)
#' # plot(as.Mclust(res))}
  as.Mclust.MoEClust      <- function(x, expert.covar = TRUE, signif = 0L, ...) {
    x             <- if(inherits(x, "MoECompare")) x$optimal else x
    if(length(signif) > 1 || !is.numeric(signif)     ||
       signif < 0 || signif     >= 1)             stop("'signif' must be a single number in the interval [0, 1)", call.=FALSE)
    uni           <- x$d  == 1
    gating        <- attr(x, "Gating")
    expert        <- attr(x, "Expert")
    x$data        <- as.matrix(x$data)
    x$loglik      <- x$loglik[length(x$loglik)]
    x$BIC         <- replace(x$BIC, !is.finite(x$BIC), NA)
    class(x$BIC)  <- "mclustBIC"
    x$uncertainty         <- if(uni)                      unname(x$uncertainty) else x$uncertainty
    x$classification      <- if(uni)                   unname(x$classification) else x$classification
    x$parameters$pro      <- if(gating)             colMeans2(x$parameters$pro) else x$parameters$pro
    x$parameters$variance <- if(isTRUE(expert.covar) &&
                                expert)  suppressWarnings(expert_covar(x, ...)) else x$parameters$variance
    x$data        <- if(signif   > 0)      apply(x$data, 2L, .trim_out, signif) else x$data
    x$modelName   <- ifelse(x$G == 0, "EII", x$modelName)
    colnames(x$z) <- NULL
    x             <- x[-which(is.element(names(x), c("ICL", "AIC", "aic", "gating", "expert", "LOGLIK", "linf", "iters", "net.covs", "resid.data", "DF", "ITERS")))]
    name.x        <- names(x)
    attributes(x) <- NULL
    names(x)      <- name.x
    x             <- x[c("call", "data", "modelName", "n", "d", "G", "BIC", "loglik", "df", "bic", "icl", "hypvol", "parameters", "z", "classification", "uncertainty")]
    class(x)      <- "Mclust"
      x
  }

#' Account for extra variability in covariance matrices with expert covariates
#'
#' In the presence of expert network covariates, this helper function modifies the component-specific covariance matrices of a \code{"MoEClust"} object, in order to account for the extra variability due to the component means, usually resulting in bigger shapes & sizes for the MVN ellipses in \code{\link{MoE_gpairs}} plots. The function also works for univariate response data.
#' @param x An object of class \code{"MoEClust"} generated by \code{\link{MoE_clust}}, or an object of class \code{"MoECompare"} generated by \code{\link{MoE_compare}}. Models with a noise component are facilitated here too.
#' @param weighted A logical indicating whether the estimated cluster membership probabilities should be used to provide a weighted estimate of the variability due to the component means. Defaults to \code{TRUE}. The option \code{weighted=FALSE} is provided only so that previous behaviour under earlier versions of this package can be recovered but is otherwise not recommended.
#' @param ... Catches unused arguments.
#' 
#' @details This function is used internally by \code{\link{MoE_gpairs}}, \code{\link{plot.MoEClust}(x, what="gpairs")}, and \code{\link[=as.Mclust.MoEClust]{as.Mclust}}, for visualisation purposes.
#' @note The \code{modelName} of the resulting \code{variance} object may not correspond to the model name of the \code{"MoEClust"} object, in particular scale, shape, &/or orientation may no longer be constrained across clusters. Usually, the \code{modelName} of the transformed \code{variance} object will be \code{"VVV"}.
#' @return The \code{variance} component only from the \code{parameters} list from the output of a call to \code{\link{MoE_clust}}, modified accordingly.
#' @seealso \code{\link{MoE_clust}}, \code{\link{MoE_gpairs}}, \code{\link{plot.MoEClust}}, \code{\link[=as.Mclust.MoEClust]{as.Mclust}}
#' @references Murphy, K. and Murphy, T. B. (2020). Gaussian parsimonious clustering models with covariates and a noise component. \emph{Advances in Data Analysis and Classification}, 14(2): 293-325. <\doi{10.1007/s11634-019-00373-8}>.
#' @author Keefe Murphy - <\email{keefe.murphy@@mu.ie}>
#' @keywords utility
#' @usage
#' expert_covar(x,
#'              weighted = TRUE,
#'              ...)
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
  expert_covar    <- function(x, weighted = TRUE, ...) {
      UseMethod("expert_covar")
  }

#' @method expert_covar MoEClust
#' @importFrom mclust "covw" "sigma2decomp"
#' @export
  expert_covar.MoEClust   <- function(x, weighted = TRUE, ...) {
    x             <- if(inherits(x, "MoECompare")) x$optimal else x
    x.sig         <- x$parameters$variance
    G             <- x$G
    if(length(weighted)   != 1   ||
       !is.logical(weighted))                     stop("'weighted' must be a single logical indicator", call.=FALSE)
    weighted      <- weighted    && G != 1
    if(attr(x, "Expert"))  {
      if(x$d      == 1)    {
        if(isTRUE(weighted))      {
          predVar <- drop(covw(fitted.MoEClust(x), x$z[,seq_len(G)], normalize=TRUE)$S)
          x.sig$sigmasq   <- x.sig$sigmasq + sqrt(predVar)
        } else     {
          n       <- x$n
          predVar <- vapply(x$expert, function(expert) stats::cov(as.matrix(expert$fitted.values)), numeric(1L)) * (n - 1L)/n
          x.sig$sigmasq   <- unname(x.sig$sigmasq + sqrt(predVar))
        }
        if(x$modelName    == "V" || length(unique(x.sig$sigmasq)) > 1) {
          x.sig$scale     <- x.sig$sigmasq
          x.sig$modelName <- "V"
        }
      }   else     {
        if(isTRUE(weighted))      {
          predVar <- covw(fitted.MoEClust(x), x$z[,seq_len(G)], normalize=TRUE)$S
        } else     {
          n       <- x$n
          predVar <- sapply(x$expert, function(expert) stats::cov(expert$fitted.values), simplify="array") * (n - 1L)/n
        }
        x.sig     <- suppressWarnings(sigma2decomp(x.sig$sigma + predVar))
      }
    }     else                                    message("No expert covariates: returning the variance object without modification\n")
      return(x.sig)
  }

#' Force diagonal elements of a triangular matrix to be positive
#'
#' This function ensures that the triangular matrix in a QR (or other) decomposition has positive values along its diagonal.
#' @param x A matrix, which must be either upper-triangular or lower-triangular.
#'
#' @return An upper or lower triangular matrix with positive diagonal entries such that the matrix is still a valid decomposition of the matrix the input \code{x} is a decomposition of.
#' @export
#' @author Keefe Murphy - <\email{keefe.murphy@@mu.ie}>
#' @keywords utility
#'
#' @examples
#' data(ais)
#' res <- MoE_clust(ais[,3:7], G=3, modelNames="EEE")
#' sig <- res$parameters$variance
#' a   <- force_posiDiag(sig$cholSigma)
#' b   <- chol(sig$Sigma)
#' all.equal(a, b)                    #TRUE
#' all.equal(crossprod(a), sig$Sigma) #TRUE
#' all.equal(crossprod(b), sig$Sigma) #TRUE
  force_posiDiag  <- function(x) {
    if(!is.matrix(x) ||
       any(x[row(x)   > col(x)] != 0) &&
       any(x[col(x)   > row(x)] != 0))            stop("'x' must be an upper or lower triangular matrix")
      provideDimnames(diag(sign(diag(x))) %*% x, base=list(rownames(x), colnames(x)))
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
  quant_clust     <- function(x, G)  {
    if(NCOL(x)     > 1  || !all(is.numeric(x)))   stop("'x' must be univariate", call.=FALSE)
    x             <- as.vector(x)
    eps           <- stats::sd(x) * sqrt(.Machine$double.eps)
    q             <- NA
    n             <- G
    while(length(q) < (G + 1L))      {
      n           <- n   + 1L
      q           <- unique(stats::quantile(x, seq(from=0L, to=1L, length=n)))
    }
    if(length(q)   > (G  + 1L))      {
      q           <- q[-order(diff(q))[seq_len(length(q) - G - 1L)]]
    }
    q[1L]         <- min(x)  - eps
    q[length(q)]  <- max(x)  + eps
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
#' @author Keefe Murphy - <\email{keefe.murphy@@mu.ie}>
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
#' Drops unseen factor levels in \code{newdata} for which predictions are required from a \code{\link[stats]{lm}} or \code{\link[nnet]{multinom}} model \code{fit}.
#' @param fit A fitted \code{\link[stats]{lm}} or \code{\link[nnet]{multinom}} model.
#' @param newdata A \code{data.frame} containing variables with which to predict.
#'
#' @return A \code{data.frame} like \code{newdata} with unseen factor levels replaced by \code{NA}.
#' @note This function is so far untested for models other than \code{\link[stats]{lm}} or \code{\link[nnet]{multinom}}, though it \emph{may} still work for other classes.
#' @author Keefe Murphy - <\email{keefe.murphy@@mu.ie}>
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
#' Computes the Mahalanobis distance between the response variable(s) and the fitted values of linear regression models with multivariate or univariate responses.
#' @param fit A fitted \code{\link[stats]{lm}} model, inheriting either the \code{"mlm"} or \code{"lm"} class.
#' @param resids The residuals. Can be residuals for observations included in the model, or residuals arising from predictions on unseen data. Must be coercible to a matrix with the number of columns being the number of response variables. Missing values are not allowed.
#' @param squared A logical. By default (\code{FALSE}), the generalized interpoint distance is computed. Set this flag to \code{TRUE} for the squared value.
#' @param identity A logical indicating whether the identity matrix is used in in place of the precision matrix in the Mahalanobis distance calculation. Defaults to \code{FALSE}; \code{TRUE} corresponds to the use of the Euclidean distance. Only relevant for multivariate response data.
#'
#' @return A vector giving the Mahalanobis distance (or squared Mahalanobis distance) between response(s) and fitted values for each observation.
#' @author Keefe Murphy - <\email{keefe.murphy@@mu.ie}>
#' @importFrom matrixStats "rowSums2"
#' @keywords utility
#' @export
#' @usage
#' MoE_mahala(fit,
#'            resids,
#'            squared = FALSE,
#'            identity = FALSE)
#' @examples
#' \dontshow{library(matrixStats)}
#' data(ais)
#' hema <- as.matrix(ais[,3:7])
#' mod  <- lm(hema ~ sex + BMI, data=ais)
#' res  <- hema - predict(mod)
#' MoE_mahala(mod, res, squared=TRUE)
#' 
#' \donttest{data(CO2data)
#' CO2  <- CO2data$CO2
#' GNP  <- CO2data$GNP
#' mod2 <- lm(CO2 ~ GNP, data=CO2data)
#' pred <- predict(mod2)
#' res2 <- CO2 - pred
#' maha <- MoE_mahala(mod2, res2)
#' 
#' # Highlight outlying observations
#' plot(GNP, CO2, type="n", ylab=expression('CO'[2]))
#' lines(GNP, pred, col="red")
#' points(GNP, CO2, cex=maha, lwd=2)
#' text(GNP, CO2, col="blue", 
#'      labels=replace(as.character(CO2data$country), maha < 1, ""))
#'      
#' # Replicate initialisation strategy using 2 randomly chosen components
#' # Repeat the random initialisation if necessary
#' # (until 'crit' at convergence is minimised)
#' G       <- 3L
#' z       <- sample(seq_len(G), nrow(CO2data), replace=TRUE)
#' old     <- Inf
#' crit    <- .Machine$double.xmax
#' while(crit < old)   {
#'   Sys.sleep(1)
#'   old   <- crit
#'   maha  <- NULL
#'   plot(GNP, CO2, type="n", ylab=expression('CO'[2]))
#'   for(g in seq_len(G)) { 
#'    ind  <- which(z == g)
#'    mod  <- lm(CO2 ~ GNP, data=CO2data, sub=ind)
#'    pred <- predict(mod, newdata=CO2data[,"CO2", drop=FALSE])
#'    maha <- cbind(maha, MoE_mahala(mod, CO2 - pred))
#'    lines(GNP, pred, col=g + 1L)
#'   }
#'   min.M <- rowMins(maha)
#'   crit  <- sum(min.M)
#'   z     <- max.col(maha == min.M)
#'   points(GNP, CO2, cex=min.M, lwd=2, col=z + 1L)
#'   text(GNP, CO2, col=z + 1L, 
#'        labels=replace(as.character(CO2data$country), which(min.M <= 1), ""))
#' }
#' crit}
  MoE_mahala      <- function(fit, resids, squared = FALSE, identity = FALSE) {
    if(!inherits(fit, "mlm") &&
       !inherits(fit, "lm"))                      stop("'fit' must inherit the class \"mlm\" or \"lm\"",  call.=FALSE)
    resids        <- tryCatch(data.matrix(as.data.frame(resids)), error=function(e) {
                                                  stop("Invalid 'resids': must be coercible to a matrix", call.=FALSE) })
    if(!is.numeric(resids)   ||
       anyNA(resids))                             stop("Invalid 'resids': must be numeric and contain no missing values", call.=FALSE)
    if(length(squared) > 1   ||
       !is.logical(squared))                      stop("'squared' must be a single logical indicator",    call.=FALSE)
    if(length(identity) > 1  ||
       !is.logical(identity))                     stop("'identity' must be a single logical indicator",   call.=FALSE)
    
    if(inherits(fit, "mlm"))  {
     if(isTRUE(identity))     {
       icov       <- diag(ncol(resids))
     } else        { 
      covar       <- crossprod(resids)/(nrow(resids) - fit$rank)
      covar[!stats::complete.cases((covar))]   <- .Machine$double.eps
      if(diff(dim(resids))   >= 0)     {
        covsvd    <- svd(covar)
        posi      <- covsvd$d > max(sqrt(.Machine$double.eps) * covsvd$d[1L], 0L)
        icov      <- if(all(posi)) covsvd$v   %*% (t(covsvd$u)/covsvd$d) else if(any(posi))
        covsvd$v[,posi, drop=FALSE]   %*% (t(covsvd$u[,posi, drop=FALSE])/covsvd$d[posi]) else array(0L, dim(covar)[2L:1L])
      } else icov <- chol2inv(.chol(covar))
     }
     res          <- rowSums2(resids  %*% icov * resids)
      return(drop(if(isTRUE(squared)) res else sqrt(res)))
    }   else       {
    #covar        <- as.numeric(crossprod(resids)/(nrow(resids) - fit$rank))
     covar        <- 1L
      return(drop(if(isTRUE(squared)) (resids  * resids)/covar else abs(resids)/covar))
    }
  }

#' Approximate Hypervolume Estimate
#'
#' Computes simple approximations to the hypervolume of univariate and multivariate data sets. Also returns the location of the centre of mass.
#' @param data A numeric vector, matrix, or data frame of observations. Categorical variables are not allowed, and covariates should not be included. If a matrix or data frame, rows correspond to observations and columns correspond to variables. There \strong{must} be more observations than variables.
#' @param method The method used to estimate the hypervolume. The default method uses the function \code{\link[mclust]{hypvol}}. The \code{"convexhull"} and \code{"ellipsoidhull"} options require loading the \code{geometry} and \code{cluster} libraries, respectively. This argument is only relevant for multivariate data; for univariate data, the range of the data is used.
#' @param reciprocal A logical variable indicating whether or not the reciprocal hypervolume is desired rather than the hypervolume itself. The default is to return the hypervolume.
#'
#' @importFrom matrixStats "colMeans2" "colRanges" "rowDiffs" "rowMeans2"
#' @importFrom mclust "hypvol"
#' @note This function is called when adding a noise component to \code{MoEClust} models via the function \code{MoE_control}, specifically it's argument \code{noise.meth}. The function internally only uses the response variables, and not the covariates. However, one can bypass the invocation of this function by specifying its \code{noise.vol} argument directly. This is explicitly necessary for models for high-dimensional data which include a noise component for which this function cannot estimate a (hyper)volume.
#' 
#' Note that supplying the volume manually to \code{\link{MoE_clust}} can affect the summary of the means in \code{parameters$mean} and by extension the location of the MVN ellipses in \code{\link{MoE_gpairs}} plots for models with \emph{both} expert network covariates and a noise component. The location cannot be estimated when the volume is supplied manually; in this case, prediction is made on the basis of renormalising the \code{z} matrix after discarding the column corresponding to the noise component. Otherwise, the mean of the noise component is accounted for. The renormalisation approach can be forced by specifying \code{noise.args$discard.noise=TRUE}, even when the mean of the noise component is available.
#' @return A list with the following two elements:
#' \describe{
#' \item{\code{vol}}{A hypervolume estimate (or its inverse). 
#' 
#' This can be used as the hypervolume parameter for the noise component when observations are designated as noise in \code{\link{MoE_clust}}.}
#' \item{\code{loc}}{A vector of length \code{ncol(data)} giving the location of the centre of mass.
#' 
#' This can help in predicting the fitted values of models fitted with noise components via \code{\link{MoE_clust}}.}}
#' @author Keefe Murphy - <\email{keefe.murphy@@mu.ie}>
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
    has.lib       <- switch(EXPR=method, hypvol=TRUE, 
                            convexhull=suppressMessages(requireNamespace("geometry",   quietly=TRUE)) && .version_above("geometry", "0.4.0"), 
                            ellipsoidhull=suppressMessages(requireNamespace("cluster", quietly=TRUE)) && .version_above("cluster",  "1.4.0"))
    if(diff(dim(data))     > 0)                   stop("The hypervolume can only be estimated when there are more observations than variables, otherwise it should be specified as an independent tuning parameter",                   call.=FALSE)
    if(!has.lib)                                  stop(paste0("Use of the ", method, " 'method' option requires loading the ", switch(EXPR=method, hypvol="'mclust'", convexhull="'geometry'", ellipsoidhull="'cluster'"), "library"), call.=FALSE)
    if(length(reciprocal) != 1 ||
       !is.logical(reciprocal))                   stop("'reciprocal' must be a single logical indicator", call.=FALSE)
    if(ncol(data) == 1)    {
      vol         <- ifelse(reciprocal, 1/abs(diff(range(data))), abs(diff(range(data))))
      loc         <- ifelse(reciprocal, 1/(2 * vol), vol/2)
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
    attr(vol,  "Inverse") <- reciprocal
    noise         <- list(vol=vol, loc=loc)
    attr(noise, "Method") <- method
    class(noise)  <- "NoiseVol"
      return(noise)
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
  
  .listof_exp     <- function(x, ...) {
    nn            <- names(x)
    ll            <- length(x)
    if(length(nn) != ll)  {
      nn          <- paste("Component", seq.int(ll))
    }
    for(i in seq_len(ll)) {
      cat(nn[i], ":\n\n")
      .listof_exp2(x[[i]], ...)
      cat("\n")
    }
  }
  
  .listof_exp2    <- function(x, ...) {
    nn            <- names(x)
    ll            <- length(x)
    if(length(nn) != ll)  {
      nn          <- paste("Component", seq.int(ll))
    }
    for(i in seq_len(ll)) {
      cat(nn[i], ":\n")
      .summ_exp(x[[i]],    ...)
      cat("\n")
    }
  }

  .mat_byrow      <- function(x, nrow, ncol) {
      matrix(x, nrow=nrow, ncol=ncol, byrow=any(dim(as.matrix(x)) == 1))
  }

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
  
  .pick_posidens  <- function(x) {
    x.ind         <- which(!is.na(x), arr.ind=TRUE)
    x.ind[,1L]    <- gsub(".*= ", "", rownames(x)[x.ind[,1L]])
    x.ind[,2L]    <- colnames(x)[as.numeric(x.ind[,2L])]
      stats::setNames(x[!is.na(x)], vapply(seq_len(nrow(x.ind)), function(p, b=x.ind[p,]) paste0(b[2L], ",", b[1L]), character(1L)))
  }
  
  #' @importFrom matrixStats "rowSums2"
  .renorm_z       <- function(z) z/rowSums2(z)
  
  #' @importFrom matrixStats "colRanges" "rowDiffs"
  .SLDC           <- function(x) sum(log(abs(rowDiffs(colRanges(x)))))

  .sq_mat         <- function(x) diag(sqrt(diag(x)))
  
  .summ_exp       <- function(x, digits = max(3L, getOption("digits") - 3L), 
                              signif.stars = getOption("show.signif.stars"),
                              symbolic.cor = x$symbolic.cor, ...) {
    resid         <- x$residuals
    df            <- x$df
    rdf           <- df[2L]
    cat("\n", if(!is.null(x$weights)   && diff(range(x$weights))) "Weighted ", "Residuals:\n", sep = "")
    if(rdf         > 5L) {
      nam         <- c("Min", "1Q", "Median", "3Q", "Max")
      rq          <- if(length(dim(resid)) == 2L) { 
        structure(apply(t(resid), 1L, stats::quantile), dimnames = list(nam, dimnames(resid)[[2L]]))
      } else       {
        zz        <- zapsmall(stats::quantile(resid), digits + 1L)
        structure(zz, names = nam)
      }
      print(rq,    digits = digits, ...)
    } else if(rdf  > 0L)  {
      print(resid, digits = digits, ...)
    } else         {
      cat("ALL", df[1L], "residuals are 0: no residual degrees of freedom!")
      cat("\n")
    }
    if(length(x$aliased) == 0L)         {
      cat("\nNo Coefficients\n")
    } else         {
      if(singular <- df[3L] - df[1L])   {
        cat("\nCoefficients: (", singular, " not defined because of singularities)\n", sep = "")
      } else cat("\nCoefficients:\n")
      coefs       <- x$coefficients
      if(any(aliased     <- x$aliased)) {
        cn        <- names(aliased)
        coefs     <- matrix(NA, length(aliased), 4L, dimnames = list(cn, colnames(coefs)))
        coefs[!aliased,] <- x$coefficients
      }
      stats::printCoefmat(coefs, digits = digits, signif.stars = signif.stars, na.print = "NA", ...)
    }
    cat("\nResidual standard error:", format(signif(x$sigma, digits)), "on", rdf, "degrees of freedom")
    cat("\n")
    if(nzchar(mess       <- stats::naprint(x$na.action))) cat("  (", mess, ")\n", sep = "")
    if(!is.null(x$fstatistic))          {
      cat("Multiple R-squared: ",    formatC(x$r.squared, digits = digits))
      cat(",\tAdjusted R-squared: ", formatC(x$adj.r.squared, digits = digits), 
          "\nF-statistic:",          formatC(x$fstatistic[1L], digits = digits), "on", 
          x$fstatistic[2L], "and", x$fstatistic[3L], "DF,  p-value:", 
          format.pval(stats::pf(x$fstatistic[1L], x$fstatistic[2L], x$fstatistic[3L], lower.tail = FALSE), digits = digits))
      cat("\n")
    }
    correl        <- x$correlation
    if(!is.null(correl))  {
      p           <- NCOL(correl)
      if(p > 1L)   {
        cat("\nCorrelation of Coefficients:\n")
        if(is.logical(symbolic.cor)    && symbolic.cor) {
          print(stats::symnum(correl, abbr.colnames = NULL))
        } else     {
          correl  <- format(round(correl, 2), nsmall = 2, digits = digits)
          correl[!lower.tri(correl)]   <- ""
          print(correl[-1L, -p, drop = FALSE], quote = FALSE)
        } 
      }
    }
    cat("\n")
      invisible(x)
  }

  .tau_noise      <- function(tau, z0)  {
    t0            <- ifelse(length(z0) == 1, z0, mean(z0))
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
  
  .unique_list    <- function(x)  {
    x             <- lapply(x,  function(x) { attributes(x) <- NULL; x} )
      sum(duplicated.default(x, nmax=1L))  == (length(x)     - 1L)
  }
  
  .version_above  <- function(pkg, version) {
    pkg           <- as.character(utils::packageVersion(pkg))
      identical(pkg, version) || (utils::compareVersion(pkg, version) >= 0)
  }
  
  .vol_ellipsoid  <- function(x)  {
   lDet           <- as.numeric(determinant(x$cov, logarithm=TRUE)$modulus/2)
   ld2pi          <- log(base::pi * x$d2)
   exp(ifelse((p2 <- length(x$loc)/2L) > 1,
              p2   * ld2pi + lDet - lgamma(p2 + 1L),
              lDet + ld2pi))
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
    if(isTRUE(attr(x, "Posdens")))                warning("Solution contains positive log-densities and may be spurious\n", call.=FALSE)
    gating        <- attr(x$gating, "Formula")
    expert        <- attr(x$expert, "Formula")
    gate.x        <- !attr(x, "Gating")
    exp.x         <- !attr(x, "Expert")
    net.x         <- !c(gate.x, exp.x)
    crit          <- round(unname(c(x$bic, x$icl, x$aic)), digits)
    hypvol        <- x$hypvol
    noise         <- !is.na(hypvol)
    equalP        <- G <= 1  || attr(x$gating,   "EqualPro")
    equalN        <- noise   && attr(x$gating, "EqualNoise") && equalP
    cat(paste0("\nBest Model", ifelse(length(x$BIC) > 1, paste0(" (according to ", toupper(attr(x, "Criterion")), "):\n"), ":\n"), 
               ifelse(G == 0, "single noise component",  paste0(mclustModelNames(name)$type, " (", name, "), with ",
               G, " component", ifelse(G   > 1, "s",   ""))), ifelse(G == 0 || !noise,   "\n\n", " (and a noise component)\n\n"),
               ifelse(!equalP |
                      G <= 1,  "",   paste0("Equal Mixing Proportions", ifelse(equalN | G <= 1 | !noise, "\n", " (with estimated noise component mixing proportion)\n"))),
               ifelse(!noise,  "",   paste0("Hypervolume of Noise Component: ", round(hypvol, digits), "\n")),
               "BIC = ", crit[1L], " | ICL = ", crit[2L], " | AIC = ",  crit[3L],
               ifelse(any(net.x),    paste0("\nIncluding",    ifelse(all(net.x), " gating and expert", ifelse(!gate.x, " gating", ifelse(!exp.x, " expert", ""))), " network covariates:\n"), "\nNo covariates\n"),
               ifelse(gate.x,  "",   paste0("\tGating: ",  gating, ifelse(exp.x,  "", "\n"))),
               ifelse(exp.x,   "",   paste0("\tExpert: ",  expert, "")), "\n"))
      invisible()
  }

#' @method summary MoEClust
#' @rdname MoE_clust
#' @usage
#' \method{summary}{MoEClust}(object,
#'         classification = TRUE,
#'         parameters = FALSE,
#'         networks = FALSE,
#'         ...)
#' @export
  summary.MoEClust        <- function(object, classification = TRUE, parameters = FALSE, networks = FALSE, ...) {
    object        <- if(inherits(object, "MoECompare")) object$optimal else object
    if(length(classification)  > 1  ||
       !is.logical(classification))               stop("'classification' must be a single logical indicator", call.=FALSE)
    if(length(parameters)  > 1      ||
       !is.logical(parameters))                   stop("'parameters' must be a single logical indicator",     call.=FALSE)
    if(length(networks)    > 1      ||
       !is.logical(networks))                     stop("'networks' must be a single logical indicator",       call.=FALSE)
    G             <- object$G
    attr(G, "range")      <- eval(object$call$G)
    params        <- object$parameters
    hypvol        <- object$hypvol
    equalPro      <- G <= 1 || attr(object$gating, "EqualPro")
    equalN        <- !is.na(hypvol) && attr(object$gating, "EqualNoise") && equalPro
    summ          <- list(data = deparse(object$call$data), n = object$n, d = object$d, G = G, modelName = object$modelName, algo=attr(object, "Algo"),
                          loglik = object$loglik[length(object$loglik)], df = object$df, iters = object$iters, gating = object$gating, expert = object$expert, 
                          bic=unname(object$bic), icl = unname(object$icl), aic = unname(object$aic), pro = params$pro, mean = params$mean, variance = params$variance$sigma, 
                          Vinv = params$Vinv, hypvol = hypvol, z = object$z, equalPro = equalPro, equalNoise = equalN, classification = object$classification, noise.gate = attr(object, "NoiseGate"),
                          expert = object$expert, gating = object$gating, printClass = classification, printParams = parameters, printNetwork = networks)
    class(summ)   <- "summary_MoEClust"
    attr(summ, "Posdens") <- attr(object, "Posdens")
      summ
 }

#' @method print summary_MoEClust
#' @importFrom mclust "mclustModelNames"
#' @export
  print.summary_MoEClust  <- function(x, digits = 3L, ...) {
    if(length(digits) > 1 || !is.numeric(digits) ||
       digits     <= 0)                           stop("Invalid 'digits'", call.=FALSE)
    tmp           <- data.frame(log.likelihood = round(x$loglik, digits), n = x$n, d = x$d, df = x$df, iters = x$iters,
                                BIC = round(x$bic, digits), ICL = round(x$icl, digits), AIC = round(x$aic, digits))
    tmp           <- if(is.na(x$hypvol))   tmp    else cbind(tmp, HypVol = x$hypvol)
    tmp           <- cbind(tmp, Algo = x$algo)
    rownames(tmp) <- NULL
    name          <- x$modelName
    G             <- x$G
    range.G       <- attr(G, "range")
    if(!is.null(range.G)  && length(range.G) > 1  &&
       G          == min(range.G))                message("Best model occurs at the min of the number of components considered\n")
    if(!is.null(range.G)  && length(range.G) > 1  &&
       G          == max(range.G))                message("Best model occurs at the max of the number of components considered\n")
    if(isTRUE(attr(x, "Posdens")))                warning("Solution contains positive log-densities and may be spurious\n", call.=FALSE)
    noise         <- !is.na(x$hypvol)
    gating        <- attr(x$gating, "Formula")
    expert        <- attr(x$expert, "Formula")
    gate.x        <- gating == "~1"
    exp.x         <- expert == "~1"
    equalP        <- x$equalPro && gate.x
    equalN        <- noise  && x$equalNoise && equalP
    title         <- "Gaussian Parsimonious Clustering Model with Covariates"
    cat(paste0("------------------------------------------------------\n", title, "\nData: ",
               x$data,"\n", "------------------------------------------------------\n\n",
               "MoEClust: ", ifelse(G == 0, "single noise component", 
                                   paste0(name, " (", mclustModelNames(name)$type, "), with ",
                                   G, " component", ifelse(G > 1, "s", ""))),
               ifelse(G == 0  || is.na(x$hypvol), "\n", " (and a noise component)\n"),
               paste0("\nGating Network Covariates:  ", ifelse(gate.x, "None", gating)),
               paste0("\nExpert Network Covariates:  ", ifelse(exp.x,  "None", expert)),
               ifelse(G  > 1  && gate.x,                paste0("\nEqual Mixing Proportions:   ",  equalP), ""),
               paste0("\nNoise Component:            ", noise, ""),
               ifelse(noise,                            paste0("\nNoise Component Estimation: ", attr(x$hypvol, "Meth")), ""),
               ifelse(G  > 1  && !gate.x && noise,      paste0("\nNoise Component Gating:     ", x$noise.gate), ""),
               ifelse(G  > 1  && noise   && equalP,     paste0("\nNoise Proportion Estimated: ", !equalN, "\n\n"), "\n\n")))
    print(tmp, row.names = FALSE)
    if(isTRUE(x$printClass))   {
      cat("\nClustering table :")
      print(table(x$classification), row.names = FALSE)
    }
    if(isTRUE(x$printParams))  {
      params     <- list("Mixing proportions"  = x$pro,
                         "Component means"     = x$mean,
                         "Component variances" = x$variance)
      if(!is.na(x$hypvol))     {
        attributes(x$hypvol)  <- NULL
        params   <- c(params, list("Hypervolume of noise component" = x$hypvol))
      }
      class(params)     <- "listof"
      cat("\n")
      print(params)
    }
    if(isTRUE(x$printNetwork) &&
       !all(gate.x, exp.x))    {
      if(isFALSE(gate.x))      {
        gating   <- list("Gating Network"      = x$gating)
        class(gating)   <- "listof"
        print(gating, call = FALSE)
        cat("\n")  
      } else                                      message("No gating network to display\n")
      if(isFALSE(exp.x))       {
        expert   <- list("Expert Network"      = x$expert)
        class(expert)   <- "listof"
        print(expert, call = FALSE)
        cat("\n") 
      } else                                      message("No expert network to display\n")
    }
    if(isTRUE(x$printParams)  && 
       isFALSE(exp.x)) {
      message("\n\n\nUsers are cautioned against interpreting the component mean parameters in the presence of expert network covariates.\nThese are in fact the posterior means of the fitted values of the expert network.\nThe observation-specific component means (i.e. the fitted values themselves) should be consulted instead.\nThese can obtained via predict(object)$mean.\n")
    }
    cat("\n")
      invisible()
  }

#' @method print MoECompare
#' @rdname MoE_compare
#' @usage
#' \method{print}{MoECompare}(x,
#'       index = seq_len(x$pick),
#'       posidens = TRUE,
#'       rerank = FALSE,
#'       digits = 3L,
#'       details = TRUE, 
#'       maxi = length(index),
#'       ...)
#' @export
  print.MoECompare       <- function(x, index=seq_len(x$pick), posidens = TRUE, rerank = FALSE, digits = 3L, details = TRUE, maxi = length(index), ...) {
    index                <- if(is.logical(index)) which(index) else index
    if(length(index) < 1 || (!is.numeric(index) &&
       (any(index    < 1  | index > x$pick))))    stop("Invalid 'index'",  call.=FALSE)
    if(length(digits)     > 1    ||
       !is.numeric(digits)       ||
       digits            <= 0)                    stop("Invalid 'digits'", call.=FALSE)
    if(length(posidens)   > 1    ||
       !is.logical(posidens))                     stop("'posidens' must be a single logical indicator", call.=FALSE)
    if(length(rerank)     > 1    ||
       !is.logical(rerank))                       stop("'rerank' must be a single logical indicator",   call.=FALSE)
    if(length(details)    > 1    ||
       !is.logical(details))                      stop("'details' must be a single logical indicator",  call.=FALSE)
    if(length(maxi)      != 1    ||
       !is.numeric(maxi) ||
       maxi        <= 0  ||
       floor(maxi) != maxi)                       stop("'maxi' must be a single strictly positive integer",  call.=FALSE)
    maxi                 <- min(maxi, length(index))
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
    title                <- "Comparison of Gaussian Parsimonious Clustering Models with Covariates"
    if(isTRUE(details))   {
      cat(paste0("---------------------------------------------------------------------\n", 
                 title, "\nData: ", x$data, "\nRanking Criterion: ", toupper(attr(x, "Crit")), "\nOptimal Only: ", attr(x, "Opt"), 
                "\n---------------------------------------------------------------------\n\n"))
    }
    comp.res             <- data.frame(do.call(cbind, x[-seq_len(3L)]))
    if(isTRUE(posidens))  {
      comp.res           <- comp.res[index,, drop=FALSE]
      comp.res$posidens  <- if(!any(comp.res$posidens == "TRUE")) NULL else comp.res$posidens 
    } else                {
      index              <- index[!x$posidens[index]]
      comp.res           <- comp.res[index,, drop=FALSE]
      comp.res$posidens  <- NULL
    }
    comp.res             <- comp.res[,c(rep(TRUE, 2L), !vapply(comp.res[-seq_len(2L)], function(x) all(x == ""), logical(1L))), drop=FALSE]
    comp.res             <- if(isTRUE(details)) cbind(rank = if(isTRUE(rerank)) seq_along(index) else index, comp.res) else comp.res[,-which(colnames(comp.res) == "MoENames"), drop=FALSE]
    rownames(comp.res)   <- NULL
    print(comp.res[seq_len(maxi),], row.names = FALSE)
    cat("\n")
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
  print.MoE_gating       <- function(x, call = FALSE, ...) {
    equalpro      <- attr(x, "EqualPro")
    formula       <- attr(x, "Formula")
    noise         <- attr(x, "Noise")
    equalNoise    <- noise   && equalpro
    gateNoise     <- noise   && !equalpro && formula != "~1"
    if(ifelse(inherits(x, "multinom"),
              x$convergence  == 1,
              isTRUE(x$converged)))               warning("Multinomial logistic regression failed to converge", call.=FALSE, immediate.=TRUE)
    class(x)      <- class(x)[class(x)    != "MoE_gating"]
    if(isTRUE(call)          &&
      !is.null(cl <- x$call)) {
      cat("Call:\n")
      dput(cl, control = NULL)
    }
    cat("\nCoefficients:\n")
    print(stats::coef(x), ...)
    cat(paste("\nFormula:", formula, "\n"))
    cat(paste("Noise:",     noise,   "\n"))
    if(gateNoise)                                 cat(paste("Noise Component Gating:", attr(x, "NoiseGate"), "\n"))
    cat(paste("EqualPro:",  equalpro, ifelse(equalNoise, "\n", "")))
    if(equalNoise)                                cat(paste("Noise Proportion Estimated:", !attr(x, "EqualNoise")))
    if(equalpro)                                  message("\n\nCoefficients set to zero as this is an equal mixing proportion model")
      invisible()
  }

#' @method print MoE_expert
#' @export
  print.MoE_expert       <- function(x, call = FALSE, ...) {
   if(all(is.na(x) | (names(x) == "Cluster0")))   stop("No expert network exists for models with only a noise component", call.=FALSE)
   formula        <- attr(x, "Formula")
   if(inherits(x, "summary_MoEexp"))      {
     attributes(x)[-1L]  <- NULL
     class(x)     <- "listof"
     if(isTRUE(call))     {
       print(x, ...)
     } else        {
       .listof_exp(x, ...)
     }
   }   else        {
     attributes(x)[-1L]  <- NULL
     class(x)     <- "listof"
     cat("\n")
     for(g in seq_along(x))               {
       cat(paste0("Cluster", g, " :\n\n"))
       if(isTRUE(call)   &&
          !is.null(cl    <- x[[g]]$call)) {
         cat("Call:\n")
         dput(cl, control = NULL)
         cat("\n")
       }
       cat("Coefficients:\n")
       print(stats::coef(x[[g]]), ...)
       cat("\n")
     }
     cat(paste("Formula:", formula))
   }
     invisible()
  }

#' @method summary MoE_gating
#' @export
  summary.MoE_gating     <- function(object, ...) {
    equalnoise    <- attr(object, "EqualNoise")
    formula       <- attr(object, "Formula")
    equalpro      <- attr(object, "EqualPro")
    noise         <- attr(object, "Noise")
    noise.gate    <- attr(object, "NoiseGate")
    class(object) <- class(object)[2L]
    summ          <- summary(object, ...)
    summ$OddsRatios           <- exp(summ$coefficients)
    class(summ)   <- "summary_MoEgate"
    attr(summ, "Class")       <- class(object)
    attr(summ, "EqualNoise")  <- equalnoise
    attr(summ, "Formula")     <- formula
    attr(summ, "EqualPro")    <- equalpro
    attr(summ, "Noise")       <- noise
    attr(summ, "NoiseGate")   <- noise.gate
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
    attr(summ, "Formula")     <- attr(object, "Formula") 
      summ
  }

#' @method print summary_MoEgate
#' @export
  print.summary_MoEgate  <- function(x, ...) {
    equalpro      <- attr(x, "EqualPro")
    formula       <- attr(x, "Formula")
    noise         <- attr(x, "Noise")
    class(x)      <- "MoE_gating"
    print(x, ...)
    cat("\n\nOddsRatios:\n")
    print(x$OddsRatios)
    cat("\nStd. Errors:\n")
    print(x$standard.errors)
    message("\n\n\nUsers are cautioned against making inferences about statistical significance from summaries of the coefficients in the gating network\n")
    class(x)      <- "summary_MoEgate"
      invisible(x)
  }

#' @method print summary_MoEexp
#' @export
  print.summary_MoEexp   <- function(x, ...) {
   if(all(is.na(x) | (names(x) == "Cluster0")))   stop("No expert network exists for models with only a noise component", call.=FALSE)
   class(x)       <- c("MoE_expert", class(x))
   print(x, ...)
   cat(paste("Formula:", attr(x, "Formula"), "\n"))
   message("\n\n\nUsers are cautioned against making inferences about statistical significance from summaries of the coefficients in the expert network\n")
   class(x)       <- "summary_MoEexp"
     invisible(x)
  }
