#' Mixtures of Experts: Model-Based Clustering with Covariates
#'
#' Fits Mixture of Experts models with \pkg{mclust}-family covariance structures. In other words, performs model-based clustering via the EM algorithm where covariates are allowed to enter neither, either, or both the mixing proportions (gating network) and/or component densities (expert network).
#' @param data A numeric vector, matrix, or data frame of observations. Categorical variables are not allowed. If a matrix or data frame, rows correspond to observations and columns correspond to variables.
#' @param G An integer vector specifying the numbers of mixture components (clusters) to fit. Defaults to \code{G=1:9}. Must be a strictly positive integer, unless a noise component is included in the estimation, in which case \code{G=0} is allowed (see \code{\link{MoE_control}}).
#' @param modelNames A vector of character strings indicating the models to be fitted in the EM phase of clustering. With \code{n} observations and \code{d} variables, the defaults are:\cr
#' \tabular{ll}{for univariate data \tab \code{c("E", "V")}\cr
#' for multivariate data \eqn{n > d}{n > d} \tab \code{mclust.options("emModelNames")}\cr
#' for high-dimensional multivariate data \eqn{n \leq d}{n <= d} \tab \code{c("EII", "VII", "EEI", "EVI", "VEI", "VVI")}
#' }
#' \cr
#' For single-component models these options reduce to:\cr
#' \tabular{ll}{
#' for univariate data \tab \code{"E"}\cr
#' for multivariate data \eqn{n > d}{n > d} \tab \code{c("EII", "EEI", "EEE")}\cr
#' for high-dimensional multivariate data \eqn{n \leq d}{n <= d}  \tab \code{c("EII", "EEI")}
#' }
#' For zero-component models with a noise component only the \code{"E"} and \code{"EII"} models will be fit for univariate and multivariate data, respectively. The help file for \code{\link[mclust]{mclustModelNames}} further describes the available models (though the \code{"X"} in the single-component models will be coerced to \code{"E"} if supplied that way).
#' @param gating A formula for determining the model matrix for the multinomial logistic regression in the gating network when covariates enter the mixing proportions. This will be ignored where \code{G=1}. Interactions etc. are permitted. The specification of the LHS of the formula is ignored.
#' @param expert A formula for determining the model matrix for the (multivariate) WLS in the expert network when covariates are included in the component densities. Interactions etc. are permitted. The specification of the LHS of the formula is ignored.
#' @param network.data An optional data frame in which to look for the covariates in the \code{gating} &/or \code{expert} network formulas, if any. If not found in \code{network.data}, any supplied \code{gating} &/or \code{expert} covariates are taken from the environment from which \code{MoE_clust} is called.
#' @param control A list of control parameters for the EM and other aspects of the algorithm. The defaults are set by a call to \code{\link{MoE_control}}.
#' @param ... An alternative means of passing control parameters directly via the named arguments of \code{\link{MoE_control}}. Do not pass the output from a call to \code{\link{MoE_control}} here! This argument is only relevant for the \code{\link{MoE_clust}} function and will be ignored for the associated \code{print} and \code{summary} functions.
#' @param x,object,digits Arguments required for the \code{print} and \code{summary} functions: \code{x} and \code{object} are objects of class \code{"MoEClust"} resulting from a call to \code{\link{MoE_clust}}, while \code{digits} gives the number of decimal places to round to for printing purposes (defaults to 2).

#' @importFrom matrixStats "rowLogSumExps"
#' @importFrom mclust "emControl" "hc" "hclass" "hcE" "hcEEE" "hcEII" "hcV" "hcVII" "hcVVV" "Mclust" "mclust.options" "mclustBIC" "mclustModelNames" "mclustVariance" "mstep" "mstepE" "mstepEEE" "mstepEEI" "mstepEEV" "mstepEII" "mstepEVE" "mstepEVI" "mstepEVV" "mstepV" "mstepVEE" "mstepVEI" "mstepVEV" "mstepVII" "mstepVVE" "mstepVVI" "mstepVVV" "nVarParams" "unmap"
#' @importFrom mvnfast "dmvn"
#' @importFrom nnet "multinom"
#' @return A list (of class \code{"MoEClust"}) with the following named entries, mostly corresponding to the chosen optimal model (as determined by the \code{criterion} within \code{\link{MoE_control}}):\cr
#' \describe{
#' \item{\code{call}}{The matched call.}
#' \item{\code{data}}{The input data, as a \code{data.frame}.}
#' \item{\code{modelName}}{A character string denoting the \pkg{mclust} model type at which the optimal \code{criterion} occurs.}
#' \item{\code{n}}{The number of observations in the \code{data}.}
#' \item{\code{d}}{The dimension of the \code{data}.}
#' \item{\code{G}}{The optimal number of mixture components.}
#' \item{\code{BIC}}{A matrix of \emph{all} BIC values with \code{length{G}} rows and \code{length(modelNames)} columns. May include missing entries: \code{NA} represents models which were not visited, \code{-Inf} represents models which were terminated due to error, for which a log-likelihood could not be estimated. Inherits the classes \code{"MoECriterion"} and \code{"mclustBIC"}, for which a dedicated plotting function exists.}
#' \item{\code{ICL}}{A matrix of \emph{all} ICL values with \code{length{G}} rows and \code{length(modelNames)} columns. May include missing entries: \code{NA} represents models which were not visited, \code{-Inf} represents models which were terminated due to error, for which a log-likelihood could not be estimated. Inherits the classes \code{"MoECriterion"} and \code{"mclustBIC"}, for which a dedicated plotting function exists.}
#' \item{\code{AIC}}{A matrix of \emph{all} AIC values with \code{length{G}} rows and \code{length(modelNames)} columns. May include missing entries: \code{NA} represents models which were not visited, \code{-Inf} represents models which were terminated due to error, for which a log-likelihood could not be estimated. Inherits the classes \code{"MoECriterion"} and \code{"mclustBIC"}, for which a dedicated plotting function exists.}
#' \item{\code{bic}}{The BIC value corresponding to the optimal model. May not necessarily be the optimal BIC.}
#' \item{\code{icl}}{The ICL value corresponding to the optimal model. May not necessarily be the optimal ICL.}
#' \item{\code{aic}}{The AIC value corresponding to the optimal model. May not necessarily be the optimal AIC.}
#' \item{\code{gating}}{An object of class \code{"MoE_gating"} and either \code{"multinom"} or \code{"glm"} (for single-component models) giving the \code{\link[nnet]{multinom}} regression coefficients of the \code{gating} network. If \code{gating} covariates were \emph{NOT} supplied (or the best model has just one component), this corresponds to a RHS of ~1, otherwise the supplied \code{gating} formula. As such, a fitted \code{gating} network is always returned even in the absence of supplied covariates. The number of parameters to penalise by for \code{\link{MoE_crit}} is given by \code{length(coef(gating))}, and the \code{gating} formula used is stored here as an attribute. If there is a noise component, its coefficients are those for the \emph{last} component. \strong{Users are cautioned against making inferences about statistical significance from summaries of the coefficients in the gating network}.}
#' \item{\code{expert}}{An object of class \code{"MoE_expert"} and \code{"lm"} giving the (multivariate) WLS regression coefficients of the \code{expert} network. If \code{expert} covariates were NOT supplied, this corresponds to a RHS of ~1, otherwise the supplied \code{expert} formula. As such, a fitted \code{expert} network is always returned even in the absence of supplied covariates. The number of parameters to penalise by for \code{\link{MoE_crit}} is given by \code{G * length(coef(expert[[1]]))}, and the \code{expert} formula used is stored here is an attribute. \strong{Users are cautioned against making inferences about statistical significance from summaries of the coefficients in the expert network}.}
#' \item{\code{loglik}}{The vector of increasing log-likelihood values for every EM iteration under the optimal model.}
#' \item{\code{df}}{The number of estimated parameters in the optimal model (i.e. the number of 'used' degrees of freedom). Subtract this number from \code{n} to get the degrees of freedom.}
#' \item{\code{hypvol}}{The hypervolume parameter for the noise component if required, otherwise set to \code{NA} (see \code{\link{MoE_control}}).}
#' \item{\code{parameters}}{A list with the following components:\cr
#' \itemize{
#' \item{\code{pro} - }{The mixing proportions: either a vector of length \code{G} or, if \code{gating} covariates were supplied, a matrix with an entry for each observation (rows) and component (columns).}
#' \item{\code{mean} - }{The means of each component. If there is more than one component, this is a matrix whose \emph{k}-th column is the mean of the \emph{k}-th component of the mixture model.}
#' \item{\code{variance} - }{A list of variance parameters of each component of the model. The components of this list depend on the model type specification. See the help file for \code{\link[mclust]{mclustVariance}} for details.}
#' \item{\code{Vinv} - }{The inverse of the hypervolume parameter for the noise component if required, otherwise set to \code{NULL} (see \code{\link{MoE_control}}).}
#' }}
#' \item{\code{z}}{The final responsibility matrix whose \code{[i,k]}-th entry is the probability that observation \emph{i} belonds to the \emph{k}-th component. If there is a noise component, its values are found in the \emph{last} column.}
#' \item{\code{classification}}{The vector of cluster labels for the chosen model corresponding to \code{z}, i.e. \code{max.col(z)}. Observations belonging the noise component will belong to component \code{0}.}
#' \item{\code{uncertainty}}{The uncertainty associated with the \code{classification}.}
#' \item{\code{net.covs}}{A data frame gathering the unique set of covariates used in the \code{gating} and \code{expert} networks, if any. Will be missing in the absence of gating or expert network covariates, and supplied gating covariates will be exluded if the optimal model has only one component.}
#' \item{\code{resid.data}}{In the presence of expert network covariates, this is the augmented data (as a data frame) actually used in the clustering at convergence consisting of the \code{(n * G) * d} matrix of (multivariate) WLS residuals. Will be missing in the absence of expert network covariates.}
#' \item{\code{DF}}{A matrix of giving numbers of estimated parameters (i.e. the number of 'used' degrees of freedom) for \emph{all} visited models, with \code{length{G}} rows and \code{length(modelNames)} columns. Subtract these numbers from \code{n} to get the degrees of freedom. May include missing entries: \code{NA} represents models which were not visited, \code{-Inf} represents models which were terminated due to error, for which parameters could not be estimated. Inherits the classes \code{"MoECriterion"} and \code{"mclustBIC"}, for which a dedicated plotting function exists.}
#' \item{\code{iters}}{A matrix giving the total number of EM iterations with \code{length{G}} rows and \code{length(modelNames)} columns. May include missing entries: \code{NA} represents models which were not visited, \code{Inf} represents models which were terminated due to singularity/error and thus would never have converged.}
#' }
#' Dedicated \code{\link[=plot.MoEClust]{plot}}, \code{print} and \code{summary} functions exist for objects of class \code{"MoEClust"}. The results can be coerced to the \code{"Mclust"} class to access other functions from the \pkg{mclust} package via \code{\link{as.Mclust}}.
#' @details The function effectively allows 4 different types of Mixture of Experts model (as well as the different models in the mclust family, for each): i) the standard finite Gaussian mixture, ii) covariates only in the gating network, iii) covariates only in the expert network, iv) the full Mixture of Experts model with covariates entering both the mixing proportions and component densities. Note that having the same covariates in both networks is allowed.\cr
#'
#' While model selection in terms of choosing the optimal number of components and the \pkg{mclust} model type is performed within \code{\link{MoE_clust}}, using one of the \code{criterion} options within \code{\link{MoE_control}}, choosing between multiple fits with different combinations of covariates or different initialisation settings can be done by supplying objects of class \code{"MoEClust"} to \code{\link{MoE_compare}}.
#' @note The EM algorithm finishes on an extra M-step once the best model has been identified.
#'
#' Where \code{BIC}, \code{ICL}, \code{AIC}, \code{DF} and \code{iters} contain \code{NA} entries, this corresponds to a model which was not run; for instance a VVV model is never run for single-component models as it is equivalent to EEE. As such, one can consider the value as not really missing, but equivalent to the EEE value. \code{BIC}, \code{ICL}, \code{AIC} and \code{DF} all inherit the class \code{"MoECriterion"}, for which a dedicated print function exists.
#' @seealso \code{\link{MoE_compare}}, \code{\link{plot.MoEClust}}, \code{\link{MoE_control}}, \code{\link{as.Mclust}}, \code{\link{MoE_crit}}, \code{\link{MoE_estep}}, \code{\link{MoE_dens}}, \code{\link[mclust]{mclustModelNames}}, \code{\link[mclust]{mclustVariance}}
#' @export
#' @references K. Murphy and T. B. Murphy (2017). Parsimonious Model-Based Clustering with Gating and Expert Network Covariates.
#'
#' C. Fraley and A. E. Raftery (2002). Model-based clustering, discriminant analysis, and density estimation. \emph{Journal of the American Statistical Association}, 97:611-631.
#' @author Keefe Murphy - \href{keefe.murphy@ucd.ie}{<keefe.murphy@ucd.ie>}
#'
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
#' # Extract the model with highest ICL
#' (comp <- MoE_compare(m1, m2, m3, m4, criterion="icl"))
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
#' splom(~ hema | z, groups=sex)
#'
#' # Convert results to the "Mclust" class for further visualisation
#' plot(as.Mclust(best), what="density")}
  MoE_clust       <- function(data, G = 1:9, modelNames = NULL, gating = NULL, expert = NULL, network.data = NULL, control = MoE_control(...), ...) {

  # Definitions and storage set-up
    call          <- match.call()
    multi         <- missing(modelNames)
    gate.x        <- !missing(gating)
    exp.x         <- !missing(expert)
    criterion     <- control$criterion
    equal.pro     <- control$equalPro
    init.z        <- control$init.z
    exp.init      <- control$exp.init
    do.joint      <- exp.init$joint
    do.mahala     <- exp.init$mahalanobis
    max.init      <- exp.init$max.init
    max.it        <- control$itmax[1]
    stopx         <- control$stopping
    tol           <- control$tol[1]
    warnit        <- control$warn.it
    noise         <- control$noise.init
    noise.meth    <- control$noise.meth
    hcName        <- control$hc.meth
    miss.init     <- control$miss.init
    miss.hc       <- control$miss.hc
    verbose       <- control$verbose
    itwarn        <- warnit > 2
    control       <- control[which(names(control) %in% c("eps", "tol", "itmax", "equalPro"))]
    if(!multi     &&
       !all(is.character(modelNames)))            stop("'modelNames' must be a vector of character strings")
    if(gate.x     &&
       !inherits(gating, "formula"))              stop("'gating' must be a formula")
    if(exp.x      &&
       !inherits(expert, "formula"))              stop("'expert' must be a formula")
    if(missing(data))                             stop("'data' must be supplied!")

    tmp.nam       <- as.character(substitute(data))
    data          <- as.data.frame(data)
    num.X         <- vapply(data, is.numeric, logical(1L))
    if(anyNA(data))    {
      if(verbose)                                 message("Rows with missing values removed from data")
      data        <- data[stats::complete.cases(data),, drop=FALSE]
    }
    if(sum(num.X) != ncol(data))    {
      if(verbose)                                 message("Non-numeric columns removed from data")
      data        <- data[,num.X, drop=FALSE]
    }
    X             <- as.matrix(data)
    n             <- nrow(X)
    d             <- ncol(X)
    noise.null    <- is.null(noise)
    if(missing(G) && !noise.null) G   <- 0:9
    if(any(G      != floor(G))        &&
       any(G < ifelse(noise.null, 1, 0)))         stop(paste0("'G' must be ", ifelse(noise.null, "strictly positive", "strictly non-negative when modelling with a noise-component")))
    if(any(G > n))   G        <- G[G  <= n]

    mod.fam       <- mclust.options("emModelNames")
    range.G       <- sort(as.integer(unique(G)))
    if(!noise.null)   {
     if(length(noise)         != n)               stop("'noise.init' must be a vector of length n")
     if(!is.logical(noise))    {
       if(any(match(noise, seq_len(n),
              nomatch = 0)    == 0))              stop("Numeric noise must correspond to row indices of data")
       noise      <- as.logical(match(seq_len(n), noise, nomatch = 0))
     }
     Vinv         <- .noise_vol(X, noise.meth)
     nnoise       <- sum(as.numeric(noise))
     noisen       <- n - nnoise
     if(any(G      > noisen)) range.G <- range.G[range.G <= noisen]
    } else   {
      noise       <- rep(FALSE, n)
      nnoise      <- 0
      noisen      <- n
      Vinv        <- NULL
    }
    Gall          <- ifelse(noise.null, all(G > 1), all(G[G > 0] > 1))
    Gany          <- ifelse(noise.null, any(G > 1), any(G[G > 0] > 1))
    allg0         <- all(G == 0)
    anyg0         <- any(G == 0)
    if(allg0)      {
      Vinv        <- ifelse(noise.null, .noise_vol(X, noise.meth), Vinv)
      noise.null  <- FALSE
      noise       <- rep(TRUE, n)
      nnoise      <- n
      noisen      <- 0
    }
    if((uni <- d  == 1))       {
      mfg         <- c("E", "V")
      mf1         <- mf0      <- "E"
      colnames(X) <- tmp.nam[length(tmp.nam)]
    } else        {
      mf0         <- "EII"
      if((low.dim <- n > d))   {
        mfg       <- mod.fam
        mf1       <- c("EII", "EEI", "EEE")
      } else      {
        mfg       <- mod.fam[1:6]
        mf1       <- c("EII", "EEI")
      }
    }
    low.dim       <- !uni && low.dim
    if(!multi)    {
      mNs         <- toupper(modelNames)
      if(any(sX   <- grepl("X",     mNs)))      {
       mNs        <- gsub("X", "E", mNs)
       if(verbose &&
          all(is.element(mNs,             mfg)))  message(paste0("'modelNames' which contain 'X' coerced to ", paste(shQuote(mNs[sX]), collapse=", ")))
      }
      if(Gany     && any(!is.element(mNs, mfg)))  stop(paste0("Invalid 'modelNames'", ifelse(uni, " for univariate data", ifelse(low.dim, "", " for high-dimensional data")), "!"))
      if(!Gall)   {
        if(any(sZ <- !is.element(mNs,     mf1))){
          mf1     <- tryCatch(unname(vapply(mNs,  function(x)  switch(EXPR=x, E=, V="E", EII=, VII="EII", EEI=, VEI=, EVI=, VVI="EEI", EEE=, EVE=, VEE=,  VVE=, EEV=, VEV=, EVV=, VVV="EEE"), character(1L))),
                              error=function(e) { e$message <- paste0("Invalid 'modelNames' for single component models", ifelse(uni, " for univariate data", ifelse(low.dim, "", " for high-dimensional data")), "!")
                                                  stop(e) } )
          if(verbose)                             message(paste0("'modelNames'", ifelse(any(sX), " further", ""), " coerced from ", paste(shQuote(mNs[sZ]), collapse=", "), " to ", paste(shQuote(mf1[sZ]), collapse=", "), " where 'G'=1"))
        }
      }
      mfg         <- mNs
    }
    mf1           <- unique(mf1)
    mfg           <- unique(mfg)
    all.mod       <- unique(c(if(any(G == 0)) mf0, if(any(G == 1)) mf1, if(any(G > 1)) mfg))
    multi         <- length(all.mod)   > 1
    BICs          <- ICLs     <-
    AICs          <- DF.x     <- it.x <- provideDimnames(matrix(NA, nrow=length(range.G), ncol=length(all.mod)), base=list(as.character(range.G), all.mod))
    crit.tx       <- crit.gx  <- -sqrt(.Machine$double.xmax)

  # Define the gating formula
    if(allg0 && gate.x)        {  if(verbose)     message("Can't include gating network covariates in a noise-only model")
      gate.x      <- FALSE
    }
    gate.G        <- ifelse((range.G   + !noise.null) > 1, gate.x, FALSE)
    if(gate.x)    {
      gating      <- tryCatch(stats::update.formula(stats::as.formula(gating), z ~ .),
                              error=function(e)   stop("Invalid 'gating' network formula supplied"))
      environment(gating)     <- environment()
      if(gating[[3]] == 1)     { if(verbose)      message("Not including gating network covariates with only intercept on gating formula RHS")
        gate.x    <- FALSE
        gate.G    <- rep(gate.x, length(range.G))
      }
      if(any(G    <= 1))       {  if(verbose)     message("Can't include gating network covariates in a single component mixture")
        gate.G[G  <= 1]       <- FALSE
      }
    } else gating <- stats::as.formula(z ~ 1)
    if(equal.pro  && gate.x)   { if(verbose)      message("Can't constrain mixing proportions to be equal when gating covariates are supplied")
      equal.pro   <- FALSE
    }
    equal.tau     <- ifelse((range.G + !noise.null) == 1, TRUE, equal.pro)

  # Define the expert formula
    if(allg0 && gate.x)        {  if(verbose)     message("Can't include expert network covariates in a noise-only model")
      exp.x       <- FALSE
    }
    if(exp.x)     {
      expert      <- tryCatch(stats::update.formula(stats::as.formula(expert), X ~ .),
                              error=function(e)   stop("Invalid 'expert' network formula supplied"))
      environment(expert)     <- environment()
      if(expert[[3]]   == 1)   { if(verbose)      message("Not including expert network covariates with only intercept on expert formula RHS")
        exp.x     <- FALSE
      }
      Nseq        <- seq_len(n)
    } else expert <- stats::as.formula(X ~ 1)
    one.G         <- all(G == 1)
    if(init.z == "mclust"     && !one.G &&
       !any(gate.x, exp.x))                       stop("Can't initialise using 'mclust' when there are no gating or expert covariates: try another 'init.z' method")

  # Tell network formulas where to look for variables
    if(!missing(network.data) &&
       !is.data.frame(network.data))              stop("'network.data' must be a data.frame if supplied")
    if(is.null(network.data))  {
      gate.covs   <- if(gate.x) stats::model.frame(gating[-2]) else matrix(0, nrow=n, ncol=0)
      expx.covs   <- if(exp.x)  stats::model.frame(expert[-2]) else matrix(0, nrow=n, ncol=0)
      netdat      <- cbind(gate.covs, expx.covs)
      netdat      <- data.frame(if(ncol(netdat) > 0) netdat[!duplicated(names(netdat))] else netdat, stringsAsFactors=TRUE)
      attr(netdat, "Gating")  <- gate.names    <- if(is.null(names(gate.covs)))      NA else names(gate.covs)
      attr(netdat, "Expert")  <- expx.names    <- if(is.null(names(expx.covs)))      NA else names(expx.covs)
      attr(netdat, "Both")    <- if(length(intersect(gate.names, expx.names)) == 0)  NA else intersect(gate.names, expx.names)
    } else netdat <- network.data
    net.cts       <- !sapply(netdat, is.factor)
    nct           <- sum(net.cts)
    init.var      <- ifelse(do.joint && !allg0, d + nct, d)
    highd         <- init.var >= n
    multv         <- init.var  > 1
    if(!multv)     {
      init.z      <- ifelse(miss.init, "quantile", init.z)
    } else if(init.z   == "quantile" && !one.G)   stop("Quantile-based initialisation of the allocations is only permitted for univariate data without continuous expert network covariates")
    if(init.z     == "hc")       {
      if(miss.hc)  {
        hcName    <- ifelse(highd, "EII", "VVV")
      } else if(init.z == "hc"  && !one.G) {
        if(multv  &&
           is.element(hcName,  c("E",     "V")))  stop("'hc.meth' can only be 'E' or 'V' for univariate data without continuous expert network covariates")
        if(highd  && !is.element(hcName, c("EII",
                                "VII",  "EEE")))  warning("Consider a diagonal 'EII' or 'VII' model (or equal volume 'EEE' model) for 'hc.meth' for initialising allocations for high-dimensional data", call.=FALSE)
        if(!multv && !is.element(hcName, c("VVV",
                                "E",      "V")))  warning("Possibly invalid 'hc.meth' for univariate data", call.=FALSE)
      }
    }

  # Loop over range of G values and initialise allocations
    for(g in range.G) {
      if(isTRUE(verbose))        cat(paste0("\n", g, " cluster model", ifelse(multi, "s", ""), " -\n"))
      x.dat       <- replicate(max(g, 1), X, FALSE)
      h           <- which(range.G == g)
      equal.pro   <- equal.tau[h]
      gate.g      <- gate.G[h]
      exp.g       <- exp.x && g > 0
      exp.init    <- exp.g && do.mahala
      Gseq        <- seq_len(g)
      gN          <- g + !noise.null
      z           <- matrix(0, n, gN)

    # Initialise Expert Network & Allocations
      XI        <- if(exp.g  && do.joint) cbind(X, netdat[,net.cts])[!noise,, drop=FALSE] else X[!noise,, drop=FALSE]
      z.tmp     <- unmap(if(g > 1) switch(init.z, hc=as.vector(hclass(tryCatch(hc(XI, modelName=hcName, minclus=g), error=function(e) stop("Error in hierarchical clustering initialisation")), g)),
                                          kmeans=stats::kmeans(XI, g)$cluster, random=sample(Gseq, noisen, replace=TRUE), quantile=MoE_qclass(XI, g),
                                          mclust=Mclust(XI, g, verbose=FALSE, control=emControl(equalPro=equal.pro))$classification) else rep(1, ifelse(noisen == 0, n, noisen)))
      if(exp.g)    {
        z.mat     <- z.alloc  <- matrix(0, nrow=n * g, ncol=g)
        muX       <- if(uni)  rep(0, g) else matrix(0, nrow=d, ncol=g)
      } else {
        exp.pen   <- g * d
      }
      if(exp.init) {
        tmp.z     <- matrix(NA, nrow=ifelse(noisen == 0, n, noisen), ncol=g)
        mahala    <- res.G    <- list()
        netnoise  <- netdat[!noise,, drop=FALSE]
        ix        <- 0
        while(!identical(tmp.z, z.tmp) && ix <= max.init) {
          tmp.z   <- z.tmp
          ix      <- ix  + 1
          for(k in Gseq) {
            sub   <- z.tmp[,k] == 1
            exp   <- tryCatch(stats::lm(expert, data=netnoise, subset=sub), error=function(e) {
                     try(stats::lm(stats::as.formula(paste("X", paste(names(netnoise[,!apply(netnoise[sub,, drop=FALSE], 2, function(x) all(x == x[1], na.rm=TRUE)), drop=FALSE]), collapse="+"), sep="~")), data=netnoise, subset=sub), silent=TRUE) })
            if(inherits(exp, "try-error")) {
              exp.init        <- FALSE
              break
            }
            res   <- as.matrix((X - tryCatch(stats::predict(exp, newdata=netdat), error=function(e) stats::predict(exp, newdata=.drop_levels(exp, netdat))))[!noise,, drop=FALSE])
            res.G[[k]]        <- res
            mahala[[k]]       <- .MoE_mahala(exp, res, uni)
          }
          if(!exp.init) {
            break
          } else if(g   > 1)   {
            maha  <- do.call(cbind, mahala)
            maha[is.na(maha)] <- Inf
            z.tmp <- maha == apply(maha, 1, min)
          }
        }
        G.res     <- if(uni) as.matrix(do.call(base::c, res.G)) else do.call(rbind, res.G)
      }

      # Account for Noise Component
      z.tmp       <- 0 + z.tmp
      if(noise.null) {
        z         <- z.init   <- z.tmp
      } else   {
        if(g   > 0)  {
          z[!noise,-gN]       <- z.tmp
          z[noise,  gN]       <- 1
        } else {
          z[]     <- 1
        }
        z.init    <- z
      }
      if(exp.init) {
        for(k in Gseq) z.alloc[(k - 1) * n + Nseq,k] <- z.init[,k]
      }
      col.z       <- colSums(z.init)
      emptyinit   <- FALSE
      if(any(col.z < 1)) {                        warning(paste0("For the ", g, " component models, one or more components was empty after initialisation"), call.=FALSE)
        emptyinit <- TRUE
      } else if(any(col.z < 2))                   warning(paste0("For the ", g, " component models, one or more components was initialised with only 1 observation"), call.=FALSE)

    # Initialise gating network
      if(gate.g)  {
        g.init    <- multinom(gating, trace=FALSE, data=netdat)
       #g.init    <- glmnet::cv.glmnet(y=z, x=model.matrix(gating)[,-1], family="multinomial", type.multinomial="grouped")
        tau       <- stats::predict(g.init, type="probs")
       #tau       <- predict(g.init, type="response", newx=model.matrix(gating)[,-1], s="lambda.1se")[,,1]
        ltau      <- log(tau)
        gate.pen  <- length(stats::coef(g.init))  + ifelse(noise.null, 0, 1)
       #gate.pen  <- sum(lengths(coef(g.init, s="lambda.1se")[-1]))
      } else      {
        tau       <- if(equal.pro)  rep(1/gN, gN) else colMeans(z.init)
        ltau      <- .mat_byrow(log(tau), nrow=n, ncol=gN)
        gate.pen  <- ifelse(equal.pro, 0, gN - 1) + ifelse(noise.null, 0, 1)
      }
      ltau.init   <- ltau

    # Loop over the mclust model type(s)
      for(modtype in if(g > 1)   mfg    else if(g == 1) mf1 else mf0)    {
        m0W       <- m0X      <- ERR <- FALSE

      # Initialise parameters from allocations
        if(isTRUE(verbose))      cat(paste0("\n\tModel: ", modtype, "\n"))
        x.df      <- ifelse(g  > 0, nVarParams(modtype, d, g), 0) + gate.pen
        if(exp.init)   {
         Mstep    <- try(mstep(modtype, G.res, z.alloc, control=control), silent=TRUE)
         exp.init <- ifelse(inherits(Mstep, "try-error"), FALSE, attr(Mstep, "returnCode") >= 0)
         mus      <- muX
        }
        if(g > 0  && !exp.init) {
         Mstep    <- try(mstep(modtype, X, if(noise.null) z.init else z.init[,-gN, drop=FALSE], control=control), silent=TRUE)
         ERR      <- inherits(Mstep, "try-error")        || attr(Mstep, "returnCode")  < 0
        }
        if(g > 0  && !ERR)      {
          mus     <- Mstep$parameters$mean
          vari    <- Mstep$parameters$variance
          sigs    <- vari$sigma
        } else mus            <- matrix(NA, nrow=n, ncol=0)

        densme    <- utils::capture.output(medensity     <- try(MoE_dens(modelName=modtype, data=if(exp.init) res.G else x.dat, mus=mus, sigs=sigs, log.tau=ltau.init, Vinv=Vinv), silent=TRUE))
        if((ERR   <- ERR || (g  > 0 && attr(Mstep, "returnCode") < 0) || inherits(medensity, "try-error"))) {
          ll      <- NA
          j       <- 1
          if(isTRUE(verbose))    cat(paste0("\t\t# Iterations: ", ifelse(ERR, "stopped at ", ""), j, "\n"))
          next
        } else    {
          Estep   <- MoE_estep(Dens=medensity)
          z       <- Estep$z
          if((ERR <- any(is.nan(z))))             next
          ll      <- c(-Inf, ifelse(g <= 1 && !exp.g, Estep$loglik, -sqrt(.Machine$double.xmax)))
          linf    <- rep(Inf, 2)
          j       <- 2
          stX     <- gN  > 1 || exp.g
        }

      # Run the EM algorithm
        while(stX)    {

        # Expert network
          if(exp.g)   {
            e.fit <- e.res    <- list()
            for(k in Gseq) {
             fitE <- stats::lm(expert,  weights=z[,k], data=netdat)
            #fitE <- glmnet::cv.glmnet(y=X, x=model.matrix(expert), weights=z[,k])
             e.fit[[k]]       <- fitE
            #e.fit[[k]]       <- coef(fitE, s="lambda.1se")
             e.res[[k]]       <- stats::residuals(fitE)
            #e.res[[k]]       <- X - predict(fitE, type="response", newx=model.matrix(expert), s="lambda.1se")[,,1]
             z.mat[(k - 1) * n + Nseq,k]       <- z[,k]
            }
            res.x <- if(uni) as.matrix(do.call(base::c, e.res)) else do.call(rbind, e.res)
          }

        # M-step
          Mstep   <- if(exp.g) mstep(modtype, res.x, z.mat, control=control) else mstep(modtype, X, if(noise.null) z else z[,-gN, drop=FALSE], control=control)
          mus     <- if(exp.g) muX else Mstep$parameters$mean
          vari    <- Mstep$parameters$variance
          sigs    <- vari$sigma

        # Gating Network
          if(gate.g)  {
            fitG  <- multinom(gating, trace=FALSE, data=netdat)
           #fitG  <- glmnet::cv.glmnet(y=z, x=model.matrix(gating)[,-1], family="multinomial", type.multinomial="grouped")
            tau   <- stats::predict(fitG, type="probs")
           #tau   <- predict(fitG, type="response", newx=model.matrix(gating)[,-1], s="lambda.1se")[,,1]
            ltau  <- log(tau)
           #gate.pen          <- sum(lengths(coef(g.init, s="lambda.1se")[-1]))
          } else  {
            tau   <- if(equal.pro)             tau else if(noise.null) Mstep$parameters$pro else colMeans(z)
            tau   <- if(!exp.g || !noise.null) tau else tau/sum(tau)
            ltau  <- if(equal.pro)            ltau else .mat_byrow(log(tau), nrow=n, ncol=gN)
          }

        # E-step & record log-likelihood
          densme  <- utils::capture.output(medensity <- try(MoE_dens(modelName=modtype, data=if(exp.g) e.res else x.dat, mus=mus, sigs=sigs, log.tau=ltau, Vinv=Vinv), silent=TRUE))
          if((ERR <- attr(Mstep, "returnCode")  < 0  || inherits(medensity, "try-error"))) {
            ll    <- c(ll, NA)
            break
          } else  {
            Estep <- MoE_estep(Dens=medensity)
            z     <- Estep$z
            ERR   <- any(is.nan(z))
            if(isTRUE(ERR))                       break
            ll    <- c(ll, Estep$loglik)
            j     <- j + 1
            if(stopx  == "aitken")  {
             ait  <- MoE_aitken(ll[seq(j - 2, j, 1)])
             linf <- c(linf[2], ait$linf)
             dX   <- ifelse(is.numeric(ait$a) && ait$a < 0, 0, abs(diff(linf)))
             dX[is.nan(dX)]   <- Inf
            } else     {
             dX   <- abs((ll[j] - ll[j - 1])/(1 + ll[j]))
            }
            stX   <- dX >= tol && j < max.it && g > 1
            if(itwarn && !m0X)  {
             m0W  <- ifelse(!m0X, warnit < j, m0X)
             if(m0W   && !m0X)  {                 tryCatch(warning("WARNIT", call.=FALSE), warning=function(w)
                                                  message(paste0("\tEM algorithm for the ", modtype, " model has yet to converge in 'warn.it'=", warnit, " iterations")))
              m0X <- TRUE
             }
            }
          }
        } # while (j)

      # Store values corresponding to the maximum BIC/ICL/AIC so far
        j2        <- max(1, j  - 2)
        if(isTRUE(verbose))      cat(paste0("\t\t# Iterations: ", ifelse(ERR, "stopped at ", ""), j2, "\n"))
        x.df      <- ifelse(exp.g, g * length(stats::coef(fitE)), exp.pen) + x.df
        ll[j]     <- ifelse(g <= 1, ll[j], switch(stopx, aitken=max(ll[j], ifelse(is.finite(linf[2]), linf[2], linf[1])), ll[j]))
        choose    <- MoE_crit(modelName=modtype, loglik=ll[j], n=n, G=g, z=z, df=x.df)
        bics      <- choose["bic",]
        icls      <- choose["icl",]
        aics      <- choose["aic",]
        crit.t    <- switch(criterion, bic=bics, icl=icls, aic=aics)
        crit.t    <- ifelse(is.na(crit.t) || ERR, -Inf, crit.t)
        if(crit.t  > crit.tx)   {
          crit.tx <- crit.t
          tau.x   <- tau
          z.x     <- z
          ll.x    <- ll
          df.x    <- x.df
          if(gate.g)   {
            gfit  <- fitG
          }
          if(exp.g)    {
            efit  <- e.fit
            eres  <- res.x
            mu.x  <- mus
            sig.x <- vari
          }
        }
        BICs[h,modtype]       <- ifelse(ERR, -Inf, bics)
        ICLs[h,modtype]       <- ifelse(ERR, -Inf, icls)
        AICs[h,modtype]       <- ifelse(ERR, -Inf, aics)
        DF.x[h,modtype]       <- ifelse(ERR, -Inf, x.df)
        it.x[h,modtype]       <- ifelse(ERR,  Inf, j2)
      } # for (modtype)

    # Pull out mclust model corresponding to highest BIC/ICL/AIC
      if(crit.tx   > crit.gx)   {
        crit.gx   <- crit.tx
        x.tau     <- tau.x
        x.z       <- z.x
        x.ll      <- ll.x
        x.DF      <- df.x
        if(gate.g)     {
          x.fitG  <- gfit
        }
        if(exp.g)      {
          x.fitE  <- efit
          x.resE  <- eres
          x.mu    <- mu.x
          x.sig   <- sig.x
        }
      }
    } # for (g)
    if(all(is.infinite(BICs[!is.na(BICs)])))      stop("All models failed!")

  # Gather results, fit extra gating & expert networks, and extra M-step
    CRITs         <- switch(criterion, bic=BICs, icl=ICLs, aic=AICs)
    best.ind      <- which(CRITs == crit.gx, arr.ind=TRUE)
    G             <- G[best.ind[1]]
    GN            <- G + !noise.null
    z             <- x.z
    rownames(z)   <- as.character(seq_len(n))
    best.mod      <- colnames(CRITs)[best.ind[2]]
    bic.fin       <- BICs[best.ind]
    icl.fin       <- ICLs[best.ind]
    aic.fin       <- AICs[best.ind]
    uncertainty   <- if(GN > 1) 1 - apply(z, 1, max)         else rep(0, n)
    exp.x         <- exp.x & G != 0
    bG            <- gate.G[which(range.G == G)]
    exp.gate      <- c(exp.x, bG)
    net.msg       <- ifelse(any(exp.gate), paste0(" (incl. ", ifelse(all(exp.gate), "gating and expert", ifelse(exp.x, "expert", ifelse(bG, "gating", ""))), paste0(" network covariates", ifelse(noise.null, ")", ", and a noise component)"))), ifelse(noise.null, "", " (and a noise component)"))
    x.ll          <- x.ll[if(G == 0 || (G == 1 && !exp.x)) 2 else if(G == 1 && exp.x) seq_len(3)[-1]      else -c(1:2)]
    x.ll          <- x.ll[!is.na(x.ll)]

    x.fitG        <- if(bG)            x.fitG  else if(GN > 1) multinom(gating, trace=FALSE, data=netdat) else suppressWarnings(stats::glm(z ~ 1, family=stats::binomial()))
    x.tau         <- if(GN > 1 && !bG) x.fitG$fitted.values[1,]                                           else x.tau
    if(!exp.x)    {
     x.fitE       <- list()
     for(k in seq_len(GN)) {
      x.fitE[[k]] <- stats::lm(expert, weights=z[,k], data=netdat)
     }
    }
    x.fitE        <- stats::setNames(x.fitE, paste0("Cluster", seq_len(G)))
    attr(x.fitG, "EqualPro")  <- equal.tau[best.ind[1]]
    attr(x.fitG, "Formula")   <- Reduce(paste, deparse(gating[-2]))
    attr(x.fitE, "Formula")   <- Reduce(paste, deparse(expert[-2]))
    class(x.fitG) <- c("MoE_gating", class(x.fitG))
    class(x.fitE) <- c("MoE_expert", class(x.fitE))
    if(g > 0) {
      extraM      <- mstep(best.mod, X, if(noise.null) z     else z[,-GN, drop=FALSE], control=control)
      mean.fin    <- extraM$parameters$mean
      vari.fin    <- if(exp.x) x.sig else extraM$parameters$variance
    } else    {
      mean.fin    <- vari.fin <- NULL
    }

    if(multi      && verbose)    cat(paste0("\n\t\tBest Model: ", mclustModelNames(best.mod)$type, " (", best.mod, "), with ", ifelse(G == 0, "only a noise component",
                                     paste0(G, " component", ifelse(G > 1, "s", ""))), ifelse(any(exp.gate) || (!noise.null && G != 0), paste0("\n\t\t\t   ", net.msg), ""), "\n\t\t",
                                     switch(criterion, bic="BIC", icl="ICL", aic="AIC"), " = ", round(switch(criterion, bic=bic.fin, icl=icl.fin, aic=aic.fin), 2), "\n"))
    if(G == 1     && gate.x)   {
      exp.names   <- attr(netdat, "Expert")
      if(attr(x.fitG, "Formula") != "None") tmpnet                    <- netdat
      netdat      <- netdat[,setdiff(gsub("[[:space:]]", ".", gsub("[[:punct:]]", ".", attr(netdat, "Gating")[!is.na(attr(netdat, "Gating"))])), colnames(netdat))]
      attr(netdat, "Gating")  <-
      attr(netdat, "Both")    <- NA
      attr(netdat, "Expert")  <- exp.names
      if(attr(x.fitG, "Formula") != "None") attr(netdat, "Discarded") <- tmpnet
    }
    if(all(x.ll   != cummax(x.ll)))               warning("Log-likelihoods are not strictly increasing", call.=FALSE)
    if(any(it.x[!is.na(it.x)] == max.it))         warning(paste0("One or more models failed to converge in the maximum number of allowed iterations (", max.it, ")"), call.=FALSE)
    class(BICs)   <- c("MoECriterion", "mclustBIC")
    class(ICLs)   <- c("MoECriterion", "mclustICL")
    class(AICs)   <- c("MoECriterion", "mclustAIC")
    class(DF.x)   <- c("MoECriterion", "mclustDF")
    attr(BICs, "G")           <-
    attr(ICLs, "G")           <-
    attr(AICs, "G")           <-
    attr(DF.x, "G")           <- rownames(BICs)
    attr(BICs, "modelNames")  <-
    attr(ICLs, "modelNames")  <-
    attr(AICs, "modelNames")  <-
    attr(DF.x, "modelNames")  <- colnames(BICs)
    if(!noise.null) {
      attr(BICs, "Vinv")      <-
      attr(ICLs, "Vinv")      <-
      attr(AICs, "Vinv")      <-
      attr(DF.x, "Vinv")      <- Vinv
    }
    attr(BICs, "control")     <-
    attr(ICLs, "control")     <-
    attr(AICs, "control")     <-
    attr(DF.x, "control")     <- control
    attr(BICs, "warn")        <-
    attr(ICLs, "warn")        <-
    attr(AICs, "warn")        <-
    attr(DF.x, "warn")        <- FALSE
    attr(BICs, "n")           <-
    attr(ICLs, "n")           <-
    attr(AICs, "n")           <-
    attr(DF.x, "n")           <- n
    attr(BICs, "d")           <-
    attr(ICLs, "d")           <-
    attr(AICs, "d")           <-
    attr(DF.x, "d")           <- d
    attr(BICs, "oneD")        <-
    attr(ICLs, "oneD")        <-
    attr(AICs, "oneD")        <-
    attr(DF.x, "oneD")        <- uni
    attr(BICs, "criterion")   <- "BIC"
    attr(ICLs, "criterion")   <- "ICL"
    attr(AICs, "criterion")   <- "AIC"
    attr(DF.x, "criterion")   <- "DF"
    attr(BICs, "returnCodes") <- provideDimnames(unname(ifelse(is.na(BICs) | is.infinite(BICs), -1, 0)), base=list(attr(BICs, "G"), colnames(BICs)))
    attr(ICLs, "returnCodes") <- provideDimnames(unname(ifelse(is.na(ICLs) | is.infinite(ICLs), -1, 0)), base=list(attr(ICLs, "G"), colnames(ICLs)))
    attr(AICs, "returnCodes") <- provideDimnames(unname(ifelse(is.na(AICs) | is.infinite(AICs), -1, 0)), base=list(attr(AICs, "G"), colnames(AICs)))
    attr(DF.x, "returnCodes") <- provideDimnames(unname(ifelse(is.na(DF.x) | is.infinite(DF.x), -1, 0)), base=list(attr(DF.x, "G"), colnames(DF.x)))
    attr(BICs, "initialization")      <-
    attr(ICLs, "initialization")      <-
    attr(AICs, "initialization")      <-
    attr(DF.x, "initialization")      <- list(hcPairs = NULL, subset = NULL, noise = noise)
    claX       <- max.col(z)
    claX[claX  == G + 1]      <- 0
    results       <- list(call = call, data = as.data.frame(X), modelName = best.mod, n = n, d = d, G = G, BIC = BICs, ICL = ICLs, AIC = AICs,
                          bic = bic.fin, icl = icl.fin, aic = aic.fin, gating = x.fitG, expert = x.fitE, loglik = x.ll, df = x.DF, hypvol = ifelse(noise.null, NA, 1/Vinv),
                          parameters = list(pro = x.tau, mean = mean.fin, variance = vari.fin, Vinv = if(!noise.null) Vinv), z = z, classification = stats::setNames(claX, seq_len(n)),
                          uncertainty = stats::setNames(uncertainty, seq_len(n)), net.covs = netdat, resid.data = if(exp.x) x.resE, DF = DF.x, iters = it.x)
    class(results)            <- "MoEClust"
    attr(results, "Details")  <- paste0(best.mod, ": ", ifelse(G == 0, "only a noise component", paste0(G, " component", ifelse(G > 1, "s", ""), net.msg)))
    attr(results, "EqualPro") <- equal.pro
    attr(results, "Expert")   <- exp.x
    attr(results, "Gating")   <- bG
      return(results)
  }

#' Density for Parameterised MVN Mixture of Experts Models
#'
#' Computes densities (or log-densities) of observations in parameterised MVN mixtures of experts.
#' @param modelName A character string indicating the model. The help file for \code{\link[mclust]{mclustModelNames}} describes the available models.
#' @param data If there are no expert network covariates, \code{data} should be a numeric matrix or data frame, wherein rows correspond to observations (n) and columns correspond to variables (d). If there are expert network covariates, this should be a list of length G containing matrices/data.frames of (multivariate) WLS residuals for each component.
#' @param mus The mean for each of G components. If there is more than one component, this is a matrix whose k-th column is the mean of the k-th component of the mixture model. For the univariate models, this is a G-vector of means. In the presence of expert network covariates, all values should be equal to zero.
#' @param sigs A list of length G of variance parameters of the model. The components of this list depend on the specification of \code{modelName}.
#' @param log.tau If covariates enter the gating network, an n times G matrix of mixing proportions, otherwise a G-vector of mixing proportions for the components of the mixture. \strong{Must} be on the log-scale in both cases. The default of \code{0} effectively means densities (or log-densities) aren't scaled by the mixing proportions.
#' @param Vinv An estimate of the reciprocal hypervolume of the data region. The default is determined by applying the function \code{\link[mclust]{hypvol}} to the data. Used only if an initial guess as to which observations are noise is supplied. Mixing proportion(s) must be included for the noise component also.
#' @param logarithm A logical value indicating whether or not the logarithm of the component densities should be returned. This defaults to \code{TRUE}, otherwise component densities are returned, obtained from the component log-densities by exponentiation. The \strong{log}-densities can be passed to \code{\link{MoE_estep}}.
#'
#' @note This function is intended for joint use with \code{\link{MoE_estep}}, using the \strong{log}-densities.
#' @importFrom mclust "mclustModelNames"
#' @importFrom mvnfast "dmvn"
#' @return A numeric matrix whose \code{[i,k]}-th entry is the density or log-density of observation \emph{i} in component \emph{k}, scaled by the mixing proportions.
#' @export
#' @author Keefe Murphy - \href{keefe.murphy@ucd.ie}{<keefe.murphy@ucd.ie>}
#'
#' @seealso \code{\link{MoE_estep}}, \code{\link{MoE_clust}}, \code{\link[mclust]{mclustModelNames}}
#'
#' @examples
#' data(ais)
#' hema  <- ais[,3:7]
#' model <- MoE_clust(hema, G=3, gating= ~ ais$BMI + ais$sex, model="EEE")
#' Dens  <- MoE_dens(modelName=model$modelName, data=hema,
#'                   mus=model$parameters$mean, sigs=model$parameters$variance$sigma,
#'                   log.tau=log(model$parameters$pro))
#'
#' # Construct the z matrix and compute the log-likelihood
#' Estep <- MoE_estep(Dens=Dens)
#' (ll   <- Estep$loglik)
#'
#' # The z matrix will be close but not exactly the same as that from the model
#' # as the EM algorithm finishes on an M-step, but the classification should be
#' identical(max.col(Estep$z), as.integer(unname(model$classification))) #TRUE
#' round(sum(Estep$z - model$z), options()$digits) == 0                  #TRUE
  MoE_dens        <- function(modelName, data, mus, sigs, log.tau = 0, Vinv = NULL, logarithm = TRUE) {
    G             <- tryCatch(ifelse(is.matrix(mus), ncol(mus), length(mus)), error=function(e) 0)
    Vnul          <- is.null(Vinv)
    Gseq          <- seq_len(G)
    if(!is.list(data)    || (is.list(data) &&
        length(data)     != max(G, 1)))     {
      data        <- replicate(G, as.matrix(data), FALSE)
    }
    dat1          <- data[[1]]
    n             <- ifelse(is.matrix(dat1), nrow(dat1), length(dat1))
    d             <- ifelse(is.matrix(dat1), ncol(dat1), 1)
    sq_mat        <- if(d > 50) function(x)  diag(sqrt(diag(x)))     else sqrt
    if(G > 0) {
      switch(EXPR=modelName, EVE=, VEE=, VVE=, EEV=, VEV=, EVV=, VVV = {
        idens     <- utils::capture.output(densi <- try(vapply(Gseq, function(k) dmvn(data[[k]], mus[,k], sigs[,,k],         log=TRUE, isChol=FALSE), numeric(n)), silent=TRUE))
      }, VII=, VEI=, EVI=, VVI = {
        idens     <- utils::capture.output(densi <- try(vapply(Gseq, function(k) dmvn(data[[k]], mus[,k], sq_mat(sigs[,,k]), log=TRUE, isChol=TRUE),  numeric(n)), silent=TRUE))
      }, EII=, EEI = {
        sigx      <- sq_mat(sigs[,,1]);
        idens     <- utils::capture.output(densi <- try(vapply(Gseq, function(k) dmvn(data[[k]], mus[,k], sigx,              log=TRUE, isChol=TRUE),  numeric(n)), silent=TRUE))
      }, EEE= {
        sigx      <- .chol(sigs[,,1]) ;
        idens     <- utils::capture.output(densi <- try(vapply(Gseq, function(k) dmvn(data[[k]], mus[,k], sigx,              log=TRUE, isChol=TRUE),  numeric(n)), silent=TRUE))
      }, E=   {
        densi     <- vapply(Gseq, function(k)       stats::dnorm(data[[k]], mus[k], sqrt(sigs), log=TRUE), numeric(n))
      }, V=   {
        sigx      <- sqrt(sigs);
        densi     <- vapply(Gseq, function(k)       stats::dnorm(data[[k]], mus[k], sigx[k],    log=TRUE), numeric(n))
      } )
      test        <- is.infinite(densi) & densi   > 0
      densi[test] <- 0
    }
    densi         <- if(Vnul) densi else if(G     > 0) cbind(densi, log(Vinv)) else matrix(log(Vinv), nrow=n, ncol=G + !Vnul)
    densi         <- densi  + if(is.matrix(log.tau) || missing(log.tau))                              log.tau else
                              if(length(log.tau)    == G + !Vnul) .mat_byrow(log.tau, nrow=n, ncol=G + !Vnul) else stop(paste0("'log.tau' must be given for every component", ifelse(Vnul, "", ", incl. the noise component if 'Vinv' is supplied")))
      if(logarithm)  densi    else exp(densi)
  }

#' Compute the responsility matrix z and log-likelihood for MoEClust models
#'
#' Softmax function to compute the responsibility matrix z and the log-likelihood for parameterised MVN mixtures of experts, with the aid of \code{\link{MoE_dens}}.
#' @inheritParams MoE_dens
#' @param Dens (Optional) A numeric matrix whose \code{[i,k]}-th entry is the \strong{log}-density of observation \emph{i} in component \emph{k}, scaled by the mixing proportions, to which the softmax function is to be applied, typically obtained by \code{\link{MoE_dens}} but this is not necessary. If this is supplied, all other arguments are ignored, otherwise \code{\link{MoE_dens}} is called according to the other supplied arguments.
#'
#' @return A list containing two elements:
#' \itemize{
#' \item{\code{z} - }{A matrix with n rows and G columns containing the probability of cluster membership for each of n observations and G clusters}
#' \item{\code{loglik} - }{The log-likelihood, computed efficiently via \code{\link[matrixStats]{rowLogSumExps}}}
#' }
#' @importFrom matrixStats "rowLogSumExps"
#' @importFrom mclust "mclustModelNames"
#' @export
#' @note This softmax function is intended for joint use with \code{\link{MoE_dens}}, using the \strong{log}-densities.
#' @author Keefe Murphy - \href{keefe.murphy@ucd.ie}{<keefe.murphy@ucd.ie>}
#'
#' @seealso \code{\link{MoE_dens}}, \code{\link{MoE_clust}}, \code{\link[matrixStats]{rowLogSumExps}}, \code{\link[mclust]{mclustModelNames}}
#'
#' @examples
#' data(ais)
#' hema   <- ais[,3:7]
#' model  <- MoE_clust(hema, G=3, gating= ~ ais$BMI + ais$sex, model="EEE")
#' Dens   <- MoE_dens(modelName=model$modelName, data=hema,
#'                    mus=model$parameters$mean, sigs=model$parameters$variance$sigma,
#'                    log.tau=log(model$parameters$pro))
#'
#' # Construct the z matrix and compute the log-likelihood
#' Estep  <- MoE_estep(Dens=Dens)
#' (ll    <- Estep$loglik)
#'
#' # The z matrix will be close but not exactly the same as that from the model
#' # as the EM algorithm finishes on an M-step, but the classification should be
#' identical(max.col(Estep$z), as.integer(unname(model$classification))) #TRUE
#' round(sum(Estep$z - model$z), options()$digits) == 0                  #TRUE
#'
#' # Call MoE_estep directly
#' Estep2 <- MoE_estep(modelName=model$modelName, data=hema,
#'                     mus=model$parameters$mean, sigs=model$parameters$variance$sigma,
#'                     log.tau=log(model$parameters$pro))
#' identical(Estep2$loglik, ll)                                          #TRUE
  MoE_estep       <- function(modelName, data, mus, sigs, log.tau = 0, Vinv = NULL, Dens = NULL) {
    if(missing(Dens)) {
      Dens        <- do.call(MoE_dens, as.list(match.call())[-1])
    } else if(!is.matrix(Dens) ||
              !is.numeric(Dens))                  stop("'Dens' must be a numeric matrix")
    norm          <- rowLogSumExps(Dens)
    z             <- exp(sweep(Dens, 1, norm, "-"))
   #ll            <- sum(z * Dens)              # complete log-likelihood
      return(list(z = z, loglik = sum(norm)))
  }

#' MoEClust BIC, ICL, and AIC Model-Selection Criteria
#'
#' Computes the BIC (Bayesian Information Criterion), ICL (Integrated Complete Likelihood), and AIC (Akaike Information Criterion) for parameterized mixture of experts models given the log-likelihood, the dimension of the data, the number of mixture components in the model, the numbers of parameters in the gating and expert networks respectively, and, for the ICL, the numbers of observations in each component.
#' @param modelName A character string indicating the model. The help file for \code{\link[mclust]{mclustModelNames}} describes the available models.
#' @param loglik The log-likelihood for a data set with respect to the Gaussian mixture model specified in the \code{modelName} argument.
#' @param n,d,G The number of observations in the data, dimension of the data, and number of components in the Gaussian mixture model, respectively, used to compute \code{loglik}.
#' @param gating.pen The number of parameters of the \emph{gating} network of the MoEClust model. Defaults to \code{G - 1}, which corresponds to no gating covariates. If covariates are included, this should be the number of regression coefficients in the fitted object. If there are no covariates and mixing proportions are further assumed to be present in equal proportion, \code{gating.pen} should be \code{0}.
#' @param expert.pen The number of parameters of the \emph{expert} network of the MoEClust model. Defaults to \code{G * d}, which corresponds to no expert covariates. If covariates are included, this should be the number of regression coefficients in the fitted object.
#' @param z The \code{n} times \code{G} responsibility matrix whose \code{[i,k]}-th entry is the probability that observation \emph{i} belonds to the \emph{k}-th component.. If supplied the ICL is also computed and returned, otherwise only the BIC and AIC.
#' @param df An alternative way to specify the number of estimated parameters (or 'used' degrees of freedom) exactly. If supplied, the arguments \code{d, gating.pen} and \code{expert.pen}, which are used to calculate the number of parameters, will be ignored. The number of parameters used in the estimation of the noise component, if any, should also be included.
#' @param delta Dirichlet hyperparameter for the prior on the mixing proportions. Defaults to 0.5. Only relevant for the ICL computation.
#'
#' @details The function is vectorized with respect to the arguments \code{modelName} and \code{loglik}.\cr
#'
#' If \code{model} is an object of class \code{"MoEClust"} with \code{G} components, the number of parameters for the \code{gating.pen} and \code{expert.pen} are \code{length(coef(model$gating))} and \code{G * length(coef(model$expert[[1]]))}, respectively.
#' @importFrom mclust "mclustModelNames" "nVarParams"
#' @return A simplified array containing the BIC, AIC, number of estimated parameters (\code{df}) and, if \code{z} is supplied, also the ICL, for each of the given input arguments.
#' @note In order to speed up repeated calls to the function inside \code{\link{MoE_clust}}, no checks take place.
#' @export
#' @author Keefe Murphy - \href{keefe.murphy@ucd.ie}{<keefe.murphy@ucd.ie>}
#'
#' @references Biernacki, C., Celeux, G., Govaert, G. (2000). Assessing a mixture model for clustering with the integrated completed likelihood. \emph{IEEE Trans. Pattern Analysis and Machine Intelligence}, 22(7): 719-725.
#'
#' @seealso \code{\link{MoE_clust}}, \code{\link[mclust]{nVarParams}}, \code{\link[mclust]{mclustModelNames}}
#'
#' @examples
#' MoE_crit(modelName=c("VVI", "VVE", "VVV"), n=120, d=8,
#'          G=3, loglik=c(-4036.99, -3987.12, -3992.45))
#'
#' data(CO2data)
#' model <- MoE_clust(CO2data$CO2, G=1:2, expert= ~ CO2data$GNP)
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
  MoE_crit        <- Vectorize(function(modelName, loglik, n, d, G, gating.pen = G - 1, expert.pen = G * d, z = NULL, df = NULL, delta = 0.5) {
    df            <- ifelse(!missing(df), df, nVarParams(modelName, d, G) + expert.pen + gating.pen)
    double.ll     <- 2 * loglik
    bic.x         <- double.ll  - df * log(n)
    aic.x         <- double.ll  - df * 2
      return(c(bic = bic.x, icl = if(!missing(z)) bic.x + 2 * sum(log(apply(z, 1, max)), na.rm = TRUE), aic = aic.x, df = df))
  }, vectorize.args = c("modelName", "loglik"), SIMPLIFY="array")

#' Set control values for use with MoEClust
#'
#' Supplies a list of arguments (with defaults) for use with \code{\link{MoE_clust}}.
#' @param criterion When either \code{G} or \code{modelNames} is a vector, \code{criterion} determines whether the "\code{bic}" (Bayesian Information Criterion), "\code{icl}" (Integrated Complete Likelihood), "\code{aic}" (Akaike Information Criterion) is used to determine the 'best' model when gathering output. Note that all criteria will be returned in any case.
#' @param stopping The criterion used to assess convergence of the EM algorithm. The default (\code{"aitken"}) uses Aitken's acceleration method via \code{\link{MoE_aitken}}, otherwise the \code{"relative"} change in log-likelihood is monitored (which may be less strict). Both stopping rules are ultimately governed by \code{tol[1]}. When the \code{"aitken"} method is employed, the final estimate of the log-likelihood is the value of \code{linf} at convergence, rather than the value of \code{ll} at convergence under the \code{"relative"} option.
#' @param exp.init A list supplying select parameters to control the initialisation routine in the presence of expert network covariates (otherwise ignored):
#' \itemize{
#' \item{\code{"joint"} - }{A logical indicating whether the initial partition is obtained on the joint distribution of the response and (continuous only) covariates (defaults to \code{TRUE}) or just the response variables (\code{FALSE}). Only relevant when \code{init.z} is not \code{"random"}. This may render the \code{"quantile"} option to \code{init.z} for univariate data unusable.}
#' \item{\code{"mahalanobis"} - }{A logical indicating whether to iteratively reallocate observations during the initialisation phase to the component corresponding to the expert network regression to which its residual is closest in terms of Mahalanobis distance (defaults to \code{TRUE}). This will ensure that each component can be well modelled by a single expert prior to running the EM algorithm.}
#' \item{\code{"max.init"} - }{The maximum number of iterations for the Mahalanobis distance-based reallocation procedure when \code{exp.init$mahalanobis} is \code{TRUE}. Defaults to \code{100}.}
#' }
#' @param init.z The method used to initialise the cluster labels. Defaults to a hierarchical clustering tree as per \code{\link[mclust]{hc}} for multivariate data, or quantile-based clustering as per \code{\link{MoE_qclass}} for univariate data (unless there are continuous expert network covariates, in which case the defaults is again \code{\link[mclust]{hc}}). The \code{"quantile"} option is only available for univariate data without continuous expert network covariates. Other options include \code{kmeans}, \code{random} initialisation, and a full run of \code{\link[mclust]{Mclust}}, although this last option is only permitted if there are \code{gating} &/or \code{expert} covariates within \code{\link{MoE_clust}}.
#' @param eps A scalar tolerance associated with deciding when to terminate computations due to computational singularity in covariances. Smaller values of eps allow computations to proceed nearer to singularity. The default is the relative machine precision \code{.Machine$double.eps}, which is approximately \emph{2e-16} on IEEE-compliant machines.
#' @param tol A vector of length two giving relative convergence tolerances for the log-likelihood and for parameter convergence in the inner loop for models with iterative M-step ("VEI", "EVE", "VEE", "VVE", "VEV"), respectively. The default is \code{c(1e-05,sqrt(.Machine$double.eps))}. If only one number is supplied, it is used as the tolerance for the outer iterations and the tolerance for the inner iterations is as in the default.
#' @param itmax A vector of length two giving integer limits on the number of EM iterations and on the number of iterations in the inner loop for models with iterative M-step ("VEI", "EVE", "VEE", "VVE", "VEV"), respectively. The default is \code{c(.Machine$integer.max, .Machine$integer.max)} allowing termination to be completely governed by \code{tol}. If only one number is supplied, it is used as the iteration limit for the outer iteration only.
#' @param equalPro Logical variable indicating whether or not the mixing proportions are equal in the model. Default: \code{equalPro = FALSE}. Only relevant when \code{gating} covariates are \emph{not} supplied within \code{\link{MoE_clust}}, otherwise ignored.
#' @param warn.it A single number giving the iteration count at which a warning will be printed if the EM algorithm has failed to converge. Defaults to \code{0}, i.e. no warning (which is true for any \code{warn.it} value less than \code{3}), otherwise the message is printed regardless of the value of \code{verbose}. If non-zero, \code{warn.it} should be moderately large, but obviously less than \code{itmax[1]}. A warning will always be printed if one of more models fail to converge in \code{itmax[1]} iterations.
#' @param noise.init A logical or numeric vector indicating an initial guess as to which observations are noise in the data. If numeric, the entries should correspond to row indices of the data. If supplied, a noise term will be added to the model in the estimation.
#' @param noise.meth The method use to estimate the volume when observations are labelled as noise via \code{noise.init}. Defaults to \code{\link[mclust]{hypvol}}. For univariate data, this argument is ignored and the range of the data is used instead. The options \code{convexhull} and \code{ellipsoidhull} require loading the \code{geometry} and \code{cluster} libraries, respectively.
#' @param hc.meth A character string indicating the model to be used when hierarchical clustering (see \code{\link[mclust]{hc}} is employed for initialisation according to \code{init.z}. Defaults to \code{"EII"} for high-dimensional data, or \code{"VVV"} otherwise.
#' @param verbose Logical indicating whether to print messages pertaining to progress to the screen during fitting. By default is \code{TRUE} if the session is interactive, and \code{FALSE} otherwise. If \code{FALSE}, warnings and error messages will still be printed to the screen, but everything else will be suppressed.
#' @param ... Catches unused arguments.
#'
#' @details \code{\link{MoE_control}} is provided for assigning values and defaults within \code{\link{MoE_clust}}.\cr
#'
#' While the \code{criterion} argument controls the choice of the optimal number of components and \pkg{mclust} model type, \code{\link{MoE_compare}} is provided for choosing between fits with different combinations of covariates or different initialisation settings.
#' @importFrom mclust "hc"
#' @note When initialising using the \code{"hc"} option for the \code{init.z} argument, the \code{"EII"} model is used for high-dimensional data, otherwise \code{"VVV"} is used.
#' @return A named list in which the names are the names of the arguments and the values are the values supplied to the arguments.
#' @export
#' @author Keefe Murphy - \href{keefe.murphy@ucd.ie}{<keefe.murphy@ucd.ie>}
#'
#' @seealso \code{\link{MoE_clust}}, \code{\link{MoE_aitken}}, \code{\link[mclust]{hc}}, \code{\link{MoE_qclass}}, \code{\link[mclust]{hypvol}}, \code{\link[geometry]{convhulln}}, \code{\link[cluster]{ellipsoidhull}}, \code{\link{MoE_compare}}
#'
#' @examples
#' ctrl <- MoE_control(criterion="icl", itmax=100, warn.it=15, init.z="random")
#'
#' data(CO2data)
#' res  <- MoE_clust(CO2data$CO2, G=2, expert = ~ CO2data$GNP, control=ctrl)
#'
#' # Alternatively, specify control arguments directly
#' res2 <- MoE_clust(CO2data$CO2, G=2, expert = ~ CO2data$GNP, stopping="relative")
#'
#' # Supplying ctrl without naming it as control throws an error,
#' # when any of {modelNames, gating, expert} are not supplied
#' \dontrun{
#' res3 <- MoE_clust(CO2data$CO2, G=2, expert = ~ CO2data$GNP, ctrl)}
  MoE_control     <- function(criterion = c("bic", "icl", "aic"), stopping = c("aitken", "relative"), exp.init = list(...),
                              init.z = c("hc", "quantile", "kmeans", "mclust", "random"), eps = .Machine$double.eps,  tol = c(1e-05, sqrt(.Machine$double.eps)),
                              itmax = c(.Machine$integer.max, .Machine$integer.max), equalPro = FALSE, warn.it = 0, noise.init = NULL,
                              noise.meth = c("hypvol", "convexhull", "ellipsoidhull"), hc.meth = NULL, verbose = interactive(), ...) {
    if(!missing(criterion) && (length(criterion) > 1 ||
       !is.character(criterion)))                 stop("'criterion' must be a single character string")
    criterion     <- match.arg(criterion)
    if(!missing(stopping)  && (length(stopping)  > 1 ||
       !is.character(stopping)))                  stop("'stopping' must be a single character string")
    stopping      <- match.arg(stopping)
    if(is.null(exp.init$joint)) {
      exp.init$joint       <- TRUE
    } else if(length(exp.init$joint)       > 1 ||
              !is.logical(exp.init$joint))        stop("'exp.init$joint' must be a single logical indicator")
    if(is.null(exp.init$mahalanobis))      {
      exp.init$mahalanobis <- TRUE
    } else if(length(exp.init$mahalanobis) > 1 ||
              !is.logical(exp.init$mahalanobis))  stop("'exp.init$mahalanobis' must be a single logical indicator")
    if(is.null(exp.init$max.init))         {
      exp.init$max.init    <- 100
    } else if(isTRUE(exp.init$mahalanobis)     &&
             (length(exp.init$max.init)    > 1 ||
             ((!is.numeric(exp.init$max.init)  ||
              exp.init$max.init    <= 0))))       stop("'exp.init$max.init' must be a single strictly positive integer when 'exp.init$mahalanobis' is TRUE")
    miss.init     <- missing(init.z)
    miss.hc       <- missing(hc.meth)
    if(!missing(init.z)    && (length(init.z)    > 1 ||
       !is.character(init.z)))                    stop("'init.z' must be a single character string")
    init.z        <- match.arg(init.z)
    if(length(eps) > 2)                           stop("'eps' can be of length at most 2")
    if(any(eps     < 0))                          stop("'eps' is negative")
    if(any(eps    >= 1))                          stop("'eps' is not less than 1")
    if((len.tol   <- length(tol))   > 2)          stop("'tol' can be of length at most 2")
    if(any(tol     < 0))                          stop("'tol' is negative")
    if(any(tol    >= 1))                          stop("'tol' is not less than 1")
    if(any(itmax   < 0))                          stop("'itmax' is negative")
    if(len.tol    == 1)        tol <- rep(tol, 2)
    if((len.itmax <- length(itmax)) > 2)          stop("'itmax' can be of length at most 2")
    if(len.itmax  == 1)      itmax <- c(itmax, .Machine$integer.max)
    inf           <- is.infinite(itmax)
    if(any(inf))        itmax[inf] <- .Machine$integer.max
    if(length(equalPro) > 1 ||
       !is.logical(equalPro))                     stop("'equalPro' must be a single logical indicator")
    if(length(warn.it)  > 1 ||
       !is.numeric(warn.it))                      stop("'warn.it' must be a numeric vector of length 1")
    noise.meth    <- match.arg(noise.meth)
    if(!is.null(noise.init)) {
     if(!missing(noise.meth) & (length(noise.meth) > 1 ||
        !is.character(noise.meth)))               stop("'noise.meth' must be a single character string")
     has.lib      <- switch(noise.meth, hypvol=TRUE, convexhull=suppressMessages(requireNamespace("geometry", quietly=TRUE)), ellipsoidhull=suppressMessages(requireNamespace("cluster", quietly=TRUE)))
     if(!has.lib)                                 stop(paste0("Use of the ", noise.meth, " 'noise.meth' option requires loading the ", switch(noise.meth, hypvol="'mclust'", convexhull="'geometry'", ellipsoidhull="'cluster'"), "library"))
    }
    if(!miss.hc   && init.z == "hc" && !is.element(hc.meth, c("E", "V",
       mclust.options("hcModelNames"))))          stop("Invalid 'hc.meth' selected for initialisation by agglomerative hierarchical clustering")
    if(length(verbose)  < 1 ||
       !is.logical(verbose))                      stop("'verbose' must be a single logical indicator")
      list(criterion = criterion, stopping = stopping, exp.init = exp.init, init.z = init.z, eps = eps, tol = tol, itmax = itmax, equalPro = equalPro,
           warn.it = warn.it, noise.init = noise.init, noise.meth = noise.meth, hc.meth = hc.meth, verbose = verbose, miss.init = miss.init, miss.hc = miss.hc)
  }

#' Aitken Acceleration
#'
#' Calculates the Aitken acceleration estimate of the final converged maximized log-likelihood.
#' @param loglik A vector of three consecutive log-likelihood values. These three values should be in ascending order, though this is not checked.
#'
#' @details The final converged maximized log-likelihood can be used to determine convergence of the EM algorithm within \code{\link{MoE_clust}}, i.e. by checking whether the absolute difference in the current and previous estimates of the final converged maximised log-likelihood is less than some tolerance.
#' @note Within \code{\link{MoE_clust}}, as specified by the \code{stopping} argument of \code{\link{MoE_control}}, \code{"aitken"} is the default method used to assess convergence. The other option monitors the \code{"relative"} change in log-likelihood against some tolerance. \cr
#'
#' When the \code{"aitken"} method is employed, the final estimate of the log-likelihood is the value of \code{linf} at convergence, rather than the value of \code{ll} at convergence under the \code{"relative"} option.
#' @return A list with the following components:
#' \itemize{
#' \item{\code{ll} - }{The most current estimate for the log-likelihood.}
#' \item{\code{linf} - }{The most current estimate of the final converged maxmised log-likelihood.}
#' \item{\code{a} - }{The Aitken acceleration value where \code{0 <= a <= 1}.}
#' }
#' @export
#' @author Keefe Murphy - \href{keefe.murphy@ucd.ie}{<keefe.murphy@ucd.ie>}
#' @references Boehning, D., Dietz, E., Schaub, R., Schlattmann, P. and Lindsay, B. G. (1994). The distribution of the likelihood ratio for mixtures of densities from the one-parameter exponential family. \emph{Annals of the Institute of Statistical Mathematics}, 46(2): 373-388.
#'
#' @seealso \code{\link{MoE_control}}
#' @examples
#' (a1 <- MoE_aitken(-c(449.61534, 442.84221, 436.58999)))
#' a2  <- MoE_aitken(-c(442.84221, 436.58999, 436.58998))
#' abs(a2$linf - a1$linf) < 1e-05 #FALSE
#' a3  <- MoE_aitken(-c(436.58998, 436.58997, 436.58997))
#' abs(a3$linf - a2$linf) < 1e-05 #TRUE
#' (ll <- a3$linf)
  MoE_aitken      <- function(loglik) {
    l1            <- loglik[1]
    l2            <- loglik[2]
    l3            <- loglik[3]
    if(any(is.infinite(loglik))) {
      linf        <- Inf
      a           <- NA
    } else {
      a           <- ifelse(l2 > l1, (l3 - l2) / (l2 - l1),     0)
      denom       <- max(1 - a,  .Machine$double.eps)
      linf        <- ifelse(a  < 1,   l3 + (l3 - l2) / denom, Inf)
    }
      list(ll = l3, linf = linf, a = a)
  }

#' Choose the best MoEClust model
#'
#' Takes one or more sets of MoEClust models fitted by \code{\link{MoE_clust}} and ranks them according to the BIC, ICL, or AIC. It's possible to respect the internal ranking within each set of models, or to discard models within each set which were already deemed sub-optimal.
#' @param ... One or more objects of class \code{"MoEClust"} outputted by \code{\link{MoE_clust}}. All models must have been fit to the same data set. A single \emph{named} list of such objects can also be supplied. This argument is only relevant for the \code{\link{MoE_compare}} function and will be ignored for the associated \code{print} function.
#' @param criterion The criterion used to determine the ranking. Defaults to \code{"bic"}.
#' @param pick The (integer) number of models to be ranked and compared. Defaults to \code{3L}. Will be constrained by the number of models within the \code{"MoEClust"} objects supplied via \code{...} if \code{optimal.only} is \code{FALSE}, otherwise constrained simply by the number of \code{"MoEClust"} objects supplied.
#' @param optimal.only Logical indicating whether to only rank models already deemed optimal within each \code{"MoEClust"} object (\code{TRUE}), or to allow models which were deemed suboptimal enter the final ranking (\code{FALSE}, the default). See \code{details}
#' @param x An object of class \code{"MoECompare"} resulting from a call to \code{\link{MoE_compare}}.
#' @param index A logical or numeric vector giving the indices of the rows of the table of ranked models to print. This defaults to the full set of ranked models. It can be useful when the table of ranked models is large to examine a subset via this \code{index} argument, for display purposes.
#'
#' @note The \code{criterion} argument here need not comply with the criterion used for model selection within each \code{"MoEClust"} object, but be aware that a mismatch in terms of \code{criterion} \emph{may} require the optimal model to be re-fit in order to be extracted, thereby slowing down \code{\link{MoE_compare}}.\cr
#'
#' A dedicated \code{print} function exists for objects of class \code{"MoECompare"}.
#'
#' @details The purpose of this function is to conduct model selection on \code{"MoEClust"} objects, fit to the same data set, with different combinations of gating/expert network covariates or different initialisation settings. \cr
#'
#' Model selection will have already been performed in terms of choosing the optimal number of components and \pkg{mclust} model type within each supplied set of results, but \code{\link{MoE_compare}} will respect the internal ranking of models when producing the final ranking if \code{optimal.only} is \code{FALSE}: otherwise only those models already deemed optimal within each \code{"MoEClust"} object will be ranked.\cr
#'
#' As such if two sets of results are supplied when \code{optimal.only} is \code{FALSE}, the 1st, 2nd and 3rd best models could all belong to the first set of results, meaning a model deemed suboptimal according to one set of covariates could be superior to one deemed optimal under another set of covariates.
#' @return A list of class \code{"MoE_compare"}, for which a dedicated print function exists, containing the following elements, each of \code{length(pick)}:
#' \describe{
#' \item{\code{optimal}}{The single optimal model (an object of class \code{"MoEClust"}) among those supplied, according to the chosen \code{criterion}.}
#' \item{\code{pick}}{The final number of ranked models. May be different (i.e. less than) the supplied \code{pick} value.}
#' \item{\code{MoENames}}{The names of the supplied \code{"MoEClust"} objects.}
#' \item{\code{modelNames}}{The \code{\link[mclust]{mclustModelNames}}.}
#' \item{\code{G}}{The optimal numbers of components.}
#' \item{\code{bic}}{The ranked BIC values.}
#' \item{\code{icl}}{The ranked ICL values.}
#' \item{\code{aic}}{The ranked AIC values.}
#' \item{\code{gating}}{The gating formulas.}
#' \item{\code{expert}}{The expert formulas.}
#' \item{\code{equalPro}}{Logical indicating whether mixing proportions were constrained to be equal across components.}
#' }
#' @export
#' @references K. Murphy and T. B. Murphy (2017). Parsimonious Model-Based Clustering with Gating and Expert Network Covariates.
#' @note \code{\link{plot.MoEClust}} and \code{\link{as.Mclust}} can both also be called on objects of class \code{"MoECompare"}.
#' @importFrom mclust "mclustModelNames"
#' @author Keefe Murphy - \href{keefe.murphy@ucd.ie}{<keefe.murphy@ucd.ie>}
#'
#' @seealso \code{\link{MoE_clust}}, \code{\link[mclust]{mclustModelNames}}, \code{\link{plot.MoEClust}}, \code{\link{as.Mclust}}
#' @examples
#' data(CO2data)
#' GNP   <- CO2data[,1]
#' CO2   <- CO2data[,2]
#' m1    <- MoE_clust(CO2, G=1:2)
#' m2    <- MoE_clust(CO2, G=1:2, gating= ~ GNP)
#' m3    <- MoE_clust(CO2, G=1:2, expert= ~ GNP)
#' m4    <- MoE_clust(CO2, G=1:2, gating= ~ GNP, expert= ~ GNP)
#' m5    <- MoE_clust(CO2, G=1:2, equalPro=TRUE)
#' m6    <- MoE_clust(CO2, G=1:2, expert= ~ GNP, equalPro=TRUE)
#'
#' # Rank only the optimal models and examine the best model
#' (comp <- MoE_compare(m1, m2, m3, m4, m5, m6, pick=6, optimal.only=TRUE))
#' (best <- comp$optimal)
#' (summ <- summary(best))
#'
#' # Examine all models visited, including those already deemed suboptimal
#' # Only print models with expert covariates & more than one component
#' comp2 <- MoE_compare(m1, m2, m3, m4, m5, m6, pick=18)
#' print(comp2, comp2$expert != "None" & comp2$G > 1)
  MoE_compare     <- function(..., criterion = c("bic", "icl", "aic"), pick = 3L, optimal.only = FALSE) {
    crit.miss     <- missing(criterion)
    if(!missing(criterion) && (length(criterion) > 1 ||
       !is.character(criterion)))                 stop("'criterion' must be a single character string")
    criterion     <- match.arg(criterion)
    num.miss      <- missing(pick)
    opt.miss      <- missing(optimal.only)
    if(length(pick)    != 1            ||
       !is.numeric(pick))                         stop("'pick' must be a single number")
    if(floor(pick)     != pick         ||
       pick        < 1)                           stop("'pick' must be a strictly positive integer")
    if(length(optimal.only)    > 1     ||
       !is.logical(optimal.only))                 stop("'optimal.only' must be a single logical indicator")
    call          <- match.call(expand.dots=TRUE)[-1]
    call          <- if(crit.miss) call else call[-which(names(call) == "criterion")]
    call          <- if(num.miss)  call else call[-which(names(call) == "pick")]
    call          <- if(opt.miss)  call else call[-which(names(call) == "optimal.only")]
    len.call      <- length(as.list(call))
    if(len.call   == 1 && is.list(...) && !inherits(..., "MoEClust")) {
      mod.names   <- names(...)
      MoEs        <- unique(...)
      if(is.null(mod.names))                      stop("When supplying models as a list, every element of the list must be named")
    } else {
      mod.names   <- unique(sapply(call, deparse))
      MoEs        <- unique(list(...))
    }
    MoEs          <- stats::setNames(MoEs, mod.names)
    if(any(sapply(MoEs, class) != "MoEClust"))    stop("All models must be of class 'MoE_clust'!")
    if(length(unique(sapply(MoEs,
       function(mod) mod$call$data)))  != 1)      stop("All models being compared must have been fit to the same data set!")
    title         <- "Comparison of Gaussian finite mixture of experts models fitted by EM algorithm"
    dat.name      <- deparse(MoEs[[1]]$call$data)
    gating        <- lapply(lapply(MoEs, "[[", "gating"), attr, "Formula")
    expert        <- lapply(lapply(MoEs, "[[", "expert"), attr, "Formula")
    hypvol        <- sapply(MoEs, "[[", "hypvol")
    BICs          <- lapply(MoEs, "[[", "BIC")
    ICLs          <- lapply(MoEs, "[[", "ICL")
    AICs          <- lapply(MoEs, "[[", "AIC")
    DFxs          <- lapply(MoEs, "[[", "DF")
    choice        <- max(lengths(switch(criterion, bic=BICs, icl=ICLs, aic=AICs)))
    bics          <- lapply(BICs, function(x) .pick_MoECrit(x, choice)$crits)
    icls          <- lapply(ICLs, function(x) .pick_MoECrit(x, choice)$crits)
    aics          <- lapply(AICs, function(x) .pick_MoECrit(x, choice)$crits)
    dfxs          <- lapply(DFxs, function(x) .pick_MoECrit(x, choice)$crits)
    if(optimal.only) {
      opt.names   <- names(.crits_names(lapply(switch(criterion, bic=bics, icl=icls, aic=aics), "[", 1)))
    }
    bics          <- .crits_names(bics)
    icls          <- .crits_names(icls)
    aics          <- .crits_names(aics)
    dfxs          <- .crits_names(dfxs)
    if(optimal.only) {
      bics        <- bics[names(bics) %in% opt.names]
      icls        <- icls[names(icls) %in% opt.names]
      aics        <- aics[names(aics) %in% opt.names]
      dfxs        <- dfxs[names(dfxs) %in% opt.names]
    }
    crits         <- switch(criterion, bic=bics, icl=icls, aic=aics)
    pick          <- min(pick, length(crits))
    max.crits     <- sort(crits, decreasing=TRUE)[seq_len(pick)]
    max.names     <- names(max.crits)
    crit.names    <- gsub("\\|.*", "",          max.names)
    G             <- as.numeric(gsub(".*,", "", max.names))
    equalPro      <- sapply(MoEs, attr, "EqualPro")
    gating        <- unname(unlist(gating[crit.names]))
    expert        <- unname(unlist(expert[crit.names]))
    modelNames    <- gsub(",.*", "", gsub(".*\\|", "", max.names))
    best.model    <- MoEs[[crit.names[1]]]
    if(best.model$modelName    != modelNames[1] || best.model$G != G[1]) {
      best.model$net.covs      <- if(best.model$G == 1) attr(best.model$net.covs, "Discarded") else best.model$net.covs
      cat("Re-fitting optimal model due to mismatched 'criterion'...\n\n")
      old.call    <- best.model$call
      old.call    <- c(as.list(old.call)[1], list(criterion=criterion), as.list(old.call)[-1])
      old.call    <- as.call(old.call[!duplicated(names(old.call))])
      best.call   <- c(list(data=best.model$data, modelNames=modelNames[1], G=G[1], verbose=FALSE, network.data=best.model$net.covs), as.list(old.call[-1]))
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
        best.model$df                  <- best.mod$df
        best.model$hypvol              <- best.mod$hypvol
        best.model$parameters          <- best.mod$parameters
        best.model$z                   <- best.mod$z
        best.model$classification      <- best.mod$classification
        attributes(best.model)         <- attributes(best.mod)
      } else best.model                <- paste0("Failed to re-fit the optimal model: ", gsub("\"", "'", deparse(old.call, width.cutoff=500L), fixed=TRUE))
    }
    gating[gating == "~1" | G  == 1]   <- "None"
    expert[expert == "~1"]             <- "None"
    comp          <- list(title = title, data = dat.name, optimal = best.model, pick = pick, MoENames = crit.names,
                          modelNames = modelNames, G = G, df = round(unname(dfxs[max.names]), 2), bic = round(unname(bics[max.names]), 2),
                          icl = round(unname(icls[max.names]), 2), aic = round(unname(aics[max.names]), 2), gating = gating, expert = expert,
                          equalPro = G == 1 | unname(equalPro[crit.names]), hypvol = unname(hypvol[crit.names]))
    class(comp)   <- c("MoECompare", "MoEClust")
      comp
  }

#' Convert MoEClust objects to the Mclust class
#'
#' Converts an object of class \code{"MoEClust"} generated by \code{\link{MoE_clust}} and converts it to an object of class \code{"Mclust"} as generated by fitting \code{\link[mclust]{Mclust}}, to facilitate use of plotting and other functions for the \code{"Mclust"} class within the \pkg{mclust} package.
#' @param x An object of class \code{"MoEClust"} generated by \code{\link{MoE_clust}} or an object of class \code{"MoECompare"} generated by \code{\link{MoE_compare}}.
#' @param resid Logical indicating whether to treat the data as the raw data (when \code{FALSE}, the default) or the augmented data comprising the residuals from the expert network (when \code{TRUE}). In the latter case, the mean and (co)variance parameters are taken to be the mean and (co)variance of the residuals. Only relevant if expert network covariates were supplied to \code{x}, otherwise coerced to \code{FALSE}.
#' @param signif Significance level for outlier removal. Must be a single number in the interval [0, 1). Corresponds to the percentage of data to be considered extreme and therefore removed (half of \code{signif} at each endpoint, on a column-wise basis). The default, \code{0}, corresponds to no outlier removal. \strong{Only} invoke this argument as an aid to visualisation via \code{\link[mclust]{plot.Mclust}}.
#' @param ... Further arguments to be passed to other methods.
#'
#' @return An object of class \code{"Mclust"}. See \code{methods(class="Mclust")} for a list of functions which can be applied to this class.
#' @details Of course, the user is always encouraged to use the dedicated \code{\link[=plot.MoEClust]{plot}} function for objects of the \code{"MoEClust"} class instead, but calling \code{plot} after converting via \code{\link{as.Mclust}} can be particularly useful for univariate mixtures.
#'
#' The \code{signif} argument is intended only to aid visualisation via \code{\link[mclust]{plot.Mclust}}, as plots therein can be sensitive to outliers, particularly with regard to axis limits. This is especially true when \code{resid} is \code{TRUE} in the presence of expert network covariates.
#' @note Of the functions which can be applied to the result of the conversion, \code{\link[mclust]{logLik.Mclust}} shouldn't be trusted in the presence of either expert network covariates, or (for more models with more than 1 component) gating network covariates.
#'
#' Mixing proportions are averaged over observations in components in the presence of gating network covariates during the coercion.
#'
#' Also note that plots may be misleading for models of univariate data with more than 1 component, in the presence of expert covariates when \code{resid} is \code{TRUE} and the \code{what} argument is either \code{"classification"} or \code{"uncertainty"} within \code{\link[mclust]{plot.Mclust}}.
#' @importFrom mclust "as.densityMclust.Mclust" "logLik.Mclust" "icl" "plot.Mclust" "plot.mclustBIC" "plot.mclustICL" "predict.Mclust" "print.Mclust" "summary.Mclust"
#' @export
#' @seealso \code{\link[mclust]{Mclust}}, \code{\link[mclust]{plot.Mclust}}, \code{\link{MoE_clust}}, \code{\link{plot.MoEClust}}
#' @author Keefe Murphy - \href{keefe.murphy@ucd.ie}{<keefe.murphy@ucd.ie>}
#' @references C. Fraley and A. E. Raftery (2002). Model-based clustering, discriminant analysis, and density estimation. \emph{Journal of the American Statistical Association}, 97:611-631.
#'
#' @examples
#' \dontrun{
#' # Fit a mixture of experts model to the ais data
#' data(ais)
#' mod <- MoE_clust(ais[,3:7], G=2, expert= ~ ais$sex, gating= ~ ais$BMI)
#'
#' # Convert to the "Mclust" class and examine the classification
#' plot(as.Mclust(mod), what="classification")
#'
#' # Examine the density using the augmented data in the expert network
#' plot(as.Mclust(mod, resid=TRUE), what="density")
#'
#' # While we could have just used plot.MoEClust above,
#' # plot.Mclust is especially useful for univariate data
#' data(CO2data)
#' res <- MoE_clust(CO2data$CO2, G=2, expert = ~ CO2data$GNP)
#' plot(as.Mclust(res))}
  as.Mclust       <- function (x, resid = FALSE, signif = 0, ...) {
    UseMethod("as.Mclust")
  }

#' @method as.Mclust MoEClust
#' @export
  as.Mclust.MoEClust      <- function(x, resid = FALSE, signif = 0, ...) {
    if(length(resid)  > 1 ||
       !is.logical(resid))                        stop("'resid' must be a single logical indicator")
    if(length(signif) > 1 || !is.numeric(signif) ||
       signif < 0 || signif >= 1)                 stop("'signif' must be a single number in the interval [0, 1)")
    x             <- if(inherits(x, "MoECompare")) x$optimal else x
    gating        <- attr(x, "Gating")
    resid         <- resid  && attr(x, "Expert")
    x$loglik      <- x$loglik[length(x$loglik)]
    x$BIC         <- replace(x$BIC, !is.finite(x$BIC), NA)
    class(x$BIC)  <- "mclustBIC"
    x$parameters$pro      <- if(gating) colMeans(x$z)                             else x$parameters$pro
    x$parameters$mean[]   <- if(resid)  0                                         else x$parameters$mean
    x$classification      <- if(resid)  unname(rep(x$classification, x$G))        else unname(x$classification)
    x$data                <- if(resid)  as.matrix(x$resid.data)                   else as.matrix(x$data)
    rownames(x$data)      <- if(resid)  seq_len(x$n * x$G)                        else rownames(x$data)
    x$data        <- if(signif > 0)     apply(x$data, 2, .trim_out, signif)       else x$data
    x$z           <- if(resid)          do.call(rbind, replicate(x$G, list(x$z))) else x$z # FIX THIS LINE
    dimnames(x$z) <- NULL
    x$uncertainty <- unname(x$uncertainty)
    x             <- x[-which(is.element(names(x), c("ICL", "icl", "AIC", "aic", "gating", "expert", "net.covs", "resid.data", "DF", "iters")))]
    name.x        <- names(x)
    attributes(x) <- NULL
    names(x)      <- name.x
    class(x)      <- "Mclust"
      x
  }

#' Quantile-Based Clustering for Univariate Data
#'
#' Returns a quantile-based clustering for univariate data.
#' @param x A vector of numeric data.
#' @param G The desired number of clusters.
#'
#' @return The vector of cluster labels.
#' @export
#'
#' @examples
#' data(CO2data)
#' MoE_qclass(CO2data$CO2, 2)
  MoE_qclass      <- function(x, G) {
    if((is.data.frame(x) || is.matrix(x)) &&
       (ncol(x)    > 1   || !all(is.numeric(x)))) stop("'x' must be univariate")
    x             <- as.vector(x)
    eps           <- stats::sd(x) * sqrt(.Machine$double.eps)
    q             <- NA
    n             <- G
    while(length(q) < (G + 1)) {
      n           <- n + 1
      q           <- unique(stats::quantile(x, seq(from=0, to=1, length=n)))
    }
    if(length(q)   > (G + 1))  {
      q           <- q[-order(diff(q))[seq_len(length(q) - G - 1)]]
    }
    q[1]          <- min(x) - eps
    q[length(q)]  <- max(x) + eps
    cl            <- rep(0, length(x))
    for(i in seq_len(G)) {
      cl[x >= q[i] & x < q[i + 1]] <- i
    }
      cl
  }

# Hidden/Print/Summary Functions
  .chol           <- function(x) tryCatch(chol(x), error=function(e) {
    d             <- nrow(x)
    eigs          <- eigen(x, symmetric = TRUE)
    eval          <- eigs$values
    evec          <- eigs$vectors
     return(chol(x + evec %*% tcrossprod(diag(pmax.int(0, 2 * max(abs(eval)) * d * .Machine$double.eps - eval), d), evec)))
    }
  )

  .crits_names    <- function(x) {
    unlist(lapply(seq_along(x), function(i) stats::setNames(x[[i]], paste0(names(x[i]), "|", names(x[[i]])))))
  }

  .drop_levels    <- function(fit, newdata) {
    factors       <- gsub("[[:space:]]", ".", gsub("[[:punct:]]", ".", gsub("[-^0-9]|as.factor|\\(|\\)", "", names(unlist(fit$xlevels)))))
    factorLevels  <- unname(unlist(fit$xlevels))
    modelFactors  <- cbind.data.frame(factors, factorLevels)
    predictors    <- names(newdata[names(newdata) %in% factors])
    for(i in seq_along(predictors))  {
      ind          <- newdata[,predictors[i]]      %in% modelFactors[modelFactors$factors == predictors[i],]$factorLevels
      if(any(!ind)) {
        newdata[!ind,predictors[i]] <-  NA
        newdata[,predictors[i]]     <- factor(newdata[,predictors[i]], levels=modelFactors[modelFactors$factors == predictors[i],]$factorLevels)
      }
    }
    newdata
  }

  .mat_byrow      <- function(x, nrow, ncol) {
    matrix(x, nrow=nrow, ncol=ncol, byrow=any(dim(as.matrix(x)) == 1))
  }

  .MoE_mahala     <- function(fit, resids, uni) {
    if(uni)  {
      resids * resids * 1/summary(fit)$sigma
    } else   {
      covar       <- stats::estVar(fit)
      inv.cov     <- try(base::solve(covar), silent=TRUE)
      if(inherits(inv.cov, "try-error"))        {
        covsvd    <- svd(covar)
        posi      <- covsvd$d > max(sqrt(.Machine$double.eps) * covsvd$d[1L], 0)
        inv.cov   <- if(all(posi)) covsvd$v %*% (t(covsvd$u)/covsvd$d) else if(any(posi))
          covsvd$v[,posi, drop=FALSE] %*% (t(covsvd$u[,posi, drop=FALSE])/covsvd$d[posi]) else array(0, dim(covar)[2L:1L])
      }
      rowSums(resids %*% inv.cov * resids)
    }
  }

#' @importFrom mclust "hypvol"
  .noise_vol      <- function(x, method = c("hypvol", "convexhull", "ellipsoidhull")) {
    ifelse(ncol(x) == 1, abs(diff(range(x))), switch(method, hypvol=hypvol(x, reciprocal = TRUE),
    convexhull=1/geometry::convhulln(x, options=c("Pp", "FA"))$vol, ellipsoidhull=1/cluster::volume(cluster::ellipsoidhull(x))))
  }

  .pick_MoECrit   <- function(x, pick = 3L) {
    if(!inherits(x, "MoECriterion"))              stop("'x' must be an object of class 'MoECriterion'")
    x             <- replace(x, !is.finite(x), NA)
    pick          <- min(pick,        length(x[!is.na(x)]))
    decrease      <- attr(x, "criterion") != "DF"
    x.sx          <- sort(x,          decreasing=decrease)[pick]
    x.crit        <- if(decrease)     x   >= x.sx else x <= x.sx
    x.ind         <- which(x.crit,    arr.ind=TRUE)
    x.val         <- sort(x[x.ind],   decreasing=decrease)
    ind.x         <- order(x[x.ind],  decreasing=decrease)
    x.ind         <- x.ind[ind.x,,    drop=FALSE]
    x.ind[,1]     <- gsub(".*= ", "", rownames(x)[x.ind[,1]])
    x.ind[,2]     <- colnames(x)[as.numeric(x.ind[,2])]
      return(list(crits = stats::setNames(x.val, vapply(seq_len(pick), function(p, b=x.ind[p,]) paste0(b[2], ",", b[1]), character(1L))), pick = pick))
  }

  .trim_out       <- function(x, signif = 0.01, na.rm = TRUE, ...) {
    qnt           <- stats::quantile(x, probs=c(signif, 2 - signif)/2, na.rm=na.rm, ...)
    H             <- 1.5    * stats::IQR(x, na.rm=na.rm)
    y             <- x
    li.qnt        <- qnt[1] - H
    ui.qnt        <- qnt[2] + H
    y[x < li.qnt] <- li.qnt
    y[x > ui.qnt] <- ui.qnt
    y
  }

#' @method print MoEClust
#' @importFrom mclust "mclustModelNames"
#' @rdname MoE_clust
#' @export
  print.MoEClust  <- function(x, digits = 2, ...) {
    cat("Call:\t");  print(x$call); cat("\n")
    if(length(digits)  > 1 || !is.numeric(digits) ||
       digits     <= 0)                           stop("Invalid 'digits'")
    name          <- x$modelName
    G             <- x$G
    equalP        <- G < 1 || attr(x, "EqualPro")
    gating        <- attr(x$gating, "Formula")
    expert        <- attr(x$expert, "Formula")
    gate.x        <- !attr(x, "Gating")
    exp.x         <- !attr(x, "Expert")
    net.x         <- !c(gate.x, exp.x)
    crit          <- round(unname(c(x$bic, x$icl, x$aic)), digits)
    Vinv          <- x$parameters$Vinv
    cat(paste0("Best Model: ",  mclustModelNames(name)$type, " (", name, "), with ",
               ifelse(G == 0, "only a noise component", paste0(G, " component", ifelse(G > 1, "s", ""))),
               ifelse(is.null(Vinv) || G == 0, "\n", " (and a noise component)\n"),
               ifelse(!equalP, "",   paste0("Equal Mixing Proportions\n")),
               ifelse(is.null(Vinv),  "", paste0("Hypervolume of Noise Component: ", round(Vinv, digits), "\n")),
               "BIC = ", crit[1], " | ICL = ", crit[2], " | AIC = ", crit[3],
               ifelse(any(net.x),    paste0("\nIncluding ", ifelse(all(net.x), "gating and expert", ifelse(!gate.x, "gating", ifelse(!exp.x, "expert", ""))), " network covariates:\n"), "\nNo covariates\n"),
               ifelse(gate.x,  "",   paste0("\tGating: ",   gating, ifelse(exp.x, "", "\n"))),
               ifelse(exp.x,   "",   paste0("\tExpert: ",   expert, ""))))
      invisible()
  }

#' @method summary MoEClust
#' @rdname MoE_clust
#' @export
  summary.MoEClust        <- function(object, ...) {
    G             <- object$G
    params        <- object$parameters
    equalPro      <- G == 1 || attr(object, "EqualPro")
    summ          <- list(title = "Gaussian finite mixture of experts model fitted by EM algorithm", data = deparse(object$call$data), n = object$n, d = object$d, G = G, modelName = object$modelName,
                          loglik = object$loglik[length(object$loglik)], df = object$df, gating = object$gating, expert = object$expert, bic=unname(object$bic), icl = unname(object$icl), aic = unname(object$aic),
                          pro = params$pro, mean = params$mean, variance = params$variance$sigma, Vinv = params$Vinv, z = object$z, equalPro = equalPro, classification = object$classification)
    class(summ)   <- "summary_MoEClust"
     summ
 }

#' @method print summary_MoEClust
#' @importFrom mclust "mclustModelNames"
#' @export
  print.summary_MoEClust  <- function(x, digits = 2, ...) {
    if(length(digits)  > 1 || !is.numeric(digits) ||
       digits     <= 0)                           stop("Invalid 'digits'")
    tmp           <- data.frame(log.likelihood = round(x$loglik, digits), n = x$n, d = x$d, df = x$df,
                                BIC = round(x$bic, digits), ICL = round(x$icl, digits), AIC = round(x$aic, digits))
    tmp           <- if(is.null(x$Vinv)) tmp else cbind(tmp, HypVol = 1/x$Vinv)
    rownames(tmp) <- NULL
    name          <- x$modelName
    G             <- x$G
    equalP        <- x$equalPro
    gating        <- attr(x$gating, "Formula")
    expert        <- attr(x$expert, "Formula")
    gate.x        <- gating == "~1"
    exp.x         <- expert == "~1"
    zs            <- stats::setNames(table(x$classification), NULL)
    cat(paste0("---------------------------------------------------------------\n", x$title, "\nData: ",
               x$data,"\n", "---------------------------------------------------------------\n\n",
               "MoEClust ",  name, " (", mclustModelNames(name)$type, "), with ",
               ifelse(G == 0, "only a noise component", paste0(G, " component", ifelse(G > 1, "s", ""))),
               ifelse(is.null(x$Vinv) || G == 0, "\n", " (and a noise component)\n"),
               ifelse(G > 1, paste0("\nEqual Mixing Proportions:  ", equalP && gate.x), ""),
               paste0("\nNoise Component:           ", !is.null(x$Vinv)),
               paste0("\nGating Network Covariates: ", ifelse(gate.x, "None", gating)),
               paste0("\nExpert Network Covariates: ", ifelse(exp.x,  "None", expert), "\n\n")))
    print(tmp, row.names = FALSE)
    cat("\nClustering table:")
    print(zs,  row.names = FALSE)
      invisible()
  }

#' @method print MoECompare
#' @rdname MoE_compare
#' @export
  print.MoECompare       <- function(x, index=seq_len(x$pick), ...) {
    index                <- if(is.logical(index)) which(index) else index
    if(length(index) < 1 || (!is.numeric(index) &&
       (any(index    < 1  | index > x$pick))))    stop("Invalid 'index'")
    x$noise              <- !is.na(x$hypvol)
    cat(paste0("------------------------------------------------------------------------------\n", x$title, "\nData: ",
               x$data,"\n", "------------------------------------------------------------------------------\n\n"))
    print(data.frame(do.call(cbind, x[-c(1:4, length(x) - 1)]))[index,])
      invisible()
  }

#' @method print MoECriterion
#' @export
  print.MoECriterion     <- function(x, pick = 3L, ...) {
    if(length(pick)      != 1    ||
       !is.numeric(pick))                         stop("'pick' must be a single number")
    if(floor(pick)       != pick ||
       pick        < 1)                           stop("'pick' be a strictly positive integer")
    crit          <- attr(x, "criterion")
    choice        <- .pick_MoECrit(x, pick)
    pick          <- choice$pick
    dim1          <- attr(x, "dim")
    dim2          <- attr(x, "dimnames")
    attributes(x) <- NULL
    attr(x, "dim")       <- dim1
    attr(x, "dimnames")  <- dim2
    cat(switch(EXPR= crit, BIC="Bayesian Information Criterion (BIC):\n", ICL="Integrated Completed Likelihood (ICL):\n",
                           AIC="Akaike Information Criterion (AIC):\n",    DF="Number of Estimated Parameters (DF):\n"))
    print(unclass(x))
    cat(paste0("\nTop ", ifelse(pick > 1, paste0(pick, " models"), "model"), " based on the ", crit, " criterion:\n"))
    print(choice$crits)
      invisible()
  }

#' @method print MoE_gating
#' @export
  print.MoE_gating       <- function(x, ...) {
    equalpro      <- attr(x, "EqualPro")
    formula       <- attr(x, "Formula")
    class(x)      <- class(x)[class(x) != "MoE_gating"]
    print(x, ...)
    cat(paste("EqualPro:", equalpro, "\n"))
    cat(paste("Formula:",  formula))
    message("\n\n\nUsers are cautioned against making inferences about statistical significance from summaries of the coefficients in the gating network")
      invisible(x)
  }

#' @method print MoE_expert
#' @export
  print.MoE_expert       <- function(x, ...) {
    formula       <- attr(x, "Formula")
    attributes(x)[-1]    <- NULL
    class(x)      <- "listof"
    print(x, ...)
    cat(paste("Formula:", formula))
    message("\n\n\nUsers are cautioned against making inferences about statistical significance from summaries of the coefficients in the expert network")
      invisible(x)
  }

#' @method summary MoE_gating
#' @export
  summary.MoE_gating     <- function(object, ...) {
    class(object) <- "multinom"
    summ          <- summary(object, ...)
    class(summ)   <- "summary_MoEgate"
      summ
  }

#' @method summary MoE_expert
#' @export
  summary.MoE_expert     <- function(object, clusters = seq_along(object), ...) {
    if(any(!is.numeric(clusters), any(clusters < 1),
           any(clusters   > length(object))))     stop("Invalid 'clusters'")
    summ          <- lapply(object[clusters], summary, ...)
    class(summ)   <- "summary_MoEexp"
      summ
  }

#' @method print summary_MoEgate
#' @export
  print.summary_MoEgate  <- function(x, ...) {
    class(x)      <- "summary.multinom"
    print(x, ...)
    message("\n\n\nUsers are cautioned against making inferences about statistical significance from summaries of the coefficients in the gating network")
      invisible(x)
  }

#' @method print summary_MoEexp
#' @export
  print.summary_MoEexp   <- function(x, ...) {
    class(x)      <- "listof"
    print(x, ...)
    message("Users are cautioned against making inferences about statistical significance from summaries of the coefficients in the expert network")
      invisible(x)
  }
