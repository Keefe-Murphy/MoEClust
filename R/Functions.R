#' Mixtures of Experts: Model-Based Clustering with Covariates
#'
#' Fits Mixture of Experts models with \pkg{mclust}-family covariance structures. In other words, performs model-based clustering via the EM algorithm where covariates are allowed to enter neither, either, or both the mixing proportions (gating network) and/or component densities (expert network).
#' @param data A numeric vector, matrix, or data frame of observations. Categorical variables are not allowed. If a matrix or data frame, rows correspond to observations and columns correspond to variables.
#' @param G An integer vector specifying the numbers of mixture components (clusters) to fit. Defaults to \code{G=1:9}.
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
#' The help file for \code{\link[mclust]{mclustModelNames}} further describes the available models (though the \code{"X"} in the single-component models will be coerced to \code{"E"} if supplied that way).
#' @param gating A formula for determining the model matrix for the multinomial logistic regression in the gating network when covariates enter the mixing proportions. This will be ignored where \code{G=1}. Interactions etc. are permitted. The specification of the LHS of the formula is ignored.
#' @param expert A formula for determining the model matrix for the multivariate WLS in the expert network when covariates are included in the component densities. Interactions etc. are permitted. The specification of the LHS of the formula is ignored.
#' @param network.data An optional data frame in which to look for the covariates in the \code{gating} &/or \code{expert} network formulas, if any. If not found in \code{network.data}, any supplied \code{gating} &/or \code{expert} covariates are taken from the environment from which \code{MoE_clust} is called.
#' @param control A list of control parameters for the EM and other aspects of the algorithm. The defaults are set by a call to \code{\link{MoE_control}}.
#' @param ... An alternative means of passing control parameters via the named arguments of \code{\link{MoE_control}}. Do not pass the output from a call to \code{\link{MoE_control}} here!
#'
#' @importFrom mclust "hc" "hclass" "hcVVV" "Mclust" "mclust.options" "mclustBIC" "mclustModelNames" "mclustVariance" "mstep" "mstepE" "mstepEEE" "mstepEEI" "mstepEEV" "mstepEII" "mstepEVE" "mstepEVI" "mstepEVV" "mstepV" "mstepVEE" "mstepVEI" "mstepVEV" "mstepVII" "mstepVVE" "mstepVVI" "mstepVVV" "nVarParams" "unmap"
#' @importFrom nnet "multinom"
#' @importFrom stats "as.formula" "binomial" "coef" "complete.cases" "glm" "kmeans" "lm" "model.frame" "predict" "residuals" "setNames" "update.formula"
#' @return A list (of class \code{"MoEClust"}) with the following named entries, mostly corresponding to the chosen optimal model (as determined by the \code{criterion} within \code{\link{MoE_control}}):\cr
#' \describe{
#' \item{\code{call}}{The matched call.}
#' \item{\code{data}}{The input data matrix.}
#' \item{\code{modelName}}{A character string denoting the \pkg{mclust} model type at which the optimal \code{criterion} occurs.}
#' \item{\code{n}}{The number of observations in the \code{data}.}
#' \item{\code{d}}{The dimension of the \code{data}.}
#' \item{\code{G}}{The optimal number of mixture components.}
#' \item{\code{BIC}}{A matrix of \emph{all} BIC values with \code{length{G}} rows and \code{length(modelNames)} columns. May include missing entries: \code{NA} represents models which were not visited, \code{-Inf} represents models which were terminated due to error, for which a log-likelihood could not be estimated.}
#' \item{\code{ICL}}{A matrix of \emph{all} ICL values with \code{length{G}} rows and \code{length(modelNames)} columns. May include missing entries: \code{NA} represents models which were not visited, \code{-Inf} represents models which were terminated due to error, for which a log-likelihood could not be estimated.}
#' \item{\code{AIC}}{A matrix of \emph{all} AIC values with \code{length{G}} rows and \code{length(modelNames)} columns. May include missing entries: \code{NA} represents models which were not visited, \code{-Inf} represents models which were terminated due to error, for which a log-likelihood could not be estimated.}
#' \item{\code{bic}}{The BIC value corresponding to the optimal model. May not necessarily be the optimal BIC.}
#' \item{\code{icl}}{The ICL value corresponding to the optimal model. May not necessarily be the optimal ICL.}
#' \item{\code{aic}}{The AIC value corresponding to the optimal model. May not necessarily be the optimal AIC.}
#' \item{\code{gating}}{The \code{\link[nnet]{multinom}} regression coefficients of the \code{gating} network. If \code{gating} covariates were \emph{NOT} supplied (or the best model has just one component), this corresponds to a RHS of ~1, otherwise the supplied \code{gating} formula. As such, a fitted \code{gating} network is always returned even in the absence of supplied covariates. The number of parameters to penalise by for \code{\link{MoE_crit}} is given by \code{gating$rank}, and the \code{gating} formula used is stored here as an attribute.}
#' \item{\code{expert}}{The multivariate WLS regression coefficients of the \code{expert} network. If \code{expert} covariates were NOT supplied, this corresponds to a RHS of ~1, otherwise the supplied \code{expert} formula. As such, a fitted \code{expert} network is always returned even in the absence of supplied covariates. The number of parameters to penalise by for \code{\link{MoE_crit}} is given by \code{G * expert[[1]]$rank}, and the \code{expert} formula used is stored here is an attribute.}
#' \item{\code{loglik}}{The vector of increasing log-likelihood values for every EM iteration under the optimal model.}
#' \item{\code{df}}{The number of estimated parameters in the optimal model.}
#' \item{\code{parameters}}{A list with the following components\cr
#' \itemize{
#' \item{\code{pro} - }{The mixing proportions: either a vector or, if \code{gating} covariates were supplied, a matrix with an entry for each observation (rows) and component (columns).}
#' \item{\code{mean} - }{The means of each component. If there are \code{expert} network covariates, this is the result of the extra M-step on the identified best model \emph{without} the \code{expert} covariates; as such it is just the cluster means according to the final allocation. The mean of the residuals used in the clustering accounting for the \code{expert} covariates can be found in \code{resid.mean}; if there are no \code{expert} covariates \code{mean} and \code{resid.mean} will be equal.}
#' \item{\code{resid.mean} - }{The mean of the residuals of each component used in the clustering accounting for the \code{expert} covariates; will be all \code{0} in the presence of \code{expert} network covariates, otherwise will be equal to \code{mean}, the cluster means according to the final allocation.}
#' \item{\code{variance} - }{A list of variance parameters of each component of the model. The components of this list depend on the model type specification. See the help file for \code{\link[mclust]{mclustVariance}} for details. If there are \code{expert} network covariates, this is the result of the extra M-step on the identified best model \emph{without} the \code{expert} covariates; as such it just the (co)variance of the clusters according to the final allocation. The (co)variance of the residuals used in the clustering accounting for the \code{expert} covariates can be found in \code{resid.variance}; if there are no \code{expert} covariates \code{variance} and \code{resid.variance} will be equal.}
#' \item{\code{resid.variance} - }{A list of variance parameters of the residuals of each component of the model used in the clustering accounting for the \code{expert} covariates. The components of this list depend on the model type specification. See the help file for \code{\link[mclust]{mclustVariance}} for details. Will be equal to \code{variance} in the absence of \code{expert} network covariates, otherwise will be the (less variable) residual (co)variance.}
#' }}
#' \item{\code{z}}{The final responsibility matrix whose \code{[i,k]}-th entry is the probability that observation \emph{i} belonds to the \emph{k}-th component.}
#' \item{\code{classification}}{The vector of cluster labels for the chosen model corresponding to \code{z}, i.e. \code{max.col(z)}.}
#' \item{\code{uncertainty}}{The uncertainty associated with the \code{classification}.}
#' \item{\code{net.covs}}{A data frame gathering the unique set of covariates used in the \code{gating} and \code{expert} networks, if any.}
#' \item{\code{DF}}{A matrix of \emph{all} Degrees of Freedom values with \code{length{G}} rows and \code{length(modelNames)} columns. May include missing entries: \code{NA} represents models which were not visited, \code{-Inf} represents models which were terminated due to error, for which parameters could not be estimated.}
#' \item{\code{iters}}{A matrix giving the total number of EM iterations with \code{length{G}} rows and \code{length(modelNames)} columns. May include missing entries: \code{NA} represents models which were not visited, \code{Inf} represents models which were terminated due to singularity/error and thus would never have converged.}
#' }
#' Dedicated \code{print} and \code{summary} functions exist for objects of class \code{"MoEClust"}.
#' @details The function effectively allows 4 different types of Mixture of Experts model (as well as the different models in the mclust family, for each): i) the standard finite Gaussian mixture, ii) covariates only in the gating network, iii) covariates only in the expert network, iv) the full Mixture of Experts model with covariates entering both the mixing proportions and component densities. Note that having the same covariates in both networks is allowed.\cr
#'
#' While model selection in terms of choosing the optimal number of components and the \pkg{mclust} model type is performed within \code{\link{MoE_clust}}, using one of the \code{criterion} options within \code{\link{MoE_control}}, choosing between multiple fits with different combinations of covariates or different initialisation settings can be done by supplying objects of class \code{"MoEClust"} to \code{\link{MoE_compare}}.
#' @note The EM algorithm finishes on an extra M-step once the best model has been identified.
#'
#' Where \code{BIC}, \code{ICL}, \code{AIC}, \code{DF} and \code{iters} contain \code{NA} entries, this corresponds to a model which was not run; for instance a VVV model is never run for single-component models as it is equivalent to EEE. As such, one can consider the value as not really missing, but equivalent to the EEE value. \code{BIC}, \code{ICL}, \code{AIC} and \code{DF} all inherit the class \code{"MoECriterion"}, for which a dedicated print function exists.
#' @seealso \code{\link{MoE_compare}}, \code{\link{MoE_control}}, \code{\link{MoE_crit}}, \code{\link{MoE_estep}}, \code{\link{MoE_dens}}, \code{\link[mclust]{mclustModelNames}}, \code{\link[mclust]{mclustVariance}}
#' @export
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
#' # Visualise the results using the 'lattice' library
#' # require("lattice")
#' # z   <- factor(best$classification, labels=paste0("Cluster", seq_len(best$G)))
#' # splom(~ hema | sex, groups=z)
#' # splom(~ hema | z, groups=sex)}
  MoE_clust       <- function(data, G = 1:9, modelNames = NULL, gating = NULL, expert = NULL, network.data = NULL, control = MoE_control(...), ...) {

  # Definitions and storage set-up
    call          <- match.call()
    multi         <- missing(modelNames)
    gate.x        <- !missing(gating)
    exp.x         <- !missing(expert)
    criterion     <- control$criterion
    equal.pro     <- control$equalPro
    init.z        <- control$init.z
    max.it        <- control$itmax[1]
    stopx         <- control$stopping
    tol           <- control$tol[1]
    verbose       <- control$verbose
    warnit        <- control$warn.it
    itwarn        <- warnit > 2
    control       <- control[-c(1:3, length(control), length(control) - 1)]
    if(!multi     &&
       !all(is.character(modelNames)))            stop("'modelNames' must be a vector of character strings")
    if(gate.x     &&
       !inherits(gating, "formula"))              stop("'gating' must be a formula")
    if(exp.x      &&
       !inherits(expert, "formula"))              stop("'expert' must be a formula")
    if(missing(data))                             stop("'data' must be supplied!")

    data          <- as.data.frame(data)
    num.X         <- vapply(data, is.numeric, logical(1L))
    if(anyNA(data))    {
      if(verbose)                                 message("Rows with missing values removed from data")
      data        <- data[complete.cases(data),, drop=FALSE]
    }
    if(sum(num.X) != ncol(data))    {
      if(verbose)                                 message("Non-numeric columns removed from data")
      data        <- data[,num.X, drop=FALSE]
    }
    X             <- as.matrix(data)
    n             <- nrow(X)
    d             <- ncol(X)
    if(!all(is.integer(G)) && any(G < 1))         stop("Invalid range of G values supplied")

    mod.fam       <- mclust.options("emModelNames")
    range.G       <- sort(unique(G))
    Gall          <- all(G  > 1)
    Gany          <- any(G  > 1)
    if((uni <- d  == 1))       {
      mfg         <- c("E", "V")
      mf1         <- "E"
    } else        {
      if(n   > d) {
        mfg       <- mod.fam
        mf1       <- c("EII", "EEI", "EEE")
      } else      {
        mfg       <- mod.fam[1:6]
        mf1       <- c("EII", "EEI")
      }
    }
    if(!multi)    {
      mNs         <- toupper(modelNames)
      sX          <- grepl("X",     mNs)
      if(any(sX))             {
       mNs        <- gsub("X", "E", mNs)
       if(verbose &&
          all(is.element(mNs,         mfg)))      message(paste0("'modelNames' which contain 'X' coerced to ", paste(shQuote(mNs[sX]), collapse=", ")))
      }
      if(Gany && any(!is.element(mNs, mfg)))      stop(paste0("Invalid 'modelNames'", ifelse(uni, " for univariate data", ifelse(n > d, "", " for high-dimensional data")), "!"))
      if(!Gall)   {
        sZ        <- !is.element(mNs, mf1)
        if(any(sZ))           {
          mf1     <- tryCatch(unname(vapply(mNs,  function(x)  switch(EXPR=x, E=, V="E", EII=, VII="EII", EEI=, VEI=, EVI=, VVI="EEI", EEE=, EVE=, VEE=,  VVE=, EEV=, VEV=, EVV=, VVV="EEE"), character(1L))),
                              error=function(e) { e$message <- paste0("Invalid 'modelNames' for single component models", ifelse(uni, " for univariate data", ifelse(n > d, "", " for high-dimensional data")), "!")
                                                  stop(e) } )
          if(verbose)                             message(paste0("'modelNames'", ifelse(any(sX), " further", ""), " coerced from ", paste(shQuote(mNs[sZ]), collapse=", "), " to ", paste(shQuote(mf1[sZ]), collapse=", "), " where 'G'=1"))
        }
      }
      mfg         <- mNs
    }
    mf1           <- unique(mf1)
    mfg           <- unique(mfg)
    all.mod       <- if(!Gall) unique(c(mf1, mfg)) else if(!Gany) mf1 else mfg
    multi         <- ifelse(Gall, length(unique(mfg)) > 1, ifelse(Gany, length(all.mod) > 1, length(unique(mf1)) > 1))
    BICs          <- ICLs     <-
    AICs          <- DF.x     <- it.x          <- provideDimnames(matrix(NA, nrow=length(range.G), ncol=length(all.mod)), base=list(paste0("G = ", as.character(range.G)), all.mod))
    crit.tx       <- crit.gx  <- -sqrt(.Machine$double.xmax)

  # Define the gating formula
    gate.G        <- rep(gate.x, length(range.G))
    if(gate.x)    {
      if(Gany)    {
        gating    <- tryCatch(update.formula(as.formula(gating), z ~ .),
                              error=function(e)   stop("Invalid 'gating' network formula supplied"))
        environment(gating)   <- environment()
        if(gating[[3]] == 1)   { if(verbose)      message("Not including gating network covariates with only intercept on gating formula RHS")
          gate.x  <- FALSE
          gate.G  <- rep(gate.x, length(range.G))
        }
      }
      if(!Gall)   {  if(verbose)                  message("Can't include gating network covariates in a single component mixture")
        gate.G[1] <- FALSE
      }
    } else gating <- as.formula(z ~ 1)
    if(equal.pro  && gate.x)   { if(verbose)      message("Can't constrain mixing proportions to be equal when gating covariates are supplied")
      equal.pro   <- FALSE
    }
    equal.pis     <- c(ifelse(!Gall, TRUE, equal.pro), rep(equal.pro, length(G) - 1))

  # Define the expert formula
    if(exp.x)     {
      expert      <- tryCatch(update.formula(as.formula(expert), X ~ .),
                              error=function(e)   stop("Invalid 'expert' network formula supplied"))
      environment(expert)     <- environment()
      if(expert[[3]]   == 1)   { if(verbose)      message("Not including expert network covariates with only intercept on expert formula RHS")
        exp.x     <- FALSE
      }
      Nseq        <- seq_len(n)
    } else expert <- as.formula(X ~ 1)
    if(init.z == "mclust"     &&
       !any(gate.x, exp.x))                       stop("Can't initialise using 'mclust' when there are no gating or expert covariates: try another 'init.z' method")

  # Tell network formulas where to look for variables
    if(!missing(network.data) &&
       !is.data.frame(network.data))              stop("'network.data' must be a data.frame if supplied")
    if(is.null(network.data))  {
     gate.covs    <- if(gate.x) model.frame(gating[-2]) else matrix(0, nrow=n, ncol=0)
     expx.covs    <- if(exp.x)  model.frame(expert[-2]) else matrix(0, nrow=n, ncol=0)
     network.data <- cbind(gate.covs, expx.covs)
     network.data <- as.data.frame(network.data[!duplicated(names(network.data))])
    }

  # Loop over range of G values and initialise allocations
    for(g in range.G)    {
      if(isTRUE(verbose))        cat(paste0("\n", g, " cluster model", ifelse(multi, "s", ""), " -\n"))
      x.dat       <- replicate(g, X, FALSE)
      h           <- which(range.G == g)
      equal.pro   <- equal.pis[h]
      gate.g      <- gate.G[h]
      Gseq        <- seq_len(g)
      z           <- z.init   <- unmap(if(g > 1) switch(init.z, hc=as.vector(hclass(hc(X, minclus=g), g)), kmeans=kmeans(X, g)$cluster,
                                       mclust=Mclust(X, g, verbose=FALSE)$classification, random=sample(Gseq, n, replace=TRUE)) else rep(1, n))

    # Initialise gating network
      lpi         <- matrix(0, nrow=n, ncol=g)
      if(gate.g)  {
        g.init    <- multinom(gating, trace=FALSE, data=network.data)
       #g.init    <- glmnet::cv.glmnet(y=z, x=model.matrix(gating)[,-1], family="multinomial", type.multinomial="grouped")
        pis       <- predict(g.init, type="probs")
       #pis       <- predict(g.init, type="response", newx=model.matrix(gating)[,-1], s="lambda.1se")[,,1]
        lpi       <- log(pis)
        gate.pen  <- g.init$rank
       #gate.pen  <- sum(lengths(coef(g.init, s="lambda.1se")[-1]))
      } else      {
        pis       <- if(equal.pro) rep(1/g, g) else 1
        lpi[,]    <- if(equal.pro)    log(pis) else 0
        gate.pen  <- ifelse(equal.pro, 0, g - 1)
      }
      if(exp.x)   {
        z.mat     <- matrix(0, nrow=n * g, ncol=g)
        muX       <- if(uni) rep(0, g) else matrix(0, nrow=d, ncol=g)
      } else expert.pen       <- g * d

    # Loop over the mclust model type(s)
      for(modtype in if(g == 1)  mf1   else mfg)  {
        m0W       <- m0X      <- FALSE

      # Initialise parameters from allocations
        if(isTRUE(verbose))      cat(paste0("\n\tModel: ", modtype, "\n"))
        x.df      <- nVarParams(modtype, d, g)  + gate.pen
        Mstep     <- mstep(modtype, X, z.init, control=control, equalPro=equal.pro)
        mus       <- Mstep$parameters$mean
        vari      <- Mstep$parameters$variance
        sigs      <- vari$sigma
        if(!gate.g)      {
          pis     <- if(equal.pro) pis  else Mstep$parameters$pro
          lpi     <- if(equal.pro) lpi  else matrix(log(pis), nrow=n, ncol=g, byrow=TRUE)
        }
        densme    <- capture.output(medensity  <- try(MoE_dens(modelName=modtype, data=x.dat, mus=mus, sigs=sigs, log.pis=lpi), silent=TRUE))
        if((ERR   <- attr(Mstep, "returnCode")  < 0 || inherits(medensity, "try-error"))) {
          ll      <- NA
          j       <- 1
          if(isTRUE(verbose))    cat(paste0("\t\t# Iterations: ", ifelse(ERR, "stopped at ", ""), j, "\n"))
          next
        } else    {
          z       <- MoE_estep(Dens=medensity)$z
          if((ERR <- any(is.nan(z))))             next
          ll      <- c(-Inf, -sqrt(.Machine$double.xmax))
          linf    <- rep(Inf, 2)
          j       <- 2
          stX     <- TRUE
        }

      # Run the EM algorithm
        while(stX)    {

        # Expert network
          if(exp.x)   {
            e.fit <- e.res    <- list()
            for(k in Gseq) {
             fitE <- lm(expert,  weights=z[,k], data=network.data)
            #fitE <- glmnet::cv.glmnet(y=X, x=model.matrix(expert), weights=z[,k])
             e.fit[[k]]       <- fitE
            #e.fit[[k]]       <- coef(fitE, s="lambda.1se")
             e.res[[k]]       <- residuals(fitE)
            #e.res[[k]]       <- X - predict(fitE, type="response", newx=model.matrix(expert), s="lambda.1se")[,,1]
             z.mat[(k - 1) * n + Nseq,k]       <- z[,k]
            }
            res.x <- if(uni) as.matrix(do.call(base::c, e.res)) else do.call(rbind, e.res)
          }
          x.df    <- ifelse(j == 2, ifelse(exp.x, g * e.fit[[1]]$rank, expert.pen) + x.df, x.df)

        # M-step
          Mstep   <- if(exp.x) mstep(modtype, res.x, z.mat, control=control) else mstep(modtype, X, z, control=control)
          mus     <- if(exp.x) muX else Mstep$parameters$mean
          vari    <- Mstep$parameters$variance
          sigs    <- vari$sigma

        # Gating Network
          if(gate.g)  {
            fitG  <- multinom(gating, trace=FALSE, data=network.data)
           #fitG  <- glmnet::cv.glmnet(y=z, x=model.matrix(gating)[,-1], family="multinomial", type.multinomial="grouped")
            pis   <- predict(fitG, type="probs")
           #pis   <- predict(fitG, type="response", newx=model.matrix(gating)[,-1], s="lambda.1se")[,,1]
            lpi   <- log(pis)
           #gate.pen          <- sum(lengths(coef(g.init, s="lambda.1se")[-1]))
          } else  {
            pis   <- if(equal.pro) pis else Mstep$parameters$pro
            pis   <- if(!exp.x)    pis else pis * g
            lpi   <- if(equal.pro) lpi else matrix(log(pis), nrow=n, ncol=g, byrow=TRUE)
          }

        # E-step & record log-likelihood
          densme  <- capture.output(medensity  <- try(MoE_dens(modelName=modtype, data=if(exp.x) e.res else x.dat, mus=mus, sigs=sigs, log.pis=lpi), silent=TRUE))
          if((ERR <- attr(Mstep, "returnCode")  < 0 || inherits(medensity, "try-error"))) {
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
              dX  <- abs((ll[j] - ll[j - 1])/(1 + ll[j]))
            }
            stX   <- dX >= tol && j < max.it
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
        if(isTRUE(verbose))      cat(paste0("\t\t# Iterations: ", ifelse(ERR, "stopped at ", ""), j, "\n"))
        ll[j]     <- switch(stopx, aitken=max(ll[j], ifelse(is.finite(linf[2]), linf[2], linf[1])), ll[j])
        choose    <- MoE_crit(modelName=modtype, loglik=ll[j], n=n, G=g, z=z, df=x.df)
        bics      <- choose["bic",]
        icls      <- choose["icl",]
        aics      <- choose["aic",]
        crit.t    <- switch(criterion, bic=bics, icl=icls, aic=aics)
        crit.t    <- ifelse(is.na(crit.t) || ERR, -Inf, crit.t)
        if(crit.t  > crit.tx)   {
          crit.tx <- crit.t
          pi.x    <- pis
          z.x     <- z
          ll.x    <- ll
          df.x    <- x.df
          if(gate.g)   {
            gfit  <- fitG
          }
          if(exp.x)    {
            efit  <- e.fit
            mu.x  <- mus
            sig.x <- vari
          }
        }
        BICs[h,modtype]       <- ifelse(ERR, -Inf, bics)
        ICLs[h,modtype]       <- ifelse(ERR, -Inf, icls)
        AICs[h,modtype]       <- ifelse(ERR, -Inf, aics)
        DF.x[h,modtype]       <- ifelse(ERR, -Inf, x.df)
        it.x[h,modtype]       <- ifelse(ERR,  Inf, j)
      } # for (modtype)

    # Pull out mclust model corresponding to highest BIC/ICL/AIC
      if(crit.tx   > crit.gx)   {
        crit.gx   <- crit.tx
        x.pis     <- pi.x
        x.z       <- z.x
        x.ll      <- ll.x
        x.df      <- df.x
        if(gate.g)     {
          x.fitG  <- gfit
        }
        if(exp.x)      {
          x.fitE  <- efit
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
    z             <- x.z
    x.ll          <- x.ll[!is.na(x.ll)][-c(1:2)]
    best.mod      <- colnames(CRITs)[best.ind[2]]
    bic.fin       <- setNames(BICs[best.ind], best.mod)
    icl.fin       <- setNames(ICLs[best.ind], best.mod)
    aic.fin       <- setNames(AICs[best.ind], best.mod)
    crit.fin      <- switch(criterion, bic=bic.fin, icl=icl.fin, aic=aic.fin)
    uncertainty   <- if(G > 1) 1 - apply(z, 1, max) else rep(0, n)
    bG            <- gate.G[which(range.G == G)]
    exp.gate      <- c(exp.x, bG)
    net.msg       <- ifelse(any(exp.gate), paste0(" (incl. ", ifelse(all(exp.gate), "gating and expert", ifelse(exp.x, "expert", ifelse(bG, "gating", ""))), " network covariates)"), "")

    x.fitG        <- if(bG)           x.fitG  else if(G > 1) multinom(gating, trace=FALSE, data=network.data) else glm(z ~ 1, family=binomial(link=logit))
    x.pis         <- if(G > 1 && !bG) x.fitG$fitted.values[1,]                                                else x.pis
    if(!exp.x)    {
     x.fitE       <- list()
     for(g in seq_len(G))  {
      x.fitE[[g]] <- lm(expert, weights=z[,g], data=network.data)
     }
    }
    x.fitE        <- setNames(x.fitE, paste0("Cluster", seq_len(G)))
    attr(x.fitG, "Equal.Pro") <- equal.pis[best.ind[1]]
    attr(x.fitG, "Formula")   <- Reduce(paste, deparse(gating[-2]))
    attr(x.fitE, "Formula")   <- Reduce(paste, deparse(expert[-2]))
    extraM        <- mstep(best.mod, X, z, control=control)
    mean.fin      <- extraM$parameters$mean
    variance.fin  <- extraM$parameters$variance
    resid.mu      <- if(exp.x) x.mu  else mean.fin
    resid.sig     <- if(exp.x) x.sig else variance.fin

    if(multi && verbose)         cat(paste0("\n\t\tBest Model: ", mclustModelNames(best.mod)$type, " (", best.mod, "), with ",
                                     G, " component", ifelse(G > 1, "s", ""), ifelse(any(exp.gate), net.msg, ""), "\n\t\t",
                                     switch(criterion, bic="BIC", icl="ICL", aic="AIC"), " = ", round(crit.fin, 2), "\n"))
    if(all(x.ll   != cummax(x.ll)))               warning("Log-likelihoods are not strictly increasing", call.=FALSE)
    if(any(it.x[!is.na(it.x)] == max.it))         warning(paste0("One or more models failed to converge in the maximum number of allowed iterations (", max.it, ")"), call.=FALSE)
    class(BICs)   <-
    class(ICLs)   <-
    class(AICs)   <-
    class(DF.x)   <- "MoECriterion"
    attr(BICs, "Criterion")   <- "BIC"
    attr(ICLs, "Criterion")   <- "ICL"
    attr(AICs, "Criterion")   <- "AIC"
    attr(DF.x, "Criterion")   <- "DF"
    results       <- list(call = call, data = X, modelName = best.mod, n = n, d = d, G = G,
                          BIC = BICs, ICL = ICLs, AIC = AICs, bic = bic.fin, icl = icl.fin, aic = aic.fin,
                          gating = x.fitG, expert = x.fitE, loglik = x.ll, df = x.df,
                          parameters = list(pro = x.pis, mean = mean.fin, resid.mean = resid.mu,
                          variance = variance.fin, resid.variance = resid.sig), z = z, classification = max.col(z), uncertainty = uncertainty,
                          net.covs = if(any(exp.x, gate.x)) network.data else "No covariates used in either network", DF = DF.x, iters = it.x)
    class(results)            <- "MoEClust"
    attr(results, "Details")  <- paste0(best.mod, ": ", G, " component", ifelse(G > 1, "s", ""), net.msg)
      return(results)
  }

#' Density for Parameterised MVN Mixture of Experts Models
#'
#' Computes densities (or log-densities) of observations in parameterised MVN mixtures of experts.
#' @param modelName A character string indicating the model. The help file for \code{\link[mclust]{mclustModelNames}} describes the available models.
#' @param data If there are no expert network covariates, \code{data} should be a numeric matrix or data frame, wherein rows correspond to observations (n) and columns correspond to variables (d). If there are expert network covariates, this should be a list of length G containing matrices/data.frames of multivariate WLS residuals for each component.
#' @param mus The mean for each of G components. If there is more than one component, this is a matrix whose k-th column is the mean of the k-th component of the mixture model. For the univariate models, this is a G-vector of means.
#' @param sigs A list of length G of variance parameters of the model. The components of this list depend on the specification of \code{modelName}.
#' @param log.pis If covariates enter the gating network, an n times G matrix of mixing proportions, otherwise a G-vector of mixing proportions for the components of the mixture. \strong{Must} be on the log-scale in both cases. The default of \code{0} effectively means densities (or log-densities) aren't scaled by the mixing proportions.
#' @param logarithm A logical value indicating whether or not the logarithm of the component densities should be returned. This defaults to \code{TRUE}, otherwise component densities are returned, obtained from the component log-densities by exponentiation. The \strong{log}-densities can be passed to \code{\link{MoE_estep}}.
#'
#' @note This function is intended for joint use with \code{\link{MoE_estep}}, using the \strong{log}-densities.
#' @importFrom mclust "mclustModelNames"
#' @importFrom mvnfast "dmvn"
#' @importFrom stats "dnorm"
#' @importFrom utils "capture.output"
#' @return A numeric matrix whose \code{[i,k]}-th entry is the density or log-density of observation \emph{i} in component \emph{k}, scaled by the mixing proportions.
#' @export
#' @author Keefe Murphy - \href{keefe.murphy@ucd.ie}{<keefe.murphy@ucd.ie>}
#'
#' @seealso \code{\link{MoE_estep}}, \code{\link{MoE_clust}}, \code{\link[mclust]{mclustModelNames}}
#'
#' @examples
#' data(ais)
#' hema  <- ais[,3:7]
#' model <- MoE_clust(hema, G=2, gating= ~ ais$BMI + ais$sex, model="EVE")
#' Dens  <- MoE_dens(modelName=model$modelName, data=hema,
#'                   mus=model$parameters$mean, sigs=model$parameters$variance$sigma,
#'                   log.pis=log(model$parameters$pro))
#'
#' # Construct the z matrix and compute the log-likelihood
#' Estep <- MoE_estep(Dens=Dens)
#' (ll   <- Estep$loglik)
#'
#' # The z matrix will be close but not exactly the same as that from the model
#' # as the EM algorithm finishes on an M-step, but the classification should be
#' identical(max.col(Estep$z), model$classification)    #TRUE
#' round(sum(Estep$z - model$z), options()$digits) == 0 #TRUE
  MoE_dens        <- function(modelName, data, mus, sigs, log.pis = 0, logarithm = TRUE) {
    G             <- ifelse(is.matrix(mus),  ncol(mus),  length(mus))
    Gseq          <- seq_len(G)
    if(!is.list(data) ||
        length(data)  != G) {
      data        <- replicate(G, as.matrix(data), FALSE)
    }
    dat1          <- data[[1]]
    n             <- ifelse(is.matrix(dat1), nrow(dat1), length(dat1))
    d             <- ifelse(is.matrix(dat1), ncol(dat1), 1)
    sq_mat        <- if(d > 50) function(x)  diag(sqrt(diag(x)))  else sqrt
    switch(EXPR=modelName, EVE=, VVE=, EEV=, EVV=, VVV = {
      idens       <- capture.output(densi      <- try(vapply(Gseq, function(k) dmvn(data[[k]], mus[,k], sigs[,,k],         log=TRUE, isChol=FALSE), numeric(n)), silent=TRUE))
    }, EEE=, VEE=, VEV =    {
      sigx        <- chol(sigs[,,1])  ;
      idens       <- capture.output(densi      <- try(vapply(Gseq, function(k) dmvn(data[[k]], mus[,k], sigx,              log=TRUE, isChol=TRUE),  numeric(n)), silent=TRUE))
    }, VII=, VEI=, VVI =    {
      idens       <- capture.output(densi      <- try(vapply(Gseq, function(k) dmvn(data[[k]], mus[,k], sq_mat(sigs[,,k]), log=TRUE, isChol=TRUE),  numeric(n)), silent=TRUE))
    }, EII=, EEI=, EVI =    {
      sigx        <- sq_mat(sigs[,,1]);
      idens       <- capture.output(densi      <- try(vapply(Gseq, function(k) dmvn(data[[k]], mus[,k], sigx,              log=TRUE, isChol=TRUE),  numeric(n)), silent=TRUE))
    }, E= {
      densi       <- vapply(Gseq, function(k)     dnorm(data[[k]], mus[k], sqrt(sigs), log=TRUE), numeric(n))
    }, V= {
      sigx        <- sqrt(sigs);
      densi       <- vapply(Gseq, function(k)     dnorm(data[[k]], mus[k], sigx[k],    log=TRUE), numeric(n))
    } )
    check         <- is.infinite(densi)         & densi > 0
    if(any(check))          {
      densi[which(check, arr.ind=TRUE)[1],]    <- rep(0, G)
    }
    densi         <- densi  + if(is.matrix(log.pis) || missing(log.pis)) log.pis else matrix(log.pis, nrow=n, ncol=G, byrow=TRUE)
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
#' model  <- MoE_clust(hema, G=2, gating= ~ ais$BMI + ais$sex, model="EVE")
#' Dens  <- MoE_dens(modelName=model$modelName, data=hema,
#'                   mus=model$parameters$mean, sigs=model$parameters$variance$sigma,
#'                   log.pis=log(model$parameters$pro))
#'
#' # Construct the z matrix and compute the log-likelihood
#' Estep  <- MoE_estep(Dens=Dens)
#' (ll    <- Estep$loglik)
#'
#' # The z matrix will be close but not exactly the same as that from the model
#' # as the EM algorithm finishes on an M-step, but the classification should be
#' identical(max.col(Estep$z), model$classification)    #TRUE
#' round(sum(Estep$z - model$z), options()$digits) == 0 #TRUE
#'
#' # Call MoE_estep directly
#' Estep2 <- MoE_estep(modelName=model$modelName, data=hema,
#'                     mus=model$parameters$mean, sigs=model$parameters$variance$sigma,
#'                     log.pis=log(model$parameters$pro))
#' identical(Estep2$loglik, ll)
  MoE_estep       <- function(modelName, data, mus, sigs, log.pis = 0, Dens = NULL) {
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
#' @param n The number of observations in the data used to compute \code{loglik}.
#' @param d The dimension of the data used to compute \code{loglik}.
#' @param G The number of components in the Gaussian mixture model used to compute \code{loglik}.
#' @param gating.pen The number of parameters of the \emph{gating} network of the MoEClust model. Defaults to \code{G - 1}, which corresponds to no gating covariates. If covariates are included, this should be the number of regression coefficients in the fitted object. If there are no covariates and mixing proportions are further assumed to be present in equal proportion, \code{gating.pen} should be \code{0}.
#' @param expert.pen The number of parameters of the \emph{expert} network of the MoEClust model. Defaults to \code{G * d}, which corresponds to no expert covariates. If covariates are included, this should be the number of regression coefficients in the fitted object.
#' @param z The \code{n} times \code{G} responsibility matrix whose \code{[i,k]}-th entry is the probability that observation \emph{i} belonds to the \emph{k}-th component.. If supplied the ICL is also computed and returned, otherwise only the BIC and AIC.
#' @param df An alternative way to specify the degrees of freedom exactly. If supplied, the arguments \code{d, gating.pen} and \code{expert.pen}, which are used to calculate the degrees of freedom, will be ignored.
#' @param delta Dirichlet hyperparameter for the prior on the mixing proportions. Defaults to 0.5. Only relevant for the ICL computation.
#'
#' @details The function is vectorized with respect to the arguments \code{modelName} and \code{loglik}.\cr
#'
#' If \code{model} is an object of class \code{"MoEClust"} with \code{G} components, the number of parameters for the \code{gating.pen} and \code{expert.pen} are \code{model$gating$rank} and \code{G * model$expert[[1]]$rank}, respectively.
#' @importFrom mclust "mclustModelNames" "nVarParams"
#' @return A simplified array containing the BIC, AIC, degrees of freedom and, if \code{z} is supplied, also the ICL, for each of the given input arguments.
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
#'                   expert.pen=G * model$expert[[1]]$rank)["bic",])
#' identical(bic2, unname(model$bic)) #TRUE
#'
#' # Make the same comparison with the known degrees of freedom
#' (bic3 <- MoE_crit(modelName=name, loglik=ll, n=n, G=G, df=model$df, z=z)["bic",])
#' identical(bic3, bic2)              #TRUE
  MoE_crit        <- Vectorize(function(modelName, loglik, n, d, G, gating.pen = G - 1, expert.pen = G * d, z = NULL, df = NULL, delta = 0.5) {
    df            <- ifelse(!missing(df), df, nVarParams(modelName, d, G) +  expert.pen + gating.pen)
    double.ll     <- 2 * loglik
    bic.x         <- double.ll  - df * log(n)
    aic.x         <- double.ll  - df * 2
      return(c(bic = bic.x, icl = if(!missing(z)) bic.x + 2 * sum((z == apply(z, 1, max)) * log(z), na.rm = TRUE), aic = aic.x, df = df))
  }, vectorize.args = c("modelName", "loglik"), SIMPLIFY="array")

#' Set control values for use with MoEClust
#'
#' Supplies a list of arguments (with defaults) for use with \code{\link{MoE_clust}}.
#' @param criterion When either \code{G} or \code{modelNames} is a vector, \code{criterion} determines whether the "\code{bic}" (Bayesian Information Criterion), "\code{icl}" (Integrated Complete Likelihood), "\code{aic}" (Akaike Information Criterion) is used to determine the 'best' model when gathering output. Note that all criteria will be returned in any case.
#' @param stopping The criterion used to assess convergence of the EM algorithm. The default (\code{"aitken"}) uses Aitken's acceleration method via \code{\link{MoE_aitken}}, otherwise the \code{"relative"} change in log-likelihood is monitored (which may be less strict). Both stopping rules are ultimately governed by \code{tol[1]}. When the \code{"aitken"} method is employed, the final estimate of the log-likelihood is the value of \code{linf} at convergence, rather than the value of \code{ll} at convergence under the \code{"relative"} option.
#' @param init.z The method used to initialise the cluster labels. Defaults to a hierarchical clustering tree as per \code{\link[mclust]{hc}}. Other options include \code{kmeans}, \code{random} initialisation, and a full run of \code{\link[mclust]{Mclust}}, although this last option is only permitted if there are \code{gating} &/or \code{expert} covariates within \code{\link{MoE_clust}}.
#' @param eps A scalar tolerance associated with deciding when to terminate computations due to computational singularity in covariances. Smaller values of eps allow computations to proceed nearer to singularity. The default is the relative machine precision \code{.Machine$double.eps}, which is approximately \emph{2e-16} on IEEE-compliant machines.
#' @param tol A vector of length two giving relative convergence tolerances for the log-likelihood and for parameter convergence in the inner loop for models with iterative M-step ("VEI", "EVE", "VEE", "VVE", "VEV"), respectively. The default is \code{c(1e-05,sqrt(.Machine$double.eps))}. If only one number is supplied, it is used as the tolerance for the outer iterations and the tolerance for the inner iterations is as in the default.
#' @param itmax A vector of length two giving integer limits on the number of EM iterations and on the number of iterations in the inner loop for models with iterative M-step ("VEI", "EVE", "VEE", "VVE", "VEV"), respectively. The default is \code{c(.Machine$integer.max, .Machine$integer.max)} allowing termination to be completely governed by \code{tol}. If only one number is supplied, it is used as the iteration limit for the outer iteration only.
#' @param equalPro Logical variable indicating whether or not the mixing proportions are equal in the model. Default: \code{equalPro = FALSE}. Only relevant when \code{gating} covariates are \emph{not} supplied within \code{\link{MoE_clust}}, otherwise ignored.
#' @param warn.it A single number giving the iteration count at which a warning will be printed if the EM algorithm has failed to converge. Defaults to \code{0}, i.e. no warning (which is true for any \code{warn.it} value less than \code{3}), otherwise the message is printed regardless of the value of \code{verbose}. If non-zero, \code{warn.it} should be moderately large, but obviously less than \code{itmax[1]}. A warning will always be printed if one of more models fail to converge in \code{itmax[1]} iterations.
#' @param verbose Logical indicating whether to print messages pertaining to progress to the screen during fitting. By default is \code{TRUE} if the session is interactive, and \code{FALSE} otherwise. If \code{FALSE}, warnings and error messages will still be printed to the screen, but everything else will be suppressed.
#'
#' @details \code{\link{MoE_control}} is provided for assigning values and defaults within \code{\link{MoE_clust}}.\cr
#'
#' While the \code{criterion} argument controls the choice of the optimal number of components and \pkg{mclust} model type, \code{\link{MoE_compare}} is provided for choosing between fits with different combinations of covariates or different initialisation settings.
#' @importFrom mclust "hc"
#' @return A named list in which the names are the names of the arguments and the values are the values supplied to the arguments.
#' @export
#' @author Keefe Murphy - \href{keefe.murphy@ucd.ie}{<keefe.murphy@ucd.ie>}
#'
#' @seealso \code{\link{MoE_clust}}, \code{\link{MoE_aitken}}, \code{\link[mclust]{hc}}, \code{\link{MoE_compare}}
#'
#' @examples
#' ctrl <- MoE_control(criterion="icl", itmax=10000, warn.it=12)
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
  MoE_control     <- function(criterion = c("bic", "icl", "aic"), stopping = c("aitken", "relative"),
                              init.z = c("hc", "kmeans", "mclust", "random"), eps = .Machine$double.eps,
                              tol = c(1e-05, sqrt(.Machine$double.eps)), itmax = c(.Machine$integer.max,
                              .Machine$integer.max), equalPro = FALSE, warn.it = 0, verbose = interactive()) {
    if(!is.character(criterion))                  stop("'criterion' must be a character vector of length 1")
    criterion     <- match.arg(criterion)
    if(!is.character(stopping))                   stop("'stopping' must be a character vector of length 1")
    stopping      <- match.arg(stopping)
    if(!is.character(init.z))                     stop("'init.z' must be a character vector of length 1")
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
    if(length(verbose)  < 1 ||
       !is.logical(verbose))                      stop("'verbose' must be a single logical indicator")
    if(length(warn.it)  > 1 ||
       !is.numeric(warn.it))                      stop("'warn.it' must be a numeric vector of length 1")
      list(criterion = criterion, stopping = stopping, init.z = init.z, eps = eps, tol = tol,
           itmax = itmax, equalPro = equalPro, warn.it=warn.it, verbose = verbose)
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
#' @references Boehning, D., Dietz, E., Schaub, R., Schlattmann, and Lindsay, B. (1994, June). The distribution of the likelihood ratio for mixtures of densities from the one-parameter exponential family. \emph{Annals of the Institute of Statistical Mathematics}, 46(2): 373-388.
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
#' Takes one or more sets of MoEClust models fitted by \code{\link{MoE_clust}} and ranks them according to the BIC, ICL, or AIC, while respecting the internal ranking within each set of models.
#' @param ... One or more objects of class \code{"MoEClust"} outputted by \code{\link{MoE_clust}}. All models must have been fit to the same data set. A single list of such objects can also be supplied.
#' @param criterion The criterion used to determine the ranking. Defaults to \code{"bic"}.
#' @param pick The (integer) number of models to be ranked and compared. Defaults to \code{3L}. Will be constrained by the number of models within the \code{"MoEClust"} objects supplied via \code{...}.
#'
#' @note The \code{criterion} argument here need not comply with the criterion used for model selection within each \code{"MoEClust"} object, but be aware that a mismatch in terms of \code{criterion} \emph{may} require the optimal model to be re-fit in order to be extracted, thereby slowing down \code{\link{MoE_compare}}.
#' @details The purpose of this function is to conduct model selection on \code{"MoEClust"} objects, fit to the same data set, with different combinations of gating/expert network covariates or different initialisation settings. \cr
#'
#' Model selection will have already been performed in terms of choosing the optimal number of components and \pkg{mclust} model type within each supplied set of results, but \code{\link{MoE_compare}} respects the internal ranking of models.\cr
#'
#' As such if two sets of results are supplied, the 1st, 2nd and 3rd best models could all belong to the first set of results, meaning a model deemed suboptimal according to one set of covariates could be superior to one deemed optimal under another set of covariates.
#' @return A list of class \code{"MoE_compare"}, for which a dedicated print function exists, containing the following elements, each of \code{length(pick)}:
#' \describe{
#' \item{\code{optimal}}{The single optimal model (an object of class \code{"MoEClust"}) among those supplied, according to the chosen \code{criterion}.}
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
#' @importFrom mclust "mclustModelNames"
#' @author Keefe Murphy - \href{keefe.murphy@ucd.ie}{<keefe.murphy@ucd.ie>}
#'
#' @seealso \code{\link{MoE_clust}}, \code{\link[mclust]{mclustModelNames}}
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
#' (comp <- MoE_compare(m1, m2, m3, m4, m5, m6, pick=5))
#' (best <- comp$optimal)
#' (summ <- summary(best))
  MoE_compare     <- function(..., criterion = c("bic", "icl", "aic"), pick = 3L) {
    crit.miss     <- missing(criterion)
    if(!is.character(criterion))                  stop("'criterion' must be a character string")
    criterion     <- match.arg(criterion)
    num.miss      <- missing(pick)
    if(length(pick)    != 1            ||
       !is.numeric(pick))                         stop("'pick' must be a single number")
    if(floor(pick)     != pick         ||
       pick        < 1)                           stop("'pick' must be a strictly positive integer")
    call          <- match.call(expand.dots=TRUE)[-1]
    call          <- if(crit.miss) call else call[-which(names(call) == "criterion")]
    call          <- if(num.miss)  call else call[-which(names(call) == "pick")]
    len.call      <- length(as.list(call))
    if(len.call   == 1 && is.list(...) && !inherits(..., "MoEClust")) {
      mod.names   <- names(...)
      MoEs        <- unique(...)
    } else {
      mod.names   <- unique(sapply(call, deparse))
      MoEs        <- unique(list(...))
    }
    MoEs          <- setNames(MoEs, mod.names)
    if(any(sapply(MoEs, class) != "MoEClust"))    stop("All models must be of class 'MoE_clust'!")
    if(length(unique(sapply(MoEs,
       function(mod) mod$call$data)))  != 1)      stop("All models being compared must have been fit to the same data set!")
    title         <- "Comparison of Gaussian finite mixture of experts models fitted by EM algorithm"
    dat.name      <- deparse(MoEs[[1]]$call$data)
    gating        <- lapply(MoEs, "[[", "gating")
    equalPro      <- sapply(gating, attr, "Equal.Pro")
    gating        <- lapply(gating, attr, "Formula")
    expert        <- lapply(lapply(MoEs, "[[", "expert"), attr, "Formula")
    BICs          <- lapply(MoEs, "[[", "BIC")
    ICLs          <- lapply(MoEs, "[[", "ICL")
    AICs          <- lapply(MoEs, "[[", "AIC")
    DFxs          <- lapply(MoEs, "[[", "DF")
    choice        <- max(lengths(switch(criterion, bic=BICs, icl=ICLs, aic=AICs)))
    bics          <- lapply(BICs, function(x) .pick_MoECrit(x, choice)$crits)
    icls          <- lapply(ICLs, function(x) .pick_MoECrit(x, choice)$crits)
    aics          <- lapply(AICs, function(x) .pick_MoECrit(x, choice)$crits)
    dfxs          <- lapply(DFxs, function(x) .pick_MoECrit(x, choice)$crits)
    bics          <- unlist(lapply(seq_along(bics), function(y) setNames(bics[[y]], paste0(names(bics[y]), "|", names(bics[[y]])))))
    icls          <- unlist(lapply(seq_along(icls), function(y) setNames(icls[[y]], paste0(names(icls[y]), "|", names(icls[[y]])))))
    aics          <- unlist(lapply(seq_along(aics), function(y) setNames(aics[[y]], paste0(names(aics[y]), "|", names(aics[[y]])))))
    dfxs          <- unlist(lapply(seq_along(dfxs), function(y) setNames(dfxs[[y]], paste0(names(dfxs[y]), "|", names(dfxs[[y]])))))
    crits         <- switch(criterion, bic=bics, icl=icls, aic=aics)
    max.crits     <- sort(crits, decreasing=TRUE)[seq_len(min(pick, length(crits)))]
    max.names     <- names(max.crits)
    crit.names    <- gsub("\\|.*", "", max.names)
    G             <- as.numeric(gsub(".*,", "", max.names))
    gating        <- unname(unlist(gating[crit.names]))
    expert        <- unname(unlist(expert[crit.names]))
    modelNames    <- gsub(",.*", "", gsub(".*\\|", "", max.names))
    best.model    <- tryCatch(get(crit.names[1]), error=function(e) MoEs[[crit.names[1]]])
    if(best.model$modelName   != modelNames[1] || best.model$G != G[1]) {
      old.call    <- best.model$call
      old.call    <- c(as.list(old.call)[1], list(criterion=criterion, modelNames=modelNames[1], G=G[1]), as.list(old.call)[-1])
      old.call    <- as.call(old.call[!duplicated(names(old.call))])
      best.call   <- c(list(data=best.model$data, modelNames=modelNames[1], G=G[1], verbose=FALSE, network.data=best.model$net.covs), as.list(old.call[-1]))
      best.call   <- best.call[!duplicated(names(best.call))]
      best.model  <- try(do.call(MoE_clust, best.call), silent=TRUE)
      if(!inherits(best.model, "try-error")) {
        best.model$call                <- old.call
      } else best.model                <- paste0("Failed to re-fit the optimal model: ", gsub("\"", "'", deparse(old.call, width.cutoff=500L), fixed=TRUE))
    }
    gating[gating == "~1" | G == 1]    <- "None"
    expert[expert == "~1"]             <- "None"
    comp          <- list(title = title, data = dat.name, optimal = best.model, MoENames = crit.names, modelNames = modelNames, G = G,
                          df = round(unname(dfxs[max.names]), 2), bic = round(unname(bics[max.names]), 2), icl = round(unname(icls[max.names]), 2),
                          aic = round(unname(aics[max.names]), 2), gating = gating, expert = expert, equalPro = unname(equalPro[crit.names]))
    class(comp)   <- "MoECompare"
      comp
  }

# Hidden/Print/Summary Functions
  .pick_MoECrit   <- function(x, pick = 3L) {
    if(!inherits(x, "MoECriterion"))              stop("'x' must be an object of class 'MoECriterion'")
    pick          <- min(pick,        length(x[!is.na(x)]))
    decrease      <- attr(x, "Criterion") != "DF"
    x.sx          <- sort(x,          decreasing=decrease)[pick]
    x.crit        <- if(decrease)     x   >= x.sx else x <= x.sx
    x.ind         <- which(x.crit,    arr.ind=TRUE)
    x.val         <- sort(x[x.ind],   decreasing=decrease)
    ind.x         <- order(x[x.ind],  decreasing=decrease)
    x.ind         <- x.ind[ind.x,,    drop=FALSE]
    x.ind[,1]     <- gsub(".*= ", "", rownames(x)[x.ind[,1]])
    x.ind[,2]     <- colnames(x)[as.numeric(x.ind[,2])]
      return(list(crits = setNames(x.val, vapply(seq_len(pick), function(p, b=x.ind[p,]) paste0(b[2], ",", b[1]), character(1L))), pick = pick))
  }

#' @method print MoEClust
#' @importFrom mclust "mclustModelNames"
#' @export
  print.MoEClust  <- function(x, ...) {
    cat("Call:\t");  print(x$call); cat("\n")
    name          <- x$modelName
    G             <- x$G
    equalP        <- attr(x$gating, "Equal.Pro") && G > 1
    gating        <- attr(x$gating, "Formula")
    expert        <- attr(x$expert, "Formula")
    gate.x        <- gating == "~1"
    exp.x         <- expert == "~1"
    net.x         <- !c(gate.x, exp.x)
    crit          <- round(unname(c(x$bic, x$icl, x$aic)), 2)
    cat(paste0("Best Model: ", mclustModelNames(name)$type, " (", name, "), with ",
               G, " component", ifelse(G > 1, "s\n", "\n"), "BIC = ", crit[1], " | ICL = ", crit[2], " | AIC = ", crit[3],
               ifelse(!equalP, "", paste0("\nEqual Mixing Proportions")),
               ifelse(any(net.x),  paste0("\nIncluding ", ifelse(all(net.x), "gating and expert", ifelse(!gate.x, "gating", ifelse(!exp.x, "expert", ""))), " network covariates:\n"), "\nNo covariates\n"),
               ifelse(gate.x,  "", paste0("\tGating: ",   gating, ifelse(exp.x, "", "\n"))),
               ifelse(exp.x,   "", paste0("\tExpert: ",   expert, ""))))
      invisible()
  }

#' @method summary MoEClust
#' @export
  summary.MoEClust       <- function(object, ...) {
    title         <- "Gaussian finite mixture of experts model fitted by EM algorithm"
    dat.name      <- deparse(object$call$data)
    params        <- object$parameters
    gating        <- object$gating
    expert        <- object$expert
    equalP        <- attr(gating, "Equal.Pro")
    summ          <- list(title = title, data = dat.name, n = object$n, d = object$d, G = object$G,
                          modelName = object$modelName, loglik = object$loglik[length(object$loglik)],
                          df = object$df, gating = gating, expert = expert, bic=unname(object$bic), icl = unname(object$icl),
                          aic = unname(object$aic), pro = params$pro, mean = params$mean, resid.mean = params$resid.mean,
                          variance = params$variance$sigma, resid.variance = params$resid.variance$sigma,
                          z = object$z, equalPro = equalP, classification = object$classification)
    class(summ)   <- "summary_MoEClust"
     summ
 }

#' @method print summary_MoEClust
#' @importFrom mclust "mclustModelNames"
#' @export
  print.summary_MoEClust <- function(x, ...) {
    tmp           <- data.frame(log.likelihood = x$loglik, n = x$n, d = x$d, df = x$df, BIC = x$bic, ICL = x$icl, AIC = x$aic)
    rownames(tmp) <- NULL
    name          <- x$modelName
    G             <- x$G
    equalP        <- x$equalPro
    gating        <- attr(x$gating, "Formula")
    expert        <- attr(x$expert, "Formula")
    gate.x        <- gating == "~1"
    exp.x         <- expert == "~1"
    zs            <- setNames(table(x$classification), NULL)
    cat(paste0("---------------------------------------------------------------\n", x$title, "\nData: ",
               x$data,"\n", "---------------------------------------------------------------\n\n",
               "MoEClust ", name, " (", mclustModelNames(name)$type, "), with ", G, " component", ifelse(G > 1, "s", ""), ":\n",
               ifelse(G > 1, paste0("\nEqual Mixing Proportions:  ", equalP && gate.x), ""),
               paste0("\nGating Network Covariates: ", ifelse(gate.x, "None", gating)),
               paste0("\nExpert Network Covariates: ", ifelse(exp.x,  "None", expert), "\n\n")))
    print(tmp, row.names = FALSE)
    cat("\nClustering table:")
    print(zs,  row.names = FALSE)
      invisible()
  }

#' @method print MoECompare
#' @export
  print.MoECompare       <- function(x, ...) {
    cat(paste0("------------------------------------------------------------------------------\n", x$title, "\nData: ",
               x$data,"\n", "------------------------------------------------------------------------------\n\n"))
    print(data.frame(do.call(cbind, x[-c(1:3)])))
      invisible()
  }

#' @method print MoECriterion
#' @export
  print.MoECriterion     <- function(x, pick = 3L, ...) {
    if(length(pick)      != 1    ||
       !is.numeric(pick))                         stop("'pick' must be a single number")
    if(floor(pick)       != pick ||
       pick        < 1)                           stop("'pick' be a strictly positive integer")
    crit          <- attr(x, "Criterion")
    choice        <- .pick_MoECrit(x, pick)
    pick          <- choice$pick
    attr(x, "Criterion") <- NULL
    cat(switch(EXPR= crit, BIC="Bayesian Information Criterion (BIC):\n", ICL="Integrated Completed Likelihood (ICL):\n",
                           AIC="Akaike Information Criterion (AIC):\n",    DF="Degrees of Freedom (DF):\n"))
    print(unclass(x))
    cat(paste0("\nTop ", ifelse(pick > 1, paste0(pick, " models"), "model"), " based on the ", crit, " criterion:\n"))
    print(choice$crits)
      invisible()
  }
#
