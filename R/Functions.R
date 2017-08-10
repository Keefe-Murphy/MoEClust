#' Mixtures of Experts: Model-Based Clustering with Covariates
#'
#' Fits Mixture of Experts models with \pkg{mclust}-family covariance structures. In other words, performs model-based clustering via the EM algorithm where covariates are allowed to enter neither, either, or both the mixing proportions (gating network) and/or component densities (expert network).
#' @param data A numeric vector, matrix, or data frame of observations. Categorical variables are not allowed. If a matrix or data frame, rows correspond to observations and columns correspond to variables.
#' @param G An integer vector specifying the numbers of mixture components (clusters) to fit. Defaults to \code{G=1:9}.
#' @param modelNames A vector of character strings indicating the models to be fitted in the EM phase of clustering. The default is:\cr
#' \tabular{ll}{for univariate data \tab \code{c("E", "V")}\cr
#' for multivariate data \eqn{N > D}{N > D} \tab \code{mclust.options("emModelNames")}\cr
#' for high-dimensional multivariate data \eqn{N \leq D}{N <= D} \tab \code{c("EII", "VII", "EEI", "EVI", "VEI", "VVI")}
#' }
#' \cr
#' For single-component models these options reduce to:\cr
#' \tabular{ll}{
#' for univariate data \tab \code{"E"}\cr
#' for multivariate data \eqn{N > D}{N > D} \tab \code{c("EII", "EEI", "EEE")}\cr
#' for high-dimensional multivariate data \eqn{N \leq D}{N <= D}  \tab \code{c("EII", "EEI")}
#' }
#' @param gating A formula for determining the model matrix for the multinomial logistic regression in the gating network when covariates enter the mixing proportions. This will be ignored where \code{G=1}. Interactions etc. are permitted. The specification of the LHS of the formula is ignored.
#' @param expert A formula for determining the model matrix for the multivariate WLS in the expert network when covariates are included in the component densities. Interactions etc. are permitted. The specification of the LHS of the formula is ignored.
#' @param control A list of control parameters for the EM and other aspects of the algorithm. The defaults are set by a call to \code{\link{MoE_control}}. Note that the argument \code{equalPro} will be ignored if \code{gating} network covariates are included.
#'
#' @details The function effectively allows 4 different types of Mixture of Expert model (as well as the different models in the mclust family, for each): i) the standard finite Gaussian mixture, ii) covariates only in the gating network, iii) covariates only in the expert network, iv) the full Mixture of Experts model with covariates entering both the mixing proportions and component densities. Note that having the same covariates in both networks is allowed.
#' @importFrom mclust "hc" "hclass" "hcVVV" "Mclust" "mclust.options" "mclustBIC" "mstep" "mstepE" "mstepEEE" "mstepEEI" "mstepEEV" "mstepEII" "mstepEVE" "mstepEVI" "mstepEVV" "mstepV" "mstepVEE" "mstepVEI" "mstepVEV" "mstepVII" "mstepVVE" "mstepVVI" "mstepVVV" "unmap"
#' @importFrom nnet "multinom"
#' @importFrom stats "as.formula" "coef" "complete.cases" "kmeans" "lm" "predict" "residuals" "setNames" "update.formula"
#' @return A list (of class \code{"MoEClust"}) with the following named entries, mostly corresponding to the chosen 'best' model (as determined by the \code{criterion} within \code{\link{MoE_control}}):\cr
#' \itemize{
#' \item{\code{mus} - }{The means of each component.}
#' \item{\code{sigmas} - }{The variance parameters of each component.}
#' \item{\code{Z} - }{The final responsibility matrix of probabilities of component membership.}
#' \item{\code{pis} - }{The mixing proportions: either a vector or, if \code{gating} covariates were supplied, a matrix with an entry for each observation (rows) and component (columns).}
#' \item{\code{classes} - }{The vector of cluster labels for the chosen model.}
#' \item{\code{BICs} - }{BIC values for every visited model. May include \code{NA} values for models which were terminated.}
#' \item{\code{ICLs} - }{ICL values for every visited model. May include \code{NA} values for models which were terminated.}
#' \item{\code{iters} - }{Total number of EM iterations for every visited model. May include \code{NA} values for models which were terminated.}
#' \item{\code{best.model} - }{A list of details pertaining to the chosen model, incl. the \pkg{mclust} model type, the number of components, as well as "\code{bic}" and "\code{icl}" values.}
#' \item{\code{log.likes} - }{The vector of increasing log-likelihood values for every EM iteration under the chosen model.}
#' \item{\code{gating} - }{The \code{\link[nnet]{multinom}} regression coefficients of the \code{gating} network if \code{gating} covariates were supplied. The number of parameters to penalise by for \code{\link{MoE_crit}} is given by \code{length(gating)}.}
#' \item{\code{expert} - }{The multivariate WLS regression coefficients of the \code{expert} network if \code{expert} covariates were supplied. The number of parameters to penalise by for \code{\link{MoE_crit}} is given by \code{sum(lengths(expert))}.}
#' }
#' @seealso \code{\link{MoE_control}}, \code{\link{MoE_crit}}, \code{\link{MoE_estep}}, \code{\link{MoE_dens}}
#' @export
#' @author Keefe Murphy - \href{keefe.murphy@ucd.ie}{<keefe.murphy@ucd.ie>}
#'
#' @examples
#' \dontrun{
#' data(ais)
#' hema <- ais[,3:7]
#'
#' # Fit a standard finite mixture model
#' m1   <- MoE_clust(hema, G=2)
#'
#' # Allow covariates to enter the mixing proportions
#' m2   <- MoE_clust(hema, G=2, gating= ~ ais$sex + ais$BMI)
#'
#' # Allow covariates to enter the component densities
#' m3   <- MoE_clust(hema, G=2, expert= ~ ais$sex)
#'
#' # Allow covariates to enter both the gating & expert network
#' m4   <- MoE_clust(hema, G=2, gating= ~ ais$BMI, expert= ~ ais$sex)
#'
#' # Extract the model with highest BIC
#' BICs <- c(m1=m1$best.model$bic, m2=m2$best.model$bic, m3=m3$best.model$bic, m4=m4$best.model$bic)
#' best <- get(gsub("\\..*", "", names(BICs)[which.max(BICs)]))
#' best$best.model
#'
#' # Visualise the results using the 'lattice' library
#' # require("lattice")
#' # z  <- factor(best$classes, labels=paste0("Cluster", seq_len(best$best.model$G)))
#' # splom(~ hema | sex, groups=z)
#' # splom(~ hema | z, groups=sex)
#' }
  MoE_clust       <- function(data, G = 1:9, modelNames = NULL, gating = NULL, expert = NULL, control = MoE_control()) {

  # Definitions and storage set-up
    data          <- as.data.frame(data)
    num.X         <- vapply(data, is.numeric, logical(1L))
    verbose       <- control$verbose
    if(anyNA(data))    {
      if(verbose)                                message("Rows with missing values removed from data")
      data        <- data[complete.cases(data),, drop=FALSE]
    }
    if(sum(num.X) != ncol(data))    {
      if(verbose)                                message("Non-numeric columns removed from data")
      data        <- data[,num.X, drop=FALSE]
    }
    X             <- as.matrix(data)
    N             <- nrow(X)
    D             <- ncol(X)

    if(!all(is.integer(G)) && any(G < 1))        stop("Invalid range of G values supplied")
    criterion     <- control$criterion
    stopping      <- control$stopping
    init.z        <- control$init.z
    multi         <- missing(modelNames)
    tol           <- control$tol[1]
    max.it        <- control$itmax[1]
    equal.pro     <- control$equalPro
    warnit        <- control$warn.it
    itwarn        <- warnit > 2
    control       <- control[-c(1:3, length(control), length(control) - 1)]
    mod.fam       <- mclust.options("emModelNames")
    range.G       <- sort(unique(G))
    Gall          <- all(G  > 1)
    Gany          <- any(G  > 1)
    uni           <- D     == 1
    if(uni)       {
      mfg         <- c("E", "V")
      mf1         <- "E"
    } else        {
      if(N > D)   {
        mfg       <- mod.fam
        mf1       <- c("EII", "EEI", "EEE")
      } else      {
        mfg       <- mod.fam[1:6]
        mf1       <- c("EII", "EEI")
      }
    }
    if(!multi)    {
      if(!all(is.element(modelNames, mfg)))      stop(paste0("Invalid 'modelNames'", ifelse(uni, ifelse(Gall, " for univariate data", "for single cluster univariate data"),
                                                      ifelse(Gany, "", " for single cluster model")), ":\n'modelNames' must be ", ifelse(all(!Gany, uni), "", "one of "), paste(shQuote(if(Gany) mfg else mf1), collapse=", "), "\n"))
      mfg         <- modelNames
      if(!Gall)   {
       mf1        <- unique(vapply(modelNames,   function(x) switch(EXPR=x, E=, V="E", EII=, VII="EII", EEI=, VEI=, EVI=, VVI="EEI", "EEE"), character(1L)))
       if(verbose &&
          any(!is.element(mf1, modelNames)))     message(paste0("'modelNames' coerced to ", paste(shQuote(mf1), collapse=", "), " as 'G' contains 1"))
      }
    }
    all.mod       <- unique(c(mf1, mfg))
    multi         <- ifelse(Gall, length(unique(mfg)) > 1, ifelse(Gany, length(all.mod) > 1, length(unique(mf1)) > 1))
    BICs          <- ICLs     <- it.x         <- provideDimnames(matrix(NA, nrow=length(range.G), ncol=length(all.mod)), base=list(paste0("G = ", as.character(range.G)), all.mod))
    crit.tx       <- crit.gx  <- -sqrt(.Machine$double.xmax)

  # Define the gating formula
    gate.x        <- !missing(gating)
    gate.G        <- rep(gate.x, length(range.G))
    if(gate.x)    {
      if(Gany)    {
        gating    <- tryCatch(update.formula(as.formula(gating), Z ~ .),
                              error=function(e)  stop("Invalid 'gating' network formula supplied"))
        environment(gating)   <- environment()
        if(gating[[3]] == 1)   { if(verbose)     message("Not including gating network covariates with only intercept on gating formula RHS")
          gate.x  <- FALSE
          gate.G  <- rep(gate.x, length(range.G))
        }
      }
      if(!Gall)   {  if(verbose)                 message("Can't include gating network covariates in a single component mixture")
        gate.G[1] <- FALSE
      }
    }
    if(gate.x     && equal.pro)                  stop("Can't constrain mixing proportions to be equal when gating covariates are supplied")

  # Define the expert formula
    exp.x         <- !missing(expert)
    if(exp.x)     {
      expert      <- tryCatch(update.formula(as.formula(expert), X ~ .),
                              error=function(e)  stop("Invalid 'expert' network formula supplied"))
      environment(expert)     <- environment()
      if(expert[[3]]   == 1)   { if(verbose)     message("Not including expert network covariates with only intercept on expert formula RHS")
        exp.x     <- FALSE
      }
      Nseq        <- seq_len(N)
    }
    if(init.z == "mclust"     &&
       !any(gate.x, exp.x))                      stop("Can't initialise using 'mclust' when there are no gating or expert covariates: try another 'init.z' method")

  # Loop over range of G values and initialise allocations
    for(g in range.G)    {
      if(isTRUE(verbose))        cat(paste0("\n", g, " cluster model", ifelse(multi, "s", ""), " -\n"))
      x.dat       <- replicate(g, X, FALSE)
      h           <- which(range.G == g)
      gate.g      <- gate.G[h]
      Gseq        <- seq_len(g)
      classes     <- if(g > 1) switch(init.z, hc=as.vector(hclass(hc(X, minclus=g), g)), kmeans=kmeans(X, g)$cluster,
                     mclust=Mclust(X, g, verbose=FALSE)$classification, random=sample(Gseq, N, replace=TRUE)) else rep(1, N)
      Z           <- Z.init   <- unmap(classes)

    # Initialise gating network
      lpi         <- matrix(0, nrow=N, ncol=g)
      if(gate.g)  {
        g.init    <- multinom(gating, trace=FALSE)
       #g.init    <- glmnet::cv.glmnet(y=Z, x=model.matrix(gating)[,-1], family="multinomial", type.multinomial="grouped")
        pis       <- predict(g.init, type="probs")
       #pis       <- predict(g.init, type="response", newx=model.matrix(gating)[,-1], s="lambda.1se")[,,1]
        lpi       <- log(pis)
        gate.pen  <- length(coef(g.init))
       #gate.pen  <- sum(lengths(coef(g.init, s="lambda.1se")[-1]))
      } else      {
        pis       <- if(equal.pro) rep(1/g, g) else 1
        lpi[,]    <- if(equal.pro)    log(pis) else 0
        gate.pen  <- ifelse(equal.pro, 0, g - 1)
      }
      if(exp.x)   {
        Z.mat     <- matrix(0, nrow=N * g, ncol=g)
      } else expert.pen       <- g * D

    # Loop over the mclust model type(s)
      modelsG1    <- if(g == 1)  mf1    else mfg
      for(modtype in modelsG1)  {
        m0W       <- m0X      <- FALSE
        m1        <- modtype  == modelsG1[1]

      # Initialise parameters from allocations
        if(isTRUE(verbose))      cat(paste0("\n\tModel: ", modtype, "\n"))
        Mstep     <- mstep(modtype, X, Z.init, control=control, equalPro=equal.pro)
        mus       <- Mstep$parameters$mean
        sigs      <- Mstep$parameters$variance$sigma
        sigs      <- if(uni)       sigs else array(unlist(sigs), dim=c(D, D, g))
        if(!gate.g)      {
          pis     <- if(equal.pro) pis  else Mstep$parameters$pro
          lpi     <- if(equal.pro) lpi  else matrix(log(pis), nrow=N, ncol=g, byrow=TRUE)
        }
        densme    <- capture.output(medensity <- try(MoE_dens(modelName=modtype, data=x.dat, mus=mus, sigs=sigs, log.pis=lpi), silent=TRUE))
        FAIL      <- attr(Mstep, "returnCode") < 0 || inherits(medensity, "try-error")
        if(FAIL)  {
          ll      <- NA
          j       <- 1
          if(isTRUE(verbose))    cat(paste0("\t\t# Iterations: ", ifelse(FAIL, "stopped at ", ""), j, "\n"))
          next
        } else    {
          Z       <- MoE_estep(Dens=medensity)$Z
          ll      <- c(-Inf, -sqrt(.Machine$double.xmax))
          j       <- 2
          stopx   <- TRUE
        }

      # Run the EM algorithm
        while(stopx)  {

        # Expert network
          if(exp.x)   {
            e.fit <- e.res    <- list()
            for(k in Gseq) {
             fitE <- lm(expert,  weights=Z[,k])
            #fitE <- glmnet::cv.glmnet(y=X, x=model.matrix(expert), weights=Z[,k])
             e.fit[[k]]       <- coef(fitE)
            #e.fit[[k]]       <- coef(fitE, s="lambda.1se")
             e.res[[k]]       <- residuals(fitE)
            #e.res[[k]]       <- X - predict(fitE, type="response", newx=model.matrix(expert), s="lambda.1se")[,,1]
             Z.mat[(k - 1) * N + Nseq,k]       <- Z[,k]
            }
            res.x <- if(uni) as.matrix(do.call(c, e.res)) else do.call(rbind, e.res)
            expert.pen        <- ifelse(m1, sum(lengths(e.fit)), expert.pen)
          }

        # M step
          Mstep   <- if(exp.x) mstep(modtype, res.x, Z.mat, control=control) else mstep(modtype, X, Z, control=control)
          mus     <- Mstep$parameters$mean
          sigs    <- Mstep$parameters$variance$sigma

        # Gating Network
          if(gate.g)  {
            fitG  <- multinom(gating, trace=FALSE)
           #fitG  <- glmnet::cv.glmnet(y=Z, x=model.matrix(gating)[,-1], family="multinomial", type.multinomial="grouped")
            pis   <- predict(fitG, type="probs")
           #pis   <- predict(fitG, type="response", newx=model.matrix(gating)[,-1], s="lambda.1se")[,,1]
            lpi   <- log(pis)
           #gate.pen          <- sum(lengths(coef(g.init, s="lambda.1se")[-1]))
          } else  {
            pis   <- if(equal.pro) pis else Mstep$parameters$pro
            pis   <- if(!exp.x)    pis else pis * g
            lpi   <- if(equal.pro) lpi else matrix(log(pis), nrow=N, ncol=g, byrow=TRUE)
          }

        # E step & record log-likelihood
          densme  <- capture.output(medensity <- try(MoE_dens(modelName=modtype, data=if(exp.x) e.res else x.dat, mus=mus, sigs=sigs, log.pis=lpi), silent=TRUE))
          FAIL    <- attr(Mstep, "returnCode") < 0 || inherits(medensity, "try-error")
          if(FAIL)     {
            ll    <- c(ll, NA)
            break
          } else  {
            Estep <- MoE_estep(Dens=medensity)
            Z     <- Estep$Z
            ll    <- c(ll, Estep$loglik)
            j     <- j + 1
            stopx <- switch(stopping, relative=abs((ll[j] - ll[j - 1])/(1 + ll[j])) >= tol,
                                      aitken=MoE_aitken(ll[seq(j - 2, j, 1)])$diff  >= tol) && j < max.it
            if(itwarn && !m0X)  {
             m0W  <- ifelse(!m0X, warnit < j, m0X)
             if(m0W   && !m0X)  {                tryCatch(warning("WARNIT", call.=FALSE), warning=function(w) message(paste0("\tEM algorithm for the ", modtype, " model has yet to converge in 'warn.it'=", warnit, " iterations")))
              m0X <- TRUE
             }
            }
          }

        } # while (j)

      # Store values corresponding to the maximum BIC/ICL so far
        if(isTRUE(verbose))      cat(paste0("\t\t# Iterations: ", ifelse(FAIL, "stopped at ", ""), j, "\n"))
        classes   <- max.col(Z)
        choose    <- MoE_crit(modelName=modtype, loglik=ll[j], N=N, D=D, G=g, gating.pen=gate.pen, expert.pen=expert.pen, nn=tabulate(classes, nbins=g))
        bics      <- choose["bic",]
        icls      <- choose["icl",]
        j.x       <- ifelse(FAIL, NA, j)
        crit.t    <- switch(criterion, bic=bics, icl=icls)
        crit.t    <- ifelse(is.na(crit.t), -Inf, crit.t)
        if(crit.t  > crit.tx)   {
          bic.x   <- bics
          icl.x   <- icls
          crit.tx <- crit.t
          mu.x    <- mus
          pi.x    <- pis
          sig.x   <- sigs
          z.x     <- Z
          log.x   <- ll
          j.x     <- j
          class.x <- classes
          if(gate.g)     {
            gfit  <- coef(fitG)
          }
          if(exp.x)      {
            efit  <- e.fit
          }
        }
        BICs[h,modtype]       <- bics
        ICLs[h,modtype]       <- icls
        it.x[h,modtype]       <- j.x
      } # for (modtype)

    # Pull out mclust model corresponding to highest BIC/ICL
      if(all(is.na(BICs)))                       stop("All models failed!")
      crit.g      <- switch(criterion, bic=bic.x, icl=icl.x)
      if(crit.g    > crit.gx)   {
        x.bic     <- bic.x
        x.icl     <- icl.x
        crit.gx   <- crit.g
        x.mu      <- mu.x
        x.pi      <- pi.x
        x.sig     <- sig.x
        x.z       <- z.x
        x.log     <- log.x
        x.class   <- class.x
        x.G       <- g
        if(gate.g)     {
          x.fitg  <- gfit
        }
        if(exp.x)      {
          x.fite  <- efit
        }
      }
    } # for (g)

  # Gather results
    CRITs         <- switch(criterion, bic=BICs, icl=ICLs)
    best.mod      <- colnames(CRITs)[which(CRITs == crit.gx, arr.ind=TRUE)[2]]
    bic.fin       <- setNames(x.bic, best.mod)
    icl.fin       <- setNames(x.icl, best.mod)
    crit.fin      <- switch(criterion, bic=bic.fin, icl=icl.fin)
    best.gate     <- gate.G[which(range.G == x.G)]
    exp.gate      <- c(exp.x, best.gate)
    net.msg       <- ifelse(any(exp.gate), paste0(" (incl. ", ifelse(all(exp.gate), "gating and expert", ifelse(exp.x, "expert", ifelse(best.gate, "gating", ""))), " network covariates)"), "")
    best.mod      <- list(modelName=best.mod, modelDetails=paste0(best.mod, ": ", x.G, " component", ifelse(x.G > 1, "s ", ""), net.msg), bic = bic.fin, icl = icl.fin, G = x.G)
    if(multi && verbose)         cat(paste0("\n\t\tBest Model: ", names(crit.fin), ", with ",
                                     x.G, " component", ifelse(x.G > 1, "s", ""), net.msg, "\n\t\t",
                                     switch(criterion, bic="BIC", "ICL"), " = ", round(crit.fin, 2), "\n"))
    if(all(x.log  != cummax(x.log)))             warning("Log-likelihoods are not strictly increasing", call.=FALSE)
    if(any(it.x[!is.na(it.x)] == max.it))        warning(paste0("One or more models failed to converge in ", max.it, " iterations"), call.=FALSE)
    results       <- c(list(mus = x.mu, sigmas = x.sig, Z = x.z, pis = x.pi, classes = x.class, BICs = BICs,
                            ICLs = ICLs, best.model = best.mod, iters = it.x, log.likes = x.log[!is.na(x.log)][-c(1:2)]),
                       if(best.gate) list(gating = x.fitg), if(exp.x) list(expert = setNames(x.fite, paste0("Cluster", seq_len(x.G)))))
    class(results)            <- "MoEClust"
      return(results)
  }

#' Density for Parameterised MVN Mixture of Experts Models
#'
#' Computes densities (or log-densities) of observations in parameterised MVN mixtures of experts.
#' @param modelName A character string indicating the model. The help file for \code{\link[mclust]{mclustModelNames}} describes the available models.
#' @param data If there are no expert network covariates, \code{data} should be a numeric matrix or data frame, wherein rows correspond to observations (N) and columns correspond to variables (D). If there are expert network covariates, this should be a list of length G containing matrices/data.frames of multivariate WLS residuals for each component.
#' @param mus The mean for each of G components. If there is more than one component, this is a matrix whose k-th column is the mean of the k-th component of the mixture model. For the univariate models, this is a G-vector of means.
#' @param sigs a list of length G of variance parameters of the model. The components of this list depend on the specification of \code{modelName}.
#' @param log.pis If covariates enter the gating network, an N times G matrix of mixing proportions, otherwise a G-vector of mixing proportions for the components of the mixture. \strong{Must} be on the log-scale in both cases. The default of \code{0} effectively means densities (or log-densities) aren't scaled by the mixing proportions.
#' @param logarithm A logical value indicating whether or not the logarithm of the component densities should be returned. This defaults to \code{TRUE}, otherwise component densities are returned, obtained from the component log-densities by exponentiation. The \strong{log}-densities can be passed to \code{\link{MoE_estep}}.
#'
#' @note This function is intended for joint use with \code{\link{MoE_estep}}, using the \strong{log}-densities.
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
#' Dens  <- MoE_dens(modelName=model$best.mod$modelName, data=hema,
#'                   mus=model$mus, sigs=model$sigmas, log.pis=log(model$pis))
#'
#' # Construct the Z matrix and compute the log-likelihood
#' Estep <- MoE_estep(Dens=Dens)
#' identical(Estep$Z, model$Z)
#' Estep$loglik
  MoE_dens        <- function(modelName, data, mus, sigs, log.pis = 0, logarithm = TRUE) {
    G             <- ifelse(is.matrix(mus),       ncol(mus),       length(mus))
    Gseq          <- seq_len(G)
    if(!is.list(data) ||
        length(data)  != G) {
      data        <- replicate(G, as.matrix(data), FALSE)
    }
    dat1          <- data[[1]]
    N             <- ifelse(is.matrix(dat1), nrow(dat1), length(dat1))
    P             <- ifelse(is.matrix(dat1), ncol(dat1), 1)
    sq_mat        <- if(P > 50) function(x)  diag(sqrt(diag(x)))  else sqrt
    switch(EXPR=modelName, EVE=, VVE=, EEV=, EVV=, VVV = {
      idens       <- capture.output(densi     <- try(vapply(Gseq, function(k) dmvn(data[[k]], mus[,k], sigs[,,k],         log=TRUE, isChol=FALSE), numeric(N)), silent=TRUE))
    }, EEE=, VEE=, VEV =    {
      sigx        <- chol(sigs[,,1])  ;
      idens       <- capture.output(densi     <- try(vapply(Gseq, function(k) dmvn(data[[k]], mus[,k], sigx,              log=TRUE, isChol=TRUE),  numeric(N)), silent=TRUE))
    }, VII=, VEI=, VVI =    {
      idens       <- capture.output(densi     <- try(vapply(Gseq, function(k) dmvn(data[[k]], mus[,k], sq_mat(sigs[,,k]), log=TRUE, isChol=TRUE),  numeric(N)), silent=TRUE))
    }, EII=, EEI=, EVI =    {
      sigx        <- sq_mat(sigs[,,1]);
      idens       <- capture.output(densi     <- try(vapply(Gseq, function(k) dmvn(data[[k]], mus[,k], sigx,              log=TRUE, isChol=TRUE),  numeric(N)), silent=TRUE))
    }, E= {
      densi       <- vapply(Gseq, function(k)    dnorm(data[[k]], mus[k], sqrt(sigs), log=TRUE), numeric(N))
    }, V= {
      sigx        <- sqrt(sigs);
      densi       <- vapply(Gseq, function(k)    dnorm(data[[k]], mus[k], sigx[k],    log=TRUE), numeric(N))
    } )
    check         <- is.infinite(densi)        & densi > 0
    if(any(check))          {
      densi[which(check, arr.ind=TRUE)[1],]   <- rep(0, G)
    }
    densi         <- densi  + if(is.matrix(log.pis) || missing(log.pis)) log.pis else matrix(log.pis, nrow=N, ncol=G, byrow=TRUE)
      if(logarithm)  densi    else exp(densi)
  }

#' Compute the Z matrix and log-likelihood for MoEClust models
#'
#' Computes the responsibility matrix Z and the log-likelihood for parameterised MVN mixtures of experts, with the aid of \code{\link{MoE_dens}}.
#' @inheritParams MoE_dens
#' @param Dens (Optional) A numeric matrix whose \code{[i,k]}-th entry is the \strong{log}-density of observation \emph{i} in component \emph{k}, scaled by the mixing proportions, typically obtained by \code{\link{MoE_dens}} but this is not necessary. If this is supplied, all other arguments are ignored, otherwise \code{\link{MoE_dens}} is called according to the other supplied arguments.
#'
#' @return A list containing two elements:
#' \itemize{
#' \item{\code{Z} - }{A matrix with N rows and G columns containing the probability of cluster membership for each of N observations and G clusters}
#' \item{\code{loglik} - }{The log-likelihood, computed efficiently via \code{\link[matrixStats]{rowLogSumExps}}}
#' }
#' @importFrom matrixStats "rowLogSumExps"
#' @export
#' @note This function is intended for joint use with \code{\link{MoE_dens}}, using the \strong{log}-densities.
#' @author Keefe Murphy - \href{keefe.murphy@ucd.ie}{<keefe.murphy@ucd.ie>}
#'
#' @seealso \code{\link{MoE_dens}}, \code{\link{MoE_clust}}, \code{\link[matrixStats]{rowLogSumExps}}, \code{\link[mclust]{mclustModelNames}}
#'
#' @examples
#' data(ais)
#' hema   <- ais[,3:7]
#' model  <- MoE_clust(hema, G=2, gating= ~ ais$BMI + ais$sex, model="EVE")
#' Dens   <- MoE_dens(modelName=model$best.mod$modelName, data=hema,
#'                    mus=model$mus, sigs=model$sigmas, log.pis=log(model$pis))
#'
#' # Construct the Z matrix and compute the log-likelihood
#' Estep  <- MoE_estep(Dens=Dens)
#' identical(Estep$Z, model$Z)
#' (ll    <- Estep$loglik)
#'
#' # Call MoE_estep directly
#' Estep2 <- MoE_estep(modelName=model$best.mod$modelName, data=hema,
#'                     mus=model$mus, sigs=model$sigmas, log.pis=log(model$pis))
#' identical(Estep2$loglik, ll)
  MoE_estep       <- function(modelName, data, mus, sigs, log.pis = 0, Dens = NULL) {
    if(missing(Dens)) {
      Dens        <- do.call(MoE_dens, as.list(match.call())[-1])
    } else if(!is.matrix(Dens) ||
              !is.numeric(Dens))                 stop("'Dens' must be a numeric matrix")
    norm          <- rowLogSumExps(Dens)
    Z             <- exp(sweep(Dens, 1, norm, "-"))
   #ll            <- sum(Z * Dens)             # complete log-likelihood
      return(list(Z = Z, loglik = sum(norm)))
  }

#' MoEClust BIC and ICL Model-Selection Criteria
#'
#' Computes the BIC (Bayesian Information Criterion) and ICL (Integrated Complete Likelihood) for parameterized mixture of experts models given the log-likelihood, the dimension of the data, the number of mixture components in the model, the numbers of parameters in the gating and expert networks respectively, and, for the ICL, the numbers of observations in each component.
#' @param modelName A character string indicating the model. The help file for \code{\link[mclust]{mclustModelNames}} describes the available models.
#' @param loglik The log-likelihood for a data set with respect to the Gaussian mixture model specified in the \code{modelName} argument.
#' @param N The number of observations in the data used to compute \code{loglik}.
#' @param D The dimension of the data used to compute \code{loglik}.
#' @param G The number of components in the Gaussian mixture model used to compute \code{loglik}.
#' @param gating.pen The number of parameters of the \emph{gating} network of the MoEClust model. Defaults to \code{G - 1}, which corresponds to no gating covariates. If covariates are included, this should be the number of regression coefficients in the fitted object. If there are no covariates and mixing proportions are further assumed to be present in equal proportion, \code{gating.pen} should be \code{0}.
#' @param expert.pen The number of parameters of the \emph{expert} network of the MoEClust model. Defaults to \code{G * D}, which corresponds to no expert covariates. If covariates are included, this should be the number of regression coefficients in the fitted object.
#' @param nn A vector of length \code{G} which sums to \code{N} giving the number of observations in each component. If supplied the ICL is also computed and returned, otherwise only the BIC.
#' @param delta Dirichlet hyperparameter for the prior on the mixing proportions. Defaults to 0.5. Only relevant for the ICL computation.
#'
#' @details The function is vectorized with respect to the arguments \code{modelName} and \code{loglik}.\cr
#'
#' If \code{model} is an object of class \code{"MoEClust"}, the number of parameters for the \code{gating.pen} and \code{expert.pen} are \code{length(model$gating)} and \code{sum(lengths(model$expert))}, respectively.
#' @importFrom mclust "mclustModelNames" "nVarParams"
#' @return A simplified array containing the BIC (and if \code{nn} is supplied, also the ICL) for each of the given input arguments.
#' @note In order to speed up repeated calls to the function inside \code{\link{MoE_clust}}, no checks take place.
#' @export
#' @author Keefe Murphy - \href{keefe.murphy@ucd.ie}{<keefe.murphy@ucd.ie>}
#'
#' @references Biernacki, C., Celeux, G., Govaert, G. (2000). Assessing a mixture model for clustering with the integrated completed likelihood. \emph{IEEE Trans. Pattern Analysis and Machine Intelligence}, 22(7): 719-725.
#'
#' @seealso \code{\link{MoE_clust}}, \code{\link[mclust]{nVarParams}}, \code{\link[mclust]{mclustModelNames}}
#'
#' @examples
#' MoE_crit(modelName=c("VVI", "VVE", "VVV"),
#'          loglik=c(-4036.99, -3987.12, -3992.45),
#'          N=120, D=8, G=3, nn=c(54, 38, 28))
#'
#' data(CO2data)
#' model <- MoE_clust(CO2data$CO2, G=2, expert= ~ CO2data$GNP)
#' (bic2 <- MoE_crit(modelName=model$best.model$modelName, N=length(CO2data$CO2),
#'                   loglik=max(model$log.likes), expert.pen=sum(lengths(model$expert)),
#'                   D=1, G=model$best.model$G, nn=tabulate(model$classes))["bic",])
#' identical(bic2, unname(model$best.model$bic))
  MoE_crit        <- Vectorize(function(modelName, loglik, N, D, G, gating.pen = G - 1, expert.pen = G * D, nn = NULL, delta = 0.5) {
    bicx          <- 2 * loglik - (nVarParams(modelName, D, G) + expert.pen + gating.pen) * log(N)
      return(c(bic = bicx, icl = if(!missing(nn)) bicx + ifelse(G == 1, 0, lgamma(G * delta) + sum(lgamma(nn + delta)) - G * lgamma(delta) - lgamma(N + G * delta))))
  }, vectorize.args = c("modelName", "loglik"), SIMPLIFY="array")

#' Set control values for use with MoEClust
#'
#' Supplies a list of arguments (with defaults) for use with \code{\link{MoE_clust}}.
#' @param criterion When either \code{G} or \code{modelNames} is a vector, \code{criterion} determines whether the "\code{bic}" (Bayesian Information Criterion) or "\code{icl}" (Integrated Complete Likelihood) is used to determine the 'best' model when gathering output. Note that both criteria will be returned in any case.
#' @param stopping The criterion used to assess convergence of the EM algorithm. The default (\code{"aitken"}) uses Aitken's acceleration method via \code{\link{MoE_aitken}}, otherwise the \code{"relative"} change in log-likelihood is monitored (which may be less strict). Both stopping rules are ultimately governed by \code{tol[1]}.
#' @param init.z The method used to initialise the cluster labels. Defaults to a hierarchical clustering tree as per \code{\link[mclust]{hc}}. Other options include \code{kmeans}, \code{random} initialisation, and a full run of \code{\link[mclust]{Mclust}}, although this last option is only permitted if there are \code{gating} &/or \code{expert} covariates within \code{\link{MoE_clust}}.
#' @param eps A scalar tolerance associated with deciding when to terminate computations due to computational singularity in covariances. Smaller values of eps allow computations to proceed nearer to singularity. The default is the relative machine precision \code{.Machine$double.eps}, which is approximately \emph{2e-16} on IEEE-compliant machines.
#' @param tol A vector of length two giving relative convergence tolerances for the log-likelihood and for parameter convergence in the inner loop for models with iterative M-step ("VEI", "EVE", "VEE", "VVE", "VEV"), respectively. The default is \code{c(1e-05,sqrt(.Machine$double.eps))}. If only one number is supplied, it is used as the tolerance for the outer iterations and the tolerance for the inner iterations is as in the default.
#' @param itmax A vector of length two giving integer limits on the number of EM iterations and on the number of iterations in the inner loop for models with iterative M-step ("VEI", "EVE", "VEE", "VVE", "VEV"), respectively. The default is \code{c(.Machine$integer.max, .Machine$integer.max)} allowing termination to be completely governed by \code{tol}. If only one number is supplied, it is used as the iteration limit for the outer iteration only.
#' @param equalPro Logical variable indicating whether or not the mixing proportions are equal in the model. Default: \code{equalPro = FALSE}. Only relevant when \code{gating} covariates are \emph{not} supplied within \code{\link{MoE_clust}}.
#' @param warn.it A single number giving the iteration count at which a warning will be printed if the EM algorithm has failed to converge. Defaults to \code{0}, i.e. no warning (which is true for any \code{warn.it} value less than \code{3}), otherwise the message is printed regardless of the value of \code{verbose}. If non-zero, \code{warn.it} should be moderately large, but obviously less than \code{itmax[1]}
#' @param verbose Logical indicating whether to print messages pertaining to progress to the screen during fitting. By default is \code{TRUE} if the session is interactive, and \code{FALSE} otherwise. If \code{FALSE}, warnings and error messages will still be printed to the screen, but everything else will be suppressed.
#'
#' @details \code{MoE_control} is provided for assigning values and defaults within \code{\link{MoE_clust}}.
#' @return A named list in which the names are the names of the arguments and the values are the values supplied to the arguments.
#' @export
#' @author Keefe Murphy - \href{keefe.murphy@ucd.ie}{<keefe.murphy@ucd.ie>}
#'
#' @seealso \code{\link{MoE_clust}}, \code{\link{MoE_aitken}}, \code{\link[mclust]{hc}}
#'
#' @examples
#' cont <- MoE_control(criterion="icl", itmax=10000, warn.it=12)
#'
#' data(CO2data)
#' res  <- MoE_clust(CO2data$CO2, G=2, expert = ~ CO2data$GNP, control=cont)
  MoE_control     <- function(criterion = c("bic", "icl"), stopping = c("aitken", "relative"),
                              init.z = c("hc", "kmeans", "mclust", "random"), eps = .Machine$double.eps,
                              tol = c(1e-05, sqrt(.Machine$double.eps)), itmax = c(.Machine$integer.max,
                              .Machine$integer.max), equalPro = FALSE, warn.it = 0, verbose = interactive()) {
    if(!is.character(criterion))                 stop("'criterion' must be a character vector of length 1")
    criterion     <- match.arg(criterion)
    if(!is.character(stopping))                  stop("'stopping' must be a character vector of length 1")
    stopping      <- match.arg(stopping)
    if(!is.character(init.z))                    stop("'init.z' must be a character vector of length 1")
    init.z        <- match.arg(init.z)
    if(length(eps) > 2)                          stop("'eps' can be of length at most 2")
    if(any(eps     < 0))                         stop("'eps' is negative")
    if(any(eps    >= 1))                         stop("'eps' is not less than 1")
    if((len.tol   <- length(tol))   > 2)         stop("'tol' can be of length at most 2")
    if(any(tol     < 0))                         stop("'tol' is negative")
    if(any(tol    >= 1))                         stop("'tol' is not less than 1")
    if(any(itmax   < 0))                         stop("'itmax' is negative")
    if(len.tol    == 1)        tol <- rep(tol, 2)
    if((len.itmax <- length(itmax)) > 2)         stop("'itmax' can be of length at most 2")
    if(len.itmax  == 1)      itmax <- c(itmax, .Machine$integer.max)
    inf           <- is.infinite(itmax)
    if(any(inf))        itmax[inf] <- .Machine$integer.max
    if(length(equalPro) > 1 ||
       !is.logical(equalPro))                    stop("'equalPro' must be a single logical indicator")
    if(length(verbose)  < 1 ||
       !is.logical(verbose))                     stop("'verbose' must be a single logical indicator")
    if(length(warn.it)  > 1 ||
       !is.numeric(warn.it))                     stop("'warn.it' must be a numeric vector of length 1")
      list(criterion = criterion, stopping = stopping, init.z = init.z, eps = eps, tol = tol,
           itmax = itmax, equalPro = equalPro, warn.it=warn.it, verbose = verbose)
  }

#' Aitken Acceleration
#'
#' Calculates the Aitken acceleration estimate of the final converged maximized log-likelihood.
#' @param log.likes A vector of three consecutive log-likelihood values. These three values should be in ascending order, though this is not checked.
#'
#' @details The final converged maximized log-likelihood can be used to determine convergence of the EM algorithm within \code{\link{MoE_clust}}, i.e. by checking whether the absolute difference in the current log-likelihood and estimated final converged maximised log-likelihood is less than some tolerance.
#' @note Within \code{\link{MoE_clust}}, as specified by the \code{stopping} argument of \code{\link{MoE_control}}, this is the default method used to assess convergence. The other option monitors the relative change in log-likelihood against some tolerance.
#' @return A list with the following components:
#' \itemize{
#' \item{\code{ll} - }{The most current estimate for the log-likelihood.}
#' \item{\code{linf} - }{An estimate of the final converged maxmised log-likelihood.}
#' \item{\code{a} - }{The Aitken acceleration value where \code{0 <= a <= 1}.}
#' \item{\code{diff} - }{The absolute difference between \code{ll} and \code{inf}. Note that if \code{a} goes to exactly \code{0} this will be also set to \code{0}.}
#' }
#' @export
#' @references Boehning, D., Dietz, E., Schaub, R., Schlattmann, and Lindsay, B. (1994, June). The distribution of the likelihood ratio for mixtures of densities from the one-parameter exponential family. \emph{Annals of the Institute of Statistical Mathematics}, 46(2): 373-388.
#'
#' @seealso \code{\link{MoE_control}}
#' @examples
#' MoE_aitken(-c(451.67, 442.84, 436.59))
#'
#' MoE_aitken(-c(359.11, 359.109, 359.109))$diff < 1e-05
  MoE_aitken      <- function(log.likes) {
    l1            <- log.likes[1]
    l2            <- log.likes[2]
    l3            <- log.likes[3]
    if(any(is.infinite(log.likes)))  {
      linf        <- diff     <- Inf
      a           <- NA
    } else {
      a           <- ifelse(l2 > l1, (l3 - l2) / (l2 - l1),       0)
      linf        <- ifelse(a  < 1,   l3 + (l3 - l2) / (1 - a), Inf)
      diff        <- abs(linf  - l3)
    }
      list(ll = l3, linf = linf, a = a, diff = ifelse(identical(a, 0), 0, diff))
  }
#
