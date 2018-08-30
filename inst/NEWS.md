__MoEClust: Gaussian Parsimonious Clustering Models -__   
=======================================================
__with Gating and Expert Network Covariates__  
=======================================================

## MoEClust v1.2.0 - (_4<sup>th</sup> release [minor update]: 2018-08-25_)
### New Features & Improvements
* New `MoE_control` arg. `algo` allows model fitting using the `"EM"` or `"CEM"` algorithm:  
    * Related new function `MoE_cstep` added.
    * Extra `algo` option `"cemEM"` allows running EM starting from convergence of CEM.
* New `MoE_control` arg. `nstarts` allows for multiple random starts when `init.z="random"`.
* If `clustMD` is invoked for initialisation, models are now run more quickly in parallel.

### Bug Fixes & Miscellaneous Edits
* Fixed bug in checking for strictly increasing log-likelihood estimates.

## MoEClust v1.2.0 - (_3<sup>rd</sup> release [minor update]: 2018-08-24_)
### New Features & Improvements
* New `predict.MoEClust` function added: predicts cluster membership probability,  
  MAP classification, and fitted response, using only new covariates or new covariates &  
  new response data, with noise components (and the `noise.gate` option) accounted for.
* New plotting function `MoE_Uncertainty` added (callable within `plot.MoEClust`):  
  visualises clustering uncertainty in the form of a barplot or an ordered profile plot,  
  allowing reference to be made to the true labels, or not, in both cases.
* Specifying `response.type="density"` to `MoE_gpairs` now works properly for models with  
  gating &/or expert network covariates. Previous approach which evaluated the density using  
  averaged gates &/or averaged means replaced by more computationally expensive but correct  
  approach, which evaluates MVN density for every observation individually and then averages.
* Added `clustMD` package to `Suggests:`. New `MoE_control` argument `exp.init$clustMD`  
  governs whether categorical/ordinal covariates are also incorporated into the initialisation  
  when `isTRUE(exp.init$joint)` & `clustMD` is loaded (defaults to `FALSE`, works with noise). 
* Added `drop.break` arg. to `MoE_control` for further control over the extra initialisation  
  step invoked in the presence of expert covariates (see Documentation for details).
* Sped-up `MoE_dens` for the `EEE` & `VVV` models by using already available Cholesky factors.
* Other new `MoE_control` arguments:  
    * `km.args` specifies `kstarts` & `kiters` when `init.z="kmeans"`.
    * Consolidated args. related to `init.z="hc"` & noise into `hc.args` & `noise.args`.
    * `hc.args` now also passed to call to `mclust` when `init.z="mclust"`.
    * `init.crit` (`"bic"`/`"icl"`) controls selection of optimal `mclust`/`clustMD`  
       model type to initialise with (if `init.z="mclust"` or `isTRUE(exp.init$clustMD)`);  
       relatedly, initialisation now sped-up when `init.z="mclust"`.

### Bug Fixes & Miscellaneous Edits
* `ITERS` replaces `iters` as the matrix of the number of EM iterations in `MoE_clust` output:  
    * `iters` now gives this number for the optimal model.  
	  * `ITERS` now behaves like `BIC`/`ICL` etc. in inheriting the `"MoECriterion"` class.  
	  * `iters` now filters down to `summary.MoEClust` and the associated printing function.  
	  * `ITERS` now filters down to `MoE_compare` and the associated printing function.
* Fixed point-size, transparency, & plotting symbols when `response.type="uncertainty"`  
  within `MoE_gpairs` to better conform to `mclust`: previously no transparency.
* `subset` arg. to `MoE_gpairs` now allows `data.ind=0` or `cov.ind=0`, allowing plotting of  
  response variables or plotting of the covariates to be suppressed entirely.
* Clarified MVN ellipses in `MoE_gpairs` plots.
* `sigs` arg. to `MoE_dens` and `MoE_estep` must now be a variance object, as per `variance`  
  in the  parameters list from `MoE_clust` & `mclust` output, the number of  clusters `G`,  
  variables `d` & `modelName` is inferred from this object: the arg. `modelName` was removed.
* `MoE_clust` no longer returns an error if `init.z="mclust"` when no gating/expert network  
   covariates are supplied; instead, `init.z="hc"` is used to better reproduce `mclust` output.
* `resid.data` now returned by `MoE_clust` as a list, to better conform to `MoE_dens`.
* Renamed functions `MoE_aitken` & `MoE_qclass` to `aitken` & `quant_clust`, respectively.
* Rows of `data` w/ missing values now dropped for gating/expert covariates too (`MoE_clust`).
* Logical covariates in gating/expert networks now coerced to factors.
* Fixed small bug calculating `linf` within `aitken` & the associated stopping criterion.
* Final `linf` estimate now returned for optimal model when `stopping="aitken"` & G > 1.
* Removed redundant extra M-step after convergence for models without expert covariates.
* Removed redundant & erroneous `resid` & `residuals` args. to `as.Mclust` & `MoE_gpairs`.
* `MoE_plotCrit`, `MoE_plotGate` & `MoE_plotLogLik` now invisibly return revelant quantities.
* Corrected degrees of freedom calculation for `G=0` models when `noise.init` is not supplied.
* Fixed `drop_levels` to handle alphanumeric variable names and ordinal variables.
* Fixed `MoE_compare` when a mix of models with and without a noise component are supplied.
* Fixed `MoE_compare` when optimal model has to be re-fit due to mismatched `criterion`.
* Fixed y-axis labelling of `MoE_Uncertainty` plots.
* `print.MoECompare` now has a `digits` arg. to control rounding of printed output.
* Better handling of tied model-selection criteria values in `MoE_clust` & `MoE_compare`.
* Interactions and higher-order terms are now accounted for within `drop_constants`.
* Replaced certain instances of `is.list(x)` with `inherits(x, "list")` for stricter checking.
* Added extra checks for invalid gating &/or expert covariates within `MoE_clust`.
* Added `mclust::clustCombi/clustCombiOptim` examples to `as.Mclust` documentation.
* Added extra precautions for empty clusters: during initialisation & during EM.
* Added utility function `MoE_news` for accessing this `NEWS` file.
* Added message if optimum `G` is at either end of the range considered.
* Tidied indentation/line-breaks for `cat`/`message`/`warning` calls for printing clarity.
* Added line-breaks to `usage` sections of multi-argument functions.
* Corrected `MoEClust-package` help file (formerly just `MoEClust`).
* Many documentation clarifications.

## MoEClust v1.1.0 - (_2<sup>nd</sup> release [minor update]: 2018-02-06_)
### New Features & Improvements
* `MoE_control` gains the `noise.gate` argument (defaults to `TRUE`): when `FALSE`,  
  the noise component's mixing proportion isn't influenced by gating network covariates.
* `x$parameters$mean` is now reported as the posterior mean of the fitted values when  
  there are expert network covariates: when there are no expert covariates, the posterior  
  mean of the response is reported, as before. This effects the centres of the MVN ellipses  
  in response vs. response panels of `MoE_gpairs` plots when there are expert covariates.
* New function `expert_covar` used to account for variability in the means, in the presence  
  of expert covariates, in order to modify shape & size of MVN ellipses in visualisations.
* `MoE_control` gains the `hcUse` argument (defaults to `"VARS"` as per old `mclust` versions).
* `MoE_mahala` gains the `squared` argument + speedup/matrix-inversion improvements.
* Speed-ups, incl. functions from `matrixStats` (on which `MoEClust` already depended).
* The `MoE_gpairs` argument `addEllipses` gains the option `"both"`.

### Bug Fixes & Miscellaneous Edits
* Fixed bug when `equalPro=TRUE` in the presence of a noise component when there are  
  no gating covariates:  now only the mixing proportions of the non-noise components  
  are constrained to be equal, after accounting for the noise component.
* `MoE_gpairs` argument `scatter.type` gains the options `lm2` & `ci2` for further control  
  over gating covariates. Fixed related bug whereby `lm` & `ci` type plots were being  
  erroneously produced for panels involving pairs of continuous covariates only.
* Fixed bugs in `MoE_mahala` and in expert network estimation with a noise component.
* `G=0` models w/ noise component only can now be fitted without having to supply `noise.init`.
* `MoE_compare` now correctly prints noise information for sub-optimal models.
* Slight edit to criterion used when `stopping="relative"`: now conforms to `mclust`.
* Added `check.margin=FALSE` to calls to `sweep()`.
* Added `call.=FALSE` to all `stop()` messages.
* Removed dependency on the `grid` library.
* Many documentation clarifications.

## MoEClust v1.0.0 - (_1<sup>st</sup> release: 2017-11-28_)
