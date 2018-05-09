__MoEClust: Gaussian Parsimonious Clustering Models -__   
=======================================================
__with Gating and Expert Network Covariates__  
=======================================================

## MoEClust v1.2.0 - (_3<sup>rd</sup> release [minor update]: 2018-02-18_)
### New Features & Improvements
* New `predict.MoEClust` function added: predicts cluster membership probability,  
  MAP classification, and fitted response, using only new covariates or new covariates &  
  new response data with noise components (and the `noise.gate` option) accounted for.
* New plotting function `MoE_Uncertainty` added (callable within `plot.MoEClust`):  
  visualises clustering uncertainty in the form of a barplot or an ordered profile plot,  
  allowing reference to the made to the true labels, or not, in both cases.
* Specifying `response.type="density"` to `MoE_gpairs` now works properly for models with  
  gating &/or expert network covariates. Previous approach which evaluated the density using  
  averaged gates &/or average means replaced by more computationally expensive but correct  
  approach, which evaluates MVN density for every observation individually and then averages.
* Added `drop.break` arg. to `MoE_control` for further control over the extra initialisation  
  step invoked in the presence of expert covariates (see Documentation for details).
* Sped-up `MoE_dens` for the `EEE` & `VVV` models by using already available Cholesky factors.
* New `MoE_control` arg. `km.args` specifies `kstarts` & `kiters` when `init.z="kmeans"`.
    * Consolidated args. related to `init.z="hc"` & noise into `hc.args` & `noise.args`.

### Bug Fixes & Miscellaneous Edits
* Fixed point-size, transparency, & plotting symbols when `response.type="uncertainty"`  
  within `MoE_gpairs` to better conform to `mclust`: previously no transparency.
* `subset` arg. to `MoE_gpairs` now allows `data.ind=0` or `cov.ind=0`, allowing plotting of  
  response variables or plotting of the covariates to be suppressed entirely.
* Clarified MVN ellipses in `MoE_gpairs` plots.
* `sigs` arg. to `MoE_dens` and `MoE_estep` must now be a variance object, as per `variance`  
  in the  parameters list from `MoE_clust` & `mclust` output, the number of  clusters `G`,  
  variables `d` & `modelName` is inferred from this object: the arg. `modelName` was removed.
* `resid.data` now returned by `MoE_clust` as a list, to better conform to `MoE_dens`.
* Renamed functions `MoE_aitken` & `MoE_qclass` to `aitken` & `quant_clust`, respectively.
* Fixed small bug calculating `linf` within `aitken` & the associated stopping criterion.
* Final `linf` estimate now returned for optimal model when `stopping="aitken"` & G > 1.
* Removed redundant extra M-step after convergence for models without expert covariates.
* Removed redundant & erroneous `resid` & `residuals` args. to `as.Mclust` & `MoE_gpairs`.
* `MoE_plotCrit`, `MoE_plotGate` & `MoE_plotLogLik` now invisibly return revelant quantities.
* Corrected degrees of freedom calculation for `G=0` models when `noise.init` is not supplied.
* Improved `drop_levels` to handle alphanumeric variable names and ordinal variables.
* Interactions and higher-order terms are now accounted for within `drop_constants`.
* Replaced certain instances of `is.list(x)` with `inherits(x, "list")` for stricter checking.
* Added extra checks for invalid gating &/or expert covariates within `MoE_clust`.
* Added `mclust::clustCombi/clustCombiOptim` examples to `as.Mclust` documentation.
* Added extra precautions for empty clusters at initialisation.
* Added utility function `MoE_news` for accessing this `NEWS` file.
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
