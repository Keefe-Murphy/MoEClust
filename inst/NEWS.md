__MoEClust: Finite Gaussian Mixtures of Experts -__   
=======================================================
__Parsimonious Model-Based Clustering with Covariates__  
=======================================================

### New Features & Improvements
* New plotting function `MoE_Uncertainty` added (callable within `plot.MoEClust`):  
  visualises clustering uncertainty in the form of a barplot or an ordered profile plot,  
  with or without reference to the true labels in both cases.
* Added `drop.break` arg. to `MoE_control` for further control over the extra initialisation  
  step invoked in the presence of expert covariates (see Documentation for details).

### Bug Fixes & Miscellaneous Edits
* Fixed point-size, transparency & plotting symbols when `response.type="uncertainty"`  
  within `MoE_gpairs` to better conform to `mclust`: previously no transparency.
* `subset` arg. to `MoE_gpairs` now allows `data.ind=0` or `cov.ind=0`, allowing plotting of  
  response variables or plotting of the covariates to be suppressed entirely.
* `sigs` arg. to `MoE_dens` and `MoE_estep` must now be a variance object,  
   as per the `variance` object in the  parameters list from `MoE_clust` & `mclust` output.
* `resid.data` now returned by `MoE_clust` as a list, to better conform to `MoE_dens`.
* Removed redundant extra M-step after convergence for models without expert covariates.
* Removed redundant & erroneous `resid` & `residuals` args. to `as.Mclust` & `MoE_gpairs`.
* Added extra checks for invalid gating &/or expert covariates within `MoE_clust`.

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
