__MoEClust: Gaussian Parsimonious Clustering Models__   
=======================================================
__with Gating and Expert Network Covariates__ 
=======================================================
__and a Noise Component__
=======================================================

## MoEClust v1.3.1 - (_9<sup>th</sup> release [patch update]: 2020-05-12_)
### New Features, Improvements, Bug Fixes & Miscellaneous Edits
* Maintenance release for compatibility with R 4.0.0 - minor edits.
* `summary.MoEClust` gains the printing-related arguments `classification=TRUE`,  
  `parameters=FALSE`, and `networks=FALSE` (thanks to a request from Prof. Kamel Gana).
* Related improvements to `print`/`summary` methods for `MoE_gating` & `MoE_expert` objects.
* Minor speed-up for G=1 models with expert network covariates.
* Improvements to `MoE_plotGate`, with new `type`, `pch`, and `xlab` defaults.
* Added informative `dimnames` to returned `parameters` from `MoE_clust()`.
* Documentation, vignette, examples, and references improvements.

## MoEClust v1.3.0 - (_8<sup>th</sup> release [minor update]: 2020-03-30_)
### New Features, Improvements, Bug Fixes & Miscellaneous Edits
* Various fixes and improvements to initialisation when there are expert network covariates:  
    * `MoE_mahala` now correctly uses the covariance of `resids` rather than the response.
    * New `MoE_mahala` arg. `identity` allow use of Euclidean distance instead:  
    this argument can also be passed via `exp.init$identity` to `MoE_control`.
    * Convergence of the initialisation procedure now explictly monitored & sped-up.
    * Values of the criterion being minimised are now returned as an attribute.
    * The number of iterations of the initialisation algorithm are also returned as an attribute.
    * `MoE_control` arg. `exp.init$max.init` now defaults to `.Machine$integer.max`.
    * Improved checks on the `resids` arg. to `MoE_mahala`.
    * Greatly expanded the `MoE_mahala` examples.
* Improvements to `predict.MoEClust`:  
    * Now returns the predicted values of the gating and expert networks.
    * Now returns the predictions from the expert network of the most probable component  
    (`MAPy`), in addition to the (aggregated) predicted responses (`y`).
    * New arg. `MAPresids` governs whether residuals are computed against `MAPy` or `y`.
    * New arg. `use.y` (see documentation for details).
    * Now properly allows empty `newdata` for models with no covariates of any kind.
    * Fixed prediction for equal mixing proportion models when `discard.noise=FALSE`.
* Fixed small `MoE_stepwise` bugs when  
    * only one of `gating` or `expert` are supplied.
    * univariate response `data` are supplied.
    * moving from G=1 to G=2 with equal mixing proportions and no covariates.
    * discarding covariates present in the response data.
* Odds ratios now returned (and printed) when calling `summary` on `x$gating`.
* `noise_vol` now returns correction location for univariate data when `reciprocal=TRUE`.
* Spell-checking of documentation and fixes to `donttest` examples.

## MoEClust v1.2.4 - (_7<sup>th</sup> release [patch update]: 2019-12-11_)
### New Features, Improvements, Bug Fixes & Miscellaneous Edits
* Fixed small bugs in `MoE_stepwise`:  
    * Improved checks on `network.data` and `data`.
    * Prevented `z.list` from being suppliable.  
    * Fixes when`equalPro="yes"` & `noise=TRUE`.
    * Fixes for supplying optional `MoE_control` arguments (also for `MoE_clust`).
    * Prevented termination if adding a component fails,  
    provided at least one other step doesn't fail.
* Fixed `discard.noise=TRUE` behaviour for `MoE_clust`, `predict.MoEClust`, &  
  `residuals.MoEClust` for models with a noise component fitted via `"CEM"`.
* Minor fixes to `noise_vol` function and handling of `noise.meth` arg. to `MoE_control`.
* Slight speed-up to E-step/C-step for models with a noise component.
* Initial allocation matrices now stored as attributes to `MoE_clust` output (see `?MoE_control`).
* Anti-aliasing of vignette images.
* Updated citation info after publication in _Advances in Data Analysis and Classification_.

## MoEClust v1.2.3 - (_6<sup>th</sup> release [patch update]: 2019-07-29_)
### New Features, Improvements, Bug Fixes & Miscellaneous Edits
* Exported function `MoE_stepwise` for conducting a greedy forward stepwise  
  search to find the optimal model in terms of the number of components, GPCM  
  covariance parameterisation, and the subsets of gating/expert network covariates.
* `MoE_control` & `predict.MoEClust` gain the arg. `discard.noise`:  
  Default of `FALSE` retains old behaviour (see documentation for details).
* `MoE_control` gains the arg. `z.list` and the `init.z` arg. gets the option `"list"`:  
  this allows manually supplying (soft or hard) initial cluster allocation matrices.
* New args. and small fixes added to `MoE_gpairs`:  
    * `uncert.cov` arg. added to control uncertainty point-size in panels with covariates.
    * `density.pars` gains arg. `label.style`.
    * `scatter.pars` & `stripplot.pars` gain args. `noise.size` & `size.noise`.
    * `barcode.pars$bar.col` slightly fixed from previous update.
    * Colours for `"violin"` type plots now accurate for MAP panels.
* Slight speed-up to `noise_vol` when `method="ellipsoidhull"`.
* Small fix to `predict.MoEClust` when `resid=TRUE` for models with expert covariates.
* Small fix related to `...` construct for `residuals.MoEClust`.
* All printing related to noise-only models no longer shows the model name (there is none!).
* Other small fixes to `print.MoEClust`, `print.summary_MoEClust`, & `print.MoECompare`.
* Cosmetic fix to returned `gating` objects for `equalPro=TRUE` models. 
* Removed `parallel` package from `Suggests:`.

## MoEClust v1.2.2 - (_5<sup>th</sup> release [patch update]: 2019-05-15_)
### New Features, Improvements, Bug Fixes & Miscellaneous Edits
* `noise_vol` now also returns the location of the centre of mass of the region  
  used to estimate the hypervolume, regardless of the method employed. This fixes:    
    * `predict.MoEClust` for any models with a noise component (see below).
    * The summary of means for models with expert covariates and a noise component.
    * The location of the MVN ellipses for such models in `MoE_gpairs` (see below).
* Furthermore, calculation of the hypervolume in `noise_vol` for data with >2 dimensions  
  is now correct when `method="ellipsoidhull"`, owing to a bug in the `cluster` package.
* Other fixes and speed-ups for the `MoE_gpairs` plotting function:  
    * Added arg. `expert.covar` (& also to `as.Mclust` function).
    * Fixed location of MVN ellipses for models with noise & expert covariates (see above).
    * Fixes when `response.type="density"` for all models with a noise component.
    * Speed-up when `response.type="density"` for models with covariates of any kind.
    * Fixes to labelling for models with a noise component.
    * Fixed handling of `subset$data.ind` & `subset$cov.ind` arguments.
    * Barcode type plots now have colour for panels involving the MAP classification.
    * Barcode type plots now respect the arg. `buffer`.
    * Use of colour in `MoE_plotGate` is now consistent with `MoE_gpairs`.
* Fixes to how `gating` & `expert` formulas are handled:  
    * Allowed specification of formulas with dropped variables of the form `~.-a-b`.
    * Allowed formulas with no intercept of the form `~c-1`.
    * Allowed interaction effects, transformations and higher-order terms using `I()`.
    * Small related fixes to `drop_levels` & `drop_constants` functions.  
* `MoE_compare` gains arg. `noise.vol` for overriding the `noise.meth` arg.:  
  this allows specifying an improper uniform density directly via the (hyper)volume,  
  & hence adding noise to models for high-dimensional data for which `noise_vol()` fails.
* Fixed bug for `equalPro` models with noise component, and also added `equalNoise` arg.  
  to `MoE_control`, further controlling `equalPro` in the presence of a noise component.
* Fixes to `predict.MoEClust` for the following special cases:  
    * Fixes for any models with a noise component (see `noise_vol` comment above).
    * Accounted for predictions of single observations for models with a noise component.
    * Accounted for models with equal mixing proportions.
* Accounted for categorical covariates in the `x.axis` arg. to `MoE_plotGate`.
* `tau0` can now also be supplied as a vector in the presence of gating covariates.
* Fix to `expert_covar` for univariate models. 
* Slight `MoE_estep` speed-up due to removal of unnecessary `sweep()`.
* Small fixes for when `clustMD` is invoked, and added `snow` package to `Suggests:`.
* The `nnet` arg. `MaxNWts` now passable to gating network `multinom` call via `MoE_control`.
* Improved printing of output and handling of ties, especially for `MoE_compare`.
* Many documentation and vignette improvements.

## MoEClust v1.2.1 - (_4<sup>th</sup> release [patch update]: 2018-12-11_)
### New Features, Improvements, Bug Fixes & Miscellaneous Edits
* New `MoE_control` arg. `algo` allows model fitting using the `"EM"` or `"CEM"` algorithm:  
    * Related new function `MoE_cstep` added.
    * Extra `algo` option `"cemEM"` allows running EM starting from convergence of CEM.
* Added `LOGLIK` to `MoE_clust` output, giving maximal log-likelihood values for all fitted models.
    * Behaves exactly as per `DF/ITERS`, etc., with associated printing/plotting functions.
    * Edited `MoE_compare`, `summary.MoEClust`, & `MoE_plotCrit` accordingly.
* New `MoE_control` arg. `nstarts` allows for multiple random starts when `init.z="random"`.
* New `MoE_control` arg. `tau0` provides another means of initialising the noise component.
* If `clustMD` is invoked for initialisation, models are now run more quickly in parallel.
* `MoE_plotGate` now allows a user-specified x-axis against which mixing proportions are plotted.
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
* `sigs` arg. to `MoE_dens` & `MoE_estep` must now be a variance object, as per `variance`  
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
