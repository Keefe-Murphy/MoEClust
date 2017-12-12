__MoEClust: Finite Gaussian Mixtures of Experts -__   
=======================================================
__Parsimonious Model-Based Clustering with Covariates__  
=======================================================

## MoEClust v1.1.0 - (_2<sup>nd</sup> release [minor update]: 2017-12-12_)
* `MoE_control` gains the `noise.gate` argument (defaults to `TRUE`): when `FALSE`,  
  the noise component's mixing proportion isn't influenced by gating network covariates.
* `MoE_control` gains the `hcUse` argument (defaults to `"VARS"` as per old `mclust` versions).
* `MoE_mahala` gains the `squared` argument + speed/matrix-inversion improvements.
* `G=0` models with noise component only can now be fit without having to supply `noise.init`.
* Speed-ups, incl. functions from `matrixStats` (on which `MoEClust` already depended).
* Fixed bugs in `MoE_mahala` and expert network estimation with a noise component.
* Removed dependency on the `grid` library.
* Many documentation improvements.

## MoEClust v1.0.0 - (_1<sup>st</sup> release: 2017-11-28_)
