# skmle 0.1.0

* Initial CRAN release.

## Survival outcomes with sparse longitudinal covariates

* `skmle()` fits the transformed hazards family (proportional hazards,
  additive hazards and the Box-Cox transformations between them) by the sieve
  maximum kernel-weighted log-likelihood estimator of Sun, Sun, Zhao and Cao
  (2025).
* `kee_cox()` and `kee_additive()` fit the proportional and additive hazards
  models by kernel-weighted estimating equations.
* `skmle_cv()` selects the bandwidth by subject-level cross-validation and
  refits on the full data.
* `sim_skmle_data()` simulates sparse longitudinal survival data.
* Half and full Epanechnikov kernels (`w_epan_half()`, `w_epan_full()`), and a
  `weight` argument that takes any function of the standardised lag.

## Asynchronous longitudinal outcomes

* `kee_async()` fits generalised linear models with time-invariant
  coefficients by the kernel-weighted estimating equations of Cao, Zeng and
  Fine (2015), with identity, log and logistic links and a half or full kernel.
* `kee_async_td()` estimates time-dependent coefficients pointwise.
* `sim_async_data()` simulates a response and a covariate on independent
  observation-time streams.

## Shared

* `coef()`, `vcov()`, `confint()`, `nobs()`, `print()`, `summary()` and
  `plot()` methods, and the broom generics `tidy()`, `glance()` and
  `augment()`.
