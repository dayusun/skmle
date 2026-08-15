# Changelog

## skmle 0.1.0

### Scope

The package covers two settings, not one. Alongside the transformed
hazards models for survival outcomes it now fits generalised linear
models for asynchronous longitudinal outcomes, where the response and
the covariate are recorded on different time grids. The `Title` drops
its “for Survival Models” restriction, and the `Description`, README,
package help page and tutorial have been rewritten to present the two
settings on equal footing rather than treating the second as an extra.

### Asynchronous longitudinal data

- [`kee_async()`](https://dayusun.github.io/skmle/reference/kee_async.md)
  and
  [`kee_async_td()`](https://dayusun.github.io/skmle/reference/kee_async_td.md)
  implement the kernel-weighted estimating equations of Cao, Zeng and
  Fine (2015) for a longitudinal response and a longitudinal covariate
  observed on **different** time grids, with time-invariant and
  time-dependent coefficients respectively. Identity, log and logistic
  links are supported.
- [`sim_async_data()`](https://dayusun.github.io/skmle/reference/sim_async_data.md)
  simulates from their Section 4 design.
- [`plot()`](https://rdrr.io/r/graphics/plot.default.html) and
  [`print()`](https://rdrr.io/r/base/print.html) methods for the
  `kee_td` coefficient curves.
- Both estimators are backed by C++: pairs are enumerated once and
  collapsed onto covariate rows, so the pair design never enters the
  Newton loop. The time-dependent weight factorises over the two
  occasion indices, so no pair is enumerated there at all.

### Half and full kernels

- [`skmle()`](https://dayusun.github.io/skmle/reference/skmle.md),
  [`kee_cox()`](https://dayusun.github.io/skmle/reference/kee_cox.md),
  [`kee_additive()`](https://dayusun.github.io/skmle/reference/kee_additive.md),
  [`skmle_cv()`](https://dayusun.github.io/skmle/reference/skmle_cv.md),
  [`kee_async()`](https://dayusun.github.io/skmle/reference/kee_async.md)
  and
  [`kee_async_td()`](https://dayusun.github.io/skmle/reference/kee_async_td.md)
  take a `one_sided` argument. `TRUE` is the half kernel: only covariate
  observations preceding a time inform it. `FALSE` is the full,
  two-sided kernel. The survival estimators default to `TRUE`, the
  estimator as published; the asynchronous ones default to `FALSE`, the
  kernel of Cao, Zeng and Fine.
- For the survival estimators the switch reaches the risk-set averages
  in the C++ backend as well as the row weights, so the two halves of an
  estimator cannot disagree.
  [`skmle_cv()`](https://dayusun.github.io/skmle/reference/skmle_cv.md)
  also uses it inside the fold loop, which previously had the half
  kernel hardcoded, so the bandwidth is now selected under the kernel
  the refit uses.
- The Epanechnikov kernel had been written out inline in three fitting
  functions; it now lives in one internal helper, so the half/full
  switch cannot be applied to two of the three.

### Initial release

- Initial CRAN release.
- [`skmle()`](https://dayusun.github.io/skmle/reference/skmle.md) fits
  transformed hazards models by sieve maximum kernel-weighted
  log-likelihood estimation (SMKLE).
- [`kee_cox()`](https://dayusun.github.io/skmle/reference/kee_cox.md)
  and
  [`kee_additive()`](https://dayusun.github.io/skmle/reference/kee_additive.md)
  fit Cox and additive hazards models via kernel-weighted estimating
  equations.
- [`skmle_cv()`](https://dayusun.github.io/skmle/reference/skmle_cv.md)
  selects the kernel bandwidth by K-fold cross-validation.
- [`sim_skmle_data()`](https://dayusun.github.io/skmle/reference/sim_skmle_data.md)
  simulates survival data with sparse, intermittently observed
  longitudinal covariates.
- [`plot()`](https://rdrr.io/r/graphics/plot.default.html),
  [`print()`](https://rdrr.io/r/base/print.html) and
  [`summary()`](https://rdrr.io/r/base/summary.html) methods for fitted
  objects.
