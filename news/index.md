# Changelog

## skmle 0.1.0

- Initial CRAN release.

### Survival outcomes with sparse longitudinal covariates

- [`skmle()`](https://www.sundayu.me/skmle/reference/skmle.md) fits the
  transformed hazards family (proportional hazards, additive hazards and
  the Box-Cox transformations between them) by the sieve maximum
  kernel-weighted log-likelihood estimator of Sun, Sun, Zhao and Cao
  (2025).
- [`kee_cox()`](https://www.sundayu.me/skmle/reference/kee_cox.md) and
  [`kee_additive()`](https://www.sundayu.me/skmle/reference/kee_additive.md)
  fit the proportional and additive hazards models by kernel-weighted
  estimating equations.
- [`skmle_cv()`](https://www.sundayu.me/skmle/reference/skmle_cv.md)
  selects the bandwidth by subject-level cross-validation and refits on
  the full data.
- [`sim_skmle_data()`](https://www.sundayu.me/skmle/reference/sim_skmle_data.md)
  simulates sparse longitudinal survival data.
- Half and full Epanechnikov kernels
  ([`w_epan_half()`](https://www.sundayu.me/skmle/reference/kernel-weights.md),
  [`w_epan_full()`](https://www.sundayu.me/skmle/reference/kernel-weights.md)),
  and a `weight` argument that takes any function of the standardised
  lag.

### Asynchronous longitudinal outcomes

- [`kee_async()`](https://www.sundayu.me/skmle/reference/kee_async.md)
  fits generalised linear models with time-invariant coefficients by the
  kernel-weighted estimating equations of Cao, Zeng and Fine (2015),
  with identity, log and logistic links and a half or full kernel.
- [`kee_async_td()`](https://www.sundayu.me/skmle/reference/kee_async_td.md)
  estimates time-dependent coefficients pointwise.
- [`sim_async_data()`](https://www.sundayu.me/skmle/reference/sim_async_data.md)
  simulates a response and a covariate on independent observation-time
  streams.

### Shared

- [`coef()`](https://rdrr.io/r/stats/coef.html),
  [`vcov()`](https://rdrr.io/r/stats/vcov.html),
  [`confint()`](https://rdrr.io/r/stats/confint.html),
  [`nobs()`](https://rdrr.io/r/stats/nobs.html),
  [`print()`](https://rdrr.io/r/base/print.html),
  [`summary()`](https://rdrr.io/r/base/summary.html) and
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) methods, and
  the broom generics
  [`tidy()`](https://generics.r-lib.org/reference/tidy.html),
  [`glance()`](https://generics.r-lib.org/reference/glance.html) and
  [`augment()`](https://generics.r-lib.org/reference/augment.html).
