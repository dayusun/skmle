# Changelog

## skmle 0.1.0

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
