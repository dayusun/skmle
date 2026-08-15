# Package index

## Package overview

- [`skmle-package`](https://dayusun.github.io/skmle/reference/skmle-package.md)
  : skmle: Sieve Kernel Maximum Likelihood Estimation

## Survival models with sparse longitudinal covariates

- [`skmle()`](https://dayusun.github.io/skmle/reference/skmle.md) : Fit
  a Transformed Hazards Model by SMKLE
- [`skmle_cv()`](https://dayusun.github.io/skmle/reference/skmle_cv.md)
  [`print(`*`<cv.skmle>`*`)`](https://dayusun.github.io/skmle/reference/skmle_cv.md)
  : Select the Bandwidth by Cross-Validation
- [`kee_cox()`](https://dayusun.github.io/skmle/reference/kee_cox.md) :
  Fit a Cox-Type KEE Model
- [`kee_additive()`](https://dayusun.github.io/skmle/reference/kee_additive.md)
  : Fit an Additive Hazards KEE Model

## Asynchronous longitudinal regression

Kernel-weighted estimating equations of Cao, Zeng and Fine (2015) for a
response and a covariate observed on different time grids.

- [`kee_async()`](https://dayusun.github.io/skmle/reference/kee_async.md)
  : Asynchronous longitudinal regression, time-invariant coefficients
- [`kee_async_td()`](https://dayusun.github.io/skmle/reference/kee_async_td.md)
  [`print(`*`<kee_td>`*`)`](https://dayusun.github.io/skmle/reference/kee_async_td.md)
  : Asynchronous longitudinal regression, time-dependent coefficients

## Simulation

- [`sim_skmle_data()`](https://dayusun.github.io/skmle/reference/sim_skmle_data.md)
  : Simulate Sparse Longitudinal Survival Data
- [`sim_async_data()`](https://dayusun.github.io/skmle/reference/sim_async_data.md)
  : Simulate asynchronous longitudinal data

## Model summaries and plots

- [`plot(`*`<skmle>`*`)`](https://dayusun.github.io/skmle/reference/plot.skmle.md)
  : Plot the estimated baseline function for skmle model
- [`plot(`*`<kee_td>`*`)`](https://dayusun.github.io/skmle/reference/plot.kee_td.md)
  : Plot estimated coefficient curves
- [`summary(`*`<skmle>`*`)`](https://dayusun.github.io/skmle/reference/summary.skmle.md)
  : Summary for skmle object
- [`summary(`*`<kee>`*`)`](https://dayusun.github.io/skmle/reference/summary.kee.md)
  : Summary for kee object
- [`print(`*`<skmle>`*`)`](https://dayusun.github.io/skmle/reference/print.skmle.md)
  : Print skmle object
- [`print(`*`<kee>`*`)`](https://dayusun.github.io/skmle/reference/print.kee.md)
  : Print kee object
- [`print(`*`<summary.skmle>`*`)`](https://dayusun.github.io/skmle/reference/print.summary.skmle.md)
  : Print summary of skmle object
- [`print(`*`<summary.kee>`*`)`](https://dayusun.github.io/skmle/reference/print.summary.kee.md)
  : Print summary of kee object
