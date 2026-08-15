# Package index

## Package overview

- [`skmle-package`](https://www.sundayu.me/skmle/reference/skmle-package.md)
  : skmle: Sieve Kernel Maximum Likelihood Estimation

## Survival models with sparse longitudinal covariates

- [`skmle()`](https://www.sundayu.me/skmle/reference/skmle.md) : Fit a
  Transformed Hazards Model by SMKLE
- [`skmle_cv()`](https://www.sundayu.me/skmle/reference/skmle_cv.md)
  [`print(`*`<cv.skmle>`*`)`](https://www.sundayu.me/skmle/reference/skmle_cv.md)
  : Select the Bandwidth by Cross-Validation
- [`kee_cox()`](https://www.sundayu.me/skmle/reference/kee_cox.md) : Fit
  a Cox-Type KEE Model
- [`kee_additive()`](https://www.sundayu.me/skmle/reference/kee_additive.md)
  : Fit an Additive Hazards KEE Model

## Asynchronous longitudinal regression

Kernel-weighted estimating equations of Cao, Zeng and Fine (2015) for a
response and a covariate observed on different time grids.

- [`kee_async()`](https://www.sundayu.me/skmle/reference/kee_async.md) :
  Asynchronous longitudinal regression with time-invariant coefficients
- [`kee_async_cv()`](https://www.sundayu.me/skmle/reference/kee_async_cv.md)
  [`print(`*`<cv.kee_async>`*`)`](https://www.sundayu.me/skmle/reference/kee_async_cv.md)
  : Choose the bandwidth for an asynchronous longitudinal fit
- [`kee_async_td()`](https://www.sundayu.me/skmle/reference/kee_async_td.md)
  [`print(`*`<kee_td>`*`)`](https://www.sundayu.me/skmle/reference/kee_async_td.md)
  : Asynchronous longitudinal regression with time-dependent
  coefficients

## Simulation

- [`sim_skmle_data()`](https://www.sundayu.me/skmle/reference/sim_skmle_data.md)
  : Simulate Sparse Longitudinal Survival Data
- [`sim_async_data()`](https://www.sundayu.me/skmle/reference/sim_async_data.md)
  : Simulate asynchronous longitudinal data

## Model summaries and plots

- [`plot(`*`<skmle>`*`)`](https://www.sundayu.me/skmle/reference/plot.skmle.md)
  : Plot the estimated baseline function for skmle model
- [`plot(`*`<kee_td>`*`)`](https://www.sundayu.me/skmle/reference/plot.kee_td.md)
  : Plot estimated coefficient curves
- [`vcov(`*`<skmle>`*`)`](https://www.sundayu.me/skmle/reference/vcov.skmle.md)
  [`vcov(`*`<kee>`*`)`](https://www.sundayu.me/skmle/reference/vcov.skmle.md)
  [`vcov(`*`<kee_td>`*`)`](https://www.sundayu.me/skmle/reference/vcov.skmle.md)
  : Extract the covariance matrix of a fitted model
- [`nobs(`*`<skmle>`*`)`](https://www.sundayu.me/skmle/reference/nobs.skmle.md)
  [`nobs(`*`<kee>`*`)`](https://www.sundayu.me/skmle/reference/nobs.skmle.md)
  [`nobs(`*`<kee_td>`*`)`](https://www.sundayu.me/skmle/reference/nobs.skmle.md)
  : Number of subjects contributing to a fit
- [`confint(`*`<kee_td>`*`)`](https://www.sundayu.me/skmle/reference/confint.kee_td.md)
  : Wald confidence intervals for a coefficient curve
- [`plot(`*`<cv.kee_async>`*`)`](https://www.sundayu.me/skmle/reference/plot.cv.skmle.md)
  [`plot(`*`<cv.skmle>`*`)`](https://www.sundayu.me/skmle/reference/plot.cv.skmle.md)
  : Plot a cross-validation curve

## Tidy output

broom generics returning tibbles, so fits compose with the rest of a
tidy workflow.

- [`tidy(`*`<skmle>`*`)`](https://www.sundayu.me/skmle/reference/tidy.skmle.md)
  [`tidy(`*`<kee>`*`)`](https://www.sundayu.me/skmle/reference/tidy.skmle.md)
  [`tidy(`*`<kee_td>`*`)`](https://www.sundayu.me/skmle/reference/tidy.skmle.md)
  : Summarise a fit as a tibble
- [`glance(`*`<skmle>`*`)`](https://www.sundayu.me/skmle/reference/glance.skmle.md)
  [`glance(`*`<kee>`*`)`](https://www.sundayu.me/skmle/reference/glance.skmle.md)
  [`glance(`*`<kee_td>`*`)`](https://www.sundayu.me/skmle/reference/glance.skmle.md)
  : One-row summary of a fit
- [`augment(`*`<kee_async>`*`)`](https://www.sundayu.me/skmle/reference/augment.kee_async.md)
  : Add fitted means to the covariate table
- [`summary(`*`<skmle>`*`)`](https://www.sundayu.me/skmle/reference/summary.skmle.md)
  : Summary for skmle object
- [`summary(`*`<kee>`*`)`](https://www.sundayu.me/skmle/reference/summary.kee.md)
  : Summary for kee object
- [`print(`*`<skmle>`*`)`](https://www.sundayu.me/skmle/reference/print.skmle.md)
  : Print skmle object
- [`print(`*`<kee>`*`)`](https://www.sundayu.me/skmle/reference/print.kee.md)
  : Print kee object
- [`print(`*`<summary.skmle>`*`)`](https://www.sundayu.me/skmle/reference/print.summary.skmle.md)
  : Print summary of skmle object
- [`print(`*`<summary.kee>`*`)`](https://www.sundayu.me/skmle/reference/print.summary.kee.md)
  : Print summary of kee object
