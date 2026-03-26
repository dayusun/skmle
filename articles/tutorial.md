# Getting Started with skmle

## Overview

`skmle` fits transformed hazards survival models when longitudinal
covariates are observed sparsely and intermittently over time.

The package currently provides three main user-facing workflows:

1.  [`skmle()`](https://dayusun.github.io/skmle/reference/skmle.md) for
    the general transformed hazards model.
2.  [`kee_cox()`](https://dayusun.github.io/skmle/reference/kee_cox.md)
    for the proportional hazards estimating-equation approach.
3.  [`kee_additive()`](https://dayusun.github.io/skmle/reference/kee_additive.md)
    for the additive hazards estimating-equation approach.

This vignette shows the usual modeling workflow: simulate data, fit a
model, inspect the summary output, plot the estimated baseline
component, and select a bandwidth by cross-validation.

``` r
library(skmle)
library(survival)
```

## Simulate Example Data

``` r
set.seed(123)

dat <- sim_skmle_data(
  n = 200,
  mu = function(tt) 8 * (0.75 + (0.5 - tt)^2),
  mu_bar = 8,
  alpha = function(tt) 0.5 * 0.75 + 0.75 * (tt * (1 - sin(2 * pi * (tt - 0.25)))),
  beta = c(1, -0.5),
  s = 0,
  cen = 0.7
)

head(dat)
#> # A tibble: 6 × 6
#>   id        X delta covariates[,1]  [,2] obs_times censoring
#>   <chr> <dbl> <lgl>          <dbl> <dbl>     <dbl>     <dbl>
#> 1 1     0.805 TRUE         -0.714      0     0.129     0.930
#> 2 1     0.805 TRUE         -0.916      0     0.217     0.930
#> 3 1     0.805 TRUE         -0.916      0     0.247     0.930
#> 4 1     0.805 TRUE         -0.387      0     0.433     0.930
#> 5 1     0.805 TRUE         -0.0782     0     0.798     0.930
#> 6 2     0.504 TRUE          0.402      1     0.288     0.865
```

The simulated data are stored in long format. Each row corresponds to
one observed longitudinal measurement time for one subject.

The key columns are:

- `id`: subject identifier
- `X`: observed event or censoring time
- `delta`: event indicator
- `covariates`: observed covariate values at that visit time
- `obs_times`: longitudinal observation time

## Fit the General Transformed Hazards Model

The `covariates` column returned by
[`sim_skmle_data()`](https://dayusun.github.io/skmle/reference/sim_skmle_data.md)
is a matrix. You can either use it directly in the formula or split it
into separate columns. Using the matrix directly is convenient for
routine work.

``` r
fit_skmle <- skmle(
  Surv(X, delta) ~ covariates,
  data = dat,
  id = id,
  obs_times = obs_times,
  s = 0,
  h = 0.5,
  nknots = 3,
  norder = 3
)

fit_skmle
```

The printed object gives the fitted coefficients. As in many R model
objects, the formatted inferential output is produced by
[`summary()`](https://rdrr.io/r/base/summary.html).

``` r
summary(fit_skmle)
```

The summary table reports:

- coefficient estimates
- standard errors
- z statistics
- p-values

## Plot the Estimated Baseline Component

``` r
plot(fit_skmle)
```

This plot visualizes the estimated nonparametric baseline component from
the sieve fit.

## Fit the Specialized Estimating-Equation Methods

When the model of interest matches one of the specialized settings, the
package also provides dedicated estimating-equation estimators.

### Cox-Type Estimator

``` r
fit_kee_cox <- kee_cox(
  Surv(X, delta) ~ covariates,
  data = dat,
  id = id,
  obs_times = obs_times,
  h = 0.5
)

summary(fit_kee_cox)
```

### Additive Hazards Estimator

For the additive hazards estimator, simulate data under `s = 1`.

``` r
set.seed(456)

dat_add <- sim_skmle_data(
  n = 200,
  mu = function(tt) 8 * (0.75 + (0.5 - tt)^2),
  mu_bar = 8,
  alpha = function(tt) 0.75 + 0.75 * (tt * (1 - sin(2 * pi * (tt - 0.25)))),
  beta = c(1, -0.5),
  s = 1,
  cen = 0.7
)

fit_kee_add <- kee_additive(
  Surv(X, delta) ~ covariates,
  data = dat_add,
  id = id,
  obs_times = obs_times,
  h = 0.5
)

summary(fit_kee_add)
```

## Select a Bandwidth by Cross-Validation

Bandwidth selection can be handled by
[`skmle_cv()`](https://dayusun.github.io/skmle/reference/skmle_cv.md).

``` r
set.seed(999)

if ("skmle_cv" %in% getNamespaceExports("skmle")) {
  cv_fit <- skmle_cv(
    Surv(X, delta) ~ covariates,
    data = dat,
    id = id,
    obs_times = obs_times,
    s = 0,
    K = 3,
    n_h = 4,
    nknots = 3,
    norder = 3,
    quiet = TRUE
  )

  cv_fit$h_cv
  cv_fit$cv_results
} else {
  cat("`skmle_cv()` is present in the source tree but is not available in the currently installed package namespace used for this render.\n")
}
```

The returned object contains:

- the selected bandwidth
- the refitted `skmle` model at that bandwidth
- the cross-validation loss table

You can then inspect the final refit in the usual way.

``` r
if (exists("cv_fit")) {
  summary(cv_fit$fit)
}
```

## Typical Workflow

For routine use, the usual sequence is:

1.  Prepare data in long format with one row per observation time.
2.  Fit [`skmle()`](https://dayusun.github.io/skmle/reference/skmle.md)
    if you want the general transformed hazards model.
3.  Use [`summary()`](https://rdrr.io/r/base/summary.html) and
    [`plot()`](https://rdrr.io/r/graphics/plot.default.html) to inspect
    the fitted model.
4.  Use
    [`skmle_cv()`](https://dayusun.github.io/skmle/reference/skmle_cv.md)
    if you want data-driven bandwidth selection.
5.  Use
    [`kee_cox()`](https://dayusun.github.io/skmle/reference/kee_cox.md)
    or
    [`kee_additive()`](https://dayusun.github.io/skmle/reference/kee_additive.md)
    when the scientific model matches those specialized estimators.

This gives a standard R model-fitting interface while keeping the core
numerical work in the Rcpp backend.
