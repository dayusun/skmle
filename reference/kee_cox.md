# Fit a Cox-Type KEE Model

Fit the proportional hazards model for sparse longitudinal covariate
data using a kernel estimating-equation approach.

## Usage

``` r
kee_cox(
  formula,
  data,
  id,
  obs_times,
  h = NULL,
  one_sided = TRUE,
  weight = NULL
)
```

## Arguments

- formula:

  A model formula with a
  [`survival::Surv()`](https://rdrr.io/pkg/survival/man/Surv.html)
  response.

- data:

  Data frame containing all variables used in the fit.

- id:

  Subject identifier aligned row-wise with `data`.

- obs_times:

  Longitudinal observation times aligned row-wise with `data`. Times may
  be on any scale, so there is no need to rescale to the unit interval
  first. `h` must be on the same scale. Observation times may be
  negative; see “Negative observation times” in `kee_cox()`.

- h:

  Positive kernel bandwidth. If omitted, one is read off the observation
  times as a rule of thumb and reported in a message. Use
  [`skmle_cv()`](https://www.sundayu.me/skmle/reference/skmle_cv.md) to
  choose it from the data.

- one_sided:

  Logical. `TRUE` (the default) uses a half kernel: only covariate
  observations strictly before the event or quadrature time inform that
  time, which is the risk-set restriction and the estimator as
  published. `FALSE` uses a full, two-sided kernel, smoothing the
  covariate path from both sides. The switch applies to the risk-set
  averages inside the C++ backend as well as to the row weights, so the
  two are always consistent.

- weight:

  Weight function of the standardised lag \\(t - r)/h\\, or `NULL` (the
  default) for the Epanechnikov kernel implied by `one_sided`. Any R
  function will do; see
  [kernel-weights](https://www.sundayu.me/skmle/reference/kernel-weights.md).

## Value

An object of class `kee` containing coefficient estimates, the estimated
variance-covariance matrix, the estimating-equation matrices `W` and
`Sigma`, the per-subject score contributions `scores` (one row per
subject, so that `Sigma` is `crossprod(scores) / n^2`), convergence
status, and the original function call.

## Details

`kee_cox()` targets the proportional hazards case without estimating a
nonparametric baseline component. It is therefore a useful specialized
alternative to
[`skmle()`](https://www.sundayu.me/skmle/reference/skmle.md) when the
scientific model is Cox-type and the main interest is in the regression
coefficients.

## Negative observation times

`obs_times` may be negative. Event times may not, and a negative event
time is an error: in
[`skmle()`](https://www.sundayu.me/skmle/reference/skmle.md) and
[`kee_additive()`](https://www.sundayu.me/skmle/reference/kee_additive.md)
the sieve basis and the cumulative-hazard quadrature are built on
`[0, max(X)]`, and every estimator takes the time origin to be zero.
Observation times carry no such constraint. They enter the score only
through the lag `X - obs_times`, and the half kernel admits a row on the
sign of that lag rather than on the sign of the time, so a covariate
observed at `-180` and an event at `30` sit a lag of `210` apart and are
handled like any other pair.

A pre-origin observation contributes only where it falls within `h` of
an event time. Use negative times for ordinary measurements taken before
entry on the same scale, such as a run-in visit or a screening draw.

## References

Sun, Dayu, Zhuowei Sun, Xingqiu Zhao, and Hongyuan Cao. "Kernel Meets
Sieve: Transformed Hazards Models with Sparse Longitudinal Covariates."
*Journal of the American Statistical Association* 120(552), 2580-2591.
[doi:10.1080/01621459.2025.2476781](https://doi.org/10.1080/01621459.2025.2476781)

Cao, Hongyuan, Mathew M. Churpek, Donglin Zeng, and Jason P. Fine.
"Analysis of the Proportional Hazards Model with Sparse Longitudinal
Covariates." *Journal of the American Statistical Association* 110(511),
1187-1196 (2015).
[doi:10.1080/01621459.2014.957289](https://doi.org/10.1080/01621459.2014.957289)

## Examples

``` r
# \donttest{
library(survival)

set.seed(123)
dat <- sim_skmle_data(
  n = 80,
  mu = function(tt) 8 * (0.75 + (0.5 - tt)^2),
  mu_bar = 8,
  alpha = function(tt) 0.5 * 0.75 + 0.75 * (tt * (1 - sin(2 * pi * (tt - 0.25)))),
  beta = c(1, -0.5),
  s = 0,
  cen = 0.7
)

fit_cox <- kee_cox(
  Surv(X, delta) ~ covariates,
  data = dat,
  id = id,
  obs_times = obs_times,
  h = 0.5
)

summary(fit_cox)
#> Call:
#> kee_cox(formula = Surv(X, delta) ~ covariates, data = dat, id = id, 
#>     obs_times = obs_times, h = 0.5)
#> 
#> Cox-type proportional hazards, kernel estimating equation (half kernel)
#>   n= 80 subjects   bandwidth h = 0.5
#> 
#>             Estimate Std. Error z value Pr(>|z|)   
#> covariates1  1.02085    0.36037  2.8328 0.004615 **
#> covariates2 -0.36751    0.29950 -1.2271 0.219797   
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# }
```
