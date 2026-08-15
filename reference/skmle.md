# Fit a Transformed Hazards Model by SMKLE

Fit a transformed hazards model for survival data with sparsely and
intermittently observed longitudinal covariates using the sieve maximum
kernel-weighted log-likelihood estimator (SMKLE).

## Usage

``` r
skmle(
  formula,
  data,
  id,
  obs_times,
  s = 0,
  h = NULL,
  nknots = 3,
  lq_nodes = 64,
  maxeval = 10000,
  xtol_rel = 1e-06,
  one_sided = TRUE
)
```

## Arguments

- formula:

  A model formula. The left-hand side must be a
  [`survival::Surv()`](https://rdrr.io/pkg/survival/man/Surv.html)
  response and the right-hand side must contain at least one covariate.

- data:

  Data frame containing the variables used in `formula`, `id`, and
  `obs_times`.

- id:

  Subject identifier. Non-numeric identifiers are allowed and are
  internally converted to integer subject codes.

- obs_times:

  Longitudinal observation times aligned row-wise with `data`.

- s:

  Box-Cox transformation parameter, defaulting to `0`. `s = 0` is the
  proportional hazards model, `s = 1` is the additive hazards model, and
  values in between interpolate. If you do not have a reason to choose
  otherwise, the default is the familiar Cox model.

- h:

  Positive kernel bandwidth. If omitted, one is read off the observation
  times as a rule of thumb and reported in a message. Use
  [`skmle_cv()`](https://www.sundayu.me/skmle/reference/skmle_cv.md) to
  choose it from the data.

- nknots:

  Number of interior knots used in the sieve approximation of the
  baseline component. Knots are placed at `(1:nknots)/(nknots + 1)`. The
  basis is a natural cubic spline
  ([`splines::ns`](https://rdrr.io/r/splines/ns.html)); its order is
  fixed, which is why there is no `norder` argument.

- lq_nodes:

  Number of Legendre-Gauss quadrature nodes used in numerical
  integration.

- maxeval:

  Maximum number of optimizer evaluations.

- xtol_rel:

  Relative convergence tolerance passed to the optimizer.

- one_sided:

  Logical. `TRUE` (the default) uses a half kernel: only covariate
  observations strictly before the event or quadrature time inform that
  time, which is the risk-set restriction and the estimator as
  published. `FALSE` uses a full, two-sided kernel, smoothing the
  covariate path from both sides. The switch applies to the risk-set
  averages inside the C++ backend as well as to the row weights, so the
  two are always consistent.

## Value

An object of class `skmle` containing:

- `coefficients`: regression coefficient estimates,

- `var`: estimated variance-covariance matrix for the regression
  coefficients,

- `gamma`: estimated spline coefficients for the baseline component,

- `loglik`: maximized log-likelihood,

- `convergence`: optimizer status code,

- model metadata such as `n`, `s`, `h`, and `call`.

## Details

`skmle()` is the main model-fitting function in the package. It
combines:

- kernel weighting to handle intermittently observed longitudinal
  covariates,

- a sieve approximation for the unknown baseline component, and

- a C++-backed numerical optimizer for the joint estimation problem.

The returned object follows the usual R model pattern: print the fitted
coefficients with [`print()`](https://rdrr.io/r/base/print.html), obtain
inferential output with
[`summary()`](https://rdrr.io/r/base/summary.html), and visualize the
estimated baseline component with
[`plot()`](https://rdrr.io/r/graphics/plot.default.html).

## References

Sun, Dayu, Zhuowei Sun, Xingqiu Zhao, and Hongyuan Cao. "Kernel Meets
Sieve: Transformed Hazards Models with Sparse Longitudinal Covariates."
*Journal of the American Statistical Association* (2025): 1-12.

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

fit <- skmle(
  Surv(X, delta) ~ covariates,
  data = dat,
  id = id,
  obs_times = obs_times,
  s = 0,
  h = 0.5,
  nknots = 3
)

summary(fit)
#> Call:
#> skmle(formula = Surv(X, delta) ~ covariates, data = dat, id = id, 
#>     obs_times = obs_times, s = 0, h = 0.5, nknots = 3)
#> 
#>   n= 80
#> 
#>             Estimate Std. Error z value Pr(>|z|)   
#> covariates1  1.10345    0.34655  3.1841 0.001452 **
#> covariates2 -0.62148    0.35735 -1.7391 0.082015 . 
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Log-likelihood: -0.1486 
# }
```
