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
  s,
  h,
  nknots = 3,
  norder = 3,
  lq_nodes = 64,
  maxeval = 10000,
  xtol_rel = 1e-06
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

  Box-Cox transformation parameter. `s = 0` gives the proportional
  hazards model and `s = 1` gives the additive hazards model.

- h:

  Positive kernel bandwidth.

- nknots:

  Number of interior knots used in the sieve approximation of the
  baseline component.

- norder:

  Order parameter supplied to the high-level interface for the spline
  approximation.

- lq_nodes:

  Number of Legendre-Gauss quadrature nodes used in numerical
  integration.

- maxeval:

  Maximum number of optimizer evaluations.

- xtol_rel:

  Relative convergence tolerance passed to the optimizer.

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

fit
#> Call:
#> skmle(formula = Surv(X, delta) ~ covariates, data = dat, id = id, 
#>     obs_times = obs_times, s = 0, h = 0.5, nknots = 3)
#> 
#> Coefficients:
#> covariates1 covariates2 
#>   0.9215573  -0.5463503 
summary(fit)
#> Call:
#> skmle(formula = Surv(X, delta) ~ covariates, data = dat, id = id, 
#>     obs_times = obs_times, s = 0, h = 0.5, nknots = 3)
#> 
#>   n= 80
#> 
#>             Estimate Std. Error z value Pr(>|z|)   
#> covariates1  0.92156    0.31342  2.9403 0.003279 **
#> covariates2 -0.54635    0.34045 -1.6048 0.108543   
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Log-likelihood: -0.0753 

# If you prefer explicit covariate names, split the matrix column first.
dat$Z1 <- dat$covariates[, 1]
dat$Z2 <- dat$covariates[, 2]
fit_named <- skmle(
  Surv(X, delta) ~ Z1 + Z2,
  data = dat,
  id = id,
  obs_times = obs_times,
  s = 0,
  h = 0.5,
  nknots = 3
)

summary(fit_named)
#> Call:
#> skmle(formula = Surv(X, delta) ~ Z1 + Z2, data = dat, id = id, 
#>     obs_times = obs_times, s = 0, h = 0.5, nknots = 3)
#> 
#>   n= 80
#> 
#>    Estimate Std. Error z value Pr(>|z|)   
#> Z1  0.92156    0.31342  2.9403 0.003279 **
#> Z2 -0.54635    0.34045 -1.6048 0.108543   
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Log-likelihood: -0.0753 
```
