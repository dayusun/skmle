# Fit an Additive Hazards KEE Model

Fit the additive hazards model for sparse longitudinal covariate data
using a kernel estimating-equation approach.

## Usage

``` r
kee_additive(formula, data, id, obs_times, h, lq_nodes = 64)
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

  Longitudinal observation times aligned row-wise with `data`.

- h:

  Positive kernel bandwidth.

- lq_nodes:

  Number of quadrature nodes used in the numerical integration step.

## Value

An object of class `kee` containing coefficient estimates, an estimated
variance-covariance matrix, intermediate matrices used for sandwich
variance estimation, and model metadata.

## Details

`kee_additive()` is the specialized additive-hazards counterpart to
[`kee_cox()`](https://dayusun.github.io/skmle/reference/kee_cox.md). It
uses a kernel-smoothed martingale estimating equation and typically runs
faster than the general `skmle(s = 1)` fit because it solves a more
specialized problem.

## References

Sun, Dayu, Zhuowei Sun, Xingqiu Zhao, and Hongyuan Cao. "Kernel Meets
Sieve: Transformed Hazards Models with Sparse Longitudinal Covariates."
*Journal of the American Statistical Association* (2025): 1-12.

Sun, Dayu, Hongyuan Cao, and Yining Chen. "Additive hazards models with
sparse longitudinal covariates." *Lifetime Data Analysis* (2022).

## Examples

``` r
# \donttest{
library(survival)

set.seed(123)
dat <- sim_skmle_data(
  n = 80,
  mu = function(tt) 8 * (0.75 + (0.5 - tt)^2),
  mu_bar = 8,
  alpha = function(tt) 0.75 + 0.75 * (tt * (1 - sin(2 * pi * (tt - 0.25)))),
  beta = c(1, -0.5),
  s = 1,
  cen = 0.7
)

fit_add <- kee_additive(
  Surv(X, delta) ~ covariates,
  data = dat,
  id = id,
  obs_times = obs_times,
  h = 0.5
)

summary(fit_add)
#> Call:
#> kee_additive(formula = Surv(X, delta) ~ covariates, data = dat, 
#>     id = id, obs_times = obs_times, h = 0.5)
#> 
#>   n= 80
#> 
#>             Estimate Std. Error z value Pr(>|z|)  
#> covariates1  1.26747    0.66610  1.9028  0.05706 .
#> covariates2 -0.54733    0.75574 -0.7242  0.46892  
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# }
```
